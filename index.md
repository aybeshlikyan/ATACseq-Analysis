* TOC {:toc}

## FastQC
```
Version (on macOS Version 10.15.5 (19F101)):
    java version "1.8.0_131"
    Java(TM) SE Runtime Environment (build 1.8.0_131-b11)
    Java HotSpot(TM) 64-Bit Server VM (build 25.131-b11, mixed mode)
    
FastQC v0.11.9 (Mac DMG image)
```

## Bowtie2 and Samtools

I followed the [ATAC-seq Guidelines from Harvard](https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html#alignment) for this portion.

### Version Information
```
Bowtie2 Version (on Hoffman2):
    /u/local/apps/bowtie2/2.2.9/bowtie2-align-s version 2.2.9
    64-bit
    Built on localhost.localdomain
    Thu Apr 21 18:36:37 EDT 2016
    Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-54)
    Options: -O3 -m64 -msse2  -funroll-loops -g3 -DPOPCNT_CAPABILITY
    Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
```
```
SAMtools Version (on Hoffman2):
    samtools 1.9
    Using htslib 1.9
Copyright (C) 2018 Genome Research Ltd.
```

### Generating Index Files with Bowtie2

This only needs to be done once at the beginning of the project.
The FASTA file can be obtained from [Ensembl](http://uswest.ensembl.org/Homo_sapiens/Info/Index).

*generate_index.sh:*
```bash
#!/bin/bash

bowtie2-build grch38_1kgmaj.fa grch38_index
```

### Alignment with Bowtie2 and Samtools

The alignment was run with the maximum DNA fragment length set at 1000 with the option `-X 1000` and the maximum number of alignments to report per read set at 10 with the option `-k 10`.

*align.sh:*
```bash
#!/bin/bash

# This is the environment variable that tells bowtie2 where the directory for the generated index file from the previous step is.
export BOWTIE2_INDEXES=<path to index directory>

bowtie2 --very-sensitive -k 10 -X 1000 -x grch38_index --end-to-end -p 8 \
        -1 <name>_1.fastq.gz \
        -2 <name>_2.fastq.gz \
| samtools view -u \
| samtools sort -o <name>_sorted.bam

samtools index <name>_sorted.bam
```

## ATACseqQC - R Package

I followed the Bioconductor [ATACseqQC Package Guide](https://bioconductor.org/packages/devel/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html) for this portion.

### Version Information

```
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.10 (Final)

ATACseqQC_1.13.6
```


### Generating the TxDb File

The TxDb file only needs to be generated once in the beginning. The GRCh38 annotation file can be found [here](https://www.gencodegenes.org/human/release_29.html). 
To build and save the database:
```R
# building the TxDb object
library(GenomicFeatures)
txdb = makeTxDbFromGFF('gencode.v29.annotation.gtf.gz')

# saving the TxDb database to be used in other files
library(AnnotationDbi)
saveDb(txdb, 'txdb.gencode29.sqlite')
```

The TxDb database can then be loaded in future files with:
```R
# loading the previously generated TxDb file
library(AnnotationDbi)
txdb = loadDb(file = 'txdb.gencode29.sqlite')
```

### Quality Control and Filtering

I have this bash script in the beginning of the ATACseqQC Rmd script to remove the mitochondrial reads, PCR duplicates, and non-unique alignments. I followed [this useful guide by Yiwei Niu](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/#alignment-and-filter).

```bash
PICARD=/u/local/apps/picard-tools/current/picard.jar # This is the path to the Picard jarfile. It might be different on another machine.
bamfile=<name>_sorted.bam

# Remove mitochondrial reads
samtools view -h $bamfile | awk  '($3 != "chrM")' | samtools sort -O bam -o <name>_sorted_rmChrM.bam

# Remove PCR duplicates
java -jar $PICARD MarkDuplicates \
  QUIET=true INPUT=<name>_sorted_rmChrM.bam \
  OUTPUT=<name>_sorted_rmChrM_marked.bam \
  METRICS_FILE=<name>_sorted_Dup.metrics \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp
samtools view -h -b -F 1024 <name>_sorted_rmChrM_marked.bam > <name>_sorted_rmChrM_rmDup.bam

# Remove non-unique alignments
samtools view -h -q 30 <name>_sorted_rmChrM_rmDup.bam > <name>_sorted_rmChrM_rmDup_rmMulti.bam

# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain properly paired reads -f 2
samtools view -h -b -F 1804 -f 2 <name>_sorted_rmChrM_rmDup_rmMulti.bam > <name>_sorted_filtered.bam

# Generate bam index file
samtools index <name>_sorted_filtered.bam
```

## Peak-calling and Blacklist Filtering
### Version Information
```
macs2 2.2.7.1
Python 3.6.3
pip 20.1.1 from /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pip (python 3.6)
```
### Installing MACS2 in a Virtual Environment

```bash
python3 -m venv MACS2-env/
source MACS2-env/bin/activate
pip install macs2
```

### Peak-calling

I originally followed [this useful guide by Yiwei Niu](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/#peak-calling-using-macs2), but I changed the parameters according to [this Twitter thread by Xi Chen](https://twitter.com/XiChenUoM/status/1336658454866325506).

```bash
#!/bin/bash

# convert bam to bed
bedtools bamtobed -i ../ATACseqQC/<name>/<name>_sorted_filtered_shifted.bam \
    > <name>_simple.bed

# start up virtual environment
source ../../MACS2-env/bin/activate

#peak calling  
macs2 callpeak -f BED -g hs --shift -100 --extsize 200 --keep-dup all --cutoff-analysis \
    -n <name> \
    -t <name>_simple.bed --outdir .\
    2> Logs/<name>_macs2.log
```

### Blacklist Filtering

The bedfile with the Encode blacklisted regions can be found [here](https://www.encodeproject.org/files/ENCFF356LFX/).
I followed [this useful guide by Yiwei Niu](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/#blacklist-filtering-for-peaks).

```bash
#!/bin/bash

PEAK=<name>_peaks.narrowPeak
FILTERED_PEAK=<name>_peaks_narrowPeak_filt.bed
BLACKLIST="ENCFF356LFX.bed.gz"

bedtools intersect -v -a ${PEAK} -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[0-9XY]+(?!_)' > ${FILTERED_PEAK}
```

## Annotation with ChIPseeker - R Package

I followed [this guide](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) to run the ChIPseeker package.

## Differential Analysis
### DiffBind - R package

I followed [this manual](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) to conduct differential analysis using the DiffBind package in 2 separate normalization methods: by library size and by reads in peaks.

### CSAW - R package

I followed the CSAW workflow proposed in [this paper](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00342-y) with 2 normalization methods: TMM and Loess.

## Browser Tracks
