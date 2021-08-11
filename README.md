# barcode_detection
Detect and count barcodes in NGS data. Can detect combinations of pre-defined 
sample barcodes and random UMIs (unique molecular barcodes) in sequencing 
reads from Illumina short read or Nanopore long-read sequencing data. Pre-defined
barcodes are detected using a user-defined list of barcode sequences. UMIs
are detected based on known anchor sequences in the reference sequence.

Note: currently barcode_detection.py can only detect pre-defined barcodes of 
length 16 and UMIs of length 12. (To be improved in future releases)


## Installation
Simply clone the repository


## Dependencies
Python 3.8 or higher

Required python packages:
  - argparse
  - Levenshtein

Required tools for preparing the input files for barcode_detection.py
   - bwa-mem
   - samtools


## Preparing the raw sequencing data for use in barcode_detection.py
The following steps have to be performed:
   - indexing of the reference sequence for bwa-mem (only once per reference)
   - alignment of reads against the reference using bwa-mem
   - Sorting of alignments (optional)
   - Indexing of BAM file
   - Quality filtering of alignments (e.g. use only Q60 alignments)
   - Extract reads from high quality alignments as Fasta file 

\
Example for Nanopore reads:
```
bwa index reference.fasta

bwa mem -M -t 4 reference.fasta nanopore_reads.fastq > example.sam

samtools view -b -F 0x900 -S -o example.bam example.sam 

samtools sort -o example.sort.bam example.bam

samtools index example.sort.bam

samtools view -b -q 60 -o example.sort.filtered.bam example.sort.bam

samtools fasta example.sort.filtered.bam > high_quality_reads.fasta
```
\
Example for Illumina reads (using only read-1 fastq file, otherwise same as for Nanopore):
```
bwa index reference.fasta

bwa mem -M -t 4 reference.fasta read1.fastq > example.sam

samtools view -b -F 0x900 -S -o example.bam example.sam 

samtools sort -o example.sort.bam example.bam

samtools index example.sort.bam

samtools view -b -q 60 -o example.sort.filtered.bam example.sort.bam

samtools fasta example.sort.filtered.bam > high_quality_reads.fasta
```


## Alignment QC (optional):
The following steps could be performed to check the fraction of high quality read alignments
   - Alignment statistics with samtools flagstat (shows raw and aligned read counts)
   - Count and compare raw and high-quality-aligned reads

\
Example:
```
samtools flagstat example.bam

wc -l high_quality_reads.fasta   # divide by 2!
```

## Input files and parameters for barcode_detection.py

Required input files:
   - high_quality_reads.fasta (see above)
   - list of fixed barcodes to be detected

Format of pre-defined barcode list: tab-delimited file (no header) with two columns for:
   - barcode-name
   - barcode-sequence

The number of pre-defined barcodes is not limited, but only barcodes of 
length 16 can be identified in the current version.

Parameters:
   - reference fasta file (-f)
   - list of pre-defined barcodes (-b) 
   - sample name (-s): this name will appear as first column of the output and can be any word
   - experiment name (-e): this name will appear as second column in the output and can be any word

Note: parameter for edit distance will be added soon (currently hard-coded in read_parser.py)


## Running barcode_detection.py
```
python barcode_detection.py -f high_quality_reads.fasta -b barcode_list.txt -s test_sample -e test_experiment > barcode_counts.tsv
```

## Output of barcode_detection.py
barcode_detection.py generates a single output file that contains 5 columns (tab-delimited with header line):
   - Sample name
   - Experiment name
   - Pre-defined barcode name
   - Pre-defined barcode sequence
   - Total number of reads detected having the respective pre-defined barcode
   - Count of unique UMI sequences detected together with the respective pre-defined barcode

