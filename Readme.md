# ExcludonFinder
![Untitled (7)](https://github.com/user-attachments/assets/07e51f4a-bec1-4305-9f3e-280db14282e8)

[![Bioconda](https://img.shields.io/conda/vn/bioconda/excludonfinder.svg)](https://bioconda.github.io/recipes/excludonfinder/README.html)

An easy to use tool for identifying and analyzing excludons in genomic data using RNA-seq data.
## Outline 
![rect1](https://github.com/user-attachments/assets/2f5b3e48-0eb2-422f-bf42-ecdc8ef19b9c)

 From a given RNA-seq data, alignment is performed against reference genome (1) and coverage per nucleotide is calculated (2).
Convergent (-> <-) and divergent (<- ->) pairs of genes are substarcted and median covergage is calculated for each of them (3). Trancriptional units (TUs) for each gene is annotated (4) based on gene coverage. A threshold of the covegare decreasing is set, gene gene expression decays under this threshold transcription start and end sites (TSS and TTS) is set. If TUs of convergen and divergent pairs overlaps, this pair is annotated as Excludon (5). 

## Features
- Fast parallel processing for large datasets
- Support for both short and long-read data
- Support for paired-end and single-end RNA-seq data
- Built-in quality checks and mapping statistics


## Installation

### Using Conda (Recommended)
```bash
conda install -c bioconda excludonfinder

### From source
```bash
git clone https://github.com/Alvarosmb/ExcludonFinder.git
cd ExcludonFinder
conda env create -f environment.yml
conda activate ExcludonFinder
```

```markdown
## Usage

### If installed with conda:
```bash
ExcludonFinder -f <reference.fasta> -1 <reads_R1.fastq> -2 <reads_R2.fastq> -g <annotation.gff>

### If installed from source
``` bash
./scripts/ExcludonFinder.sh -f <reference.fasta> -1 <reads_R1.fastq> -2 <reads_R2.fastq> -g <annotation.gff>

### Options

- `-f`: Reference genome in FASTA format
- `-1`: Input FASTQ file for Read 1
- `-2`: Input FASTQ file for Read 2
- `-g`: Annotation file in GFF format
- `-t`: Coverage threshold (default: 0.5)
- `-j`: Number of threads (default: 8)
- `-l`: Long-read data

## Example

```bash
./scripts/ExcludonFinder.sh \
  -f data/example/E.coli_K12_MG1655.fasta \
  -1 data/example/test_R1.fastq \
  -2 data/example/test_R2.fastq \
  -g data/example/E.coli_K12_MG1655.gff \
  -t 0.5 \
  -j 4
```

## Examples 
The data/examples directory contains test RNA-seq data from E. coli K12 MG1655. For faster testing and analysis, the dataset is reduced to reads mapping only to the first 50 genes. Expected results can be found in data/examples/output/.

## Citation
If you found this tool useful, please cite:
```text
Alvaro Sanmartin, Pablo Iturbe, Jeronimo Rodriguez-Beltran, Iñigo Lasa. ExcludonFinder: Mapping Transcriptional Overlaps Between Neighboring Genes

```
