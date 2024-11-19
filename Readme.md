# ExcludonFinder

A powerful tool for identifying and analyzing excludons in genomic data using RNA-seq data.

## Features

- Fast parallel processing for large datasets
- Support for paired-end and single-end RNA-seq data
- Built-in quality checks and mapping statistics
- Support for both short and long-read data

## Installation

```bash
git clone https://github.com/yourusername/ExcludonFinder.git
cd ExcludonFinder
conda env create -f environment.yml
conda activate Excludon_finder_env
```

## Usage

Basic command:
```bash
./scripts/ExcludonFinder.sh -f <reference.fasta> -1 <reads_R1.fastq> -2 <reads_R2.fastq> -g <annotation.gff>
```

### Options

- `-f`: Reference genome in FASTA format
- `-1`: Input FASTQ file for Read 1
- `-2`: Input FASTQ file for Read 2
- `-g`: Annotation file in GFF format
- `-t`: Coverage threshold (default: 0.5)
- `-j`: Number of threads (default: 8)
- `-l`: Use minimap2 for long-read data

## Example

```bash
./scripts/ExcludonFinder.sh \
  -f data/example/EColi_K12_substr._MG1655.fasta \
  -1 data/example/test_R1.fastq \
  -2 data/example/test_R2.fastq \
  -g data/example/test_ecoli.gff \
  -t 0.6 \
  -j 4
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
