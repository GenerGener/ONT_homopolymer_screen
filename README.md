# ONT_homopolymer_screen
Identify regions in nucleic acid sequences which may coincide with sequencing artifacts and contribute toward incorrect interpretations. This script was originally written with nanopore sequencing in mind, but other platforms have historically had issues with homopolymers. Additionally, repetitive sequences may be more likely to change between replication cycles.

Related tool: [telomere-finder.py](https://github.com/GenerGener/telomere)

# ONT Problematic Regions Detection Tutorial

This tutorial demonstrates how to use the ONT problematic regions detection script "ont-problems.py" to identify potentially problematic sequences for Oxford Nanopore Technologies (ONT) sequencing, including homopolymers, repeats, and user-defined k-mers. ONT wet and dry lab development has advanced at a rapid pace. This repo is designed to be used to help inform analyses on legacy ONT datasets, and to provide users with regions to prioritize during standard bioinformatics workflows.

Note as this tool may be under continued development, exact naming and usage might change. Usage should be based on user preference/local running environment.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/ont-problems
cd ont-problems

# Install required packages
pip install biopython pandas
```

## Input Data

For this tutorial, we'll use two example files:
1. A demo genome with telomeric repeats (`demo-genome.fasta`)
2. The HIV-1 HXB2 reference genome

### Demo Genome Structure
The demo genome contains three chromosomes:
```fasta
>chr1
ACTGACTGACTGACTGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTAGGGTTAGGGACT
GACTGACTGACTGACTGACTGACTGACTGACTGACTGTTAGGGTTAGGGTTAGGG

>chr2
ACTGACTGACTGACTGTTAGGGTTAGGGTTAGGGACTGACTGACTGACTGTTAGGGTTAGGG
CCCTAACCCTAACCCTAACCCTAA

>chr3_reverse
CCCTAACCCTAACCCTAACCCTAACCCTAAGCTGACTGACTGACTGACTGCCCTAACCCTAA
CCCTAACCCTAACCCTAA
```

This file contains both forward ("TTAGGG") and reverse ("CCCTAA") telomeric repeats, making it perfect for testing our k-mer detection.

## Basic Usage

The script can be run with minimal parameters:

```bash
python ont_problems.py --input demo-genome.fasta --output-prefix results/demo
```

This will:
1. Identify homopolymers and repeats
2. Find default k-mers (6-mers and 7-mers)
3. Generate BED*, CSV, and metadata files *TODO check BED/GTF/GFF; check basing as in telomere-finder

## Advanced Usage

### Custom K-mer Lengths

```bash
python ont_problems.py --input demo-genome.fasta \
    --output-prefix results/demo \
    --kmer-lengths 6 7 8
```

### Different Output Formats

```bash
# GFF output (1-based coordinates)
python ont_problems.py --input demo-genome.fasta \
    --output-prefix results/demo \
    --output-format gff

# GTF output (1-based coordinates)
python ont_problems.py --input demo-genome.fasta \
    --output-prefix results/demo \
    --output-format gtf
```

### Compressed Input

The script automatically handles gzipped FASTA files:

```bash
gzip demo-genome.fasta
python ont_problems.py --input demo-genome.fasta.gz \
    --output-prefix results/demo
```

## Example Analysis: HIV-1 HXB2

Let's analyze the HIV-1 HXB2 reference genome for potential ONT sequencing issues:

```bash
python ont_problems.py --input HXB2.fasta \
    --output-prefix results/hiv \
    --kmer-lengths 6 7
```

Expected output for HIV-1:
```
Top k-mers at contig ends:
K03455.1:
  TTAGCC: 12
  AAGCTT: 10
  TTTGCC: 8
  GCCTGT: 8
  GACTGG: 7

Problem region counts by type:
  homopolymer: 42
  dimer_repeat: 28
  trimer_repeat: 15
  kmer: 1834
```

## Understanding the Output

### 1. Region Files (BED/GFF/GTF)
- BED format (0-based):
  ```
  chr1    0    4    homopolymer_AAAA    4    +
  chr1    10   16   kmer_TTAGGG    6    +
  ```

- GFF format (1-based):
  ```
  ##gff-version 3
  chr1    ONT_problems    homopolymer    1    4    .    +    .    ID=homopolymer_0;Sequence=AAAA
  ```

### 2. Problems CSV
Contains detailed information about all identified regions:
- Contig name
- Start/end positions
- Sequence
- Type (homopolymer/repeat/kmer)
- Length
- Strand

### 3. Metadata File
Includes:
- Run timestamp
- Total regions found
- Counts by problem type
- Most common k-mers
- K-mers near contig ends

## Best Practices

1. **K-mer Selection**:
   - Use k=6 for telomeric repeats (TTAGGG)
   - Use k=7 for other common repeats
   - Larger k
