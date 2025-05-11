# Genome Assembly Algorithms

This repository contains implementations of genome assembly algorithms as part of the Bioinformatics Algorithms Assignment 2: Genome Assembly and Evaluation.

## Overview

The project consists of two main tasks:

### Task 1: Basic Genome Assembly Algorithms
1. De Bruijn Graph (DBG) assembly algorithm
2. Overlap-Layout-Consensus (OLC) assembly algorithm

Both implementations are designed to process genomic reads from FASTQ files and output assembled contigs in FASTA format. The code includes evaluation metrics for assembly quality assessment and visualization capabilities.

### Task 2: Lizard Assembly Pipeline
A pipeline for assembling and evaluating the Scincus mitranus (sandfish lizard) genome using high-performance computing resources on the Ibex cluster. The pipeline leverages professional assembly tools like Hifiasm and includes comprehensive evaluation of assembly quality.

## Prerequisites

### For Task 1 (Basic Assembly Algorithms)
- Python 3.7+
- BioPython
- NetworkX
- matplotlib
- numpy
- QUAST (http://quast.sourceforge.net/quast) for assembly evaluation
- Bandage (https://rrwick.github.io/Bandage/) for assembly graph visualization
- SPAdes (https://github.com/ablab/spades) for comparison with professional DBG assembler
- Canu (https://github.com/marbl/canu) for comparison with professional OLC assembler
- SRA Toolkit (https://github.com/ncbi/sra-tools) for downloading SRA data
- seqtk (https://github.com/lh3/seqtk) for subsampling FASTQ files

### For Task 2 (Lizard Assembly)
- Access to the Ibex cluster or equivalent high-performance computing resources
- Hifiasm - for high-quality long-read assembly
- YAHS - Hi-C scaffolding tool
- Samtools - for handling alignment files
- Minimap2 - for read alignment
- R with ggplot2 - for visualization
- BUSCO - for gene completeness assessment
- Asmgene - alternative gene completeness tool
- Merqury - for k-mer distribution and QV score
- QUAST - for assembly evaluation
- Inspector - for misassembly identification

## Project Structure

```
.
├── data/                               # From assignment
│   ├── toy_dataset/                    # Small datasets for testing
│   │   ├── reads_r.fastq
│   │   ├── reads_b.fastq
│   │   ├── reference_r.fasta
│   │   └── reference_b.fasta
│   └── synthetic_dataset/              # MERS-CoV datasets
│       ├── GCF_000901155.1_ViralProj183710_genomic.fna
│       └── reads/
│           ├── no_error_reads_hiseq_5k.fastq
│           ├── reads_hiseq_5k.fastq
│           ├── no_error_ont_hq_50x.fastq
│           └── ont_hq_50x.fastq
├── dbg_assembly/                       # De Bruijn Graph implementation
│   ├── main.py
│   └── results/
├── olc_assembly/                       # Overlap-Layout-Consensus implementation
│   ├── main.py
│   └── results/
├── lizard_assembly_pipeline.py         # Script for Task 2 genome assembly
├── lizard_assembly_results/            # Results from lizard genome assembly (on Ibex)
│   ├── assembly/
│   ├── evaluation/
│   │   ├── quast_results/
│   │   ├── busco_results/
│   │   ├── merqury_results/
│   │   └── misassembly_results/
```

## De Bruijn Graph (DBG) Assembly

The DBG assembly algorithm creates a graph where nodes represent k-mers and edges represent overlaps between k-mers, then identifies Eulerian paths to construct contigs.

### Usage

```bash
python dbg_assembly/main.py \
  --fastq <path_to_fastq> \
  --kmer <kmer_size> \
  --output <output_fasta> \
  --gfa <output_gfa> \
  --min-length <min_contig_length> \
  --verbose
```

Parameters:
- `--fastq`: Path to the input FASTQ file
- `--kmer`: k-mer size for graph construction (default: 31)
- `--output`: Path to the output FASTA file
- `--gfa`: Path to output GFA format file for graph visualization
- `--min-length`: Minimum contig length to report (default: 100)
- `--verbose`: Enable verbose output

## Overlap-Layout-Consensus (OLC) Assembly

The OLC assembly algorithm computes all-vs-all read overlaps, constructs an overlap graph, identifies non-branching paths, and generates a consensus sequence for each contig.

### Usage

```bash
python olc_assembly/main.py \
  --fastq <path_to_fastq> \
  --min-overlap <min_overlap_length> \
  --output <output_fasta> \
  --min-length <min_contig_length> \
  --error-rate <max_error_rate> \
  --verbose
```

Parameters:
- `--fastq`: Path to the input FASTQ file
- `--min-overlap`: Minimum overlap length for graph construction (default: 50)
- `--output`: Path to the output FASTA file
- `--min-length`: Minimum contig length to report (default: 100)
- `--error-rate`: Maximum error rate for overlap detection (default: 0.05)
- `--verbose`: Enable verbose output

## Example Commands

### Task 1.3.1: Visualizing De Bruijn Graph

```bash
# Generate assembly and graph file for dataset B with k=40
python dbg_assembly/main.py \
  --fastq data/toy_dataset/reads_b.fastq \
  --kmer 40 \
  --output dbg_assembly/results/toy_b_k40.fasta \
  --gfa dbg_assembly/results/toy_b_k40.gfa \
  --verbose

# Visualize the graph using Bandage
Bandage.app/Contents/MacOS/Bandage load dbg_assembly/results/toy_b_k40.gfa
```

### Task 1.3.2: Comparing Different k-mer Sizes

```bash
# Generate assembly for dataset R with k=35
python dbg_assembly/main.py \
  --fastq data/toy_dataset/reads_r.fastq \
  --kmer 35 \
  --output dbg_assembly/results/toy_r_k35.fasta \
  --gfa dbg_assembly/results/toy_r_k35.gfa \
  --verbose

# Generate assembly for dataset R with k=45
python dbg_assembly/main.py \
  --fastq data/toy_dataset/reads_r.fastq \
  --kmer 45 \
  --output dbg_assembly/results/toy_r_k45.fasta \
  --gfa dbg_assembly/results/toy_r_k45.gfa \
  --verbose

# Visualize both graphs
Bandage.app/Contents/MacOS/Bandage load dbg_assembly/results/toy_r_k35.gfa
Bandage.app/Contents/MacOS/Bandage load dbg_assembly/results/toy_r_k45.gfa

# Compare the assemblies using QUAST
quast.py \
  dbg_assembly/results/toy_r_k35.fasta \
  dbg_assembly/results/toy_r_k45.fasta \
  -r data/toy_dataset/reference_r.fasta \
  -o quast_toy_r \
  --threads 1
```

### Task 1.3.3: Assembling MERS Virus Reads

```bash
# DBG Assembly for error-free Illumina reads
python dbg_assembly/main.py \
  --fastq data/synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq \
  --kmer 45 \
  --output dbg_assembly/results/dbg_mers_errorless_illumina.fasta \
  --min-length 100 \
  --verbose

# DBG Assembly for Illumina reads with errors
python dbg_assembly/main.py \
  --fastq data/synthetic_dataset/reads/reads_hiseq_5k.fastq \
  --kmer 45 \
  --output dbg_assembly/results/dbg_mers_error_illumina.fasta \
  --min-length 100 \
  --verbose

# DBG Assembly for error-free ONT reads
python dbg_assembly/main.py \
  --fastq data/synthetic_dataset/reads/no_error_ont_hq_50x.fastq \
  --kmer 45 \
  --output dbg_assembly/results/dbg_mers_errorless_ont.fasta \
  --min-length 100 \
  --verbose

# DBG Assembly for ONT reads with errors
python dbg_assembly/main.py \
  --fastq data/synthetic_dataset/reads/ont_hq_50x.fastq \
  --kmer 45 \
  --output dbg_assembly/results/dbg_mers_error_ont.fasta \
  --min-length 100 \
  --verbose

# OLC Assembly for error-free Illumina reads
python olc_assembly/main.py \
  --fastq data/synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq \
  --min-overlap 50 \
  --output olc_assembly/results/olc_mers_errorless_illumina.fasta \
  --min-length 100 \
  --verbose

# OLC Assembly for Illumina reads with errors
python olc_assembly/main.py \
  --fastq data/synthetic_dataset/reads/reads_hiseq_5k.fastq \
  --min-overlap 50 \
  --output olc_assembly/results/olc_mers_error_illumina.fasta \
  --min-length 100 \
  --error-rate 0.1 \
  --verbose

# OLC Assembly for error-free ONT reads
python olc_assembly/main.py \
  --fastq data/synthetic_dataset/reads/no_error_ont_hq_50x.fastq \
  --min-overlap 500 \
  --output olc_assembly/results/olc_mers_errorless_ont.fasta \
  --min-length 100 \
  --verbose

# OLC Assembly for ONT reads with errors
python olc_assembly/main.py \
  --fastq data/synthetic_dataset/reads/ont_hq_50x.fastq \
  --min-overlap 500 \
  --output olc_assembly/results/olc_mers_error_ont.fasta \
  --min-length 100 \
  --error-rate 0.25 \
  --verbose

# Evaluate all assemblies using QUAST
quast.py \
  dbg_assembly/results/dbg_mers_*.fasta \
  olc_assembly/results/olc_mers_*.fasta \
  -r data/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna \
  -o quast_mers \
  -t 1
```

### Task 1.3.4: Comparison with SPAdes Assembler

```bash
# SPAdes Assembly for error-free Illumina reads
spades.py \
  -s data/synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq \
  -o results/spades_mers_errorless_illumina \
  --isolate \
  -m 500

# SPAdes Assembly for Illumina reads with errors
spades.py \
  -s data/synthetic_dataset/reads/reads_hiseq_5k.fastq \
  -o results/spades_mers_error_illumina \
  --isolate \
  -m 500

# SPAdes Assembly for error-free ONT reads
spades.py \
  -s data/synthetic_dataset/reads/no_error_ont_hq_50x.fastq \
  -o results/spades_mers_errorless_ont \
  --isolate \
  -m 500

# SPAdes Assembly for ONT reads with errors
spades.py \
  -s data/synthetic_dataset/reads/ont_hq_50x.fastq \
  -o results/spades_mers_error_ont \
  --isolate \
  -m 500

# Compare SPAdes assemblies using QUAST
quast.py \
  results/spades_mers_errorless_illumina/contigs.fasta \
  results/spades_mers_error_illumina/contigs.fasta \
  results/spades_mers_errorless_ont/contigs.fasta \
  results/spades_mers_error_ont/contigs.fasta \
  -r data/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna \
  -o quast_mers_comparison_spades \
  --threads 1
```

## Assembly Evaluation Metrics

For all assemblies, the following metrics are calculated:

- **Sequence length**: Total assembly length
- **Number of contigs**: Total number of contigs
- **GC content (%)**: GC percentage
- **Genome fraction (%)**: Percentage of reference covered by assembly
- **Duplication ratio**: Average number of times a reference base is covered
- **Largest contig**: Length of the largest contig
- **N50**: Length where 50% of assembly is in contigs of this size or larger
- **N90**: Length where 90% of assembly is in contigs of this size or larger
- **L50**: Number of contigs to reach 50% of total assembly length
- **Misassemblies**: Number of positions with breakpoints relative to reference
- **Mismatches per 100 kbp**: Number of mismatches per 100,000 bases
- **Indels per 100 kbp**: Number of insertions/deletions per 100,000 bases

## Implementation Details

### De Bruijn Graph Assembly

The DBG implementation:
1. Constructs a de Bruijn graph from k-mers extracted from input reads
2. Simplifies the graph by merging non-branching paths
3. Identifies Eulerian paths to generate contigs
4. Outputs the assembled contigs as a FASTA file
5. Optionally outputs the graph structure in GFA format for visualization

### Overlap-Layout-Consensus Assembly

The OLC implementation:
1. Computes all-vs-all read overlaps using a prefix-suffix comparison
2. Constructs an overlap graph with reads as nodes and overlaps as edges
3. Identifies non-branching paths in the graph to form contigs
4. Generates a consensus sequence for each contig
5. Outputs the assembled contigs as a FASTA file

## Results

The repository includes example outputs from running the assembly algorithms on various datasets, including:
- Toy synthetic datasets with and without errors
- MERS-CoV (Middle East Respiratory Syndrome Coronavirus) genome
- Comparison with professional assemblers like SPAdes

## Task 2: Lizard Assembly Pipeline

For the second task, we implement a realistic de novo assembly pipeline for the Scincus mitranus (sandfish lizard) genome using HPC resources.

### Prerequisites

The pipeline requires access to the Ibex cluster and the following modules:
- hifiasm
- yahs
- samtools
- minimap2
- R
- busco
- merqury
- quast

### SLURM Job Submission

Create a SLURM job script (`submit_lizard_assembly.sh`) with the following content:

```bash
#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --constraint=rome 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=475G
#SBATCH --job-name=cs249-liz
#SBATCH --mail-type=ALL
#SBATCH --output=bin/logs/%x-%j-slurm.out
#SBATCH --error=bin/logs/%x-%j-slurm.err
#SBATCH --account=cs249

# entire script fails if single command fails
set -e

# load required modules
module purge
module avail
module load hifiasm yahs samtools minimap2 R busco merqury quast

# Run the assembly pipeline
python lizard_assembly_pipeline.py \
  --input-dir /ibex/reference/course/cs249/lizard/input/ncbi_upload \
  --outdir ./lizard_assembly_results \
  --threads 128  \
  --skip-assembly \
  --assembly-path ./lizard_assembly_results/assembly/sandfish.hic.p_ctg.fasta \
  --force
```

Submit the job using:

```bash
sbatch submit_lizard_assembly.sh
```

### Pipeline Components

The `lizard_assembly_pipeline.py` script automates the following steps:

1. **Genome Assembly**:
   - Uses Hifiasm assembler optimized for PacBio HiFi reads
   - Incorporates Hi-C data for chromosome-level scaffolding
   - Produces high-quality primary contigs

2. **Assembly Evaluation**:
   - QUAST for basic assembly metrics
   - BUSCO for gene completeness assessment
   - Merqury for k-mer distribution and QV score
   - Misassembly identification using alignment-based methods

### Usage Options

```bash
python lizard_assembly_pipeline.py \
  --input-dir <input_directory> \
  --outdir <output_directory> \
  --threads <number_of_threads> \
  [--skip-assembly] \
  [--assembly-path <path_to_existing_assembly>] \
  [--force]
```

Parameters:
- `--input-dir`: Directory containing the input sequencing files
- `--outdir`: Directory to store assembly results
- `--threads`: Number of CPU threads to use
- `--skip-assembly`: Skip the assembly step (if you already have an assembly)
- `--assembly-path`: Path to an existing assembly (used with --skip-assembly)
- `--force`: Overwrite existing output files

### Expected Results

After running the pipeline, the following metrics will be generated:

1. **Basic Assembly Metrics** (QUAST):
   - Total assembly length
   - Number of contigs/scaffolds
   - N50/L50 statistics
   - GC content

2. **Gene Completeness** (BUSCO):
   - Percentage of complete, fragmented, and missing genes
   - Comparison to vertebrate lineage database

3. **Assembly Quality** (Merqury):
   - QV (Quality Value) score
   - k-mer completeness
   - k-mer spectrum analysis

4. **Misassembly Detection**:
   - Identification of potential structural errors
   - Alignment discrepancies against reference databases

