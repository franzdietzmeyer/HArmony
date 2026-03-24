```text
    __  _____                                 
   / / / /   |  _________ ___  ____  ____  __  __
  / /_/ / /| | / ___/ __ `__ \/ __ \/ __ \/ / / /
 / __  / ___ |/ /  / / / / / / /_/ / / / / /_/ / 
/_/ /_/_/  |_/_/  /_/ /_/ /_/\____/_/ /_/\__, /  
                                        /____/
```

![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Release](https://img.shields.io/badge/release-v1.1.0--alpha-purple.svg)

**HArmony: Universal Influenza HA Numbering, Antigenic Annotation, and Comparative Evolution Engine.**

## Overview

Mapping diverse Influenza A Hemagglutinin (HA) subtypes (for example H5 or H7) onto a universal coordinate framework is difficult because insertions and deletions shift apparent positions between strains and clades.

HArmony addresses this problem by using a pre-trained Profile Hidden Markov Model (HMM) to align each sequence onto a common structural grid (H3/Aichi/68-compatible numbering). This provides robust, alignment-independent coordinate mapping that remains stable across subtype diversity.

## Key Features

- **Universal Grid**: Maps any HA sequence to H3-compatible coordinates with HMMER-derived confidence scores.
- **Antigenic Intelligence**: Annotates residues with Burke and Smith (1994) antigenic sites (A-E) and structural domains (HA1/HA2).
- **Glycan Tracking**: Detects N-linked glycosylation motifs (`N-X-[S/T]`) and reports de novo gain/loss behavior.
- **Pathogenicity Diagnostics**: Evaluates HA1/HA2 cleavage region polybasic signatures and furin-like motifs to classify HPAI-like vs LPAI-like patterns.
- **Evolutionary Compare Mode**: Performs pairwise coordinate-level comparison to track structural and antigenic drift.
- **Publication Plots**: Generates linear protein feature maps and comparison visualizations using matplotlib.

## Installation (Conda Recommended)

```bash
git clone https://github.com/yourusername/harmony.git
cd harmony
conda env create -f environment.yml
conda activate harmony-env
```

Conda installs the required HMMER binaries (including `hmmsearch`) through Bioconda, so no manual binary path configuration is needed.

## Usage Guide

### 1) `map` - Batch Processing

```bash
harmony map --input my_sequences.fasta --output-dir ./results --ref H3 --plot
```

This command generates:
- detailed per-residue mapping output (CSV/JSON depending on options),
- a run-level summary CSV,
- linear structural map PNG files when `--plot` is enabled.

### 2) `compare` - Pairwise Analysis

```bash
harmony compare --seq1 old_h5.fa --seq2 new_h5.fa --ref H3 --plot
```

This produces:
- a coordinate-level comparison report (`comparison_report.csv`),
- a mirrored stacked comparison figure highlighting mutations and glycan shifts.

## Understanding the Output

Example mapping columns:

| Column | Meaning |
|---|---|
| `H3_Coord` | Reference coordinate assigned by HArmony |
| `Residue` | Amino acid at mapped position |
| `Site` | Antigenic site annotation (A-E or None) |
| `Glycan` | Glycosylation status annotation |
| `Confidence` | HMM alignment confidence (0.0-1.0) |

## Citation

If you use HArmony in your research, please cite:

**Dietzmeyer, F. (2026). HArmony: Universal HA Numbering and Antigenic Site Mapping.**
