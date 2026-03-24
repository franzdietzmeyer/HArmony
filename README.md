```text
██   ██  █████  ██████  ███    ███  ██████  ███    ██ ██    ██ 
██   ██ ██   ██ ██   ██ ████  ████ ██    ██ ████   ██  ██  ██  
███████ ███████ ██████  ██ ████ ██ ██    ██ ██ ██  ██   ████   
██   ██ ██   ██ ██   ██ ██  ██  ██ ██    ██ ██  ██ ██    ██    
██   ██ ██   ██ ██   ██ ██      ██  ██████  ██   ████    ██    
                                                                                                                 
```

Universal Influenza Hemagglutinin Numbering, Antigenic Annotation, and Comparative Evolution Engine.

## Overview

Mapping highly divergent Influenza A Hemagglutinin (HA) subtypes (e.g., H5, H7, H9) onto a universal coordinate framework is inherently challenging. Lineage-specific insertions and deletions (indels) continuously shift the apparent topological positions of critical amino acids, breaking traditional sequence alignments.

HArmony resolves this sequence-to-structure alignment problem by deploying pre-trained Profile Hidden Markov Models (pHMMs). This architecture anchors highly variable input sequences to a standardized structural grid based on the mature H3 (Aichi/68) nomenclature. This approach yields robust, alignment-independent coordinate mapping, enabling high-resolution comparative immunology and evolutionary tracking across the entire viral phylogeny.

## Key Capabilities

### Universal Structural Grid:
Transposes any input HA sequence into H3-standardized coordinates, leveraging HMMER-derived posterior probabilities to assign per-residue confidence scores.

### Antigenic Topology Mapping:
Automatically annotates structural domains (HA1/HA2) and maps residues against the canonical Burke and Smith (1994) antigenic sites (A-E) on the globular head.

### Glycan Shield Profiling:
Detects classical N-linked glycosylation sequons (N-X-[S/T], where X != P) and computes de novo gain/loss dynamics relative to the reference state or a paired sequence.

### Pathogenicity Motif Diagnostics:
Scans the HA1/HA2 junction for polybasic cleavage site (PBCS) expansions and furin-recognition motifs to robustly stratify Highly Pathogenic (HPAI) versus Low Pathogenic (LPAI) signatures.

### Evolutionary Compare Module:
Executes coordinate-level pairwise comparative analyses to precisely quantify structural drift, receptor-binding site (RBS) mutations, and shifts in glycan shielding.

### Publication-Grade Visualizations:
Automates the generation of linear protein feature maps and stacked mirror alignments using matplotlib, detailing mutational hotspots directly on the structural topology.

## Installation

HArmony relies on the HMMER suite (`hmmsearch`) for probabilistic modeling. We highly recommend using Conda/Mamba for installation, as it automatically resolves all Python and C-binary dependencies via Bioconda.

```bash
# 1. Clone the repository
git clone https://github.com/franzdietzmeyer/harmony.git
cd HArmony

# 2. Create the isolated Conda environment
conda env create -f environment.yml

# 3. Activate the environment
conda activate harmony-env
```

Note: The environment installation handles the HMMER binaries automatically. No manual `$PATH` configuration is required.

## Usage Guide

### 1. Batch Mapping & Annotation (`map`)

Process single sequences or large FASTA libraries to generate standardized coordinate grids, antigenic labels, and diagnostic statistics.

```bash
harmony map --input data/my_sequences.fasta --output-dir ./results --ref H3 --plot
```

Generates:

- Detailed per-sequence coordinate matrices (CSV/JSON).
- Run-level summary report (aggregating HPAI/LPAI status and glycan counts).
- High-resolution linear structural maps (`.png`).

### 2. Evolutionary Pairwise Analysis (`compare`)

Perform a coordinate-locked comparison between two sequences (e.g., an ancestral strain vs. a newly emerged outbreak isolate) to track molecular drift.

```bash
harmony compare --seq1 data/ancestral_h5.fa --seq2 data/emerged_h5.fa --ref H3 --plot
```

Generates:

- A `comparison_report.csv` detailing the specific residue substitutions, topological locations, and glycan shifts.
- A mirrored, stacked visualization highlighting mutational connectors between the two viral states.

## Understanding the Output

The core output of HArmony is the coordinate matrix. Below is an example of the generated data structure:

```text
H3_Coord  Residue  Domain  Antigenic_Site  Glycan_Status      Pathogenicity      Confidence
156       K        HA1     Site B          None               -                 0.99
158       N        HA1     Site B          N-Glyc (Gained)    -                 0.98
329a      R        HA1     None            None               HPAI (Polybasic)  0.95
```

## Citation

If HArmony facilitates your research, please cite:

Dietzmeyer, F. (2026). HArmony: Universal Influenza Hemagglutinin Numbering and Antigenic Site Mapping. (Manuscript in preparation). GitHub repository: https://github.com/franzdietzmeyer/harmony

## License

This project is licensed under the **Open Science Protector (GNU GPLv3)** license. See `LICENSE` for the full text and terms.