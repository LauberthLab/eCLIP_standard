# eCLIP-seq Analysis Pipeline

A Snakemake pipeline for end-to-end processing of single-end enhanced CLIP (eCLIP) sequencing data, from raw FASTQ files through differential enrichment analysis, strand-specific peak calling, and motif discovery. Designed for execution on the Northwestern Quest HPC cluster (SLURM, allocation `p32170`).

---

## Overview

The pipeline performs 16 sequential steps:

| Step | Tool | Description |
|------|------|-------------|
| 1 | FastQC | Quality check on raw reads |
| 2 | Trim Galore | Illumina adapter trimming + quality filtering (Q ≥ 28) |
| 3 | STAR | Alignment to hg38 (EndToEnd, unique mappers only) |
| 4 | samtools | BAM indexing |
| 5 | featureCounts | Strand-specific gene-level read counting (all samples jointly) |
| 6 | process_fc.py | Tidy the raw counts matrix (clean column names, drop annotation columns) |
| 7 | DESeq2 | Differential enrichment: IP vs IgG (default) or IP vs Input |
| 8 | samtools merge | Merge replicate BAMs per condition × sample type |
| 9 | bamCompare | Log2-ratio bigWig generation (individual replicates + merged) |
| 10 | plotFingerprint | Enrichment QC (IP vs control cumulative read distribution) |
| 11 | MultiQC | Aggregated QC report across all steps |
| 12 | samtools | Split aligned BAMs into plus- and minus-strand |
| 13 | samtools merge | Merge IgG replicates per strand per condition (peak calling control) |
| 14 | MACS2 | Peak calling — individual IP replicates vs merged IgG, per strand |
| 15 | MACS2 | Peak calling — pooled IP replicates vs merged IgG, per strand |
| 16 | HOMER | De novo motif discovery on every narrowPeak file (`-rna` mode) |

---

## Repository Structure

```
.
├── Snakefile                 # Main workflow definition (all 16 rules)
├── config.yaml               # Pipeline parameters (paths, thresholds, tool options)
├── config.yaml (profile)     # SLURM profile config (sbatch template, default resources)
├── metadata.tsv              # Sample sheet (sample ↔ condition ↔ type ↔ replicate ↔ FASTQ)
└── scripts/
    ├── process_fc.py         # Tidy featureCounts output
    └── de_analysis.R         # DESeq2 differential enrichment analysis
```

---

## Requirements

### Software

The following tools must be available in `$PATH` (e.g., via `module load` on Quest):

- **Snakemake** ≥ 7.32 (recommended: install via mamba into a dedicated conda env)
- **FastQC**
- **Trim Galore** (includes Cutadapt)
- **STAR** (with a pre-built hg38 genome index)
- **samtools**
- **Subread** (featureCounts)
- **deepTools** (bamCompare, plotFingerprint)
- **MACS2** (or MACS3; set `macs_cmd` in config)
- **HOMER** (findMotifsGenome.pl; hg38 genome must be installed)
- **MultiQC**
- **R** with packages: `DESeq2`, `edgeR`, `optparse`
- **Python 3** with `pandas`

### Reference Files

| File | Config Key | Description |
|------|-----------|-------------|
| STAR genome index | `genome_dir` | Pre-built hg38 STAR index directory |
| Gencode GTF | `gtf` | Gene annotation (e.g., gencode.v45.annotation.gtf) |

---

## Metadata Format

The pipeline is driven by a tab-separated `metadata.tsv` file with the following required columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier (used as wildcard throughout) |
| `condition` | Experimental condition grouping (e.g., `TOP1`) |
| `type` | One of `IP`, `IgG`, or `Input` |
| `replicate` | Integer replicate number |
| `fastq` | Absolute path to the raw `.fastq.gz` file |

Example:

```
sample              condition  type   replicate  fastq
21161FL-30-01_S1    TOP1       IP     1          /path/to/sample_R1.fastq.gz
21161FL-30-05_S5    TOP1       IgG    1          /path/to/igg_R1.fastq.gz
21161FL-30-03_S3    TOP1       Input  1          /path/to/input_R1.fastq.gz
```

The pipeline automatically derives all pairwise comparisons (IP vs IgG, IP vs Input) and merging groups from this file. Peak calling (steps 12–16) uses IgG as the control and is only run for conditions where both IP and IgG samples exist.

---

## Configuration

All tunable parameters live in `config.yaml`:

### Paths

- `metadata` — path to the sample sheet
- `genome_dir` — STAR genome index directory
- `gtf` — Gencode annotation GTF
- `results_dir` — top-level output directory (created automatically)

### Differential Enrichment (DESeq2)

- `de_control` — which type to use as the baseline: `"IgG"` (default) or `"Input"`
- `min_cpm` — CPM filter threshold (genes with CPM > this value in ≥ 2 samples are retained; default `0.3`)
- `lfc_threshold` — minimum log2 fold-change for significance (default `0.5849`, i.e., ~1.5-fold)
- `padj_threshold` — adjusted p-value cutoff (default `0.05`)

### BigWig Generation (bamCompare)

- `bin_size` — genomic bin size in bp (default `1`)
- `normalize_using` — normalization method: `RPKM`, `CPM`, `BPM`, `RPGC`, or `None`
- `pseudocount` — pseudocount added before log2 ratio (default `1`)

### Peak Calling (MACS2)

- `macs_cmd` — executable name (`macs2` or `macs3`)
- `genome_size` — effective genome size (`hs` for human, `mm` for mouse)
- `extsize` — fragment extension size in `--nomodel` mode (default `30`)

### Motif Finding (HOMER)

- `homer_genome` — HOMER genome build (default `hg38`)
- `homer_preparsed_dir` — shared cache directory for HOMER preparsed files
- `homer_len` — comma-separated list of motif lengths to search (default `"5,6,7,8,9,10"`)

---

## Usage

### 1. Set up the SLURM profile

Place the cluster `config.yaml` (the SLURM profile) in a directory accessible to Snakemake — for example:

```bash
mkdir -p ~/.config/snakemake/slurm_profile
cp config.yaml ~/.config/snakemake/slurm_profile/config.yaml  # the SLURM profile version
```

### 2. Activate the environment

```bash
mamba activate snake   # or your Snakemake conda env
```

### 3. Dry run

```bash
snakemake -n --profile slurm_profile
```

### 4. Execute

```bash
snakemake --profile slurm_profile --cores all
```

The profile sets `jobs: 50`, `latency-wait: 60`, `keep-going: true`, and `rerun-incomplete: true` by default. Two rules (`all` and `process_featurecounts`) run locally.

---

## Output Directory Structure

```
results/
├── qc/
│   ├── fastqc/                       # Per-sample FastQC reports
│   ├── fingerprint/                  # plotFingerprint PNGs (IP vs control)
│   └── multiqc_report.html           # Aggregated QC report
├── trimmed/                          # Adapter-trimmed FASTQs + trimming reports
├── aligned/                          # STAR-aligned BAMs + indices + STAR logs
├── counts/
│   ├── raw_counts.txt                # featureCounts raw output
│   ├── raw_counts.txt.summary        # Assignment summary
│   └── processed_counts.txt          # Tidied counts matrix
├── de/
│   ├── {condition}_IP_vs_{ctrl}_all.tsv          # Full DESeq2 results
│   └── {condition}_IP_vs_{ctrl}_significant.tsv  # Filtered significant hits
├── merged/                           # Replicate-merged BAMs per condition × type
├── bigwigs/                          # Log2-ratio bigWig tracks
│   ├── {condition}_rep{N}_IP_vs_{ctrl}_bs1.bw    # Per-replicate
│   └── {condition}_merged_IP_vs_{ctrl}_bs1.bw    # Merged
├── stranded/                         # Strand-split BAMs + merged IgG strand BAMs
├── peaks/{condition}/                # MACS2 narrowPeak, summits, XLS
│   ├── {sample}_{strand}_peaks.narrowPeak         # Per-replicate
│   └── combined_{strand}_peaks.narrowPeak         # Pooled IP
├── motifs/{condition}/               # HOMER de novo motif results
│   ├── {sample}_{strand}/homerResults.html
│   └── combined_{strand}/homerResults.html
└── logs/                             # Per-rule log files
    └── slurm/                        # SLURM stdout/stderr
```

---

## Key Design Decisions

**Single-end, strand-specific.** The STAR alignment uses `--alignEndsType EndToEnd` with stringent filtering (`--outFilterMultimapNmax 1`, `--outFilterMatchNmin 24`), and featureCounts runs with `-s 1` (stranded counting). This reflects the single-end eCLIP library geometry.

**Strand-aware peak calling.** BAMs are split by strand (SAM flag 0x10) before MACS2 peak calling, with merged IgG replicates serving as the control for each strand independently. This avoids artifacts from antisense signal contamination.

**DESeq2 enrichment model.** The R script filters low-abundance genes by CPM, then runs a standard DESeq2 model comparing IP vs control (IgG or Input). Significance is defined by both adjusted p-value and a minimum log2 fold-change threshold, which defaults to ~1.5-fold enrichment.

**HOMER in RNA mode.** Motif discovery uses the `-rna` flag, appropriate for identifying RNA-binding motifs from eCLIP peak sequences.

---

## Troubleshooting

**SLURM exit code 127 (command not found).** If running Snakemake inside a Singularity container, the container may mask SLURM binaries. Install Snakemake directly into a conda environment (e.g., `mamba create -n snake snakemake=7.32.4`) instead.

**HOMER preparsed directory errors.** HOMER caches genome parsing results in `homer_preparsed_dir`. Ensure this directory is writable and shared across runs to avoid redundant computation. On a shared HPC filesystem, set this to a location within your project allocation.

**featureCounts column name mismatches.** The `process_fc.py` script strips BAM path prefixes from featureCounts column headers. If sample names in the tidied output don't match `metadata.tsv`, check that your FASTQ filenames follow the expected `{sample}_S*` or `{sample}.bam` naming convention.