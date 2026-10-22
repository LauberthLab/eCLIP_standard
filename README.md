# eCLIP-seq Analysis Pipeline

Snakemake v7.32.4 and Python==3.10


A Snakemake pipeline for end-to-end processing of single-end enhanced CLIP (eCLIP) sequencing data, from raw FASTQ files through UMI-aware deduplication, differential enrichment analysis, strand-specific peak calling, and motif discovery.

---

## Overview

The pipeline performs the following sequential steps:

| Step | Tool | Description |
|------|------|-------------|
| 0 | UMI-tools `extract` | Move the 5′ UMI into the read name (before trimming) |
| 1 | FastQC | Quality check on raw reads (pre-extraction) |
| 2 | Trim Galore | Illumina adapter trimming + quality filtering (Q ≥ 28), on UMI-extracted reads |
| 3 | STAR | Alignment to hg38 (EndToEnd, unique mappers only) |
| 4 | samtools | BAM indexing |
| 4b | UMI-tools `dedup` | UMI-aware PCR-duplicate removal (5′ coordinate + UMI) |
| 5 | featureCounts | Strand-specific gene-level read counting on deduplicated BAMs (all samples jointly) |
| 6 | process counts | Tidy the raw counts matrix (clean column names, drop annotation columns) |
| 7 | DESeq2 | Differential enrichment: IP vs IgG |
| 8 | samtools merge | Merge replicate deduplicated BAMs per condition × sample type |
| 9 | bamCompare | Log2-ratio bigWig generation (individual replicates + merged) |
| 10 | plotFingerprint | Enrichment QC (IP vs control cumulative read distribution) |
| 11 | MultiQC | Aggregated QC report across all steps (incl. UMI-tools dedup stats) |
| 12 | samtools | Split deduplicated BAMs into plus- and minus-strand |
| 13 | samtools merge | Merge IgG replicates per strand per condition (peak-calling control) |
| 14 | MACS3 | Peak calling — individual IP replicates vs merged IgG, per strand |
| 15 | MACS3 | Peak calling — pooled IP replicates vs merged IgG, per strand |
| 16 | BEDTools | Reproducible peak set: pooled-IP peaks overlapping a peak in every replicate |
| 17 | HOMER | De novo motif discovery on the final reproducible peak set (`-rna` mode) |

**Why UMI dedup matters here.** UMIs are extracted at the 5′ end before trimming and used after alignment to collapse reads sharing both a 5′ mapping coordinate and a UMI — i.e. genuine PCR duplicates. Because duplicates are removed this way, `--keep-dup all` in the MACS3 rules is **correct**: reads still sharing a coordinate after dedup carry distinct UMIs and represent independent crosslink events that should be retained.

---

## Repository Structure

```
.
├── Snakefile                 # Main workflow definition
├── config.yaml               # Pipeline parameters (paths, thresholds, tool options)
├── envs/
│   └── eclip.yaml            # Conda environment used by --use-conda rules
├── metadata.tsv              # Sample sheet (sample ↔ condition ↔ type ↔ replicate ↔ batch ↔ FASTQ)
├── slurm_profile/            # SLURM profile (account, partition, default resources)
│   └── config.yaml
└── scripts/
    └── de_analysis.R         # DESeq2 differential enrichment analysis
```

---

## Requirements

Most tools are provisioned automatically per-rule via `--use-conda` from `envs/eclip.yaml`. Two rules deliberately use Quest **modules** instead of conda, to avoid heavyweight conda solves:

- **DESeq2** (`deseq2` rule) uses `module load R/4.2.0` with `DESeq2`, `edgeR`, `optparse`.
- **HOMER** (`homer_motif` rule) uses `module load homer`, because the module ships the pre-installed hg38 genome package (a fresh conda HOMER does not).

### Conda environment (`envs/eclip.yaml`)

Provides: UMI-tools (1.1.4), FastQC, Trim Galore (+ Cutadapt), STAR (2.7.5a), samtools (1.6), Subread (featureCounts), deepTools (bamCompare, plotFingerprint), MACS3 (3.0.2), BEDTools, MultiQC, pandas.

Snakemake itself should be installed separately (recommended: `mamba create -n snake snakemake`). STAR requires a pre-built hg38 index; if using the module HOMER, the hg38 genome package must already be installed system-wide (it is, on Quest).

### Reference Files

| File | Config Key | Description |
|------|-----------|-------------|
| STAR genome index | `genome_dir` | Pre-built hg38 STAR index directory |
| Gencode GTF | `gtf` | Gene annotation (e.g., gencode.v45.annotation.gtf) |

---

## Metadata Format

Tab-separated `metadata.tsv` with the following required columns (note: `fastq` is the **sixth** column — a `batch` column precedes it):

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier (used as wildcard throughout) |
| `condition` | Experimental condition grouping (e.g., `TOP1`) |
| `type` | One of `IP`, `IgG`, or `Input` |
| `replicate` | Integer replicate number |
| `batch` | Prep/sequencing batch label |
| `fastq` | Absolute path to the raw `.fastq.gz` file |

Example:

```
sample              condition  type   replicate  batch     fastq
21161FL-30-01_S1    TOP1       IP     1          batch_30  /path/to/sample_R1.fastq.gz
21161FL-30-05_S5    TOP1       IgG    1          batch_30  /path/to/igg_R1.fastq.gz
21161FL-30-03_S3    TOP1       Input  1          batch_30  /path/to/input_R1.fastq.gz
```

**Input samples are excluded from all processing.** Any row with `type == Input` is filtered out at load time and is never extracted, aligned, counted, or used as a control. Peak calling uses IgG as the control and runs only for conditions where both IP and IgG samples exist.

---

## Configuration

All tunable parameters live in `config.yaml`:

### Paths

- `metadata` — path to the sample sheet
- `genome_dir` — STAR genome index directory
- `gtf` — Gencode annotation GTF
- `results_dir` — top-level output directory (created automatically)

### UMI Deduplication (UMI-tools)

- `umi_pattern` — barcode pattern passed to `extract --bc-pattern` (default `"NNNNNNNNNN"`, a 10-nt 5′ UMI)
- `umi_method` — dedup collapse method (`unique` used)

### Differential Enrichment (DESeq2)

- `de_control` — baseline type: `"IgG"` (default)
- `min_cpm` — CPM filter threshold (default `0.3`)
- `lfc_threshold` — minimum log2 fold-change for significance (default `0.5849`, ~1.5-fold)
- `padj_threshold` — adjusted p-value cutoff (default `0.05`)

### BigWig Generation (bamCompare)

- `bin_size` — genomic bin size in bp (default `1`, single-base resolution)
- `normalize_using` — `RPKM`, `CPM`, `BPM`, `RPGC`, or `None`
- `pseudocount` — pseudocount added before log2 ratio (default `1`)

### Peak Calling (MACS3)

- `macs_cmd` — executable name (`macs3`)
- `genome_size` — effective genome size (`hs` for human)
- `extsize` — read extension size in `--nomodel` mode (default `30`)

### Motif Finding (HOMER)

- `homer_genome` — HOMER genome build (default `hg38`; requires the module HOMER with hg38 installed)
- `homer_preparsed_dir` — shared cache directory for HOMER preparsed files
- `homer_len` — comma-separated motif lengths (default `"5,6,7,8,9,10"`)

---

## Usage

The SLURM account and partition live in `slurm_profile/config.yaml`, not in the Snakefile.

### 1. Activate the Snakemake environment

```bash
mamba activate snake   # environment containing Snakemake
```

### 2. Pre-build the conda env once (so it is not on the clock during the run)

```bash
snakemake --use-conda --conda-create-envs-only --cores 1
```

### 3. Dry run

```bash
snakemake --profile slurm_profile --use-conda -n
```

### 4. Execute

```bash
snakemake --profile slurm_profile --use-conda --cores <N>
```

Set `--cores` at least as high as the largest rule's `threads` (STAR uses 50), or per-rule thread requests will be clamped. The profile controls concurrent SLURM jobs (`jobs: 50`), `latency-wait`, `keep-going`, and `rerun-incomplete`. Rules `all` and `process_featurecounts` run locally. Because each rule submits its own SLURM job with its own resources, samples are processed in parallel across the cluster.

---

## Output Directory Structure

```
results/
├── umi/                              # UMI-extracted FASTQs
├── qc/
│   ├── fastqc/                       # Per-sample FastQC reports
│   ├── fingerprint/                  # plotFingerprint PNGs (IP vs control)
│   └── multiqc_report.html           # Aggregated QC report (incl. dedup stats)
├── trimmed/                          # Adapter-trimmed FASTQs + trimming reports
├── aligned/                          # STAR-aligned BAMs + indices + STAR logs
├── dedup/                            # UMI-deduplicated BAMs + indices
├── counts/
│   ├── raw_counts.txt                # featureCounts raw output
│   ├── raw_counts.txt.summary        # Assignment summary
│   └── processed_counts.txt          # Tidied counts matrix
├── de/
│   ├── {condition}_IP_vs_IgG_all.tsv          # Full DESeq2 results
│   └── {condition}_IP_vs_IgG_significant.tsv  # Filtered significant hits
├── merged/                           # Replicate-merged deduplicated BAMs
├── bigwigs/                          # Log2-ratio bigWig tracks (bs1)
│   ├── {condition}_rep{N}_IP_vs_IgG_bs1.bw    # Per-replicate
│   └── {condition}_merged_IP_vs_IgG_bs1.bw    # Merged
├── stranded/                         # Strand-split BAMs + merged IgG strand BAMs
├── peaks/{condition}/                # MACS3 narrowPeak, summits, XLS
│   ├── {sample}_{strand}_peaks.narrowPeak      # Per-replicate
│   ├── combined_{strand}_peaks.narrowPeak      # Pooled IP
│   └── final_peaks.narrowPeak                  # Reproducible set (all-replicate overlap)
├── motifs/{condition}/final/         # HOMER de novo motif results on final peaks
│   └── homerResults.html
└── logs/                             # Per-rule log files (incl. umi_extract, umi_dedup)
```

---

## Key Design Decisions

**UMI-aware deduplication.** A 10-nt UMI is extracted from the 5′ end before trimming and used after alignment (`dedup --method`) to collapse PCR duplicates by 5′ coordinate + UMI. This is what makes `--keep-dup all` correct at the MACS3 step. The `umi_pattern` should be confirmed against the library's actual read structure before a production run; `directional` is the recommended dedup method (error-correcting), whereas `unique` under-collapses.

**Single-end, strand-specific.** STAR uses `--alignEndsType EndToEnd` with stringent filtering (`--outFilterMultimapNmax 1`, `--outFilterMatchNmin 24`), and featureCounts runs with `-s 1`. This reflects single-end eCLIP library geometry. Note that unique-mappers-only alignment under-represents multi-copy loci (rDNA, some snoRNAs); handling of those loci should be described explicitly if they are central to the analysis.

**Strand-aware peak calling.** Deduplicated BAMs are split by strand (SAM flag 0x10) before MACS3, with merged IgG replicates as the per-strand control, avoiding antisense contamination artifacts.

**Reproducible peak set.** The final peak set is defined as pooled-IP peaks that overlap a called peak in every individual replicate (BEDTools `intersect`), after merging plus/minus strands. This is an overlap-based reproducibility filter rather than IDR.

**HOMER in RNA mode.** Motif discovery uses `-rna`, appropriate for RNA-binding motifs, and runs on the final reproducible peak set.

---

## Troubleshooting

**STAR out-of-memory.** The hg38 index needs ~32–40 GB; the rule requests more for safety. If the profile caps memory via `default-resources`, that overrides the rule — check `slurm_profile/config.yaml`.

**HOMER "Genome hg38 not found."** A conda-installed HOMER ships without genome packages. Use the module HOMER (which has hg38 installed), or point `homer_genome` at an hg38 FASTA whose contig names match the peak file (the pipeline prepends `chr`).

**HOMER preparsed directory.** Set `homer_preparsed_dir` to a writable, persistent location within your project allocation so preparsing is cached across runs.

**Conda activation aborts with an unbound-variable error.** Some conda env deactivate hooks reference unset variables that are fatal under `set -u`. `--use-conda` avoids this by letting Snakemake manage activation per job; do not manually switch conda envs inside rule shells.
