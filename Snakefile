"""
eCLIP-seq Analysis Pipeline
============================
Steps:
  1.  FastQC          – quality check on raw reads
  2.  Trim Galore      – adapter trimming + quality filter (--illumina, -q 28)
  3.  STAR             – alignment to hg38 (EndToEnd, unique mappers)
  4.  samtools index   – index BAMs
  5.  featureCounts    – gene-level counting across all samples
  6.  process_fc.py    – tidy the counts matrix
  7.  DESeq2           – differential enrichment: IP vs IgG (default) or Input
  8.  samtools merge   – merge replicate BAMs per condition × type
  9.  bamCompare       – log2 bigwigs, individual replicates + merged
  10. plotFingerprint  – enrichment QC plots
  11. MultiQC          – aggregate QC report
  12. split_strand     – split BAMs into plus/minus strand
  13. merge_igg_stranded – merge IgG replicates per strand per condition
  14. macs_individual  – MACS2 per IP replicate × strand vs merged IgG
  15. macs_combined    – MACS2 all IP reps pooled × strand vs merged IgG
  16. intersect_peaks  – merge strands + intersect reps → final peak set
  17. homer_motif      – findMotifsGenome.pl on final reproducible peaks

Metadata TSV columns: sample, condition, type (IP|IgG|Input), replicate, fastq

Usage:
    snakemake --profile slurm_profile --cores all
"""

import re
import pandas as pd
from pathlib import Path

configfile: "config.yaml"

# ─────────────────────────────────────────────────────────────────────────────
# Environment: load HPC modules for every shell job
# Adjust module names to match your cluster (run `module avail <tool>`)
# ─────────────────────────────────────────────────────────────────────────────
shell.prefix(
    "set -euo pipefail; "
    "module purge; "
    "module load fastqc; "
    "module load cutadapt/4.2; "
    "module load TrimGalore; "
    "module load STAR; "
    "module load samtools; "
    "module load subread; "
    "module load deeptools; "
    "module load MACS3; "
    "module load homer; "
    "module load multiqc; "
    "module load R/4.2.0; "
    "module load bedtools; "
)

RESULTS = config["results_dir"]
GENOME  = config["genome_dir"]
GTF     = config["gtf"]
DE_CTRL = config.get("de_control", "IgG")
BIN_SIZE    = config.get("bin_size", 1)
NORMALIZE   = config.get("normalize_using", "RPKM")
PSEUDOCOUNT = config.get("pseudocount", 1)

# Peak calling / motif config
MACS_CMD        = config.get("macs_cmd",           "macs3")
GENOME_SIZE     = config.get("genome_size",         "hs")
EXTSIZE         = config.get("extsize",             30)
HOMER_GENOME    = config.get("homer_genome",        "hg38")
HOMER_PREPARSED = config.get("homer_preparsed_dir", f"{config.get('results_dir','results')}/homer_preparsed")
HOMER_LEN       = config.get("homer_len",           "5,6,7,8,9,10")

meta = pd.read_csv(config["metadata"], sep="\t").set_index("sample", drop=False)
meta["replicate"] = meta["replicate"].astype(int)

# ─────────────────────────────────────────────────────────────────────────────
# Helper lookups
# ─────────────────────────────────────────────────────────────────────────────
def samples_of(type_=None, condition=None):
    df = meta
    if type_:      df = df[df["type"]      == type_]
    if condition:  df = df[df["condition"] == condition]
    return df["sample"].tolist()

def get_raw_fastq(wildcards):
    return meta.loc[wildcards.sample, "fastq"]

def get_all_bams(_):
    return [f"{RESULTS}/aligned/{s}Aligned.sortedByCoord.out.bam" for s in all_samples]

def get_all_bam_indices(_):
    return [f"{RESULTS}/aligned/{s}Aligned.sortedByCoord.out.bam.bai" for s in all_samples]

def get_bams_for_merge(wildcards):
    return [f"{RESULTS}/aligned/{s}Aligned.sortedByCoord.out.bam"
            for s in samples_of(type_=wildcards.stype, condition=wildcards.condition)]

def get_bam_indices_for_merge(wildcards):
    return [f"{RESULTS}/aligned/{s}Aligned.sortedByCoord.out.bam.bai"
            for s in samples_of(type_=wildcards.stype, condition=wildcards.condition)]

def _bam_for(condition, stype, rep):
    rows = meta[(meta["condition"] == condition) &
                (meta["type"]      == stype) &
                (meta["replicate"] == int(rep))]
    if rows.empty:
        raise ValueError(f"No sample found: condition={condition} type={stype} rep={rep}")
    return f"{RESULTS}/aligned/{rows.iloc[0]['sample']}Aligned.sortedByCoord.out.bam"

def get_ip_bam(wildcards):   return _bam_for(wildcards.condition, "IP",            wildcards.rep)
def get_ctrl_bam(wildcards): return _bam_for(wildcards.condition, wildcards.ctrl,  wildcards.rep)
def get_ip_bai(wildcards):   return get_ip_bam(wildcards)   + ".bai"
def get_ctrl_bai(wildcards): return get_ctrl_bam(wildcards) + ".bai"

# ─────────────────────────────────────────────────────────────────────────────
# Derived target lists
# ─────────────────────────────────────────────────────────────────────────────
conditions   = meta["condition"].unique().tolist()
all_samples  = meta["sample"].tolist()
ip_samples   = samples_of(type_="IP")
avail_ctrls  = [t for t in ["IgG", "Input"] if t in meta["type"].values]

wildcard_constraints:
    sample = "|".join(re.escape(s) for s in all_samples)

de_pairs = [
    (c, DE_CTRL) for c in conditions
    if samples_of(type_="IP", condition=c) and samples_of(type_=DE_CTRL, condition=c)
]

indiv_bw_targets, fingerprint_targets = [], []
for cond in conditions:
    ip_reps = set(meta[(meta["condition"] == cond) & (meta["type"] == "IP")]["replicate"])
    for ctrl in avail_ctrls:
        ctrl_reps = set(meta[(meta["condition"] == cond) & (meta["type"] == ctrl)]["replicate"])
        for rep in ip_reps & ctrl_reps:
            indiv_bw_targets.append(f"{RESULTS}/bigwigs/{cond}_rep{rep}_IP_vs_{ctrl}_bs{BIN_SIZE}.bw")
            fingerprint_targets.append(f"{RESULTS}/qc/fingerprint/{cond}_rep{rep}_IP_vs_{ctrl}.png")

merged_bw_targets = [
    f"{RESULTS}/bigwigs/{cond}_merged_IP_vs_{ctrl}_bs{BIN_SIZE}.bw"
    for cond in conditions for ctrl in avail_ctrls
    if samples_of(type_="IP", condition=cond) and samples_of(type_=ctrl, condition=cond)
]

de_all_targets = [f"{RESULTS}/de/{c}_IP_vs_{ctrl}_all.tsv"         for c, ctrl in de_pairs]
de_sig_targets = [f"{RESULTS}/de/{c}_IP_vs_{ctrl}_significant.tsv" for c, ctrl in de_pairs]

# Peak calling targets (strand-specific, IgG control only)
STRANDS = ["plus", "minus"]
indiv_peak_targets, combined_peak_targets = [], []
final_peak_targets, final_homer_targets = [], []
for cond in conditions:
    ip_samps  = samples_of(type_="IP",  condition=cond)
    igg_samps = samples_of(type_="IgG", condition=cond)
    if not (ip_samps and igg_samps):
        continue
    for strand in STRANDS:
        for samp in ip_samps:
            indiv_peak_targets.append(f"{RESULTS}/peaks/{cond}/{samp}_{strand}_peaks.narrowPeak")
        combined_peak_targets.append(f"{RESULTS}/peaks/{cond}/combined_{strand}_peaks.narrowPeak")
    final_peak_targets.append(f"{RESULTS}/peaks/{cond}/final_peaks.narrowPeak")
    final_homer_targets.append(f"{RESULTS}/motifs/{cond}/final/homerResults.html")

# ─────────────────────────────────────────────────────────────────────────────
# rule all
# ─────────────────────────────────────────────────────────────────────────────
localrules: all, process_featurecounts

rule all:
    input:
        expand(f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.html", sample=all_samples),
        expand(f"{RESULTS}/trimmed/{{sample}}_trimming_report.txt", sample=all_samples),
        expand(f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam.bai", sample=all_samples),
        f"{RESULTS}/counts/raw_counts.txt",
        f"{RESULTS}/counts/processed_counts.txt",
        de_all_targets, de_sig_targets,
        indiv_bw_targets, merged_bw_targets,
        fingerprint_targets,
        f"{RESULTS}/qc/multiqc_report.html",
        indiv_peak_targets,
        combined_peak_targets,
        final_peak_targets,
        final_homer_targets,

# ─────────────────────────────────────────────────────────────────────────────
# 1. FastQC – raw reads
# ─────────────────────────────────────────────────────────────────────────────
rule fastqc_raw:
    input:  get_raw_fastq
    output:
        html = f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.html",
        zip  = f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.zip",
    params: outdir = f"{RESULTS}/qc/fastqc"
    log:    f"{RESULTS}/logs/fastqc/{{sample}}.log"
    resources: ntasks=4, mem="8gb", time="0:30:00"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {resources.ntasks} -o {params.outdir} {input} 2> {log}
        raw_base=$(basename {input} .fastq.gz); raw_base=$(basename $raw_base .fq.gz)
        if [ "$raw_base" != "{wildcards.sample}" ]; then
            mv {params.outdir}/${{raw_base}}_fastqc.html {output.html}
            mv {params.outdir}/${{raw_base}}_fastqc.zip  {output.zip}
        fi
        """

# ─────────────────────────────────────────────────────────────────────────────
# 2. Trim Galore
# ─────────────────────────────────────────────────────────────────────────────
rule trim_galore:
    input:  get_raw_fastq
    output:
        trimmed = f"{RESULTS}/trimmed/{{sample}}_trimmed.fq.gz",
        report  = f"{RESULTS}/trimmed/{{sample}}_trimming_report.txt",
    params: outdir = f"{RESULTS}/trimmed"
    log:    f"{RESULTS}/logs/trim_galore/{{sample}}.log"
    resources: ntasks=4, mem="30gb", time="1:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        ln -sf $(realpath {input}) {params.outdir}/{wildcards.sample}.fastq.gz
        trim_galore --illumina -q 28 --fastqc -j {resources.ntasks} \
            -o {params.outdir} {params.outdir}/{wildcards.sample}.fastq.gz 2> {log}
        mv {params.outdir}/{wildcards.sample}.fastq.gz_trimming_report.txt {output.report} \
            2>/dev/null || mv {params.outdir}/{wildcards.sample}_trimming_report.txt {output.report}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 3. STAR alignment
# ─────────────────────────────────────────────────────────────────────────────
rule star_align:
    input:  f"{RESULTS}/trimmed/{{sample}}_trimmed.fq.gz"
    output:
        bam = f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam",
        log = f"{RESULTS}/aligned/{{sample}}Log.final.out",
    params: genome=GENOME, prefix=f"{RESULTS}/aligned/{{sample}}"
    log:    f"{RESULTS}/logs/star/{{sample}}.log"
    resources: ntasks=20, mem="100gb", time="2:00:00"
    shell:
        """
        mkdir -p {RESULTS}/aligned
        STAR --runMode alignReads --runThreadN {resources.ntasks} \
            --readFilesCommand zcat --genomeDir {params.genome} \
            --alignEndsType EndToEnd --genomeLoad NoSharedMemory \
            --alignMatesGapMax 15000 --readFilesIn {input} \
            --outFilterMultimapNmax 1 --outFileNamePrefix {params.prefix} \
            --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout --outReadsUnmapped Fastx \
            --outFilterScoreMin 10 --outFilterMatchNmin 24 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 4. samtools index
# ─────────────────────────────────────────────────────────────────────────────
rule samtools_index:
    input:  f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam"
    output: f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam.bai"
    log:    f"{RESULTS}/logs/samtools_index/{{sample}}.log"
    resources: ntasks=8, mem="8gb", time="0:30:00"
    shell:  "samtools index -@ {resources.ntasks} {input} 2> {log}"

# ─────────────────────────────────────────────────────────────────────────────
# 5. featureCounts – all samples at once
# ─────────────────────────────────────────────────────────────────────────────
rule featurecounts:
    input:  bams=get_all_bams, indices=get_all_bam_indices
    output:
        counts  = f"{RESULTS}/counts/raw_counts.txt",
        summary = f"{RESULTS}/counts/raw_counts.txt.summary",
    params: gtf=GTF, outdir=f"{RESULTS}/counts"
    log:    f"{RESULTS}/logs/featurecounts/featurecounts.log"
    resources: ntasks=6, mem="16gb", time="1:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        featureCounts -T {resources.ntasks} -t exon -g gene_name -s 1 \
            -a {params.gtf} -o {output.counts} {input.bams} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 6. Process featureCounts output
# ─────────────────────────────────────────────────────────────────────────────
rule process_featurecounts:
    input:  f"{RESULTS}/counts/raw_counts.txt"
    output: f"{RESULTS}/counts/processed_counts.txt"
    log:    f"{RESULTS}/logs/process_fc/process_fc.log"
    run:
        import os
        fc = pd.read_csv(str(input), sep="\t", header=1)
        fc = fc.drop(columns=["Chr","Start","End","Strand"], errors="ignore")
        fc.columns = [
            c if c in ("Geneid","Length")
            else (os.path.basename(c).replace("Aligned.sortedByCoord.out.bam","")
                  if "Aligned.sortedByCoord.out.bam" in os.path.basename(c)
                  else os.path.basename(c).replace(".bam","").split("_S")[0]
                       if "_S" in os.path.basename(c)
                       else os.path.basename(c).replace(".bam",""))
            for c in fc.columns
        ]
        os.makedirs(os.path.dirname(os.path.abspath(str(output))), exist_ok=True)
        fc.to_csv(str(output), sep="\t", index=False)
        with open(str(log), "w") as lf:
            lf.write(f"Processed {len(fc)} genes x {len(fc.columns)-2} samples -> {output}\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7. DESeq2 differential enrichment
# ─────────────────────────────────────────────────────────────────────────────
rule deseq2:
    input:  counts=f"{RESULTS}/counts/processed_counts.txt", metadata=config["metadata"]
    output:
        all_res = f"{RESULTS}/de/{{condition}}_IP_vs_{{control}}_all.tsv",
        sig_res = f"{RESULTS}/de/{{condition}}_IP_vs_{{control}}_significant.tsv",
    params:
        min_cpm=config.get("min_cpm",0.3), lfc_thresh=config.get("lfc_threshold",0.5849),
        padj_thresh=config.get("padj_threshold",0.05)
    log:    f"{RESULTS}/logs/deseq2/{{condition}}_vs_{{control}}.log"
    resources: ntasks=2, mem="16gb", time="1:00:00"
    shell:
        """
        mkdir -p {RESULTS}/de
        Rscript scripts/de_analysis.R \
            --counts {input.counts} --metadata {input.metadata} \
            --condition {wildcards.condition} --control {wildcards.control} \
            --min_cpm {params.min_cpm} --lfc {params.lfc_thresh} --padj {params.padj_thresh} \
            --out_all {output.all_res} --out_sig {output.sig_res} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 8. Merge replicate BAMs (per condition × type)
# ─────────────────────────────────────────────────────────────────────────────
rule merge_bams:
    input:  bams=get_bams_for_merge, indices=get_bam_indices_for_merge
    output: bam=f"{RESULTS}/merged/{{condition}}_{{stype}}_merged.bam"
    log:    f"{RESULTS}/logs/merge/{{condition}}_{{stype}}.log"
    resources: ntasks=16, mem="16gb", time="1:00:00"
    shell:
        """
        mkdir -p {RESULTS}/merged
        if [ $(echo {input.bams} | wc -w) -eq 1 ]; then
            ln -sf $(realpath {input.bams}) {output.bam}
        else
            samtools merge -f -@ {resources.ntasks} {output.bam} {input.bams} 2> {log}
        fi
        """

rule index_merged_bam:
    input:  f"{RESULTS}/merged/{{condition}}_{{stype}}_merged.bam"
    output: f"{RESULTS}/merged/{{condition}}_{{stype}}_merged.bam.bai"
    log:    f"{RESULTS}/logs/merge/{{condition}}_{{stype}}_index.log"
    resources: ntasks=8, mem="8gb", time="0:30:00"
    shell:  "samtools index -@ {resources.ntasks} {input} 2> {log}"

# ─────────────────────────────────────────────────────────────────────────────
# 9a. bamCompare – individual replicates
# ─────────────────────────────────────────────────────────────────────────────
rule bamcompare_individual:
    input:  ip_bam=get_ip_bam, ctrl_bam=get_ctrl_bam, ip_bai=get_ip_bai, ctrl_bai=get_ctrl_bai
    output: f"{RESULTS}/bigwigs/{{condition}}_rep{{rep}}_IP_vs_{{ctrl}}_bs{BIN_SIZE}.bw"
    log:    f"{RESULTS}/logs/bigwigs/{{condition}}_rep{{rep}}_vs_{{ctrl}}.log"
    resources: ntasks=16, mem="16gb", time="2:00:00"
    shell:
        """
        mkdir -p {RESULTS}/bigwigs
        bamCompare -b1 {input.ip_bam} -b2 {input.ctrl_bam} -o {output} \
            -p {resources.ntasks} --normalizeUsing {NORMALIZE} \
            --scaleFactorsMethod None --pseudocount {PSEUDOCOUNT} --binSize {BIN_SIZE} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 9b. bamCompare – merged replicates
# ─────────────────────────────────────────────────────────────────────────────
rule bamcompare_merged:
    input:
        ip_bam   = f"{RESULTS}/merged/{{condition}}_IP_merged.bam",
        ctrl_bam = f"{RESULTS}/merged/{{condition}}_{{ctrl}}_merged.bam",
        ip_bai   = f"{RESULTS}/merged/{{condition}}_IP_merged.bam.bai",
        ctrl_bai = f"{RESULTS}/merged/{{condition}}_{{ctrl}}_merged.bam.bai",
    output: f"{RESULTS}/bigwigs/{{condition}}_merged_IP_vs_{{ctrl}}_bs{BIN_SIZE}.bw"
    log:    f"{RESULTS}/logs/bigwigs/{{condition}}_merged_vs_{{ctrl}}.log"
    resources: ntasks=16, mem="16gb", time="2:00:00"
    shell:
        """
        mkdir -p {RESULTS}/bigwigs
        bamCompare -b1 {input.ip_bam} -b2 {input.ctrl_bam} -o {output} \
            -p {resources.ntasks} --normalizeUsing {NORMALIZE} \
            --scaleFactorsMethod None --pseudocount {PSEUDOCOUNT} --binSize {BIN_SIZE} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 10. plotFingerprint
# ─────────────────────────────────────────────────────────────────────────────
rule fingerprint:
    input:  ip_bam=get_ip_bam, ctrl_bam=get_ctrl_bam, ip_bai=get_ip_bai, ctrl_bai=get_ctrl_bai
    output: f"{RESULTS}/qc/fingerprint/{{condition}}_rep{{rep}}_IP_vs_{{ctrl}}.png"
    log:    f"{RESULTS}/logs/fingerprint/{{condition}}_rep{{rep}}_vs_{{ctrl}}.log"
    resources: ntasks=8, mem="8gb", time="1:00:00"
    shell:
        """
        mkdir -p {RESULTS}/qc/fingerprint
        plotFingerprint -b {input.ip_bam} {input.ctrl_bam} \
            --labels "IP_{wildcards.condition}_rep{wildcards.rep}" \
                     "{wildcards.ctrl}_{wildcards.condition}_rep{wildcards.rep}" \
            --plotTitle "{wildcards.condition} rep{wildcards.rep} IP vs {wildcards.ctrl}" \
            -p {resources.ntasks} -plot {output} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 11. MultiQC
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc:
    input:
        expand(f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.zip", sample=all_samples),
        expand(f"{RESULTS}/aligned/{{sample}}Log.final.out",  sample=all_samples),
        f"{RESULTS}/counts/raw_counts.txt.summary",
    output: f"{RESULTS}/qc/multiqc_report.html"
    params: outdir=f"{RESULTS}/qc", scandir=RESULTS
    log:    f"{RESULTS}/logs/multiqc/multiqc.log"
    resources: ntasks=2, mem="8gb", time="0:30:00"
    shell:  "multiqc {params.scandir} -o {params.outdir} --force 2> {log}"

# ─────────────────────────────────────────────────────────────────────────────
# 12. Split BAM by strand
#     plus  → -F 16  (exclude reverse-strand flag)
#     minus → -f 16  (require reverse-strand flag)
# ─────────────────────────────────────────────────────────────────────────────
rule split_strand:
    input:
        bam = f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam",
        bai = f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam.bai",
    output:
        plus  = f"{RESULTS}/stranded/{{sample}}_plus.bam",
        minus = f"{RESULTS}/stranded/{{sample}}_minus.bam",
    log:    f"{RESULTS}/logs/stranded/{{sample}}.log"
    resources: ntasks=8, mem="8gb", time="0:30:00"
    shell:
        """
        mkdir -p {RESULTS}/stranded
        samtools view -@ {resources.ntasks} -F 16 -b {input.bam} -o {output.plus}  2>  {log}
        samtools view -@ {resources.ntasks} -f 16 -b {input.bam} -o {output.minus} 2>> {log}
        samtools index -@ {resources.ntasks} {output.plus}  2>> {log}
        samtools index -@ {resources.ntasks} {output.minus} 2>> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 13. Merge IgG strand BAMs across replicates (per condition)
# ─────────────────────────────────────────────────────────────────────────────
def get_igg_stranded_bams(wildcards):
    return [f"{RESULTS}/stranded/{s}_{wildcards.strand}.bam"
            for s in samples_of(type_="IgG", condition=wildcards.condition)]

rule merge_igg_stranded:
    input:  bams=get_igg_stranded_bams
    output:
        bam = f"{RESULTS}/stranded/{{condition}}_IgG_merged_{{strand}}.bam",
        bai = f"{RESULTS}/stranded/{{condition}}_IgG_merged_{{strand}}.bam.bai",
    log:    f"{RESULTS}/logs/stranded/{{condition}}_IgG_merged_{{strand}}.log"
    resources: ntasks=16, mem="16gb", time="0:30:00"
    shell:
        """
        mkdir -p {RESULTS}/stranded
        if [ $(echo {input.bams} | wc -w) -eq 1 ]; then
            ln -sf $(realpath {input.bams}) {output.bam}
        else
            samtools merge -f -@ {resources.ntasks} {output.bam} {input.bams} 2> {log}
        fi
        samtools index -@ {resources.ntasks} {output.bam} 2>> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 14. MACS – individual IP replicate vs merged IgG (per strand)
# ─────────────────────────────────────────────────────────────────────────────
def _ip_condition(wildcards):
    return meta.loc[wildcards.sample, "condition"]

rule macs_individual:
    input:
        ip_bam  = f"{RESULTS}/stranded/{{sample}}_{{strand}}.bam",
        igg_bam = lambda wc: f"{RESULTS}/stranded/{_ip_condition(wc)}_IgG_merged_{wc.strand}.bam",
    output:
        narrowpeak = f"{RESULTS}/peaks/{{condition}}/{{sample}}_{{strand}}_peaks.narrowPeak",
        xls        = f"{RESULTS}/peaks/{{condition}}/{{sample}}_{{strand}}_peaks.xls",
        summits    = f"{RESULTS}/peaks/{{condition}}/{{sample}}_{{strand}}_summits.bed",
    params: outdir=f"{RESULTS}/peaks/{{condition}}", name="{sample}_{strand}",
            gsize=GENOME_SIZE, ext=EXTSIZE, macs=MACS_CMD
    log:    f"{RESULTS}/logs/peaks/{{condition}}/{{sample}}_{{strand}}.log"
    resources: ntasks=4, mem="16gb", time="1:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        {params.macs} callpeak \
            -t {input.ip_bam} -c {input.igg_bam} \
            -f BAM -n {params.name} -g {params.gsize} \
            -B --keep-dup all \
            --outdir {params.outdir} \
            --nomodel --extsize {params.ext} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 15. MACS – all IP replicates combined vs merged IgG (per strand)
# ─────────────────────────────────────────────────────────────────────────────
def get_all_ip_stranded_bams(wildcards):
    return [f"{RESULTS}/stranded/{s}_{wildcards.strand}.bam"
            for s in samples_of(type_="IP", condition=wildcards.condition)]

rule macs_combined:
    input:
        ip_bams = get_all_ip_stranded_bams,
        igg_bam = f"{RESULTS}/stranded/{{condition}}_IgG_merged_{{strand}}.bam",
    output:
        narrowpeak = f"{RESULTS}/peaks/{{condition}}/combined_{{strand}}_peaks.narrowPeak",
        xls        = f"{RESULTS}/peaks/{{condition}}/combined_{{strand}}_peaks.xls",
        summits    = f"{RESULTS}/peaks/{{condition}}/combined_{{strand}}_summits.bed",
    params: outdir=f"{RESULTS}/peaks/{{condition}}", name="combined_{strand}",
            gsize=GENOME_SIZE, ext=EXTSIZE, macs=MACS_CMD
    log:    f"{RESULTS}/logs/peaks/{{condition}}/combined_{{strand}}.log"
    resources: ntasks=4, mem="16gb", time="1:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        {params.macs} callpeak \
            -t {input.ip_bams} -c {input.igg_bam} \
            -f BAM -n {params.name} -g {params.gsize} \
            -B --keep-dup all \
            --outdir {params.outdir} \
            --nomodel --extsize {params.ext} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 16. Merge strands + intersect reps → final reproducible peak set
# ─────────────────────────────────────────────────────────────────────────────
def _ip_samples_for_condition(wildcards):
    return samples_of(type_="IP", condition=wildcards.condition)

def _get_indiv_strand_peaks(wildcards):
    """All per-rep, per-strand narrowPeaks for this condition."""
    samps = _ip_samples_for_condition(wildcards)
    return [f"{RESULTS}/peaks/{wildcards.condition}/{s}_{strand}_peaks.narrowPeak"
            for s in samps for strand in STRANDS]

rule intersect_peaks:
    input:
        indiv      = _get_indiv_strand_peaks,
        comb_plus  = f"{RESULTS}/peaks/{{condition}}/combined_plus_peaks.narrowPeak",
        comb_minus = f"{RESULTS}/peaks/{{condition}}/combined_minus_peaks.narrowPeak",
    output:
        final = f"{RESULTS}/peaks/{{condition}}/final_peaks.narrowPeak",
    params:
        outdir     = f"{RESULTS}/peaks/{{condition}}",
        ip_samples = _ip_samples_for_condition,
    log: f"{RESULTS}/logs/peaks/{{condition}}/intersect.log"
    resources: ntasks=2, mem="8gb", time="0:30:00"
    shell:
        """
        set -euo pipefail
        # 1. Merge strands per replicate
        rep_beds=()
        for samp in {params.ip_samples}; do
            cat {params.outdir}/${{samp}}_plus_peaks.narrowPeak \
                {params.outdir}/${{samp}}_minus_peaks.narrowPeak \
              | sort -k1,1 -k2,2n > {params.outdir}/${{samp}}_all.narrowPeak
            rep_beds+=( {params.outdir}/${{samp}}_all.narrowPeak )
        done

        # 2. Merge strands for combined
        cat {input.comb_plus} {input.comb_minus} \
          | sort -k1,1 -k2,2n > {params.outdir}/combined_all.narrowPeak

        # 3. Intersect: combined peaks supported by every replicate
        cp {params.outdir}/combined_all.narrowPeak {params.outdir}/_intersect_tmp.narrowPeak
        for bed in "${{rep_beds[@]}}"; do
            bedtools intersect \
                -a {params.outdir}/_intersect_tmp.narrowPeak \
                -b "$bed" -u \
              > {params.outdir}/_intersect_next.narrowPeak
            mv {params.outdir}/_intersect_next.narrowPeak \
               {params.outdir}/_intersect_tmp.narrowPeak
        done
        mv {params.outdir}/_intersect_tmp.narrowPeak {output.final}

        # Clean up intermediates
        rm -f {params.outdir}/*_all.narrowPeak

        wc -l {output.final} > {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# 17. HOMER motif finding on final reproducible peak set
# ─────────────────────────────────────────────────────────────────────────────
rule homer_motif:
    input:  peaks=f"{RESULTS}/peaks/{{condition}}/final_peaks.narrowPeak"
    output: html=f"{RESULTS}/motifs/{{condition}}/final/homerResults.html"
    params:
        outdir=f"{RESULTS}/motifs/{{condition}}/final",
        genome=HOMER_GENOME, preparsed=HOMER_PREPARSED, lengths=HOMER_LEN
    log:    f"{RESULTS}/logs/motifs/{{condition}}/final.log"
    resources: ntasks=16, mem="32gb", time="4:00:00", homer_lock=1
    shell:
        """
        mkdir -p {params.outdir} {params.preparsed}
        awk 'BEGIN{{OFS="\\t"}} {{if($1 !~ /^chr/) $1="chr"$1; print}}' \
            {input.peaks} > {params.outdir}/peaks_chrFixed.narrowPeak
        findMotifsGenome.pl \
            {params.outdir}/peaks_chrFixed.narrowPeak {params.genome} {params.outdir} \
            -p {resources.ntasks} -rna -S 10 \
            -len {params.lengths} \
            -preparsedDir {params.preparsed} 2> {log}
        """
