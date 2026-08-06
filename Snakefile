"""
eCLIP-seq pipeline (Snakemake, --use-conda)                        UMI-aware

Run:
    snakemake --profile slurm_profile --use-conda --cores all

- Every rule uses envs/eclip.yaml EXCEPT deseq2, which uses the R/4.2.0 module
  (keeping a heavy R/Bioconductor stack out of the conda solve).
- account/partition come from your slurm_profile (not set here).
- Input-type samples are ignored everywhere (see all_samples / used types).
- UMI dedup upstream makes MACS `--keep-dup all` correct.

metadata.tsv columns: sample  condition  type(IP|IgG|Input)  replicate  batch  fastq
"""

import re
import pandas as pd

configfile: "config.yaml"

RESULTS = config["results_dir"]
GENOME  = config["genome_dir"]
GTF     = config["gtf"]
DE_CTRL = config.get("de_control", "IgG")

UMI_PATTERN = config.get("umi_pattern", "NNNNNNNNNN")
UMI_METHOD  = config.get("umi_method",  "unique")

BIN_SIZE    = config.get("bin_size", 1)
NORMALIZE   = config.get("normalize_using", "RPKM")
PSEUDOCOUNT = config.get("pseudocount", 1)

MACS_CMD     = config.get("macs_cmd", "macs3")
GENOME_SIZE  = config.get("genome_size", "hs")
EXTSIZE      = config.get("extsize", 30)
HOMER_GENOME = config.get("homer_genome", "hg38")
HOMER_PREPARSED = config.get("homer_preparsed_dir", f"{RESULTS}/homer_preparsed")
HOMER_LEN    = config.get("homer_len", "5,6,7,8,9,10")

CONDA_ENV = "envs/eclip.yaml"

# ── metadata (fastq is column 6; Input excluded from processing) ──────────────
meta = pd.read_csv(config["metadata"], sep="\t").set_index("sample", drop=False)
meta["replicate"] = meta["replicate"].astype(int)
meta_use = meta[meta["type"] != "Input"]          #  Input ignored everywhere

def samples_of(type_=None, condition=None):
    df = meta_use
    if type_:     df = df[df["type"] == type_]
    if condition: df = df[df["condition"] == condition]
    return df["sample"].tolist()

def get_raw_fastq(wildcards):
    return meta.loc[wildcards.sample, "fastq"]

def get_all_bams(_):
    return [f"{RESULTS}/dedup/{s}.dedup.bam" for s in all_samples]
def get_all_bam_indices(_):
    return [f"{RESULTS}/dedup/{s}.dedup.bam.bai" for s in all_samples]
def get_bams_for_merge(wildcards):
    return [f"{RESULTS}/dedup/{s}.dedup.bam"
            for s in samples_of(type_=wildcards.stype, condition=wildcards.condition)]
def get_bam_indices_for_merge(wildcards):
    return [f"{RESULTS}/dedup/{s}.dedup.bam.bai"
            for s in samples_of(type_=wildcards.stype, condition=wildcards.condition)]

def _bam_for(condition, stype, rep):
    rows = meta_use[(meta_use["condition"] == condition) &
                    (meta_use["type"] == stype) &
                    (meta_use["replicate"] == int(rep))]
    if rows.empty:
        raise ValueError(f"No sample: {condition}/{stype}/rep{rep}")
    return f"{RESULTS}/dedup/{rows.iloc[0]['sample']}.dedup.bam"

def get_ip_bam(w):   return _bam_for(w.condition, "IP",     w.rep)
def get_ctrl_bam(w): return _bam_for(w.condition, w.ctrl,   w.rep)
def get_ip_bai(w):   return get_ip_bam(w)   + ".bai"
def get_ctrl_bai(w): return get_ctrl_bam(w) + ".bai"

conditions  = meta_use["condition"].unique().tolist()
all_samples = meta_use["sample"].tolist()             # excludes Input
avail_ctrls = [t for t in ["IgG"] if t in meta_use["type"].values]   # Input never a control

wildcard_constraints:
    sample    = "|".join(re.escape(s) for s in all_samples),
    condition = "|".join(re.escape(c) for c in conditions),
    stype     = "IP|IgG",
    ctrl      = "IgG",
    rep       = r"\d+"

de_pairs = [(c, DE_CTRL) for c in conditions
            if samples_of("IP", c) and samples_of(DE_CTRL, c)]

indiv_bw_targets, fingerprint_targets = [], []
for cond in conditions:
    ip_reps = set(meta_use[(meta_use.condition==cond)&(meta_use.type=="IP")]["replicate"])
    for ctrl in avail_ctrls:
        ctrl_reps = set(meta_use[(meta_use.condition==cond)&(meta_use.type==ctrl)]["replicate"])
        for rep in ip_reps & ctrl_reps:
            indiv_bw_targets.append(f"{RESULTS}/bigwigs/{cond}_rep{rep}_IP_vs_{ctrl}_bs{BIN_SIZE}.bw")
            fingerprint_targets.append(f"{RESULTS}/qc/fingerprint/{cond}_rep{rep}_IP_vs_{ctrl}.png")

merged_bw_targets = [f"{RESULTS}/bigwigs/{cond}_merged_IP_vs_{ctrl}_bs{BIN_SIZE}.bw"
    for cond in conditions for ctrl in avail_ctrls
    if samples_of("IP", cond) and samples_of(ctrl, cond)]

de_all_targets = [f"{RESULTS}/de/{c}_IP_vs_{ctrl}_all.tsv"         for c, ctrl in de_pairs]
de_sig_targets = [f"{RESULTS}/de/{c}_IP_vs_{ctrl}_significant.tsv" for c, ctrl in de_pairs]

STRANDS = ["plus", "minus"]
indiv_peak_targets, combined_peak_targets = [], []
final_peak_targets, final_homer_targets = [], []
for cond in conditions:
    if not (samples_of("IP", cond) and samples_of("IgG", cond)):
        continue
    for strand in STRANDS:
        for samp in samples_of("IP", cond):
            indiv_peak_targets.append(f"{RESULTS}/peaks/{cond}/{samp}_{strand}_peaks.narrowPeak")
        combined_peak_targets.append(f"{RESULTS}/peaks/{cond}/combined_{strand}_peaks.narrowPeak")
    final_peak_targets.append(f"{RESULTS}/peaks/{cond}/final_peaks.narrowPeak")
    final_homer_targets.append(f"{RESULTS}/motifs/{cond}/final/homerResults.html")

localrules: all, process_featurecounts

rule all:
    input:
        expand(f"{RESULTS}/trimmed/{{sample}}_trimming_report.txt", sample=all_samples),
        expand(f"{RESULTS}/dedup/{{sample}}.dedup.bam.bai", sample=all_samples),
        f"{RESULTS}/counts/processed_counts.txt",
        de_all_targets, de_sig_targets,
        indiv_bw_targets, merged_bw_targets, fingerprint_targets,
        f"{RESULTS}/qc/multiqc_report.html",
        indiv_peak_targets, combined_peak_targets,
        final_peak_targets, final_homer_targets,

# ── 0. UMI extract ────────────────────────────────────────────────────────────
rule umi_extract:
    input:  get_raw_fastq
    output: f"{RESULTS}/umi/{{sample}}_umi.fq.gz"
    log:    f"{RESULTS}/logs/umi_extract/{{sample}}.log"
    conda:  CONDA_ENV
    resources: mem_mb=8000, runtime=60
    shell:
        "umi_tools extract --extract-method=string --bc-pattern={UMI_PATTERN} "
        "--stdin={input} --stdout={output} --log={log}"

# ── 1. FastQC (raw) ───────────────────────────────────────────────────────────
rule fastqc_raw:
    input:  get_raw_fastq
    output: html=f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.html",
            zip =f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.zip"
    params: outdir=f"{RESULTS}/qc/fastqc"
    log:    f"{RESULTS}/logs/fastqc/{{sample}}.log"
    conda:  CONDA_ENV
    threads: 4
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input} 2> {log}
        b=$(basename {input} .fastq.gz); b=$(basename $b .fq.gz)
        if [ "$b" != "{wildcards.sample}" ]; then
            mv {params.outdir}/${{b}}_fastqc.html {output.html}
            mv {params.outdir}/${{b}}_fastqc.zip  {output.zip}
        fi
        """

# ── 2. Trim Galore (on UMI-extracted reads) ───────────────────────────────────
rule trim_galore:
    input:  f"{RESULTS}/umi/{{sample}}_umi.fq.gz"
    output: trimmed=f"{RESULTS}/trimmed/{{sample}}_trimmed.fq.gz",
            report =f"{RESULTS}/trimmed/{{sample}}_trimming_report.txt"
    params: outdir=f"{RESULTS}/trimmed"
    log:    f"{RESULTS}/logs/trim_galore/{{sample}}.log"
    conda:  CONDA_ENV
    threads: 4
    resources: mem_mb=30000, runtime=60
    shell:
        r"""
        mkdir -p {params.outdir}
        ln -sf $(realpath {input}) {params.outdir}/{wildcards.sample}.fastq.gz
        trim_galore --illumina -q 28 --fastqc -j {threads} \
            -o {params.outdir} {params.outdir}/{wildcards.sample}.fastq.gz 2> {log}
        mv {params.outdir}/{wildcards.sample}.fastq.gz_trimming_report.txt {output.report} \
            2>/dev/null || mv {params.outdir}/{wildcards.sample}_trimming_report.txt {output.report}
        """

# ── 3. STAR ───────────────────────────────────────────────────────────────────
rule star_align:
    input:  f"{RESULTS}/trimmed/{{sample}}_trimmed.fq.gz"
    output: bam=f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam",
            log=f"{RESULTS}/aligned/{{sample}}Log.final.out"
    params: genome=GENOME, prefix=f"{RESULTS}/aligned/{{sample}}"
    log:    f"{RESULTS}/logs/star/{{sample}}.log"
    conda:  CONDA_ENV
    threads: 50
    resources: ntasks=20, mem="100gb", time="2:00:00"
    shell:
        r"""
        mkdir -p {RESULTS}/aligned
        STAR --runMode alignReads --runThreadN {threads} \
            --readFilesCommand zcat --genomeDir {params.genome} \
            --alignEndsType EndToEnd --genomeLoad NoSharedMemory \
            --alignMatesGapMax 15000 --readFilesIn {input} \
            --outFilterMultimapNmax 1 --outFileNamePrefix {params.prefix} \
            --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout --outReadsUnmapped Fastx \
            --outFilterScoreMin 10 --outFilterMatchNmin 24 2> {log}
        """

# ── 4. index raw BAM ──────────────────────────────────────────────────────────
rule samtools_index:
    input:  f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam"
    output: f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam.bai"
    conda:  CONDA_ENV
    threads: 8
    resources: mem_mb=8000, runtime=30
    shell: "samtools index -@ {threads} {input}"

# ── 4b. UMI dedup ─────────────────────────────────────────────────────────────
rule umi_dedup:
    input:  bam=f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam",
            bai=f"{RESULTS}/aligned/{{sample}}Aligned.sortedByCoord.out.bam.bai"
    output: bam=f"{RESULTS}/dedup/{{sample}}.dedup.bam",
            bai=f"{RESULTS}/dedup/{{sample}}.dedup.bam.bai"
    log:    f"{RESULTS}/logs/umi_dedup/{{sample}}.log"
    conda:  CONDA_ENV
    threads: 4
    resources: mem_mb=32000, runtime=120
    shell:
        "umi_tools dedup --stdin={input.bam} --stdout={output.bam} "
        "--method={UMI_METHOD} --log={log} && samtools index -@ {threads} {output.bam}"

# ── 5. featureCounts ──────────────────────────────────────────────────────────
rule featurecounts:
    input:  bams=get_all_bams, indices=get_all_bam_indices
    output: counts=f"{RESULTS}/counts/raw_counts.txt",
            summary=f"{RESULTS}/counts/raw_counts.txt.summary"
    params: gtf=GTF
    log:    f"{RESULTS}/logs/featurecounts/featurecounts.log"
    conda:  CONDA_ENV
    threads: 6
    resources: mem_mb=16000, runtime=60
    shell:
        "featureCounts -T {threads} -t exon -g gene_name -s 1 "
        "-a {params.gtf} -o {output.counts} {input.bams} 2> {log}"

# ── 6. process counts (localrule, no conda needed but harmless) ────────────────
rule process_featurecounts:
    input:  f"{RESULTS}/counts/raw_counts.txt"
    output: f"{RESULTS}/counts/processed_counts.txt"
    run:
        import os
        fc = pd.read_csv(str(input), sep="\t", header=1)
        fc = fc.drop(columns=["Chr","Start","End","Strand"], errors="ignore")
        def clean(c):
            if c in ("Geneid","Length"): return c
            b = os.path.basename(c)
            return b.replace(".dedup.bam","").replace("Aligned.sortedByCoord.out.bam","").replace(".bam","")
        fc.columns = [clean(c) for c in fc.columns]
        os.makedirs(os.path.dirname(os.path.abspath(str(output))), exist_ok=True)
        fc.to_csv(str(output), sep="\t", index=False)

# ── 7. DESeq2 (module R, NOT conda) ───────────────────────────────────────────
rule deseq2:
    input:  counts=f"{RESULTS}/counts/processed_counts.txt", metadata=config["metadata"]
    output: all_res=f"{RESULTS}/de/{{condition}}_IP_vs_{{control}}_all.tsv",
            sig_res=f"{RESULTS}/de/{{condition}}_IP_vs_{{control}}_significant.tsv"
    params: min_cpm=config.get("min_cpm",0.3), lfc=config.get("lfc_threshold",0.5849),
            padj=config.get("padj_threshold",0.05)
    log:    f"{RESULTS}/logs/deseq2/{{condition}}_vs_{{control}}.log"
    resources: mem_mb=16000, runtime=60
    shell:
        r"""
        module load R/4.2.0
        Rscript scripts/de_analysis.R \
            --counts {input.counts} --metadata {input.metadata} \
            --condition {wildcards.condition} --control {wildcards.control} \
            --min_cpm {params.min_cpm} --lfc {params.lfc} --padj {params.padj} \
            --out_all {output.all_res} --out_sig {output.sig_res} 2> {log}
        """

# ── 8. merge BAMs ─────────────────────────────────────────────────────────────
rule merge_bams:
    input:  bams=get_bams_for_merge, indices=get_bam_indices_for_merge
    output: f"{RESULTS}/merged/{{condition}}_{{stype}}_merged.bam"
    conda:  CONDA_ENV
    threads: 16
    resources: mem_mb=16000, runtime=60
    shell:
        r"""
        mkdir -p {RESULTS}/merged
        if [ $(echo {input.bams} | wc -w) -eq 1 ]; then ln -sf $(realpath {input.bams}) {output}
        else samtools merge -f -@ {threads} {output} {input.bams}; fi
        """

rule index_merged_bam:
    input:  f"{RESULTS}/merged/{{condition}}_{{stype}}_merged.bam"
    output: f"{RESULTS}/merged/{{condition}}_{{stype}}_merged.bam.bai"
    conda:  CONDA_ENV
    threads: 8
    resources: mem_mb=8000, runtime=30
    shell: "samtools index -@ {threads} {input}"

# ── 9. bamCompare ─────────────────────────────────────────────────────────────
rule bamcompare_individual:
    input:  ip_bam=get_ip_bam, ctrl_bam=get_ctrl_bam, ip_bai=get_ip_bai, ctrl_bai=get_ctrl_bai
    output: f"{RESULTS}/bigwigs/{{condition}}_rep{{rep}}_IP_vs_{{ctrl}}_bs{BIN_SIZE}.bw"
    log:    f"{RESULTS}/logs/bigwigs/{{condition}}_rep{{rep}}_vs_{{ctrl}}.log"
    conda:  CONDA_ENV
    threads: 16
    resources: mem_mb=16000, runtime=120
    shell:
        "bamCompare -b1 {input.ip_bam} -b2 {input.ctrl_bam} -o {output} -p {threads} "
        "--normalizeUsing {NORMALIZE} --scaleFactorsMethod None "
        "--pseudocount {PSEUDOCOUNT} --binSize {BIN_SIZE} 2> {log}"

rule bamcompare_merged:
    input:  ip_bam=f"{RESULTS}/merged/{{condition}}_IP_merged.bam",
            ctrl_bam=f"{RESULTS}/merged/{{condition}}_{{ctrl}}_merged.bam",
            ip_bai=f"{RESULTS}/merged/{{condition}}_IP_merged.bam.bai",
            ctrl_bai=f"{RESULTS}/merged/{{condition}}_{{ctrl}}_merged.bam.bai"
    output: f"{RESULTS}/bigwigs/{{condition}}_merged_IP_vs_{{ctrl}}_bs{BIN_SIZE}.bw"
    log:    f"{RESULTS}/logs/bigwigs/{{condition}}_merged_vs_{{ctrl}}.log"
    conda:  CONDA_ENV
    threads: 16
    resources: mem_mb=16000, runtime=120
    shell:
        "bamCompare -b1 {input.ip_bam} -b2 {input.ctrl_bam} -o {output} -p {threads} "
        "--normalizeUsing {NORMALIZE} --scaleFactorsMethod None "
        "--pseudocount {PSEUDOCOUNT} --binSize {BIN_SIZE} 2> {log}"

# ── 10. fingerprint ───────────────────────────────────────────────────────────
rule fingerprint:
    input:  ip_bam=get_ip_bam, ctrl_bam=get_ctrl_bam, ip_bai=get_ip_bai, ctrl_bai=get_ctrl_bai
    output: f"{RESULTS}/qc/fingerprint/{{condition}}_rep{{rep}}_IP_vs_{{ctrl}}.png"
    log:    f"{RESULTS}/logs/fingerprint/{{condition}}_rep{{rep}}_vs_{{ctrl}}.log"
    conda:  CONDA_ENV
    threads: 8
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        mkdir -p {RESULTS}/qc/fingerprint
        plotFingerprint -b {input.ip_bam} {input.ctrl_bam} \
            --labels "IP_{wildcards.condition}_rep{wildcards.rep}" "{wildcards.ctrl}_{wildcards.condition}_rep{wildcards.rep}" \
            --plotTitle "{wildcards.condition} rep{wildcards.rep} IP vs {wildcards.ctrl}" \
            -p {threads} -plot {output} 2> {log}
        """

# ── 11. MultiQC ───────────────────────────────────────────────────────────────
rule multiqc:
    input:  expand(f"{RESULTS}/qc/fastqc/{{sample}}_fastqc.zip", sample=all_samples),
            expand(f"{RESULTS}/aligned/{{sample}}Log.final.out", sample=all_samples),
            expand(f"{RESULTS}/logs/umi_dedup/{{sample}}.log", sample=all_samples),
            f"{RESULTS}/counts/raw_counts.txt.summary"
    output: f"{RESULTS}/qc/multiqc_report.html"
    params: outdir=f"{RESULTS}/qc", scandir=RESULTS
    log:    f"{RESULTS}/logs/multiqc/multiqc.log"
    conda:  CONDA_ENV
    resources: mem_mb=8000, runtime=30
    shell: "multiqc {params.scandir} -o {params.outdir} -n multiqc_report.html --force 2> {log}"

# ── 12. split strand ──────────────────────────────────────────────────────────
rule split_strand:
    input:  bam=f"{RESULTS}/dedup/{{sample}}.dedup.bam",
            bai=f"{RESULTS}/dedup/{{sample}}.dedup.bam.bai"
    output: plus=f"{RESULTS}/stranded/{{sample}}_plus.bam",
            minus=f"{RESULTS}/stranded/{{sample}}_minus.bam"
    conda:  CONDA_ENV
    threads: 8
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        mkdir -p {RESULTS}/stranded
        samtools view -@ {threads} -F 16 -b {input.bam} -o {output.plus}
        samtools view -@ {threads} -f 16 -b {input.bam} -o {output.minus}
        samtools index -@ {threads} {output.plus}
        samtools index -@ {threads} {output.minus}
        """

# ── 13. merge IgG stranded ────────────────────────────────────────────────────
def get_igg_stranded_bams(w):
    return [f"{RESULTS}/stranded/{s}_{w.strand}.bam" for s in samples_of("IgG", w.condition)]

rule merge_igg_stranded:
    input:  bams=get_igg_stranded_bams
    output: bam=f"{RESULTS}/stranded/{{condition}}_IgG_merged_{{strand}}.bam",
            bai=f"{RESULTS}/stranded/{{condition}}_IgG_merged_{{strand}}.bam.bai"
    conda:  CONDA_ENV
    threads: 16
    resources: mem_mb=16000, runtime=30
    shell:
        r"""
        mkdir -p {RESULTS}/stranded
        if [ $(echo {input.bams} | wc -w) -eq 1 ]; then ln -sf $(realpath {input.bams}) {output.bam}
        else samtools merge -f -@ {threads} {output.bam} {input.bams}; fi
        samtools index -@ {threads} {output.bam}
        """

# ── 14. MACS per rep ──────────────────────────────────────────────────────────
def _ip_condition(w): return meta_use.loc[w.sample, "condition"]

rule macs_individual:
    input:  ip_bam=f"{RESULTS}/stranded/{{sample}}_{{strand}}.bam",
            igg_bam=lambda w: f"{RESULTS}/stranded/{_ip_condition(w)}_IgG_merged_{w.strand}.bam"
    output: narrowpeak=f"{RESULTS}/peaks/{{condition}}/{{sample}}_{{strand}}_peaks.narrowPeak",
            xls=f"{RESULTS}/peaks/{{condition}}/{{sample}}_{{strand}}_peaks.xls",
            summits=f"{RESULTS}/peaks/{{condition}}/{{sample}}_{{strand}}_summits.bed"
    params: outdir=f"{RESULTS}/peaks/{{condition}}", name="{sample}_{strand}"
    log:    f"{RESULTS}/logs/peaks/{{condition}}/{{sample}}_{{strand}}.log"
    conda:  CONDA_ENV
    threads: 4
    resources: mem_mb=16000, runtime=60
    shell:
        "mkdir -p {params.outdir} && {MACS_CMD} callpeak -t {input.ip_bam} -c {input.igg_bam} "
        "-f BAM -n {params.name} -g {GENOME_SIZE} -B --keep-dup all "
        "--outdir {params.outdir} --nomodel --extsize {EXTSIZE} 2> {log}"

# ── 15. MACS combined ─────────────────────────────────────────────────────────
def get_all_ip_stranded_bams(w):
    return [f"{RESULTS}/stranded/{s}_{w.strand}.bam" for s in samples_of("IP", w.condition)]

rule macs_combined:
    input:  ip_bams=get_all_ip_stranded_bams,
            igg_bam=f"{RESULTS}/stranded/{{condition}}_IgG_merged_{{strand}}.bam"
    output: narrowpeak=f"{RESULTS}/peaks/{{condition}}/combined_{{strand}}_peaks.narrowPeak",
            xls=f"{RESULTS}/peaks/{{condition}}/combined_{{strand}}_peaks.xls",
            summits=f"{RESULTS}/peaks/{{condition}}/combined_{{strand}}_summits.bed"
    params: outdir=f"{RESULTS}/peaks/{{condition}}", name="combined_{strand}"
    log:    f"{RESULTS}/logs/peaks/{{condition}}/combined_{{strand}}.log"
    conda:  CONDA_ENV
    threads: 4
    resources: mem_mb=16000, runtime=60
    shell:
        "mkdir -p {params.outdir} && {MACS_CMD} callpeak -t {input.ip_bams} -c {input.igg_bam} "
        "-f BAM -n {params.name} -g {GENOME_SIZE} -B --keep-dup all "
        "--outdir {params.outdir} --nomodel --extsize {EXTSIZE} 2> {log}"

# ── 16. intersect ─────────────────────────────────────────────────────────────
def _get_indiv_strand_peaks(w):
    return [f"{RESULTS}/peaks/{w.condition}/{s}_{strand}_peaks.narrowPeak"
            for s in samples_of("IP", w.condition) for strand in STRANDS]

rule intersect_peaks:
    input:  indiv=_get_indiv_strand_peaks,
            comb_plus=f"{RESULTS}/peaks/{{condition}}/combined_plus_peaks.narrowPeak",
            comb_minus=f"{RESULTS}/peaks/{{condition}}/combined_minus_peaks.narrowPeak"
    output: final=f"{RESULTS}/peaks/{{condition}}/final_peaks.narrowPeak"
    params: outdir=f"{RESULTS}/peaks/{{condition}}",
            ip_samples=lambda w: samples_of("IP", w.condition)
    log:    f"{RESULTS}/logs/peaks/{{condition}}/intersect.log"
    conda:  CONDA_ENV
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        set -euo pipefail
        rep_beds=()
        for samp in {params.ip_samples}; do
            cat {params.outdir}/${{samp}}_plus_peaks.narrowPeak {params.outdir}/${{samp}}_minus_peaks.narrowPeak \
              | sort -k1,1 -k2,2n > {params.outdir}/${{samp}}_all.narrowPeak
            rep_beds+=( {params.outdir}/${{samp}}_all.narrowPeak )
        done
        cat {input.comb_plus} {input.comb_minus} | sort -k1,1 -k2,2n > {params.outdir}/combined_all.narrowPeak
        cp {params.outdir}/combined_all.narrowPeak {params.outdir}/_tmp.narrowPeak
        for bed in "${{rep_beds[@]}}"; do
            bedtools intersect -a {params.outdir}/_tmp.narrowPeak -b "$bed" -u > {params.outdir}/_next.narrowPeak
            mv {params.outdir}/_next.narrowPeak {params.outdir}/_tmp.narrowPeak
        done
        mv {params.outdir}/_tmp.narrowPeak {output.final}
        rm -f {params.outdir}/*_all.narrowPeak
        wc -l {output.final} > {log}
        """

# ── 17. HOMER ─────────────────────────────────────────────────────────────────
rule homer_motif:
    input:  f"{RESULTS}/peaks/{{condition}}/final_peaks.narrowPeak"
    output: f"{RESULTS}/motifs/{{condition}}/final/homerResults.html"
    params: outdir=f"{RESULTS}/motifs/{{condition}}/final", preparsed=HOMER_PREPARSED
    log:    f"{RESULTS}/logs/motifs/{{condition}}/final.log"
    threads: 16
    resources: mem_mb=32000, runtime=240
    shell:
        r"""
        module load homer
        mkdir -p {params.outdir} {params.preparsed}
        awk 'BEGIN{{OFS="\t"}} {{if($1 !~ /^chr/) $1="chr"$1; print}}' {input} > {params.outdir}/peaks_chrFixed.narrowPeak
        findMotifsGenome.pl {params.outdir}/peaks_chrFixed.narrowPeak {HOMER_GENOME} {params.outdir} \
            -p {threads} -rna -S 10 -len {HOMER_LEN} -preparsedDir {params.preparsed} 2> {log}
        """
