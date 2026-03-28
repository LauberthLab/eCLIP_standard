suppressMessages({ library(DESeq2); library(edgeR); library(optparse) })

option_list <- list(
  make_option("--counts",    type="character"),
  make_option("--metadata",  type="character"),
  make_option("--condition", type="character"),
  make_option("--control",   type="character", default="IgG"),
  make_option("--min_cpm",   type="double",    default=0.3),
  make_option("--lfc",       type="double",    default=0.5849),
  make_option("--padj",      type="double",    default=0.05),
  make_option("--out_all",   type="character"),
  make_option("--out_sig",   type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Load counts, drop Length column
counts_raw <- read.delim(opt$counts, row.names=1, check.names=FALSE)
counts_raw <- counts_raw[, !colnames(counts_raw) %in% "Length", drop=FALSE]

# Filter metadata to this condition + control type
meta     <- read.delim(opt$metadata)
meta_sub <- meta[meta$condition==opt$condition & meta$type %in% c("IP",opt$control),]
meta_sub <- meta_sub[meta_sub$sample %in% colnames(counts_raw),]
counts   <- counts_raw[, meta_sub$sample, drop=FALSE]

group <- factor(meta_sub$type, levels=c(opt$control, "IP"))

# CPM filter: >min_cpm in >=2 samples
keep <- rowSums(edgeR::cpm(counts, normalized.lib.sizes=TRUE) > opt$min_cpm) >= 2
counts_filt <- counts[keep,, drop=FALSE]

# DESeq2
dds <- DESeqDataSetFromMatrix(round(counts_filt), DataFrame(group=group), ~group)
colData(dds)$group <- relevel(colData(dds)$group, ref=opt$control)
dds    <- DESeq(dds, quiet=TRUE)
result <- na.omit(as.data.frame(results(dds)))
result <- result[order(result$padj),]

dir.create(dirname(opt$out_all), showWarnings=FALSE, recursive=TRUE)
write.table(result, opt$out_all, sep="\t", quote=FALSE, col.names=NA)

sig <- result[result$padj <= opt$padj & result$log2FoldChange >= opt$lfc,]
write.table(sig, opt$out_sig, sep="\t", quote=FALSE, col.names=NA)
cat(sprintf("Significant genes: %d\n", nrow(sig)))
