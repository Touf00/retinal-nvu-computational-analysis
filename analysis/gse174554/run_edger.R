suppressPackageStartupMessages(library(edgeR))

counts <- read.csv(gzfile("combined/pseudobulk_counts.csv.gz"),
                   row.names=1, check.names=FALSE)
meta <- read.csv("combined/sample_metadata.csv", stringsAsFactors=FALSE)

stopifnot(nrow(meta)==28, length(unique(meta$pair))==14)
stopifnot(all(meta$sample %in% colnames(counts)))
counts <- counts[, meta$sample, drop=FALSE]

meta$pair <- factor(meta$pair)
meta$progression <- factor(meta$progression, levels=c("Primary","Recurrent"))
design <- model.matrix(~ pair + progression, data=meta)
stopifnot(qr(design)$rank == ncol(design))

y0 <- DGEList(counts=counts)
keep <- filterByExpr(y0, design=design)
y <- y0[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
coef_name <- "progressionRecurrent"
stopifnot(coef_name %in% colnames(design))
qlf <- glmQLFTest(fit, coef=coef_name)

tt <- topTags(qlf, n=Inf, sort.by="PValue")$table
tt$gene <- rownames(tt)
tt <- tt[, c("gene","logFC","logCPM","F","PValue","FDR")]
write.csv(tt, "combined/edgeR_paired_recurrent_vs_primary_all.csv", row.names=FALSE)

sig <- tt[tt$FDR < 0.05, , drop=FALSE]
write.csv(sig, "combined/edgeR_paired_recurrent_vs_primary_FDR05.csv", row.names=FALSE)

cpm_mat <- cpm(y, log=TRUE, prior.count=2)
write.csv(data.frame(gene=rownames(cpm_mat), cpm_mat, check.names=FALSE),
          "combined/logCPM_filtered.csv", row.names=FALSE)

mds <- plotMDS(y, plot=FALSE)
mds_df <- data.frame(sample=colnames(y), pair=meta$pair,
                     progression=meta$progression,
                     MDS1=mds$x, MDS2=mds$y)
write.csv(mds_df, "combined/MDS_coordinates.csv", row.names=FALSE)

qc <- data.frame(
  sample=colnames(y),
  pair=meta$pair,
  progression=meta$progression,
  raw_library_size=y0$samples$lib.size,
  filtered_library_size=y$samples$lib.size,
  norm_factor=y$samples$norm.factors
)
write.csv(qc, "combined/edgeR_sample_QC.csv", row.names=FALSE)

audit <- data.frame(
  metric=c("genes_input","genes_filterByExpr","genes_FDR05",
           "design_rank","design_columns","common_dispersion"),
  value=c(nrow(y0), nrow(y), nrow(sig), qr(design)$rank,
          ncol(design), y$common.dispersion)
)
write.csv(audit, "combined/edgeR_audit.csv", row.names=FALSE)

cat("GENES INPUT:", nrow(y0), "\n")
cat("GENES RETAINED:", nrow(y), "\n")
cat("FDR<0.05:", nrow(sig), "\n")
cat("COMMON DISPERSION:", y$common.dispersion, "\n")
cat("TOP 20\n")
print(head(tt,20), row.names=FALSE)
