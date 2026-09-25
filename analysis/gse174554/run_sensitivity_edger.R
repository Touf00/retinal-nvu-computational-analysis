suppressPackageStartupMessages(library(edgeR))

defs <- c("EC95","EC99_ge2","EC99_broad")
for (d in defs) {
  base <- file.path("sensitivity", d)
  counts <- read.csv(gzfile(file.path(base,"pseudobulk_counts.csv.gz")),
                     row.names=1, check.names=FALSE)
  meta <- read.csv(file.path(base,"sample_metadata.csv"), stringsAsFactors=FALSE)
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
  q <- glmQLFTest(fit, coef="progressionRecurrent")
  tt <- topTags(q, n=Inf, sort.by="PValue")$table
  tt$gene <- rownames(tt)
  tt <- tt[,c("gene","logFC","logCPM","F","PValue","FDR")]
  write.csv(tt,file.path(base,"edgeR_all.csv"),row.names=FALSE)
  write.csv(tt[tt$FDR<0.05,,drop=FALSE],file.path(base,"edgeR_FDR05.csv"),row.names=FALSE)
  audit <- data.frame(metric=c("pairs","samples","genes_input","genes_filterByExpr","genes_FDR05","design_rank","common_dispersion"),
                      value=c(length(unique(meta$pair)),nrow(meta),nrow(y0),nrow(y),
                              sum(tt$FDR<0.05),qr(design)$rank,y$common.dispersion))
  write.csv(audit,file.path(base,"edgeR_audit.csv"),row.names=FALSE)
  cat("\n",d,"\n")
  print(audit,row.names=FALSE)
  print(head(tt,10),row.names=FALSE)
}
