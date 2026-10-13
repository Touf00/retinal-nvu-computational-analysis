#!/usr/bin/env Rscript
# Extract *original author-assigned spatial region labels* from published Seurat RDS.
# Source: Zucha et al. PNAS 2024, Mendeley DOI 10.17632/gnb2dsjms2.1.
# Example: Rscript IMACS_Stroke_PET/scripts/export_zucha_region_annotations.R \
#   /path/to/spatial_seurat_1DP_13_nygen.rds \
#   ./zucha_author_region_metadata.csv
#
# Need R + SeuratObject to load RDS class. This exports metadata ONLY (small CSV);
# there is no count matrix redistribution or inferred lesion label.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: Rscript export_zucha_region_annotations.R <published_final.rds> <output.csv>")
rds <- args[1]; output <- args[2]
if(!file.exists(rds)) stop("Input RDS missing: ",rds)
if(!requireNamespace("SeuratObject",quietly=TRUE))
  stop("Install the R SeuratObject package: install.packages('SeuratObject')")
obj <- readRDS(rds)
md <- obj@meta.data
if(is.null(md) || nrow(md)==0) stop("RDS has no metadata")
required <- c("Condition","DetailedRegionAnnoShort")
missing <- setdiff(required,colnames(md))
if(length(missing)>0) stop("Author annotations absent: ",paste(missing,collapse=", "),
                           ". Do not infer core/periphery without original metadata.")
out <- data.frame(
  original_seurat_spot_id=rownames(md),
  condition=as.character(md$Condition),
  DetailedRegionAnnoShort=as.character(md$DetailedRegionAnnoShort),
  stringsAsFactors=FALSE
)
if("DetailedIntAnno" %in% colnames(md))
  out$DetailedIntAnno <- as.character(md$DetailedIntAnno)
if("IntegratedAnno" %in% colnames(md))
  out$IntegratedAnno <- as.character(md$IntegratedAnno)
if("barcode" %in% colnames(md))
  out$author_barcode <- as.character(md$barcode)
if(anyDuplicated(out$original_seurat_spot_id)) stop("Nonunique author Seurat spot IDs")
if(anyNA(out$DetailedRegionAnnoShort)) warning("Some author regions are missing; keep them NA")
if(anyNA(out$condition)) stop("Missing section condition; investigate before join")
write.csv(out,file=output,row.names=FALSE,na="",quote=TRUE)
message("Exported ",nrow(out)," rows to ",normalizePath(output))
message("Conditions: "); print(table(out$condition,useNA="ifany"))
message("Regions: "); print(table(out$DetailedRegionAnnoShort,useNA="ifany"))
message("Next: inspect section-specific barcode formatting before mapping to GEO spots.")
