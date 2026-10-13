#!/usr/bin/env Rscript
# Export original author-provided spatial annotations from Zucha et al. PNAS 2024.
# DOI 10.17632/gnb2dsjms2.1; preserves source Seurat spot IDs and actual
# metadata-defined sample keys. It does NOT infer lesion region or merge day7
# sections. A separate, manually verified section-to-GSM map is REQUIRED for
# downstream spatial source alignment.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop("Usage: Rscript export_zucha_region_annotations.R <published_final.rds> <output.csv>")
rds <- args[[1]]; output <- args[[2]]
if (!file.exists(rds)) stop("Missing published Seurat RDS: ", rds)
if (!requireNamespace("SeuratObject",quietly=TRUE))
  stop("Install SeuratObject: install.packages('SeuratObject')")
obj <- readRDS(rds)
if (!inherits(obj, "Seurat")) stop("Expected author Seurat object, got ", paste(class(obj),collapse=","))
md <- obj@meta.data
if (is.null(md) || nrow(md)==0 || is.null(rownames(md))) stop("No author cell/spot metadata")
required <- c("Condition","DetailedRegionAnnoShort")
if (!all(required %in% colnames(md)))
  stop("Required published author fields absent: ", paste(setdiff(required,colnames(md)),collapse=", "))
out <- data.frame(
  original_seurat_spot_id=rownames(md),
  condition=as.character(md$Condition),
  DetailedRegionAnnoShort=as.character(md$DetailedRegionAnnoShort),
  stringsAsFactors=FALSE
)
if (anyDuplicated(out$original_seurat_spot_id))
  stop("Non-unique Seurat spot IDs. Never join with duplicates.")
# The original codebase uses Barcodes and Condition, but Condition may have
# only four time labels for five GEO sections. Do not use Condition as sample ID.
barcode_field <- c("Barcodes","barcode","Barcode","barcodes")
barcode_field <- barcode_field[barcode_field %in% colnames(md)]
out$author_barcode <- if (length(barcode_field))
  as.character(md[[barcode_field[[1]]]]) else rep("",nrow(md))
# Keep raw candidate sample/slide identifiers for manual verification.
sample_fields <- c("orig.ident","sample_id","sample","Sample","slide","Slide",
                   "image","Image","image_id","section","Section")
present_fields <- sample_fields[sample_fields %in% colnames(md)]
for (field in present_fields)
  out[[paste0("author_meta_", gsub("[^A-Za-z0-9]","_",field))]] <- as.character(md[[field]])
if (length(present_fields)) {
  out$author_sample <- as.character(md[[present_fields[[1]]]])
  out$author_sample_field <- present_fields[[1]]
} else {
  out$author_sample <- rep("",nrow(md))
  out$author_sample_field <- rep("",nrow(md))
  warning("No explicit sample identifier. A 5-section barcode mapping CANNOT yet be assumed.")
}
if ("DetailedIntAnno" %in% names(md)) out$DetailedIntAnno <- as.character(md$DetailedIntAnno)
if ("IntegratedAnno" %in% names(md)) out$IntegratedAnno <- as.character(md$IntegratedAnno)
if (anyNA(out$condition) || any(!nzchar(out$condition)))
  stop("Missing Condition metadata")
dir.create(dirname(output), recursive=TRUE,showWarnings=FALSE)
write.csv(out,file=output,row.names=FALSE,na="",quote=TRUE)
message("Exported ",nrow(out)," original author spot rows.")
message("Sample-field candidates: ",paste(present_fields,collapse=", "))
message("Unique chosen author sample IDs: ",paste(unique(out$author_sample),collapse=", "))
message("Condition counts: "); print(table(out$condition,useNA="ifany"))
message("Author-region counts: "); print(table(out$DetailedRegionAnnoShort,useNA="ifany"))
message("Next: verify a distinct author_sample identifier for each of the FIVE")
message("GSE233814 sections and produce a manually audited mapping CSV.")
message("No region or cell-type biological inference was performed here.")
