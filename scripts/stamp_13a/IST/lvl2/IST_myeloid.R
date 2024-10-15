# Unsupervised clustering - InSituType

# Libraries
library(SingleCellExperiment)
library(dplyr)
library(here)
library(glue)
library(qs)
library(InSituType)

# Data loading
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_13a"
res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)

sub <- "myeloid"
sce <- sce[,sce$lvl0 == sub]

##############################################################################
# IST function
##############################################################################
ist <- function(sce, nk, gs=TRUE, pbs=NULL, bkg=TRUE) {
  # dependencies
  library(InSituType)
  library(SingleCellExperiment)
  # load counts
  mtx <- counts(sce[gs, ])
  mtx <- as(t(mtx), "dgCMatrix")
  # cohorting based on IF data
  j <- names(cd <- colData(sce))
  i <- grep("^Mean", j, value=TRUE)
  i <- setdiff(i, "Mean.G")
  i <- c("Area", "AspectRatio", i)
  coh <- fastCohorting(as.matrix(cd[i]))
  # background estimation
  neg <- grep("^neg", altExpNames(sce), value=TRUE)
  neg <- sce$nCount_negprobes/nrow(altExp(sce, neg))
  # update reference profiles
  pbs <- if (!is.null(pbs)) {
    bkg <- if (bkg) {
      rna <- sce$nCount_RNA
      rna*mean(neg)/mean(rna) 
    }
    updateReferenceProfiles(
      reference_profiles=pbs, counts=mtx,
      neg=neg, bg=bkg)$updated_profiles
  }
  # clustering
  insitutype(mtx, 
             reference_profiles=pbs,
             update_reference_profiles=FALSE,
             neg=neg, cohort=coh, n_clusts=nk)
}

##############################################################################
# Feature selection with IST profiles
##############################################################################
# Load IST profiles
dir <- glue("{proj_dir}/data/{stamp}/Ist")
unsup <- qread(glue("{dir}/unsup.qs"))

# Average the profiles
norm <- unsup$profiles /  Matrix::colSums(unsup$profiles) # average
norm <- norm + 1e-6

#Calculate logFC
result <- norm 
for (i in 1:nrow(norm)) {
  for (j in 1:ncol(norm)) {
    result[i, j] <- log2(norm[i, j] / mean(norm[i, -j]))
  }
}

#Take top features for each cluster by logFC
feat <- apply(result, 2, function(column) {
  names(sort(column, decreasing = TRUE))[1:80]
})

hvg <- c(feat[,"f"],feat[,"g"],feat[,"j"],feat[,"b"])

# Run IST
sce <- sce[hvg,]
unsup <- ist(sce, 2:7)

# Save
dir <- glue("{proj_dir}/data/{stamp}/Ist")
dir.create(dir, showWarnings = F, recursive = T)
qsave(unsup, file = glue("{dir}/unsup_{sub}.qs"))
