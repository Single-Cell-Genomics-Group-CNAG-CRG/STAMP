# Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(scater)
  library(glue)
  library(qs)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
  library(data.table)
  library(spatstat)
  library(InSituType)
  library(EnhancedVolcano)
  library(DESeq2)
  
})
# Data loading
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_13a"
res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)
####################################################################################
# 1. Activated monocytes : ctrl vs lps at 4h
####################################################################################
sub <- sce[,sce$lvl1 == "act. mono." & sce$experiment != "aCD3aCD28" & sce$timepoint == "4h"] # sub

sub <- logNormCounts(sub) # lognorm
# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ replicate + experiment) # indicate replicates and exp columns

dds <- DESeq(dds, betaPrior=FALSE)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="experiment_LPS_vs_ctrl")
# Shrink log fold changes association with condition: only when few cells
#res <- lfcShrink(dds, coef="experiment_LPS_vs_ctrl", type="apeglm")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.25,
                title = "LPS vs Ctrl - 4h",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
            #    selectLab = c('CXCL8','CCL5','INSIG1',
            #                   'CCL3/L1/L3','ITGB8','G0S2','INHBA','PTGS2','IL6',"STAT4"),
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                labSize = 3,
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-1,1),
                ylim = c(0,25))

####################################################################################
# 2. Classical monocytes : ctrl vs lps at 4h
####################################################################################
sub <- sce[,sce$lvl1 == "class. mono." & sce$experiment != "aCD3aCD28" & sce$timepoint == "4h"] # sub
sub <- logNormCounts(sub) # lognorm

table(sub$experiment)

# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ replicate + experiment) # indicate replicates and exp. columns

dds <- DESeq(dds, betaPrior=FALSE)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="experiment_LPS_vs_ctrl")
# Shrink log fold changes association with condition: only when few cells
#res <- lfcShrink(dds, coef="experiment_LPS_vs_ctrl", type="apeglm")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.10,
                title = "LPS vs Ctrl - 4h",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
                #    selectLab = c('CXCL8','CCL5','INSIG1',
                #                   'CCL3/L1/L3','ITGB8','G0S2','INHBA','PTGS2','IL6',"STAT4"),
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                labSize = 3,
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-1,1),
                ylim = c(0,25))


####################################################################################
# 3. Classical monocytes : 4h vs 24h in LPS
####################################################################################
sub <- sce[,sce$lvl1 == "class. mono." & sce$experiment == "LPS"] # sub
sub <- logNormCounts(sub) # lognorm

table(sub$timepoint,sub$replicate)

#sub <- sub[,sample(colnames(sub),1000)]
# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ replicate + timepoint) # indicate replicates and exp. columns

dds <- DESeq(dds, betaPrior=FALSE, parallel = TRUE, BPPARAM = bp)
#dds <- DESeq(dds, betaPrior=FALSE)

resultsNames(dds) # lists the coefficients
res <- results(dds, name="timepoint_24h_vs_4h")
# Shrink log fold changes association with condition: only when few cells
#res <- lfcShrink(dds, coef="timepoint_24h_vs_4h", type="apeglm")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.25,
                title = "24h vs 4h - LPS",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
                #    selectLab = c('CXCL8','CCL5','INSIG1',
                #                   'CCL3/L1/L3','ITGB8','G0S2','INHBA','PTGS2','IL6',"STAT4"),
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                labSize = 3,
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-1,1),
                ylim = c(0,25))


####################################################################################
# 4. CD8 T effector: 4h vs 24 in aCD3aCD28
####################################################################################
sub <- sce[,sce$lvl1 == "CD4 TN" & sce$experiment == "aCD3aCD28"] # sub
sub <- logNormCounts(sub) # lognorm

table(sub$timepoint,sub$replicate)

#sub <- sub[,sample(colnames(sub),1000)]
# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ replicate + timepoint) # indicate replicates and exp. columns

dds <- DESeq(dds, betaPrior=FALSE, parallel = TRUE, BPPARAM = bp)

resultsNames(dds) # lists the coefficients
res <- results(dds, name="timepoint_24h_vs_4h")
# Shrink log fold changes association with condition: only when few cells
#res <- lfcShrink(dds, coef="timepoint_24h_vs_4h", type="apeglm")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.25,
                title = "24h vs 4h - aCD3aCD28",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
                #    selectLab = c('CXCL8','CCL5','INSIG1',
                #                   'CCL3/L1/L3','ITGB8','G0S2','INHBA','PTGS2','IL6',"STAT4"),
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                labSize = 3,
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-1,1),
                ylim = c(0,25))
