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
sce$id <- paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)

agg <- aggregateAcrossCells(sce, ids = sce$id, use.assay.type = "counts", statistics = "sum")

sizeFactors(agg) <- NULL

sub <- agg[,agg$lvl1 == "act. mono." & agg$experiment != "aCD3aCD28" & agg$timepoint == "4h"] # sub

# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ experiment) # indicate replicates and exp columns

dds <- DESeq(dds, BPPARAM = bp)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="experiment_LPS_vs_ctrl")
# Shrink log fold changes association with condition: only when few cells
#res <- lfcShrink(dds, coef="experiment_LPS_vs_ctrl", type="apeglm")
act_mono_lps_ctrl_4h <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                title = "LPS vs Ctrl - 4h",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
                selectLab = c('CXCL8','CCL5','CSF3',
                               'CCL3/L1/L3','G0S2','PTGS2','IL6','IL1B','CDKN1A'),
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                labSize = 3,
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-1,1),
                ylim = c(0,15)) 
act_mono_lps_ctrl_4h
####################################################################################
# 2. Activated monocytes : 4h vs 24h in LPS
####################################################################################
sub <- agg[,agg$lvl1 == "act. mono." & agg$experiment == "LPS"] # sub

table(sub$timepoint,sub$replicate)

# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ timepoint) # indicate replicates and exp. columns

dds <- DESeq(dds, BPPARAM = bp)

resultsNames(dds) # lists the coefficients
res <- results(dds, name="timepoint_24h_vs_4h")

act_mono_4h_24h_lps <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 0.25,
                title = "24h vs 4h - LPS",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                labSize = 3,
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-1,1),
                ylim = c(0,5)) 
act_mono_4h_24h_lps
####################################################################################
# 3. Classical monocytes : ctrl vs lps at 4h
####################################################################################
sub <- agg[,agg$lvl1 == "class. mono." & agg$experiment != "aCD3aCD28" & agg$timepoint == "4h"] # sub

table(sub$experiment)

# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ experiment) # indicate replicates and exp. columns

dds <- DESeq(dds, BPPARAM = bp)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="experiment_LPS_vs_ctrl")

class_mono_lps_ctrl_4h <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 0.10,
                pCutoffCol = 'padj',
                title = "LPS vs Ctrl - 4h",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
                    selectLab = c('CXCL8','CCL5','INSIG1',
                                   'CCL3/L1/L3','ITGB8','G0S2','INHBA','PTGS2','IL6',"STAT4"),
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
class_mono_lps_ctrl_4h
####################################################################################
# 4. Classical monocytes : 4h vs 24h in LPS
####################################################################################
sub <- agg[,agg$lvl1 == "class. mono." & agg$experiment == "LPS"] # sub

table(sub$timepoint,sub$replicate)

# Create obj
dds <- DESeqDataSetFromMatrix(countData = counts(sub), # take matrix
                              colData = colData(sub), # take coldata
                              design= ~ timepoint) # indicate replicates and exp. columns

dds <- DESeq(dds, BPPARAM = bp)

resultsNames(dds) # lists the coefficients
res <- results(dds, name="timepoint_24h_vs_4h")
# Shrink log fold changes association with condition: only when few cells
#res <- lfcShrink(dds, coef="timepoint_24h_vs_4h", type="apeglm")
class_mono_4h_24h_lps <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                title = "24h vs 4h - LPS",
                subtitle = glue("subset: {unique(sub$lvl1)}"),
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
class_mono_4h_24h_lps

####################################################################################
pdf("/Users/emanuelepitino/Desktop/stamp_13a/volcanos.pdf", width = 5, height = 10)
wrap_plots(act_mono_lps_ctrl_4h,act_mono_4h_24h_lps, ncol = 1) + 
  plot_layout(guides = "collect") &
  theme_bw(18) &
  theme(panel.grid = element_blank(), legend.position = "none") & 
  scale_color_manual(values = c("#cccccc","#77dd77","#ff6961"))
dev.off()

