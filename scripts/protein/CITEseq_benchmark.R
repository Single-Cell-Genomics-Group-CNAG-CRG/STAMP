## correlation between STAMP12c and STAMP11c for the shared proteins

#STAMP11c
cp_counts <- data.table::fread('/mnt/scratch2/STAMP/stamp_11c/STAMPCP_exprMat_file.csv.gz')
cp_counts <- cp_counts[, -c(1:2)] # remove fov and cell id
cp_counts <- cp_counts[, -c(31:35,68:69)] # remove membrane markers and control iGGs
cp_counts <- as.data.frame(t(cp_counts))
# for standarization
rownames(cp_counts) <- toupper(rownames(cp_counts))


# STAMP 12c
stamp12c<- as.data.frame(data.table::fread('/mnt/scratch2/STAMP/STAMP_PCF/24066_segmentation.tsv'))
w <- grep('Cell: Mean', colnames(stamp12c))
stamp12c <- stamp12c[, w]
colnames(stamp12c) <- gsub(': Cell: Mean', '', colnames(stamp12c))
colnames(stamp12c) <- toupper(colnames(stamp12c))
stamp12c <- t(stamp12c)

shared <- intersect(rownames(cp_counts), rownames(stamp12c))

avg_cor <- data.frame(
  protein = shared,
  STAMP11c = rowMeans(cp_counts[shared,]),
  STAMP12c =  rowMeans(stamp12c[shared,])
)

pp <- ggplot(avg_cor, aes(STAMP11c, STAMP12c)) + geom_point() +
  geom_smooth(method = "lm", color = "darkred", se = FALSE) +  # Regression line
  stat_cor(method = "pearson", label.x = min(avg_cor$STAMP11c), label.y = max(avg_cor$STAMP12c)) +
  theme_bw() + xlab('Mean intensity/protein STAMP11c (CosMx)') + ylab('Mean intensity/protein STAMP12c (PhF)') +
  theme(text=element_text(size = 15))
ggsave(plot = pp, filename = 'CorCP-PCF-Fig6.pdf', width = 8, height = 6)

## CITEseq
pbmc_cite <- Matrix::readMM('/mnt/scratch2/STAMP/pbmc_CITEseq/GSM5008738_ADT_3P-matrix.mtx.gz')
pbmc_cite_features <- read.table('/mnt/scratch2/STAMP/pbmc_CITEseq/GSM5008738_ADT_3P-features.tsv.gz')
pbmc_cite_barcodes <- read.table('/mnt/scratch2/STAMP/pbmc_CITEseq/GSM5008738_ADT_3P-barcodes.tsv.gz')
rownames(pbmc_cite) <- pbmc_cite_features$V1
colnames(pbmc_cite) <- pbmc_cite_barcodes$V1
pbmc_cite_met <- read.csv('/mnt/scratch2/STAMP/pbmc_CITEseq/GSE164378_sc.meta.data_3P.csv.gz')
pbmc_cite_met <- pbmc_cite_met[pbmc_cite_met$X %in% colnames(pbmc_cite), ]
rownames(pbmc_cite_met) <- pbmc_cite_met$X

pbmc_cite <- CreateSeuratObject(counts = pbmc_cite, meta.data = pbmc_cite_met, assay = 'ADT')
pbmc_cite <- NormalizeData(pbmc_cite)
pbmc_cite$Tech <- "CITE-Seq"

stamp12_pcf_merged <- readRDS("/mnt/scratch1/STAMP-X/STAMP-Analysis/stamp12_pcf_merged.rds")
stamp12_pcf_merged$Tech <- stamp12_pcf_merged$Sample
STAMP11_b_c_merged_annotated <- readRDS("/mnt/scratch1/STAMP-X/STAMP-Analysis/STAMP11_b_c_merged_annotated.rds")
STAMP11_b_c_merged_annotated$Tech <- STAMP11_b_c_merged_annotated$STAMP

# annotate stamp12b
tmp <- subset(stamp12_pcf_merged, Tech == 'STAMP-PhF')
tmp <- tmp %>% NormalizeData(normalization.method = 'CLR') %>% FindVariableFeatures() %>% ScaleData()

anchors <- FindTransferAnchors(reference = stamp12b_pcf_pbmc,
                               query = tmp,
                               dims = 1:20,  # Adjust based on your data
                               reduction = "rpca")
labels <- TransferData(anchorset = anchors,
                            refdata = stamp12b_pcf_pbmc$Lineage,  # Replace with your metadata column
                            dims = 1:20)
tmp$Lineage <- labels$predicted.id

stamp12b_pcf_pbmc$Tech <- 'STAMP-XPhF'
stamp12_pcf_merged <- merge(stamp12b_pcf_pbmc, tmp)
# At this point I should have the information from all 4 STAMP proteome samples + CITE-seq





# Create data frame with nCount, nFeature, CellType
meta_all_prot <-
  data.frame(nCount = c(pbmc_cite$nCount_ADT, stamp12_pcf_merged$nCount_PCF, STAMP11_b_c_merged_annotated$nCount_CP),
             nFeature = c(pbmc_cite$nFeature_ADT, stamp12_pcf_merged$nFeature_PCF, STAMP11_b_c_merged_annotated$nFeature_CP),
             Sample = c(pbmc_cite$Tech,stamp12_pcf_merged$Tech, STAMP11_b_c_merged_annotated$Tech),
             CellType = c(pbmc_cite$celltype.l1, stamp12_pcf_merged$Lineage,STAMP11_b_c_merged_annotated$CellType_pred ))

meta_all_prot$PanelSize <- dim(STAMP11_b_c_merged_annotated)[1]
meta_all_prot$PanelSize[meta_all_prot$Sample %in% 'CITE-Seq'] <- dim(pbmc_cite)[1]
meta_all_prot$PanelSize[meta_all_prot$Sample %in% c('STAMP-PhF', 'STAMP-XPhF')] <- dim(stamp12_pcf_merged)[1]

# Change names of cell types to match my annotation
cite_cellt <- as.factor(meta_all_prot$CellType[meta_all_prot$Sample %in% 'CITE-Seq'])
levels(cite_cellt) <- c('B', 'T','T', 'Myeloid', 'Myeloid', 'NK', 'other', 'T')
cite_cellt <- as.character(cite_cellt)
meta_all_prot$CellType[meta_all_prot$Sample %in% 'CITE-Seq'] <- cite_cellt

cp_cellt <- as.factor(meta_all_prot$CellType[meta_all_prot$Sample %in% c('STAMP11b', 'STAMP11c')])
levels(cp_cellt) <- c('B', 'T', 'T', 'Myeloid', 'Myeloid', 'NK', 'other', 'other')
cp_cellt <- as.character(cp_cellt)
meta_all_prot$CellType[meta_all_prot$Sample %in% c('STAMP11b', 'STAMP11c')] <- cp_cellt

#remove lowQ and other
meta_all_prot <- meta_all_prot[!meta_all_prot$CellType %in% c('other', 'LowQ'), ]
#rename samples
meta_all_prot$Sample <- as.factor(meta_all_prot$Sample)
levels(meta_all_prot$Sample)[2:5] <- c('STAMP-PCF', 'STAMP-X-PCF', 'STAMP-X-CP', 'STAMP-CP')
meta_all_prot$Sample <- as.character(meta_all_prot$Sample)

p1 <- ggplot(meta_all_prot, aes(Sample, nCount, col = Sample)) + geom_boxplot() + scale_y_log10() + facet_wrap(~ CellType) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),  legend.position = "top",
                     legend.title = element_blank() , text=element_text(size = 16))  +
  scale_color_manual(values = awtools::a_pal()(5)) + xlab('')

#separate CITE seq and STAMP
cols = awtools::a_pal()(5)


p1 <- ggplot(meta_all_prot[meta_all_prot$Sample == 'CITE-Seq',], aes(Sample, nCount, col = Sample)) + geom_boxplot() + scale_y_log10() + facet_wrap(~ CellType) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),  legend.position = "top",
                     legend.title = element_blank() , text=element_text(size = 16))  +
  scale_color_manual(values = cols[1]) + xlab('')

p1.1 <- ggplot(meta_all_prot[meta_all_prot$Sample != 'CITE-Seq',],
               aes(x = factor(Sample, levels = c("STAMP-X-CP", "STAMP-CP", "STAMP-X-PCF", "STAMP-PCF")),
                   y = nCount,
                   col = Sample)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ CellType) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(size = 16)
  ) +
  scale_color_manual(values = cols[2:5]) +
  guides(color = guide_legend(nrow = 2)) +  # Arrange legend in 2 rows
  xlab('') +
  ylab('Fluorescence Intensity')


p2 <- ggplot(meta_all_prot, aes(Sample, nFeature, col = Sample)) + geom_boxplot() + facet_wrap(~ CellType) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),  legend.position = "none",
                     legend.title = element_blank() , text=element_text(size = 16))  +
  scale_color_manual(values = awtools::a_pal()(5)) + xlab('')

p3 <- ggplot(meta_all_prot, aes(Sample, nFeature / PanelSize, col = Sample)) + geom_boxplot() + facet_wrap(~ CellType) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),  legend.position = "none",
                     legend.title = element_blank() , text=element_text(size = 16))  +
  scale_color_manual(values = awtools::a_pal()(5)) + xlab('')

p1 / p2 / p3


# prepare and plot genes of interest as barplots of avg expression

#normalize the citeseq the same way
pbmc_cite <- NormalizeData(pbmc_cite, normalization.method = 'CLR')


genes_of_interest <- c("CD4", "CD4-1", "CD8", "CD8-1", 'CD14', 'CD20')

# Function to get average expression from a Seurat object
get_avg_expr <- function(seurat_obj, sample_name, assay = 'PCF', group = '') {
  avg_expr <- AverageExpression(seurat_obj, features = genes_of_interest, assays = assay, group.by = group, slot = 'counts')
  avg_expr <- avg_expr[[assay]]
  avg_expr <- as.data.frame(avg_expr)
  avg_expr$Gene <- rownames(avg_expr)
  avg_expr$Sample <- sample_name
  return(avg_expr)
}

# Compute for each Seurat object
pbmc_cite$Lineage <- meta_all_prot$CellType[meta_all_prot$Sample == 'CITE-Seq']

avg_expr_1 <- get_avg_expr(pbmc_cite, "CITE-Seq", assay = 'ADT', group = 'Lineage')
avg_expr_1$Gene[1] <- 'CD4'
stamp12_pcf_merged <- JoinLayers(stamp12_pcf_merged)
avg_expr_2 <- get_avg_expr(subset(stamp12_pcf_merged, Tech == 'STAMP-XPhF') , sample_name = "STAMP-X-PCF",
                           group = 'Lineage', assay = 'PCF')
avg_expr_3 <- get_avg_expr(subset(stamp12_pcf_merged, Tech == 'STAMP-PhF') , sample_name = "STAMP-PCF",
                           group = 'Lineage', assay = 'PCF')
#remove LowQ
avg_expr_2 <- avg_expr_2[, -2]; avg_expr_3 <- avg_expr_3[, -2];

STAMP11_b_c_merged_annotated$Lineage <-  meta_all_prot$CellType[meta_all_prot$Sample %in% c('STAMP-X-CP', 'STAMP-CP"')]

avg_expr_4 <- get_avg_expr(subset(STAMP11_b_c_merged_annotated, Tech == 'STAMP11b') , sample_name = "STAMP-X-CP",
                           group = 'Lineage', assay = 'CP')
avg_expr_5 <- get_avg_expr(subset(STAMP11_b_c_merged_annotated, Tech == 'STAMP11c') , sample_name = "STAMP-CP",
                           group = 'Lineage', assay = 'CP')




# Combine into one data frame
avg_expr_combined <- bind_rows(avg_expr_1, avg_expr_2, avg_expr_3, avg_expr_4,avg_expr_5) %>%
  pivot_longer(cols = -c(Gene, Sample), names_to = "CellType", values_to = "AvgExpression")

p4 <- ggplot(avg_expr_combined[avg_expr_combined$Sample == 'CITE-Seq',], aes(Sample, AvgExpression, fill = Sample)) + geom_bar(stat='identity') + facet_wrap(~Gene) +
  scale_y_log10() + scale_fill_manual(values = cols[1]) +  theme_bw() + xlab('') + ylab('log10 counts') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),  legend.position = "none",
        legend.title = element_blank() , text=element_text(size = 16))

p4.1 <- ggplot(avg_expr_combined[avg_expr_combined$Sample != 'CITE-Seq',], aes(Sample, AvgExpression, fill = Sample)) + geom_bar(stat='identity') + facet_wrap(~Gene) +
  scale_y_log10() + scale_fill_manual(values = cols[2:5]) +  theme_bw() + xlab('') + ylab('log10 Fluorescence Intensity') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),  legend.position = "none",
        legend.title = element_blank() , text=element_text(size = 16))





# Combine p1, p2, and p3 in one column
left_panel <- ((p1 + p1.1) / p2 / p3)

((p1 + p1.1) / (p2 / p3)+ (p4+p4.1)   )





p1 <- ggplot(meta_all_prot[meta_all_prot$Sample != 'CITE-Seq',],
             aes(CellType, nCount, col = Sample)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "top"  # Move legend to the top
  ) +
  scale_color_manual(values = cols[2:5]) +
  xlab('') +
  guides(color = guide_legend(nrow = 2))
p1.1 <- ggplot(meta_all_prot[meta_all_prot$Sample %in% 'CITE-Seq',], aes(CellType, nCount, col = Sample)) + geom_boxplot() +
  theme_bw() + theme( text=element_text(size = 16), legend.position = "top")  +
  scale_color_manual(values = awtools::a_pal()(5)) + xlab('')




p2 <- ggplot(meta_all_prot[meta_all_prot$Sample !='CITE-Seq',], aes(CellType, nFeature, col = Sample)) + geom_boxplot() +
  theme_bw() + theme( text=element_text(size = 16), legend.position = 'none')  +
  scale_color_manual(values = cols[2:5]) + xlab('')

p2.1 <- ggplot(meta_all_prot[meta_all_prot$Sample %in% 'CITE-Seq',], aes(CellType, nFeature, col = Sample)) + geom_boxplot() +
  theme_bw() + theme( text=element_text(size = 16), legend.position = 'none')  +
  scale_color_manual(values = awtools::a_pal()(5)) + xlab('')


p3 <- ggplot(meta_all_prot, aes(CellType, nFeature / PanelSize, col = Sample)) + geom_boxplot()  +
  theme_bw() + theme( text=element_text(size = 16), legend.position = 'none')  +
  scale_color_manual(values = awtools::a_pal()(5)) + xlab('')

left_panel <- ((p1 + p1.1) / (p2+p2.1) / p3)

ggsave(plot= p1, filename = 'test.png', width = 15, height = 10)


# Combine left_panel and p4 in a row, with widths 80% and 20%
combined_plot <- ( left_panel | (p4 + p4.1) ) + plot_layout(widths = c(5, 3)) +
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = "")
combined_plot
ggsave(plot = combined_plot, filename  = 'Benchmark_proteome_3.png', width = 16, height = 10, dpi = 300 )








