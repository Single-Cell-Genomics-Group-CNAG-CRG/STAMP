library(dplyr)
library(patchwork)
library(ggplot2)
library(qs)
library(here)

stamp <- "stamp_7b"
sample <- "iPSCs"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# CosMx data
res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
cosmx <- qread(glue("{res_dir}/anno_sce.qs"))
cosmx$tech <- "cosmx"
# Flex data
sample <- "combined"
outdir <- glue("{proj_dir}/data/{stamp}/processed/flex/iPSC/{sample}")
flex <- qread(glue("{outdir}/anno_sce.qs"), nthreads = 8)
flex$tech <- "flex"


##  marker genes
mrk <- c("INHBB", "NPR1","CMKLR1", "IL1A", "CD69", "FAP",
  "TNFSF8", "ITK", "CD40LG","C5AR2",
  "NRG1", "PTGS1", "EPHA4", "BMP3", "NR2F2",
  "IGFBP7", "APOA1", "PPARG", "TTR", "CACNA1C",
  "WNT9A", "LAMA4", "COL5A3","LY75", "LDHA",
  "PDGFB","DHRS2", "SPOCK2", "MST1R", "ACTG2",
  "HLA-DPB1", "TFEB", "INS", "KLRF1", "CASR",
  "CCL3/L1/L3","LIF", "FGR", "IL34", "IL18R1")

### CORR SPLIT BY CELL TYPE
# Take matrices
c <- cosmx
f <- flex
cs <- as(counts(c),"dgCMatrix")
fl <- as(counts(f), "dgCMatrix")

feat <- intersect(rownames(cs),rownames(fl)) # take features intersection
fl <- fl[feat, ] # subset 
cs <- cs[feat,] # subset

# Aggregate the flex  matrix
df <- as.data.frame(colData(f))
colnames(fl) <- df[colnames(fl), "cluster"]
#fl <- t(rowsum(t(fl), group = colnames(fl)))
fl <- as(fl, "dgCMatrix")
# Aggregate the cosmx matrix
df <- as.data.frame(colData(c))
colnames(cs) <- df[colnames(cs), "cluster"]
#cs <- t(rowsum(t(cs), group = colnames(cs)))
cs <- as(cs, "dgCMatrix")

# function to aggregate and normalize
aggr_norm <- function(mat){
  res <- do.call(cbind, lapply(unique(colnames(mat)), function(col_name) {
    sub <- mat[, colnames(mat) %in% col_name, drop = FALSE]
    agg <- rowSums(sub) / ncol(sub)
    agg <- rowSums(sub)
    agg
  }))
  colnames(res) <- unique(colnames(mat))
  return(res)
}

agg_fl <- as.data.frame(aggr_norm(fl))
agg_fl$gene <- rownames(agg_fl)
agg_cs <- as.data.frame(aggr_norm(cs))
agg_cs$gene <- rownames(agg_cs)


agg_cs <- melt(agg_cs, id.vars = "gene")
names(agg_cs)[names(agg_cs) == "value"] <- "cosmx"

agg_fl <- melt(agg_fl, id.vars = "gene")
names(agg_fl)[names(agg_fl) == "value"] <- "flex"



df <- merge(agg_cs, agg_fl, by = c("gene", "variable"), all = TRUE)


corr_plot <- function(df,subset){
  
  genes_p1 <- c("INHBB", "NPR1","CMKLR1", "IL1A", "CD69")
  genes_p2 <- c("FAP","TNFSF8", "ITK", "CD40LG","C5AR2")
  genes_ec1 <- c("NRG1", "PTGS1", "EPHA4", "BMP3", "NR2F2")
  genes_m1 <- c("IGFBP7", "APOA1", "PPARG", "TTR", "CACNA1C")
  genes_m2 <- c("WNT9A", "LAMA4", "COL5A3","LY75", "LDHA")
  genes_n1 <- c("PDGFB","DHRS2", "SPOCK2", "MST1R", "ACTG2")
  genes_n2 <- c("HLA-DPB1", "TFEB", "INS", "KLRF1", "CASR")
  genes_n3 <- c("CCL3/L1/L3","LIF", "FGR", "IL34", "IL18R1")
  
  
  if(subset == "P1"){ feat_color = genes_p1}
  if(subset == "P2"){ feat_color = genes_p2}
  if(subset == "EC1"){ feat_color = genes_ec1}
  if(subset == "M1"){ feat_color = genes_m1}
  if(subset == "M2"){ feat_color = genes_m2}
  if(subset == "EN1"){ feat_color = genes_n1}
  if(subset == "EN2"){ feat_color = genes_n2}
  if(subset == "EN3"){ feat_color = genes_n3}
  
  sub <- df[df$variable == subset,]
  # Create a new variable to indicate if a gene is a marker
  sub$Is_Marker <- ifelse(sub$gene %in% feat_color, "Yes", "No")
  
  # Separate the data
  sub_no <- sub[sub$Is_Marker == "No", ]
  sub_yes <- sub[sub$Is_Marker == "Yes", ]
  # Plot
  plt <- ggplot(sub, aes(x = flex, y = cosmx)) +
    # Plot "No" points first
    ggrastr::rasterise(geom_point(data = sub_no, aes(color = Is_Marker),
                                  shape = 16, size = 2, alpha = 0.8), dpi = 800) +
    ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + 
    # Plot "Yes" points on top
    geom_point(data = sub_yes, aes(color = Is_Marker), shape = 16, size = 2, alpha = 0.8) +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(slope = 1, intercept = 0, color = "red4", linetype = "dashed") +
    labs(title = subset, color = "Is Marker") +
    theme_bw() +
    theme(
      text = element_text(size = 10, color = "black"), 
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      panel.grid = element_blank(),
      legend.position = "none"  # Remove legend for combining later
    ) +
    labs(x = "Mean counts/gene - Flex", y= "Mean counts/gene - CosMx") +
    scale_color_manual(values = c("Yes" = "red", "No" = "grey"))
  
  return(plt)
}


gg_corr <- wrap_plots(
  corr_plot(df,"P1"),
  corr_plot(df,"P2"),
  corr_plot(df,"EC1"),
  corr_plot(df,"EN1"),
  corr_plot(df,"EN2"),
  corr_plot(df,"EN3"),
  ncol = 3) +
  plot_layout(axis_titles = "collect")
#gg_corr


pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/gex_corr.pdf", width = 6, height = 4)
gg_corr
dev.off()



# CosMx data
res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
cosmx <- qread(glue("{res_dir}/anno_sce.qs"))
cosmx$tech <- "cosmx"

sample <- "combined"
outdir <- glue("{proj_dir}/data/{stamp}/processed/flex/iPSC/{sample}")
flex <- qread(glue("{outdir}/anno_sce.qs"), nthreads = 8)
flex$tech <- "flex"


cd_f <- as.data.frame(colData(flex))
cd_c <- as.data.frame(colData(cosmx))

cd_c <- cd_c %>% select(tech,cluster)
cd_f <- cd_f %>% select(tech,cluster)

cd <- rbind(cd_c,cd_f)

df <- as.data.frame(table(cd$tech,cd$cluster))
df <- df %>% group_by(Var1) %>% mutate(perc = Freq/sum(Freq)) %>% ungroup()

cnumb <- ggplot(df, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 12,color = "black"),
    axis.text.y = element_text(size = 12, color ="black"),
    axis.title = element_text(size = 18),
    legend.position = "top") +
  labs(x = "Cluster", y = "# Cells", fill = "Tech")

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/c_d.pdf", width = 14, height = 6)
wrap_plots(cnumb,gg_corr, ncol = 2) + plot_layout(widths = c(1,2.5))
dev.off()
