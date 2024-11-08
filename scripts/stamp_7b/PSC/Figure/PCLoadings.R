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
sce <- qread(glue("{res_dir}/PreProcNew.qs"))

pc <- reducedDim(sce,"PCA")

loadings <- attr(pc,"rotation")

# Define the list of PCs
pcs <- lapply(1:3, function(i) {
  pc <- sort(loadings[, i], decreasing = TRUE)
  pc <- c(head(pc, 5), tail(pc, 5))
  pc_df <- as.data.frame(pc)
  pc_df$gene <- rownames(pc_df)
  pc_df$val <- pc_df$pc
  pc_df$gene <- factor(pc_df$gene, levels = rev(pc_df$gene))
  pc_df$pc <- paste0("PC",i)
  pc_df
})


pc_plt <- \(pc){
p <- ggplot(pcs[[pc]], aes(x = val, y = gene, fill = val > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#fe6346", "FALSE" = "#6395eb")) +
  theme_bw(10) +
  labs(x = "Loading", y = "", subtitle = paste0("PC",pc)) +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 2 / 1,
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    text = element_text(color = "black"), 
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 15),
    legend.position = "none",
  )

if (pc == 1) {p <- p + scale_x_continuous(breaks = c(-0.3,0,0.2))}
if (pc == 2) {p <- p + scale_x_continuous(breaks = c(-0.2,0,0.1))}
if (pc == 3) {p <- p + scale_x_continuous(breaks = c(-0.2,0,0.1))}
return(p)
}

pc_plots <- wrap_plots(
  pc_plt(1) + theme(axis.title.x = element_blank()),
  pc_plt(2),
  pc_plt(3) + theme(axis.title.x = element_blank()),
  nrow = 1)

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC_V2/pc_plots.pdf", width = 6, height =4)
pc_plots
dev.off()
