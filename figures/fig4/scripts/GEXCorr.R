## Dependencies
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(glue)
  library(here)
  library(scales)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(dplyr)
})
## Data
cd <- qread("./../obj/cd.qs")
sce <- qread("./../obj/sce.qs")
source(glue("{here()}/scripts/misc/paths.R"))
source(glue("{here()}/scripts/misc/BIN.R"))

cd <- cd[cd$stamp != "stamp_9",]
sce <- sce[,sce$stamp != "stamp_9"]
cd$replicate[cd$stamp == "stamp_11"] <- "replicate 1"
cd$replicate[cd$stamp == "stamp_12"] <- "replicate 2"

sce$replicate[sce$stamp == "stamp_11"] <- "replicate 1"
sce$replicate[sce$stamp == "stamp_12"] <- "replicate 2"

pal <- Polychrome::createPalette(26, c("#A3B8CC", "#550A46"))
val <- unique(cd$replicate)
names(pal) <- val

df <- as.data.frame(table(cd$sub,cd$replicate))

gg_cnumb <- ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.5) +
  geom_text(aes(label = Freq), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, 
            size = 4) +  
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(color = "black"),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(color = "black", size = 16),
        legend.position = "top") +
  labs(x = "", y = "# Cells", fill = "Replicate")


create_boxplot <- function(df, var_name, title, ylab, log) {
  # Dynamically calculate the median and max for each sample
stats <- df %>%
  group_by(replicate, sub) %>%
  summarize(median_value = median(!!sym(var_name)),
            max_value = max(!!sym(var_name)))

# Create the plot
p <- ggplot(df, aes(x = sub, y = !!sym(var_name), color = replicate)) +
  ggrastr::rasterise(geom_boxplot(aes(fill = replicate), alpha = 0.5, outlier.size = 0.1), dpi = 800) +  # Set alpha for fill only
  geom_text(data = stats, aes(x = sub, y = max_value + 0.05 * max_value, label = round(median_value, 2), group = replicate),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 5, color = "black") +  # Add median values on top
  theme_bw() +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(y = ylab, x = "", color = "Median") +
  guides(color = "none", fill = "none") +  # Remove the legend
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 16),
        axis.text = element_text(color = "black", size = 16),
        text = element_text(color = "black", size = 20))
  
  if(log == TRUE) {p <- p + scale_y_log10()}
  return(p)
}




# Example usage
df <- cd
gg_sum <- create_boxplot(df, "sum", "Sample", "nCount", log = TRUE)
gg_feat <- create_boxplot(df, "detected", "Sample", "nFeature", log = TRUE)
gg_area <- create_boxplot(df, "cell_area", "Sample", "Cell area (um2)", log = FALSE)


# Take count matrices
c <- sce[,sce$replicate == "replicate 1"]
n <- sce[,sce$replicate == "replicate 2"]
c <- as(counts(c), "dgCMatrix")
n <- as(counts(n), "dgCMatrix")

# norm by total n cells
c <- rowSums(c) / ncol(n)
n <- rowSums(n) / ncol(n)
c <- as.data.frame(as.matrix(c))
n <- as.data.frame(as.matrix(n))

names(c)[1] <- "replicate_1"
names(n)[1] <- "replicate_2"

c$gene <- rownames(c)
n$gene <- rownames(n)
df <- merge(c, n, by = "gene", all = TRUE)
#df <- df[df$gene != "EEF1G",]


gg_corr <- ggplot(df, aes(x = replicate_1, y = replicate_2)) + 
  ggrastr::rasterise(geom_point(shape = 16, size = 3, alpha = 0.8), dpi = 800) +
  ggpubr::stat_cor(method = "pearson", label.x = -1, label.y = 25, size = 7) + 
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(x ="Mean counts/gene - replicate 1", y ="Mean counts/gene - replicate 2") +
  theme(panel.grid =  element_blank(),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 18, color = "black"))


pdf("/Users/emanuelepitino/Desktop/fig2/h_i.pdf", width = 23, height = 7)
wrap_plots(
  gg_cnumb,
  gg_sum,
  gg_feat,
  gg_area,
  gg_corr,
  nrow = 1
) +
  plot_layout(widths = c(1.5,1,1,1,2))

dev.off()


##

