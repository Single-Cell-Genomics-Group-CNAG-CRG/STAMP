# dependencies
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(glue)
  library(qs)
  library(ggpubr)
})

# setup
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_3"
dir <- glue("{here()}")


# read cd
dir <- glue("{proj_dir}/data/{stamp}/processed/combined")
fl <- list.files(dir, full.names = T)
cd <- lapply(fl, qread)

# filter
columns_to_keep <- c("sample", "sum", "detected", "fov", "Area.um2", "label", "CenterX_global_px","CenterY_global_px")
cd <- lapply(cd, function(x) {
  x[, columns_to_keep, drop = FALSE]
})
cd <- do.call(rbind,cd)

#  General palette

cd$label[cd$label == "SKBR3"] <- "SK-BR-3"
cd$label[cd$label == "MCF7"] <- "MCF-7"

pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
val <- unique(cd$label)
names(pal) <- val

cd$sample <- recode(cd$sample, "c100" = "100 C", "c250" = "250 C", "c500" = "500 C", "c1000" = "1000 C",
                    "c20k" = "20K C", "n20k" = "20K N")



psamp <- Polychrome::createPalette(10, c("#a6dba0", "#c2a5cf", "#f4a582"))
names(psamp) <- unique(cd$sample)

cdfull <- cd
###### #### ###### #### ###### #### ###### #### ###### #### ###### #### ###### ##


#cd <- cd[!cd$sample %in% c("20K C", "20K N"), ]



# CELL position
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #
df <- cd

# Invert the levels of the sample factor
#df$sample <- factor(df$sample, levels = rev(levels(factor(df$sample))))
df$sample <- factor(df$sample, levels = c("20K C","20K N","1000 C","500 C","250 C","100 C"))
# Custom labels
# Plot
gg_pos <- ggplot(df, aes(x = CenterX_global_px, y = CenterY_global_px, color = sample)) + 
  ggrastr::rasterise(geom_point(shape = 16, size = 0.8), dpi = 800) +
  scale_color_manual(values = psamp) +
  theme_bw() +
  labs(color = "sub-STAMP", x = "x (px)", y = "y (px)") +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme(text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          panel.grid = element_blank()) +
    coord_equal()

dir <- glue("/Users/emanuelepitino/Desktop/fig2")

pdf(glue("{dir}/gg_pos.pdf"), width = 5, height = 4)
gg_pos
dev.off()

# CELL NUMBERS
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
cd$sample <- factor(cd$sample, levels = c("100 C","250 C","500 C","1000 C", "20K C", "20K N"))
df <- as.data.frame(table(cd$sample, cd$label)) %>%
  group_by(Var1) %>%
  mutate(Proportion = round((Freq / sum(Freq)) * 100, 2))

gg_prop <- ggplot(df, aes(x = Var1, y = Proportion)) + 
  geom_col(aes(fill = Var2), alpha = 0.9) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size = 25, color = "black"),
        axis.text = element_text(size = 18,color = "black"),
        axis.text.x = element_text( angle = 45, hjust = 1), 
        legend.position = "right",
        panel.grid = element_blank()) +
  labs(fill = "Cluster ID", y = "Proportions", x= "")  +
  scale_y_continuous(expand = c(0, 0))
  
  pdf(glue("{dir}/gg_prop.pdf"), width = 5, height = 6)
  gg_prop
  dev.off()
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

# CORRELATION OF CELL NUMBERS
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
df <- as.data.frame(table(cd$sample)) %>%
      mutate(InputCells = c(100,225,500,1000, 20000,20000))
avg <- round(mean(df$pct_diff),2)


df <- df %>%
  mutate(pct_diff = (round((Freq/InputCells),2))*100)


gg_cnumb <- ggplot(df, aes(x = Var1, y = pct_diff, fill = Var1)) +
  geom_col() +
  scale_fill_manual(values = psamp) +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3, color = "black") +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x = "",y = "% recovered cells") + 
  theme(legend.position = "none", 
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 18, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  scale_y_continuous(expand = c(0,5))




pdf(glue("{dir}/gg_cnumb.pdf"), width = 3, height = 6)
gg_cnumb
dev.off()

pdf(glue("{dir}/b_c.pdf"), width = 10, height = 6)
wrap_plots(
  gg_cnumb,
  gg_prop,
  nrow = 1
) 
dev.off()

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

# GENE COUNTS & Features
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
df <- cd
# plot function
create_boxplot <- function(df, var_name, title, ylab, log) {
  # Dynamically calculate the median and max for each sample
  stats <- df %>%
    group_by(sample) %>%
    summarize(median_value = median(!!sym(var_name)),
              max_value = max(!!sym(var_name)))
  
  # Create the plot
  p <- ggplot(df, aes(x = sample, y = !!sym(var_name), color = sample)) +
    ggrastr::rasterise(geom_boxplot(aes(fill = sample), alpha = 0.5, outlier.size = 0.1), dpi = 800) +  # Set alpha for fill only
    theme_bw() +
    scale_color_manual(values = psamp) +
    scale_fill_manual(values = psamp) +
    labs(y = ylab, x = "", color = "Median") +
    guides(color = "none", fill = "none") +  # Remove the legend
    geom_text(data = stats, aes(x = sample, y = max_value, label = round(median_value, 0)),
              vjust = -0.5, size = 5, color = "black") +  # Place median above the max value
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
          axis.text = element_text(color ="black", size = 16),
          text = element_text(color = "black", size = 20))
  
  if(log == TRUE) {p <- p + scale_y_log10()}
  return(p)
}

# Example usage
gg_sum <- create_boxplot(df, "sum", "Sample", "nCount", log = TRUE)
gg_feat <- create_boxplot(df, "detected", "Sample", "nFeature", log = TRUE)
gg_area <- create_boxplot(df, "Area.um2", "Sample", "Cell area (um2)", log = FALSE)


pdf(glue("{dir}/gg_cfa.pdf"), width = 12, height = 6)
wrap_plots(
  gg_sum,
  gg_feat,
  gg_area,
  nrow = 1)
dev.off()

####### ####### ####### ####### ####### ####### ####### ####### ####### ######## 
####### ####### ####### ####### ####### ####### ####### ####### ####### ########
######## ####### ####### ####### ## PART 2 ########### ####### ####### #########  
####### ####### ####### Visualization Cell/nuclei ####### ####### ####### ######
####### ####### ####### ####### ####### ####### ####### ####### ####### ########

# GEX correlation
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
res_dir <- glue("{proj_dir}/data/Stamp_3/processed/")
n <- qread(glue("{res_dir}/n20k/qc_n20k.qs"), nthreads = 8)
c <- qread(glue("{res_dir}/c20k/qc_c20k.qs"), nthreads = 8)

# Take count matrices
c <- as(counts(c), "dgCMatrix")
n <- as(counts(n), "dgCMatrix")

# norm by total n cells
c <- rowSums(c) / ncol(n)
n <- rowSums(n) / ncol(n)
c <- as.data.frame(as.matrix(c))
n <- as.data.frame(as.matrix(n))

names(c)[1] <- "Cells"
names(n)[1] <- "Nuclei"

c$gene <- rownames(c)
n$gene <- rownames(n)
df <- merge(c, n, by = "gene", all = TRUE)


gg_corr <- ggplot(df, aes(x = Cells, y = Nuclei)) + 
ggrastr::rasterise(geom_point(shape = 16, size = 3, alpha = 0.8), dpi = 800) +
  ggpubr::stat_cor(method = "pearson", label.x = -3, label.y = 50, size = 7) + 
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(x ="Mean counts/gene - 20K cells", y ="Mean counts/gene - 20K nuclei") +
  theme(panel.grid =  element_blank(),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 18, color = "black"))
  

pdf(glue("{dir}/gg_corr.pdf"), width = 5, height = 5)
gg_corr
dev.off()

pdf(glue("{dir}/d_e.pdf"), width = 18, height = 6)
wrap_plots(
  gg_sum,
  gg_feat,
  gg_area,
  gg_corr,
  ncol = 4) + 
  plot_layout(widths = c(1,1,1,1.5))
dev.off()
 
