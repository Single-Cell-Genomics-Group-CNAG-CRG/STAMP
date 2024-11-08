# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(glue)
library(qs)
library(here)

set.seed(123)
pal <- Polychrome::createPalette(31, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V",
                "W","X","Y","MX1","Q","R",
                "stamp_15","stamp_16","stamp_17","stamp_18")
pal[28] <- "#A6CEE3"
pal[29] <- "#1F78B4" 
pal[30] <- "#B2DF8A" 
pal[31] <- "#33A02C"

stamps <- c("stamp_17","stamp_18","stamp_15","stamp_16")
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/clust_sce.qs"), nthreads = 8)
s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/clust_sce.qs"), nthreads = 8)
s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/qc_sce.qs"), nthreads = 8)
s16 <- qread(glue("{proj_dir}/data/{stamps[4]}/processed/qc_sce.qs"), nthreads = 8)

qcmet_s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/qcmet_df.qs"))
qcmet_s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/qcmet_df.qs"))
qcmet_s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/qcmet_df.qs"))
qcmet_s16 <- qread(glue("{proj_dir}/data/{stamps[4]}/processed/qcmet_df.qs"))

qcmet_s17$stamp <- stamps[1]
qcmet_s18$stamp <- stamps[2]
qcmet_s15$stamp <- stamps[3]
qcmet_s16$stamp <- stamps[4]

qcmet_df <- rbind(qcmet_s17,qcmet_s18,qcmet_s15,qcmet_s16)
s17$stamp <- stamps[1]
s18$stamp <- stamps[2]
s15$stamp <- stamps[3]
s16$stamp <- stamps[4]


cd17 <- as.data.frame(colData(s17)) %>% select(stamp,sample,sum,detected,cell_area)
cd18 <- as.data.frame(colData(s18)) %>% select(stamp,sample,sum,detected,cell_area)
cd15 <- as.data.frame(colData(s15)) %>% select(stamp,sample,sum,detected,Area.um2)
cd16 <- as.data.frame(colData(s16)) %>% select(stamp,sample,sum,detected,Area.um2)

colnames(cd15)[colnames(cd15) == "Area.um2"] <- "cell_area"
colnames(cd16)[colnames(cd16) == "Area.um2"] <- "cell_area"

cd <- rbind(cd17,cd18,cd15,cd16)


mean_df <- cd %>%
  group_by(sample, stamp) %>%
  summarise(
    mean_sum = mean(sum),
    mean_detected = mean(detected),
    mean_carea = mean(cell_area)) %>%
  ungroup()
mean_df <- as.data.frame(mean_df)

mean_df$tech <- ifelse(mean_df$stamp %in% c("stamp_15","stamp_16"),"CosMx","Xenium")

# qcmet corr function
qcmet_corr <- \(x, y,qcmet_var) {
  # subset
  plot_data <- mean_df %>%
    filter(stamp %in% c(x, y)) %>%
    select(qcmet_var, sample, stamp) %>%
    pivot_wider(names_from = stamp, values_from = qcmet_var)
  
  # Plot
  ggplot(plot_data, aes(x = .data[[x]], y = .data[[y]])) +
        geom_smooth(method = "lm", color = "black", size = 0.8) +
        geom_point(aes(color = sample)) +
        scale_color_manual(values = pal) +
        theme_bw() +
        coord_equal() +
        ggpubr::stat_cor(
          method = "pearson", 
          label.x.npc = "left", 
          label.y.npc = "top", 
          size = 3
        ) +
        theme(legend.position = "none",
              panel.grid = element_blank(),
              text = element_text(size = 15, color = "black"),
              axis.text = element_text(color ="black", size = 10),
              axis.title  = element_text(size = 13, color = "black"),
              aspect.ratio = 1) 
  
}
# CosMx vs Xenium qcmet corr function
agg_mean_df <- mean_df %>%
  group_by(sample, tech) %>%
  summarize(
    mean_sum = mean(mean_sum, na.rm = TRUE),
    mean_detected = mean(mean_detected, na.rm = TRUE),
    mean_carea = mean(mean_carea, na.rm = TRUE)
  ) %>%
  ungroup()

cosmx_xenium_corr <- function(value_col) {
  # Aggregate the specified value column by sample and tech
  agg_mean_df <- mean_df %>%
    group_by(sample, tech) %>%
    summarize(value = mean(.data[[value_col]], na.rm = TRUE)) %>%
    ungroup()
  
  # Reshape for plotting
  plot_data <- agg_mean_df %>%
    select(value, sample, tech) %>%
    pivot_wider(names_from = tech, values_from = value)
  
  # Plot
  ggplot(plot_data, aes(x = CosMx, y = Xenium)) +
    geom_smooth(method = "lm", color = "black", size = 0.8) +
    geom_point(aes(color = sample)) +
    scale_color_manual(values = pal) +
    theme_bw() +
    ggpubr::stat_cor(
      method = "pearson", 
      label.x.npc = "left", 
      label.y.npc = "top", 
      size = 3
    ) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          text = element_text(size = 15, color = "black"),
          axis.text = element_text(color ="black", size = 10),
          axis.title  = element_text(size = 13, color = "black"),
          aspect.ratio = 1) 
}


# Plot
counts <- wrap_plots(
  qcmet_corr("stamp_15", "stamp_16","mean_sum") + labs(x = "STAMP-C #15", y = "STAMP-C #16"),
  qcmet_corr("stamp_17", "stamp_18","mean_sum") + labs(x = "STAMP-X #17", y = "STAMP-X #18"),
  cosmx_xenium_corr("mean_sum"),
  ncol = 1) +
  plot_annotation(title = "nCount") & 
  theme(plot.title = element_text(hjust = 0.5, size = 18))
counts

# Run the function
features <- wrap_plots(
  qcmet_corr("stamp_15", "stamp_16","mean_detected") + labs(x = "STAMP-C #15", y = "STAMP-C #16"),
  qcmet_corr("stamp_17", "stamp_18","mean_detected") + labs(x = "STAMP-X #17", y = "STAMP-X #18"),
  cosmx_xenium_corr("mean_detected"),
  ncol = 1) +
  plot_annotation(title = "nFeatures") & 
  theme(plot.title = element_text(hjust = 0.5, size = 18))
features

carea <- wrap_plots(
  qcmet_corr("stamp_15", "stamp_16","mean_carea") + labs(x = "STAMP-C #15", y = "STAMP-C #16"),
  qcmet_corr("stamp_17", "stamp_18","mean_carea") + labs(x = "STAMP-X #17", y = "STAMP-X #18"),
  cosmx_xenium_corr("mean_carea"),
  ncol = 1) +
  plot_annotation(title = "Cell Area (um2)") & 
  theme(plot.title = element_text(hjust = 0.5, size = 18))
carea


pdf("/Users/emanuelepitino/Desktop/high_multi/counts_corr.pdf", width = 4, height = 7)
counts
dev.off()

pdf("/Users/emanuelepitino/Desktop/high_multi/feats_corr.pdf", width = 4, height = 7)
features
dev.off()

pdf("/Users/emanuelepitino/Desktop/high_multi/carea_corr.pdf", width = 4, height = 7)
carea
dev.off()


