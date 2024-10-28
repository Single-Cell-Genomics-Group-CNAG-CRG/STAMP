# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

set.seed(123)
pal <- Polychrome::createPalette(31, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V",
                "W","X","Y","MX2","Q","R",
                "stamp_15","stamp_16","stamp_17","stamp_18")
pal[28] <- "#A6CEE3"
pal[29] <- "#1F78B4" 
pal[30] <- "#B2DF8A" 
pal[31] <- "#33A02C"

stamps <- c("stamp_17","stamp_18","stamp_15","stamp_16")
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/clust_sce.qs"))
s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/clust_sce.qs"))
s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/qc_sce.qs"))

s16 <- s15 # mock until s16 is not available
qcmet_s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/qcmet_df.qs"))
qcmet_s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/qcmet_df.qs"))
qcmet_s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/qcmet_df.qs"))
qcmet_s16 <- qcmet_s15 # mock until s16 is not available

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
    mean_carea = mean(cell_area))

mean_bxp <- \(var){
  
set.seed(124)
  if (var == "mean_sum") {ylab = "nCount"}
  if (var == "mean_detected") {ylab = "nFeature"}
  if (var == "mean_carea") {ylab = "Cell Area (um2)"}
  
  p <- ggplot(mean_df, aes(x = stamp, y = !!sym(var))) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(aes(color = sample), width = 0.25, size = 1, alpha = 0.7) +
    scale_color_manual(values = pal) +
    theme_bw() +
    labs(x = "Stamp", y = ylab) +
    geom_signif(comparisons = list(c("stamp_17", "stamp_18"), c("stamp_15", "stamp_16")),
                map_signif_level = TRUE, 
                test = "wilcox.test") +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          aspect.ratio = 1.5/1,
          legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    scale_x_discrete(labels = c("stamp_15" = "C1","stamp_16" = "C2",
                                "stamp_17" = "X1", "stamp_18" = "X2"))

if(var == "mean_sum") {p <- p + 
  scale_y_continuous(transform = "log10",
                     breaks = c(300,1000,2000),
                     expand = expansion(mult = c(0.05, 0.15)))}
if(var == "mean_detected") {p <- p + 
  scale_y_continuous(transform = "log10",
                     breaks = c(200,500,1000),
                     expand = expansion(mult = c(0.05, 0.15)))}
if(var == "mean_carea") {p <- p + 
  scale_y_continuous(breaks = c(100,200,300),
                     expand = expansion(mult = c(0.05, 0.15)))}

return(p)
}

qcmet <- wrap_plots(mean_bxp("mean_sum") + theme(legend.position = "none"),
           mean_bxp("mean_detected") + theme(legend.position = "none"),
           mean_bxp("mean_carea"),
           ncol =3) +
  plot_layout(axis_titles = "collect")

pdf(glue("{outdir}/qcmet.pdf"),width = 10,height = 5)
qcmet
dev.off()
################################################################################
# percentage of cells kept

pct_c <- ggplot(qcmet_df, aes(y = Var1, x = pct_kept, fill = stamp)) + 
  geom_col(position = "dodge", width = 0.8) +
  theme_bw() + 
  scale_fill_manual(values = pal, 
                    labels = c("stamp_15" = "C1",
                               "stamp_16" = "C2",
                               "stamp_17" = "X1",
                               "stamp_18" = "X2")) +
  geom_vline(xintercept = mean(qcmet_df$pct_kept),
             color = "black", linetype = "dashed", size =1) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 2/1,
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        plot.subtitle = element_text(hjust = 1, face = "bold"),
        plot.title = element_text(size = 10,color = "black", hjust = 0),
        legend.position = "none") +
  scale_x_continuous(breaks = c(0.00,0.50,1.00)) +
  labs(fill = "Replicate", x = "% cells kept", y = "Cell line IDs",
       subtitle = round(mean(qcmet_df$pct_kept),2)) 
pct_c
################################################################################
# number of cells kept
n_c <- ggplot(qcmet_df, aes( y = Var1, x = Freq_filt, fill = stamp)) +
  geom_col(position = "dodge", width = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = pal, 
                    labels = c("stamp_15" = "C1",
                               "stamp_16" = "C2",
                               "stamp_17" = "X1",
                               "stamp_18" = "X2")) +
  labs(x = "# Cells", y = "", fill = "Replicate") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 2/1,
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        plot.subtitle = element_text(hjust = 1, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y =  element_blank(),
        plot.title = element_text(size = 10,color = "black", hjust = 0),
        legend.position = "right") +
  scale_x_continuous(breaks = c(0,15000,30000),labels = c(0,15000,30000)) +
  geom_vline(xintercept = mean(qcmet_df$Freq_filt),
             color = "black", linetype = "dashed", size = 1) + 
  labs(title = glue("N = {sum(qcmet_df$Freq_filt)}"),
       subtitle = round(mean(qcmet_df$Freq_filt),0))
n_c

#### save
stamp <- "high_multi"
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/pct_nc.pdf"),width = 6)
wrap_plots(pct_c,n_c)
dev.off()
