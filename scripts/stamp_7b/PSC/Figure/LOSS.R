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
qcmet <- qread(glue("{res_dir}/qcmet_df.qs"))


pal <- c("ectoderm" = "#8DF3F2", "iPSC_parental" = "#FE98FC", "mesoderm" = "#EBE682")

qcmet$pct_kept <-  round(qcmet$filt / qcmet$raw,2)

qcmet$loss <- 1 - qcmet$pct_kept
qcmet$Var1 <- factor(qcmet$Var1, 
                            levels = qcmet$Var1[order(-qcmet$loss)])


loss <- ggplot(qcmet, aes(x = Var1, y = loss, fill = Var1)) + 
  geom_col(position = "dodge", width = 0.8) +
  theme_bw() + 
  scale_fill_manual(
    values = pal
  ) +
  geom_hline(
    yintercept = mean(qcmet$loss),
    color = "black", linetype = "dashed", size = 1
  ) +
  # Add text annotation in the top-right corner inside the plot
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = round(mean(qcmet$loss),2),
    hjust = 1.2,  
    vjust = 1.2,
    size = 4,
    color = "black",
    fontface = "bold"
  ) +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 2.5/1,
    axis.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10,angle = 30, hjust = 1, vjust = 1),
    axis.title = element_text(size = 15, color = "black"),
    plot.subtitle = element_text(hjust = 0),  # Align subtitle to the left
    plot.title = element_text(size = 10, color = "black", hjust = 0),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0.00,0.20),breaks = c(0.00,0.10,0.20)) +
  labs(
    fill = "Sample",
    y = "LOSS",
    x = "Sub-stamp"
  )

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC_V2/loss.pdf", width = 3, height =4)
loss
dev.off()
