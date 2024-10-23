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
    axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1,color = "black"),
    axis.text.y = element_text(size = 12, color ="black"),
    axis.title = element_text(size = 18),
    legend.position = "right",
    aspect.ratio = 1.2/1
  ) +
  labs(x = "Cluster", y = "# Cells", fill = "Tech")

pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(df$Var2)

gg_prop <- ggplot(df, aes(x = Var1, y = perc, fill = Var2)) + 
  geom_col() + 
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        aspect.ratio = 2/1) +
  labs(x = "Tech", y = "Proportions", fill = "Cluster")
                                                     
pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/cnumbers.pdf", width = 6, height = 4)
cnumb
dev.off()

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/prop.pdf", width = 5, height = 6)
gg_prop
dev.off()