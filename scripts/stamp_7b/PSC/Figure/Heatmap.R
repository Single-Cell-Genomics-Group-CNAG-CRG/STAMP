library(dplyr)
library(patchwork)
library(ggplot2)
library(qs)
library(here)
library(glue)

dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_7b"
sample <- "iPSCs"

dir <- glue("{proj_dir}/data/{stamp}/{sample}/Ist")
unsup <- qread(glue("{dir}/unsup_3_10.qs"))


## InSituType Profiles
#Normalize & log-transform the profiles
norm <- unsup$profiles /  Matrix::colSums(unsup$profiles) # average
norm <- norm + 1e-6

#Calculate logFC
result <- norm 
for (i in 1:nrow(norm)) {
  for (j in 1:ncol(norm)) {
    result[i, j] <- log2(norm[i, j] / mean(norm[i, -j]))
  }
}

feat <- unique(as.vector(unlist(apply(result, 2, function(column) {
  names(sort(column, decreasing = TRUE))[1:5]
}))))

#Remove outliers at 2.5 sd for better visualization
x <- result[feat,]
z <- \(x, th=2.5) {
  if (is.null(dim(x))) {
    x[x < 0] <- 0
    sd <- sd(x, na.rm=TRUE)
    x <- x-mean(x, na.rm=TRUE)
    if (sd != 0) x <- x/sd
  } else {
    mus <- colMeans(x)
    sds <- colSds(x)
    x <- sweep(x, 2, mus, `-`)
    x <- sweep(x, 2, sds, `/`)
  }
  x[x > +th] <- +th
  x[x < -th] <- -th
  return(x)
}
mtx_scaled <- z(x)

#Plot
library(reshape2)
library(ggdendro)

# Perform hierarchical clustering on columns
column_dist_euclidean <- dist(t(mtx_scaled), method = "euclidean")  # Use transpose to cluster columns
column_hclust_complete <- hclust(column_dist_euclidean, method = "average")

# Reorder the columns based on clustering
ordered_columns <- column_hclust_complete$order
mtx_scaled_ordered <- mtx_scaled[, ordered_columns]
colnames(mtx_scaled_ordered) <- colnames(mtx_scaled)[ordered_columns]

# Change annotation
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "a"] <- "amnion-precursors"
#colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "a"] <- "mesoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "b"] <- "pluripotent"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "c"] <- "pluripotent"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "d"] <- "BMP-induced Prog."
#colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "d"] <- "mesoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "e"] <- "endoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "f"] <- "mesoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "g"] <- "late meso."
#colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "g"] <- "mesoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "h"] <- "endoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "i"] <- "endoderm"
colnames(mtx_scaled_ordered)[colnames(mtx_scaled_ordered) == "j"] <- "ectoderm"

# Convert the matrix into a format suitable for ggplot
mtx_melted <- melt(t(mtx_scaled_ordered))

# Reorder annotation factors
mtx_melted$Var1 <- factor(mtx_melted$Var1, levels = rev(c("pluripotent","endoderm",
                                                      "ectoderm",
                                                      "mesoderm")))

# Reorder genes factors
var2 <- c("IGFBP7", "APOA1", "PPARG", "TTR", "CACNA1C", "INHBB", "NPR1", 
                  "CMKLR1", "IL1A", "CD69", "FAP", "TNFSF8", "ITK", "CD40LG", 
                  "C5AR2", "CYP1B1", "TACSTD2", "RAMP1", "PLAC8", "LTBR", "PDGFB", 
                  "DHRS2", "SPOCK2", "MST1R", "ACTG2", "WNT9A", "LAMA4", "COL5A3", 
                  "LY75", "LDHA", "LGR5", "FOXF1", "LUM", "KLF2", "ITGA8", 
                  "HLA-DPB1", "TFEB", "INS", "KLRF1", "CASR", "CCL3/L1/L3", 
                  "LIF", "FGR", "IL34", "IL18R1", "NRG1", "PTGS1", "EPHA4", 
                  "BMP3", "NR2F2")

mtx_melted$Var2 <- factor(mtx_melted$Var2, 
                          levels = c("INHBB", "NPR1","CMKLR1", "IL1A", "CD69", "FAP",
                                     "TNFSF8", "ITK", "CD40LG","C5AR2",
                                     "CCL3/L1/L3","LIF", "FGR", "IL34", "IL18R1",
                                     "NRG1", "PTGS1", "EPHA4", "BMP3", "NR2F2",
                                     "WNT9A", "LAMA4", "COL5A3","LY75", "LDHA",
                                     "IGFBP7", "APOA1", "PPARG", "TTR", "CACNA1C",
                                     "PDGFB","DHRS2", "SPOCK2", "MST1R", "ACTG2",
                                     "HLA-DPB1", "TFEB", "INS", "KLRF1", "CASR"))
mtx_melted <- na.omit(mtx_melted)
# Create the heatmap using ggplot2
hm <- ggplot(mtx_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, name = "expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    aspect.ratio = 1/3
  ) +
  coord_flip() +
  labs(x = "", y = "")
hm

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/hm.pdf", width = 10)
hm
dev.off()
