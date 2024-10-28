# Libraries
library(SingleCellExperiment)
library(here)
library(glue)
library(data.table)
library(qs)
library(DropletUtils)
# paths
dir <- glue("{here()}")
# Parameters and paths
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_18"
data_dir <- glue("{proj_dir}/data/{stamp}/raw")

# Define the samples to process
samples <- c("A", "B", "C", "D", "E", "F", "G")

# Initialize a list to store SingleCellExperiment objects
sce_list <- list()

# Loop over each sample
for (sample in samples) {
#  # Construct file paths using glue
  h5_file <- glue("{data_dir}/{sample}.h5")
  csv_file <- glue("{data_dir}/MD_{sample}.csv")
  
  # Read the data
  sce <- read10xCounts(h5_file)
  rownames(sce) <- rowData(sce)$Symbol
  rowd <- rowData(sce)
  cd <- fread(csv_file, sep =",", header = T)
  
  # Split matrix
  counts_mat <- counts(sce)
  gs <- rownames(counts_mat)
  ncp <- grep("NegControlProbe", gs)
  ncc <- grep("NegControlCodeword", gs)
  uc <- grep("UnassignedCodeword", gs)
  dep <- grep("Deprecated", gs)
  int <- grep("Intergenic", gs)
  
  # Alternative experiment
  as <- list(counts = counts_mat[-c(ncp, ncc, uc, dep, int), ])
  ae <- list(
    NegControlProbe = SingleCellExperiment(list(counts = counts_mat[ncp, ])),
    NegControlCodeword = SingleCellExperiment(list(counts = counts_mat[ncc, ])),
    DeprecatedCodeword = SingleCellExperiment(list(counts = counts_mat[dep, ])),
    IntergenicRegion = SingleCellExperiment(list(counts = counts_mat[int, ])),
    UnassignedCodeword = SingleCellExperiment(list(counts = counts_mat[uc, ]))
  )
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(as, colData = cd, altExps = ae)
  sce$region <- sample
  
  # filter rowData
  rowd <- rowd[rowd$Symbol %in% rownames(sce),]
  rowData(sce) <- rowd # add rowdata back to sce
  # Store the SCE object in the list
  sce_list[[sample]] <- sce
  
  # Optional: Plotting for each sample
  df <- as.data.frame(colData(sce))
  p <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
    geom_point(shape = 16, size = 0.01) +
    coord_equal() +
    ggtitle(glue("Sample {sample}")) +
    theme_minimal()
  
  # Display the plot
  print(p)
}

sce <- do.call(cbind,sce_list)

df <- as.data.frame(colData(sce))

ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point(shape = 16, size = 0.01) +
  coord_equal() +
  ggtitle(glue("Sample {sample}")) +
  theme_minimal() +
  facet_wrap(~region)

# save
res_dir <- glue("{proj_dir}/data/{stamp}/processed")
dir.create(res_dir, showWarnings = F)
qsave(sce, file = glue("{res_dir}/raw_sce.qs"))
