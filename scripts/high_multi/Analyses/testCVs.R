# Load necessary libraries
library(ggplot2)
library(dplyr)

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

# Set parameters
set.seed(42)
n_sub_samples <- 21      
n_datasets <- 4         
n_measurements <- 1000   
desired_cv <- 5       

# Generate sub-sample labels (A, B, C, ...)
sub_sample_labels <- LETTERS[1:n_sub_samples]

# Initialize an empty data frame to store the simulated data
data <- data.frame()

# Simulate data
for (dataset_num in 1:n_datasets) {
  dataset_name <- paste0("stamp-",dataset_num)
  for (sub_sample in sub_sample_labels) {
    # Set mean for each dataset (can vary if needed)
    mean_value <- 10  # Keeping mean constant for simplicity
    # Calculate standard deviation to achieve the desired CV
    sd_value <- (desired_cv / 100) * mean_value
    # Generate measurements for each sub-sample in each dataset
    values <- rnorm(n_measurements, mean = mean_value, sd = sd_value)
    temp_df <- data.frame(
      sub_stamp = sub_sample,
      stamp = dataset_name,
      nCount = values
    )
    data <- rbind(data, temp_df)
  }
}

# Calculate CV for each sub-sample in each dataset
cv_data <- data %>%
  group_by(sub_stamp, stamp) %>%
  summarise(
    mean_value = mean(nCount),
    sd_value = sd(nCount),
    cv = (sd_value / mean_value) * 100
  )

names(pal) <- unique(cv_data$sub_stamp)

# Plot CV for each sub-sample across datasets
ggplot(cv_data, aes(x = stamp, y = cv, color = sub_stamp, group = sub_stamp)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Stamp",
    y = "Coefficient of Variation (%)",
    title = "Coefficient of Variation Across Stamps for Each Sub-stamp",
    color = "Sub-stamp"
  ) +
  theme_bw() +
  scale_color_manual(values = pal)
  
