## SNR for CP and PhF
getSNR <- function(obj, markers = rownames(obj), data_type = 'counts') {
  # Prepare an empty list to store the results
  results <- list()

  # Loop over each marker to compute mean expressions
  for (marker in markers) {
    # Extract expression for the marker
    expression_values <- FetchData(obj, vars = marker, layer = data_type)[[marker]]

    # Identify top 20 highest-expressing cells
    top_20_values <- sort(expression_values, decreasing = TRUE)[1:20]
    top_20_mean <- mean(top_20_values)

    # Identify bottom 10% lowest-expressing cells
    bottom_10_percent_values <- sort(expression_values, decreasing = FALSE)[1:round(0.1 * length(expression_values))]
    bottom_10_percent_mean <- mean(bottom_10_percent_values)

    # Calculate the ratio (signal-to-background)
    ratio_top_bottom <- top_20_mean / ( bottom_10_percent_mean + 1 )

    # Store results
    results[[marker]] <- data.frame(
      Marker = marker,
      Top_20_Mean = top_20_mean,
      Bottom_10_Percent_Mean = bottom_10_percent_mean,
      Top_Bottom_Ratio = ratio_top_bottom
    )
  }

  # Combine results into a data frame
  results_df <- do.call(rbind, results)
}

