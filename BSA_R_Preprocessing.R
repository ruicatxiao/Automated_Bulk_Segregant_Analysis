#!/usr/bin/env Rscript

# Load necessary library
library(ggplot2)
library(readr)

# Define the format_genomic function
format_genomic <- function(...) {
  function(x) {
    limits <- c(1e0, 1e3, 1e6)
    i <- findInterval(abs(x), limits)
    i <- ifelse(i == 0, which(limits == 1e0), i)
    paste(
      format(round(x / limits[i], 1), trim = TRUE, scientific = FALSE, ...)
    )
  }
}

# Read input file from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input file provided.")
}
input_file <- args[1]

# Get the base name of the input file without extension
input_file_base <- sub("\\.[^.]*$", "", basename(input_file))
input_file_base <- sub("[.]$", "", input_file_base)

# Read in the data
BSA <- read.table(input_file, sep = '\t', header = TRUE)
BSAdim <- dim(BSA)

# Extract AD and DP columns
AD <- BSA[, seq(5, BSAdim[2], 4)]
DP <- BSA[, seq(6, BSAdim[2], 4)]

# Function to extract reference allele depth
ref.DP <- function(X) {
  as.numeric(strsplit(as.character(X), ",")[[1]])[1]
}

ADdim <- dim(AD)
refFre.AD <- matrix(ncol = ADdim[2], nrow = ADdim[1])

# Calculate reference allele frequency
for (j in 1:ADdim[2]) {
  refFre.AD[, j] <- sapply(AD[, j], ref.DP, simplify = "array") / as.numeric(DP[, j])
}

# Filter out low DP values
refFre.AD[DP < 10] <- NA
colnames(refFre.AD) <- colnames(AD)

# Combine with original data
refFre.AD <- cbind(BSA[, 1:4], refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

# Display head of the data and unique chromosomes
head(refFre.AD)
unique(refFre.AD$CHROM)

# Optional: Remove specific chromosomes (uncomment if needed)
# refFre.AD <- refFre.AD[!(refFre.AD$CHROM %in% c("Ctyz_00_1", "Ctyz_00_2", "Ctyz_00_3")), ]

RefFreADdim <- dim(refFre.AD)

# Write the processed data to a TSV file if needed
write.table(
  cbind(refFre.AD, DP[1:RefFreADdim[1], ]),
  file = paste0(input_file_base, ".BSA.refFre.AD.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Remove rows with missing values
refFre.AD <- refFre.AD[complete.cases(refFre.AD), ]

# Check if there's any data left after filtering
if (nrow(refFre.AD) == 0) {
  cat("Warning: No data points left after filtering for plotting.\n")
  quit(status = 0) # Exit gracefully
}

# Generate plots for each sample column starting from column 5
for (i in 5:ncol(refFre.AD)) {
  colname_i <- colnames(refFre.AD)[i]
  cat("Processing column:", colname_i, "\n") # Print current column being processed

  # Check if the current sample column has any non-NA values
  sample_data <- refFre.AD[[colname_i]]
  if (all(is.na(sample_data))) {
    cat("  Skipping", colname_i, "as all values are NA.\n")
    next # Move to the next sample
  }

  # Subset data to remove NA values for the current sample
  plot_data <- refFre.AD[!is.na(sample_data), ]

  # Check if there's sufficient data for plotting after removing NAs
  if (nrow(plot_data) == 0) {
    cat("  Skipping", colname_i, "as no non-NA data points remain after filtering.\n")
    next
  }

  # Check if there are enough unique POS values for the smoothing polynomial
  unique_pos_count <- length(unique(plot_data$POS))
  required_points_for_poly <- 8 + 1 # degree 8 polynomial needs at least 9 unique points
  if (unique_pos_count < required_points_for_poly) {
    cat("  Warning:", colname_i, "has only", unique_pos_count, "unique POS values, insufficient for 8th degree polynomial smoothing.\n")
    # Optionally, skip the geom_smooth part or use a lower degree polynomial
    # For now, we'll proceed without geom_smooth if there are insufficient points
    use_smooth <- FALSE
  } else {
    use_smooth <- TRUE
  }

  # Prepare the base plot
  p <- ggplot(data = plot_data) +
    facet_grid(~CHROM, scales = "free_x", space = "free_x") +
    ylim(0, 1) +
    ggtitle(paste(colname_i, "raw allele frequency plot")) +
    geom_hline(yintercept=0.5, linetype="dashed",  color = "black") +
    geom_point(aes_string(x = "POS", y = colname_i), color = "#3933ff", size = 0.5, alpha = 0.9)

  # Add smoothing line only if conditions are met
  if (use_smooth) {
    p <- p + geom_smooth(aes_string(x = "POS", y = colname_i), method = "lm", formula = y ~ poly(x,8), se = TRUE, color = "red")
  } else {
    cat("  Skipping polynomial smoothing for", colname_i, "due to insufficient data points.\n")
  }

  # Try to set x-axis breaks and labels, handle potential errors
  tryCatch({
    # Calculate breaks - ensure max(POS) > 0 to avoid log10(0) or log10(negative)
    max_pos <- max(plot_data$POS)
    if (max_pos > 0) {
      calculated_breaks <- seq(
        from = 0,
        to = max_pos,
        by = 10^(floor(log10(max_pos)))
      )
      # Only apply breaks if the sequence is not empty or trivially short
      if (length(calculated_breaks) > 1) {
        p <- p + scale_x_continuous(
          breaks = calculated_breaks,
          labels = format_genomic()
        )
      } else {
        # If calculated breaks are too few, let ggplot handle it automatically
        cat("  Using default x-axis breaks for", colname_i, "due to low max POS.\n")
      }
    } else {
      # If max POS is 0 or negative, let ggplot handle it automatically
      cat("  Using default x-axis breaks for", colname_i, "due to non-positive max POS.\n")
    }
  }, error = function(e) {
    # If scale_x_continuous fails for any reason, warn and proceed with defaults
    cat("  Warning: Could not set custom x-axis breaks for", colname_i, ":", conditionMessage(e), "\n")
    # ggplot will use defaults if scale_x_continuous is not added or fails within tryCatch
  })

  # Save each plot to a PDF file
  plot_filename <- paste0(input_file_base, ".", colname_i, "raw.pdf")
  cat("  Saving plot to:", plot_filename, "\n")
  ggsave(filename = plot_filename, plot = p, width=30, height=6)

  # Optional: Print a message if the plot was saved successfully
  # (ggsave usually doesn't print an error message if it fails within Rscript, but the tryCatch above helps catch scale_x_continuous issues)
}

cat("R script completed successfully.\n")