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

# Write the processed data to a CSV file if needed
# write.csv(
#  cbind(refFre.AD, DP[1:RefFreADdim[1], ]),
#  file = paste0(input_file_base, ".BSA.refFre.AD.csv"),
#  row.names = FALSE
# )


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

# Generate plots for each sample column starting from column 5, which is the start of bulk sample data
for (i in 5:ncol(refFre.AD)) {
  colname_i <- colnames(refFre.AD)[i]
  p <- ggplot(data = refFre.AD) +
    scale_x_continuous(
      breaks = seq(
        from = 0,
        to = max(refFre.AD$POS),
        by = 10^(floor(log10(max(refFre.AD$POS))))
      ),
      labels = format_genomic()
    ) +
    facet_grid(~CHROM, scales = "free_x", space = "free_x") +
    ylim(0, 1) +
    ggtitle(paste(colname_i, "raw allele frequency plot")) +
    geom_hline(yintercept=0.5, linetype="dashed",  color = "black") +
    geom_point(aes_string(x = "POS", y = colname_i), color = "#3933ff", size = 0.5, alpha = 0.9) +
    geom_smooth(aes_string(x = "POS", y = colname_i), method = "lm", formula = y ~ poly(x,8), se = TRUE, color = "red")
  
  # Save each plot to a PDF file
  ggsave(filename = paste0(input_file_base, ".", colname_i, "raw.pdf"), plot = p, width=30, height=6)
}
