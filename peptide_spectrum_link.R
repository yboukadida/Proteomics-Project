# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("mzR", quietly = TRUE)) BiocManager::install("mzR")
library(mzR)

# Path to the mzML file
mzml_file <- "C:/Users/yesmi/OneDrive/Desktop/TMT.mzML"
msdata <- openMSfile(mzml_file)

# Extract MS2 metadata
hdr <- header(msdata)
ms2_scans <- hdr[hdr$msLevel == 2, ]

# Check if required data objects are loaded
if (!exists("res") || !exists("qnt")) {
  stop("Objects 'res' and/or 'qnt' not found. Please load them first.")
}

# Identify significant peptides
sig_peptides <- rownames(res[res$sig == "yes", ])
if (length(sig_peptides) == 0) {
  stop("No significant peptides found in 'res'.")
}

# Select the first significant peptide and its sequence
peptide_id <- sig_peptides[1]
sequence <- fData(qnt)[peptide_id, "sequence"]
cat("Significant peptide:", peptide_id, "\nSequence:", sequence, "\n")

# Search for a valid MS2 scan with peaks
valid_scan_found <- FALSE
for (i in 1:min(10, nrow(ms2_scans))) {
  scan_id <- ms2_scans$seqNum[i]
  spec <- peaks(msdata, scan_id)
  if (!is.null(spec) && nrow(spec) > 0) {
    valid_scan_found <- TRUE
    break
  }
}

if (!valid_scan_found) {
  stop("No valid MS2 scan found among the first 10 scans.")
}

# Plot and save the MS2 spectrum
png("MS2_spectrum_peptide.png", width = 800, height = 600)
plot(spec[,1], spec[,2], type = "h",
     main = paste("MS2 Spectrum - Peptide:", peptide_id, "\nScan", scan_id),
     xlab = "m/z", ylab = "Intensity", col = "darkblue")
dev.off()

# Close the mzML connection
close(msdata)

cat("Spectrum saved to: MS2_spectrum_peptide.png\n")

