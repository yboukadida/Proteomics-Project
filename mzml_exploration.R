# 05_ms2_spectrum_explorer.R

# Load mzR
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("mzR", quietly = TRUE)) BiocManager::install("mzR")
library(mzR)

# Local .mzML file path
mzml_file <- "C:/Users/yesmi/OneDrive/Desktop/TMT.mzML"

# Open file with error handling
msdata <- tryCatch({
  openMSfile(mzml_file)
}, error = function(e) {
  stop("Failed to open mzML file: ", e$message,
       "\nCheck the path and file existence.")
})

# Get scan metadata
hdr <- header(msdata)
print(head(hdr))

# Plot Total Ion Current (TIC) chromatogram
plot(hdr$retentionTime, hdr$totIonCurrent, type = "l", lwd = 2,
     col = "blue", main = "TIC Chromatogram",
     xlab = "Retention Time (s)", ylab = "Total Ion Current")

# Get MS2 scans
ms2_scans <- hdr[hdr$msLevel == 2, ]

# Select scan to plot (10th or first)
if (nrow(ms2_scans) < 10) {
  warning("Less than 10 MS2 scans found. Using the first one.")
  scan_id <- ms2_scans$seqNum[1]
} else {
  scan_id <- ms2_scans$seqNum[10]
}

# Get MS2 peaks
spec <- peaks(msdata, scan_id)

# Plot MS2 spectrum
plot(spec[, 1], spec[, 2], type = "h", col = "darkred",
     main = paste("MS2 Spectrum – Scan", scan_id),
     xlab = "m/z", ylab = "Intensity")

# Close connection
close(msdata)

