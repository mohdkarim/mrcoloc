library(openxlsx)
library(tidyverse)

wb <- loadWorkbook("output/mrcoloc_supplement.xlsx")
sheets <- getSheetNames("output/mrcoloc_supplement.xlsx")

cat("\n========================================\n")
cat("SUPPLEMENTARY TABLES SUMMARY REPORT\n")
cat("========================================\n\n")

for (sheet in sheets) {
  cat(sprintf("\n--- %s ---\n", sheet))
  
  df <- read.xlsx(wb, sheet)
  
  # Basic info
  cat(sprintf("  Dimensions: %d rows × %d columns\n", nrow(df), ncol(df)))
  cat(sprintf("  Columns: %s\n", paste(names(df)[1:min(5, ncol(df))], collapse=", ")))
  if (ncol(df) > 5) cat("          ... and", ncol(df)-5, "more\n")
  
  # Show first 2 rows as example
  cat("\n  First 2 rows (key columns only):\n")
  print(head(df[, 1:min(4, ncol(df))], 2))
  cat("\n")
}

cat("\n========================================\n")
cat("END OF REPORT\n")
cat("========================================\n")
