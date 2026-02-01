# Check for installed packages, if missing, install first
list.of.packages <- c("bigrquery", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
library(bigrquery)
library(dplyr)

# Authenticate with BigQuery
bq_auth(email = "mohd@variantbio.com")

# SQL query
sql <- '
SELECT
  id,
  approvedSymbol
FROM
  `open-targets-prod.platform.target`
'

# Run query
tb <- bq_project_query("starlit-vim-382123", sql)

# Download results to R
genes <- bq_table_download(tb)
