# Download association_by_datasource_indirect from OT
# Check for installed packages, if missing, install first
list.of.packages <- c("bigrquery", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load libraries
library(bigrquery)
library(dplyr)
bq_auth(email = "mohd@variantbio.com")
# SQL
ensid2 <- toString(sprintf("'%s'", ensid))
efo2 <- toString(sprintf("'%s'", efo))
sql <- paste0('SELECT
  targetId,
  diseaseId,
  STRING_AGG(DISTINCT datatypeId, ", ") AS concatenatedDatatypeIds,
  STRING_AGG(DISTINCT datasourceId, ", ") AS concatenatedDatasourceIds,
  MAX(CASE WHEN datatypeId = "genetic_association" THEN score END) AS harmonic_genetic_score
FROM
  `open-targets-prod.platform.association_by_datasource_indirect`
WHERE
  targetId IN (', ensid2, ')
  AND diseaseId IN (', efo2, ')
GROUP BY
  targetId,
  diseaseId
')
tb <- bq_project_query("starlit-vim-382123", sql)
tr <- bq_table_download(tb)
tr$key <- with(tr, paste0(targetId, "_", diseaseId))