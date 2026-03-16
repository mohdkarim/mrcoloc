# Open Targets evidence triangulation
# Loads from local cache (data/ot_triangulation.rds) if available.
# Falls back to BigQuery if cache is missing and bigrquery is installed.
# Cache is generated on first successful run; delete to force refresh.
#
# Expects `ensid` and `efo` vectors in calling environment.

library(dplyr)

cache_path <- file.path(
  Sys.getenv("MRCOLOC_ROOT", getwd()),
  "data", "ot_triangulation.rds"
)

if (file.exists(cache_path)) {
  tr <- readRDS(cache_path)
  message("   Loaded ", nrow(tr), " triangulation rows from cache: ", cache_path)
} else {
  message("   Cache not found, querying BigQuery...")
  if (!requireNamespace("bigrquery", quietly = TRUE)) {
    stop("Cache file ", cache_path, " not found and bigrquery is not installed.\n",
         "Run the supplement pipeline once with BigQuery access to generate the cache.")
  }
  library(bigrquery)
  bq_project <- Sys.getenv("BQ_PROJECT_ID", "")
  bq_email <- Sys.getenv("BQ_EMAIL", "")
  if (bq_project == "" || bq_email == "") {
    stop("Set BQ_PROJECT_ID and BQ_EMAIL environment variables for BigQuery access.\n",
         "Or run the supplement pipeline once with cached data.")
  }
  bq_auth(email = bq_email)

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
  tb <- bq_project_query(bq_project, sql)
  tr <- bq_table_download(tb)
  tr$key <- with(tr, paste0(targetId, "_", diseaseId))

  # Save cache for future runs
  saveRDS(tr, cache_path)
  message("   Queried ", nrow(tr), " rows and saved cache to: ", cache_path)
}
