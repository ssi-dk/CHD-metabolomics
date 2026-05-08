# ---------------------------------------------------------------------------
# Creates an anonymized feature table for submission to GNPS by:
#   1. Replacing sample IDs with generic names (Sample_N.mzML)
#   2. Zeroing out all peak area values
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Paths — edit these before running
# ---------------------------------------------------------------------------

input_path  <- "path/to/input.csv"
output_path <- "path/to/output.csv"

# ---------------------------------------------------------------------------
# Column names that are metadata, not sample names
# ---------------------------------------------------------------------------

cids <- c(
  "row ID", "row m/z", "row retention time",
  "correlation group ID", "annotation network number",
  "best ion", "auto MS2 verify", "identified by n=",
  "partners", "neutral M mass"
)

# ---------------------------------------------------------------------------
# Load feature table
# ---------------------------------------------------------------------------

ft <- read.table(opt$input, sep = ",", header = TRUE, check.names = FALSE)
cat(sprintf("Loaded feature table: %d rows x %d columns\n", nrow(ft), ncol(ft)))

# ---------------------------------------------------------------------------
# Remove columns that are entirely NA (excluding metadata columns)
# ---------------------------------------------------------------------------

all_na_cols <- which((colSums(is.na(ft))==nrow(ft)) == T)
all_na_cols <- all_na_cols[-which(names(all_na_cols) %in% cids)]

if (length(all_na_cols) > 0) {
  ft <- ft[, -all_na_cols]
  cat(sprintf("Removed %d all-NA columns\n", length(all_na_cols)))
}

# ---------------------------------------------------------------------------
# Replace sample column names with generic identifiers
# ---------------------------------------------------------------------------

idxs <- which(colnames(ft)%in% cids)[length(cids)]+1
colnames(ft)[idxs:ncol(ft)] <- paste(paste('Sample',idxs:ncol(ft), sep = '', '.mzML'), ' filtered Peak area', sep = '')

# ---------------------------------------------------------------------------
# Zero out all peak area values
# ---------------------------------------------------------------------------

ft[,idxs:ncol(ft)] <- 0

# ---------------------------------------------------------------------------
# Save anonymized feature table
# ---------------------------------------------------------------------------

write.csv(ft, output_path, quote = FALSE, row.names = FALSE)
cat(sprintf("Anonymized feature table saved to: %s\n", output_path))
cat(sprintf("Final dimensions: %d rows x %d columns\n", nrow(ft), ncol(ft)))