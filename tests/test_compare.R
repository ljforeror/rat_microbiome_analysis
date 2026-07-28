# ============================================================
# tests/test_script.R
#
# Verifies that the code in this repository runs correctly:
# 1. Loads the example data (example_data/)
# 2. Recomputes alpha diversity (Observed, Chao1, Shannon, Simpson)
#    the same way as in the main analysis (via phyloseq::estimate_richness)
# 3. Compares the new result against tests/expected_output.csv
#
# Run from the root of the repository:
#   Rscript tests/test_script.R
# ============================================================

library(phyloseq)


#LOAD THE EXAMPLE DATA
# ------------------------------------------------------------

metadata_example <- read.csv("example_data/metadata_example.csv", row.names = "Sample")

# asv_counts_example: samples (rows) x ASVs (columns)
asv_counts_example <- as.matrix(read.csv("example_data/asv_counts_example.csv", row.names = 1))
asv_counts_example <- round(asv_counts_example)
storage.mode(asv_counts_example) <- "integer"

taxa_example <- as.matrix(read.csv("example_data/taxa_example.csv", row.names = 1))


# phyloseq OBJECT AND COMPUTE ALPHA DIVERSITY
#    taxa_are_rows = FALSE because in asv_counts_example
#    samples are in rows and ASVs are in columns
# ------------------------------------------------------------

ps_example <- phyloseq(
  otu_table(asv_counts_example, taxa_are_rows = FALSE),
  sample_data(metadata_example),
  tax_table(taxa_example)
)

new_result <- estimate_richness(
  ps_example,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)
new_result$Sample <- rownames(new_result)
new_result <- new_result[order(new_result$Sample), ]
rownames(new_result) <- NULL


#LOAD THE EXPECTED OUTPUT AND COMPARE
# ------------------------------------------------------------

expected_result <- read.csv("tests/expected_output.csv")
expected_result <- expected_result[order(expected_result$Sample), ]
rownames(expected_result) <- NULL

comparison <- all.equal(new_result, expected_result,
                        tolerance = 1e-6, check.attributes = FALSE)

if (isTRUE(comparison)) {
  cat("PASSED: the result matches tests/expected_output.csv\n")
} else {
  cat("FAILED: the result does not match the expected output.\n")
  cat("Difference details:\n")
  print(comparison)
  quit(status = 1)
}
##ENDDD