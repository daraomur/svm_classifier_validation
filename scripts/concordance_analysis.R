# =============================================================================
# Sinonasal Methylation Classifier Validation — Concordance Analysis
# SPELCASTER WP1b | RCSI
# Author: Dara O'Murchu
# Date: June 2026
# =============================================================================

# --- 0. Libraries ------------------------------------------------------------
library(dplyr)
library(readxl)
library(stringr)

# --- 1. Paths ----------------------------------------------------------

base_dir <- "~/Desktop/Research Project/My Project/Dara_Validation"

clf_path <- file.path(base_dir, "results/classifier/classifier_results_with_qc.csv")
cohort_path <- file.path(base_dir, "results/pathology_concordance/pathology_diagnoses_cohort.xlsx")
out_results <- file.path(base_dir, "results/pathology_concordance")

dir.create(out_results, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load classifier results -----------------------------------------

clf <- read.csv(clf_path, stringsAsFactors = FALSE) %>%
  mutate(
    sample_id = str_trim(sample_id),
    sample_id = str_remove(sample_id, "\\s+T$")  # strip stray " T" suffix
  )

# --- 3. Load curated pathology cohort (with Matched Label) --------------

cohort <- read_excel(cohort_path) %>%
  setNames(c("sample_id", "tumour_type", "pathology_note", "matched_label")) %>%
  mutate(
    sample_id = str_trim(sample_id),
    tumour_type = str_trim(tumour_type),
    matched_label = str_trim(matched_label)
  )

# --- 4. Reference-class status lookup (for sens/spec eligibility only) --
# Does NOT drive the concordance label (Matched Label already does that).
# Flags why a tumour_type lacks a 1:1 classifier-class match:
#   "No Ref"    - entity not represented in the 18-class taxonomy
#   "Ambiguous" - conventional diagnosis covers multiple possible classes
#                 (e.g. old "SNUC" wastebasket diagnosis)

no_ref_entities <- c(
  "Schneiderian Papilloma",
  "Papilloma",
  "Polymorphous ADCA",
  "Chondrosarcoma",
  "Biphenotypic Sinonasal Sarcoma",
  "Nasopharyngheal Carcinoma"
)

ambiguous_candidates <- list(
  "SNUC" = "SMARCB1|NEC-like IDH2|NEC-like SMARCA4 ARID1A|ACC",
  "SNUC Sinonasal Undifferentiated Carcinoma" = "SMARCB1|NEC-like IDH2|NEC-like SMARCA4 ARID1A|ACC",
  "LCC" = "SMARCB1|NEC-like IDH2|NEC-like SMARCA4 ARID1A|ACC",
  "Poorly differentiated high grade neuroendocrine CA" = "NEC-like IDH2|NEC-like SMARCA4 ARID1A|SMARCB1",
  "Sinonasal Neuroendocrine Carcinoma" = "NEC-like IDH2|NEC-like SMARCA4 ARID1A|SMARCB1",
  "Round Blue Cell Neoplasm" = "ONB|ALV RMS|EMB RMS|EWS|MCC",
  "Myoepithelial carcinoma" = "SMARCB1"
)

lookup_notes <- c(
  "Schneiderian Papilloma" = "Benign; no malignant-class equivalent",
  "Papilloma" = "Benign; no malignant-class equivalent",
  "Polymorphous ADCA" = "Salivary-gland-type lineage, distinct from ACC/ADC",
  "Chondrosarcoma" = "Mesenchymal, not represented in taxonomy",
  "Biphenotypic Sinonasal Sarcoma" = "Mesenchymal (PAX3-fusion), not represented in taxonomy",
  "Nasopharyngheal Carcinoma" = "Diagnosis TBC - site/differential in report inconsistent with NPC label",
  "SNUC" = "Jurmeister SNUC reclassification group",
  "SNUC Sinonasal Undifferentiated Carcinoma" = "Jurmeister SNUC reclassification group",
  "LCC" = "Poorly-differentiated wastebasket, same group as SNUC",
  "Poorly differentiated high grade neuroendocrine CA" = "NEC differential group",
  "Sinonasal Neuroendocrine Carcinoma" = "NEC differential group",
  "Round Blue Cell Neoplasm" = "Morphology-only differential, no further IHC/molecular detail in report",
  "Myoepithelial carcinoma" = "SMARCB1-deficient carcinoma frequently misdiagnosed as myoepithelial CA on H&E",
  "Well differentiated Adenocarcinoma" = "Secondary diagnosis - original biopsy inaccessible"
)

get_status <- function(tumour_type) {
  if (tumour_type %in% no_ref_entities) return("No Ref")
  if (tumour_type %in% names(ambiguous_candidates)) return("Ambiguous")
  return("Direct")
}

# --- 5. Merge classifier output with pathology cohort --------------------

merged <- clf %>%
  left_join(cohort, by = "sample_id") %>%
  rowwise() %>%
  mutate(
    status = get_status(tumour_type),
    lookup_note = ifelse(tumour_type %in% names(lookup_notes), lookup_notes[[tumour_type]], NA)
  ) %>%
  ungroup()

unmatched <- merged %>% filter(is.na(matched_label))
if (nrow(unmatched) > 0) {
  warning(
    nrow(unmatched), " sample(s) have no Matched Label - check sample_id join: ",
    paste(unmatched$sample_id, collapse = "; ")
  )
}

# --- 6. Derive 3-category concordance label -------------------------------
# Per Jurmeister et al. outcome framework, there are only three outcomes.
# "Confirmed" = predicted_class exactly matches the curated Matched Label.
# Everything else with a confident call is "Reclassified" - true even where
# the conventional diagnosis has no reference class of its own (e.g.
# Papilloma -> SCC, Polymorphous ADCA -> ACC), per Jurmeister's own usage
# of "reclassified" for exactly this scenario.

merged <- merged %>%
  mutate(
    concordance = case_when(
      predicted_class == "Unknown" ~ "Non-classifiable",
      predicted_class == matched_label ~ "Diagnosis Confirmed",
      TRUE ~ "Diagnosis Reclassified"
    )
  )

# --- 7. Sens/spec eligibility flag for downstream per-class stats --------

merged <- merged %>%
  mutate(
    sens_spec_eligible = case_when(
      status == "Direct" ~ "True class - include",
      status == "No Ref" & predicted_class != "Unknown" ~ "No ref class, but confident call - review as possible FP",
      status == "No Ref" ~ "No reference class - exclude",
      status == "Ambiguous" ~ "Ambiguous differential - exclude from per-class stats",
      TRUE ~ NA_character_
    )
  )

# --- 8. Output table -------------------------------------------------------

results <- merged %>%
  select(
    sample_id, tumour_type, matched_label, predicted_class, confidence_score,
    confidence_tier, qc_flag, concordance, sens_spec_eligible,
    lookup_note, pathology_note
  ) %>%
  arrange(sample_id)

cat("\n=== Concordance summary ===\n")
print(table(results$concordance))

cat("\n=== Sens/spec eligibility summary ===\n")
print(table(results$sens_spec_eligible))

write.csv(results, file.path(out_results, "concordance_results.csv"), row.names = FALSE)

cat("\nSaved to:", file.path(out_results, "concordance_results.csv"), "\n")
