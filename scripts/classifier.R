# =============================================================================
# Sinonasal Methylation Classifier Validation — Jurmeister et al. SNT Application
# SPELCASTER WP1b | RCSI
# Author: Dara O'Murchu
# Date: May 2026
# =============================================================================

# --- 0. Libraries ------------------------------------------------------------

library(e1071)
library(glmnet)
library(dplyr)
library(tibble)

# --- 1. Paths ----------------------------------------------------------------

base_dir <- "~/Desktop/Research Project/My Project/Dara_Validation"
classifier_dir <- file.path(base_dir, "data/collaborator_files")

out_results <- file.path(base_dir, "results/classifier")
dir.create(out_results, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load inputs ----------------------------------------------------------

betas <- readRDS(file.path(base_dir, "results/beta_matrices/batch01_betas_QCDPB_lifted.rds"))
cat("Beta matrix loaded:", nrow(betas), "probes x", ncol(betas), "samples\n")

qc <- read.csv(file.path(base_dir, "results/qc/batch01_qc_summary.csv"))

# --- 3. Load classifier models -----------------------------------------------

pred_model <- readRDS(file.path(classifier_dir, "SVM_model.rds"))
cal_model <- readRDS(file.path(classifier_dir, "SVM_calibration_model.rds"))

features <- scan(file.path(classifier_dir, "sel_CpGs_20000.txt"),
                 what = character(), quiet = TRUE)
cat("Classifier feature set:", length(features), "CpGs\n")

# --- 4. Probe alignment ------------------------------------------------------

probes_in_betas <- rownames(betas)
probes_found <- intersect(features, probes_in_betas)
probes_missing <- setdiff(features, probes_in_betas)

cat("Classifier CpGs found in beta matrix:", length(probes_found), "\n")
cat("Classifier CpGs missing:", length(probes_missing), "\n")
cat("Coverage:", round(length(probes_found) / length(features) * 100, 1), "%\n")

# --- 5. Feature subsetting and imputation ------------------------------------

n_samples <- ncol(betas)
sample_ids <- colnames(betas)

# Initialise feature matrix with NA
model_feature_betas <- matrix(
  NA_real_,
  nrow = length(features),
  ncol = n_samples,
  dimnames = list(features, sample_ids)
)

# Fill in features present in the beta matrix
model_feature_betas[probes_found, ] <- betas[probes_found, ]

cat("NAs before imputation:", sum(is.na(model_feature_betas)), "\n")

frac_imputed <- colMeans(is.na(model_feature_betas))

# Impute missing values per sample with sample mean (if-else: assign 0.5)
for (j in seq_len(ncol(model_feature_betas))) {
  miss <- is.na(model_feature_betas[, j])
  if (!any(miss)) next
  fill <- mean(model_feature_betas[, j], na.rm = TRUE)
  if (!is.finite(fill)) fill <- 0.5
  model_feature_betas[miss, j] <- fill
}

cat("NAs remaining after imputation:", sum(is.na(model_feature_betas)), "\n")

# Transpose: classifier expects samples x features
x <- t(model_feature_betas)

# --- 6. Helper functions -----------------------------------------------------

predict_prob_SVM <- function(model, x) {
  pred_obj <- predict(model, x, probability = TRUE)
  probs <- attr(pred_obj, "probabilities")
  probs[, attr(pred_obj, "levels"), drop = FALSE]
}

predict_cal_probs <- function(cal_model, raw_probs) {
  probs <- predict(cal_model, raw_probs, type = "response", s = "lambda.min")[, , 1]
  if (is.vector(probs)) {
    probs <- matrix(probs, nrow = 1,
                    dimnames = list(rownames(raw_probs)[1], names(probs)))
  }
  probs
}

normalize_class_names <- function(pred) {
  pred <- ifelse(pred == "APA", "PIT AD", pred)
  pred <- ifelse(pred == "NEC IDH2", "NEC-like IDH2", pred)
  pred <- ifelse(pred == "NEC SWI SNF", "NEC-like SMARCA4 ARID1A", pred)
  pred <- ifelse(pred == "unknown", "Unknown", pred)
  pred
}

apply_target_class_names <- function(prob_mat) {
  target <- c("SMARCB1", "NEC-like SMARCA4 ARID1A", "ALV RMS", "ACC", "ADC", 
              "PIT AD", "CPH", "EWS", "ONB", "SCC", "NEC-like IDH2", "EMB RMS", 
              "NUT", "MCC", "LECA", "MELA", "GPC", "CTRL", "Unknown")
  if (ncol(prob_mat) == length(target)) colnames(prob_mat) <- target
  prob_mat
}

compute_svm_feature_screen_scores <- function(pred_model, feature_names) {
  svm_model <- pred_model$model
  if (is.null(svm_model$SV) || is.null(svm_model$coefs))
    stop("SVM model does not expose SV/coefs for feature screening.")
  sv_mat <- svm_model$SV
  sv_names <- colnames(sv_mat)
  if (is.null(sv_names)) stop("Support vector matrix has no feature names.")
  sv_weights <- rowSums(abs(svm_model$coefs))
  if (!all(is.finite(sv_weights)) || sum(sv_weights) <= 0)
    stop("Invalid support-vector weights.")
  sw <- sv_weights / sum(sv_weights)
  sv_mean <- as.numeric(crossprod(sw, sv_mat))
  centered <- sweep(sv_mat, 2, sv_mean, "-")
  sv_var <- as.numeric(crossprod(sw, centered^2))
  names(sv_var) <- sv_names
  scores <- rep(0, length(feature_names))
  names(scores) <- feature_names
  shared <- intersect(feature_names, names(sv_var))
  scores[shared] <- sv_var[shared]
  scores
}

select_explain_feature_indices <- function(pred_model, feature_names, n_candidates = 300) {
  p <- length(feature_names)
  n_candidates <- as.integer(n_candidates)
  if (is.na(n_candidates) || n_candidates <= 0) n_candidates <- min(300L, p)
  if (n_candidates >= p) return(seq_len(p))
  scores <- compute_svm_feature_screen_scores(pred_model, feature_names)
  ord <- order(scores, decreasing = TRUE)
  ord[seq_len(n_candidates)]
}

explain_sample_local_sensitivity <- function(
    sample_name, x_row, predicted_class_raw, predicted_class_display,
    top_prob, pred_model, cal_model, candidate_idx,
    chunk_size = 50, eps = 0.02
) {
  feature_names <- names(x_row)
  x_row <- as.numeric(x_row)
  names(x_row) <- feature_names
  if (is.null(feature_names)) stop("x_row must be a named vector.")
  if (length(candidate_idx) == 0) return(tibble())
  chunk_size <- max(1L, as.integer(chunk_size))
  eps <- as.numeric(eps)
  if (!is.finite(eps) || eps <= 0) stop("eps must be > 0.")
  dp_dx <- numeric(length(candidate_idx))
  p_plus <- numeric(length(candidate_idx))
  p_minus <- numeric(length(candidate_idx))
  step <- numeric(length(candidate_idx))
  for (start in seq(1L, length(candidate_idx), by = chunk_size)) {
    stop_i <- min(start + chunk_size - 1L, length(candidate_idx))
    idx_chunk <- candidate_idx[start:stop_i]
    n_chunk <- length(idx_chunk)
    x_both <- matrix(
      rep(x_row, times = n_chunk),
      nrow = n_chunk, byrow = TRUE,
      dimnames = list(NULL, feature_names)
    )
    x_both <- rbind(x_both, x_both)
    step_chunk <- numeric(n_chunk)
    for (k in seq_len(n_chunk)) {
      j <- idx_chunk[[k]]
      plus_val <- min(1, x_row[[j]] + eps)
      minus_val <- max(0, x_row[[j]] - eps)
      x_both[k, j] <- plus_val
      x_both[k + n_chunk, j] <- minus_val
      step_chunk[k] <- plus_val - minus_val
    }
    probs_both <- predict_cal_probs(
      cal_model,
      predict_prob_SVM(pred_model$model, x_both)
    )
    cls_col <- predicted_class_raw
    if (!cls_col %in% colnames(probs_both)) cls_col <- predicted_class_display
    if (!cls_col %in% colnames(probs_both)) {
      dp_dx[start:stop_i] <- NA
      p_plus[start:stop_i] <- NA
      p_minus[start:stop_i] <- NA
      step[start:stop_i] <- step_chunk
      next
    }
    p_p <- probs_both[seq_len(n_chunk), cls_col]
    p_m <- probs_both[seq_len(n_chunk) + n_chunk, cls_col]
    dp_dx[start:stop_i] <- ifelse(step_chunk > 0, (p_p - p_m) / step_chunk, 0)
    p_plus[start:stop_i] <- p_p
    p_minus[start:stop_i] <- p_m
    step[start:stop_i] <- step_chunk
  }
  tibble(
    sample = sample_name,
    predicted_class = predicted_class_display,
    top_class_probability = top_prob,
    cpg = feature_names[candidate_idx],
    sample_beta = x_row[candidate_idx],
    eps = eps,
    p_plus = p_plus,
    p_minus = p_minus,
    dp_dx = dp_dx,
    abs_dp_dx = abs(dp_dx),
    local_effect_for_eps = abs(dp_dx * eps),
    direction = ifelse(dp_dx >= 0,
                       "increasing_beta_increases_top_prob",
                       "increasing_beta_decreases_top_prob")
  ) %>%
    arrange(desc(abs_dp_dx))
}

# --- 7. Apply SVM classifier -------------------------------------------------

cat("\nApplying SVM classifier...\n")

raw_probs <- predict_prob_SVM(pred_model$model, x)

# Raw (native model label space) and display (renamed) calibrated probabilities
cal_probs_raw <- predict_cal_probs(cal_model, raw_probs)
cal_probs <- apply_target_class_names(cal_probs_raw)

if (is.null(rownames(cal_probs_raw)) || any(!nzchar(rownames(cal_probs_raw)))) {
  rownames(cal_probs_raw) <- rownames(x)
}
if (is.null(rownames(cal_probs)) || any(!nzchar(rownames(cal_probs)))) {
  rownames(cal_probs) <- rownames(x)
}

top_score <- apply(cal_probs, 1, max)
top_class_raw <- colnames(cal_probs_raw)[apply(cal_probs_raw, 1, which.max)]
top_class <- normalize_class_names(colnames(cal_probs)[apply(cal_probs, 1, which.max)])

# --- 8. Compile results ------------------------------------------------------

confidence_tier <- case_when(
  top_score >= 0.90 ~ "Strong (>=0.9)",
  top_score >= 0.84 ~ "Acceptable (0.84-0.9)",
  top_score >= 0.50 ~ "Suggestive (0.5-0.84)",
  TRUE ~ "Non-informative (<0.5)"
)

results_summary <- tibble(
  sample_id = rownames(cal_probs),
  predicted_class = top_class,
  confidence_score = round(top_score, 4),
  confidence_tier = confidence_tier
)

cat("\n=== Classification Results ===\n")
print(results_summary)

cat("\nConfidence tier distribution:\n")
print(table(results_summary$confidence_tier))

cat("\nPredicted class distribution:\n")
print(table(results_summary$predicted_class))

# --- 9. Merge with QC -------------------------------------------------------

results_with_qc <- results_summary %>%
  left_join(
    qc %>% select(sample_id, frac_na, frac_unmeth_ch, mean_intensity, qc_flag),
    by = "sample_id"
  ) %>%
  left_join(
    tibble(sample_id = names(frac_imputed), frac_imputed = unname(frac_imputed)),
    by = "sample_id"
  ) %>%
  arrange(desc(confidence_score))

cat("\n=== Results with QC metrics ===\n")
print(results_with_qc)

cat("\nMean confidence score by QC flag:\n")
print(results_with_qc %>%
        group_by(qc_flag) %>%
        summarise(
          n = n(),
          mean_confidence = round(mean(confidence_score, na.rm = TRUE), 3),
          n_unknown = sum(predicted_class == "Unknown"),
          n_strong = sum(confidence_tier == "Strong (>=0.9)"),
          n_acceptable = sum(confidence_tier == "Acceptable (0.84-0.9)"),
          n_suggestive = sum(confidence_tier == "Suggestive (0.5-0.84)"),
          n_non_informative = sum(confidence_tier == "Non-informative (<0.5)")
        ))

# --- 10. CpG sensitivity analysis --------------------------------------------

# Analysis parameters
explain_contributions <- TRUE
explain_candidate_cpgs <- 300L
explain_top_cpgs_per_sample <- 25L
explain_chunk_size <- 50L
explain_eps <- 0.02

contrib_full_tbl <- tibble()
contrib_tbl <- tibble()
contrib_status <- "disabled"
candidate_cpgs_used <- 0L
estimated_counterfactual_rows <- 0L

if (isTRUE(explain_contributions)) {
  
  # Screen candidate CpGs for sensitivity analysis
  candidate_idx <- select_explain_feature_indices(
    pred_model = pred_model,
    feature_names = colnames(x),
    n_candidates = explain_candidate_cpgs
  )
  candidate_cpgs_used <- length(candidate_idx)
  estimated_counterfactual_rows <- as.integer(2L * candidate_cpgs_used * nrow(x))
  contrib_list <- vector("list", nrow(x))
  
  pred_lookup <- setNames(top_class_raw, rownames(cal_probs_raw))
  pred_display_lookup <- setNames(results_summary$predicted_class, results_summary$sample_id)
  score_lookup <- setNames(results_summary$confidence_score, results_summary$sample_id)
  
  # Per-sample local sensitivity (+/- eps perturbation per candidate CpG)
  for (i in seq_len(nrow(x))) {
    sid <- rownames(x)[[i]]
    if (is.null(sid) || !nzchar(sid)) sid <- sample_ids[[i]]
    
    cls_raw <- unname(pred_lookup[[sid]])
    cls_display <- unname(pred_display_lookup[[sid]])
    p_top <- unname(score_lookup[[sid]])
    
    if (is.null(cls_raw) || is.null(p_top)) {
      cls_raw <- as.character(top_class_raw[[i]])
      cls_display <- as.character(top_class[[i]])
      p_top <- as.numeric(top_score[[i]])
    }
    
    contrib_list[[i]] <- explain_sample_local_sensitivity(
      sample_name = sid,
      x_row = x[i, ],
      predicted_class_raw = cls_raw,
      predicted_class_display = cls_display,
      top_prob = p_top,
      pred_model = pred_model,
      cal_model = cal_model,
      candidate_idx = candidate_idx,
      chunk_size = explain_chunk_size,
      eps = explain_eps
    )
  }
  
  contrib_full_tbl <- bind_rows(contrib_list)
  contrib_tbl <- contrib_full_tbl %>%
    group_by(sample) %>%
    slice_max(order_by = abs_dp_dx, n = as.integer(explain_top_cpgs_per_sample), with_ties = FALSE) %>%
    ungroup()
  contrib_status <- "computed"
}

tibble(
  contribution_analysis = contrib_status,
  explain_eps = as.numeric(explain_eps),
  explain_candidate_cpgs = as.integer(explain_candidate_cpgs),
  candidate_cpgs_used = as.integer(candidate_cpgs_used),
  estimated_counterfactual_rows = estimated_counterfactual_rows,
  explain_top_cpgs_per_sample = as.integer(explain_top_cpgs_per_sample),
  explain_chunk_size = as.integer(explain_chunk_size)
)

contrib_tbl

# --- 11. Save outputs --------------------------------------------------------

prob_tbl <- as.data.frame(cal_probs) %>%
  rownames_to_column("sample_id") %>%
  as_tibble()

saveRDS(model_feature_betas,file.path(out_results, "model_feature_betas_imputed.rds"))
write.csv(results_summary, file.path(out_results, "classifier_predictions_summary.csv"), row.names = FALSE)
write.csv(prob_tbl, file.path(out_results, "classifier_probabilities_full.csv"), row.names = FALSE)
write.csv(results_with_qc, file.path(out_results, "classifier_results_with_qc.csv"), row.names = FALSE)
writeLines(probes_missing, file.path(out_results, "classifier_cpgs_missing.txt"))
if (nrow(contrib_tbl) > 0) {
  write.csv(contrib_tbl, file.path(out_results, "top_class_cpg_contributions_top.csv"), row.names = FALSE)
}
if (nrow(contrib_full_tbl) > 0) {
  write.csv(contrib_full_tbl, file.path(out_results, "top_class_cpg_contributions_full.csv"), row.names = FALSE)
}

cat("\nOutputs saved to:", out_results, "\n")

# --- 12. Session info --------------------------------------------------------

sink(file.path(out_results, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Classifier complete ===\n")
