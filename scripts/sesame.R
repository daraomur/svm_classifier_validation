concordance_summary <- data.frame(
  category = c("Diagnoses confirmed", "Diagnoses reclassified", "Non-classifiable (Unknown)"),
  count = c(12, 4, 22)
)

concordance_summary$category <- factor(concordance_summary$category,
                                       levels = c("Diagnoses confirmed", "Diagnoses reclassified", "Non-classifiable (Unknown)"))

p_pie <- ggplot(concordance_summary, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, colour = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(count, "\n(", round(100 * count / sum(count), 1), "%)")),
            position = position_stack(vjust = 0.5), size = 6, colour = "white") +
  scale_fill_manual(values = c(
    "Diagnoses confirmed" = "#154069",
    "Diagnoses reclassified" = "#3B8BD4",
    "Non-classifiable (Unknown)" = "#B4B2A9"
  )) +
  labs(title = "Concordance Outcomes", fill = "Category") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_pie)
ggsave(file.path(out_qc, "plot_concordance_pie.pdf"), p_pie, width = 6, height = 6)

# =============================================================================
# Sinonasal Methylation Classifier Validation — SeSAMe Processing Pipeline
# SPELCASTER WP1b | RCSI
# Author: Dara O'Murchu
# Date: May 2026
# =============================================================================

# --- 0. Libraries ------------------------------------------------------------

library(sesame)
library(sesameData)
library(readxl)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pals)
library(KernSmooth)

# --- 1. Paths ----------------------------------------------------------------

base_dir <- "~/Desktop/Research Project/My Project/Dara_Validation"
idat_dir <- file.path(base_dir, "data/SinonasalMethyArray_Batch01_2026")
sample_sheet_path <- file.path(base_dir, "data/SinonasalMethyArray_Batch01_2026/Sample IDs.xlsx")

out_betas <- file.path(base_dir, "results/beta_matrices")
out_qc <- file.path(base_dir, "results/qc")

dir.create(out_betas, recursive = TRUE, showWarnings = FALSE)
dir.create(out_qc, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load sample sheet ----------------------------------------------------

sample_sheet <- read_excel(sample_sheet_path) %>%
  dplyr::select(1, 2) %>%
  setNames(c("sample_id", "file_name")) %>%
  mutate(
    sample_id = trimws(sample_id),
    chip_id = sub("_R[0-9]+C[0-9]+$", "", file_name),
    position = sub(".*_(R[0-9]+C[0-9]+)$", "\\1", file_name)
  )

# Typo correction for SPL-DU-00057
sample_sheet <- sample_sheet %>%
  mutate(file_name = ifelse(file_name == "2102522410105_R08C01", # Removed additional 1
                            "210252240105_R08C01",
                            file_name))

print(sample_sheet)
cat("Total samples in sheet:", nrow(sample_sheet), "\n")

# --- 3. Discover IDAT files --------------------------------------------------

idat_prefixes <- searchIDATprefixes(idat_dir)
cat("IDAT prefixes found:", length(idat_prefixes), "\n")

idat_basenames <- basename(idat_prefixes)
missing_idats <- sample_sheet$file_name[!sample_sheet$file_name %in% idat_basenames]
if (length(missing_idats) > 0) {
  warning("No IDAT found for:\n", paste0(" ", missing_idats, collapse = "\n"))
}

# --- 4. Helper functions -----------------------------------------------------

# Safely extract metrics from sesameQC S4 object or plain list
qc_stats_to_vector <- function(qc_obj) {
  if (methods::is(qc_obj, "sesameQC")) {
    return(unlist(qc_obj@stat, use.names = TRUE))
  }
  unlist(qc_obj, use.names = TRUE)
}

# Collapse EPICv2 replicate cg-prefixed probes to their mean
collapse_epic_v2_probe_prefix <- function(beta_vec) {
  if (is.null(names(beta_vec))) stop("beta_vec must be a named numeric vector.")
  
  probe_ids <- names(beta_vec)
  probe_prefix <- ifelse(
    grepl("^cg", probe_ids),
    sub("_.*$", "", probe_ids),
    probe_ids
  )
  
  if (!anyDuplicated(probe_prefix)) return(beta_vec)
  
  split_vec <- split(as.numeric(beta_vec), probe_prefix)
  vapply(split_vec, function(x) {
    m <- mean(x, na.rm = TRUE)
    if (is.nan(m)) NA_real_ else m
  }, numeric(1))
}

# --- 5. Per-sample processing loop -------------------------------------------
# Each sample processed individually:
# QC on raw SigDF, QCDPB preprocessing, EPICv2 -> EPICv1 liftover, probe collapse.

mft <- sesameAnno_buildManifestGRanges(
  sesameAnno_download("EPICv2.hg38.manifest.tsv.gz"),
  columns = "nextBase"
)
extR <- names(mft)[!is.na(mft$nextBase) & mft$nextBase == "R"]
extA <- names(mft)[!is.na(mft$nextBase) & mft$nextBase == "A"]

cat("\nProcessing samples individually...\n")

betas_list <- vector("list", nrow(sample_sheet))
qc_list <- vector("list", nrow(sample_sheet))
bis_conv_list <- vector("list", nrow(sample_sheet))
names(betas_list) <- sample_sheet$sample_id

for (i in seq_len(nrow(sample_sheet))) {
  
  sid <- sample_sheet$sample_id[i]
  file_name <- sample_sheet$file_name[i]
  prefix <- idat_prefixes[basename(idat_prefixes) == file_name]
  
  if (length(prefix) == 0) {
    warning("No IDAT prefix found for: ", sid, " (", file_name, ")")
    next
  }
  
  cat(sprintf("  [%d/%d] %s\n", i, nrow(sample_sheet), sid))
  
  sdf <- readIDATpair(prefix)
  
  # QC metrics extracted from raw SigDF before preprocessing
  qc_vals <- qc_stats_to_vector(sesame::sesameQC_calcStats(sdf))
  qc_list[[i]] <- as_tibble(as.list(qc_vals)) %>%
    mutate(
      sample_id = sid,
      file_name = file_name,
      idat_prefix = prefix,
      array_version = tryCatch(sdfPlatform(sdf), error = function(e) NA_character_),
      .before = 1
    )
  
  # GCT bisulfite conversion score from raw EPICv2 SigDF
  bis_conv_list[[i]] <- tryCatch(
    bisConversionControl(sdf, extR, extA),
    error = function(e) NA_real_
  )
  
  tryCatch({
    pdf(file.path(out_qc, paste0(sid, "_intens_vs_betas.pdf")), width = 7, height = 6)
    sesameQC_plotIntensVsBetas(sdf)
    dev.off()
  }, error = function(e) {
    dev.off()
    warning("sesameQC_plotIntensVsBetas failed for ", sid, ": ", conditionMessage(e))
  })
  
  # Preprocess:
  # Q=qualityMask, C=inferInfiniumIChannel, D=dyeBiasNL, P=pOOBAH, B=noob
  sdf_prep <- prepSesame(sdf, prep = "QCDPB")
  beta_vec <- getBetas(sdf_prep)
  
  # Liftover EPICv2 -> EPICv1 probe space (hard stop on failure)
  beta_vec <- tryCatch(
    sesame::mLiftOver(beta_vec, target = "EPIC", impute = FALSE),
    error = function(e) {
      stop("mLiftOver failed for ", sid, ": ", conditionMessage(e),
           "\nCache the resource: sesameData::sesameDataCache('anno.hg19.EPIC')")
    }
  )
  
  beta_vec <- collapse_epic_v2_probe_prefix(beta_vec)
  
  betas_list[[i]] <- beta_vec
  
  rm(sdf, sdf_prep, qc_vals, beta_vec)
  gc(verbose = FALSE)
}

# --- 6. Post-loop liftover verification --------------------------------------

probe_counts <- sapply(betas_list[!sapply(betas_list, is.null)], length)
cat("\nProbe counts per sample after liftover:\n")
print(table(probe_counts))

liftover_failures <- names(probe_counts)[probe_counts > 900000]
if (length(liftover_failures) > 0) {
  stop("Liftover failed for: ", paste(liftover_failures, collapse = ", "))
} else {
  cat("All samples successfully lifted over to EPICv1 probe space.\n")
}

# --- 7. Assemble beta matrix -------------------------------------------------

cat("\nAssembling beta matrix...\n")

all_probes <- Reduce(intersect,
                     lapply(betas_list[!sapply(betas_list, is.null)], names))
cat("Probes present in all samples:", length(all_probes), "\n")

betas <- do.call(cbind, lapply(betas_list, function(bv) {
  if (is.null(bv)) return(rep(NA_real_, length(all_probes)))
  bv[all_probes]
}))

cat("Beta matrix dimensions:", dim(betas), "\n")

# --- 8. Assemble QC table ----------------------------------------------------

qc_annotated <- bind_rows(qc_list) %>%
  mutate(
    bis_conv_score = as.numeric(bis_conv_list),
    qc_flag = case_when(
      is.na(frac_na) | is.na(bis_conv_score) ~ "FLAGGED",
      frac_na > 0.2 | bis_conv_score > 1.5 ~ "FLAGGED",
      TRUE ~ "PASS"
    ),
    qc_reason = case_when(
      is.na(frac_na) | is.na(bis_conv_score) ~ "missing metric",
      frac_na > 0.2 & bis_conv_score > 1.5 ~ "high NA + elevated GCT",
      frac_na > 0.2 ~ "high NA",
      bis_conv_score > 1.5 ~ "elevated GCT",
      TRUE ~ "pass"
    )
  )

cat("\nQC Summary:\n")
print(qc_annotated %>% select(sample_id, array_version, frac_na, bis_conv_score, qc_flag, qc_reason))
cat("\nQC flag counts:\n")
print(table(qc_annotated$qc_flag))

# --- 9. QC visualisation -----------------------------------------------------

# Fraction NA per sample
p_na <- ggplot(qc_annotated,
               aes(x = reorder(sample_id, frac_na),
                   y = frac_na,
                   fill = qc_flag)) +
  geom_col() +
  scale_fill_manual(values = c("PASS" = "#66BB6A", "FLAGGED" = "#EF5350")) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "grey40") +
  coord_flip() +
  labs(title = "Fraction of NA beta values per preprocessed sample",
       subtitle = "Dashed line = 20% threshold",
       x = NULL, y = "Fraction NA", fill = "QC status") +
  theme_minimal(base_size = 10)

print(p_na)
ggsave(file.path(out_qc, "plot_frac_na.pdf"), p_na, width = 8, height = 10)

# Bisulfite conversion control (GCT) per sample
p_conv <- ggplot(qc_annotated,
                 aes(x = reorder(sample_id, bis_conv_score),
                     y = bis_conv_score,
                     fill = qc_flag)) +
  geom_col() +
  scale_fill_manual(values = c("PASS" = "#66BB6A", "FLAGGED" = "#EF5350")) +
  geom_hline(yintercept = 1.5, linetype = "dashed", colour = "grey40") +
  coord_flip() +
  labs(title = "GCT bisulfite-conversion control score per sample",
       subtitle = "Dashed line = GCT > 1.5 flagging threshold",
       x = NULL,
       y = "GCT bisulfite-conversion score",
       fill = "QC status") +
  theme_minimal(base_size = 10)

print(p_conv)
ggsave(file.path(out_qc, "plot_bisulfite_conversion.pdf"), p_conv, width = 8, height = 10)

# Joint QC scatter plot: probe detection vs bisulfite conversion
p_joint <- ggplot(qc_annotated,
                  aes(x = bis_conv_score,
                      y = frac_na,
                      colour = qc_flag,
                      label = sample_id)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_manual(values = c("PASS" = "#66BB6A", "FLAGGED" = "#EF5350")) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey40") +
  labs(title = "Sample QC: probe detection vs bisulfite conversion",
       subtitle = "Dashed lines = fraction NA > 20% or GCT > 1.5",
       x = "GCT bisulfite-conversion score",
       y = "Fraction NA probes",
       colour = "QC status") +
  theme_minimal(base_size = 10)

print(p_joint)
ggsave(file.path(out_qc, "plot_qc_joint.pdf"), p_joint, width = 9, height = 7)


# Beta value density plot
set.seed(42)
probe_sample <- sample(rownames(betas), 100000)

betas_long <- betas[probe_sample, ] %>%
  as.data.frame() %>%
  rownames_to_column("probe_id") %>%
  pivot_longer(-probe_id, names_to = "sample_id", values_to = "beta") %>%
  left_join(qc_annotated %>% select(sample_id, qc_flag), by = "sample_id") %>%
  filter(!is.na(beta))

## Highlight 8 poorest quality samples by frac_na
samples_of_interest <- qc_annotated %>%
  arrange(desc(frac_na)) %>%
  slice_head(n = 8) %>%
  pull(sample_id)

interest_colours <- setNames(
  scales::hue_pal()(length(samples_of_interest)),
  samples_of_interest
)

betas_long <- betas_long %>%
  mutate(highlight = sample_id %in% samples_of_interest)

p_density <- ggplot() +
  geom_density(
    data = betas_long %>% filter(!highlight),
    aes(x = beta, group = sample_id),
    colour = "grey80",
    linewidth = 0.3
  ) +
  geom_density(
    data = betas_long %>% filter(highlight),
    aes(x = beta, group = sample_id, colour = sample_id),
    linewidth = 0.6
  ) +
  scale_colour_manual(values = interest_colours) +
  labs(
    title = "Beta value density per sample",
    subtitle = "Highlighted samples exhibit highest fraction of NA probes",
    x = "Beta value",
    y = "Density",
    colour = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.width = unit(1.5, "cm")
  ) +
  guides(colour = guide_legend(nrow = 2))

print(p_density)
ggsave(file.path(out_qc, "plot_beta_density.pdf"), p_density, width = 9, height = 7)
cat("\nQC plots saved.\n")

unknown_samples <- c(
  "SPL-DU-00003", "SPL-DU-00006", "SPL-DU-00009", "SPL-DU-00011",
  "SPL-DU-00012", "SPL-DU-00015", "SPL-DU-00022", "SPL-DU-00023",
  "SPL-DU-00024", "SPL-DU-00031", "SPL-DU-00034", "SPL-DU-00038",
  "SPL-DU-00039", "SPL-DU-00041", "SPL-DU-00042", "SPL-DU-00043",
  "SPL-DU-00045", "SPL-DU-00048", "SPL-DU-00050", "SPL-DU-00052",
  "SPL-DU-00053", "SPL-DU-00056"
)

qc_annotated <- qc_annotated %>%
  mutate(classifier_call = ifelse(sample_id %in% unknown_samples,
                                  "Unknown", "Classified"))

p_joint <- ggplot(qc_annotated,
                  aes(x = bis_conv_score,
                      y = frac_na,
                      colour = qc_flag,
                      shape = classifier_call,
                      label = sample_id)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_manual(values = c("PASS" = "#66BB6A", "FLAGGED" = "#EF5350")) +
  scale_shape_manual(values = c("Classified" = 16, "Unknown" = 17)) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey40") +
  labs(title = "Sample QC: probe detection vs bisulfite conversion",
       subtitle = "Dashed lines = fraction NA > 20% or GCT > 1.5, triangle = Unknown",
       x = "GCT bisulfite-conversion score",
       y = "Fraction NA probes",
       colour = "QC status",
       shape = "Classifier result") +
  theme_minimal(base_size = 10)
print(p_joint)
ggsave(file.path(out_qc, "plot_qc_joint.pdf"), p_joint, width = 9, height = 7)

# --- 10. SNP heatmap ---------------------------------------------------------

sdfs <- vector("list", nrow(sample_sheet))
names(sdfs) <- sample_sheet$sample_id

for (i in seq_len(nrow(sample_sheet))) {
  sid <- sample_sheet$sample_id[i]
  prefix <- idat_prefixes[basename(idat_prefixes) == sample_sheet$file_name[i]]
  if (length(prefix) == 0) next
  sdfs[[i]] <- readIDATpair(prefix)
}

sdfs <- Filter(Negate(is.null), sdfs)

pdf(file.path(out_qc, "plot_snp_heatmap.pdf"), width = 10, height = 8)
sesameQC_plotHeatSNPs(sdfs)
dev.off()

rm(sdfs)
gc(verbose = FALSE)

# --- 11. Save outputs --------------------------------------------------------

saveRDS(betas, file.path(out_betas, "batch01_betas_QCDPB_lifted.rds"))
write.csv(betas, file.path(out_betas, "batch01_betas_QCDPB_lifted.csv"), quote = FALSE)
write.csv(qc_annotated, file.path(out_qc, "batch01_qc_summary.csv"), row.names = FALSE)

cat("\nOutputs saved to:", base_dir, "\n")

# --- 12. Session info --------------------------------------------------------

sink(file.path(out_qc, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Pipeline complete ===\n")
cat("Beta matrix:", nrow(betas), "probes x", ncol(betas), "samples\n")
cat("Samples passing QC:", sum(qc_annotated$qc_flag == "PASS", na.rm = TRUE), "\n")
cat("Samples flagged:   ", sum(qc_annotated$qc_flag == "FLAGGED", na.rm = TRUE), "\n")
cat("Samples with no IDAT match:", nrow(sample_sheet) - nrow(qc_annotated), "\n")