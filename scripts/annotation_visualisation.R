# =============================================================================
# Sinonasal Methylation Classifier Validation — Sample Annotation & Visualisation
# SPELCASTER WP1b | RCSI
# Author: Dara O'Murchu
# Date: May 2026
# =============================================================================

# --- 0. Libraries ------------------------------------------------------------

library(ExperimentHub)
library(umap)
library(Rtsne)
library(HiTIMED)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(data.table)
library(sesame)
library(sesameData)
library(readxl)
library(GenomicRanges)

# --- 1. Paths ----------------------------------------------------------------

base_dir <- "~/Desktop/Research Project/My Project/Dara_Validation"
jur_dir <- "~/Desktop/Research Project/My Project/Dara_Validation/data/jurmeister_2022"

idat_dir <- file.path(base_dir, "data/SinonasalMethyArray_Batch01_2026")
sample_sheet_path <- file.path(base_dir, "data/SinonasalMethyArray_Batch01_2026/Sample IDs.xlsx")

out_dimred <- file.path(base_dir, "results/dimensionality_reduction")
out_purity <- file.path(base_dir, "results/tumour_purity")
out_cnv <- file.path(base_dir, "results/cnv")

dir.create(out_dimred, recursive = TRUE, showWarnings = FALSE)
dir.create(out_purity, recursive = TRUE, showWarnings = FALSE)
dir.create(out_cnv, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load inputs ----------------------------------------------------------

betas_full <- readRDS(file.path(base_dir, "results/beta_matrices/batch01_betas_QCDPB_lifted.rds"))
cat("Full beta matrix loaded:", nrow(betas_full), "probes x", ncol(betas_full), "samples\n")

results <- read.csv(file.path(base_dir, "results/classifier/classifier_results_with_qc.csv"))

# --- 3. Load sample sheet and IDAT prefixes (for CNV) ------------------------

sample_sheet <- read_excel(sample_sheet_path) %>%
  select(1, 2) %>%
  setNames(c("sample_id", "file_name")) %>%
  mutate(
    sample_id = trimws(sample_id),
    file_name = ifelse(file_name == "2102522410105_R08C01", "210252240105_R08C01", file_name),
    chip_id = sub("_R[0-9]+C[0-9]+$", "", file_name),
    position = sub(".*_(R[0-9]+C[0-9]+)$", "\\1", file_name)
  )

idat_prefixes <- searchIDATprefixes(idat_dir)
cat("IDAT prefixes found:", length(idat_prefixes), "\n")

# --- 4. Load Jurmeister reference data ---------------------------------------

cat("Loading Jurmeister beta matrix...\n")
jur_path <- path.expand(file.path(jur_dir, "jurmeister_full_betas_round3.csv"))
jur_betas_raw <- fread(file = jur_path, data.table = FALSE)
rownames(jur_betas_raw) <- jur_betas_raw$probe
jur_betas_raw$probe <- NULL
jur_betas_mat <- as.matrix(jur_betas_raw)
cat("Jurmeister beta matrix loaded:", nrow(jur_betas_mat), "probes x", ncol(jur_betas_mat), "samples\n")

jur_samples <- read.csv(file.path(jur_dir, "jurmeister_tumor_betas_v2_sample_list.csv"))
cat("Jurmeister sample annotation loaded:", nrow(jur_samples), "samples\n")

# Remove Test_set samples (no class annotation)
test_set_samples <- colnames(jur_betas_mat)[grepl("^Test_set", colnames(jur_betas_mat))]
cat("Removing", length(test_set_samples), "unannotated Test_set samples\n")
jur_betas_mat_cleaned <- jur_betas_mat[, !colnames(jur_betas_mat) %in% test_set_samples]
cat("Jurmeister matrix after removal:", ncol(jur_betas_mat_cleaned), "samples\n")

excluded_sample_titles <- jur_samples$sample_title[jur_samples$excluded_supp %in% c("Yes", "Noise")]
cat("Excluding", length(excluded_sample_titles), "Jurmeister reference samples flagged in excluded_supp\n")

jur_betas_mat_cleaned <- jur_betas_mat_cleaned[, !colnames(jur_betas_mat_cleaned) %in% excluded_sample_titles]
cat("Jurmeister matrix after exclusion filtering:", ncol(jur_betas_mat_cleaned), "samples\n")

# --- 5. Tumour purity estimation via HiTIMED ---------------------------------
# tumor_type = "HNSC" - head and neck squamous cell carcinoma
# closest available reference to sinonasal tumours in HiTIMED

cat("\nRunning HiTIMED deconvolution...\n")

hitimed_result <- HiTIMED_deconvolution(
  tumor_beta = betas_full,
  tumor_type = "HNSC",
  h = 6,
  tissue_type = "tumor"
)

cat("HiTIMED output dimensions:", dim(hitimed_result), "\n")
cat("Cell types estimated:\n")
print(colnames(hitimed_result))

tumour_purity <- as.data.frame(hitimed_result) %>%
  rownames_to_column("sample_id") %>%
  select(sample_id, tumour_purity_pct = Tumor) %>%
  mutate(tumour_purity_pct = round(tumour_purity_pct, 1))

cat("\nTumour purity estimates:\n")
print(tumour_purity %>% arrange(desc(tumour_purity_pct)))

# Tumour purity visualisation 

purity_plot_data <- tumour_purity %>%
  left_join(
    results %>% select(sample_id, confidence_score, confidence_tier, qc_flag, predicted_class),
    by = "sample_id"
  )

qc_colours <- c("PASS" = "#2196F3", "FLAGGED" = "#F44336")

# Purity distribution by QC flag 

p_purity_dist <- ggplot(purity_plot_data, aes(x = tumour_purity_pct, fill = qc_flag)) +
  geom_histogram(binwidth = 8, colour = "white", linewidth = 0.5, alpha = 0.8,
                 position = "stack") +
  geom_vline(
    xintercept = median(purity_plot_data$tumour_purity_pct, na.rm = TRUE),
    linetype = "dashed", colour = "black", linewidth = 0.8
  ) +
  annotate(
    "text",
    x = median(purity_plot_data$tumour_purity_pct, na.rm = TRUE) + 2,
    y = Inf, vjust = 1.5, hjust = 0, size = 3.2, colour = "black",
    label = paste0("Median: ",
                   round(median(purity_plot_data$tumour_purity_pct, na.rm = TRUE), 0), "%")
  ) +
  scale_fill_manual(values = qc_colours) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Tumour Purity Distribution",
    subtitle = "HiTIMED deconvolution (HNSC reference) | n = 38",
    x = "Tumour Purity (%)",
    y = "Number of Samples",
    fill = "QC Flag"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
    )

ggsave(file.path(out_purity, "plot_purity_distribution.pdf"),
       p_purity_dist, width = 7, height = 5)
cat("Purity distribution plot saved.\n")

# Purity vs confidence score scatter 

# Samples to annotate — key reclassified cases and high-imputation outliers
annotate_sids <- c("SPL-DU-00010", "SPL-DU-00027", "SPL-DU-00023", "SPL-DU-00034 T")

p_purity_conf <- ggplot(purity_plot_data,
                        aes(x = tumour_purity_pct, y = confidence_score,
                            colour = qc_flag)) +
  # Confidence tier threshold lines
  geom_hline(yintercept = c(0.5, 0.84, 0.9),
             linetype = "dotted", colour = "grey60", linewidth = 0.6) +
  annotate("text", x = 101, y = 0.95, label = "Strong",
           size = 2.8, colour = "grey50", hjust = 0) +
  annotate("text", x = 101, y = 0.87, label = "Acceptable",
           size = 2.8, colour = "grey50", hjust = 0) +
  annotate("text", x = 101, y = 0.67, label = "Suggestive",
           size = 2.8, colour = "grey50", hjust = 0) +
  annotate("text", x = 101, y = 0.25, label = "Non-informative",
           size = 2.8, colour = "grey50", hjust = 0) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_text_repel(
    data = purity_plot_data %>% filter(sample_id %in% annotate_sids),
    aes(label = paste0(sample_id, "\n(", predicted_class, ")")),
    size = 2.5, colour = "grey20", show.legend = FALSE,
    box.padding = 0.5, max.overlaps = 20
  ) +
  scale_colour_manual(values = qc_colours) +
  scale_x_continuous(limits = c(0, 115), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Tumour Purity vs Classifier Confidence Score",
    subtitle = "Coloured by QC flag | Key samples annotated",
    x = "Tumour Purity (%)",
    y = "Classifier Confidence Score",
    colour = "QC Flag"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(out_purity, "plot_purity_vs_confidence.pdf"),
       p_purity_conf, width = 8, height = 6)
cat("Purity vs confidence scatter saved.\n")

# --- 6. Build annotation table -----------------------------------------------

annot <- results %>%
  select(sample_id, predicted_class, confidence_score, confidence_tier, frac_na, qc_flag) %>%
  mutate(sample_id = trimws(sample_id)) %>%
  left_join(tumour_purity, by = "sample_id")

# --- 7. RCSI-only UMAP and t-SNE ---------------------------------------------
# Confidence score, tumour purity, and predicted class.
# QC overlay covered by combined plots below.

probe_vars <- apply(betas_full, 1, var, na.rm = TRUE)
top_var_probes <- names(sort(probe_vars, decreasing = TRUE))[1:20000]
cat("Most variable probes selected:", length(top_var_probes), "\n")

x <- t(betas_full[top_var_probes, ])
cat("Input matrix for UMAP/t-SNE:", nrow(x), "samples x", ncol(x), "features\n")

x_imp <- apply(x, 2, function(probe) {
  na_idx <- is.na(probe)
  if (!any(na_idx)) return(probe)
  fill <- mean(probe, na.rm = TRUE)
  if (!is.finite(fill)) fill <- 0.5
  probe[na_idx] <- fill
  probe
})
cat("NAs remaining after imputation:", sum(is.na(x_imp)), "\n")

# --- 7a. RCSI UMAP -----------------------------------------------------------

set.seed(42)
umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"

cat("\nRunning UMAP on 20,000 most variable probes (RCSI only)...\n")
umap_result <- umap(x_imp, config = umap_config)

umap_df <- as.data.frame(umap_result$layout) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(sample_id = rownames(x_imp)) %>%
  left_join(annot, by = "sample_id")

# Plot - confidence score
p_umap_conf <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = confidence_score, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_gradient(low = "#EF5350", high = "#66BB6A") +
  labs(
    title = "UMAP of sinonasal tumour methylation profiles",
    subtitle = "20,000 most variable probes | Coloured by classifier confidence score",
    colour = "Confidence score"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_umap_confidence.pdf"), p_umap_conf, width = 9, height = 7)

# Plot - tumour purity
p_umap_purity <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = tumour_purity_pct, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_gradient(low = "#EF5350", high = "#66BB6A") +
  labs(
    title = "UMAP of sinonasal tumour methylation profiles",
    subtitle = "20,000 most variable probes | Coloured by tumour purity (HiTIMED, HNSC reference)",
    colour = "Tumour purity (%)"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_umap_tumour_purity.pdf"), p_umap_purity, width = 9, height = 7)

# Plot - predicted class
p_umap_class <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = predicted_class, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  labs(
    title = "UMAP of sinonasal tumour methylation profiles",
    subtitle = "20,000 most variable probes | Coloured by predicted methylation class",
    colour = "Predicted class"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_umap_class.pdf"), p_umap_class, width = 9, height = 7)

cat("\nRCSI UMAP plots saved.\n")

# --- 7b. RCSI t-SNE ----------------------------------------------------------

set.seed(42)
cat("\nRunning t-SNE on 20,000 most variable probes (RCSI only)...\n")

x_unique <- unique(x_imp)
cat("Unique samples for t-SNE:", nrow(x_unique), "\n")

tsne_result <- Rtsne(
  x_unique,
  dims = 2,
  perplexity = min(10, floor((nrow(x_unique) - 1) / 3)),
  max_iter = 1000,
  check_duplicates = FALSE
)

tsne_df <- as.data.frame(tsne_result$Y) %>%
  setNames(c("tSNE1", "tSNE2")) %>%
  mutate(sample_id = rownames(x_unique)) %>%
  left_join(annot, by = "sample_id")

# Plot - confidence score
p_tsne_conf <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, colour = confidence_score, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_gradient(low = "#EF5350", high = "#66BB6A") +
  labs(
    title = "t-SNE of sinonasal tumour methylation profiles",
    subtitle = "20,000 most variable probes | Coloured by classifier confidence score",
    colour = "Confidence score"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_tsne_confidence.pdf"), p_tsne_conf, width = 9, height = 7)

# Plot - tumour purity
p_tsne_purity <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, colour = tumour_purity_pct, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_gradient(low = "#EF5350", high = "#66BB6A") +
  labs(
    title = "t-SNE of sinonasal tumour methylation profiles",
    subtitle = "20,000 most variable probes | Coloured by tumour purity (HiTIMED, HNSC reference)",
    colour = "Tumour purity (%)"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_tsne_tumour_purity.pdf"), p_tsne_purity, width = 9, height = 7)

# Plot - predicted class
p_tsne_class <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, colour = predicted_class, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  labs(
    title = "t-SNE of sinonasal tumour methylation profiles",
    subtitle = "20,000 most variable probes | Coloured by predicted methylation class",
    colour = "Predicted class"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_tsne_class.pdf"), p_tsne_class, width = 9, height = 7)

cat("\nRCSI t-SNE plots saved.\n")

# --- 8. Combined UMAP and t-SNE (RCSI + Jurmeister) --------------------------
# Probe selection uses variance from YOUR samples only across the intersection
# of probes common to both datasets.

common_probes <- intersect(rownames(betas_full), rownames(jur_betas_mat_cleaned))
cat("\nProbes in your data:", nrow(betas_full), "\n")
cat("Probes in Jurmeister data:", nrow(jur_betas_mat_cleaned), "\n")
cat("Common probes (intersection):", length(common_probes), "\n")

probe_vars_common <- apply(betas_full[common_probes, ], 1, var, na.rm = TRUE)
top_var_probes_common <- names(sort(probe_vars_common, decreasing = TRUE))[1:20000]
cat("Most variable common probes selected:", length(top_var_probes_common), "\n")

# Merge beta matrices
betas_cohort <- betas_full[top_var_probes_common, ]
betas_ref <- jur_betas_mat_cleaned[top_var_probes_common, ]
betas_merged <- cbind(betas_cohort, betas_ref)
cat("Merged matrix:", nrow(betas_merged), "probes x", ncol(betas_merged), "samples\n")
cat("Cohort samples:", ncol(betas_cohort), "\n")
cat("Jurmeister samples:", ncol(betas_ref), "\n")

x_merged <- t(betas_merged)

# Build combined annotation
yours_annot <- results %>%
  select(sample_id, predicted_class, confidence_score, confidence_tier, frac_na, qc_flag) %>%
  mutate(
    sample_id = trimws(sample_id),
    dataset = "RCSI"
  )

normalize_jur_class <- function(raw) {
  raw <- gsub("_", " ", raw)
  raw <- ifelse(raw == "APA", "PIT AD", raw)
  raw <- ifelse(raw == "NEC IDH2", "NEC-like IDH2", raw)
  raw <- ifelse(raw == "NEC SWI SNF", "NEC-like SMARCA4 ARID1A", raw)
  raw
}

jur_annot <- data.frame(
  sample_id = colnames(jur_betas_mat_cleaned),
  predicted_class = normalize_jur_class(gsub("[0-9]+$", "", colnames(jur_betas_mat_cleaned))),
  confidence_score = NA_real_,
  confidence_tier = NA_character_,
  frac_na = NA_real_,
  qc_flag = NA_character_,
  dataset = "Jurmeister",
  stringsAsFactors = FALSE
)

annot_combined <- bind_rows(yours_annot, jur_annot)

non_target_classes <- c("ATRT", "Others", "PDCA", "SCNEC")

rcsi_classes <- sort(unique(annot_combined$predicted_class[
  annot_combined$dataset == "RCSI" &
    !is.na(annot_combined$predicted_class) &
    !annot_combined$predicted_class %in% non_target_classes
]))

rcsi_pass_classes <- sort(unique(annot_combined$predicted_class[
  annot_combined$dataset == "RCSI" &
    annot_combined$qc_flag == "PASS" &
    !is.na(annot_combined$predicted_class) &
    !annot_combined$predicted_class %in% non_target_classes
]))

jurmeister_classes <- sort(unique(annot_combined$predicted_class[
  annot_combined$dataset == "Jurmeister" &
    !is.na(annot_combined$predicted_class) &
    !annot_combined$predicted_class %in% non_target_classes
]))

jurmeister_only_classes <- setdiff(jurmeister_classes, rcsi_classes)
jurmeister_only_pass_classes <- setdiff(jurmeister_classes, rcsi_pass_classes)
all_display_classes <- sort(union(rcsi_classes, jurmeister_classes))

class_palette <- setNames(
  scales::hue_pal()(length(all_display_classes)),
  all_display_classes
)

cat("\nRCSI classes represented in first legend:\n")
print(rcsi_classes)
cat("\nAdditional Jurmeister classes represented in second legend:\n")
print(jurmeister_only_classes)
cat("\nQC PASS RCSI classes represented in first legend:\n")
print(rcsi_pass_classes)
cat("\nAdditional Jurmeister classes in QC PASS plots:\n")
print(jurmeister_only_pass_classes)

# Impute NAs
x_imp_merged <- apply(x_merged, 2, function(probe) {
  na_idx <- is.na(probe)
  if (!any(na_idx)) return(probe)
  fill <- mean(probe, na.rm = TRUE)
  if (!is.finite(fill)) fill <- 0.5
  probe[na_idx] <- fill
  probe
})
cat("NAs after imputation:", sum(is.na(x_imp_merged)), "\n")

# --- 8a. Combined UMAP -------------------------------------------------------

set.seed(42)
umap_config_combined <- umap.defaults
umap_config_combined$n_neighbors <- 15
umap_config_combined$min_dist <- 0.1
umap_config_combined$metric <- "euclidean"

cat("\nRunning UMAP on merged dataset...\n")
umap_result_combined <- umap(x_imp_merged, config = umap_config_combined)

umap_combined_df <- as.data.frame(umap_result_combined$layout) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(sample_id = rownames(x_imp_merged)) %>%
  left_join(annot_combined, by = "sample_id")

# Plot 1a - methylation class (All)
p_umap_combined_class <- ggplot() +
  geom_point(
    data = umap_combined_df %>% filter(dataset == "Jurmeister", !predicted_class %in% non_target_classes),
    aes(x = UMAP1, y = UMAP2, colour = predicted_class),
    shape = 16, size = 2, alpha = 0.8
  ) +
  geom_point(
    data = umap_combined_df %>% filter(dataset == "RCSI", !predicted_class %in% non_target_classes),
    aes(x = UMAP1, y = UMAP2, fill = predicted_class),
    shape = 24, size = 2, colour = "black", stroke = 0.5
  ) +
  scale_fill_manual(values = class_palette, breaks = rcsi_classes, name = "RCSI predicted classes") +
  scale_colour_manual(values = class_palette, breaks = jurmeister_only_classes, name = "Additional Jurmeister classes") +
  geom_text_repel(
    data = filter(umap_combined_df, dataset == "RCSI"),
    aes(x = UMAP1, y = UMAP2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black", show.legend = FALSE
  ) +
  labs(
    title = "UMAP - methylation class overlay",
    subtitle = "Triangles = RCSI samples, Circles = Jurmeister reference"
  ) +
  theme_minimal(base_size = 10) +
  guides(
    fill = guide_legend(order = 1, ncol = 2, override.aes = list(shape = 24, colour = "black", size = 3)),
    colour = guide_legend(order = 2, ncol = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
  )

ggsave(file.path(out_dimred, "plot_umap_combined_class.pdf"), p_umap_combined_class, width = 12, height = 8)

# Plot 1b - methylation class (QC = PASS)
p_umap_combined_pass <- ggplot() +
  geom_point(
    data = umap_combined_df %>% filter(dataset == "Jurmeister", !predicted_class %in% non_target_classes),
    aes(x = UMAP1, y = UMAP2, colour = predicted_class),
    shape = 16, size = 2, alpha = 0.8
  ) +
  geom_point(
    data = umap_combined_df %>% filter(dataset == "RCSI", qc_flag == "PASS", !predicted_class %in% non_target_classes),
    aes(x = UMAP1, y = UMAP2, fill = predicted_class),
    shape = 24, size = 2, colour = "black", stroke = 0.5
  ) +
  scale_fill_manual(values = class_palette, breaks = rcsi_pass_classes, name = "RCSI predicted classes") +
  scale_colour_manual(values = class_palette, breaks = jurmeister_only_pass_classes, name = "Additional Jurmeister classes") +
  geom_text_repel(
    data = umap_combined_df %>% filter(dataset == "RCSI", qc_flag == "PASS"),
    aes(x = UMAP1, y = UMAP2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black", show.legend = FALSE
  ) +
  labs(
    title = "UMAP - methylation class overlay (QC PASS samples only)",
    subtitle = "Triangles = RCSI samples, Circles = Jurmeister reference"
  ) +
  theme_minimal(base_size = 10) +
  guides(
    fill = guide_legend(order = 1, ncol = 2, override.aes = list(shape = 24, colour = "black", size = 3)),
    colour = guide_legend(order = 2, ncol = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
  )

ggsave(file.path(out_dimred, "plot_umap_combined_class_pass.pdf"), p_umap_combined_pass, width = 12, height = 8)

# Plot 2 - QC status (RCSI samples only, Jurmeister greyed)
p_umap_combined_qc <- ggplot() +
  geom_point(
    data = filter(umap_combined_df, dataset == "Jurmeister"),
    aes(x = UMAP1, y = UMAP2),
    colour = "#D0D0D0", size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = filter(umap_combined_df, dataset == "RCSI"),
    aes(x = UMAP1, y = UMAP2, colour = qc_flag),
    size = 2
  ) +
  geom_text_repel(
    data = filter(umap_combined_df, dataset == "RCSI"),
    aes(x = UMAP1, y = UMAP2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black"
  ) +
  scale_colour_manual(values = c("PASS" = "#66BB6A", "FLAGGED" = "#EF5350")) +
  labs(
    title = "UMAP - RCSI cohort within Jurmeister reference space",
    subtitle = "20,000 most variable probes | Coloured by QC status",
    colour = "QC status"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_umap_combined_qc.pdf"), p_umap_combined_qc, width = 10, height = 8)

# Plot 3 - confidence score (RCSI samples only, Jurmeister greyed)
p_umap_combined_conf <- ggplot() +
  geom_point(
    data = filter(umap_combined_df, dataset == "Jurmeister"),
    aes(x = UMAP1, y = UMAP2),
    colour = "#D0D0D0", size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = filter(umap_combined_df, dataset == "RCSI"),
    aes(x = UMAP1, y = UMAP2, colour = confidence_score),
    size = 2
  ) +
  geom_text_repel(
    data = filter(umap_combined_df, dataset == "RCSI"),
    aes(x = UMAP1, y = UMAP2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black"
  ) +
  scale_colour_gradient(low = "#EF5350", high = "#66BB6A") +
  labs(
    title = "UMAP - RCSI cohort within Jurmeister reference space",
    subtitle = "20,000 most variable probes | Coloured by classifier confidence score",
    colour = "Confidence score"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_umap_combined_confidence.pdf"), p_umap_combined_conf, width = 10, height = 8)

cat("\nCombined UMAP plots saved.\n")

# --- 8b. Combined t-SNE ------------------------------------------------------

set.seed(42)
cat("\nRunning t-SNE on merged dataset...\n")

x_unique_merged <- unique(x_imp_merged)
n_unique <- nrow(x_unique_merged)
perp <- max(5, min(50, floor(sqrt(n_unique))))
cat("t-SNE perplexity:", perp, "(n unique samples:", n_unique, ")\n")

tsne_result_combined <- Rtsne(
  x_unique_merged,
  dims = 2,
  perplexity = perp,
  max_iter = 1000,
  check_duplicates = FALSE,
  pca = TRUE,
  pca_center = TRUE,
  pca_scale = FALSE,
  initial_dims = 50
)

tsne_combined_df <- as.data.frame(tsne_result_combined$Y) %>%
  setNames(c("tSNE1", "tSNE2")) %>%
  mutate(sample_id = rownames(x_unique_merged)) %>%
  left_join(annot_combined, by = "sample_id")

# Plot 1a - methylation class (All)
p_tsne_combined_class <- ggplot() +
  geom_point(
    data = tsne_combined_df %>% filter(dataset == "Jurmeister", !predicted_class %in% non_target_classes),
    aes(x = tSNE1, y = tSNE2, colour = predicted_class),
    shape = 16, size = 2, alpha = 0.8
  ) +
  geom_point(
    data = tsne_combined_df %>% filter(dataset == "RCSI", !predicted_class %in% non_target_classes),
    aes(x = tSNE1, y = tSNE2, fill = predicted_class),
    shape = 24, size = 2, colour = "black", stroke = 0.5
  ) +
  scale_fill_manual(values = class_palette, breaks = rcsi_classes, name = "RCSI predicted classes") +
  scale_colour_manual(values = class_palette, breaks = jurmeister_only_classes, name = "Additional Jurmeister classes") +
  geom_text_repel(
    data = filter(tsne_combined_df, dataset == "RCSI"),
    aes(x = tSNE1, y = tSNE2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black", show.legend = FALSE
  ) +
  labs(
    title = "t-SNE - methylation class overlay",
    subtitle = "Triangles = RCSI samples, Circles = Jurmeister reference"
  ) +
  theme_minimal(base_size = 10) +
  guides(
    fill = guide_legend(order = 1, ncol = 2, override.aes = list(shape = 24, colour = "black", size = 3)),
    colour = guide_legend(order = 2, ncol = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
  )

ggsave(file.path(out_dimred, "plot_tsne_combined_class.pdf"), p_tsne_combined_class, width = 12, height = 8)

# Plot 1b - methylation class (QC = PASS)
p_tsne_combined_pass <- ggplot() +
  geom_point(
    data = tsne_combined_df %>% filter(dataset == "Jurmeister", !predicted_class %in% non_target_classes),
    aes(x = tSNE1, y = tSNE2, colour = predicted_class),
    shape = 16, size = 2, alpha = 0.8
  ) +
  geom_point(
    data = tsne_combined_df %>% filter(dataset == "RCSI", qc_flag == "PASS", !predicted_class %in% non_target_classes),
    aes(x = tSNE1, y = tSNE2, fill = predicted_class),
    shape = 24, size = 2, colour = "black", stroke = 0.5
  ) +
  scale_fill_manual(values = class_palette, breaks = rcsi_pass_classes, name = "RCSI predicted classes") +
  scale_colour_manual(values = class_palette, breaks = jurmeister_only_pass_classes, name = "Additional Jurmeister classes") +
  geom_text_repel(
    data = tsne_combined_df %>% filter(dataset == "RCSI", qc_flag == "PASS"),
    aes(x = tSNE1, y = tSNE2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black", show.legend = FALSE
  ) +
  labs(
    title = "t-SNE - methylation class overlay (QC PASS samples only)",
    subtitle = "Triangles = RCSI samples, Circles = Jurmeister reference"
  ) +
  theme_minimal(base_size = 10) +
  guides(
    fill = guide_legend(order = 1, ncol = 2, override.aes = list(shape = 24, colour = "black", size = 3)),
    colour = guide_legend(order = 2, ncol = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
  )

ggsave(file.path(out_dimred, "plot_tsne_combined_class_pass.pdf"), p_tsne_combined_pass, width = 12, height = 8)

# Plot 2 - QC status (RCSI samples only, Jurmeister greyed)
p_tsne_combined_qc <- ggplot() +
  geom_point(
    data = filter(tsne_combined_df, dataset == "Jurmeister"),
    aes(x = tSNE1, y = tSNE2),
    colour = "#D0D0D0", size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = filter(tsne_combined_df, dataset == "RCSI"),
    aes(x = tSNE1, y = tSNE2, colour = qc_flag),
    size = 2
  ) +
  geom_text_repel(
    data = filter(tsne_combined_df, dataset == "RCSI"),
    aes(x = tSNE1, y = tSNE2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black"
  ) +
  scale_colour_manual(values = c("PASS" = "#66BB6A", "FLAGGED" = "#EF5350")) +
  labs(
    title = "t-SNE - RCSI cohort within Jurmeister reference space",
    subtitle = "20,000 most variable probes | Coloured by QC status",
    colour = "QC status"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_tsne_combined_qc.pdf"), p_tsne_combined_qc, width = 10, height = 8)

# Plot 3 - confidence score (RCSI samples only, Jurmeister greyed)
p_tsne_combined_conf <- ggplot() +
  geom_point(
    data = filter(tsne_combined_df, dataset == "Jurmeister"),
    aes(x = tSNE1, y = tSNE2),
    colour = "#D0D0D0", size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data = filter(tsne_combined_df, dataset == "RCSI"),
    aes(x = tSNE1, y = tSNE2, colour = confidence_score),
    size = 2
  ) +
  geom_text_repel(
    data = filter(tsne_combined_df, dataset == "RCSI"),
    aes(x = tSNE1, y = tSNE2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black"
  ) +
  scale_colour_gradient(low = "#EF5350", high = "#66BB6A") +
  labs(
    title = "t-SNE - RCSI cohort within Jurmeister reference space",
    subtitle = "20,000 most variable probes | Coloured by classifier confidence score",
    colour = "Confidence score"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dimred, "plot_tsne_combined_confidence.pdf"), p_tsne_combined_conf, width = 10, height = 8)

cat("\nCombined t-SNE plots saved.\n")

# Plot 1a-ii - methylation class (classified RCSI samples only, unknowns excluded)
rcsi_classified_classes <- setdiff(rcsi_classes, "Unknown")

p_tsne_combined_class_classified <- ggplot() +
  geom_point(
    data = tsne_combined_df %>% filter(dataset == "Jurmeister", !predicted_class %in% non_target_classes),
    aes(x = tSNE1, y = tSNE2, colour = predicted_class),
    shape = 16, size = 2, alpha = 0.8
  ) +
  geom_point(
    data = tsne_combined_df %>% filter(dataset == "RCSI", !predicted_class %in% non_target_classes, predicted_class != "Unknown"),
    aes(x = tSNE1, y = tSNE2, fill = predicted_class),
    shape = 24, size = 2, colour = "black", stroke = 0.5
  ) +
  scale_fill_manual(values = class_palette, breaks = rcsi_classified_classes, name = "RCSI predicted classes") +
  scale_colour_manual(values = class_palette, breaks = jurmeister_only_classes, name = "Additional Jurmeister classes") +
  geom_text_repel(
    data = tsne_combined_df %>% filter(dataset == "RCSI", predicted_class != "Unknown"),
    aes(x = tSNE1, y = tSNE2, label = sample_id),
    size = 2.5, max.overlaps = 20, colour = "black", show.legend = FALSE
  ) +
  labs(
    title = "t-SNE - methylation class overlay",
    subtitle = "Triangles = RCSI samples, Circles = Jurmeister reference"
  ) +
  theme_minimal(base_size = 10) +
  guides(
    fill = guide_legend(order = 1, ncol = 2, override.aes = list(shape = 24, colour = "black", size = 3)),
    colour = guide_legend(order = 2, ncol = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
  )

ggsave(file.path(out_dimred, "plot_tsne_combined_class_classified.pdf"), p_tsne_combined_class_classified, width = 12, height = 8)
# --- 9. Save dimensionality reduction outputs --------------------------------

write.csv(umap_df, file.path(out_dimred, "umap_coordinates.csv"), row.names = FALSE)
write.csv(tsne_df, file.path(out_dimred, "tsne_coordinates.csv"), row.names = FALSE)
write.csv(umap_combined_df, file.path(out_dimred, "umap_combined_coordinates.csv"), row.names = FALSE)
write.csv(tsne_combined_df, file.path(out_dimred, "tsne_combined_coordinates.csv"), row.names = FALSE)
write.csv(tumour_purity, file.path(out_purity, "hitimed_tumour_purity.csv"), row.names = FALSE)

cat("\nAll dimensionality reduction outputs saved.\n")

# --- 10. CNV profiling -------------------------------------------------------

cat("\nGenerating CNV profiles...\n")

# Build lookup from classifier results (already loaded)
cnv_lookup <- results %>%
  select(sample_id, predicted_class, qc_flag) %>%
  mutate(sample_id = trimws(sample_id))

# Load X/Y probe list once outside loop
sesameDataCache("EPIC.probeInfo")
probe_info <- sesameDataGet("EPIC.probeInfo")
xy_probes <- names(probe_info$mapped.probes.hg19)[
  as.character(GenomicRanges::seqnames(probe_info$mapped.probes.hg19)) %in% c("chrX", "chrY")
]
cat("X/Y probes to exclude:", length(xy_probes), "\n")

pdf(file.path(out_cnv, "all_samples_cnv.pdf"), width = 12, height = 4)

for (i in seq_len(nrow(sample_sheet))) {
  sid <- sample_sheet$sample_id[i]
  file_name <- sample_sheet$file_name[i]
  prefix <- idat_prefixes[basename(idat_prefixes) == file_name]
  
  if (length(prefix) == 0) {
    warning("No IDAT found for: ", sid)
    next
  }
  
  cat(sprintf("  [%d/%d] %s\n", i, nrow(sample_sheet), sid))
  
  meta <- cnv_lookup[cnv_lookup$sample_id == sid, ]
  pred_class <- if (nrow(meta) > 0 && !is.na(meta$predicted_class)) meta$predicted_class else "Unknown"
  qc_flag_label <- if (nrow(meta) > 0 && !is.na(meta$qc_flag)) meta$qc_flag else "Unknown"
  
  tryCatch({
    sdf <- readIDATpair(prefix)
    
    # Filter X/Y probes from SigDF
    sdf_probe_base <- sub("_.*$", "", sdf$Probe_ID)
    sdf <- sdf[!sdf_probe_base %in% xy_probes, ]
    
    seg <- cnSegmentation(sdf)
    
    p <- visualizeSegments(seg, to.plot = paste0("chr", 1:22)) +
      labs(
        title = sid,
        subtitle = paste0("Predicted class: ", pred_class, "  |  QC: ", qc_flag_label)
      ) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(
          size = 9,
          colour = ifelse(qc_flag_label == "PASS", "#2E7D32", "#C62828")
        )
      )
    
    print(p)
    
  }, error = function(e) {
    warning("CNV failed for ", sid, ": ", conditionMessage(e))
  })
  
  gc(verbose = FALSE)
}

dev.off()
cat("CNV profiles saved to:", out_cnv, "\n")
rm(sdf, seg)
gc()

# --- 11. Session info --------------------------------------------------------

sink(file.path(out_dimred, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Annotation and visualisation complete ===\n")
cat("Dimensionality reduction outputs saved to:", out_dimred, "\n")
cat("Tumour purity saved to:", out_purity, "\n")
cat("CNV profiles saved to:", out_cnv, "\n")