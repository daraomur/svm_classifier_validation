# =============================================================================
# Sinonasal Methylation Classifier — Enrichment Analysis (Sensitivity)
# SPELCASTER WP1b | RCSI
# Author: Dara O'Murchu
# Date: May 2026
# =============================================================================

# --- 0. Libraries ------------------------------------------------------------

library(limma)
library(missMethyl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(ggplot2)

# --- 1. Paths ----------------------------------------------------------------

base_dir <- "~/Desktop/Research Project/My Project/Dara_Validation"

out_dm <- file.path(base_dir, "results/enrichment_analysis/n35")
dir.create(out_dm, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load inputs ----------------------------------------------------------

cat("Loading inputs...\n")
betas   <- readRDS(file.path(base_dir, "results/beta_matrices/batch01_betas_QCDPB_lifted.rds"))
results <- read.csv(file.path(base_dir, "results/classifier/classifier_results_with_qc.csv"))

cat("Beta matrix:", nrow(betas), "probes x", ncol(betas), "samples\n")

# --- 3. Exclude poor-quality samples -----------------------------------------

exclude_samples <- c(
  "SPL-DU-00052",
  "SPL-DU-00053 T",
  "SPL-DU-00056"
)

cat("Excluding", length(exclude_samples), "poor-quality samples.\n")

results <- results %>%
  filter(!trimws(sample_id) %in% exclude_samples)

cat("Samples retained:", nrow(results), "\n")

# --- 4. Prepare sample groups ------------------------------------------------

annot <- results %>%
  filter(sample_id %in% colnames(betas)) %>%
  mutate(sample_id = trimws(sample_id))

betas <- betas[, annot$sample_id]

group <- factor(annot$qc_flag, levels = c("PASS", "FLAGGED"))
cat("PASS samples:", sum(group == "PASS"), "\n")
cat("FLAGGED samples:", sum(group == "FLAGGED"), "\n")

# --- 5. Filter and impute ----------------------------------------------------

na_frac    <- rowMeans(is.na(betas))
betas_filt <- betas[na_frac <= 0.2, ]
cat("Probes after NA filtering:", nrow(betas_filt), "\n")

betas_filt <- t(apply(betas_filt, 1, function(x) {
  nas <- is.na(x)
  if (any(nas)) x[nas] <- mean(x, na.rm = TRUE)
  x
}))

betas_clip <- pmax(pmin(betas_filt, 0.99), 0.01)
mvals      <- log2(betas_clip / (1 - betas_clip))
cat("M-value matrix:", nrow(mvals), "probes x", ncol(mvals), "samples\n")

# --- 6. Differential methylation with limma ----------------------------------

cat("\nRunning limma differential methylation analysis...\n")

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "FLAGGED_vs_PASS")

fit <- lmFit(mvals, design)
fit <- eBayes(fit)

dm_results <- topTable(fit,
                       coef    = "FLAGGED_vs_PASS",
                       number  = Inf,
                       sort.by = "p")

dm_results <- dm_results %>%
  rownames_to_column("CpG") %>%
  mutate(
    direction   = ifelse(logFC > 0,
                         "Hypermethylated_in_FLAGGED",
                         "Hypomethylated_in_FLAGGED"),
    significant = P.Value < 0.001 & abs(logFC) > 0.5
  )

write.csv(dm_results,
          file.path(out_dm, "dm_results_flagged_vs_pass.csv"),
          row.names = FALSE)

cat("Total CpGs tested:", nrow(dm_results), "\n")
cat("Significant (FDR < 0.05, |logFC| > 0.5):",
    sum(dm_results$adj.P.Val < 0.05 & abs(dm_results$logFC) > 0.5), "\n")
cat("Nominally significant (uncorrected p < 0.001, |logFC| > 0.5):",
    sum(dm_results$significant), "\n")
cat("Direction breakdown:\n")
print(table(dm_results$direction[dm_results$significant]))

# --- 7. Volcano plot ---------------------------------------------------------

fdr_threshold_pvalue <- max(dm_results$P.Value[dm_results$adj.P.Val < 0.05])

p_volcano <- ggplot(dm_results,
                    aes(x = logFC,
                        y = -log10(P.Value),
                        colour = adj.P.Val < 0.05 & abs(logFC) > 0.5)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#E53935"),
                      labels = c("Not significant", "FDR < 0.05, |logFC| > 0.5")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(fdr_threshold_pvalue), linetype = "dashed", colour = "grey40") +
  labs(title    = "Differential methylation: FLAGGED vs PASS (n=31, 7 excluded)",
       subtitle = "Dashed lines = FDR < 0.05 (Benjamini-Hochberg) threshold and |logFC| > 0.5",
       x        = "Log2 fold change (M-value)",
       y        = "-log10(p-value)",
       colour   = NULL) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dm, "plot_volcano.pdf"), p_volcano, width = 8, height = 6)
cat("Volcano plot saved.\n")

# --- 8. Exploratory enrichment on nominally significant CpGs -----------------

sig_cpgs <- dm_results %>%
  filter(significant)

cat("\nNominally significant CpGs:", nrow(sig_cpgs), "\n")
cat("Direction breakdown:\n")
print(table(sig_cpgs$direction))

if (nrow(sig_cpgs) > 0) {
  
  cat("\nRunning gometh GO enrichment...\n")
  
  gometh_go <- gometh(
    sig.cpg    = sig_cpgs$CpG,
    all.cpg    = dm_results$CpG,
    collection = "GO",
    array.type = "EPIC",
    plot.bias  = FALSE
  )
  
  gometh_go_sig <- gometh_go %>%
    arrange(FDR) %>%
    filter(FDR < 0.05)
  
  cat("Significant GO terms (FDR < 0.05):", nrow(gometh_go_sig), "\n")
  
  write.csv(gometh_go_sig,
            file.path(out_dm, "gometh_go_results.csv"),
            row.names = FALSE)
  
  cat("\nRunning gometh KEGG enrichment...\n")
  
  gometh_kegg <- gometh(
    sig.cpg    = sig_cpgs$CpG,
    all.cpg    = dm_results$CpG,
    collection = "KEGG",
    array.type = "EPIC",
    plot.bias  = FALSE
  )
  
  gometh_kegg_sig <- gometh_kegg %>%
    arrange(FDR) %>%
    filter(FDR < 0.05)
  
  cat("Significant KEGG pathways (FDR < 0.05):", nrow(gometh_kegg_sig), "\n")
  
  write.csv(gometh_kegg_sig,
            file.path(out_dm, "gometh_kegg_results.csv"),
            row.names = FALSE)
  
  # --- 9. clusterProfiler ----------------------------------------------------
  
  cat("\nRunning clusterProfiler enrichment...\n")
  
  cpg_genes <- getMappedEntrezIDs(
    sig.cpg    = sig_cpgs$CpG,
    all.cpg    = dm_results$CpG,
    array.type = "EPIC"
  )
  
  sig_entrez <- unique(cpg_genes$sig.eg)
  all_entrez <- unique(cpg_genes$universe)
  
  cat("Significant genes (Entrez):", length(sig_entrez), "\n")
  cat("Background genes (Entrez):", length(all_entrez), "\n")
  
  ego <- enrichGO(
    gene          = sig_entrez,
    universe      = all_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  cat("Significant GO BP terms:", nrow(as.data.frame(ego)), "\n")
  
  if (nrow(as.data.frame(ego)) > 0) {
    
    write.csv(as.data.frame(ego),
              file.path(out_dm, "clusterProfiler_go_bp.csv"),
              row.names = FALSE)
    
    p_dot <- dotplot(ego, showCategory = 20) +
      labs(title    = "GO Biological Process enrichment",
           subtitle = "Nominally significant DMCs (uncorrected p < 0.001): FLAGGED vs PASS (n=31)") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(out_dm, "plot_go_dotplot.pdf"), p_dot, width = 10, height = 8)
    cat("GO dotplot saved.\n")
    
    p_emap <- emapplot(enrichplot::pairwise_termsim(ego), showCategory = 30) +
      labs(title = "GO term enrichment network (n=31)") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(out_dm, "plot_go_emapplot.pdf"), p_emap, width = 12, height = 10)
    cat("GO enrichment map saved.\n")
  }
  
  ekegg <- enrichKEGG(
    gene          = sig_entrez,
    universe      = all_entrez,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  cat("Significant KEGG pathways:", nrow(as.data.frame(ekegg)), "\n")
  
  if (nrow(as.data.frame(ekegg)) > 0) {
    
    write.csv(as.data.frame(ekegg),
              file.path(out_dm, "clusterProfiler_kegg.csv"),
              row.names = FALSE)
    
    p_kegg <- dotplot(ekegg, showCategory = 20) +
      labs(title    = "KEGG pathway enrichment (n=31)",
           subtitle = "Nominally significant DMCs (uncorrected p < 0.001)") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(out_dm, "plot_kegg_dotplot.pdf"), p_kegg, width = 10, height = 8)
    cat("KEGG dotplot saved.\n")
  }
  
} else {
  cat("No nominally significant CpGs — enrichment analysis skipped.\n")
}

# --- 10. Session info --------------------------------------------------------

sink(file.path(out_dm, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Sensitivity enrichment analysis complete ===\n")
cat("Outputs saved to:", out_dm, "\n")
