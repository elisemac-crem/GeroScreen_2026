## =============================================================================
## ARI GeroScreen — Supplementary figure: validating the brain aging references

## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
setwd("/restricted/projectnb/crem-trainees/Murphy_Lab/Elise/GeroScreen_Paper")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(readxl)
})

logcpm_raw <- read.csv("26_01_29_logcpm.csv", check.names = FALSE,
                       stringsAsFactors = FALSE)
tage_raw   <- read.csv("tAge_GeroScreen_vehicle_cell_line.csv",
                       stringsAsFactors = FALSE, encoding = "latin1")

annot_cols  <- c("gene_key", "gene_id", "symbol", "feature_name",
                 "Length", "gene_biotype")
sample_cols <- setdiff(colnames(logcpm_raw), annot_cols)

out_dir <- "FigS_reference_validation"
dir.create(out_dir, showWarnings = FALSE)

base_font <- "sans"
fdr_cut   <- 0.05
set.seed(2026)

## ---- 1. Sample metadata, colours, expression -----------------------------------
group_levels <- c("CTRL", "RUT", "SIT", "GEN", "EBS", "CAL", "PI", "TRA", "RAP")
treatments   <- setdiff(group_levels, "CTRL")
group_cols <- setNames(c("grey55", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                         "#0072B2", "#D55E00", "#CC79A7", "#000000"),
                       group_levels)
label_cols <- group_cols
label_cols["EBS"] <- "#B59A00"     # darker yellow, readable as text

meta <- data.frame(sample = sample_cols, stringsAsFactors = FALSE) %>%
  mutate(line  = factor(sub("_.*$", "", sample)),
         group = factor(sub("[0-9]+$", "", sub("^[^_]+_", "", sample)),
                        levels = group_levels))
stopifnot(!any(is.na(meta$group)))

expr_df <- logcpm_raw %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(row_mean = rowMeans(across(all_of(sample_cols)))) %>%
  group_by(symbol) %>%
  slice_max(row_mean, n = 1, with_ties = FALSE) %>%
  ungroup()
expr <- as.matrix(expr_df[, meta$sample])
rownames(expr) <- expr_df$symbol
ens2sym <- setNames(logcpm_raw$symbol, sub("\\..*$", "", logcpm_raw$gene_id))

## ---- 2. The three references ------------------------------------------------------
c2 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2"),
               error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2"))
lu_tab <- bind_rows(
  data.frame(gene = unique(c2$gene_symbol[c2$gs_name == "LU_AGING_BRAIN_UP"]), dir = 1),
  data.frame(gene = unique(c2$gene_symbol[c2$gs_name == "LU_AGING_BRAIN_DN"]), dir = -1))

cat_tab <- read_excel("TableS3.Significant_genes_and_peaks_age.xlsx",
                      sheet = "Significant Genes") %>%
  mutate(sym = ifelse(grepl("^ENSG", gene), ens2sym[gene], gene)) %>%
  filter(cell_type == "ExN", !is.na(sym)) %>%
  group_by(sym) %>% slice_min(p_value, n = 1, with_ties = FALSE) %>% ungroup()

psy_tab <- read.csv("PsychAD_EN/psychad_EN_age_signature.csv", stringsAsFactors = FALSE)
psy_all <- read.csv("PsychAD_EN/psychad_EN_age_genes_all.csv", stringsAsFactors = FALSE)

sets <- list(
  `Lu 2004 up`       = lu_tab$gene[lu_tab$dir > 0],
  `Lu 2004 down`     = lu_tab$gene[lu_tab$dir < 0],
  `Catching up`      = cat_tab$sym[cat_tab$slope > 0],
  `Catching down`    = cat_tab$sym[cat_tab$slope < 0],
  `PsychAD up`       = psy_tab$gene[psy_tab$direction == "up"],
  `PsychAD down`     = psy_tab$gene[psy_tab$direction == "down"]
)

## ---- 3. Shared theme and save helper ------------------------------------------------
theme_fig <- theme_bw(base_size = 10, base_family = base_font) +
  theme(panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
        axis.text         = element_text(colour = "black"),
        legend.key.size   = unit(0.32, "cm"),
        plot.title        = element_text(face = "bold", size = 10),
        plot.subtitle     = element_text(size = 7.5),
        plot.title.position = "plot")

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = 600, device = png, type = "cairo")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf)
}

## ---- 4. Panel a: PsychAD vs Catching, gene by gene --------------------------------
## Independent cohorts, same cell type. Agreement here is what justifies using
## the PsychAD signature derived in this study.
rep_df <- cat_tab %>%
  inner_join(psy_all, by = c("sym" = "gene")) %>%
  transmute(gene = sym, catching = slope, psychad = logFC_per_decade,
            agree = sign(catching) == sign(psychad))
n_agree <- sum(rep_df$agree)
bt      <- binom.test(n_agree, nrow(rep_df))
rho_rep <- cor(rep_df$catching, rep_df$psychad, method = "spearman")
write.csv(rep_df, file.path(out_dir, "PsychAD_vs_Catching_genes.csv"), row.names = FALSE)

p_a <- ggplot(rep_df, aes(catching, psychad)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "grey85", alpha = 0.45) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey85", alpha = 0.45) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_point(size = 1.4, alpha = 0.75, colour = "#0072B2") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.35, size = 2.6,
           family = base_font,
           label = sprintf("%d / %d agree in direction\n(p = %.0e)\nrho = %.2f",
                           n_agree, nrow(rep_df), bt$p.value, rho_rep)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks = c(-30, -10, -3, 0, 3, 10)) +
  labs(x = "Catching 2026 age slope (pseudo-log scale)",
       y = expression(PsychAD~age~effect~(log[2]*FC/decade)),
       title = "Two neuron references agree",
       subtitle = "Age genes in excitatory neurons,\nshared between cohorts") +
  theme_fig

save_fig(p_a, "FigS_a_PsychAD_vs_Catching", 3.6, 3.6)

## ---- 5. Panel b: how little the references share ------------------------------------
ov <- expand.grid(set1 = names(sets), set2 = names(sets),
                  stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(n_shared = length(intersect(sets[[set1]], sets[[set2]])),
         n_set1   = length(sets[[set1]]),
         pct      = 100 * n_shared / n_set1) %>%
  ungroup() %>%
  mutate(set1 = factor(set1, levels = rev(names(sets))),
         set2 = factor(set2, levels = names(sets)),
         diag = as.character(set1) == as.character(set2))
write.csv(ov %>% select(set1, set2, n_shared, n_set1, pct),
          file.path(out_dir, "reference_overlap.csv"), row.names = FALSE)

p_b <- ggplot(ov, aes(set2, set1)) +
  geom_tile(aes(fill = ifelse(diag, NA, pct)), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n_shared,
                fontface = ifelse(diag, "bold", "plain")), size = 2.6,
            family = base_font) +
  scale_fill_gradient(low = "white", high = "#B2182B", na.value = "grey90",
                      name = "% of row\nshared") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "The three references barely overlap",
       subtitle = "Genes shared between signature arms (bold = set size)") +
  theme_fig +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 0, size = 7.5),
        axis.text.y = element_text(size = 7.5))

save_fig(p_b, "FigS_b_reference_overlap", 4.2, 3.4)

## ---- 6. Panel c: transcriptional response size vs tAge shift -----------------------
design <- model.matrix(~ 0 + group + line, data = meta)
colnames(design) <- make.names(sub("^group", "", colnames(design)))
contr <- makeContrasts(contrasts = paste0(make.names(treatments), " - CTRL"),
                       levels = design)
colnames(contr) <- treatments
fit <- eBayes(contrasts.fit(lmFit(expr, design), contr), trend = TRUE, robust = TRUE)

deg_n <- sapply(treatments, function(tx)
  sum(topTable(fit, coef = tx, number = Inf, sort.by = "none")$adj.P.Val < fdr_cut))

tage <- tage_raw %>%
  mutate(line  = factor(cell_line),
         group = factor(sub("[0-9]+$", "", sub("^[^_]+_", "", sample_name)),
                        levels = group_levels))
tage_est <- function(col) {
  m <- summary(lm(tage[[col]] ~ group + line, data = tage))$coefficients
  setNames(m[paste0("group", treatments), "Estimate"], treatments)
}

resp_df <- data.frame(treatment = treatments,
                      n_deg  = deg_n,
                      chrono = tage_est("tAge_Chronoage"),
                      mort   = tage_est("tAge_Mortality")) %>%
  mutate(group = factor(treatment, levels = group_levels))
write.csv(resp_df, file.path(out_dir, "response_vs_tAge.csv"), row.names = FALSE)

rho_c <- cor(resp_df$n_deg, abs(resp_df$chrono), method = "spearman")

lab_layer_c <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = treatment, colour = group), size = 2.6,
                           fontface = "bold", family = base_font,
                           box.padding = 0.4, min.segment.length = 0,
                           segment.size = 0.25, seed = 1)
} else {
  geom_text(aes(label = treatment, colour = group), size = 2.6, vjust = -1.2,
            fontface = "bold", family = base_font)
}

p_c <- ggplot(resp_df, aes(n_deg + 1, chrono)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(fill = group), shape = 21, size = 3.2, stroke = 0.6,
             colour = "black") +
  lab_layer_c +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_colour_manual(values = label_cols, guide = "none") +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                labels = c("0", "10", "100", "1,000", "10,000")) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.35, size = 2.6,
           family = base_font,
           label = sprintf("rho = %.2f (n = 8)", rho_c)) +
  labs(x = "Differentially expressed genes vs vehicle",
       y = "Chronological tAge shift (years)",
       title = "Clock shifts without a transcriptional response",
       subtitle = sprintf("DEGs at FDR < %.2f; sitagliptin and calpeptin move the clock\nwith almost no detectable response", fdr_cut)) +
  theme_fig

save_fig(p_c, "FigS_c_response_vs_tAge", 3.8, 3.6)

## ---- 7. Panel d: technical-replicate sensitivity -----------------------------------
## The two replicates per line are technical, so treating all four samples per
## compound as independent overstates confidence. Here the same tAge models are
## refitted two ways: all samples (n = 4 per compound) and replicate means
## (n = 2 lines per compound). Nothing survives correction the conservative way,
## which is why the tAge results are framed as exploratory.
tage_avg <- tage %>%
  group_by(group, line) %>%
  summarise(tAge_Chronoage = mean(tAge_Chronoage),
            tAge_Mortality = mean(tAge_Mortality), .groups = "drop")

fit_tage <- function(dat, col, model_label) {
  m  <- summary(lm(dat[[col]] ~ group + line, data = dat))$coefficients
  rn <- paste0("group", treatments)
  data.frame(treatment = treatments,
             estimate  = m[rn, "Estimate"],
             P.Value   = m[rn, "Pr(>|t|)"],
             readout   = ifelse(col == "tAge_Chronoage",
                                "Chronological tAge", "Mortality tAge"),
             model     = model_label)
}

sens_df <- bind_rows(
  fit_tage(tage,     "tAge_Chronoage", "All samples\n(n = 4)"),
  fit_tage(tage,     "tAge_Mortality", "All samples\n(n = 4)"),
  fit_tage(tage_avg, "tAge_Chronoage", "Replicate means\n(n = 2 lines)"),
  fit_tage(tage_avg, "tAge_Mortality", "Replicate means\n(n = 2 lines)")
) %>%
  group_by(readout, model) %>%
  mutate(adj.P.Val = p.adjust(P.Value, method = "BH")) %>%
  ungroup() %>%
  mutate(model   = factor(model, levels = c("All samples\n(n = 4)",
                                            "Replicate means\n(n = 2 lines)")),
         readout = factor(readout, levels = c("Chronological tAge", "Mortality tAge")),
         group   = factor(treatment, levels = rev(treatments)),
         sig     = adj.P.Val < fdr_cut)
write.csv(sens_df %>% mutate(model = gsub("\n", " ", model)),
          file.path(out_dir, "technical_replicate_sensitivity.csv"), row.names = FALSE)

p_d <- ggplot(sens_df, aes(adj.P.Val, group)) +
  geom_vline(xintercept = fdr_cut, linetype = "dashed", colour = "#B2182B",
             linewidth = 0.4) +
  geom_line(aes(group = group), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(fill = group, shape = model), size = 2.4, stroke = 0.5,
             colour = "black") +
  facet_wrap(~ readout, nrow = 1) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_x_log10(breaks = c(0.01, 0.05, 0.2, 1),
                labels = c("0.01", "0.05", "0.2", "1")) +
  labs(x = "BH-adjusted p vs vehicle", y = NULL,
       title = "Technical replicates inflate confidence",
       subtitle = "Same models refitted on replicate means; dashed line = 0.05") +
  theme_fig +
  theme(axis.text.y = element_text(face = "bold",
                                   colour = label_cols[levels(sens_df$group)]),
        legend.position = "bottom", legend.margin = margin(t = -4),
        panel.grid.major.y = element_blank())

save_fig(p_d, "FigS_d_technical_replicates", 4.4, 3.6)

## ---- 8. Combined layout ---------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- (p_a | p_b) / (p_c | p_d) +
    plot_layout(widths = c(1, 1.15)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 11, family = base_font))
  save_fig(combined, "FigS_reference_validation", 7.8, 7.2)
}

## ---- 9. Session info ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs written to: ", normalizePath(out_dir))