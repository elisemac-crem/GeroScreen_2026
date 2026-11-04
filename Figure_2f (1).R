## =============================================================================
## ARI GeroScreen — Figure 2f: age-up genes suppressed by trametinib
## Inputs (all in the working directory):
##   26_01_29_logcpm.csv                      log2 CPM, genes x samples
##   PsychAD_EN/psychad_EN_age_signature.csv  from psychad_EN_age_signature.R
##   PsychAD_EN/psychad_EN_age_genes_all.csv  (same script; all genes tested)

## ---- 0. Setup ---------------------------------------------------------------
setwd("/restricted/projectnb/crem-trainees/Murphy_Lab/Elise/GeroScreen_Paper")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(fgsea)
})

logcpm_raw <- read.csv("26_01_29_logcpm.csv", check.names = FALSE,
                       stringsAsFactors = FALSE)

annot_cols  <- c("gene_key", "gene_id", "symbol", "feature_name",
                 "Length", "gene_biotype")
sample_cols <- setdiff(colnames(logcpm_raw), annot_cols)

out_dir <- "Fig2f_leading_edge"
dir.create(out_dir, showWarnings = FALSE)

## ---- Figure settings (change these to fit the layout) ------------------------
fig_width  <- 4.0    # inches — narrow, to sit beside panel d/e
fig_height <- 5.0    # inches — tall
n_le       <- 20     # genes shown; raise for a taller panel, lower for shorter
gene_size  <- 7      # gene label size (pt)

base_font <- "sans"
tra_col   <- "#CC79A7"   # trametinib colour used in every panel
set.seed(2026)

## ---- 1. Sample metadata and expression ------------------------------------------
group_levels <- c("CTRL", "RUT", "SIT", "GEN", "EBS", "CAL", "PI", "TRA", "RAP")
treatments   <- setdiff(group_levels, "CTRL")

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

## ---- 2. PsychAD "up with age" genes ------------------------------------------------
psy_tab <- read.csv("PsychAD_EN/psychad_EN_age_signature.csv", stringsAsFactors = FALSE)
psy_all <- read.csv("PsychAD_EN/psychad_EN_age_genes_all.csv", stringsAsFactors = FALSE)
age_up  <- intersect(psy_tab$gene[psy_tab$direction == "up"], rownames(expr))
message(sprintf("PsychAD age-up genes expressed in these neurons: %d", length(age_up)))

## ---- 3. limma: trametinib vs vehicle, pooled and per line --------------------------
design <- model.matrix(~ 0 + group + line, data = meta)
colnames(design) <- make.names(sub("^group", "", colnames(design)))
contr <- makeContrasts(contrasts = paste0(make.names(treatments), " - CTRL"),
                       levels = design)
colnames(contr) <- treatments
fit <- eBayes(contrasts.fit(lmFit(expr, design), contr), trend = TRUE, robust = TRUE)

meta$line_group <- factor(paste(meta$line, meta$group, sep = "_"))
design_l <- model.matrix(~ 0 + line_group, data = meta)
colnames(design_l) <- sub("^line_group", "", colnames(design_l))
line_contrasts <- unlist(lapply(levels(meta$line), function(l)
  paste0(l, "_", treatments, " - ", l, "_CTRL")))
contr_l <- makeContrasts(contrasts = line_contrasts, levels = design_l)
colnames(contr_l) <- sub(" - .*$", "", line_contrasts)
fit_l <- eBayes(contrasts.fit(lmFit(expr, design_l), contr_l),
                trend = TRUE, robust = TRUE)

## ---- 4. Leading edge of the age-up set for trametinib -------------------------------
gsea_tra <- fgsea(list(age_up = age_up), sort(fit$t[, "TRA"], decreasing = TRUE),
                  minSize = 5, maxSize = 1000)
le_genes <- unique(unlist(gsea_tra$leadingEdge))
message(sprintf("Leading edge: %d genes (NES %.2f, p = %.3g)",
                length(le_genes), gsea_tra$NES, gsea_tra$pval))

le_df <- data.frame(gene = le_genes, stringsAsFactors = FALSE) %>%
  mutate(pooled     = fit$coefficients[gene, "TRA"],
         BU6        = fit_l$coefficients[gene, "BU6_TRA"],
         BU8        = fit_l$coefficients[gene, "BU8_TRA"],
         adj.P.Val  = fit$p.value[gene, "TRA"],
         age_effect = psy_all$logFC_per_decade[match(gene, psy_all$gene)]) %>%
  arrange(pooled) %>%
  slice_head(n = n_le) %>%
  mutate(gene = factor(gene, levels = rev(gene)))   # strongest at the top
write.csv(le_df %>% mutate(gene = as.character(gene)),
          file.path(out_dir, "Fig2f_leading_edge_genes.csv"), row.names = FALSE)

## ---- 5. Figure ----------------------------------------------------------------------
theme_fig <- theme_bw(base_size = 10, base_family = base_font) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.x        = element_text(colour = "black"),
        axis.text.y        = element_text(face = "italic", size = gene_size,
                                          colour = "black"),
        axis.ticks.y       = element_blank(),
        legend.key.size    = unit(0.32, "cm"),
        legend.position    = "bottom",
        legend.margin      = margin(t = -4),
        plot.title         = element_text(face = "bold", size = 10),
        plot.subtitle      = element_text(size = 7.5),
        plot.title.position = "plot")
p_f <- ggplot(le_df, aes(pooled, gene)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_col(fill = tra_col, width = 0.7) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.02))) +
  labs(x = expression(Trametinib~log[2]*FC~vs~vehicle), y = NULL,
       title = "Age-up genes suppressed by trametinib",
       subtitle = sprintf("Leading edge of the PsychAD \"up with age\" set (top %d);\npooled estimate, blocked on cell line",
                          n_le)) +
  theme_fig

p_f
## ---- 6. Save (cairo for PNG and PDF) ---------------------------------------------
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = 600, device = png, type = "cairo")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf)
}
save_fig(p_f, "Fig2f_TRA_leading_edge", fig_width, fig_height)

## ---- 7. Session info ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs written to: ", normalizePath(out_dir))