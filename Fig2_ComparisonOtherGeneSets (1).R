## =============================================================================
## ARI GeroScreen — brain aging signatures across three references
## Inputs (all in the working directory):
##   26_01_29_logcpm.csv                         log2 CPM, genes x samples
##   TableS3_Significant_genes_and_peaks_age.xlsx  Catching et al. 2026 (Zenodo 20834804)
##   PsychAD_EN/psychad_EN_age_signature.csv     from psychad_EN_age_signature.R
##   Lu et al. 2004 sets come from MSigDB C2 via msigdbr

## ---- 0. Setup ---------------------------------------------------------------
setwd("/restricted/projectnb/crem-trainees/Murphy_Lab/Elise/GeroScreen_Paper")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(fgsea)
  library(readxl)
})

logcpm_raw <- read.csv("26_01_29_logcpm.csv", check.names = FALSE,
                       stringsAsFactors = FALSE)

annot_cols  <- c("gene_key", "gene_id", "symbol", "feature_name",
                 "Length", "gene_biotype")
sample_cols <- setdiff(colnames(logcpm_raw), annot_cols)

out_dir <- "Fig2_brain_aging_refs"
dir.create(out_dir, showWarnings = FALSE)

base_font <- "sans"
fdr_cut   <- 0.05
set.seed(2026)

## ---- 1. Sample metadata and expression ------------------------------------------
group_levels <- c("CTRL", "RUT", "SIT", "GEN", "EBS", "CAL", "PI", "TRA", "RAP")
treatments   <- setdiff(group_levels, "CTRL")
group_cols <- setNames(c("grey55", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                         "#0072B2", "#D55E00", "#CC79A7", "#000000"),
                       group_levels)

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

## Ensembl -> symbol, for the Catching table (its gene column mixes both)
ens2sym <- setNames(logcpm_raw$symbol, sub("\\..*$", "", logcpm_raw$gene_id))

## ---- 2. Three brain aging references -----------------------------------------------
## Lu et al. 2004: bulk frontal cortex, ~30 donors (MSigDB C2)
c2 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2"),
               error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2"))

## Catching et al. 2026: excitatory neurons (ExN), snMultiome DLPFC
cat_tab <- read_excel("TableS3.Significant_genes_and_peaks_age.xlsx",
                      sheet = "Significant Genes") %>%
  mutate(sym = ifelse(grepl("^ENSG", gene), ens2sym[gene], gene)) %>%
  filter(cell_type == "ExN", !is.na(sym))

## PsychAD: excitatory-neuron pseudobulk age signature (derived in this project)
psy_tab <- read.csv("PsychAD_EN/psychad_EN_age_signature.csv", stringsAsFactors = FALSE)

sets <- list(
  `Lu 2004.up`       = unique(c2$gene_symbol[c2$gs_name == "LU_AGING_BRAIN_UP"]),
  `Lu 2004.down`     = unique(c2$gene_symbol[c2$gs_name == "LU_AGING_BRAIN_DN"]),
  `Catching 2026.up`   = unique(cat_tab$sym[cat_tab$slope > 0]),
  `Catching 2026.down` = unique(cat_tab$sym[cat_tab$slope < 0]),
  `PsychAD.up`       = unique(psy_tab$gene[psy_tab$direction == "up"]),
  `PsychAD.down`     = unique(psy_tab$gene[psy_tab$direction == "down"])
)
sets <- lapply(sets, function(g) intersect(g, rownames(expr)))
print(sapply(sets, length))

## ---- 3. limma, then fGSEA on each arm ----------------------------------------------
design <- model.matrix(~ 0 + group + line, data = meta)
colnames(design) <- make.names(sub("^group", "", colnames(design)))
contr <- makeContrasts(contrasts = paste0(make.names(treatments), " - CTRL"),
                       levels = design)
colnames(contr) <- treatments
fit <- eBayes(contrasts.fit(lmFit(expr, design), contr), trend = TRUE, robust = TRUE)

gsea <- bind_rows(lapply(treatments, function(tx) {
  res <- fgsea(sets, sort(fit$t[, tx], decreasing = TRUE),
               minSize = 5, maxSize = 1000)
  data.frame(treatment = tx, set = res$pathway, NES = res$NES,
             pval = res$pval, size = res$size)
})) %>%
  group_by(set) %>%                                   # BH within each arm
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  separate(set, into = c("reference", "arm"), sep = "\\.", remove = FALSE) %>%
  mutate(reference = factor(reference, levels = c("Lu 2004", "Catching 2026", "PsychAD")),
         arm       = factor(ifelse(arm == "up", "Genes up\nwith age",
                                   "Genes down\nwith age"),
                            levels = c("Genes up\nwith age", "Genes down\nwith age")),
         compound  = factor(treatment, levels = rev(treatments)),
         sig       = padj < fdr_cut,
         star      = as.character(cut(padj, c(-Inf, 0.001, 0.01, 0.05, Inf),
                                      labels = c("***", "**", "*", ""), right = FALSE)))
write.csv(gsea %>% select(treatment, reference, arm = set, NES, pval, padj, size),
          file.path(out_dir, "brain_aging_NES_all_refs.csv"), row.names = FALSE)

## Gene counts per arm, for the axis labels
set_n <- sapply(sets, length)
ref_lab <- sprintf("%s\n(%d up / %d down)",
                   levels(gsea$reference),
                   set_n[paste0(levels(gsea$reference), ".up")],
                   set_n[paste0(levels(gsea$reference), ".down")])
names(ref_lab) <- levels(gsea$reference)

## ---- 4. Figure ------------------------------------------------------------------------
nes_lim <- ceiling(max(abs(gsea$NES)) * 2) / 2

## Compound labels are coloured to match every other panel. Yellow and grey read
## poorly as text, so those get a darker version of the same hue.
label_cols <- group_cols
label_cols["EBS"] <- "#B59A00"   # darker ebselen yellow
theme_fig <- theme_bw(base_size = 10, base_family = base_font) +
  theme(panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
        strip.background  = element_rect(fill = "grey92"),
        strip.text        = element_text(face = "bold", size = 9),
        axis.text.y       = element_text(face = "bold",
                                         colour = label_cols[levels(gsea$compound)]),
        axis.text.x       = element_text(colour = "black"),
        legend.key.size   = unit(0.35, "cm"),
        plot.title        = element_text(face = "bold", size = 11),
        plot.subtitle     = element_text(size = 8),
        plot.title.position = "plot")

p_nes <- ggplot(gsea, aes(arm, compound)) +
  geom_point(aes(size = abs(NES), fill = NES, colour = sig), shape = 21, stroke = 0.7) +
  geom_text(aes(label = star), size = 3, vjust = 0.72, family = base_font) +
  facet_grid(~ reference, labeller = labeller(reference = ref_lab)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-nes_lim, nes_lim), name = "NES") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"), guide = "none") +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Brain aging signatures across three references",
       subtitle = paste("fGSEA vs vehicle (limma, blocked on line); black outline and stars = BH adj. p < 0.05.",
                        "Blue on \"up with age\" and red on \"down with age\" are younger-like.",
                        sep = "\n")) +
  theme_fig

p_nes

## ---- 5. Save (cairo for PNG and PDF) ---------------------------------------------
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = 600, device = png, type = "cairo")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf)
}
save_fig(p_nes, "Fig2_brain_aging_NES", 7, 4.4)

## ---- 6. Session info ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs written to: ", normalizePath(out_dir))