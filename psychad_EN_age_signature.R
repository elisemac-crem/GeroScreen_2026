## =============================================================================
## ARI GeroScreen — PsychAD excitatory-neuron aging signature
## Input : PsychAD_EN/psychad_EN_pseudobulk_counts.csv.gz  (from psychad_pseudobulk.py)
##         PsychAD_EN/psychad_EN_donor_meta.csv
## Output: PsychAD_EN/psychad_EN_age_genes_all.csv   — every gene tested
##         PsychAD_EN/psychad_EN_age_signature.csv    — significant genes, up/down
##         PsychAD_EN/psychad_EN_age_QC.png           — donor age distribution + MDS
##
## Model: adult donors only (age >= min_age), voom pseudobulk,
##        expression ~ age (per decade) + sex + log(cells per donor)
## This mirrors how the Catching et al. ExN genes were defined (linear age in
## adults), so the two references can be compared directly.
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
setwd("/restricted/projectnb/crem-trainees/Murphy_Lab/Elise/GeroScreen_Paper")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(edgeR)
})

in_dir   <- "PsychAD_EN"
min_age  <- 20     # adults only; the developmental period is a different biology
min_cells <- 50    # donors with fewer excitatory neurons are dropped
fdr_cut  <- 0.05
max_genes_per_direction <- 250   # cap so the signature is comparable in size to Lu 2004

## Some gene names appear more than once (different Ensembl IDs with the same
## symbol), so read the gene column as data and sum counts across duplicates
counts_raw <- read.csv(file.path(in_dir, "psychad_EN_pseudobulk_counts.csv.gz"),
                       check.names = FALSE, stringsAsFactors = FALSE)
gene_names <- as.character(counts_raw[[1]])
counts <- rowsum(as.matrix(counts_raw[, -1]), group = gene_names, reorder = TRUE)
message(sprintf("%d rows read; %d unique genes after merging duplicate names",
                length(gene_names), nrow(counts)))
rm(counts_raw)

meta <- read.csv(file.path(in_dir, "psychad_EN_donor_meta.csv"),
                 stringsAsFactors = FALSE)
meta$donor_id <- as.character(meta$donor_id)
stopifnot(all(meta$donor_id %in% colnames(counts)))
counts <- counts[, meta$donor_id]

## ---- 1. Donor filtering ----------------------------------------------------------
meta <- meta %>%
  mutate(keep = !is.na(age_years) & age_years >= min_age & n_cells >= min_cells)
message(sprintf("Donors: %d total, %d adults with >= %d excitatory neurons",
                nrow(meta), sum(meta$keep), min_cells))
meta   <- filter(meta, keep)
counts <- counts[, meta$donor_id]

meta <- meta %>%
  mutate(age_decades = age_years / 10,
         sex         = factor(sex),
         log_cells   = log(n_cells))

## ---- 2. voom pseudobulk ----------------------------------------------------------
dge <- DGEList(counts)
design <- if (nlevels(meta$sex) > 1) {
  model.matrix(~ age_decades + sex + log_cells, data = meta)
} else {
  model.matrix(~ age_decades + log_cells, data = meta)
}
keep_genes <- filterByExpr(dge, design)
dge <- calcNormFactors(dge[keep_genes, , keep.lib.sizes = FALSE])
message(sprintf("Genes tested: %d", nrow(dge)))

v   <- voom(dge, design)
fit <- eBayes(lmFit(v, design), robust = TRUE)

age_res <- topTable(fit, coef = "age_decades", number = Inf, sort.by = "P") %>%
  mutate(gene = rownames(.)) %>%
  select(gene, logFC_per_decade = logFC, AveExpr, t, P.Value, adj.P.Val)
write.csv(age_res, file.path(in_dir, "psychad_EN_age_genes_all.csv"), row.names = FALSE)

## ---- 3. Signature: significant genes, capped per direction -----------------------
sig <- age_res %>%
  filter(adj.P.Val < fdr_cut) %>%
  mutate(direction = ifelse(logFC_per_decade > 0, "up", "down")) %>%
  group_by(direction) %>%
  slice_min(P.Value, n = max_genes_per_direction, with_ties = FALSE) %>%
  ungroup()
write.csv(sig, file.path(in_dir, "psychad_EN_age_signature.csv"), row.names = FALSE)
message(sprintf("Signature: %d genes up with age, %d down with age (FDR < %.2f)",
                sum(sig$direction == "up"), sum(sig$direction == "down"), fdr_cut))

## ---- 4. QC plot ------------------------------------------------------------------
mds <- plotMDS(v, plot = FALSE)
qc  <- meta %>% mutate(MDS1 = mds$x, MDS2 = mds$y)
p_qc <- ggplot(qc, aes(MDS1, MDS2, colour = age_years, size = n_cells)) +
  geom_point(alpha = 0.8) +
  scale_colour_viridis_c(name = "Age (years)") +
  scale_size_continuous(range = c(1, 4), name = "Excitatory\nneurons") +
  labs(title = "PsychAD excitatory-neuron pseudobulk (adult donors)") +
  theme_bw(base_size = 10)
ggsave(file.path(in_dir, "psychad_EN_age_QC.png"), p_qc, width = 6, height = 4.5,
       dpi = 300, device = png, type = "cairo")

writeLines(capture.output(sessionInfo()), file.path(in_dir, "sessionInfo_signature.txt"))
message("Done.")