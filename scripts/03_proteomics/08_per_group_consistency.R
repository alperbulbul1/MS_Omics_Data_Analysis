#!/usr/bin/env Rscript
## 08_per_group_consistency.R  —  generated from notebook spec
## Run: Rscript 08_per_group_consistency.R


## ============================================================
## # 08 — Per-group cross-study consistency (29 tissue × cell-type groups)
## 
## R port of Python `02_per_group_consistency.py`. Reads
## `per_group_results/all_DEPs_with_group.csv` and for each group
## computes per-gene `(n_studies_hit, n_up, n_down, signed_consistency)`.
## 
## **Outputs**
## - `per_group_results/<group>/consistency_R.tsv` × 29
## - `figures/all_groups_grid_R.png`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(ggrepel)
  library(gridExtra)
})
source("helpers.R")
PG <- file.path(PROT_ROOT, "per_group_results")


rows <- fread(file.path(PG, "all_DEPs_with_group.csv"))
rows <- rows[!cell_type %in% c("REVIEW","META","UNKNOWN") & gene != ""]
cat(sprintf("Loaded %d DEPs across %d groups · %d studies\n",
            nrow(rows), uniqueN(rows$group), uniqueN(rows$study_id)))

per_gene <- rows[, .(n_studies_hit = uniqueN(study_id),
                     n_up   = sum(direction == "Up"),
                     n_down = sum(direction == "Down"),
                     n_other = .N - sum(direction %in% c("Up","Down")),
                     studies = paste(sort(unique(study_id)), collapse = ";"),
                     directions = paste(direction, collapse = ";")),
                 by = .(group, gene)]
per_gene[, signed_consistency := (n_up - n_down) / pmax(n_up + n_down, 1)]
per_gene[, is_cross_omics := gene %in% CROSS_OMICS]
per_gene[, is_recurring   := gene %in% RECURRING]
per_gene[, is_ECM         := gene %in% ECM_FAMILY]

# Save per-group tables
for (g in unique(per_gene$group)) {
  out_dir <- file.path(PG, g)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sub <- per_gene[group == g][order(-n_studies_hit, -abs(signed_consistency))]
  fwrite(sub, file.path(out_dir, "consistency_R.tsv"), sep = "\t")
}
cat(sprintf("Wrote per-group consistency_R.tsv for %d groups\n",
            uniqueN(per_gene$group)))


# ---- group summary ----
group_summary <- per_gene[, .(
  n_studies = uniqueN(rows[group == .BY$group]$study_id),
  n_unique_genes = .N,
  n_recurring_in_group_2plus = sum(n_studies_hit >= 2),
  n_cross_omics_hit = sum(is_cross_omics),
  n_recurring_signature_hit = sum(is_recurring),
  n_ECM_family_hit = sum(is_ECM)
), by = group][order(-n_studies, -n_unique_genes)]

print(group_summary)
fwrite(group_summary, file.path(PG, "_group_summary_R.tsv"), sep = "\t")


# ---- composite grid of all groups ----
make_panel <- function(grp_name) {
  d <- per_gene[group == grp_name]
  if (nrow(d) == 0) return(NULL)
  n_studies <- group_summary[group == grp_name]$n_studies
  ggplot(d, aes(signed_consistency, n_studies_hit)) +
    geom_jitter(width = 0.02, height = 0.1,
                aes(colour = factor(pmin(n_studies_hit, 3))),
                size = 1.2, alpha = 0.7) +
    scale_colour_manual(values = c("1"="#CCCCCC","2"="#3E92CC","3"="#1F4E79"),
                        guide = "none") +
    geom_point(data = subset(d, is_ECM),
               shape = 22, fill = "#F4A261", colour = "black",
               size = 2.6, stroke = 0.3, alpha = 0.9) +
    geom_point(data = subset(d, is_recurring),
               shape = 23, fill = "#7B3FA0", colour = "black",
               size = 3.0, stroke = 0.3, alpha = 0.9) +
    geom_point(data = subset(d, is_cross_omics),
               shape = 8, colour = "#D62828", size = 5.5, stroke = 1.4) +
    geom_text_repel(data = subset(d, is_cross_omics),
                    aes(label = gene), size = 2.8, fontface = "bold",
                    colour = "#D62828") +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    scale_x_continuous(limits = c(-1.15, 1.15)) +
    ggtitle(sprintf("%s\n(%d studies, %d prot)",
                    gsub("__", "·", grp_name),
                    n_studies, nrow(d))) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(size = 7.5, face = "bold"),
          axis.title = element_blank())
}

groups_ord <- group_summary$group
panels <- Filter(Negate(is.null), lapply(groups_ord, make_panel))
ga <- gridExtra::arrangeGrob(grobs = panels, ncol = 4)
ggsave(file.path(FIG_DIR, "all_groups_grid_R.png"), ga,
       width = 22, height = 4 * ceiling(length(panels)/4), dpi = 150,
       limitsize = FALSE)
cat(sprintf("Wrote all_groups_grid_R.png (%d panels)\n", length(panels)))

