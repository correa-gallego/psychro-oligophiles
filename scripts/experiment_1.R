rm(list = ls())
set.seed(42)

req_pkgs <- c(
  "tidyverse", "janitor", "vegan", "patchwork",
  "ggtext", "cowplot", "glue", "scales", "grid", "gridExtra",
  "car", "effectsize", "iNEXT"
)

to_install <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(vegan)
  library(patchwork)
  library(ggtext)
  library(cowplot)
  library(glue)
  library(scales)
  library(grid)
  library(gridExtra)
  library(car)
  library(effectsize)
  library(iNEXT)
})

script_dir <- file.path(
  "/Users/scorrea/Library/CloudStorage",
  "GoogleDrive-correagsebastian2204@gmail.com",
  "My Drive/Correa-Gallego2026/Scripts"
)

input_dir <- file.path(script_dir, "pilot", "data")
out_dir   <- file.path(script_dir, "output", "pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_counts <- file.path(input_dir, "counts_dilution_e2.csv")
path_morph  <- file.path(input_dir, "morphotype_counts_pilot.csv")

plated_volume_mL <- 0.1
sediment_wet_mL_est <- 1.5
mother_fraction <- 0.2 / 1.0

make_sample_id <- function(zone, subzone, replicate) {
  zone <- toupper(zone)
  subzone <- toupper(subzone)
  replicate <- toupper(replicate)
  
  rep_num <- readr::parse_number(replicate)
  subzone_offset <- ifelse(subzone == "A", 0L, 3L)
  final_rep <- rep_num + subzone_offset
  zone_prefix <- c("Z1" = "L", "Z2" = "D")[zone]
  
  paste0(zone_prefix, "-R", final_rep)
}

extract_final_replicate <- function(sample_id) {
  paste0("R", readr::parse_number(sample_id))
}

theme_nature <- function(base_size = 11.4, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, color = "grey30"),
      axis.title = element_text(size = base_size + 0.3, face = "bold"),
      axis.text = element_text(size = base_size - 0.7),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.major = element_line(linewidth = 0.22, color = "grey94"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      legend.title = element_text(face = "bold", size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.9),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      strip.text = element_text(face = "bold", size = base_size - 0.3),
      strip.background = element_rect(color = NA, fill = "grey95"),
      plot.margin = margin(t = 6, r = 12, b = 6, l = 6)
    )
}

extract_es <- function(es_obj) {
  df <- as.data.frame(es_obj)
  val <- as.numeric(df[1, 1])
  ci_cols <- grep("CI_", names(df), value = TRUE)
  if (length(ci_cols) >= 2) {
    ci_vals <- sort(as.numeric(df[1, ci_cols]))
  } else {
    ci_cols <- grep("CI", names(df), value = TRUE)
    ci_cols <- ci_cols[ci_cols != "CI"]
    ci_vals <- sort(as.numeric(df[1, ci_cols]))
  }
  list(g = val, ci_low = ci_vals[1], ci_high = ci_vals[2])
}

palette_zones <- c("Light" = "#029E73", "Dark" = "#DE8F05")

fig_sizes <- list(
  qc = list(width = 175, height = 105),
  diagnostics = list(width = 205, height = 105),
  rarefaction = list(width = 175, height = 118),
  ordination = list(width = 128, height = 142),
  heatmap = list(width = 178, height = 168),
  simper = list(width = 198, height = 130),
  composition = list(width = 225, height = 138),
  rank = list(width = 145, height = 108),
  cfu_zone = list(width = 118, height = 108),
  alpha_panel = list(width = 235, height = 98),
  density_diversity_composite = list(width = 235, height = 190),
  composition_composite = list(width = 235, height = 205)
)

save_svg <- function(plot, filename, size, path = out_dir, show_legend = TRUE) {
  s <- fig_sizes[[size]]
  if (is.null(s)) stop("Unknown figure size key: ", size)
  if (!show_legend) plot <- plot + theme(legend.position = "none")
  ggsave(
    filename = file.path(path, filename),
    plot = plot,
    width = s$width,
    height = s$height,
    units = "mm",
    dpi = 300,
    bg = "white"
  )
}

save_pdf <- function(plot, filename, size, path = out_dir, show_legend = TRUE) {
  s <- fig_sizes[[size]]
  if (is.null(s)) stop("Unknown figure size key: ", size)
  if (!show_legend) plot <- plot + theme(legend.position = "none")
  ggsave(
    filename = file.path(path, filename),
    plot = plot,
    width = s$width,
    height = s$height,
    units = "mm",
    device = "pdf",
    bg = "white"
  )
}

counts <- read_csv(path_counts, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    zone = toupper(zone),
    subzone = toupper(subzone),
    replicate_raw = toupper(replicate),
    sample_id = make_sample_id(zone, subzone, replicate_raw),
    replicate = extract_final_replicate(sample_id),
    dilution = as.numeric(dilution),
    zone_label = c("Z1" = "Light", "Z2" = "Dark")[zone],
    zone_label = factor(zone_label, levels = c("Light", "Dark"))
  )

stopifnot(nrow(counts) == 12)

morph_raw <- read_csv(path_morph, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    zone = toupper(zone),
    subzone = toupper(subzone),
    replicate_raw = toupper(replicate),
    sample_id = make_sample_id(zone, subzone, replicate_raw),
    replicate = extract_final_replicate(sample_id)
  )

morph_cols <- names(morph_raw)[str_detect(names(morph_raw), "^org_sa")]

cat("\nMorphotype columns:", length(morph_cols), "\n")
cat("Codes:", paste(morph_cols, collapse = ", "), "\n")

qc_check <- morph_raw %>%
  mutate(
    sum_morphotypes = rowSums(across(all_of(morph_cols)), na.rm = TRUE),
    delta = colony_count - sum_morphotypes,
    pct_discrepancy = round(100 * delta / colony_count, 1)
  ) %>%
  select(sample_id, colony_count, sum_morphotypes, delta, pct_discrepancy)

cat("\nQC summary\n")
print(as.data.frame(qc_check))
cat("\nMean absolute discrepancy:", round(mean(abs(qc_check$delta)), 1), "\n")
cat("Max absolute discrepancy:", max(abs(qc_check$delta)), "\n")

write_csv(qc_check, file.path(out_dir, "QC_morphotype_sum_vs_count.csv"))

p_qc <- ggplot(qc_check, aes(x = reorder(sample_id, delta), y = delta)) +
  geom_col(aes(fill = delta > 0), width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#4575b4", "FALSE" = "#d73027")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Colony count − Σ morphotypes",
    title = NULL
  ) +
  theme_nature()

save_svg(p_qc, "QC_discrepancy.svg", "qc")

cfu_df <- counts %>%
  mutate(
    dilution_total = dilution * mother_fraction,
    cfu_per_mL_slurry = colony_count / (dilution_total * plated_volume_mL),
    cfu_per_mL_wet_sed = cfu_per_mL_slurry / sediment_wet_mL_est,
    log10_cfu = log10(cfu_per_mL_wet_sed)
  )

cfu_summary <- cfu_df %>%
  group_by(zone_label) %>%
  summarise(
    n = n(),
    mean = mean(log10_cfu),
    sd = sd(log10_cfu),
    se = sd / sqrt(n),
    min = min(log10_cfu),
    max = max(log10_cfu),
    .groups = "drop"
  )

cat("\nCFU summary by zone\n")
print(as.data.frame(cfu_summary))

write_csv(cfu_summary, file.path(out_dir, "CFU_summary_by_zone.csv"))

m_cfu <- lm(log10_cfu ~ zone_label, data = cfu_df)

sw_cfu <- shapiro.test(residuals(m_cfu))
lev_cfu <- leveneTest(log10_cfu ~ zone_label, data = cfu_df)

cat("\nCFU diagnostics\n")
cat("Shapiro-Wilk W =", round(sw_cfu$statistic, 4), " p =", round(sw_cfu$p.value, 4), "\n")
cat("Levene F =", round(lev_cfu$`F value`[1], 4), " p =", round(lev_cfu$`Pr(>F)`[1], 4), "\n")

diag_df <- tibble(
  fitted = fitted(m_cfu),
  residual = residuals(m_cfu)
)

p_qq <- ggplot(diag_df, aes(sample = residual)) +
  geom_qq(size = 2.35, alpha = 0.85) +
  geom_qq_line(color = "red", linewidth = 0.65) +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles",
    title = NULL
  ) +
  theme_nature()

p_resfit <- ggplot(diag_df, aes(x = fitted, y = residual)) +
  geom_point(size = 2.35, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.65) +
  labs(
    x = "Fitted values",
    y = "Residuals",
    title = NULL
  ) +
  theme_nature()

p_diag <- p_qq + p_resfit +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_svg(p_diag, "CFU_diagnostics.svg", "diagnostics")

t_cfu <- t.test(log10_cfu ~ zone_label, data = cfu_df, var.equal = FALSE)
w_cfu <- wilcox.test(log10_cfu ~ zone_label, data = cfu_df, exact = FALSE)
d_cfu_raw <- cohens_d(log10_cfu ~ zone_label, data = cfu_df, hedges.correction = TRUE)
d_cfu <- extract_es(d_cfu_raw)

cat("\nCFU tests\n")
cat("Welch t =", round(t_cfu$statistic, 3), " df =", round(t_cfu$parameter, 1), " p =", format.pval(t_cfu$p.value, digits = 4), "\n")
cat("Wilcoxon W =", w_cfu$statistic, " p =", format.pval(w_cfu$p.value, digits = 4), "\n")
cat("Hedges g =", round(d_cfu$g, 3), " [", round(d_cfu$ci_low, 3), ",", round(d_cfu$ci_high, 3), "]\n")

comm_wide <- morph_raw %>%
  select(sample_id, all_of(morph_cols)) %>%
  column_to_rownames("sample_id")

comm <- as.matrix(comm_wide)

stopifnot(all(comm >= 0))
stopifnot(all(comm == round(comm)))

cat("\nCommunity matrix\n")
cat("Samples:", nrow(comm), "\n")
cat("Morphotypes:", ncol(comm), "\n")
cat("Total colonies:", sum(comm), "\n")
cat("Sample total range:", min(rowSums(comm)), "to", max(rowSums(comm)), "\n")

comm_hell <- decostand(comm, method = "hellinger")

alpha_df <- tibble(sample_id = rownames(comm)) %>%
  mutate(
    richness = specnumber(comm),
    shannon = diversity(comm, index = "shannon"),
    simpson = diversity(comm, index = "simpson")
  ) %>%
  left_join(
    cfu_df %>% select(sample_id, zone_label, zone, replicate),
    by = "sample_id"
  )

cat("\nAlpha diversity summary\n")
alpha_df %>%
  group_by(zone_label) %>%
  summarise(
    across(c(richness, shannon, simpson), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  print(width = Inf)

m_rich <- lm(richness ~ zone_label, data = alpha_df)
sw_rich <- shapiro.test(residuals(m_rich))
t_rich <- t.test(richness ~ zone_label, data = alpha_df, var.equal = FALSE)
w_rich <- wilcox.test(richness ~ zone_label, data = alpha_df, exact = FALSE)
d_rich_raw <- cohens_d(richness ~ zone_label, data = alpha_df, hedges.correction = TRUE)
d_rich <- extract_es(d_rich_raw)

m_sha <- lm(shannon ~ zone_label, data = alpha_df)
sw_sha <- shapiro.test(residuals(m_sha))
t_sha <- t.test(shannon ~ zone_label, data = alpha_df, var.equal = FALSE)
w_sha <- wilcox.test(shannon ~ zone_label, data = alpha_df, exact = FALSE)
d_sha_raw <- cohens_d(shannon ~ zone_label, data = alpha_df, hedges.correction = TRUE)
d_sha <- extract_es(d_sha_raw)

m_sim <- lm(simpson ~ zone_label, data = alpha_df)
sw_sim <- shapiro.test(residuals(m_sim))
t_sim <- t.test(simpson ~ zone_label, data = alpha_df, var.equal = FALSE)
w_sim <- wilcox.test(simpson ~ zone_label, data = alpha_df, exact = FALSE)
d_sim_raw <- cohens_d(simpson ~ zone_label, data = alpha_df, hedges.correction = TRUE)
d_sim <- extract_es(d_sim_raw)

cat("\nAlpha diversity tests\n")
cat("Richness Welch p =", format.pval(t_rich$p.value, digits = 4), " Wilcoxon p =", format.pval(w_rich$p.value, digits = 4), "\n")
cat("Shannon Welch p =", format.pval(t_sha$p.value, digits = 4), " Wilcoxon p =", format.pval(w_sha$p.value, digits = 4), "\n")
cat("Simpson Welch p =", format.pval(t_sim$p.value, digits = 4), " Wilcoxon p =", format.pval(w_sim$p.value, digits = 4), "\n")

inext_list <- as.list(as.data.frame(t(comm)))
names(inext_list) <- rownames(comm)

inext_out <- iNEXT(inext_list, q = 0, datatype = "abundance", conf = 0.95, nboot = 200)

if (is.data.frame(inext_out$iNextEst$size_based)) {
  rare_df <- inext_out$iNextEst$size_based %>%
    as_tibble() %>%
    rename(sample_id = Assemblage, x = m, y = qD, y.lwr = qD.LCL, y.upr = qD.UCL)
} else {
  rare_df <- bind_rows(
    lapply(names(inext_out$iNextEst$size_based), function(nm) {
      inext_out$iNextEst$size_based[[nm]] %>%
        as_tibble() %>%
        mutate(sample_id = nm)
    })
  ) %>%
    rename(x = m, y = qD, y.lwr = qD.LCL, y.upr = qD.UCL)
}

rare_df <- rare_df %>%
  mutate(
    sample_id = as.character(sample_id),
    zone_label = case_when(
      str_detect(sample_id, "^L-") ~ "Light",
      str_detect(sample_id, "^D-") ~ "Dark"
    ),
    zone_label = factor(zone_label, levels = c("Light", "Dark")),
    Method = as.character(Method)
  )

p_rare <- ggplot(rare_df, aes(x = x, y = y, color = zone_label)) +
  geom_ribbon(
    aes(ymin = y.lwr, ymax = y.upr, fill = zone_label, group = sample_id),
    alpha = 0.08,
    color = NA
  ) +
  geom_line(
    data = rare_df %>% filter(Method == "Rarefaction"),
    aes(group = sample_id),
    linewidth = 0.58,
    linetype = "solid"
  ) +
  geom_line(
    data = rare_df %>% filter(Method == "Extrapolation"),
    aes(group = sample_id),
    linewidth = 0.58,
    linetype = "dashed"
  ) +
  geom_point(
    data = rare_df %>% filter(Method == "Observed"),
    aes(group = sample_id),
    size = 2.35,
    shape = 16
  ) +
  scale_color_manual(values = palette_zones) +
  scale_fill_manual(values = palette_zones, guide = "none") +
  labs(
    x = "Number of individuals",
    y = "Morphotype richness (q = 0)",
    color = "Zone",
    title = NULL
  ) +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  )

save_svg(p_rare, "rarefaction_curves.svg", "rarefaction")

bc <- vegdist(comm_hell, method = "bray")

nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 200, autotransform = FALSE, trace = 0)

cat("\nOrdination fit\n")
cat("NMDS stress =", round(nmds$stress, 4), "\n")

scores_nmds <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("sample_id") %>%
  left_join(cfu_df %>% select(sample_id, zone_label), by = "sample_id")

p_nmds <- ggplot(scores_nmds, aes(x = NMDS1, y = NMDS2, color = zone_label)) +
  geom_point(size = 3.35, alpha = 0.9) +
  stat_ellipse(
    aes(group = zone_label, fill = zone_label),
    geom = "polygon",
    alpha = 0.10,
    level = 0.68,
    linetype = "dashed",
    linewidth = 0.42
  ) +
  scale_color_manual(values = palette_zones) +
  scale_fill_manual(values = palette_zones, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Zone",
    title = NULL
  ) +
  coord_equal(clip = "off") +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  )

save_svg(p_nmds, "NMDS_braycurtis.svg", "ordination")

pcoa <- cmdscale(bc, eig = TRUE, k = 2)
eig_pos <- pcoa$eig[pcoa$eig > 0]
var_exp <- 100 * pcoa$eig[1:2] / sum(eig_pos)

scores_pcoa <- tibble(
  sample_id = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(cfu_df %>% select(sample_id, zone_label), by = "sample_id")

p_pcoa <- ggplot(scores_pcoa, aes(x = PC1, y = PC2, color = zone_label)) +
  geom_point(size = 3.35, alpha = 0.9) +
  stat_ellipse(
    aes(group = zone_label, fill = zone_label),
    geom = "polygon",
    alpha = 0.10,
    level = 0.68,
    linetype = "dashed",
    linewidth = 0.42
  ) +
  scale_color_manual(values = palette_zones) +
  scale_fill_manual(values = palette_zones, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = glue("PC1 ({round(var_exp[1], 1)}%)"),
    y = glue("PC2 ({round(var_exp[2], 1)}%)"),
    color = "Zone",
    title = NULL
  ) +
  coord_equal(clip = "off") +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  )

save_svg(p_pcoa, "PCoA_braycurtis.svg", "ordination")

bc_mat <- as.matrix(bc)
hc <- hclust(as.dist(bc_mat), method = "average")
ord <- hc$order

bc_long <- as.data.frame(as.table(bc_mat[ord, ord])) %>%
  rename(sample_i = Var1, sample_j = Var2, bray = Freq) %>%
  mutate(
    sample_i = factor(sample_i, levels = rownames(bc_mat)[ord]),
    sample_j = factor(sample_j, levels = rownames(bc_mat)[ord])
  )

p_heatmap <- ggplot(bc_long, aes(x = sample_i, y = sample_j, fill = bray)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(bray, 2)), size = 2.45, color = "black") +
  scale_fill_viridis_c(option = "D", direction = -1, name = "Bray-Curtis") +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 1
  )

save_svg(p_heatmap, "BrayCurtis_heatmap.svg", "heatmap")

perm_meta <- tibble(sample_id = rownames(comm_hell)) %>%
  left_join(cfu_df %>% select(sample_id, zone_label), by = "sample_id")

set.seed(42)
perm_zone <- adonis2(
  comm_hell ~ zone_label,
  data = perm_meta,
  permutations = 999,
  method = "bray",
  by = "margin"
)

bd_zone <- betadisper(bc, perm_meta$zone_label)
bd_anova <- anova(bd_zone)

set.seed(42)
bd_perm <- permutest(bd_zone, permutations = 999)
bd_perm_F <- bd_perm$tab[1, "F"]
bd_perm_p <- bd_perm$tab[1, "Pr(>F)"]

cat("\nPERMANOVA and dispersion\n")
print(perm_zone)
cat("\nBetadisper ANOVA F =", round(bd_anova$`F value`[1], 3), " p =", round(bd_anova$`Pr(>F)`[1], 4), "\n")
cat("Betadisper permutation F =", round(bd_perm_F, 3), " p =", round(bd_perm_p, 4), "\n")

set.seed(42)
sim_result <- simper(comm, perm_meta$zone_label, permutations = 999)

sim_summary_list <- summary(sim_result)
comparison_name <- names(sim_summary_list)[1]
sim_summary <- sim_summary_list[[1]]

comp_groups <- strsplit(comparison_name, "_")[[1]]
group_a <- comp_groups[1]
group_b <- comp_groups[2]

simper_df <- tibble(
  morphotype = toupper(str_replace(rownames(sim_summary), "org_sa", "ORG-SA")),
  contribution = sim_summary$average,
  sd = sim_summary$sd,
  ratio = sim_summary$ratio,
  cumulative = sim_summary$cumsum,
  p_value = sim_summary$p
)

simper_df[[paste0("mean_", group_a)]] <- sim_summary$ava
simper_df[[paste0("mean_", group_b)]] <- sim_summary$avb

simper_df <- simper_df %>%
  arrange(desc(contribution))

write_csv(simper_df, file.path(out_dir, "SIMPER_results.csv"))

simper_plot_df <- simper_df %>%
  slice_head(n = 10) %>%
  mutate(morphotype = factor(morphotype, levels = rev(morphotype)))

p_simper <- ggplot(simper_plot_df, aes(x = morphotype, y = contribution)) +
  geom_col(aes(fill = cumulative <= 0.70), width = 0.7, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = pmax(contribution - sd, 0), ymax = contribution + sd),
    width = 0.2,
    linewidth = 0.3
  ) +
  geom_text(aes(label = paste0(round(100 * cumulative, 1), "%")), hjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean contribution to Bray-Curtis dissimilarity",
    title = NULL
  ) +
  theme_nature()

save_svg(p_simper, "SIMPER_contributions.svg", "simper")

sample_order <- morph_raw %>%
  mutate(
    zone_label = c("Z1" = "Light", "Z2" = "Dark")[toupper(zone)],
    total = rowSums(across(all_of(morph_cols)))
  ) %>%
  arrange(zone_label, total) %>%
  pull(sample_id)

morph_long <- morph_raw %>%
  select(sample_id, all_of(morph_cols)) %>%
  pivot_longer(cols = all_of(morph_cols), names_to = "morphotype", values_to = "count") %>%
  mutate(
    morphotype = toupper(str_replace(morphotype, "org_sa", "ORG-SA")),
    zone_label = case_when(
      str_detect(sample_id, "^L-") ~ "Light",
      str_detect(sample_id, "^D-") ~ "Dark"
    ),
    zone_label = factor(zone_label, levels = c("Light", "Dark")),
    sample_id = factor(sample_id, levels = sample_order)
  ) %>%
  group_by(sample_id) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

morph_order <- morph_long %>%
  group_by(morphotype) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(morphotype)

morph_long <- morph_long %>%
  mutate(morphotype = factor(morphotype, levels = rev(morph_order)))

morph_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"
)
names(morph_palette) <- morph_order

p_comp <- ggplot(morph_long, aes(x = sample_id, y = proportion, fill = morphotype)) +
  geom_col(width = 0.8) +
  facet_wrap(~ zone_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = morph_palette, name = "Morphotype") +
  labs(
    x = NULL,
    y = "Relative abundance",
    title = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.8),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7.2)
  )

save_svg(p_comp, "composition_per_replicate.svg", "composition")

rank_df <- morph_long %>%
  group_by(zone_label, morphotype) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop") %>%
  group_by(zone_label) %>%
  arrange(desc(mean_prop)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p_rank <- ggplot(rank_df, aes(x = rank, y = mean_prop, color = zone_label)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.25) +
  scale_color_manual(values = palette_zones) +
  scale_y_log10(labels = scales::label_number()) +
  labs(
    x = "Morphotype rank",
    y = "Mean relative abundance (log scale)",
    color = "Zone",
    title = NULL
  ) +
  theme_nature()

save_svg(p_rank, "rank_abundance.svg", "rank")

p_cfu_zone <- ggplot(cfu_df, aes(x = zone_label, y = log10_cfu, color = zone_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
  geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
  scale_color_manual(values = palette_zones) +
  labs(
    x = NULL,
    y = expression(log[10]~"(CFU mL"^{-1}~"wet sediment)"),
    title = NULL
  ) +
  theme_nature() +
  theme(legend.position = "none")

save_svg(p_cfu_zone, "CFU_by_zone.svg", "cfu_zone")

p_alpha_rich <- ggplot(alpha_df, aes(x = zone_label, y = richness, color = zone_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
  geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
  scale_color_manual(values = palette_zones) +
  labs(x = NULL, y = "Morphotype richness", title = NULL) +
  theme_nature() +
  theme(legend.position = "none")

p_alpha_sha <- ggplot(alpha_df, aes(x = zone_label, y = shannon, color = zone_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
  geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
  scale_color_manual(values = palette_zones) +
  labs(x = NULL, y = "Shannon (H')", title = NULL) +
  theme_nature() +
  theme(legend.position = "none")

p_alpha_sim <- ggplot(alpha_df, aes(x = zone_label, y = simpson, color = zone_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
  geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
  scale_color_manual(values = palette_zones) +
  labs(x = NULL, y = "Simpson (1 − D)", title = NULL) +
  theme_nature() +
  theme(legend.position = "none")

p_alpha_panel <- p_alpha_rich + p_alpha_sha + p_alpha_sim +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_svg(p_alpha_panel, "alpha_diversity_panel.svg", "alpha_panel")

density_diversity_composite <- (
  p_cfu_zone + p_alpha_rich +
    p_alpha_sha + p_alpha_sim
) +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_pdf(density_diversity_composite, "exp1_density_diversity.pdf", "density_diversity_composite")

composition_composite <- (
  (p_nmds + p_pcoa) /
    p_comp
) +
  plot_layout(heights = c(1, 1.35)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_pdf(composition_composite, "exp1_composition.pdf", "composition_composite")

within_bc <- function(zone_lbl) {
  idx <- which(perm_meta$zone_label == zone_lbl)
  d <- as.matrix(bc)[idx, idx]
  d[upper.tri(d)]
}

within_light <- within_bc("Light")
within_dark <- within_bc("Dark")

cat("\nWithin-zone Bray-Curtis distances\n")
cat("Light mean =", round(mean(within_light), 3), " range =", round(min(within_light), 3), "to", round(max(within_light), 3), "\n")
cat("Dark mean =", round(mean(within_dark), 3), " range =", round(min(within_dark), 3), "to", round(max(within_dark), 3), "\n")

report_file <- file.path(out_dir, "pilot_stats_report.txt")
sink(report_file)

cat("Organal San Antonio — Pilot Experiment\n")
cat("Statistical report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("Design\n")
cat("Zones: Light and Dark\n")
cat("Replicates: 6 per zone, 12 samples total\n")
cat("Single dilution: 10^-2\n")
cat("Morphotypes observed: 14 of 19 described\n\n")

cat("Inference\n")
cat("The pilot experiment was analyzed as a two-zone design with six biological replicates per zone.\n")
cat("Within-zone heterogeneity is descriptive and exploratory.\n\n")

cat("CFU summary\n")
print(as.data.frame(cfu_summary), row.names = FALSE)
cat("\n")

cat("CFU diagnostics\n")
cat("Shapiro-Wilk W =", round(sw_cfu$statistic, 4), " p =", round(sw_cfu$p.value, 4), "\n")
cat("Levene F =", round(lev_cfu$`F value`[1], 4), " p =", round(lev_cfu$`Pr(>F)`[1], 4), "\n\n")

cat("CFU tests\n")
cat("Welch t =", round(t_cfu$statistic, 3), " df =", round(t_cfu$parameter, 1), " p =", format.pval(t_cfu$p.value, digits = 4), "\n")
cat("95% CI =", round(t_cfu$conf.int[1], 3), "to", round(t_cfu$conf.int[2], 3), "\n")
cat("Wilcoxon W =", w_cfu$statistic, " p =", format.pval(w_cfu$p.value, digits = 4), "\n")
cat("Hedges g =", round(d_cfu$g, 3), " [", round(d_cfu$ci_low, 3), ",", round(d_cfu$ci_high, 3), "]\n\n")

cat("Alpha diversity tests\n")
cat("Richness: Welch p =", format.pval(t_rich$p.value, digits = 4), " Wilcoxon p =", format.pval(w_rich$p.value, digits = 4), " Hedges g =", round(d_rich$g, 3), "\n")
cat("Shannon: Welch p =", format.pval(t_sha$p.value, digits = 4), " Wilcoxon p =", format.pval(w_sha$p.value, digits = 4), " Hedges g =", round(d_sha$g, 3), "\n")
cat("Simpson: Welch p =", format.pval(t_sim$p.value, digits = 4), " Wilcoxon p =", format.pval(w_sim$p.value, digits = 4), " Hedges g =", round(d_sim$g, 3), "\n\n")

cat("Beta diversity\n")
cat("NMDS stress =", round(nmds$stress, 4), "\n\n")
print(perm_zone)
cat("\n")
cat("Betadisper ANOVA F =", round(bd_anova$`F value`[1], 3), " p =", round(bd_anova$`Pr(>F)`[1], 4), "\n")
cat("Betadisper permutation F =", round(bd_perm_F, 3), " p =", round(bd_perm_p, 4), "\n\n")

cat("SIMPER top contributors\n")
print(as.data.frame(simper_df %>% head(10)), row.names = FALSE)
cat("\n")

cat("QC summary\n")
print(as.data.frame(qc_check), row.names = FALSE)
cat("\n")

sink()

analysis_df <- cfu_df %>%
  select(sample_id, zone, replicate, zone_label, colony_count, dilution, cfu_per_mL_wet_sed, log10_cfu) %>%
  left_join(alpha_df %>% select(sample_id, richness, shannon, simpson), by = "sample_id")

write_csv(analysis_df, file.path(out_dir, "pilot_analysis_table.csv"))

cat("\nOutput directory\n")
cat(out_dir, "\n\n")
cat("Files generated\n")
list.files(out_dir) %>% walk(~ cat(.x, "\n"))

sink(file.path(out_dir, "session_info.txt"))
cat("Session info — Pilot analysis\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
sessionInfo()
sink()

cat("\nPilot analysis complete\n")
cat("Outputs saved in:", out_dir, "\n")