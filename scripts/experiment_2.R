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

data_dir <- file.path(script_dir, "current", "data")
out_dir  <- file.path(script_dir, "output", "current")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plated_volume_mL <- 0.1
sediment_total_volume_mL <- 1.5
sediment_aliquot_mL <- 0.2

zone_order <- c("Entrance", "Transition", "Dark")
zone_codes <- c("E" = "Entrance", "T" = "Transition", "D" = "Dark")
palette_3z <- c("Entrance" = "#029E73", "Transition" = "#0173B2", "Dark" = "#DE8F05")

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

pairwise_adonis2 <- function(dist_obj, group, permutations = 999,
                             p_adjust_method = "BH", seed = 42) {
  set.seed(seed)
  stopifnot(length(group) == attr(dist_obj, "Size"))
  levs <- levels(factor(group))
  res  <- vector("list", 0)
  
  for (i in seq_along(levs)[-length(levs)]) {
    for (j in (i + 1):length(levs)) {
      idx   <- which(group %in% c(levs[i], levs[j]))
      d_sub <- as.dist(as.matrix(dist_obj)[idx, idx])
      meta  <- tibble(g = droplevels(factor(group[idx])))
      a     <- adonis2(d_sub ~ g, data = meta, permutations = permutations)
      res[[length(res) + 1]] <- tibble(
        contrast = paste(levs[i], "vs", levs[j]),
        R2       = a[1, "R2"],
        F_stat   = a[1, "F"],
        p_value  = a[1, "Pr(>F)"]
      )
    }
  }
  
  bind_rows(res) %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
}

fig_sizes <- list(
  qc = list(width = 175, height = 105),
  diagnostics = list(width = 205, height = 105),
  rarefaction = list(width = 180, height = 120),
  ordination = list(width = 140, height = 148),
  heatmap = list(width = 182, height = 170),
  simper = list(width = 198, height = 150),
  composition_zone = list(width = 220, height = 138),
  composition_cross = list(width = 225, height = 145),
  rank = list(width = 150, height = 110),
  cfu_zone = list(width = 125, height = 110),
  alpha_panel = list(width = 245, height = 100),
  mass = list(width = 125, height = 110),
  crosswalk = list(width = 225, height = 185),
  density_diversity_composite = list(width = 255, height = 175),
  composition_composite = list(width = 260, height = 215)
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

save_pdf <- function(plot, filename, size, path = out_dir, use_dingbats = FALSE) {
  s <- fig_sizes[[size]]
  if (is.null(s)) stop("Unknown figure size key: ", size)
  
  grDevices::pdf(
    file = file.path(path, filename),
    width = s$width / 25.4,
    height = s$height / 25.4,
    family = "sans",
    useDingbats = use_dingbats
  )
  print(plot)
  dev.off()
}

counts_raw <- read_csv(
  file.path(data_dir, "colony_counts_day7.csv"),
  col_types = cols(zone = col_character()),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    zone_label = factor(zone_codes[zone], levels = zone_order),
    replicate  = as.integer(replicate)
  )

cat("\nTotal plates:", nrow(counts_raw), "\n")
cat("By zone and category\n")
print(counts_raw %>% count(zone_label, category))
cat("\n")

mass <- read_csv(
  file.path(data_dir, "sediment_mass.csv"),
  col_types = cols(zone = col_character()),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    zone_label = factor(zone_codes[zone], levels = zone_order),
    replicate  = as.integer(replicate)
  )

cat("Sediment mass summary (g)\n")
mass %>%
  group_by(zone_label) %>%
  summarise(mean = mean(mass_g), sd = sd(mass_g), min = min(mass_g), max = max(mass_g), .groups = "drop") %>%
  print()
cat("\n")

morph_E <- read_csv(
  file.path(data_dir, "morphotype_matrix_E.csv"),
  col_types = cols(zone = col_character()),
  show_col_types = FALSE
) %>% clean_names()

morph_T <- read_csv(
  file.path(data_dir, "morphotype_matrix_T.csv"),
  col_types = cols(zone = col_character()),
  show_col_types = FALSE
) %>% clean_names()

morph_D <- read_csv(
  file.path(data_dir, "morphotype_matrix_D.csv"),
  col_types = cols(zone = col_character()),
  show_col_types = FALSE
) %>% clean_names()

control <- read_csv(
  file.path(data_dir, "cave_water_control.csv"),
  show_col_types = FALSE
) %>% clean_names()

crosswalk <- read_csv(
  file.path(data_dir, "morphotype_crosswalk_v3.csv"),
  show_col_types = FALSE
) %>% clean_names()

qc_all <- bind_rows(
  morph_E %>% mutate(zone = "E"),
  morph_T %>% mutate(zone = "T"),
  morph_D %>% mutate(zone = "D")
) %>%
  select(plate_id, zone, replicate, dilution, is_tech_rep, n_original, n_effective, delta)

cat("QC discrepancy summary\n")
cat("Plates with |delta| > 10\n")
qc_all %>% filter(abs(delta) > 10) %>% print()
cat("\nDelta summary\n")
print(summary(qc_all$delta))
cat("\n")

write_csv(qc_all, file.path(out_dir, "QC_N_discrepancy.csv"))

p_qc_delta <- ggplot(qc_all, aes(x = reorder(plate_id, delta), y = delta, fill = zone)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(
    values = c("E" = "#029E73", "T" = "#0173B2", "D" = "#DE8F05"),
    labels = c("E" = "Entrance", "T" = "Transition", "D" = "Dark")
  ) +
  coord_flip() +
  labs(x = NULL, y = "N_original − N_effective", fill = "Zone", title = NULL) +
  theme_nature()

save_svg(p_qc_delta, "QC_N_discrepancy.svg", "qc")

tech_pairs <- bind_rows(morph_E, morph_T, morph_D) %>%
  select(plate_id, zone, replicate, dilution, is_tech_rep, n_effective) %>%
  group_by(zone, replicate, dilution) %>%
  filter(n() > 1) %>%
  summarise(
    n_primary  = n_effective[!is_tech_rep][1],
    n_tech_rep = n_effective[is_tech_rep][1],
    ratio      = n_tech_rep / n_primary,
    diff       = abs(n_primary - n_tech_rep),
    .groups    = "drop"
  )

cat("Technical replicate concordance\n")
print(as.data.frame(tech_pairs))
cat("\nMean |diff|:", round(mean(tech_pairs$diff, na.rm = TRUE), 1), "\n")
cat("Max |diff|:", max(tech_pairs$diff, na.rm = TRUE), "\n\n")

write_csv(tech_pairs, file.path(out_dir, "QC_tech_rep_concordance.csv"))

cfu_plates <- counts_raw %>%
  filter(!is_tech_rep, category %in% c("Countable", "Low-count")) %>%
  group_by(zone, replicate) %>%
  arrange(factor(category, levels = c("Countable", "Low-count")), desc(dilution)) %>%
  slice_head(n = 1) %>%
  ungroup()

cfu_df <- cfu_plates %>%
  left_join(mass %>% select(zone, replicate, mass_g), by = c("zone", "replicate")) %>%
  mutate(
    zone_label = factor(zone_codes[zone], levels = zone_order),
    mass_aliquot_g = mass_g * sediment_aliquot_mL / sediment_total_volume_mL,
    cfu_per_g = colony_count / (dilution * plated_volume_mL * mass_aliquot_g),
    log10_cfu = log10(cfu_per_g)
  )

cat("CFU plates selected\n")
cfu_df %>%
  select(
    plate_id, zone_label, replicate, dilution, colony_count,
    mass_g, mass_aliquot_g, cfu_per_g, log10_cfu
  ) %>%
  print(width = Inf)
cat("\n")

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

cat("CFU summary (log10 CFU/g sediment)\n")
print(as.data.frame(cfu_summary))
cat("\n")

write_csv(cfu_summary, file.path(out_dir, "CFU_summary_by_zone.csv"))

m_cfu <- lm(log10_cfu ~ zone_label, data = cfu_df)

sw_cfu <- shapiro.test(residuals(m_cfu))
lev_cfu <- leveneTest(log10_cfu ~ zone_label, data = cfu_df)

cat("CFU diagnostics\n")
cat("Shapiro-Wilk W =", round(sw_cfu$statistic, 4), " p =", round(sw_cfu$p.value, 4), "\n")
cat("Levene F =", round(lev_cfu$`F value`[1], 4), " p =", round(lev_cfu$`Pr(>F)`[1], 4), "\n\n")

diag_df <- tibble(fitted = fitted(m_cfu), residual = residuals(m_cfu))

p_qq <- ggplot(diag_df, aes(sample = residual)) +
  geom_qq(size = 2.35, alpha = 0.85) +
  geom_qq_line(color = "red", linewidth = 0.65) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", title = NULL) +
  theme_nature()

p_resfit <- ggplot(diag_df, aes(x = fitted, y = residual)) +
  geom_point(size = 2.35, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.65) +
  labs(x = "Fitted values", y = "Residuals", title = NULL) +
  theme_nature()

p_diag <- p_qq + p_resfit +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_svg(p_diag, "CFU_diagnostics.svg", "diagnostics")

anova_cfu <- anova(m_cfu)
eta_cfu_val <- anova_cfu$`Sum Sq`[1] / sum(anova_cfu$`Sum Sq`)
kw_cfu <- kruskal.test(log10_cfu ~ zone_label, data = cfu_df)

cat("CFU tests\n")
cat("ANOVA F =", round(anova_cfu$`F value`[1], 3), " p =", round(anova_cfu$`Pr(>F)`[1], 4), "\n")
cat("Eta2 =", round(eta_cfu_val, 3), "\n")
cat("Kruskal-Wallis H =", round(kw_cfu$statistic, 3), " p =", format.pval(kw_cfu$p.value, digits = 4), "\n\n")

pw_cfu <- pairwise.t.test(cfu_df$log10_cfu, cfu_df$zone_label, p.adjust.method = "BH", pool.sd = FALSE)

cat("Pairwise Welch t-tests (BH-adjusted)\n")
print(pw_cfu$p.value)
cat("\n")

pairs <- combn(zone_order, 2, simplify = FALSE)
es_cfu <- map_dfr(pairs, function(p) {
  sub <- cfu_df %>%
    filter(zone_label %in% p) %>%
    mutate(zone_label = droplevels(zone_label))
  es <- extract_es(cohens_d(log10_cfu ~ zone_label, data = sub, hedges.correction = TRUE))
  tibble(contrast = paste(p[1], "vs", p[2]), hedges_g = es$g, ci_low = es$ci_low, ci_high = es$ci_high)
})

cat("Pairwise Hedges' g\n")
print(as.data.frame(es_cfu))
cat("\n")

aggregate_zone_matrix <- function(morph_df, zone_code) {
  mt_cols <- names(morph_df)[str_detect(names(morph_df), "^[etd]_mt")]
  morph_df %>%
    group_by(replicate) %>%
    summarise(across(all_of(mt_cols), sum), .groups = "drop") %>%
    mutate(
      zone = zone_code,
      sample_id = paste0(zone_code, "-R", replicate)
    ) %>%
    select(sample_id, zone, replicate, all_of(mt_cols))
}

comm_E <- aggregate_zone_matrix(morph_E, "E")
comm_T <- aggregate_zone_matrix(morph_T, "T")
comm_D <- aggregate_zone_matrix(morph_D, "D")

cat("Community matrices (biological replicate level)\n")
cat("E:", nrow(comm_E), "replicates x", sum(str_detect(names(comm_E), "^e_mt")), "morphotypes\n")
cat("T:", nrow(comm_T), "replicates x", sum(str_detect(names(comm_T), "^t_mt")), "morphotypes\n")
cat("D:", nrow(comm_D), "replicates x", sum(str_detect(names(comm_D), "^d_mt")), "morphotypes\n\n")

compute_alpha <- function(comm_df, zone_code) {
  mt_cols <- names(comm_df)[str_detect(names(comm_df), "^[etd]_mt")]
  mat <- as.matrix(comm_df[, mt_cols])
  rownames(mat) <- comm_df$sample_id
  tibble(
    sample_id = comm_df$sample_id,
    zone = zone_code,
    replicate = comm_df$replicate,
    richness = specnumber(mat),
    shannon = diversity(mat, index = "shannon"),
    simpson = diversity(mat, index = "simpson"),
    N_total = rowSums(mat)
  )
}

alpha_df <- bind_rows(
  compute_alpha(comm_E, "E"),
  compute_alpha(comm_T, "T"),
  compute_alpha(comm_D, "D")
) %>%
  mutate(zone_label = factor(zone_codes[zone], levels = zone_order))

cat("Alpha diversity summary\n")
alpha_df %>%
  group_by(zone_label) %>%
  summarise(
    across(c(richness, shannon, simpson, N_total), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  print(width = Inf)
cat("\n")

test_alpha <- function(metric_name, alpha_data) {
  f <- reformulate("zone_label", response = metric_name)
  m <- lm(f, data = alpha_data)
  sw <- shapiro.test(residuals(m))
  lev <- leveneTest(f, data = alpha_data)
  aov_result <- anova(m)
  kw <- kruskal.test(f, data = alpha_data)
  
  cat(metric_name, "\n")
  cat("  Shapiro-Wilk W =", round(sw$statistic, 4), " p =", round(sw$p.value, 4), "\n")
  cat("  Levene F =", round(lev$`F value`[1], 3), " p =", round(lev$`Pr(>F)`[1], 4), "\n")
  cat("  ANOVA F =", round(aov_result$`F value`[1], 3), " p =", round(aov_result$`Pr(>F)`[1], 4), "\n")
  cat("  Kruskal-Wallis H =", round(kw$statistic, 3), " p =", format.pval(kw$p.value, digits = 4), "\n\n")
  
  list(sw = sw, lev = lev, aov = aov_result, kw = kw)
}

tests_rich <- test_alpha("richness", alpha_df)
tests_sha  <- test_alpha("shannon", alpha_df)
tests_sim  <- test_alpha("simpson", alpha_df)

run_rarefaction <- function(comm_df) {
  mt_cols <- names(comm_df)[str_detect(names(comm_df), "^[etd]_mt")]
  mat <- as.matrix(comm_df[, mt_cols])
  rownames(mat) <- comm_df$sample_id
  inext_list <- as.list(as.data.frame(t(mat)))
  names(inext_list) <- rownames(mat)
  iNEXT(inext_list, q = 0, datatype = "abundance", conf = 0.95, nboot = 200)
}

inext_E <- run_rarefaction(comm_E)
inext_T <- run_rarefaction(comm_T)
inext_D <- run_rarefaction(comm_D)

extract_inext <- function(inext_out, zone_code) {
  if (is.data.frame(inext_out$iNextEst$size_based)) {
    df <- inext_out$iNextEst$size_based %>%
      as_tibble() %>%
      rename(sample_id = Assemblage, x = m, y = qD, y.lwr = qD.LCL, y.upr = qD.UCL)
  } else {
    df <- bind_rows(
      lapply(names(inext_out$iNextEst$size_based), function(nm) {
        inext_out$iNextEst$size_based[[nm]] %>%
          as_tibble() %>%
          mutate(sample_id = nm)
      })
    ) %>%
      rename(x = m, y = qD, y.lwr = qD.LCL, y.upr = qD.UCL)
  }
  
  df %>%
    mutate(
      zone = zone_code,
      zone_label = factor(zone_codes[zone_code], levels = zone_order)
    )
}

rare_df <- bind_rows(
  extract_inext(inext_E, "E"),
  extract_inext(inext_T, "T"),
  extract_inext(inext_D, "D")
) %>%
  mutate(Method = as.character(Method), sample_id = as.character(sample_id))

p_rare <- ggplot(rare_df, aes(x = x, y = y, color = zone_label)) +
  geom_ribbon(
    aes(ymin = y.lwr, ymax = y.upr, fill = zone_label, group = sample_id),
    alpha = 0.06,
    color = NA
  ) +
  geom_line(
    data = rare_df %>% filter(Method == "Rarefaction"),
    aes(group = sample_id),
    linewidth = 0.56,
    linetype = "solid"
  ) +
  geom_line(
    data = rare_df %>% filter(Method == "Extrapolation"),
    aes(group = sample_id),
    linewidth = 0.56,
    linetype = "dashed"
  ) +
  geom_point(
    data = rare_df %>% filter(Method == "Observed"),
    aes(group = sample_id),
    size = 2.35,
    shape = 16
  ) +
  scale_color_manual(values = palette_3z) +
  scale_fill_manual(values = palette_3z, guide = "none") +
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

build_crosszone_matrix <- function(comm_zone, crosswalk_df, zone_col) {
  clusters_use <- crosswalk_df %>%
    filter(!is.na(.data[[zone_col]]) & .data[[zone_col]] != "")
  
  mt_cols <- names(comm_zone)[str_detect(names(comm_zone), "^[etd]_mt")]
  result <- tibble(sample_id = comm_zone$sample_id)
  
  for (i in seq_len(nrow(clusters_use))) {
    cl_id <- clusters_use$cluster_id[i]
    mt_codes <- str_trim(unlist(str_split(clusters_use[[zone_col]][i], ",")))
    col_names <- tolower(str_replace_all(mt_codes, "-", "_"))
    valid_cols <- intersect(col_names, mt_cols)
    if (length(valid_cols) > 0) {
      result[[cl_id]] <- rowSums(comm_zone[, valid_cols, drop = FALSE])
    }
  }
  
  result
}

cross_E <- build_crosszone_matrix(comm_E, crosswalk, "e_code")
cross_T <- build_crosszone_matrix(comm_T, crosswalk, "t_code")
cross_D <- build_crosszone_matrix(comm_D, crosswalk, "d_code")

all_clusters <- sort(unique(c(
  setdiff(names(cross_E), "sample_id"),
  setdiff(names(cross_T), "sample_id"),
  setdiff(names(cross_D), "sample_id")
)))

pad_clusters <- function(df, clusters) {
  for (cl in clusters) {
    if (!cl %in% names(df)) df[[cl]] <- 0L
  }
  df %>% select(sample_id, all_of(clusters))
}

cross_E <- pad_clusters(cross_E, all_clusters)
cross_T <- pad_clusters(cross_T, all_clusters)
cross_D <- pad_clusters(cross_D, all_clusters)

cross_comm <- bind_rows(cross_E, cross_T, cross_D) %>%
  column_to_rownames("sample_id")

cross_mat <- as.matrix(cross_comm)
stopifnot(all(cross_mat >= 0))

cat("Cross-zone community matrix\n")
cat("Dimensions:", nrow(cross_mat), "samples x", ncol(cross_mat), "clusters\n")
cat("Cluster IDs:", paste(colnames(cross_mat), collapse = ", "), "\n\n")

cross_hell <- decostand(cross_mat, method = "hellinger")

cross_meta <- tibble(
  sample_id  = rownames(cross_hell),
  zone       = str_extract(rownames(cross_hell), "^[ETD]"),
  zone_label = factor(zone_codes[zone], levels = zone_order)
)

bc_cross <- vegdist(cross_hell, method = "bray")

set.seed(42)
nmds_cross <- metaMDS(cross_hell, distance = "bray", k = 2, trymax = 200, autotransform = FALSE, trace = 0)

cat("Cross-zone NMDS stress =", round(nmds_cross$stress, 4), "\n\n")

scores_nmds <- as.data.frame(scores(nmds_cross, display = "sites")) %>%
  rownames_to_column("sample_id") %>%
  left_join(cross_meta, by = "sample_id")

p_nmds <- ggplot(scores_nmds, aes(x = NMDS1, y = NMDS2, color = zone_label)) +
  geom_point(size = 3.5, alpha = 0.9) +
  stat_ellipse(
    aes(fill = zone_label),
    geom = "polygon",
    alpha = 0.10,
    level = 0.68,
    linetype = "dashed",
    linewidth = 0.42
  ) +
  scale_color_manual(values = palette_3z) +
  scale_fill_manual(values = palette_3z, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
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

save_svg(p_nmds, "NMDS_crosszone.svg", "ordination")

pcoa_cross <- cmdscale(bc_cross, eig = TRUE, k = 2)
eig_pos <- pcoa_cross$eig[pcoa_cross$eig > 0]
var_exp <- 100 * pcoa_cross$eig[1:2] / sum(eig_pos)

scores_pcoa <- tibble(
  sample_id = rownames(pcoa_cross$points),
  PC1 = pcoa_cross$points[, 1],
  PC2 = pcoa_cross$points[, 2]
) %>%
  left_join(cross_meta, by = "sample_id")

p_pcoa <- ggplot(scores_pcoa, aes(x = PC1, y = PC2, color = zone_label)) +
  geom_point(size = 3.5, alpha = 0.9) +
  stat_ellipse(
    aes(fill = zone_label),
    geom = "polygon",
    alpha = 0.10,
    level = 0.68,
    linetype = "dashed",
    linewidth = 0.42
  ) +
  scale_color_manual(values = palette_3z) +
  scale_fill_manual(values = palette_3z, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
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

save_svg(p_pcoa, "PCoA_crosszone.svg", "ordination")

bc_mat <- as.matrix(bc_cross)
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
  geom_text(aes(label = round(bray, 2)), size = 2.35, color = "black") +
  scale_fill_viridis_c(option = "D", direction = -1, name = "Bray-Curtis") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 1
  )

save_svg(p_heatmap, "BrayCurtis_heatmap_crosszone.svg", "heatmap")

cat("PERMANOVA: cross-zone community ~ zone\n")
set.seed(42)
perm_cross <- adonis2(
  cross_hell ~ zone_label,
  data = cross_meta,
  permutations = 999,
  method = "bray",
  by = "margin"
)
print(perm_cross)
cat("\n")

cat("Pairwise PERMANOVA (BH-adjusted)\n")
pw_perm <- pairwise_adonis2(bc_cross, cross_meta$zone_label, permutations = 999, p_adjust_method = "BH")
print(as.data.frame(pw_perm))
cat("\n")

cat("Betadisper\n")
bd_cross <- betadisper(bc_cross, cross_meta$zone_label)
print(round(tapply(bd_cross$distances, cross_meta$zone_label, mean), 4))
cat("\n")

bd_anova <- anova(bd_cross)
set.seed(42)
bd_perm <- permutest(bd_cross, permutations = 999)
bd_perm_F <- bd_perm$tab[1, "F"]
bd_perm_p <- bd_perm$tab[1, "Pr(>F)"]

cat("Betadisper ANOVA F =", round(bd_anova$`F value`[1], 3), " p =", round(bd_anova$`Pr(>F)`[1], 4), "\n")
cat("Betadisper permutation F =", round(bd_perm_F, 3), " p =", round(bd_perm_p, 4), "\n\n")

cat("SIMPER: cluster contributions to between-zone dissimilarity\n")
set.seed(42)
sim_cross <- simper(cross_mat, cross_meta$zone_label, permutations = 999)

sim_summary_list <- summary(sim_cross)
simper_all <- map_dfr(names(sim_summary_list), function(comp_name) {
  ss <- sim_summary_list[[comp_name]]
  groups <- strsplit(comp_name, "_")[[1]]
  tibble(
    comparison = comp_name,
    cluster = rownames(ss),
    contribution = ss$average,
    sd = ss$sd,
    cumulative = ss$cumsum,
    p_value = ss$p,
    mean_A = ss$ava,
    mean_B = ss$avb,
    group_A = groups[1],
    group_B = groups[2]
  )
})

cat("Top contributors per comparison\n")
simper_all %>%
  group_by(comparison) %>%
  slice_head(n = 5) %>%
  print(n = 20, width = Inf)
cat("\n")

write_csv(simper_all, file.path(out_dir, "SIMPER_crosszone.csv"))

simper_plot_df <- simper_all %>%
  group_by(comparison) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    comparison = factor(comparison, levels = unique(comparison)),
    cluster_plot = factor(
      paste(comparison, cluster, sep = "___"),
      levels = rev(paste(comparison, cluster, sep = "___"))
    )
  )

p_simper <- ggplot(simper_plot_df, aes(x = cluster_plot, y = contribution)) +
  geom_col(aes(fill = cumulative <= 0.70), width = 0.7, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = pmax(contribution - sd, 0), ymax = contribution + sd),
    width = 0.2,
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = paste0(round(100 * cumulative, 1), "%")),
    hjust = -0.25,
    size = 2.6
  ) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B")) +
  scale_x_discrete(labels = function(x) str_replace(x, ".*___", "")) +
  coord_flip() +
  facet_wrap(~ comparison, scales = "free_y", ncol = 1) +
  labs(
    x = NULL,
    y = "Mean contribution to Bray-Curtis dissimilarity",
    title = NULL
  ) +
  theme_nature() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

save_svg(p_simper, "SIMPER_crosszone_contributions.svg", "simper")

make_comp_long <- function(comm_df, zone_code) {
  mt_cols <- names(comm_df)[str_detect(names(comm_df), "^[etd]_mt")]
  comm_df %>%
    select(sample_id, all_of(mt_cols)) %>%
    pivot_longer(cols = all_of(mt_cols), names_to = "morphotype", values_to = "count") %>%
    mutate(
      morphotype = toupper(str_replace_all(morphotype, "_", "-")),
      zone = zone_code,
      zone_label = factor(zone_codes[zone_code], levels = zone_order)
    ) %>%
    group_by(sample_id) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
}

comp_long <- bind_rows(
  make_comp_long(comm_E, "E"),
  make_comp_long(comm_T, "T"),
  make_comp_long(comm_D, "D")
)

sample_order <- alpha_df %>% arrange(zone_label, N_total) %>% pull(sample_id)
comp_long <- comp_long %>%
  mutate(sample_id = factor(sample_id, levels = sample_order))

make_comp_plot <- function(data, zone_name) {
  zone_data <- data %>% filter(zone_label == zone_name)
  n_mt <- n_distinct(zone_data$morphotype)
  
  mt_order <- zone_data %>%
    group_by(morphotype) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(morphotype)
  
  zone_data <- zone_data %>%
    mutate(morphotype = factor(morphotype, levels = rev(mt_order)))
  
  pal <- scales::hue_pal()(n_mt)
  names(pal) <- mt_order
  
  ggplot(zone_data, aes(x = sample_id, y = proportion, fill = morphotype)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = pal, name = "Morphotype") +
    labs(x = NULL, y = "Relative abundance", title = NULL) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.8),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 7.0)
    )
}

p_comp_E <- make_comp_plot(comp_long, "Entrance")
p_comp_T <- make_comp_plot(comp_long, "Transition")
p_comp_D <- make_comp_plot(comp_long, "Dark")

save_svg(p_comp_E, "composition_Entrance.svg", "composition_zone")
save_svg(p_comp_T, "composition_Transition.svg", "composition_zone")
save_svg(p_comp_D, "composition_Dark.svg", "composition_zone")

cross_comp_long <- as.data.frame(cross_mat) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "cluster", values_to = "count") %>%
  mutate(
    zone = str_extract(sample_id, "^[ETD]"),
    zone_label = factor(zone_codes[zone], levels = zone_order)
  ) %>%
  group_by(sample_id) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

cross_sample_order <- cross_meta %>%
  left_join(
    tibble(sample_id = rownames(cross_mat), N_total = rowSums(cross_mat)),
    by = "sample_id"
  ) %>%
  arrange(zone_label, N_total) %>%
  pull(sample_id)

cross_cluster_order <- cross_comp_long %>%
  group_by(cluster) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(cluster)

cross_comp_long <- cross_comp_long %>%
  mutate(
    sample_id = factor(sample_id, levels = cross_sample_order),
    cluster = factor(cluster, levels = rev(cross_cluster_order))
  )

cross_palette <- scales::hue_pal()(length(cross_cluster_order))
names(cross_palette) <- cross_cluster_order

p_comp_cross <- ggplot(cross_comp_long, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_col(width = 0.8) +
  facet_wrap(~ zone_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = cross_palette, name = "Cluster") +
  labs(x = NULL, y = "Relative abundance", title = NULL) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.6),
    legend.key.size = unit(0.38, "cm"),
    legend.text = element_text(size = 6.8)
  )

save_svg(p_comp_cross, "composition_crosswalk_clusters.svg", "composition_cross")

rank_df <- comp_long %>%
  group_by(zone_label, morphotype) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop") %>%
  group_by(zone_label) %>%
  arrange(desc(mean_prop)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p_rank <- ggplot(rank_df, aes(x = rank, y = mean_prop, color = zone_label)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.25) +
  scale_color_manual(values = palette_3z) +
  scale_y_log10(labels = scales::label_number()) +
  labs(x = "Morphotype rank", y = "Mean relative abundance (log scale)", color = "Zone", title = NULL) +
  theme_nature()

save_svg(p_rank, "rank_abundance.svg", "rank")

p_cfu <- ggplot(cfu_df, aes(x = zone_label, y = log10_cfu, color = zone_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
  geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
  scale_color_manual(values = palette_3z) +
  labs(x = NULL, y = expression(log[10]~"(CFU g"^{-1}~"sediment)"), title = NULL) +
  theme_nature() +
  theme(legend.position = "none")

save_svg(p_cfu, "CFU_by_zone.svg", "cfu_zone")

make_alpha_plot <- function(data, metric, ylab) {
  ggplot(data, aes(x = zone_label, y = .data[[metric]], color = zone_label)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
    geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
    scale_color_manual(values = palette_3z) +
    labs(x = NULL, y = ylab, title = NULL) +
    theme_nature() +
    theme(legend.position = "none")
}

p_alpha_rich <- make_alpha_plot(alpha_df, "richness", "Morphotype richness")
p_alpha_sha  <- make_alpha_plot(alpha_df, "shannon", "Shannon (H')")
p_alpha_sim  <- make_alpha_plot(alpha_df, "simpson", "Simpson (1 − D)")

p_alpha <- p_alpha_rich + p_alpha_sha + p_alpha_sim +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_svg(p_alpha, "alpha_diversity_panel.svg", "alpha_panel")

p_mass <- ggplot(mass, aes(x = zone_label, y = mass_g, color = zone_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.40) +
  geom_jitter(width = 0.1, size = 2.8, alpha = 0.82) +
  scale_color_manual(values = palette_3z) +
  labs(x = NULL, y = "Sediment mass (g)", title = NULL) +
  theme_nature() +
  theme(legend.position = "none")

save_svg(p_mass, "sediment_mass_by_zone.svg", "mass")

crosswalk_viz <- crosswalk %>%
  filter(!str_detect(cluster_name, "(?i)fungal")) %>%
  select(cluster_id, cluster_name, confidence, e_code, t_code, d_code) %>%
  pivot_longer(cols = c(e_code, t_code, d_code), names_to = "zone_col", values_to = "mt_code") %>%
  mutate(
    zone_label = case_when(
      zone_col == "e_code" ~ "Entrance",
      zone_col == "t_code" ~ "Transition",
      zone_col == "d_code" ~ "Dark"
    ),
    zone_label = factor(zone_label, levels = zone_order),
    present = !is.na(mt_code) & mt_code != "",
    cluster_label = paste0(cluster_id, " · ", cluster_name)
  )

cl_order <- crosswalk_viz %>%
  distinct(cluster_label, confidence) %>%
  mutate(conf_rank = case_when(
    confidence == "High" ~ 1,
    confidence == "Moderate" ~ 2,
    TRUE ~ 3
  )) %>%
  arrange(conf_rank, cluster_label) %>%
  pull(cluster_label)

crosswalk_viz <- crosswalk_viz %>%
  mutate(cluster_label = factor(cluster_label, levels = rev(cl_order)))

p_crosswalk <- ggplot(crosswalk_viz, aes(x = zone_label, y = cluster_label)) +
  geom_tile(
    aes(fill = interaction(present, confidence)),
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    data = crosswalk_viz %>% filter(present),
    aes(label = mt_code),
    size = 2.4,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE.High" = "#B5D4F4",
      "TRUE.Moderate" = "#D6E8F7",
      "TRUE.Low" = "#EBF2F9",
      "FALSE.High" = "grey97",
      "FALSE.Moderate" = "grey97",
      "FALSE.Low" = "grey97"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_nature() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7)
  )

save_svg(p_crosswalk, "crosswalk_visualization.svg", "crosswalk")

fungal_cols_E <- c("e_mt13", "e_mt14", "e_mt15")
fungal_cols_T <- c("t_mt13", "t_mt14")
fungal_cols_D <- c("d_mt13")

compute_alpha_nofungi <- function(comm_df, fungal_cols) {
  mt_cols <- names(comm_df)[str_detect(names(comm_df), "^[etd]_mt")]
  bact_cols <- setdiff(mt_cols, fungal_cols)
  mat <- as.matrix(comm_df[, bact_cols])
  rownames(mat) <- comm_df$sample_id
  tibble(
    sample_id = comm_df$sample_id,
    richness_nofungi = specnumber(mat),
    shannon_nofungi = diversity(mat, index = "shannon")
  )
}

alpha_nf <- bind_rows(
  compute_alpha_nofungi(comm_E, fungal_cols_E) %>% mutate(zone = "E"),
  compute_alpha_nofungi(comm_T, fungal_cols_T) %>% mutate(zone = "T"),
  compute_alpha_nofungi(comm_D, fungal_cols_D) %>% mutate(zone = "D")
) %>%
  left_join(alpha_df %>% select(sample_id, richness, shannon), by = "sample_id") %>%
  mutate(
    delta_richness = richness - richness_nofungi,
    delta_shannon = shannon - shannon_nofungi
  )

cat("Sensitivity: fungal exclusion effect\n")
alpha_nf %>%
  group_by(zone) %>%
  summarise(mean_delta_rich = mean(delta_richness), mean_delta_sha = round(mean(delta_shannon), 3), .groups = "drop") %>%
  print()
cat("\n")

write_csv(alpha_nf, file.path(out_dir, "sensitivity_fungal_exclusion.csv"))

cat("Cave water control\n")
print(control)

cfu_control <- control %>%
  mutate(
    cfu_per_mL = colony_count / (dilution_factor * (volume_plated_u_l / 1000))
  )

cat("\nEstimated CFU/mL of cave water\n")
print(cfu_control %>% select(dilution, colony_count, cfu_per_mL))
cat("\n")

analysis_df <- cfu_df %>%
  select(
    sample_id = plate_id, zone, zone_label, replicate,
    dilution, colony_count, mass_g, mass_aliquot_g, cfu_per_g, log10_cfu
  ) %>%
  mutate(sample_id = paste0(zone, "-R", replicate)) %>%
  left_join(alpha_df %>% select(sample_id, richness, shannon, simpson, N_total), by = "sample_id")

write_csv(analysis_df, file.path(out_dir, "current_analysis_table.csv"))

report_file <- file.path(out_dir, "current_stats_report.txt")
sink(report_file)

cat("Organal San Antonio — Current Experiment\n")
cat("Statistical report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("Design\n")
cat("Zones: Entrance (E), Transition (T), Dark (D)\n")
cat("Replicates: 5 per zone, 15 total\n")
cat("Medium: 25% R2A agar + NaHCO3\n")
cat("Incubation: 15 C, darkness, aerobic\n")
cat("Colony counts: Day 7\n")
cat("Morphotypes: E = 15 (3 fungal), T = 14 (2 fungal), D = 13 (1 fungal)\n")
cat("CFU normalisation: estimated field-moist mass of the 200 uL sediment aliquot\n\n")

cat("CFU summary\n")
print(as.data.frame(cfu_summary), row.names = FALSE)
cat("\n")

cat("CFU diagnostics\n")
cat("Shapiro-Wilk W =", round(sw_cfu$statistic, 4), " p =", round(sw_cfu$p.value, 4), "\n")
cat("Levene F =", round(lev_cfu$`F value`[1], 4), " p =", round(lev_cfu$`Pr(>F)`[1], 4), "\n")
cat("ANOVA F =", round(anova_cfu$`F value`[1], 3), " p =", round(anova_cfu$`Pr(>F)`[1], 4), "\n")
cat("Eta2 =", round(eta_cfu_val, 3), "\n")
cat("Kruskal-Wallis H =", round(kw_cfu$statistic, 3), " p =", format.pval(kw_cfu$p.value, digits = 4), "\n\n")

cat("Pairwise Welch t-tests (BH-adjusted)\n")
print(pw_cfu$p.value)
cat("\n")

cat("Pairwise Hedges' g\n")
print(as.data.frame(es_cfu), row.names = FALSE)
cat("\n")

cat("Alpha diversity\n")
for (metric_info in list(
  list(name = "richness", t = tests_rich),
  list(name = "shannon",  t = tests_sha),
  list(name = "simpson",  t = tests_sim)
)) {
  tt <- metric_info$t
  cat(metric_info$name, ":\n")
  cat("  Shapiro-Wilk: p =", round(tt$sw$p.value, 4), "\n")
  cat("  Levene: F =", round(tt$lev$`F value`[1], 3), " p =", round(tt$lev$`Pr(>F)`[1], 4), "\n")
  cat("  ANOVA: F =", round(tt$aov$`F value`[1], 3), " p =", round(tt$aov$`Pr(>F)`[1], 4), "\n")
  cat("  Kruskal-Wallis: H =", round(tt$kw$statistic, 3), " p =", format.pval(tt$kw$p.value, digits = 4), "\n\n")
}

cat("Beta diversity\n")
cat("NMDS stress =", round(nmds_cross$stress, 4), "\n\n")
print(perm_cross)
cat("\n")
print(as.data.frame(pw_perm), row.names = FALSE)
cat("\n")
cat("Betadisper ANOVA F =", round(bd_anova$`F value`[1], 3), " p =", round(bd_anova$`Pr(>F)`[1], 4), "\n")
cat("Betadisper permutation F =", round(bd_perm_F, 3), " p =", round(bd_perm_p, 4), "\n\n")

cat("SIMPER top contributors per comparison\n")
print(
  as.data.frame(
    simper_all %>%
      group_by(comparison) %>%
      slice_head(n = 5) %>%
      ungroup()
  ),
  row.names = FALSE
)
cat("\n")

cat("QC\n")
cat("Mean delta:", round(mean(qc_all$delta), 1), "\n")
cat("Delta range:", min(qc_all$delta), "to", max(qc_all$delta), "\n")
cat("Plates with |delta| > 10:", sum(abs(qc_all$delta) > 10), "\n")
cat("Mean |tech diff|:", round(mean(tech_pairs$diff, na.rm = TRUE), 1), "\n")
cat("Max |tech diff|:", max(tech_pairs$diff, na.rm = TRUE), "\n\n")

cat("Cave water control CFU/mL\n")
print(as.data.frame(cfu_control %>% select(dilution, colony_count, cfu_per_mL)), row.names = FALSE)
cat("\n")

sink()

cat("\nOutput directory\n")
cat(out_dir, "\n\n")
cat("Files generated\n")
list.files(out_dir) %>% walk(~ cat(.x, "\n"))

sink(file.path(out_dir, "session_info.txt"))
cat("Session info — Current experiment analysis\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
sessionInfo()
sink()

cat("\nCurrent experiment analysis complete\n")
cat("Outputs saved in:", out_dir, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# Composite figures for manuscript
# ──────────────────────────────────────────────────────────────────────────────

density_diversity_composite <- (
  p_cfu + p_alpha_rich + p_alpha_sha + p_alpha_sim
) +
  plot_layout(ncol = 2, widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_svg(density_diversity_composite, "density_diversity_composite.svg", "density_diversity_composite")
save_pdf(density_diversity_composite, "density_diversity_composite.pdf", "density_diversity_composite")

composition_composite <- (
  (p_nmds | p_pcoa) / p_comp_cross
) +
  plot_layout(heights = c(1.0, 1.45)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_svg(composition_composite, "composition_composite.svg", "composition_composite")
save_pdf(composition_composite, "composition_composite.pdf", "composition_composite")