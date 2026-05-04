<p align="left">
  <img src="assets/eafit_logo.svg" alt="Universidad EAFIT" height="48"/>
</p>

# Cultivable Microbial Community Structure Along a Light Gradient in a Tropical Volcaniclastic Cave

**Sebastian Correa-Gallego** · B.Sc. Thesis in Biology · Universidad EAFIT · Medellín, Colombia · 2026

Advisor: **Dr. Nicolás Pinel Peláez**

---

## Overview

This repository contains the data, R scripts, statistical outputs, and figures associated with the
undergraduate thesis *Cultivable Microbial Community Structure Along a Light Gradient in a
Tropical Volcaniclastic Cave* (Correa-Gallego, 2026).

The study characterizes the aerobic, low-temperature, dilute-R2A cultivable microbial fraction
of sediments collected along the light-defined ecological gradient of the Organal San Antonio —
a tropical volcaniclastic pseudokarstic cave at approximately 2,350 m a.s.l. in Támesis,
Antioquia, Colombia. Two independent cultivation experiments were conducted: an initial
binary Light–Dark contrast and a primary three-sector experiment resolving Entrance,
Transition, and Dark zones. Community structure was assessed through cultivable density,
morphotype-level alpha diversity, and multivariate compositional analysis using a 23-cluster
operational phenotypic crosswalk.

The primary results show that sector identity along the cave gradient is associated with a
roughly sixfold reduction in cultivable density in the aphotic Dark sector relative to illuminated
sectors (*η*² = 0.652), strong compositional turnover across zones (PERMANOVA *R*² = 0.573,
*p* = 0.001) under balanced multivariate dispersion, and no detectable sector-level differences
in local morphotype-cluster richness or evenness. The Transition sector supported a
distinctive set of zone-restricted phenotypic clusters, indicating that it is not a simple
intermediate between Entrance and Dark but an ecologically distinct component of the system.

To our knowledge, this constitutes the first spatially resolved cultivable microbial baseline
for an organal-type tropical volcaniclastic cave.

---

## Repository Structure

```
psychro-oligophiles/
│
├── data/
│   ├── raw/                 # Original, unmodified field and laboratory records
│   └── processed/           # Analysis-ready tables derived by the scripts
│
├── scripts/                 # R scripts for Experiment I and Experiment II
│
├── outputs/
│   ├── figures/
│   │   ├── main/            # Main-text figures (Figures 1–6)
│   │   └── supplementary/   # Supplementary figures (Figures S1–S4)
│   └── stats/               # Statistical summary reports
│
├── docs/                    # Compiled thesis PDF
│
└── assets/                  # Repository assets (logo)
```

---

## Data

| File | Description |
|------|-------------|
| `data/raw/colony_counts_exp1.csv` | Day-7 colony counts for all plates, Experiment I |
| `data/raw/colony_counts_exp2.csv` | Day-7 colony counts for all plates, Experiment II |
| `data/raw/sediment_mass.csv` | Field-moist sediment mass per replicate, Experiment II |
| `data/raw/morphotype_proportions_exp1.csv` | Per-plate morphotype abundance matrix, Experiment I |
| `data/raw/morphotype_proportions_exp2.csv` | Per-plate morphotype abundance matrix, Experiment II |
| `data/raw/metadata_env.csv` | *In situ* abiotic conditions recorded during sampling |
| `data/processed/pilot_analysis_table.csv` | Integrated density and diversity table, Experiment I |
| `data/processed/current_analysis_table.csv` | Integrated density and diversity table, Experiment II |
| `data/processed/morphotype_crosswalk_v3.csv` | Final 23-cluster cross-sector phenotypic crosswalk |

---

## Scripts

All analyses were conducted in **R v4.5.2**. Key packages: `vegan` v2.7-3, `iNEXT` v3.0.2,
`effectsize` v1.0.2, `tidyverse`. Full session information is recorded in
`scripts/session_info.txt`.

| Script | Description |
|--------|-------------|
| `scripts/experiment_1.R` | CFU estimation, alpha and beta diversity, PERMANOVA, and figures for Experiment I |
| `scripts/experiment_2.R` | CFU estimation, crosswalk integration, alpha and beta diversity, PERMANOVA, SIMPER, rarefaction, and figures for Experiment II |

Scripts are self-contained and can be run independently. Both scripts source data from
`data/raw/` and `data/processed/` and write figures to `outputs/figures/`.

---

## Methods Summary

- **Study site:** Organal San Antonio, Támesis, Antioquia, Colombia (5.664°N, 75.729°W; ~2,350 m a.s.l.)
- **Medium:** 25% R2A agar + 0.1 g/L NaHCO₃; 15 g/L bacteriological agar
- **Incubation:** Experiment I, ~20 °C ambient; Experiment II, 14.5–15 °C in total darkness
- **Colony census:** Day 7
- **Density normalization:** Experiment I, volumetric (CFU mL⁻¹ wet sediment); Experiment II, aliquot-level gravimetric (CFU g⁻¹ field-moist sediment)
- **Community analysis:** Hellinger-transformed Bray–Curtis dissimilarities; NMDS; PERMANOVA (*adonis2*, 999 permutations); betadisper; SIMPER
- **Diversity:** Morphotype richness, Shannon *H*′, Simpson 1 − *D*; rarefaction via iNEXT
- **Effect sizes:** Hedges' *g* (pairwise), *η*² (ANOVA), *R*² (PERMANOVA); 95% confidence intervals reported throughout

---

## Permits

Field collection and mobilization of biological material were conducted under the collection
permit granted to Universidad EAFIT by the Autoridad Nacional de Licencias Ambientales
(ANLA), Resolution No. 000548 of 28 March 2025.

---

## Citation

> Correa-Gallego, S. (2026). *Cultivable Microbial Community Structure Along a Light Gradient
> in a Tropical Volcaniclastic Cave*. B.Sc. Thesis in Biology, Universidad EAFIT, Medellín,
> Colombia.

---

## License

© 2026 Sebastian Correa-Gallego.
This work is licensed under a
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/).