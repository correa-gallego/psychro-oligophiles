<p align="left">
  <img src="assets/eafit_logo.svg" alt="Universidad EAFIT" height="40"/>
</p>

# psychro-oligophiles

**Cultivable Microbial Community Structure Along a Light Gradient in a Tropical Volcaniclastic Cave**

Sebastian Correa-Gallego — B.Sc. Thesis in Biology · Universidad EAFIT · 2026
Advisor: Dr. Nicolás Pinel Peláez

---

## About

This repository contains the data, R scripts, and analytical outputs for the undergraduate thesis *Cultivable Microbial Community Structure Along a Light Gradient in a Tropical Volcaniclastic Cave* (Correa-Gallego, 2026).

The study characterizes the aerobic, low-temperature, dilute-R2A cultivable microbial fraction of sediments collected along the light-defined ecological gradient of the Organal San Antonio — a tropical volcaniclastic pseudokarstic cave at approximately 2,350 m a.s.l. in Támesis, Antioquia, Colombia. Two independent cultivation experiments were conducted: an initial binary Light–Dark contrast (Experiment I) and a primary three-sector experiment resolving Entrance, Transition, and Dark zones (Experiment II). Community structure was assessed through cultivable density, morphotype-level alpha diversity, and multivariate compositional analysis using a 23-cluster operational phenotypic crosswalk.

Primary results: sector identity is associated with a roughly sixfold reduction in cultivable density in the aphotic Dark sector (*η*² = 0.652), strong compositional turnover across zones (PERMANOVA *R*² = 0.573, *p* = 0.001) under balanced multivariate dispersion, and no detectable sector-level differences in local morphotype-cluster richness or evenness. The Transition sector supported a distinctive assemblage of zone-restricted phenotypic clusters, indicating it is not a compositional intermediate between Entrance and Dark but an ecologically distinct component of the gradient.

---

## Repository Structure

```
psychro-oligophiles/
├── data/
│   ├── raw/                          # Original, unmodified field and laboratory records
│   └── processed/                    # Analysis-ready tables derived by the scripts
├── scripts/                          # R scripts for Experiment I and Experiment II
├── outputs/
│   ├── figures/
│   │   ├── main/                     # Main-text figures (Figures 1–6)
│   │   └── supplementary/            # Supplementary figures (Figures S1–S4)
│   └── stats/                        # Statistical summary reports
└── assets/                           # Repository assets
```

The compiled thesis will be deposited in the Universidad EAFIT institutional repository upon approval and linked here.

---

## Data

**Raw**

| File | Description |
|------|-------------|
| `data/raw/colony_counts_day7.csv` | Day-7 colony counts for all plates, Experiment II |
| `data/raw/counts_dilution_e2.csv` | Dilution-series counts, Experiment I |
| `data/raw/sediment_mass.csv` | Field-moist sediment mass per replicate, Experiment II |
| `data/raw/morphotype_matrix_E.csv` | Per-plate morphotype abundance matrix, Entrance sector |
| `data/raw/morphotype_matrix_T.csv` | Per-plate morphotype abundance matrix, Transition sector |
| `data/raw/morphotype_matrix_D.csv` | Per-plate morphotype abundance matrix, Dark sector |
| `data/raw/morphotype_counts_pilot.csv` | Morphotype abundance matrix, Experiment I |
| `data/raw/morphotypes.csv` | Morphotype catalogue, Experiment I |
| `data/raw/cave_water_control.csv` | Fresh cave-water environmental control counts |
| `data/raw/metadata_env.csv` | *In situ* abiotic conditions recorded during sampling |

**Processed**

| File | Description |
|------|-------------|
| `data/processed/current_analysis_table.csv` | Integrated density and diversity table, Experiment II |
| `data/processed/pilot_analysis_table.csv` | Integrated density and diversity table, Experiment I |
| `data/processed/morphotype_crosswalk_v3.csv` | Final 23-cluster cross-sector phenotypic crosswalk |

---

## Scripts

All analyses were conducted in R v4.5.2. Key packages: `vegan` v2.7-3, `iNEXT` v3.0.2, `effectsize` v1.0.2, `tidyverse`. Full session information is in `scripts/session_info.txt`.

| Script | Description |
|--------|-------------|
| `scripts/experiment_1.R` | CFU estimation, alpha and beta diversity, PERMANOVA, and figures for Experiment I |
| `scripts/experiment_2.R` | CFU estimation, crosswalk integration, alpha and beta diversity, PERMANOVA, SIMPER, rarefaction, and figures for Experiment II |

Scripts are self-contained and can be run independently. Both source data from `data/` and write figures to `outputs/figures/`.

---

## Methods

| Parameter | Experiment I | Experiment II |
|-----------|-------------|---------------|
| Design | Light / Dark, 6 replicates per zone | Entrance / Transition / Dark, 5 replicates per zone |
| Medium | 25% R2A + 0.1 g/L NaHCO₃ | 25% R2A + 0.1 g/L NaHCO₃ |
| Incubation | ~20 °C, ambient | 14.5–15 °C, total darkness |
| Colony census | Day 7 | Day 7 |
| Density normalization | Volumetric (CFU mL⁻¹ wet sediment) | Aliquot-level gravimetric (CFU g⁻¹ field-moist sediment) |
| Community units | 14 global morphotypes | 23-cluster phenotypic crosswalk |

Community analysis: Hellinger-transformed Bray–Curtis dissimilarities; NMDS; PERMANOVA (`adonis2`, 999 permutations); homogeneity of multivariate dispersion (`betadisper`); SIMPER. Effect sizes (*η*², Hedges' *g*, PERMANOVA *R*²) with 95% confidence intervals are reported as primary inferential evidence; *p*-values are auxiliary.

---

## Permits

Field collection and mobilization of biological material were conducted under the collection permit granted to Universidad EAFIT by the Autoridad Nacional de Licencias Ambientales (ANLA), Resolution No. 000548 of 28 March 2025.

---

## Citation

Correa-Gallego, S. (2026). *Cultivable microbial community structure along a light gradient in a tropical volcaniclastic cave* [Undergraduate thesis, Universidad EAFIT]. Repositorio institucional EAFIT.

---

## License

**Code** (`scripts/`): [MIT License](LICENSE)

**Data** (`data/`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

**Thesis manuscript**: © 2026 Sebastian Correa-Gallego. Licensed under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/). The manuscript is a separate work from this repository.