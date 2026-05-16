<p align="left">
  <img src="assets/eafit_logo.svg" alt="Universidad EAFIT" height="36"/>
</p>

# psychro-oligophiles

**Cultivable Microbial Community Structure Along a Light Gradient in a Tropical Volcaniclastic Cave**

Sebastian Correa-Gallego
B.Sc. Thesis in Biology, Universidad EAFIT, 2026  
Advisor: Dr. Nicolás Pinel Peláez

---

Sediments from the Organal San Antonio cave (Támesis, Antioquia, ~2,350 m a.s.l.) were cultured on dilute R2A at low temperature across three light-defined sectors: Entrance, Transition, and Dark. Two cultivation experiments characterized cultivable density, morphotype diversity, and community composition via a 23-cluster phenotypic crosswalk.

Cultivable density declined ~sixfold into the aphotic zone (η² = 0.652). Strong compositional turnover was detected across sectors (PERMANOVA *R*² = 0.573, *p* = 0.001), with the Transition sector forming a distinct assemblage rather than a gradient intermediate. No sector-level differences in alpha diversity were found.

---

## Repository

```
psychro-oligophiles/
├── data/
│   ├── raw/          # Field and laboratory records
│   └── processed/    # Analysis-ready tables
├── scripts/          # R scripts (Experiments I & II)
├── outputs/
│   ├── figures/      # SVG and PDF plots
│   └── stats/        # Summary statistics
└── docs/figures/     # Publication figures
```

Analyses were run in R v4.5.2 (`vegan`, `iNEXT`, `effectsize`, `tidyverse`). Scripts are self-contained; both source data from `data/` and write to `outputs/`. See `scripts/session_info.txt` for full session details.

Field collection was authorized under ANLA Resolution No. 000548 (28 March 2025), granted to Universidad EAFIT.

---

## Citation

Correa-Gallego, S. (2026). *Cultivable microbial community structure along a light gradient in a tropical volcaniclastic cave* [Undergraduate thesis, Universidad EAFIT]. Repositorio institucional EAFIT.

---

## License

Code (`scripts/`): [MIT](LICENSE) · Data (`data/`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) · Manuscript: [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/)
