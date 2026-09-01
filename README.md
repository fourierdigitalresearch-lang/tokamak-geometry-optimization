**Tokamak Geometry Optimization - High-βN Study**

## Project information

Preprint title:

**“Optimized Plasma Geometry for a High-βN Tokamak: 2.79 GW Fusion Power Under Simultaneous Engineering and Stability Constraints”**

Author: **Jean Lauro Muller**
Affiliation: **Fourier Digital Research**
Contact: **[fourierdigitalresearch@gmail.com](mailto:fourierdigitalresearch@gmail.com)**

Zenodo DOI:

**10.5281/zenodo.22231424**

Include the standard Zenodo DOI badge:

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22231424.svg)](https://doi.org/10.5281/zenodo.22231424)
```

## Key optimized configuration

Present the main parameters clearly, preferably in a compact Markdown table:

* Aspect Ratio: **A = 2.0**
* Elongation: **κ = 1.8**
* Triangularity: **δ = -0.10**
* Major Radius: **R₀ = 6.2 m**
* Minor Radius: **a = 3.1 m**
* Toroidal Magnetic Field: **B₀ = 5.3 T**

Explain briefly that the optimum occurs at the boundary/vertex of the imposed engineering and stability constraints:

* **A = 2.0**: lowest scanned aspect ratio satisfying the imposed current constraint at the working resolution.
* **κ = 1.8**: saturates the imposed vertical-stability ceiling.
* **δ = -0.10**: negative triangularity gives the highest fusion power within the scanned parameter range.

Do not describe this as a globally proven optimum beyond the scanned parameter space.

## Performance results

Present the following values explicitly as **working-resolution estimates from the 180×220 mesh**:

* Fusion Power: **P_fus = 2.79 GW**
* Plasma Current: **I_p = 20.22 MA**
* Safety Factor: **q₉₅ = 5.015**
* Wall Loading: **P_wall = 2.53 MW/m²**
* Central Ion Temperature: **T₀ = 16 keV**
* Normalized Beta: **β_N = 2.785**

Make clear that these are numerical estimates at the stated working resolution and should not be interpreted as asymptotically mesh-converged values.

## Repository structure

Use exactly this repository structure:

```text
tokamak-geometry-optimization/
├── scan_extended_grid.py
├── temperature_scan_final.py
├── mesh_convergence.py
├── README.md
├── LICENSE
└── data/
    ├── scan_with_restrictions.csv
    ├── ranking_valid.csv
    ├── pareto_valid.csv
    ├── temperature_scan_final.csv
    └── mesh_convergence.csv
```

Briefly describe each script:

### `scan_extended_grid.py`

Performs the extended geometry parameter scan over:

* **A ∈ [1.8, 3.4]**
* **κ ∈ [1.5, 2.4]**
* **δ ∈ [-0.10, 0.35]**

The scan contains **80 geometries**.

Outputs:

* `scan_with_restrictions.csv`: complete scan and constraint information.
* `ranking_valid.csv`: valid configurations satisfying all imposed constraints.
* `pareto_valid.csv`: Pareto-valid configurations.

### `temperature_scan_final.py`

Scans the central ion temperature **T₀** for the selected geometry.

Output:

* `temperature_scan_final.csv`

### `mesh_convergence.py`

Evaluates the numerical behavior of the solution across multiple mesh resolutions.

The documented resolutions range from **60×80 to 180×220**, with an additional **240×300** run discussed separately in the convergence caveat.

Output:

* `mesh_convergence.csv`

## Installation and requirements

State:

* Python **3.8 or newer**
* NumPy
* SciPy
* Matplotlib, optional for plotting

Provide:

```bash
pip install numpy scipy matplotlib
```

Then provide concise commands:

```bash
python scan_extended_grid.py
python temperature_scan_final.py
python mesh_convergence.py
```

## Physical and engineering constraints

Clearly document the three principal constraints:

1. **Vertical stability:** κ ≤ 1.8
2. **Engineering current limit:** I_p ≤ 20 MA
3. **MHD stability:** q₉₅ ≥ 3.0

Explain that the optimization is performed within the scanned geometry space subject to these constraints.

Avoid implying that these constraints represent a complete reactor engineering design assessment.

## Numerical method

Include a concise technical description:

### Grad-Shafranov solver

* Second-order centered finite differences
* SciPy sparse linear solver
* Successive under-relaxation with **ω = 0.6**
* Convergence tolerance **10⁻⁷**
* Mask-based rectangular computational grid
* Miller parametrization for plasma boundary geometry

### Physics models

Include:

* D-T fusion reactivity using the **Bosch-Hale parametrization**
* Pressure profile:

```text
p(ψ_N) = p₀(1 − ψ_N)^0.5
```

* Temperature profile:

```text
T(ψ_N) = T₀(1 − ψ_N)^0.5
```

* Vacuum toroidal-field approximation:

```text
F(ψ) = R₀B₀ = constant
```

Keep the description factual and do not add physical models that are not specified.

## IMPORTANT: Numerical convergence caveat

Create a prominent section titled:

**Numerical Convergence Caveat**

State clearly:

The **180×220 mesh is the working baseline**, but convergence is not clean enough to claim an asymptotically converged solution.

The relative changes in fusion power between successive mesh resolutions are:

* **+9.3%**
* **+14.0%**
* **+4.1%**

These changes do not form a monotonically decreasing sequence.

An additional **240×300** calculation was attempted. The solver itself converged within the prescribed tolerance, but boundary-contour quantities became numerically unreliable. In particular:

* **q₉₅ increased by approximately 31%**
* **l_i decreased by approximately 64%**

The likely cause is the staircase representation of the plasma boundary associated with the mask-based rectangular grid.

Therefore:

* The **240×300 run is not used to revise the reported baseline results**.
* **P_fus = 2.79 GW** and **I_p = 20.22 MA** must be described as **working-resolution estimates**.
* They should **not** be presented as fully mesh-converged asymptotic values.
* Readers should consult the associated preprint, particularly **Sections 3.2 and 4.3**, for the complete discussion.

This caveat is important and must not be hidden or minimized.

## Interpretation of the optimization

Include a short section explaining that the scan identifies a high-performance configuration at the intersection of the imposed constraints.

Use careful scientific language:

* The result is an optimum **within the explored parameter space and model assumptions**.
* It is not a proof of a globally optimal tokamak configuration.
* The numerical values depend on the adopted profiles, physics models, boundary representation, constraints, and mesh resolution.

Do not add claims about commercial viability, net electric power, economic competitiveness, or reactor feasibility unless explicitly supported by the supplied information.

## Citation

Include the following BibTeX entry exactly:

```bibtex
@misc{muller2026tokamak,
  author = {Muller, Jean Lauro},
  title = {Optimized Plasma Geometry for a High-$\beta_N$ Tokamak: 2.79 GW Fusion Power Under Simultaneous Engineering and Stability Constraints},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.22231424},
  url = {https://doi.org/10.5281/zenodo.22231424}
}
```

Also include a short sentence asking users who employ the code or data in research to cite the work.

## License

State:

**MIT License**

and link to:

```markdown
[LICENSE](LICENSE)
```

## Acknowledgments

Include:

**Computational Resources: Fourier Digital Research**

## Links

Include a concise Links section containing:

* Zenodo DOI
* Preprint
* GitHub repository

Do not invent a preprint URL or GitHub URL. If the exact URL is not supplied, use a neutral placeholder such as:

```text
Preprint: Zenodo record
GitHub Repository: this repository
```

## Style requirements

The final README must:

* Be written entirely in **English**.
* Use professional scientific English.
* Be concise rather than verbose.
* Use Markdown headings and tables effectively.
* Use code blocks for commands, formulas, and BibTeX.
* Avoid excessive emojis.
* Avoid marketing language.
* Avoid unsupported claims.
* Clearly distinguish **numerical results**, **model assumptions**, and **engineering constraints**.
* Preserve the convergence caveat prominently.
* Be suitable for direct publication as `README.md`.

At the end include:

**Last Updated: September 2026**

Return **only the complete Markdown content of `README.md`**, beginning with the repository title.
