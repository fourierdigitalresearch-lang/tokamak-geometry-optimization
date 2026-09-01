# Tokamak Geometry Optimization - High-βN Study

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22231424.svg)](https://doi.org/10.5281/zenodo.22231424)

## Project Information

- **Preprint Title**: Optimized Plasma Geometry for a High-βN Tokamak: 2.79 GW Fusion Power Under Simultaneous Engineering and Stability Constraints
- **Author**: Jean Lauro Muller
- **Affiliation**: Fourier Digital Research
- **Contact**: fourierdigitalresearch@gmail.com
- **Zenodo DOI**: 10.5281/zenodo.22231424

## Key Optimized Configuration

The optimization identifies a high-performance configuration at the boundary/vertex of the imposed engineering and stability constraints.

| Parameter | Value |
| :--- | :--- |
| Aspect Ratio ($A$) | 2.0 |
| Elongation ($\kappa$) | 1.8 |
| Triangularity ($\delta$) | -0.10 |
| Major Radius ($R_0$) | 6.2 m |
| Minor Radius ($a$) | 3.1 m |
| Toroidal Magnetic Field ($B_0$) | 5.3 T |

**Context of the Optimum:**
- $A = 2.0$: The lowest scanned aspect ratio satisfying the imposed current constraint at the working resolution.
- $\kappa = 1.8$: Saturates the imposed vertical-stability ceiling.
- $\delta = -0.10$: Negative triangularity yields the highest fusion power within the scanned parameter range.

*Note: This result represents an optimum within the explored parameter space and model assumptions. It is not a proof of a globally optimal tokamak configuration.*

## Performance Results

The following values are explicitly presented as **working-resolution estimates** derived from the $180 \times 220$ mesh. They should not be interpreted as asymptotically mesh-converged values.

| Metric | Value |
| :--- | :--- |
| Fusion Power ($P_{\text{fus}}$) | 2.79 GW |
| Plasma Current ($I_p$) | 20.22 MA |
| Safety Factor ($q_{95}$) | 5.015 |
| Wall Loading ($P_{\text{wall}}$) | 2.53 MW/m² |
| Central Ion Temperature ($T_0$) | 16 keV |
| Normalized Beta ($\beta_N$) | 2.785 |

## Repository Structure

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

### Script Descriptions

- **`scan_extended_grid.py`**: Performs the extended geometry parameter scan over $A \in [1.8, 3.4]$, $\kappa \in [1.5, 2.4]$, and $\delta \in [-0.10, 0.35]$. The scan contains 80 geometries.
  - *Outputs*: `scan_with_restrictions.csv` (complete scan and constraint information), `ranking_valid.csv` (valid configurations satisfying all imposed constraints), `pareto_valid.csv` (Pareto-valid configurations).
- **`temperature_scan_final.py`**: Scans the central ion temperature $T_0$ for the selected optimal geometry.
  - *Output*: `temperature_scan_final.csv`.
- **`mesh_convergence.py`**: Evaluates the numerical behavior of the solution across multiple mesh resolutions. The documented resolutions range from $60 \times 80$ to $180 \times 220$, with an additional $240 \times 300$ run discussed separately in the convergence caveat.
  - *Output*: `mesh_convergence.csv`.

## Installation and Requirements

- Python 3.8 or newer
- NumPy
- SciPy
- Matplotlib (optional, for plotting)

Install dependencies via:
```bash
pip install numpy scipy matplotlib
```

Execute the scripts using the following commands:
```bash
python scan_extended_grid.py
python temperature_scan_final.py
python mesh_convergence.py
```

## Physical and Engineering Constraints

The optimization is performed within the scanned geometry space subject to three principal constraints. These constraints are baseline viability criteria and do not represent a complete reactor engineering design assessment:

1. **Vertical stability**: $\kappa \le 1.8$
2. **Engineering current limit**: $I_p \le 20$ MA
3. **MHD stability**: $q_{95} \ge 3.0$

## Numerical Method

- Grad-Shafranov solver
- Second-order centered finite differences
- SciPy sparse linear solver
- Successive under-relaxation with $\omega = 0.6$
- Convergence tolerance $10^{-7}$
- Mask-based rectangular computational grid
- Miller parametrization for plasma boundary geometry

## Physics Models

- D-T fusion reactivity using the Bosch-Hale parametrization
- Pressure profile: $p(\psi_N) = p_0(1 - \psi_N)^{0.5}$
- Temperature profile: $T(\psi_N) = T_0(1 - \psi_N)^{0.5}$
- Vacuum toroidal-field approximation: $F(\psi) = R_0B_0 = \text{constant}$

*The description above is strictly factual and does not include physical models beyond those specified in the associated preprint.*

## ⚠️ IMPORTANT: Numerical Convergence Caveat

The $180 \times 220$ mesh serves as the working baseline; however, convergence is not clean enough to claim an asymptotically converged solution. 

The relative changes in fusion power between successive mesh resolutions are:
- **+9.3%**
- **+14.0%**
- **+4.1%**

These changes do not form a monotonically decreasing sequence. An additional $240 \times 300$ calculation was attempted. While the solver itself converged within the prescribed tolerance, boundary-contour quantities became numerically unreliable. In particular:
- $q_{95}$ increased by approximately **31%**
- $l_i$ decreased by approximately **64%**

The likely cause is the staircase representation of the plasma boundary associated with the mask-based rectangular grid, which perturbs volume integrals mildly but affects contour-tracing quantities severely.

**Therefore:**
1. The $240 \times 300$ run is not used to revise the reported baseline results.
2. $P_{\text{fus}} = 2.79$ GW and $I_p = 20.22$ MA must be described strictly as working-resolution estimates.
3. They should not be presented as fully mesh-converged asymptotic values.

Readers should consult the associated preprint, particularly Sections 3.2 and 4.3, for the complete discussion. This caveat is fundamental to the interpretation of the results and must not be hidden or minimized.

## Interpretation of the Optimization

The numerical scan identifies a high-performance configuration at the intersection of the imposed constraints. It is critical to note that the numerical values depend directly on the adopted profiles, physics models, boundary representation, constraints, and mesh resolution. No claims are made regarding commercial viability, net electric power, economic competitiveness, or absolute reactor feasibility, as these require comprehensive systems-level analysis beyond the scope of this geometric scoping study.

## Citation

If you employ this code or data in your research, please cite this work using the following BibTeX entry:

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

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Computational Resources**: Fourier Digital Research

## Links

- **Zenodo DOI**: [10.5281/zenodo.22231424](https://doi.org/10.5281/zenodo.22231424)
- **Preprint**: Zenodo record
- **GitHub Repository**: this repository

---
Last Updated: September 2026
