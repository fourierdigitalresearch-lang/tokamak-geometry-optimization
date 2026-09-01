[README (1).md](https://github.com/user-attachments/files/31690730/README.1.md)[Uploading# Tokamak Geometry Optimization - High-βN Study

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

<!-- Replace zenodo.XXXXXXX above and the badge/DOI below once the Zenodo record is minted -->

## Description

This repository contains the numerical scripts and data supporting the preprint:

**"Optimized Plasma Geometry for a High-βN Tokamak: 2.79 GW Fusion Power Under Simultaneous Engineering and Stability Constraints"**

**Author:** Jean Lauro Muller
**Affiliation:** Fourier Digital Research
**Contact:** contact@fourierresearch.org

## Key Results

### Optimal Configuration
- **Aspect Ratio (A):** 2.0
- **Elongation (κ):** 1.8
- **Triangularity (δ):** -0.10 (negative)
- **Major Radius (R₀):** 6.2 m
- **Minor Radius (a):** 3.1 m
- **Toroidal Field (B₀):** 5.3 T

### Performance (working-resolution estimate, 180×220 mesh -- see mesh-convergence caveat below)
- **Fusion Power (P_fus):** 2.79 GW
- **Plasma Current (I_p):** 20.22 MA
- **Safety Factor (q₉₅):** 5.015
- **Wall Loading (P_wall):** 2.53 MW/m²
- **Ion Temperature (T₀):** 16 keV
- **Normalized Beta (β_N):** 2.785

## Repository Structure

```
tokamak-geometry-optimization/
├── scan_extended_grid.py          # Full parameter scan (80 geometries)
├── temperature_scan_final.py      # Temperature scan (7 simulations)
├── mesh_convergence.py            # Mesh convergence study (4 resolutions)
├── README.md                      # This file
├── LICENSE                        # MIT License
└── data/
    ├── scan_with_restrictions.csv
    ├── ranking_valid.csv
    ├── pareto_valid.csv
    ├── temperature_scan_final.csv
    └── mesh_convergence.csv
```

## Requirements

- **Python:** 3.8 or higher
- **NumPy:** `pip install numpy`
- **SciPy:** `pip install scipy`
- **Matplotlib:** `pip install matplotlib` (optional, for plotting)

## Usage

### 1. Extended Grid Scan

Performs a comprehensive scan over aspect ratio (A ∈ [1.8, 3.4]), elongation (κ ∈ [1.5, 2.4]), and triangularity (δ ∈ [-0.1, 0.35]):

```bash
python scan_extended_grid.py
```

**Outputs:**
- `scan_with_restrictions.csv` - All 80 geometries with constraints
- `ranking_valid.csv` - 6 valid geometries satisfying all constraints

### 2. Temperature Scan

Scans central ion temperature (T₀) for the optimal geometry:

```bash
python temperature_scan_final.py
```

**Outputs:**
- `temperature_scan_final.csv` - Fusion power vs. temperature data

### 3. Mesh Convergence Study

Validates numerical convergence across multiple mesh resolutions:

```bash
python mesh_convergence.py
```

**Outputs:**
- `mesh_convergence.csv` - Convergence data for 4 mesh resolutions (60×80 to 180×220)

## Main Findings

### Constraints Applied
1. **Vertical Stability:** κ ≤ 1.8
2. **Engineering Current Limit:** I_p ≤ 20 MA (Nb₃Sn TF coils)
3. **MHD Stability:** q₉₅ ≥ 3.0

### Optimal Geometry
The optimization reveals that the optimal configuration lies at the **vertex of the constraint space**:
- **A = 2.0** (smallest aspect ratio respecting I_p ≤ 20 MA at this mesh resolution)
- **κ = 1.8** (saturates vertical stability ceiling)
- **δ = -0.10** (negative triangularity maximizes fusion power within the scanned range)

### Mesh Convergence -- read before citing P_fus or I_p
- The 180×220 mesh is used as the working baseline, but convergence is **not** clean:
  relative changes in P_fus between successive meshes are +9.3%, +14.0%, +4.1% -- not a
  monotonically decreasing sequence.
- An additional 240×300 run was attempted: the solver converged well within tolerance,
  but boundary-contour quantities (q₉₅, l_i) became numerically unreliable (q₉₅ jumped
  +31%, l_i collapsed -64%), most likely due to the staircase boundary representation of
  the mask-based rectangular grid. That run is not used to revise the reported baseline.
- Treat P_fus = 2.79 GW and I_p = 20.22 MA as working-resolution estimates, not converged
  asymptotic values. See the preprint, Section 3.2/4.3, for full discussion.

## Numerical Method

### Grad-Shafranov Solver
- **Discretization:** Second-order centered finite differences
- **Solver:** SciPy sparse linear solver
- **Acceleration:** Successive under-relaxation (ω = 0.6)
- **Convergence tolerance:** 10⁻⁷
- **Boundary representation:** Mask-based rectangular grid (Miller parametrization)

### Physics Models
- **D-T Reactivity:** Bosch-Hale parametrization
- **Pressure profile:** p(ψ_N) = p₀(1 - ψ_N)^0.5
- **Temperature profile:** T(ψ_N) = T₀(1 - ψ_N)^0.5
- **Vacuum toroidal field approximation:** F(ψ) = R₀B₀ = const

## Citation

If you use this code or data in your research, please cite:

```bibtex
@misc{muller2026tokamak,
  author = {Muller, Jean Lauro},
  title = {Optimized Plasma Geometry for a High-$\beta_N$ Tokamak: 2.79 GW Fusion Power Under Simultaneous Engineering and Stability Constraints},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.XXXXXXX},
  url = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Computational Resources:** Fourier Digital Research

## Links

- **Preprint:** Zenodo record *(update link when available)*
- **GitHub Repository:** *(confirm exact organization/URL before publishing)*

---

**Last Updated:** September 2026
 README (1).md…]()
