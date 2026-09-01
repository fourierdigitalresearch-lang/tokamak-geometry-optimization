"""
Extended Grid Scan for Tokamak Geometry Optimization
Performs comprehensive parameter scan over A, κ, δ (80 combinations)
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from matplotlib.path import Path
import csv
import time

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================
MU0 = 4.0e-7 * np.pi
KEV_J = 1.602176634e-16
E_DT = 17.59e6 * 1.602176634e-19

# ============================================================
# DEFAULT PARAMETERS
# ============================================================
R0 = 6.2
B0 = 5.3
TARGET_BETAN = 2.8
ALPHA_P = 0.5
ALPHA_T = 0.5
T0_FIXED = 16.0

# Mesh settings
NR, NZ = 80, 100
MAX_ITER = 150
RELAX = 0.6
TOL = 1e-7
MAX_BISECT = 10
EPS = 1e-8


def boundary_curve(R0, a, kappa, delta, n=500):
    """Generate plasma boundary using Miller parametrization"""
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    Rb = R0 + a * np.cos(theta + delta * np.sin(theta))
    Zb = kappa * a * np.sin(theta)
    return Rb, Zb


def build_operator(R0, a, kappa, delta, nr, nz):
    """Build Grad-Shafranov operator with mask-based boundary"""
    Rb, Zb = boundary_curve(R0, a, kappa, delta)
    R = np.linspace(Rb.min(), Rb.max(), nr)
    Z = np.linspace(Zb.min(), Zb.max(), nz)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    
    polygon = Path(np.column_stack([Rb, Zb]))
    points = np.column_stack([RR.ravel(), ZZ.ravel()])
    mask = polygon.contains_points(points).reshape(RR.shape)
    
    idx = -np.ones(mask.shape, dtype=int)
    idx[mask] = np.arange(mask.sum())
    
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]
    
    rows, cols, data = [], [], []
    for i, j in zip(*np.where(mask)):
        k = idx[i, j]
        r = R[i]
        rows.append(k); cols.append(k)
        data.append(-2.0 / dR**2 - 2.0 / dZ**2)
        
        for ii, coeff in [(i + 1, 1.0 / dR**2 - 1.0 / (2.0 * r * dR)),
                          (i - 1, 1.0 / dR**2 + 1.0 / (2.0 * r * dR))]:
            if 0 <= ii < nr and mask[ii, j]:
                rows.append(k); cols.append(idx[ii, j]); data.append(coeff)
        
        for jj, coeff in [(j + 1, 1.0 / dZ**2), (j - 1, 1.0 / dZ**2)]:
            if 0 <= jj < nz and mask[i, jj]:
                rows.append(k); cols.append(idx[i, jj]); data.append(coeff)
    
    operator = sp.csr_matrix((data, (rows, cols)), shape=(mask.sum(), mask.sum()))
    return R, Z, RR, ZZ, mask, idx, operator, Rb, Zb


def solve_equilibrium(p0, R0, a, kappa, delta, B0, alpha_p=ALPHA_P,
                       nr=NR, nz=NZ, max_iter=MAX_ITER, relax=RELAX,
                       tol=TOL, eps=EPS):
    """Solve Grad-Shafranov equilibrium"""
    R, Z, RR, ZZ, mask, idx, operator, Rb, Zb = build_operator(R0, a, kappa, delta, nr, nz)
    F0 = R0 * B0
    psi_vec = np.zeros(mask.sum())
    
    for iteration in range(max_iter):
        psi = np.zeros(mask.shape)
        psi[mask] = psi_vec
        psi_axis = psi[mask].min()
        dpsi = max(-psi_axis, 1e-12)
        psi_n = np.clip((psi - psi_axis) / dpsi, 0.0, 1.0)
        psi_n_safe = np.clip(psi_n, 0.0, 1.0 - eps)
        
        dpdpsi = -p0 * alpha_p / dpsi * np.maximum(1.0 - psi_n_safe, eps) ** (alpha_p - 1.0)
        rhs = (-MU0 * RR**2 * dpdpsi)[mask]
        
        new_psi = spla.spsolve(operator, rhs)
        error = np.linalg.norm(new_psi - psi_vec) / max(np.linalg.norm(new_psi), 1e-14)
        
        psi_vec = relax * new_psi + (1.0 - relax) * psi_vec
        if error < tol:
            break
    
    psi = np.zeros(mask.shape)
    psi[mask] = psi_vec
    psi_axis = psi[mask].min()
    dpsi = max(-psi_axis, 1e-12)
    psi_n = np.clip((psi - psi_axis) / dpsi, 0.0, 1.0)
    psi_n_safe = np.clip(psi_n, 0.0, 1.0 - eps)
    
    p = p0 * np.maximum(1.0 - psi_n_safe, eps) ** alpha_p
    dpsi_dR = np.gradient(psi, R, axis=0)
    dpsi_dZ = np.gradient(psi, Z, axis=1)
    Bp2 = (dpsi_dR**2 + dpsi_dZ**2) / RR**2
    
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]
    volume_weight = 2.0 * np.pi * RR
    dV = dR * dZ
    
    V = np.sum(volume_weight[mask]) * dV
    p_avg = np.sum(p[mask] * volume_weight[mask]) * dV / V
    beta = 2.0 * MU0 * p_avg / B0**2
    
    dpdpsi_final = -p0 * alpha_p / dpsi * np.maximum(1.0 - psi_n_safe, eps) ** (alpha_p - 1.0)
    jphi = RR * dpdpsi_final
    Ip = np.sum(np.abs(jphi[mask])) * dR * dZ / 1e6
    
    beta_N = 100.0 * beta * a * B0 / Ip if Ip > 0 else np.inf
    bp_avg = np.sum(Bp2[mask] * volume_weight[mask]) * dV / V
    
    dRb = np.diff(Rb, append=Rb[0])
    dZb = np.diff(Zb, append=Zb[0])
    dl_wall = np.sqrt(dRb**2 + dZb**2)
    A_wall = np.sum(2.0 * np.pi * Rb * dl_wall)
    
    return {
        "R": R, "Z": Z, "RR": RR, "ZZ": ZZ, "mask": mask,
        "psi": psi, "psi_n": psi_n, "p": p, "Bp2": Bp2,
        "Ip": Ip, "beta_N": beta_N, "beta": beta, "bp_avg": bp_avg,
        "V": V, "A_wall": A_wall, "p0": p0, "F0": F0,
        "iterations": iteration + 1, "Rb": Rb, "Zb": Zb
    }


def calibrate_pressure(target_betaN, R0, a, kappa, delta, B0,
                        alpha_p=ALPHA_P, p_lo=2e4, p_hi=3e6,
                        max_bisect=MAX_BISECT, **kwargs):
    """Calibrate p0 via bisection to match target beta_N"""
    low = solve_equilibrium(p_lo, R0, a, kappa, delta, B0, alpha_p=alpha_p, **kwargs)
    high = solve_equilibrium(p_hi, R0, a, kappa, delta, B0, alpha_p=alpha_p, **kwargs)
    
    if not (low["beta_N"] <= target_betaN <= high["beta_N"]):
        return min([low, high], key=lambda x: abs(x["beta_N"] - target_betaN))
    
    best = None
    lo, hi = p_lo, p_hi
    
    for _ in range(max_bisect):
        mid = (lo + hi) / 2.0
        eq = solve_equilibrium(mid, R0, a, kappa, delta, B0, alpha_p=alpha_p, **kwargs)
        
        if best is None or abs(eq["beta_N"] - target_betaN) < abs(best["beta_N"] - target_betaN):
            best = eq
        
        if eq["beta_N"] < target_betaN:
            lo = mid
        else:
            hi = mid
    
    return best


def compute_q95(equilibrium):
    """Compute safety factor at psi_N = 0.95"""
    from scipy.interpolate import RectBivariateSpline
    
    R, Z = equilibrium["R"], equilibrium["Z"]
    F0 = equilibrium["F0"]
    Bp2 = equilibrium["Bp2"]
    interp_Bp2 = RectBivariateSpline(R, Z, Bp2, kx=3, ky=3)
    
    psi = equilibrium["psi"]
    mask = equilibrium["mask"]
    interp_psi = RectBivariateSpline(R, Z, psi, kx=3, ky=3)
    
    psi_masked = psi.copy()
    psi_masked[~mask] = np.nan
    idx_axis = np.nanargmin(psi_masked)
    i_axis, j_axis = np.unravel_index(idx_axis, psi.shape)
    R_axis, Z_axis = R[i_axis], Z[j_axis]
    psi_axis = psi[i_axis, j_axis]
    psi_target = psi_axis * (1.0 - 0.95)
    
    n_theta = 120
    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    rr, zz = [], []
    
    for th in theta:
        dr, dz = np.cos(th), np.sin(th)
        rho = 0.0
        step = 0.005 * max(np.ptp(R), np.ptp(Z))
        
        for _ in range(2000):
            rho += step
            R_test = R_axis + rho * dr
            Z_test = Z_axis + rho * dz
            
            if not (R.min() <= R_test <= R.max() and Z.min() <= Z_test <= Z.max()):
                break
            
            psi_val = interp_psi.ev(R_test, Z_test).item()
            if psi_val >= psi_target:
                R_prev = R_axis + (rho - step) * dr
                Z_prev = Z_axis + (rho - step) * dz
                psi_prev = interp_psi.ev(R_prev, Z_prev).item()
                
                if psi_prev != psi_val:
                    frac = (psi_target - psi_prev) / (psi_val - psi_prev)
                else:
                    frac = 0.0
                
                rho_exact = (rho - step) + frac * step
                rr.append(R_axis + rho_exact * dr)
                zz.append(Z_axis + rho_exact * dz)
                break
    
    if len(rr) < 8:
        return np.nan
    
    rr, zz = np.array(rr), np.array(zz)
    idx = np.argsort(np.arctan2(zz - Z_axis, rr - R_axis))
    rr, zz = rr[idx], zz[idx]
    
    bp2_vals = interp_Bp2.ev(rr, zz)
    bp = np.sqrt(np.maximum(bp2_vals, 1e-30))
    btor = F0 / rr
    dl = np.sqrt(np.diff(rr, append=rr[0])**2 + np.diff(zz, append=zz[0])**2)
    
    return np.sum(btor / (rr * bp) * dl) / (2.0 * np.pi)


def sigma_v_DT(T_keV):
    """Bosch-Hale D-T reactivity"""
    T = np.asarray(T_keV, dtype=float)
    BG = 34.3827
    MRC2 = 1124656.0
    C1, C2, C3, C4, C5, C6, C7 = (1.17302e-9, 1.51361e-2, 7.51886e-2,
                                  4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5)
    theta = T / (1.0 - T * (C2 + T * (C4 + T * C6)) / (1.0 + T * (C3 + T * (C5 + T * C7))))
    xi = (BG**2 / (4.0 * theta)) ** (1.0 / 3.0)
    sv_cm3_s = C1 * theta * np.sqrt(xi / (MRC2 * T**3)) * np.exp(-3.0 * xi)
    return sv_cm3_s * 1e-6


def fusion_power(equilibrium, T0_keV, alpha_T=ALPHA_T):
    """Compute fusion power"""
    psi_n = equilibrium["psi_n"]
    p = equilibrium["p"]
    R, Z = equilibrium["R"], equilibrium["Z"]
    RR, mask = equilibrium["RR"], equilibrium["mask"]
    
    fT = np.maximum(1.0 - psi_n, 0.0) ** alpha_T
    T_keV = T0_keV * fT
    T_J = np.maximum(T_keV, 1e-6) * KEV_J
    
    n_i = p / (2.0 * T_J)
    n_D = n_T = n_i / 2.0
    
    rate = n_D * n_T * sigma_v_DT(np.maximum(T_keV, 0.2))
    power_density = rate * E_DT
    
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]
    dV = 2.0 * np.pi * RR * dR * dZ
    
    P_fus = np.sum(power_density[mask] * dV[mask])
    P_wall = P_fus / equilibrium["A_wall"]
    
    return P_fus / 1e6, P_wall / 1e6


def main():
    """Main scan loop"""
    print("=" * 70)
    print("EXTENDED GRID SCAN - 80 GEOMETRIES")
    print("=" * 70)
    
    # CORREÇÃO: 5 * 4 * 4 = 80 combinações exatas
    A_values = np.linspace(1.8, 3.4, 5)
    kappa_values = np.linspace(1.5, 2.4, 4)
    delta_values = np.linspace(-0.1, 0.35, 4)
    
    results = []
    
    for A in A_values:
        a = R0 / A
        for kappa in kappa_values:
            for delta in delta_values:
                print(f"\nScanning: A={A:.2f}, κ={kappa:.2f}, δ={delta:.2f}")
                
                try:
                    eq = calibrate_pressure(TARGET_BETAN, R0, a, kappa, delta, B0,
                                           alpha_p=ALPHA_P, nr=NR, nz=NZ)
                    
                    q95 = compute_q95(eq)
                    P_fus, P_wall = fusion_power(eq, T0_FIXED, alpha_T=ALPHA_T)
                    
                    valid = (kappa <= 1.8 and eq["Ip"] <= 20.0 and 
                             q95 >= 3.0 and not np.isnan(q95))
                    
                    results.append({
                        "A": A, "kappa": kappa, "delta": delta,
                        "q95": q95, "Ip": eq["Ip"], "Pfus": P_fus,
                        "Pwall": P_wall, "beta_N": eq["beta_N"],
                        "valid": valid
                    })
                    
                    status = "VALID" if valid else "invalid"
                    print(f"  → Pfus={P_fus:.1f} MW, Ip={eq['Ip']:.2f} MA, q95={q95:.3f} [{status}]")
                    
                except Exception as e:
                    print(f"  → FAILED: {e}")
                    continue
    
    with open("scan_with_restrictions.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["A", "kappa", "delta", "q95", 
                                                "Ip", "Pfus", "Pwall", "beta_N", "valid"])
        writer.writeheader()
        writer.writerows(results)
    
    valid_results = [r for r in results if r["valid"]]
    with open("ranking_valid.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["A", "kappa", "delta", "q95", 
                                                "Ip", "Pfus", "Pwall", "beta_N"])
        writer.writeheader()
        writer.writerows(valid_results)
    
    print("\n" + "=" * 70)
    print(f"Scan complete: {len(results)} geometries, {len(valid_results)} valid")
    print("Results saved to: scan_with_restrictions.csv, ranking_valid.csv")
    print("=" * 70)


if __name__ == "__main__":
    main()