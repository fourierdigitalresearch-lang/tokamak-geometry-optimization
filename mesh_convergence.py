"""
Mesh Convergence Study
Tests resolutions: 60×80, 80×100, 120×150, 180×220, 240×300
"""

import numpy as np
import csv
import time

from scan_extended_grid import (
    R0, B0, TARGET_BETAN, ALPHA_P, ALPHA_T, T0_FIXED,
    calibrate_pressure, compute_q95, fusion_power
)
from scipy.interpolate import RectBivariateSpline

# ============================================================
# PARAMETERS
# ============================================================
A = 2.0
KAPPA = 1.8
DELTA = -0.10
a = R0 / A

MESHES = [
    (60, 80),
    (80, 100),
    (120, 150),
    (180, 220),
    (240, 300)  # Diagnostic run
]


def compute_li(equilibrium):
    """Compute internal inductance"""
    R, Z = equilibrium["R"], equilibrium["Z"]
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
    psi_target = psi_axis * (1.0 - 0.98)
    
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
    bp2_edge_vals = interp_Bp2.ev(rr, zz)
    bp_edge = np.mean(np.maximum(bp2_edge_vals, 0.0))
    
    if bp_edge <= 0:
        return np.nan
    
    return equilibrium["bp_avg"] / bp_edge


def main():
    print("=" * 70)
    print("MESH CONVERGENCE STUDY")
    print(f"Geometry: A={A}, κ={KAPPA}, δ={DELTA}, T0={T0_FIXED} keV")
    print("=" * 70)
    
    results = []
    
    for nr, nz in MESHES:
        print(f"\nMesh: {nr}×{nz}")
        t0 = time.time()
        
        try:
            eq = calibrate_pressure(TARGET_BETAN, R0, a, KAPPA, DELTA, B0,
                                   alpha_p=ALPHA_P, nr=nr, nz=nz)
            
            q95 = compute_q95(eq)
            li = compute_li(eq)
            P_fus, P_wall = fusion_power(eq, T0_FIXED, alpha_T=ALPHA_T)
            
            dt = time.time() - t0
            
            results.append({
                "nr": nr,
                "nz": nz,
                "Pfus": P_fus,
                "q95": q95,
                "Ip": eq["Ip"],
                "beta_N": eq["beta_N"],
                "li": li,
                "iterations": eq["iterations"],
                "time": dt
            })
            
            print(f"  Pfus = {P_fus:.1f} MW")
            print(f"  Ip = {eq['Ip']:.2f} MA")
            print(f"  q95 = {q95:.3f}")
            # CORREÇÃO: f-string válida. np.nan formata naturalmente como 'nan'
            print(f"  li = {li:.3f}")
            print(f"  βN = {eq['beta_N']:.4f}")
            print(f"  Iterations: {eq['iterations']}")
            print(f"  Time: {dt:.1f} s")
            
        except Exception as e:
            print(f"  FAILED: {e}")
            continue
    
    with open("mesh_convergence.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["nr", "nz", "Pfus", "q95", 
                                                "Ip", "beta_N", "li", "iterations", "time"])
        writer.writeheader()
        writer.writerows(results)
    
    print("\n" + "=" * 70)
    print("CONVERGENCE ANALYSIS")
    print("=" * 70)
    
    if len(results) >= 2:
        for i in range(1, len(results)):
            prev = results[i-1]
            curr = results[i]
            
            delta_P = (curr["Pfus"] - prev["Pfus"]) / prev["Pfus"] * 100
            delta_Ip = (curr["Ip"] - prev["Ip"]) / prev["Ip"] * 100
            
            print(f"\n{prev['nr']}×{prev['nz']} → {curr['nr']}×{curr['nz']}:")
            print(f"  ΔPfus = {delta_P:+.1f}%")
            print(f"  ΔIp = {delta_Ip:+.1f}%")
    
    print(f"\nResults saved to: mesh_convergence.csv")
    print("=" * 70)


if __name__ == "__main__":
    main()