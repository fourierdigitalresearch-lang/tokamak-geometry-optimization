"""
Temperature Scan for Optimal Geometry
Scans T0 from 8 to 25 keV for A=2.0, κ=1.8, δ=-0.10
"""

import numpy as np
import csv
import time

from scan_extended_grid import (
    R0, B0, TARGET_BETAN, ALPHA_P, ALPHA_T,
    calibrate_pressure, compute_q95, fusion_power
)

# ============================================================
# PARAMETERS
# ============================================================
A = 2.0
KAPPA = 1.8
DELTA = -0.10
a = R0 / A

T0_VALUES = [8, 10, 12, 14, 16, 20, 25]
NR, NZ = 180, 220


def main():
    print("=" * 70)
    print("TEMPERATURE SCAN - OPTIMAL GEOMETRY")
    print(f"Geometry: A={A}, κ={KAPPA}, δ={DELTA}")
    print(f"Mesh: {NR}×{NZ}")
    print("=" * 70)
    
    print("\nCalibrating equilibrium...")
    eq = calibrate_pressure(TARGET_BETAN, R0, a, KAPPA, DELTA, B0,
                           alpha_p=ALPHA_P, nr=NR, nz=NZ)
    
    print(f"  p0 = {eq['p0']/1e6:.4f} MPa")
    print(f"  Ip = {eq['Ip']:.2f} MA")
    print(f"  βN = {eq['beta_N']:.4f}")
    
    results = []
    
    print("\nScanning temperatures:")
    for T0 in T0_VALUES:
        print(f"\nT0 = {T0} keV")
        
        P_fus, P_wall = fusion_power(eq, T0, alpha_T=ALPHA_T)
        q95 = compute_q95(eq)
        
        results.append({
            "T0": T0,
            "Pfus": P_fus,
            "Pwall": P_wall,
            "q95": q95,
            "Ip": eq["Ip"],
            "beta_N": eq["beta_N"]
        })
        
        print(f"  → Pfus = {P_fus:.1f} MW, Pwall = {P_wall:.3f} MW/m²")
    
    peak = max(results, key=lambda x: x["Pfus"])
    print("\n" + "=" * 70)
    print(f"PEAK: T0 = {peak['T0']} keV, Pfus = {peak['Pfus']:.1f} MW")
    print("=" * 70)
    
    with open("temperature_scan_final.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["T0", "Pfus", "Pwall", "q95", "Ip", "beta_N"])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults saved to: temperature_scan_final.csv")


if __name__ == "__main__":
    main()