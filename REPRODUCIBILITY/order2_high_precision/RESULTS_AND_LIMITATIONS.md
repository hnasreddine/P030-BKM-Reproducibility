# P030 — Direct microscopic verification of the second BKM horizon for transfer order m=2

**Date:** 21 November 2026. **Model:** the actual Hermitian one-excitation continuum discretization `bath_poly(2, gs=1, gf=3)` from the P030 V0.7 reproducibility package. **No replacement of the microscopic propagation with a two-pole evolution.**

## Result

For the physical test bath, gamma_s=1, gamma_f=3, m=2, the proposed threshold is **c_BKM=4**, not 8 (8 corresponds to rates 1 and 2). Recomputed BKM defects at macroscopic time tau=c log(1/lambda), microscopic physical time T=tau/lambda²:

| lambda | c | Bath modes N | BKM defect, direct microscopic | Existing corrected-pole prediction | Relative discrepancy |
|---|---:|---:|---:|---:|---:|
| 0.05 | 4.5 | 8,500 | 8.626903479e-39 | 8.626903650e-39 | 0.000% |
| 0.05 | 5 | 8,500 | 3.969035535e-41 | 3.968736229e-41 | 0.008% |
| 0.05 | 5 | 10,500 | 3.969035535e-41 | 3.968736229e-41 | 0.008% |
| 0.05 | 5.5 | 8,500 | 2.074247545e-42 | 2.074084554e-42 | 0.008% |
| 0.035 | 5 | 17,000 | 2.988999640e-45 | 2.989237080e-45 | 0.008% |

**Bath-grid convergence:** for lambda=0.05, c=5, N=8,500 versus 10,500, the direct BKM defect changes by relative 9.9872758e-15.

**Mechanism of the earlier apparent failure.** The earlier V0.7 direct calculations at lambda=0.05, c=5 reported approximately 6.33e-39 from a long-double propagation, compared with a corrected-pole prediction 3.97e-41. Raising arithmetic precision alone was not sufficient: at these late times, the discretized bath-frequency nodes and form factors must also be generated at >30-digit precision, or small errors in a rapidly decaying fast-channel matrix element dominate the BKM defect. The recomputation applies arbitrary precision to the bath inputs, compensated double-double arithmetic to a *direct Chebyshev propagation of the full Hermitian discretized Hamiltonian*, and 100-digit arithmetic to the entropy and BKM evaluation. Bessel coefficients and truncation are computed with at least 55-digit precision. The Chebyshev cutoff was selected so that the last coefficient is below approximately 1e-95.

**Interpretation.** The values at c=4.5, 5.0, 5.5, across and beyond the finite-c crossover, agree with the previously computed corrected-resonance predictions to better than 0.01%. The coefficient-weighted slow/fast ratio is not 1 exactly at the nominal c_BKM=4: the crossing of terms at finite coupling is displaced by logarithms of the prefactors. The data validate the **second threshold regime for one explicit order-two fixed positive-frequency bath**, not every m>=2 and not a uniform numerical theorem. The gap-1 example gamma_s=1, gamma_f=2, whose nominal c_BKM=8, was also tried at lambda=0.05, c=9 but did not have quantitatively convergent asymptotic amplitudes (direct/leading ≈4.24); that case is not certified here.

**Mathematical status.** Direct computational agreement is not proof of the analytic theorem. The compared corrected-pole values are from the existing P030 package and are not newly evaluated as independent high-precision contour integrals. A true independent continuum-limit certificate would require rigorous quadrature and Chebyshev truncation error bounds at the 1e-42 level, beyond empirical N convergence.

**Method reproduction.** `prep_hp.py` generates an N-mode positive semidefinite density and Hermitian-system-to-bath couplings using mpmath (55 decimal digits), with R=2. `cheb_dd.cpp` applies the discretized *Hermitian Hamiltonian directly*, advancing Chebyshev vectors by hardware-FMA-based compensated double-double arithmetic (about 31 digits) and summing high-precision Bessel coefficients. `evaluate_dd.py` reconstructs the 2x2 microscopic survival matrix and evaluates entropy/BKM quantities by independent 100-digit Hermitian diagonalizations. The supplied source and raw outputs reproduce the reported numbers. The original manuscript and V0.7 archive remain unchanged.
