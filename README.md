# P030 — BKM entropy and nonnormal quantum decay: code and data

Reproducibility materials for the manuscript *Memory-Induced Entropic Crossovers
in Nonnormal Quantum Decay: Uniform Reconstruction and BKM Defects*, including
its Supplementary Numerics.

Author: Hassan Nasreddine ([ORCID 0009-0001-0339-0978](https://orcid.org/0009-0001-0339-0978)).

## Numerical result

For the exactly orthogonal order-two bath with λ = 0.05, c = 5, N = 8500,
(γs,γf) = (1,3), T = 2000 ln(20), C0 = 0.2 I, b = (0,0.1) and a0 = 0.6,
the computer-assisted continuum enclosure is

```text
3.95299e-41 < Delta_BKM,continuum < 3.98513e-41
```

The certificate derives this one-point bound using an exact-moment polynomial,
a finite bath calculation, continuum discretization estimates and an explicit
arithmetic error budget. Its computational assumptions include IEEE 754
round-to-nearest-even, disabled floating-point contraction and the specified
exact/interval-arithmetic libraries. See `docs/PROOF_MAP.md`.

## Contents

- `REPRODUCIBILITY/code/` — simulation and plotting code for the paper and
  Supplementary Numerics.
- `REPRODUCIBILITY/data/` — numerical results underlying their tables and figures.
- `REPRODUCIBILITY/order2_high_precision/` — separately labeled numerical
  calculations for the order-two reservoir discussed in the supplement.
- `REPRODUCIBILITY/exact_m2_certificate/` — exact-moment model and inputs,
  recorded numerical output, interval and analytical certificates, C++ and
  Python reproduction programs.
- `docs/PROOF_MAP.md` — direct map from mathematical obligations to sources.
- `SHA256SUMS.txt` and `RELEASE_MANIFEST.json` — file-level integrity inventory.

## Reproduction

Test environment: Linux x86-64, Python 3.13, C++17 with libquadmath.

```bash
python -m pip install -r requirements-tested.txt
python tools/verify_release.py
bash RUN_PRECOMPUTED.sh
bash RUN_EXACT.sh
python tools/verify_release.py
```

`RUN_PRECOMPUTED.sh` verifies the supplied recorded values and error budget.
`RUN_EXACT.sh` additionally rebuilds the microscopic propagation, recomputes
interval enclosures and reruns the instrumented double-double calculation.
Read `REPRODUCIBILITY/exact_m2_certificate/README_REPRODUCE.md` for exact
compiler flags, mathematical assumptions and output interpretation.

## Citation and licensing

The software and author-owned data are licensed BSD-3-Clause (`LICENSE` and
`DATA_LICENSE.md`). The citation metadata are in `CITATION.cff`; cite the
version-specific archived DOI once that version is deposited on Zenodo.
