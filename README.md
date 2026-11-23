# P030 — BKM entropy and nonnormal quantum decay: reproduction code and data

Version **1.0.0**. Author: Hassan Nasreddine (ORCID
[0009-0001-0339-0978](https://orcid.org/0009-0001-0339-0978)).

This **standalone public code-and-data package** supports the numerical
experiments and computer-assisted mathematical enclosure in *Memory-Induced
Entropic Crossovers in Nonnormal Quantum Decay: Uniform Reconstruction and BKM
Defects*. It does not include the unpublished manuscript, private project
history, AI referee reports or conversational records.

## Certified numerical scope

For the **exactly orthogonal order-2 bath** at `lambda=0.05`, `c=5`,
`N=8500`, `(gamma_s,gamma_f)=(1,3)`, `T=2000 ln(20)`, and initial state
`C0=0.2 I`, `b=(0,0.1)`, `a0=0.6`, the supplied computer-assisted bounds give:

```text
3.95299e-41 < Delta_BKM,continuum < 3.98513e-41
```

This enclosure has explicitly stated mathematical and arithmetic hypotheses;
it is **not** machine-checked in a proof assistant. In particular the
specified C++ run uses IEEE 754 binary64 round-to-nearest-even with
`-fno-fast-math -ffp-contract=off`; the calculations also rely on the named
exact-rational and interval-arithmetic libraries. No comparable bound is
asserted for different parameters or historical table entries.

## Archive structure

- `REPRODUCIBILITY/code/` — simulation and plotting scripts supporting the
  manuscript's historical figures and tables.
- `REPRODUCIBILITY/data/` — numerical tables and parameters.
- `REPRODUCIBILITY/order2_high_precision/` — separate historical rounded-SVD
  direct-bath experiments (not the exact-polynomial certificate).
- `REPRODUCIBILITY/exact_m2_certificate/` — exactly orthogonal 8500-mode
  inputs, recorded outputs, analytic proof notes, double-double C++ code,
  numerical interval scripts and the continuum-discretization certificate.
- `tools/verify_release.py` — checks content, SHA-256 inventory, mandatory
  inputs, code invariants and that no internal referee report is present.
- `tools/reproduce_certified_point.sh` — recompiles and replays the exact bath.
- `SHA256SUMS.txt` and `RELEASE_MANIFEST.json` — checksums and archive inventory.

## Reproduce the one-point enclosure

```sh
python tools/verify_release.py
bash tools/reproduce_certified_point.sh
bash tools/reproduce_certified_point.sh full  # slower mathematical checks
```

See `REPRODUCIBILITY/exact_m2_certificate/README_REPRODUCE.md` for packages,
compiler flags, proof assumptions and numerical limitations.

## Citation and license

Please cite the archived **version DOI** assigned by Zenodo to the GitHub
release when it exists. No DOI is invented or pre-assigned in this package.
Code and author-owned data are distributed under `LICENSE` / `DATA_LICENSE.md`
(BSD-3-Clause). This repository and the archive are independent of the
unpublished manuscript.

## Release engineering and GitHub Actions

To mirror the verification discipline of the separate P014 archive, this repository
provides `RUN_PRECOMPUTED.sh`, `RUN_EXACT.sh`, `requirements-tested.txt`,
`docs/PROOF_MAP.md`, `.github/workflows/verify.yml`, a public/private release
boundary and `PUBLISHING_CHECKLIST.md`. Reproduction is conditional on the
mathematical proof obligations and compiler/environment assumptions stated above.

No GitHub public release or Zenodo DOI is claimed by this source archive.
