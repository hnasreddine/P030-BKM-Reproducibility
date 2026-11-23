# Continuum discretization and BKM transfer

The main derivations are in `CONTINUUM_BOUND_AND_TRANSFER.md`.
`python full_certified_operator_bound.py` computes the derivative enclosures,
Poisson alias bound and resulting 8500-mode continuum error. The phase of the
shifted grid is `(-1)^k` and cancels under absolute-value bounds.

The complete single-point BKM enclosure also requires the C++ propagation,
Bessel interval bounds and positive-kernel integration described in
`../CERTIFICATE_IDEAL_ORDER2_2026-11-22.md`.
`estimate_actual_Cr.py` is exploratory and does not form part of the bound.
