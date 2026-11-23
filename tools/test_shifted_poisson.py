#!/usr/bin/env python3
"""Independent Gaussian diagnostic for the shifted-grid Poisson identity.

This tests lattice phase/parity, not the certified continuum bound. The latter
is obtained from analytic derivative estimates and is unchanged by this fix.
"""
import cmath
import math


def check(n: int) -> tuple[float, float]:
    h = 2.0 / (n - 1)
    period = 2.0 * math.pi / h
    a, center, s = 100.0, -0.17, 1.3

    def transform(t: float) -> complex:
        return math.sqrt(math.pi / a) * math.exp(-t*t / (4*a)) * cmath.exp(-1j*center*t)

    lattice = sum(
        h * math.exp(-a*(x-center)**2) * cmath.exp(-1j*x*s)
        for j in range(-40, 41)
        for x in [-1.0+j*h]
    )
    aliases = sum(
        ((-1)**(k*(n-1))) * transform(s+k*period)
        for k in range(-15, 16)
    )
    wrong = sum(transform(s+k*period) for k in range(-15, 16))
    error = abs(lattice-aliases)
    if error >= 1e-12:
        raise AssertionError(f"N={n} shifted Poisson residual {error:.3g}")
    if n % 2 == 0 and abs(lattice-wrong) <= 1e-3:
        raise AssertionError("Unshifted aliasing spuriously passes for even N")
    return error, abs(lattice-wrong)


if __name__ == '__main__':
    for n in (6, 7):
        e, bad = check(n)
        print(f"PASS: N={n}, alias residual={e:.3e}, missing-phase discrepancy={bad:.3e}")
    print("PASS: shifted-grid parity regression; no continuum or BKM computation rerun")
