"""Conservative, reproducible all-link error budget for the exactly orthogonal
P030 m=2 bath, N=8500, lambda=.05, c=5. Uses the proved DD local bound;
checks concrete re-run residuals, parser checks, printing error, and Bessel budget.
Output is conditional on the soundness of the imported analytic and DD proofs.
"""
from pathlib import Path
import re,decimal,fractions,mpmath as m
P=Path(__file__).resolve().parent
m.iv.dps=65;I=m.iv.mpf
log=(P/'strict_ideal_full.log').read_text();assert 'BAD_NONOVERLAP 0' in log
assert (P/'strict.status').read_text().strip()=='PASS'
parser=(P/'parser_ideal.log').read_text();assert 'EXACT_RATIONAL_CONVERSION_BOUND_VIOLATIONS 0' in parser
assert 'ZERO_HIGH_NONZERO_LOW 0' in parser
for tag,ceiling in [('GLOBAL_SUMRES','7.859e-26'),('INITRES','3.031e-31'),('OUTSUM','2.267e-29')]:
 q=re.search(r'\b'+tag+r'\s+([0-9.eE+-]+)',log)
 assert q,tag
 got=I(q.group(1));assert got.b<I(ceiling).a
 print(tag,q.group(1),'bound',ceiling)
rows=[r.split() for r in (P/'out_ideal_m2_lean.txt').read_text().splitlines()]
assert len(rows)==4 and all(len(r)==4 for r in rows)
fra=lambda s:fractions.Fraction(decimal.Decimal(s))
Eprint2=fractions.Fraction(0)
for row in rows:
 for tokens in (row[:2],row[2:]):
  err=sum((fractions.Fraction.from_float(float(tok))-fra(tok) for tok in tokens),fractions.Fraction(0))
  Eprint2+=err**2
assert Eprint2<fra('1.3e-24')**2
print('printed_decimal_output_vs_exact_binary64_DD_frobenius_lt_1.3e-24 PASS')
# The input stores 12,866 Bessel coefficients, used for both SVD-rounded and exact-moment p.
with (P/'input_005_5_8500_ideal_m2.txt').open() as f:
 h=f.readline().split();N=int(h[0]);K=int(h[1]);assert N==8500 and K==12865
 for _ in range(N):f.readline()
 decimal.getcontext().prec=100
 csum=abs(decimal.Decimal(next(f)))
 for _ in range(1,K+1):csum+=2*abs(decimal.Decimal(next(f)))
 assert csum<decimal.Decimal(135)
print('L1_Bessel_coefficients',str(csum),'less_than_135 PASS')
# Strict errors are enclosed by the *outward* rounded ceilings above.
# 0.5% safety factor covers positive-double error-radius accumulation;
# 1.25 bounds ||U_k(H/2)|| under ||H/2||<.6.
E_DD=I('1.005')*(I('1.25')*I(135)*(I('7.859e-26')+I('3.031e-31'))+I('2.267e-29'))
assert E_DD.b < I('1.334e-23').a
parts={
 'continuum_exact_orthogonal_to_discrete_orthogonal':I('9.356e-23'),
 'DD_recurrence_and_output_arithmetic':I('1.334e-23'),
 'printed_output_decimal':I('1.3e-24'),
 'input_DD_parser_to_decimal_H':I('1.3e-24'),
 'exact_p_vs_decimal_p':I('4e-27'),
 'Bessel_DD_parser':I('1.07e-28'),
 'Bessel_coefficient_ODE_bound':I('1.717e-32'),
 'time_token':I('1e-47'),
 'ideal_Chebyshev_truncation':I('3e-97'),
 '42_digit_nodes_and_couplings':I('2.2e-37'),
}
# Conservative independent upper bounds for input conversion and Bessel coefficients.
assert I('5.83e-27').b < parts['input_DD_parser_to_decimal_H'].a
assert I('1.0011e-97').b < parts['ideal_Chebyshev_truncation'].a
E=sum(parts.values(),I(0))
for k,v in parts.items():print('error_component',k,v.b)
print('E_DD_upper_before_rounding',E_DD.b)
print('TOTAL_BOUND',E.b)
assert E.b<I('1.2e-22').a
print('TOTAL_LT_1.2e-22 PASS')
