"""Certified Duhamel comparison of exact orthogonal p2 with its high-precision
recorded 60-digit decimal coefficients for P030 N=8500.
The mpmath.iv moment certificate bounds sup |delta p| < 5e-30.
Derivation: where theta>0, both J matrices have lambda_min >=.19*theta.
For the positive matrix square roots, Sylvester inequality
||sqrtJ - sqrtJ'||_F <= ||J-J'||_F /(sqrt(lambda_min J)+sqrt(lambda_min J')).
The Hamiltonian H has off-diagonal lambda * sqrt(h)*sqrt(J).
"""
import mpmath as m
m.iv.dps=70;I=m.iv.mpf
N=8500;epsilon=I('.05');T=I(2000)*m.iv.log(I(20));dx=I(2)/(N-1)
# Certified precision threshold from certify_polynomial_moments.py, rounded outwards
p_error=I('5e-30')
# 0 <= theta <= 1, support width .8, grid sum h*theta <= .8+h <.801
sumbound=I('.801')
# For J= diag(gamma phi/pi)+theta [[2,1.8p],[1.8p,2]], |p|<=1+5e-30,
# lambda_min J >=(2-1.8*(1+5e-30))*theta>.19theta.
# ΔJ has off-diag entries 1.8 theta Δp; Frobenius =sqrt(2)*1.8*theta*|Δp|.
sqrt_factor=I('1.8')*m.iv.sqrt(I(2))/(I(2)*m.iv.sqrt(I('.19')))
G_diff=sqrt_factor*m.iv.sqrt(sumbound)*p_error
# The *operator* norm of the Hermitian offdiagonal block ΔH=[0,λΔG;λΔG*,0] is λ||ΔG||.
H_diff=epsilon*G_diff
V_diff=T*H_diff
assert H_diff.b < I('7e-31').a
assert V_diff.b < I('4e-27').a
print('p_coeff_sum_less_than',p_error)
print('sqrt_matrix_lipschitz_factor',sqrt_factor.b)
print('||G_ideal - G_decimal||_F_upper',G_diff.b)
print('||H_ideal - H_decimal||_operator_upper',H_diff.b)
print('T',T)
print('||V_ideal - V_decimal||_operator_upper',V_diff.b)
print('certificate_roundup_less_than_4e-27: PASS')
