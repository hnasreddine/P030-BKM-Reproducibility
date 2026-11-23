"""Reproduce the ten original supplementary tables from recorded data.

Includes two separately derived order-two tables published later. All twelve
tables are compared to the final JMP supplementary LaTeX numerical entries.
No prose or historical narrative output is distributed by this generator.
"""
import json, math
from pathlib import Path
import argparse, ast, re
from decimal import Decimal
D = str((Path(__file__).resolve().parent.parent/'data').resolve()) + '/'
A = json.load(open(D+'testA_m1.json')); B1 = json.load(open(D+'testB_v03_m2.json')); B2 = json.load(open(D+'testB_poly_m2.json'))
C = json.load(open(D+'testC_EP.json')); Dc = json.load(open(D+'testD_clusters.json')); Fx = json.load(open(D+'testF_fixed_contour.json'))
Hh = json.load(open(D+'testH_pole_asymptotics.json')); E = json.load(open(D+'testE_convergence.json')); S = json.load(open(D+'testS_static.json'))
def sci(x, d=1):
    if x == 0: return '0'
    e = math.floor(math.log10(abs(x))); m = x/10**e
    if round(m, d) >= 10: m /= 10; e += 1
    return r'%.*f\times10^{%d}' % (d, m, e)
def get(rows, **kw):
    r = [x for x in rows if all(abs(x[k]-v) < 1e-9 if isinstance(v, float) else x[k] == v for k, v in kw.items())]
    assert len(r) == 1, (kw, len(r)); return r[0]
out = []
w = out.append
w(r'''\section{Numerical tests}\label{sec:numerics}
The tests below compare three objects: (i) the exact propagator $V_\eps(T)$ of the microscopic one-excitation Hamiltonian, (ii) the reconstruction $Z=e^{TA_\eps}S_\eps$, and (iii) the asymptotic formulas of Theorem~\ref{thm:thresholds}. They are numerical evidence, not proofs; in particular they are carried out at finite coupling, where the truncation length is chosen by hand.

\subsection{Methods}
\emph{Direct propagation.} The continuum $x\in[-1,1]$ is replaced by $N$ equally spaced modes $x_n$ with couplings $\sqrt{\Delta x}\,J(x_n)^{1/2}$. For densities that vanish to infinite order at the band edges this trapezoidal discretization reproduces the continuum propagator up to an aliasing error governed by the correlation function at the recurrence time $2\pi/\Delta x$; we take $2\pi/\Delta x\ge T+45/\eps$ and verify convergence in $N$. The $(d+dN)$-dimensional Hermitian matrix $H_1$ is applied through its arrowhead structure, and $V_\eps(T)=\sum_k(2-\delta_{k0})(-i)^kJ_k(RT)\,P_{\mathcal E}T_k(H_1/R)P_{\mathcal E}$ is summed with Chebyshev moments computed once in 80-bit floating point and Bessel coefficients computed by downward recurrence in multiprecision. No time stepping is involved. Grid and couplings are generated in 80-bit arithmetic; with a double-precision grid the result carries an $N$-dependent error of $10^{-14}$.

\emph{Reconstruction.} $F(s)$ is computed by the trapezoidal rule in $x$ and $K(z)$ by composite Gauss--Legendre quadrature in $s$, in 80-bit arithmetic, with a $C^\infty$ cutoff $\chi(s/L)$. The contour is a circle of radius $O(\eps)$ around the resonances (or a thin rectangle for well separated groups), deformable to $\Gamma_\eps$ without crossing zeros; $S,P$ are trapezoidal sums, and $A=PS^{-1}$, $e^{TA}S$ are evaluated in multiprecision. At finite $\eps$ the cutoff $L$ must satisfy two constraints: the tail $\eps T\int_L^\infty\norm F$ must be negligible, and the contour must stay within a few multiples of $1/L$ of the imaginary axis. We use $L=\min(300,2/\eps)$ unless stated otherwise.

\emph{Entropies.} Given a $2\times2$ matrix $V$, $\cD$ and $\cB$ are evaluated from eigenvalues with 80--120 significant digits; the defect is their difference.

\subsection{Reconstruction accuracy, exceptional point, several frequency groups}
Table~\ref{tab:recon} shows $\norm{V_\eps(T)-Z}$ for the bath of Example~\ref{ex:m1} with $\gamma_s=1$, $\gamma_f=2$, $\xi=0.4$, together with the error of the Davies semigroup. The reconstruction error decreases faster than any fixed power of $\lambda$ and reaches the accuracy of the arithmetic ($\sim10^{-18}$) at $\lambda\le0.07$. Where it is above that floor it is independent of $L$ once $L\gtrsim150$ and of the contour, and is therefore the genuine non-exponential part of the microscopic dynamics, of the order of the reservoir correlation function at time $T$.
''')
lams = [0.2, 0.14, 0.1, 0.07, 0.05, 0.035]
w(r'''\begin{table}[h]
\caption{Bath of Example~\ref{ex:m1}. Operator-norm errors of the Davies semigroup $e^{-\tau G_0}$ and of the reconstruction $Z$ against direct propagation, at $\tau=c\log(1/\lambda)$.}\label{tab:recon}
\begin{ruledtabular}
\begin{tabular}{ccccccc}
 & \multicolumn{2}{c}{$c=1$} & \multicolumn{2}{c}{$c=2$} & \multicolumn{2}{c}{$c=4$}\\
$\lambda$ & Davies & $Z$ & Davies & $Z$ & Davies & $Z$\\ \hline''')
for lam in lams:
    r = [get(A, lam=lam, c=c) for c in (1.0, 2.0, 4.0)]
    w('$%.3g$ & ' % lam + ' & '.join('$%s$ & $%s$' % (sci(x['err_davies']), sci(x['err_pole'])) for x in r) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
# EP
w(r'''\emph{Exactly at an exceptional point.} For $J(x)=\phi(x)\bigl(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\bigr)/2\pi$ and $D=\eps\,\mathrm{diag}(\delta/2,-\delta/2)$ the Davies matrix has an exceptional point at $\delta=1$. For each $\eps$ we determine by Newton iteration the detuning $\delta_*(\eps)$ at which the two zeros of the truncated pencil coalesce on the real axis; the eigenvalues of $A_\eps/\eps$ then differ by about $10^{-8}$ and the condition number of its eigenvector matrix is about $10^8$. Table~\ref{tab:ep} shows that the reconstruction error at $\delta_*$ is the same as at $0.9\,\delta_*$ and $1.1\,\delta_*$, and that the entropy and the defect are reproduced with small relative errors there as well.
\begin{table}[h]
\caption{Exceptional-point test at $\tau=2.2\log(1/\lambda)$: coalescence detuning, eigenvalue splitting of $A_\eps/\eps$ at $\delta_*$, reconstruction errors below, at and above the exceptional point, and relative errors of entropy and defect at $\delta_*$ ($C_0=0.2\,\Id$, $b=(0,0.1)$).}\label{tab:ep}
\begin{ruledtabular}
\begin{tabular}{cccccccc}
$\lambda$ & $\delta_*$ & splitting & $0.9\,\delta_*$ & $\delta_*$ & $1.1\,\delta_*$ & rel.\ err.\ $\cD$ & rel.\ err.\ $\cD-\cB$\\ \hline''')
for lam in [0.2, 0.14, 0.1, 0.07, 0.05]:
    b_, e_, a_ = [get(C, lam=lam, c=2.2, tag=t) for t in ('below', 'EP', 'above')]
    w('$%.3g$ & $%.5f$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$' % (lam, e_['delta_star'], sci(e_['eig_split']), sci(b_['err_pole']), sci(e_['err_pole']), sci(a_['err_pole']), sci(e_['relD']), sci(e_['relF'])) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
# clusters
gaps = ['0.60', '0.10', '0.02', '3eps', '0.3eps', '0']
w(r'''\emph{Several Bohr-frequency groups.} For $d=4$ we take $D=\mathrm{diag}(-\tfrac{\mathsf g}2,-\tfrac{\mathsf g}2+0.7\eps,\tfrac{\mathsf g}2,\tfrac{\mathsf g}2-1.3\eps)$, two doublets with internal splittings of order $\eps$ separated by a gap $\mathsf g$, coupled through one bath $J(x)=0.1\,\phi(x)B(x)B(x)^{\mathsf T}/\pi$ with a fixed real $4\times4$ matrix polynomial $B(x)=B_0+xB_1$ (the Davies rates range from $0.004$ to $1.93$). A single thin rectangular contour encloses all four resonances. Table~\ref{tab:clusters} shows that the error does not deteriorate as $\mathsf g$ goes from $0.6$ to $0$; in all cases $\norm{S_\eps-\Id}/\eps\in[2.6,2.8]$ and $\norm{A_\eps+iD}/\eps\in[1.8,2.0]$.
\begin{table}[h]
\caption{Four levels in two groups separated by a gap $\mathsf g$: $\norm{V_\eps(T)-Z}$ at $\tau=2.2\log(1/\lambda)$ (direct propagation in double precision, floor $\sim10^{-14}$; $L=\min(260,1.5/\text{left edge})$).}\label{tab:clusters}
\begin{ruledtabular}
\begin{tabular}{ccccccc}
$\lambda$ & $\mathsf g=0.6$ & $0.1$ & $0.02$ & $3\eps$ & $0.3\eps$ & $0$\\ \hline''')
for lam in [0.1, 0.07, 0.05, 0.035]:
    w('$%.3g$ & ' % lam + ' & '.join('$%s$' % sci(get(Dc, lam=lam, gap=g)['err_pole']) for g in gaps) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
# fixed contour
w(r'''\emph{Fixed contour.} Table~\ref{tab:fixedcontour} documents Remark~\ref{rem:fixedcontour}.
\begin{table}[h]
\caption{Truncated pencil of Example~\ref{ex:m1} with $L=\eps^{-\nu}$: number of zeros of $\det\fT$ (argument principle) and deviation of the zeroth moment from the identity, for a rectangle with left edge at $\re z=-1/L$ and for the fixed rectangle $[-0.3,0.5]\times[-1.3,1.3]$; $\max\eps\norm{K}$ is taken over the left edge of the fixed rectangle.}\label{tab:fixedcontour}
\begin{ruledtabular}
\begin{tabular}{cccccccc}
$\lambda$ & $\nu$ & $L$ & zeros (thin) & $\norm{S-\Id}/\eps$ (thin) & zeros (fixed) & $\norm{S-\Id}$ (fixed) & $\max\eps\norm K$ (fixed)\\ \hline''')
for x in Fx:
    w('$%.3g$ & $%.2f$ & $%.0f$ & $%d$ & $%.2f$ & $%d$ & $%s$ & $%s$' % (x['lam'], x['alpha'], x['L'], round(x['zeros_in_thin']), x['S_thin_minus_I_over_eps'], round(x['zeros_in_fixed']), sci(x['S_fixed_minus_I']), sci(x['eps_K_fixed'])) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
# thresholds m=1
cs1 = [1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0]
w(r'''\subsection{Transfer order $m=1$: both thresholds}
For the bath of Example~\ref{ex:m1}, $\delta\gamma=1$, $\kappa_1=-i\pi\xi a_\phi$, $\abs{h_1}=%.5f$, so $c_{\cD}=2$ and $c_{\mathrm{BKM}}=4$. Table~\ref{tab:m1} compares the microscopic entropy and defect, computed from the directly propagated $V_\eps$, with the two-term formulas \eqref{eq:thD}--\eqref{eq:thF} (with $\mathsf e_j$ taken from the computed resonances), for $C_0=0.2\,\Id$, $b=(0,0.1)$. The ratios tend to $1$ with a deviation proportional to $\eps$ (about $5\,\eps$ where the fast term dominates and $10\,\eps$ to $18\,\eps$ where the slow term does), as predicted by Remark~\ref{rem:hyp}(3); the same holds for a non-diagonal $C_0$ (data in the supplementary material). The last columns show the failure of the Davies semigroup in the relative sense and the relative accuracy of the reconstruction. Because of the finite-$\lambda$ shift of Remark~\ref{rem:hyp}(4), the slow term overtakes the fast one in the defect near $c\approx5$ at these couplings; at $c=6$ it carries more than $98\%%$ of the formula, and the microscopic defect follows it. Figure~\ref{fig:thresholds}(a,b) displays the data for $\lambda=0.035$.
''' % abs(-math.pi*0.4*0.15811363626379824))
w(r'''\begin{table}[h]
\caption{$m=1$. Ratios of microscopic entropy and defect to the two-term formulas \eqref{eq:thD}, \eqref{eq:thF}; share of the slow term in the formula for the defect; ratio of the microscopic entropy to that of the Davies semigroup; relative deviation of the reconstruction from the microscopic values.}\label{tab:m1}
\begin{ruledtabular}
\begin{tabular}{cccccccc}
$\lambda$ & $c$ & $\cD/\cD_{\mathrm{form}}$ & $(\cD-\cB)/(\cD-\cB)_{\mathrm{form}}$ & slow share & $\cD/\cD_{\mathrm{Davies}}$ & rel.\ err.\ $\cD(Z)$ & rel.\ err.\ $[\cD-\cB](Z)$\\ \hline''')
for lam in [0.1, 0.07, 0.05, 0.035]:
    for c in cs1:
        x = get(A, lam=lam, c=c); s = x['iso']
        if math.exp(-2.0*x['tau']) < 5e-17:
            w('$%.3g$ & $%.1f$ & $%.4f$ & --- & $%.3f$ & $%s$ & $%s$ & ---' % (lam, c, s['D_micro']/s['D_ren'], s['F_slow']/(s['F_slow']+s['F_fast']), sci(s['D_micro']/s['D_davies']), sci(s['relD_pole'])) + r'\\')
        else:
            w('$%.3g$ & $%.1f$ & $%.4f$ & $%.4f$ & $%.3f$ & $%s$ & $%s$ & $%s$' % (lam, c, s['D_micro']/s['D_ren'], s['F_micro']/s['F_ren'], s['F_slow']/(s['F_slow']+s['F_fast']), sci(s['D_micro']/s['D_davies']), sci(s['relD_pole']), sci(s['relF_pole'])) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
Entries marked --- are those for which the fast amplitude $e^{-\gamma_f\tau}$ is below $5\times10^{-17}$ and is not resolved by 80-bit arithmetic.
''')
w(r'''\begin{figure}[h]
\includegraphics[width=\textwidth]{fig_thresholds.pdf}
\caption{Entropy and BKM defect at $\tau=c\log(1/\lambda)$, $\lambda=0.035$, from direct propagation of the microscopic Hamiltonian (circles), with the fast and slow terms of \eqref{eq:thDlead}, \eqref{eq:thFlead} and the Davies values (dotted, coinciding with the fast term on this scale). (a,b) $m=1$, Example~\ref{ex:m1}; (c,d) $m=2$, bath \eqref{eq:polybath}. Vertical lines mark the asymptotic horizons $c_{\cD}$ and $c_{\mathrm{BKM}}$; at finite $\lambda$ the crossover is displaced as described in Remark~\ref{rem:hyp}(4). In (b) and (d) only points with $e^{-\gamma_f\tau}>5\times10^{-17}$ are shown.}\label{fig:thresholds}
\end{figure}
''')
# m=2
w(r'''\subsection{Transfer order $m=2$}
\emph{A reservoir for which the asymptotic regime is remote.} Take $\gamma_s=1$, $\gamma_f=2$ and an off-diagonal density made of three bumps of half-width $0.075$ with amplitude $0.025$, for which $\kappa_2=%s\,i$. At $\lambda\in[0.1,0.2]$ and $\tau=2\log(1/\lambda)$ the directly computed $\abs{V_{sf}}$ exceeds the leading pole term by factors between $6$ and $139$ (Table~\ref{tab:v03}); we report this failure of the leading-order formula at moderate coupling and extend the test to smaller couplings. The discrepancy is the non-exponential memory tail: to first order in $\eps$ the non-pole part of $V_{sf}(T)$ is $-\eps\int_T^\infty(s-T)F_{sf}(s)\,\dd s$, and the narrow bumps give a correlation function that decays only like $\exp(-\sqrt{0.15\,s})$. The residual $V_{sf}-V_{sf}^{\mathrm{pole}}$ agrees with this tail integral where the residual is large, and the ratio to the pole term converges to $1$ for $\lambda\le0.05$. The asymptotic statement is therefore not contradicted, but for this bath the asymptotic regime begins only below $\lambda\approx0.05$.
\begin{table}[h]
\caption{The three-bump $m=2$ bath at $\tau=2\log(1/\lambda)$: directly computed $\abs{V_{sf}}$, leading pole term $\eps^2\abs{\kappa_2}\abs{\gamma_f^2e^{-\gamma_f\tau}-\gamma_s^2e^{-\gamma_s\tau}}/\delta\gamma$, their ratio, first-Born tail integral, and ratio of the residual to the tail.}\label{tab:v03}
\begin{ruledtabular}
\begin{tabular}{ccccccc}
$\lambda$ & $T$ & $\abs{V_{sf}}$ & pole term & ratio & Born tail & residual/tail\\ \hline''' % ('%.5f' % 0.0075736587900006565))
for x in B1:
    w('$%.3g$ & $%.0f$ & $%s$ & $%s$ & $%.3f$ & $%s$ & $%.2f$' % (x['lam'], x['T'], sci(x['Vsf_abs'], 2), sci(x['pole_abs'], 2), x['ratio_pole'], sci(x['born_tail_abs'], 2), x['resid_over_born']) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
x0 = get(B2, lam=0.035, c=1.0)
w(r'''\emph{A bath with an order-one transfer coefficient.} We use \eqref{eq:polybath} with $\gamma_s=1$, $\gamma_f=3$, $A=2$, $\eta=1.8$, $\vartheta(x)=\phi((x-0.5)/0.4)$ and $m=2$, for which $\kappa_2=0.31845\,i$, $\abs{h_2}=0.15923$, $c_{\cD}=2$, $c_{\mathrm{BKM}}=4$. Table~\ref{tab:m2} shows the comparison with \eqref{eq:thD}--\eqref{eq:thF}. The first threshold is crossed: for $c\ge2.5$ the entropy is carried by the slow term (at $\lambda=0.035$, $c=3$, it exceeds the Davies value by a factor $%s$) and agrees with the formula up to a deviation that decreases proportionally to $\eps$, from $26\%%$ at $\lambda=0.07$ to $7\%%$ at $\lambda=0.035$. The constant is large (about $57\,\eps$) because the off-resonant weight $A\vartheta$ makes the first diagonal jet large. The defect agrees with its fast term as long as $\abs{V_{ff}}$ is resolved. The second threshold could not be reached by direct propagation: at $\lambda\le0.05$ the fast amplitude at the crossover is below $10^{-17}$, and at $\lambda\ge0.07$ the memory tail and the $O(\eps)$ corrections are not small. For $m=2$ the second threshold is thus supported only by the test of the next subsection.
''' % sci(get(B2, lam=0.035, c=3.0)['D_micro']/get(B2, lam=0.035, c=3.0)['D_davies']))
w(r'''\begin{table}[h]
\caption{$m=2$, bath \eqref{eq:polybath}. Columns as in Table~\ref{tab:m1}, with the share of the slow term in the entropy formula.}\label{tab:m2}
\begin{ruledtabular}
\begin{tabular}{ccccccc}
$\lambda$ & $c$ & $\abs{V_{ff}}$ & slow share in $\cD_{\mathrm{form}}$ & $\cD/\cD_{\mathrm{form}}$ & $(\cD-\cB)/(\cD-\cB)_{\mathrm{form}}$ & $\cD/\cD_{\mathrm{Davies}}$\\ \hline''')
for lam in [0.07, 0.05, 0.035]:
    for c in [1.0, 2.0, 2.5, 3.0, 3.5, 4.0]:
        x = get(B2, lam=lam, c=c)
        f = ('$%.4f$' % (x['F_micro']/x['F_ren'])) if math.exp(-3.0*x['tau']) > 5e-17 else '---'
        w('$%.3g$ & $%.1f$ & $%s$ & $%.3f$ & $%.4f$ & %s & $%s$' % (lam, c, sci(x['V_ff_abs']), x['D_slow']/(x['D_slow']+x['D_fast']), x['D_micro']/x['D_ren'], f, sci(x['D_micro']/x['D_davies'])) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
# H
w(r'''\subsection{The two-pole dynamics deep in the asymptotic regime}
Theorem~\ref{thm:thresholds} combines the reconstruction (tested above against the microscopic dynamics) with an asymptotic analysis of the two-pole expression \eqref{eq:twopole}. The second ingredient can be tested at much smaller coupling and across both thresholds without floating-point limitations: we compute $\zeta_j$ and $Q_j$ from the truncated pencil in 80-bit arithmetic ($L=\min(1500,3/\eps)$) and evaluate \eqref{eq:twopole} and the entropies with 120 digits; the entries of $Q_j$ of order $\eps^{2m}$, which are below 80-bit resolution, are fixed by the exact rank-one structure of Proposition~\ref{prop:poles}. Table~\ref{tab:H} gives, for the baths \eqref{eq:polybath} with $m=1,2,3$ ($\gamma_s=1$, $\gamma_f=3$, $A=2$, $\eta=1.8$) and a non-diagonal $C_0$, the residue ratio $(Q_s)_{sf}/h_m\eps^m$ of Proposition~\ref{prop:poles} and the ratios of entropy and defect to \eqref{eq:thD}--\eqref{eq:thF} at five horizons. All ratios converge to $1$ linearly in $\eps$. This test does not involve the microscopic Hamiltonian.
\begin{table}[h]
\caption{Two-pole model \eqref{eq:twopole} for $m=1,2,3$: $\abs{(Q_s)_{sf}/h_m\eps^m}$ and ratios (exact two-pole value)/(formula) for $\cD$ and $\cD-\cB$ at $c/c_{\cD}=0.5,1$, at the midpoint $c=\frac12(c_{\cD}+c_{\mathrm{BKM}})$, and at $c/c_{\mathrm{BKM}}=1,1.5$.}\label{tab:H}
\begin{ruledtabular}
\begin{tabular}{cccccccccccc}
 & & & \multicolumn{4}{c}{$\cD/\cD_{\mathrm{form}}$} & \multicolumn{5}{c}{$(\cD-\cB)/(\cD-\cB)_{\mathrm{form}}$}\\
$m$ & $\eps$ & residue & $0.5c_{\cD}$ & $c_{\cD}$ & mid & $1.5c_{\mathrm{BKM}}$ & $0.5c_{\cD}$ & $c_{\cD}$ & mid & $c_{\mathrm{BKM}}$ & $1.5c_{\mathrm{BKM}}$\\ \hline''')
for m in (1, 2, 3):
    for eps in (1e-2, 1e-3, 1e-4):
        rr = [x for x in Hh['residues'] if x['m'] == m and abs(x['eps']-eps) < 1e-12]
        if not rr: continue
        res = math.hypot(*rr[0]['Qsf_over_h_eps_m'])
        rows = [x for x in Hh['asymptotics'] if x['m'] == m and abs(x['eps']-eps) < 1e-12]
        def pick(f): return min(rows, key=lambda x: abs(f(x)))
        r1 = pick(lambda x: x['c_over_cD']-0.5); r2 = pick(lambda x: x['c_over_cD']-1.0); r3 = pick(lambda x: x['c_over_cD']-1.5); r4 = pick(lambda x: x['c_over_cB']-1.0); r5 = pick(lambda x: x['c_over_cB']-1.5)
        w('$%d$ & $10^{%d}$ & $%.4f$ & ' % (m, round(math.log10(eps)), res) + ' & '.join('$%.4f$' % x['D_over_Dren'] for x in (r1, r2, r3, r5)) + ' & ' + ' & '.join('$%.4f$' % x['F_over_Fren'] for x in (r1, r2, r3, r4, r5)) + r'\\')
w(r'''\end{tabular}
\end{ruledtabular}
\end{table}
''')
# ---- nonnormal Davies generator (v0.5) ----
N1 = json.load(open(D+'testN1_effective.json')); N2 = json.load(open(D+'testN2_micro.json')); N3 = json.load(open(D+'testN3_generic.json')); N4 = json.load(open(D+'testN4_convergence.json'))
me = N2['meta']; kp = me['consts']['pI']; kg = me['consts']['general']
def resolved(r):
    p = r['pI']; return (r['err_pole']/r['fast_amp'])**2*p['F_fast']/p['F_slow'] < 1e-3
resrat = {lam: [x for x in N2['rows'] if abs(x['lam']-lam) < 1e-9][0]['Qs_w_over_eps_alpha'] for lam in (0.2, 0.1, 0.07, 0.05, 0.035)}
dav = [x for x in N1['a_davies_semigroup'] if x['C0'] == 'pI']
nd = N1['b_orthogonal_but_not_dark']
w(r"""\subsection{Nonnormal Davies generator: direct microscopic test of Theorem~\ref{thm:thresholds}}\label{sec:num-nonnormal}
We take Example~\ref{ex:nonnormal} with $q=\tfrac12$, $b_0=0.1$ and $\phi(x)=\exp(1-(1-x^2)^{-1})$. The quantities entering the theorem are computed independently of the propagation: $\mathsf b_1=%.6f$ by quadrature, $\alpha=%.6f$ from \eqref{eq:alpha1}, which agrees with \eqref{eq:alphaex} to machine precision, $w_\parallel=-0.05$, $w_\perp^2=0.01$, and for $C_0=0.2\,\Id$: $c_s=%.4g$, $c_f=%.4g$, $\varrho_0=%.4g$, $d_0=%.4g$. The spectral projector $P_s$ has norm $\sqrt{1+q^2}>1$.

\emph{Residue.} The slow residue $Q_s(\eps)$ is obtained from the contour moments without using \eqref{eq:Q1}. The ratio $\norm{Q_s(\eps)w}/(\eps\abs\alpha)$ equals $%.4f$, $%.4f$, $%.4f$, $%.4f$, $%.4f$ at $\lambda=0.2$, $0.1$, $0.07$, $0.05$, $0.035$: it tends to one with a deviation proportional to $\eps$.

\emph{Entropy and defect.} Table~\ref{tab:nonnormal} and Fig.~\ref{fig:nonnormal} compare the entropy and the BKM defect of the directly propagated microscopic state with \eqref{eq:thD} and \eqref{eq:thF} (with the numerically computed resonances $\zeta_s,\zeta_f$ and no fitted parameter). Both horizons $c_{\cD}=2$ and $c_{\mathrm{BKM}}=4$ are crossed: the slow term accounts for a fraction between $10^{-8}$ and more than $0.999$ of the defect formula. The ratios approach one from above with deviations proportional to $\eps$ (about $2.5\,\eps$ for the entropy and at most $7\,\eps$ for the defect beyond $c=2$). The Davies semigroup, for which the slow channel stays dark, underestimates the entropy by up to ten orders of magnitude. The column $\cD/\cD_{\mathrm{nt}}$ shows the formula of the normal case transplanted to the eigenbasis of $M_0$, $\norm w^2\mathsf e_f\log(1/\mathsf e_f)+\abs\alpha^2\eps^2\mathsf e_s\log(1/\mathsf e_s)$: it is off by $4$--$30\%%$ in the table, and below $c_{\cD}$ its ratio to the true entropy tends to $(2q^2+2)/(q^2+2)\neq1$ (see the Davies-level test below). A second initial state with a non-scalar excited block, $C_0=\bigl(\begin{smallmatrix}0.25&0.05\\0.05&0.15\end{smallmatrix}\bigr)$ ($c_s=%.4g$, $c_f=%.4g$, $\varrho_0=%.4g$), gives the same ratios within $7\times10^{-3}$ at all resolved points.

Two entries of the table are \emph{not} resolved and are marked: at $\lambda=0.035$, $c\ge6$ the fast amplitude $e^{-\gamma_f\tau}\approx3\times10^{-18}$ is below the absolute accuracy ($\sim10^{-16}$) of the 80-bit propagation, the smallest eigenvalue of the excited block is then wrong, and the computed defect is contaminated at the $10\%%$ level although the entropy is not. We keep these entries as a record of the numerical limit; no claim is based on them.
\begin{table}[t]
\caption{Nonnormal generator of Example~\ref{ex:nonnormal}, $q=\frac12$, $b_0=0.1$, $C_0=0.2\,\Id$, direct microscopic propagation at $\tau=c\log(1/\lambda)$. Ratios $\cD/\cD_{\mathrm{form}}$ and $(\cD-\cB)/(\cD-\cB)_{\mathrm{form}}$ to the formulas \eqref{eq:thD}, \eqref{eq:thF}; share of the slow term in the defect formula; ratio to the transplanted normal-case formula; ratio to the Davies semigroup; distance between direct propagation and reconstruction. $\dagger$: not resolved (see text).}\label{tab:nonnormal}
\begin{ruledtabular}
\begin{tabular}{llllllll}
$\lambda$ & $c$ & $\cD$ ratio & defect ratio & slow share & $\cD/\cD_{\mathrm{nt}}$ & $\cD/\cD_{\mathrm{Davies}}$ & $\norm{V_\eps-Z}$\\ \hline""" % (me['b1'], me['alpha'][0], kp['cs'], kp['cf'], kp['rho0'], kp['d0'], resrat[0.2], resrat[0.1], resrat[0.07], resrat[0.05], resrat[0.035], kg['cs'], kg['cf'], kg['rho0']))
maxdiff = 0.0
for lam in (0.1, 0.07, 0.05, 0.035):
    for c in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        r = get(N2['rows'], lam=lam, c=c); p = r['pI']; g = r['general']; ok = resolved(r)
        if ok: maxdiff = max(maxdiff, abs(p['D_micro']/p['D_form']-g['D_micro']/g['D_form']), abs(p['F_micro']/p['F_form']-g['F_micro']/g['F_form']))
        sh = p['F_slow']/(p['F_slow']+p['F_fast'])
        w('$%.3g$ & $%.0f$ & $%.4f$ & $%.4f$%s & $%s$ & $%.3f$ & $%s$ & $%s$' % (lam, c, p['D_micro']/p['D_form'], p['F_micro']/p['F_form'], '' if ok else r'$^\dagger$', ('%.4f' % sh) if sh > 1e-3 else sci(sh), p['D_micro']/p['D_naive'], ('%.3f' % (p['D_micro']/p['D_davies'])) if p['D_micro']/p['D_davies'] < 100 else sci(p['D_micro']/p['D_davies']), sci(r['err_pole'])) + r'\\')
assert maxdiff < 7e-3, maxdiff
w(r"""\end{tabular}
\end{ruledtabular}
\end{table}
\begin{figure}[t]
\includegraphics[width=\textwidth]{fig_nonnormal}
\caption{Nonnormal generator of Example~\ref{ex:nonnormal}, $\lambda=0.035$: entropy (left) and BKM defect (right) from direct microscopic propagation (circles; unresolved points omitted), the two terms of \eqref{eq:thD} and \eqref{eq:thF}, their sum, and the Davies semigroup. Vertical lines: $c_{\cD}=2$ and $c_{\mathrm{BKM}}=4$.}\label{fig:nonnormal}
\end{figure}

\emph{Davies level.} For the semigroup $e^{\tau M_0}$ alone (no memory, $\alpha=0$) formula \eqref{eq:thD} reduces to its first bracket plus $\abs{w_\parallel}^2\mathsf e_f[\log(1/\mathsf e_s)-\log c_s]$. At $\tau=5,10,20,40$ the ratio of the exact entropy to it is $%.5f$, $%.7f$, $%.9f$, $%.9f$, while the ratio to the transplanted normal-case expression is $%.4f$, $%.4f$, $%.4f$, $%.4f$, converging to $(q^2+2)/(2q^2+2)=0.9$ and not to one. For the vector $w=b_0(0,1)^{\mathsf T}$, which is orthogonal to $\hat u$ but has $l_sw=%.2f\neq0$, $\tau^{-1}\log\cD=%.3f$, $%.3f$, $%.3f$ at $\tau=5,10,20$, tending to $-2\gamma_s=-2$ and not to $-2\gamma_f=-4$: orthogonality to the slow eigenvector does not make the slow channel dark (Remark~\ref{rem:hyp}).
""" % tuple([x['D_over_formula'] for x in dav] + [x['D_over_naive_diagonal'] for x in dav] + [nd[0]['l_s_w']] + [x['logD_over_tau'] for x in nd]))
# generic complex nonnormal baths at pole level
w(r"""\emph{Generic nonnormal reservoirs with complex resonances.} Three reservoirs $J(x)=\phi(x)(B_0+xB_1)$ with randomly drawn complex Hermitian $B_0>0$, $B_1$ and Hermitian $\widetilde D$ (screened only for positivity of $J$ and $\delta\gamma\ge0.4$) have Davies matrices with complex, non-orthogonal eigenvectors (overlaps $\abs{\langle\hat u,\hat v_f\rangle}=%.2f$, $%.2f$, $%.2f$). Here $\omega_s\neq\omega_f$, the transfer coefficient \eqref{eq:alpha1} is complex, and the two-pole dynamics \eqref{eq:twopole} can be evaluated far below the couplings accessible to direct propagation. Table~\ref{tab:generic} shows the convergence of the residue and of the two formulas, and the failure of the formula in which the interference term of \eqref{eq:thD} is dropped.
\begin{table}[t]
\caption{Generic nonnormal reservoirs, two-pole dynamics \eqref{eq:twopole} with numerically computed residues, general $C_0$. $\abs{r_\eps w/(\alpha\eps)}$; ratios (exact two-pole value)/(formula) at $c/c_{\cD}=0.5,1,2,3$ ($c=2c_{\cD}$ is $c_{\mathrm{BKM}}$); last column: entropy divided by the formula without the interference term at $c=c_{\cD}$.}\label{tab:generic}
\begin{ruledtabular}
\begin{tabular}{llllllllllll}
 & & & \multicolumn{4}{c}{$\cD/\cD_{\mathrm{form}}$} & \multicolumn{4}{c}{$(\cD-\cB)/(\cD-\cB)_{\mathrm{form}}$} & no interf.\\
case & $\eps$ & res. & $0.5$ & $1$ & $2$ & $3$ & $0.5$ & $1$ & $2$ & $3$ & $1$\\ \hline""" % tuple(sorted(set(round(x['eigvec_overlap'], 6) for x in N3['residues']), key=lambda v: [y['case'] for y in N3['residues'] if abs(y['eigvec_overlap']-v) < 1e-5][0])))
for case in (0, 1, 2):
    for eps in (1e-2, 1e-3, 1e-4, 1e-5):
        rr = [x for x in N3['residues'] if x['case'] == case and abs(x['eps']/eps-1) < 1e-9][0]
        rows = [x for x in N3['asymptotics'] if x['case'] == case and abs(x['eps']/eps-1) < 1e-9]
        pk = lambda v: min(rows, key=lambda x: abs(x['c_over_cD']-v))
        cs_ = [pk(v) for v in (0.5, 1.0, 2.0, 3.0)]
        w('$%d$ & $10^{%d}$ & $%.4f$ & ' % (case+1, round(math.log10(eps)), math.hypot(*rr['a_num_over_alpha_eps'])) + ' & '.join('$%.3f$' % x['D_over_formula'] for x in cs_) + ' & ' + ' & '.join('$%.3f$' % x['F_over_formula'] for x in cs_) + ' & $%.2f$' % pk(1.0)['D_over_incoherent'] + r'\\')
w(r"""\end{tabular}
\end{ruledtabular}
\end{table}
""")
f1 = N4['0.1']; f7 = N4['0.07']
mx = lambda key, sub: max(abs(N4[l][c][key][sub]) if isinstance(N4[l][c][key], dict) else 0 for l in N4 for c in N4[l])
bd = lambda i: max(N4[l][c]['bath_discretization'][i]['dV'] for l in N4 for c in N4[l])
cn = lambda i: max(N4[l][c]['contour_nodes'][i]['dZ'] for l in N4 for c in N4[l])
bad = N4['0.1']['c=4']['bath_discretization'][0]
w(r"""\emph{Convergence of the nonnormal test.} At $\lambda=0.1$ and $0.07$, $c=2$ and $4$, relative to the finest setting: (i) reservoir discretization with recurrence time $1.5\,T$, $2\,T$, $3\,T$ changes $V$ by at most $%s$, $%s$, $%s$; with recurrence time $1.1\,T$ the propagator is wrong by up to $%s$ and the defect at $\lambda=0.1$, $c=4$ by a factor $%.0f$, so that a recurrence time of at least $2\,T$ is necessary (the runs of Table~\ref{tab:nonnormal} use at least $3\,T_{\max}$); (ii) double versus 80-bit propagation: $%s$ in $V$ but up to $%s$ relative in the defect; (iii) contour nodes $n=16,24,32,48$: $%s$, $%s$, $%s$, $%s$; (iv) memory cutoff $L\ge150$: converged to the level of the direct propagation. With the final parameters the reconstruction and the direct propagation agree to a relative $%s$ in the entropy and $%s$ in the defect.

In case~2 of Table~\ref{tab:generic} the convergence at $c=0.5\,c_{\cD}$ is slow. This is consistent with the theorem: the error terms of \eqref{eq:thD} contain a relative correction of order $\mathsf e_s=\lambda^{2c\gamma_s}$ from $a_V=1-\Tr C$ and of order $\lambda^{c\delta\gamma}\log(1/\lambda)$ from the fast--slow overlap, which decay slowly when $c\gamma_s$ and $c\,\delta\gamma$ are small.
""" % (sci(bd(1)), sci(bd(2)), sci(bd(3)), sci(bd(0)), bad['relF']+1, sci(mx('precision_direct_double_vs_80bit', 'dV')), sci(mx('precision_direct_double_vs_80bit', 'relF')), sci(cn(0)), sci(cn(1)), sci(cn(2)), sci(cn(3)), sci(mx('final', 'relD')), sci(mx('final', 'relF'))))
# convergence
e1 = E['0.1']; e7 = E['0.07']
def fd(lst, key, val, out_key): return [r for r in lst if r.get(key) == val][0][out_key]
w(r'''\subsection{Convergence checks}
For Example~\ref{ex:m1} at $\tau=2.2\log(1/\lambda)$, $\lambda=0.1$ and $0.07$ (reference: 80-bit propagation on the finest grid): (i) \emph{reservoir discretization}: with recurrence time $2\pi/\Delta x$ equal to $1.1\,T$, $1.3\,T$, $1.6\,T$, $2\,T$, $3\,T$ the propagator differs from the reference by $%s$, $%s$, $%s$, $%s$, $%s$ at $\lambda=0.1$; (ii) \emph{frequency window}: extending the grid beyond the support of $J$ to $\abs x\le1.6$ changes the result by less than $10^{-14}$ (double-precision runs); (iii) \emph{precision}: double versus 80-bit direct propagation differ by $%s$ and $%s$; (iv) \emph{long-time propagation}: increasing the Chebyshev order by $50\%%$ or the spectral bound by $30\%%$ changes the result by less than $10^{-17}$, and a dense diagonalization of the discretized Hamiltonian agrees with the Chebyshev sum to $%s$; (v) \emph{contour quadrature}: $16$, $24$, $32$, $48$ nodes give differences $%s$, $%s$, $%s$, $%s$ from $320$ nodes; four different circles agree to $10^{-19}$; (vi) \emph{memory cutoff} at $\lambda=0.07$: $L=25$, $50$, $100$, $150$, $200$, $300$ give errors $%s$, $%s$, $%s$, $%s$, $%s$, $%s$ against direct propagation. The bound $2\eps Tt_L$ of Lemma~\ref{lem:truncation} overestimates the actual truncation error by several orders of magnitude.
''' % tuple([sci(e1['bath_discretization'][i]['diff_vs_ref']) for i in (1, 2, 3, 4, 5)] + [sci(e1['precision_direct']['double_vs_longdouble']), sci(e7['precision_direct']['double_vs_longdouble']), sci(e1['long_time_propagation'][3]['diff'])] + [sci(fd(e1['contour_quadrature'], 'n', n, 'diff_vs_n320')) for n in (16, 24, 32, 48)] + [sci(fd(e7['memory_cutoff'], 'L', Lc, 'err_vs_direct')) for Lc in (25, 50, 100, 150, 200, 300)]))
w(r'''The static statements of Sec.~\ref{sec:static} were tested on $%d$ random states with $d\le4$ and excited blocks with eigenvalues down to $10^{-30}$: the identities \eqref{eq:Dkernel}, \eqref{eq:defectkernel}, \eqref{eq:DplusIntegral} hold to at least $15$ digits, and no violation of \eqref{eq:schurbound}--\eqref{eq:mainuniform} was found.
''' % S['random_tests']['ntests'])
# The original data-generation expressions above render ten tables. Retain
# their numeric provenance but distribute only the actual supplementary rows.
rendered='\n'.join(out)
labels=['recon','ep','clusters','fixedcontour','m1','v03','m2','H','nonnormal','generic']

def find_rows(tex,label):
    matches=[t for t in re.findall(r'\\begin\{table\}(?:\[[^\]]*\])?.*?\\end\{table\}',tex,re.S)
             if r'\label{tab:'+label+'}' in t]
    if len(matches)!=1:raise ValueError(f'Expected exactly one table {label}, found {len(matches)}')
    body=matches[0].split(r'\hline',1)[1].split(r'\end{tabular}',1)[0]
    return [' '.join(line.split()) for line in body.splitlines() if '&' in line]

all_rows={label:find_rows(rendered,label) for label in labels}
# The high-precision microscopic bath is a separate, later published dataset.
hp=json.loads((Path(__file__).resolve().parent.parent/'order2_high_precision'/'direct_m2_results.json').read_text(),parse_float=Decimal)
def select_hp(lam,c,N):
    matches=[r for r in hp if Decimal(str(r['lam']))==Decimal(str(lam)) and Decimal(str(r['c']))==Decimal(str(c)) and r['N']==N]
    if len(matches)!=1:raise ValueError((lam,c,N,len(matches)))
    return matches[0]

def scientific(x,places):
    x=Decimal(str(x))
    a,b=f'{x:.{places}e}'.split('e')
    return a+r'\times10^{'+str(int(b))+'}'

rows=[]
for lam,c,N in [('0.050','4.5',8500),('0.050','5.0',8500),('0.050','5.0',10500),('0.050','5.5',8500),('0.035','5.0',17000)]:
    p=select_hp(lam,c,N)
    rows.append(' & '.join([f'${lam}$',f'${c}$',f'${N}$',f'${scientific(p["F_micro"],9)}$',f'${scientific(p["F_pole"],9)}$',f'${scientific(p["relative_discrepancy"],2)}$'])+' \\\\')
all_rows['m2directhp']=rows

# The last table comes from the actual certified error-budget source code.
# Read only the literal numeric expressions I('...') in the validated parts
# dictionary; do not eval the module or hardcode a second numerical budget.
final=(Path(__file__).resolve().parent.parent/'exact_m2_certificate'/'final_error_budget.py').read_text()
parsed=ast.parse(final)
assign=[node for node in parsed.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='parts' for t in node.targets)]
assert len(assign)==1 and isinstance(assign[0].value,ast.Dict)
parts={}
for key,call in zip(assign[0].value.keys,assign[0].value.values):
    assert isinstance(key,ast.Constant) and isinstance(call,ast.Call) and len(call.args)==1 and isinstance(call.args[0],ast.Constant)
    parts[str(key.value)]=Decimal(str(call.args[0].value))
labels_and_keys=[
 ('Continuum to exact $N$-mode evolution','continuum_exact_orthogonal_to_discrete_orthogonal',3),
 ('Double-double propagation/accumulation','DD_recurrence_and_output_arithmetic',3),
 ('Printed-output representation','printed_output_decimal',1),
 ('Decimal input / Hamiltonian parsing','input_DD_parser_to_decimal_H',1),
 ('Exact vs. stored order-two polynomial','exact_p_vs_decimal_p',0),
 ('Bessel-input parsing','Bessel_DD_parser',2),
]
budget=[]
for label,key,places in labels_and_keys:
    budget.append(label+' & $'+scientific(parts[key],places)+'$ \\\\')
small=sum((parts[k] for k in ['Bessel_coefficient_ODE_bound','time_token','ideal_Chebyshev_truncation','42_digit_nodes_and_couplings']),Decimal(0))
assert small<Decimal('1.72e-32'),small
budget.append(r'Bessel residual, series tail, time token, 42-digit input printing & $<1.72\times10^{-32}$ \\')
assert sum(parts.values(),Decimal(0))<Decimal('1.2e-22')
budget.append(r'Total after outward rounding & $<1.2\times10^{-22}$ \\')
all_rows['m2certbudget']=budget

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent.parent/'generated'/'supplementary_tables')
parser.add_argument('--check',action='store_true')
parser.add_argument('--supplement-section',type=Path,help='optional author-supplied supplement LaTeX for numerical row validation')
args=parser.parse_args()
reference=Path(__file__).resolve().parent.parent/'reference'/'published_supplement_table_rows.json'
if args.check:
    expected=json.loads(reference.read_text())
    assert all_rows==expected, 'Supplementary tables diverge from the submitted values'
    print('PASS: all 12 supplementary tables, including exact certificate budget')
if args.supplement_section:
    source=args.supplement_section.read_text()
    for name,lines in all_rows.items():
        actual=find_rows(source,name)
        assert actual==lines,(name,actual[:2],lines[:2])
    print('PASS: all 12 tables match the supplied supplementary manuscript')
args.output_dir.mkdir(parents=True,exist_ok=True)
for label,lines in all_rows.items():
    (args.output_dir/(label+'_rows.tex')).write_text('\n'.join(lines)+'\n')
print('Generated',len(all_rows),'supplementary table-row files')
