"""Figure 2: nonnormal Davies generator, entropy and BKM defect versus c."""
import argparse
from pathlib import Path
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=ROOT/'generated'/'figures'/'fig_nonnormal.pdf')
args=parser.parse_args()
args.output.parent.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'axes.linewidth': 0.6, 'lines.linewidth': 1.2, 'mathtext.fontset': 'cm'})
N2 = json.loads((ROOT/'data'/'testN2_micro.json').read_text()); lam = 0.035
INK, BLUE, ORANGE, GRAY = '#111111', '#0072B2', '#D55E00', '#888888'
rows = [x for x in N2['rows'] if abs(x['lam']-lam) < 1e-9]
c = np.array([x['c'] for x in rows]); g = lambda k: np.array([x['pI'][k] for x in rows])
ok = np.array([(x['err_pole']/x['fast_amp'])**2*x['pI']['F_fast']/x['pI']['F_slow'] < 1e-3 for x in rows])
fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.9))
for a, mic, t1, t2, tot, dav, cth, title, okk in (
        (ax[0], g('D_micro'), g('D_fastdir'), g('D_slowdir'), g('D_form'), g('D_davies'), 2.0, r'(a) entropy $\mathrm{D}$   ($c_{\mathrm{D}}=2$)', np.ones(len(c), bool)),
        (ax[1], g('F_micro'), g('F_fast'), g('F_slow'), g('F_form'), g('F_davies'), 4.0, r'(b) defect $\mathrm{D}-\mathrm{B}$   ($c_{\rm BKM}=4$)', ok)):
    a.semilogy(c, t1, color=BLUE, ls='--', label=r'fast term ($w_\perp$)')
    a.semilogy(c, t2, color=ORANGE, ls='-.', label=r'slow-direction term')
    a.semilogy(c, tot, color=INK, lw=0.7, label='two-term formula')
    a.semilogy(c, dav, color=GRAY, ls=':', label='Davies semigroup')
    a.semilogy(c[okk], mic[okk], 'o', ms=3.2, mfc='white', mec=INK, mew=0.8, label='microscopic (direct)')
    a.axvline(cth, color=INK, lw=0.5, alpha=0.6); a.set_title(title, fontsize=8)
    a.set_ylim(min(mic[okk].min(), tot.min())*1e-3, mic.max()*10)
    a.grid(True, which='major', lw=0.3, color='#dddddd'); a.set_axisbelow(True)
    a.set_xlabel(r'horizon parameter $c$   ($\tau=c\,\log(1/\lambda)$)')
ax[0].legend(frameon=False, fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(args.output); plt.close(fig)
print('ok', int((~ok).sum()), 'unresolved defect points omitted')
