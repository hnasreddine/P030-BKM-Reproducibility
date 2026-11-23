"""Figure 1: entropy and BKM defect versus horizon c (m=1 and m=2)."""
import argparse
from pathlib import Path
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=ROOT/'generated'/'figures'/'fig_thresholds.pdf')
args=parser.parse_args()
args.output.parent.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'axes.linewidth': 0.6, 'lines.linewidth': 1.2, 'mathtext.fontset': 'cm'})
A = json.loads((ROOT/'data'/'testA_m1.json').read_text()); B = json.loads((ROOT/'data'/'testB_poly_m2.json').read_text())
INK, BLUE, ORANGE, GRAY = '#111111', '#0072B2', '#D55E00', '#888888'
fig, ax = plt.subplots(2, 2, figsize=(6.8, 5.0), sharex='col')
def panel(a, rows, key, lam, title, cth, floor=None):
    r = [x for x in rows if abs(x['lam']-lam) < 1e-9]
    c = np.array([x['c'] for x in r])
    g = (lambda x: x[key]) if key in r[0] else (lambda x: x['iso'][key])
    pre = 'D' if key.startswith('D') else 'F'
    mic = np.array([g(x) for x in r]); fast = np.array([(x if pre+'_fast' in x else x['iso'])[pre+'_fast'] for x in r]); slow = np.array([(x if pre+'_slow' in x else x['iso'])[pre+'_slow'] for x in r])
    dav = np.array([(x if pre+'_davies' in x else x['iso'])[pre+'_davies'] for x in r])
    ok = np.ones(len(r), bool) if floor is None else np.array([np.exp(-floor*x['tau']) > 5e-17 for x in r])
    a.semilogy(c, fast, color=BLUE, ls='--', label='fast-resonance term')
    a.semilogy(c, slow, color=ORANGE, ls='-.', label='slow-resonance term')
    a.semilogy(c, dav, color=GRAY, ls=':', label='Davies semigroup')
    a.semilogy(c[ok], mic[ok], 'o', ms=3.2, mfc='white', mec=INK, mew=0.8, label='microscopic (direct)')
    a.axvline(cth, color=INK, lw=0.5, alpha=0.6); a.set_title(title, fontsize=8)
    lo = min(mic[ok].min(), max(fast.min(), slow.min())); a.set_ylim(lo*1e-2, max(mic.max(), fast.max())*10)
    a.grid(True, which='major', lw=0.3, color='#dddddd'); a.set_axisbelow(True)
panel(ax[0, 0], A, 'D_micro', 0.035, r'(a) $m=1$: entropy $\mathrm{D}$, $\lambda=0.035$   ($c_{\mathrm{D}}=2$)', 2.0)
panel(ax[1, 0], A, 'F_micro', 0.035, r'(b) $m=1$: defect $\mathrm{D}-\mathrm{B}$, $\lambda=0.035$   ($c_{\rm BKM}=4$)', 4.0, floor=2.0)
panel(ax[0, 1], B, 'D_micro', 0.035, r'(c) $m=2$: entropy $\mathrm{D}$, $\lambda=0.035$   ($c_{\mathrm{D}}=2$)', 2.0)
panel(ax[1, 1], B, 'F_micro', 0.035, r'(d) $m=2$: defect $\mathrm{D}-\mathrm{B}$, $\lambda=0.035$   ($c_{\rm BKM}=4$)', 4.0, floor=3.0)
for a in ax[1]: a.set_xlabel(r'horizon parameter $c$   ($\tau=c\,\log(1/\lambda)$)')
ax[0, 0].legend(frameon=False, fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(args.output); plt.close(fig)
print('ok')
