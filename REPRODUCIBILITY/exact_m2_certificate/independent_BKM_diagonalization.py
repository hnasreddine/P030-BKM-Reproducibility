"""Independent high-precision numerical diagnostic: diagonalize the 3x3 density
matrix and its dephased 2x2 excited block. Not a proof certificate."""
from pathlib import Path
import sys,mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import bkm_evaluation as L
mp.mp.dps=160
p=Path(__file__).resolve().parent/'out_ideal_m2_lean.txt'
rows=[[mp.mpf(v) for v in line.split()] for line in p.read_text().splitlines()]
z=lambda r:mp.mpc(r[0]+r[1],r[2]+r[3])
V=mp.matrix([[z(rows[0]),z(rows[2])],[z(rows[1]),z(rows[3])]])
C0=mp.diag([mp.mpf('.2'),mp.mpf('.2')]);b=mp.matrix([[0,mp.mpf('.1')]])
D,B,F=L.entropy_and_defect(V,C0,b,dps=160)
print('D=',mp.nstr(D,55))
print('B=',mp.nstr(B,55))
print('D-B=',mp.nstr(F,55))
assert mp.mpf('3.9670526841476009e-41') <F<mp.mpf('3.9710197867922034e-41')
print('INDEPENDENT_DIAGONALIZATION_WITHIN_INTERVAL PASS')
