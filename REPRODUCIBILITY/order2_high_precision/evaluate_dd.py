import sys,json,os,mpmath as mp
sys.path.insert(0,os.path.dirname(__file__)); import p030lib as L
mp.mp.dps=100
path=sys.argv[1]; lam=float(sys.argv[2]); c=float(sys.argv[3]);
lines=open(path).read().splitlines(); raw=[[mp.mpf(x) for x in line.split()] for line in lines]; rows=[[r[0]+r[1],r[2]+r[3]] for r in raw]
V=mp.matrix([[mp.mpc(*rows[0]),mp.mpc(*rows[2])],[mp.mpc(*rows[1]),mp.mpc(*rows[3])]])
C0=mp.matrix([[mp.mpf('.2'),0],[0,mp.mpf('.2')]]);b=mp.matrix([[0,mp.mpf('.1')]])
D,B,F=L.entropy_and_defect(V,C0,b,dps=110)
prev=json.load(open(os.path.join(os.path.dirname(__file__),'testB_poly_m2.json')))
old=next((r for r in prev if r['lam']==lam and r['c']==c),None)
print('V11',mp.nstr(V[0,0],24),'V12',mp.nstr(V[0,1],24),'V22',mp.nstr(V[1,1],24))
print('D',mp.nstr(D,40),'B',mp.nstr(B,40),'F',mp.nstr(F,40))
if old:
 print('old_F',old['F_micro'],'pole',old['F_pole'],'asym',old['F_asym'],'quad/oldF',mp.nstr(F/old['F_micro'],15),'quad/pole',mp.nstr(F/old['F_pole'],15))
