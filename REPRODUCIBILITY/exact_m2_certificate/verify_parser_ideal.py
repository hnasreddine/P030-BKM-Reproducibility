"""Exact-rational audit of all recorded decimal->DD input conversions.
Requires the same libquadmath C++ parser as the original source; no MPFR assumptions.
"""
from pathlib import Path
from fractions import Fraction
import decimal, subprocess, re
P=Path(__file__).resolve().parent
SOURCE=(Path(__file__).resolve().parent / 'cheb_dd_lambda_corrected.cpp')
src=SOURCE.read_text()
prefix=src[:src.index('int main(')]
program=prefix+r'''
#include <cstring>
int main(int argc,char**argv){
 FILE*f=fopen(argv[1],"r");
 char R[256], T[256],lam[256]; int N,K;
 fscanf(f,"%d %d %255s %255s %255s",&N,&K,R,T,lam);
 auto emit=[](const char* v){DD a=ddparse(v);printf("%s %a %a\n",v,a.h,a.l);};
 emit(lam);
 char token[256];
 for (int j=0;j<4*N+K+1;j++){if(fscanf(f,"%255s",token)!=1)return 8;emit(token);}
 fclose(f);
}
'''
(P/'verify_input_parser.cpp').write_text(program)
subprocess.run(['g++','-O3','-std=c++17','-fno-fast-math','-ffp-contract=off',str(P/'verify_input_parser.cpp'),'-lquadmath','-o',str(P/'verify_input_parser')],check=True)
p=subprocess.run([str(P/'verify_input_parser'),str(P/'input_005_5_8500_ideal_m2.txt')],text=True,capture_output=True,check=True)
u=Fraction(1,2**53); eps100=Fraction(1,2**100); eta=Fraction(1,10**300)
mxr=Fraction(0);maxrline='';bad=0; zero_high=0;maxoverlap=Fraction(0); skipped_tiny=0
for row in p.stdout.splitlines():
    tok,hs,ls=row.split(); h=float.fromhex(hs); l=float.fromhex(ls)
    m=re.search(r'[eE]([-+]?\d+)$',tok)
    if m and m.group(1).startswith('-') and len(m.group(1))>4:
        assert h==0 and l==0, ('expected binary zero for ultratiny token',tok)
        skipped_tiny+=1
        continue
    exponent=int(m.group(1)) if m else 0
    if exponent < -1000:
        assert h==0 and l==0, ('expected binary zero for tiny token',tok)
        skipped_tiny+=1
        continue
    exact=Fraction(decimal.Decimal(tok))
    fracH=Fraction.from_float(h);fracL=Fraction.from_float(l)
    normerr=abs(fracH+fracL-exact)
    allowance=eps100*abs(exact)+eta
    if normerr>allowance:bad+=1
    if h==0 and l!=0:zero_high+=1
    if h!=0:
        ratio=abs(fracL)/abs(fracH)
        if ratio>maxoverlap:maxoverlap=ratio
    if exact!=0 and abs(exact)>=eta:
        rel=normerr/(eps100*abs(exact))
        if rel>mxr:mxr=rel;maxrline=tok
print('PARSED_TOKENS',len(p.stdout.splitlines()))
print('TINY_TOKENS_BELOW_1e-300',skipped_tiny)
print('EXACT_RATIONAL_CONVERSION_BOUND_VIOLATIONS',bad)
print('ZERO_HIGH_NONZERO_LOW',zero_high)
print('MAX_REL_ERROR_OVER_2^-100',float(mxr))
print('MAX_OVERLAP_OVER_u',float(maxoverlap/u))
print('WORST_INPUT',maxrline)
assert bad==0 and maxoverlap<=u and zero_high==0
