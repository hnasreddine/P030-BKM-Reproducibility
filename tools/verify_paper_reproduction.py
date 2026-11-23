#!/usr/bin/env python3
"""Clean-environment regression for the 4 published numerical tables and 2 figures.

This test does not recertify the separate exact-m=2 continuum theorem. It
reconstructs publication outputs from independently stored public datasets.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/'REPRODUCIBILITY'/'code'
DATA=ROOT/'REPRODUCIBILITY'/'data'
HP=ROOT/'REPRODUCIBILITY'/'order2_high_precision'/'direct_m2_results.json'
REF=ROOT/'REPRODUCIBILITY'/'reference'/'published_table_rows.json'


def run(cmd:list[str],cwd:Path)->str:
    p=subprocess.run(cmd,cwd=cwd,text=True,capture_output=True)
    if p.returncode:
        raise AssertionError(f'Failed {cmd!r} in {cwd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}')
    return p.stdout


def datum(rows: list[dict],**kw):
    matches=[x for x in rows if all(Decimal(str(x[k]))==Decimal(str(v)) for k,v in kw.items())]
    assert len(matches)==1,kw
    return matches[0]


def validate_science():
    ep=json.loads((DATA/'testC_EP.json').read_text())
    n2=json.loads((DATA/'testN2_micro.json').read_text())
    m2=json.loads(HP.read_text())
    assert len([x for x in ep if x['c']==2.2 and x['tag']=='EP'])==5
    # Historical rounded-SVD bath, deliberately not the exact-bath certificate.
    exact_grid1=datum(m2,lam=0.05,c=5,N=8500)
    exact_grid2=datum(m2,lam=0.05,c=5,N=10500)
    a=Decimal(exact_grid1['F_micro']);b=Decimal(exact_grid2['F_micro'])
    grid_rel=abs(a-b)/a
    assert grid_rel<Decimal('2e-14'),grid_rel
    # Historical m=2 point differs from separately certified exact orthogonal bath.
    assert abs(a-Decimal('3.96903554965836e-41'))/a>Decimal('1e-10')
    rows=n2['rows']
    assert len([x for x in rows if x['lam']==0.035])>=3
    # The four unresolved nonnormal-defect points omitted from the figure.
    sel=[x for x in rows if x['lam']==0.035]
    unresolved=[x for x in sel if not ((x['err_pole']/x['fast_amp'])**2*x['pI']['F_fast']/x['pI']['F_slow']<1e-3)]
    assert len(unresolved)==4, len(unresolved)
    print('PASS: source-data point mapping, historical-vs-exact distinction, grid convergence, unresolved-point masking')


def main():
    validate_science()
    with tempfile.TemporaryDirectory(prefix='p030-paper-repro-') as tmp:
        temp=Path(tmp)
        output=temp/'newly-generated'
        print(run([sys.executable,str(CODE/'make_tables.py'),'--check','--output-dir',str(output/'tables')],cwd=temp).strip())
        # Explicitly test the generated content against an adversarially corrupted
        # source copy. This ensures no formatting-only test can pass after data drift.
        work=temp/'mutation'
        code=work/'code';data=work/'data';hp=work/'order2_high_precision';ref=work/'reference'
        for directory in [code,data,hp,ref]:directory.mkdir(parents=True)
        shutil.copy2(CODE/'make_tables.py',code/'make_tables.py')
        shutil.copy2(DATA/'testA_m1.json',data/'testA_m1.json')
        shutil.copy2(DATA/'testB_poly_m2.json',data/'testB_poly_m2.json')
        shutil.copy2(DATA/'testC_EP.json',data/'testC_EP.json')
        shutil.copy2(DATA/'testN2_micro.json',data/'testN2_micro.json')
        shutil.copy2(HP,hp/'direct_m2_results.json')
        shutil.copy2(REF,ref/'published_table_rows.json')
        f=data/'testN2_micro.json';d=json.loads(f.read_text())
        point=next(x for x in d['rows'] if x['lam']==0.035 and x['c']==1.0)
        point['pI']['F_micro']*=1.5
        f.write_text(json.dumps(d))
        bad=subprocess.run([sys.executable,str(code/'make_tables.py'),'--check','--output-dir',str(output/'bad')],cwd=temp,capture_output=True,text=True)
        assert bad.returncode!=0 and 'diverge' in bad.stderr,'A corrupted source was not detected'
        print('PASS: injected numerical corruption is rejected by the publication-data gate')
        for script,name in [('make_figures.py','fig_thresholds.pdf'),('make_figure_nonnormal.py','fig_nonnormal.pdf')]:
            dest=output/'figures'/name
            print(run([sys.executable,str(CODE/script),'--output',str(dest)],cwd=temp).strip())
            content=dest.read_bytes()
            assert content.startswith(b'%PDF') and len(content)>15000,(name,len(content))
            print('PASS:',name,'regenerated',len(content),'bytes from clean outside-repository cwd')
        assert len(list((output/'tables').glob('*.tex')))==4
    print('PASS: complete four-table/two-figure reproduction regression')

if __name__=='__main__':main()
