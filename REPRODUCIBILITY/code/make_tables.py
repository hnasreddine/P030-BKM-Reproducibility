#!/usr/bin/env python3
"""Reconstruct four independent numerical table datasets printed in the P030 JMP manuscript.

All table entries come from the recorded simulation datasets, not a copy of
pre-existing LaTeX rows. The rounded-SVD m=2 runs are kept separate
from the exactly orthogonal single-point continuum certificate.

Inputs are archived simulation results, not new microscopic propagation.

Use:
 python REPRODUCIBILITY/code/make_tables.py --check
 python REPRODUCIBILITY/code/make_tables.py --output-dir /tmp/p030-tables
"""
from __future__ import annotations
import argparse
from decimal import Decimal
import difflib
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parent.parent
DATA=ROOT/'data'
REF=ROOT/'reference'/'published_table_rows.json'
TAB_END=r'\\'

def load(name: str):
    return json.loads((DATA/name).read_text(encoding='utf-8'),parse_float=Decimal)

def match(rows,**kw):
    out=[x for x in rows if all((Decimal(str(x[k]))==Decimal(str(v))) if isinstance(x[k],(Decimal,int,float)) else (x[k]==v) for k,v in kw.items())]
    if len(out)!=1:raise ValueError(f'Expected unique row {kw}, got {len(out)}')
    return out[0]

def sci(x,digits=1):
    if not isinstance(x,Decimal):x=Decimal(str(x))
    if x==0:return '0'
    norm=f'{x:.{digits}e}'
    mant,exp=norm.split('e')
    # Display e.g. 3.3e-09 as LaTeX 3.3 x 10^-9
    return mant+r'\times10^{'+str(int(exp))+'}'

def math_sci(x,digits=1):return '$'+sci(x,digits)+'$'

def latex_rows():
    A=load('testA_m1.json');B=load('testB_poly_m2.json')
    C=load('testC_EP.json');N2=load('testN2_micro.json')['rows']
    HP=json.loads((ROOT/'order2_high_precision'/'direct_m2_results.json').read_text(),parse_float=Decimal)
    REN=json.loads((ROOT/'order2_high_precision'/'testB_poly_m2.json').read_text(),parse_float=Decimal)
    result={}
    ep=[]
    for lam in ('0.2','0.14','0.1','0.07','0.05'):
        below=match(C,lam=lam,tag='below',c='2.2')
        at=match(C,lam=lam,tag='EP',c='2.2')
        above=match(C,lam=lam,tag='above',c='2.2')
        values=[f'${lam}$',f'${at["delta_star"]:.5f}$',math_sci(at['eig_split'])]
        values.extend(math_sci(x['err_pole']) for x in [below,at,above])
        values.extend([math_sci(at['relD']),math_sci(at['relF'])])
        ep.append(' & '.join(values)+TAB_END)
    result['ep']=ep
    summ=[]
    for m,lam,c in [(1,'0.035','3'),(1,'0.035','5'),(2,'0.070','3'),(2,'0.050','3'),(2,'0.035','3')]:
        z=match(A if m==1 else B,lam=lam,c=c)
        d=z['iso'] if m==1 else z
        ratioD=d['D_micro']/d['D_ren'];ratioF=d['F_micro']/d['F_ren']
        summ.append(f'${m}$ & ${lam}$ & ${c}$ & ${ratioD:.4f}$ & ${ratioF:.4f}$'+TAB_END)
    result['thresholdsummary']=summ
    hp=[]
    for lam,c,n in [('0.050','4.5',8500),('0.050','5.0',8500),('0.050','5.5',8500),('0.035','5.0',17000)]:
        z=match(HP,lam=lam,c=c,N=n)
        ren=match(REN,lam=lam,c=c)
        # F_ren: two-term threshold formula with fitted resonances (VI.14).
        # F_leading: the bare-rate asymptotic specialization (VI.16).
        hp.append(' & '.join((f'${lam}$',f'${c}$',math_sci(z['F_micro'],8),math_sci(z['F_pole'],8),math_sci(z['relative_discrepancy'],2),f'${Decimal(str(z["F_micro"]))/Decimal(str(ren["F_ren"])):.4f}$',f'${Decimal(str(z["F_micro"]))/Decimal(str(z["F_leading"])):.3f}$'))+' '+TAB_END)
    result['m2hp']=hp
    non=[]
    for c in ('1','3','5'):
        r=match(N2,lam='0.035',c=c)['pI']
        ratD=r['D_micro']/r['D_form'];ratF=r['F_micro']/r['F_form'];ratDav=r['D_micro']/r['D_davies']
        dav=math_sci(ratDav,1) if ratDav>=Decimal('100000') else f'${ratDav:.3f}$'
        non.append(f'${c}$ & ${ratD:.4f}$ & ${ratF:.4f}$ & {dav}'+TAB_END)
    result['nonnormal']=non
    return result

def replace_section_rows(section: str, rows: dict[str,list[str]]) -> str:
    """Replace only numerical rows inside four identified TeX tables.

    Scientific prose, captions and disclosure in the author's LaTeX file are
    preserved exactly. The table data, however, are generated from JSON.
    """
    table_pattern=re.compile(r'\\begin\{table\}\[[^\]]*\].*?\\end\{table\}',re.S)
    replaced=[]
    def update(matchobj):
        table=matchobj.group()
        labels=[label for label in rows if r'\label{tab:'+label+'}' in table]
        if not labels:return table
        if len(labels)!=1:raise ValueError(f'Ambiguous table labels: {labels}')
        label=labels[0]
        first=table.index(r'\hline')+len(r'\hline')
        last=table.index(r'\end{tabular}',first)
        existing=table[first:last]
        if existing.count('\\\\')!=len(rows[label]):
            raise ValueError(f'Table row count changed for {label}')
        replaced.append(label)
        return table[:first]+'\n'+'\n'.join(rows[label])+'\n'+table[last:]
    regenerated=table_pattern.sub(update,section)
    if set(replaced)!=set(rows):raise ValueError(f'Missing table labels: {set(rows)-set(replaced)}')
    return regenerated

def validate_section(section_path: Path,rows: dict[str,list[str]]):
    source=section_path.read_text(encoding='utf-8')
    reconstructed=replace_section_rows(source,rows)
    if reconstructed!=source:
        diffs=list(difflib.unified_diff(source.splitlines(),reconstructed.splitlines(),fromfile='source',tofile='data-regenerated',lineterm=''))
        raise ValueError('Manuscript numerical tables do not match public JSON:\n'+'\n'.join(diffs[:120]))
    print('PASS: manuscript numerical section contains all four data-reconstructed tables unchanged')

def validate(rows):
    ref=json.loads(REF.read_text(encoding='utf-8'))
    if rows==ref:return
    report=[]
    for label in sorted(set(ref)|set(rows)):
        if rows.get(label)==ref.get(label):continue
        report.extend(difflib.unified_diff(ref.get(label,[]),rows.get(label,[]),fromfile='submitted/'+label,tofile='generated/'+label,lineterm=''))
    raise ValueError('Table data diverge from the submitted numerical section:\n'+'\n'.join(report))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'generated'/'tables')
    parser.add_argument('--check',action='store_true',help='assert exact submitted table values')
    parser.add_argument('--manuscript-section',type=Path,help='optional author LaTeX section whose rows must match the data')
    parser.add_argument('--rebuild-section',type=Path,help='optional destination to rewrite the four numeric tables into the supplied manuscript section')
    args=parser.parse_args()
    rows=latex_rows()
    if args.check:validate(rows)
    if args.manuscript_section:
        if args.check:validate_section(args.manuscript_section,rows)
        if args.rebuild_section:
            args.rebuild_section.parent.mkdir(parents=True,exist_ok=True)
            args.rebuild_section.write_text(replace_section_rows(args.manuscript_section.read_text(encoding='utf-8'),rows),encoding='utf-8')
            print('Rebuilt four numeric tables in',args.rebuild_section)
    elif args.rebuild_section:
        parser.error('--rebuild-section requires --manuscript-section')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for label,lines in rows.items():
        dest=args.output_dir/(label+'_rows.tex')
        dest.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        print(f'{label}: {len(lines)} rows, {dest}')
    if args.check:print('PASS: all 4 data-backed tables match the submitted manuscript')

if __name__=='__main__':main()
