#!/usr/bin/env python3
"""Validate contents of a self-contained public P030 scientific release."""
import hashlib,json,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
def checksum(path):return hashlib.sha256(path.read_bytes()).hexdigest()
known_inputs={
 'REPRODUCIBILITY/exact_m2_certificate/input_005_5_8500.txt':'bebc1a3b3c2de88c469835d86f3596c96dedcf8f483a09e87cfa55732fea9b96',
 'REPRODUCIBILITY/exact_m2_certificate/input_005_5_8500_ideal_m2.txt':'708414eef62d3546257d1e0202d220246837e32c401f63574ffe38157a1726c1',
}
errors=[]
for rel,expected in known_inputs.items():
 path=R/rel
 if not path.is_file():errors.append('missing input: '+rel)
 elif checksum(path)!=expected:errors.append('SHA256 mismatch: '+rel)
required=['LICENSE','DATA_LICENSE.md','.zenodo.json','CITATION.cff','README.md',
 'RELEASE_MANIFEST.json','SHA256SUMS.txt',
 'REPRODUCIBILITY/exact_m2_certificate/CERTIFICATE_IDEAL_ORDER2_2026-11-22.md',
 'REPRODUCIBILITY/exact_m2_certificate/PROOF_DD_AND_NUMERICAL_CERTIFICATE_2026-11-22.md',
 'REPRODUCIBILITY/exact_m2_certificate/continuum_proof/CONTINUUM_BOUND_AND_TRANSFER.md',
 'REPRODUCIBILITY/exact_m2_certificate/final_error_budget.py',
 'REPRODUCIBILITY/exact_m2_certificate/certified_BKM_tight.py',
 'REPRODUCIBILITY/exact_m2_certificate/strict_dd_audit.cpp',
 'REPRODUCIBILITY/exact_m2_certificate/out_ideal_m2_lean.txt',
 'tools/reproduce_certified_point.sh']
for rel in required:
 if not (R/rel).is_file():errors.append('missing: '+rel)
manifest_path=R/'RELEASE_MANIFEST.json'
if manifest_path.is_file():
 m=json.loads(manifest_path.read_text())
 listed={}
 for entry in m.get('files',[]):
  rel=entry['path']; listed[rel]=entry
  path=R/rel
  if not path.is_file():errors.append('manifest missing: '+rel)
  elif path.stat().st_size != entry['bytes'] or checksum(path)!=entry['sha256']:
   errors.append('manifest mismatch: '+rel)
 actual={p.relative_to(R).as_posix() for p in R.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc'}
 expected=set(listed)|{'RELEASE_MANIFEST.json','SHA256SUMS.txt'}
 for rel in sorted(actual-expected):errors.append('unmanifested: '+rel)
 for rel in sorted(expected-actual):errors.append('absent: '+rel)
 if m.get('version')!='1.0.0':errors.append('version mismatch')
 sha_path=R/'SHA256SUMS.txt'
 if sha_path.is_file():
  lines=sha_path.read_text().splitlines()
  pairs={v.split('  ',1)[1]:v.split('  ',1)[0] for v in lines if '  ' in v}
  if set(pairs)!=set(listed):errors.append('SHA256SUMS inventory mismatch')
  else:
   for rel,digest in pairs.items():
    if digest!=listed[rel]['sha256']:errors.append('SHA256SUMS hash mismatch: '+rel)
# Publication hygiene: distinguish mathematical verification scripts from IA reports.
banned={'AUDIT','history','P030_final_audit'}
for path in R.rglob('*'):
 if not path.is_file():continue
 if set(path.relative_to(R).parts)&banned:errors.append('internal file in public package: '+str(path.relative_to(R)))
 if path.name in ('AUDIT_FINAL.md','PROMPT_CLAUDE_CODE_P030_AUDIT_FINAL.md',
                  'README_PUBLISH_FR.md','DATA_AVAILABILITY_FOR_MANUSCRIPT.md'):
  errors.append('non-public document included: '+path.name)
for file in ['REPRODUCIBILITY/exact_m2_certificate/certified_BKM_tight.py',
             'REPRODUCIBILITY/exact_m2_certificate/final_error_budget.py']:
 path=R/file
 if path.is_file():
  src=path.read_text()
  if file.endswith('certified_BKM_tight.py') and "'out_ideal_m2_lean.txt'" not in src:
   errors.append('incorrect BKM output default')
  if file.endswith('final_error_budget.py') and '42_digit_nodes_and_couplings' not in src:
   errors.append('omitted grid-printing contribution')
try:
 meta=json.loads((R/'.zenodo.json').read_text())
 if meta.get('upload_type')!='software' or meta.get('version')!='1.0.0':errors.append('incorrect Zenodo metadata')
except Exception as e:errors.append('unreadable Zenodo metadata: '+str(e))
if errors:
 print('FAILED PUBLIC RELEASE CHECK:')
 for error in errors[:50]:print(' -',error)
 sys.exit(1)
print('PASS: clean public scientific content; required proofs and code present; two original inputs match SHA256')
print('PASS: all',len(m['files']),'manifest entries and SHA256SUMS checked')
print('NOTE: this is an offline integrity check, not an actual GitHub release or Zenodo deposition')
