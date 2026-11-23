"""Rebuild a stable SHA256 manifest excluding generated binaries and itself."""
import hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent
allow={'.py','.cpp','.md','.txt','.log','.status'}
files={}
for f in sorted(p.rglob('*')):
 if f.is_file() and f.suffix in allow and f.name!='SHA256_MANIFEST.json':
  files[str(f.relative_to(p))]={"sha256":hashlib.sha256(f.read_bytes()).hexdigest(),"bytes":f.stat().st_size}
(p/'SHA256_MANIFEST.json').write_text(json.dumps(files,indent=2,ensure_ascii=False)+'\n')
print('manifest entries',len(files),'total bytes',sum(x['bytes'] for x in files.values()))
