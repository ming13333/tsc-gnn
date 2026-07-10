import os, glob, subprocess

base = 'data/GSE174574'
gz = [f for f in glob.glob(base + '/**/*.gz', recursive=True)
      if '.part' not in f and '.tmp' not in f and '.hdr' not in f]
ok = 0
for f in gz:
    if subprocess.run(['gzip', '-t', f], capture_output=True).returncode == 0:
        ok += 1
print('完整gz文件:', len(gz), '| gzip OK:', ok, '/ 18')
for d in sorted(glob.glob(base + '/GSM*')):
    name = os.path.basename(d)
    m = glob.glob(d + '/*matrix*.gz')
    if m:
        mm = m[0]
        sz = os.path.getsize(mm)
        g = subprocess.run(['gzip', '-t', mm], capture_output=True).returncode == 0
        print(f'  {name}: matrix={sz}/{35916235} {"OK" if g else "partial"}')
    else:
        print(f'  {name}: matrix MISSING')
