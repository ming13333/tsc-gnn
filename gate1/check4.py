import os, glob, subprocess, time

base = 'data/GSE174574'
for d in sorted(glob.glob(base + '/GSM*')):
    name = os.path.basename(d)
    m = glob.glob(d + '/*matrix*.gz')
    if m:
        mm = m[0]
        sz = os.path.getsize(mm)
        g = subprocess.run(['gzip', '-t', mm], capture_output=True).returncode == 0
        print(f'{name}: matrix={sz} {"OK" if g else "partial"}')
