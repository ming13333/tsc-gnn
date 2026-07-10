import os, glob, subprocess

base = 'data/GSE174574'
tmp = [f for f in glob.glob(base + '/**/*', recursive=True)
       if f.endswith('.tmp') or f.endswith('.hdr')]
gz = [f for f in glob.glob(base + '/**/*.gz', recursive=True) if '.part' not in f]
print('正在下载的临时文件数:', len(tmp))
for t in tmp[:6]:
    print('  ', t, os.path.getsize(t) if os.path.exists(t) else 'gone')
ok = sum(1 for f in gz if subprocess.run(['gzip', '-t', f], capture_output=True).returncode == 0)
print('完整gz OK:', ok, '/ 18')
for d in sorted(glob.glob(base + '/GSM*')):
    name = os.path.basename(d)
    ms = glob.glob(d + '/*matrix*.gz')
    if ms:
        m = ms[0]
        g = subprocess.run(['gzip', '-t', m], capture_output=True).returncode == 0
        print(f'  {name} matrix={os.path.getsize(m)} {"OK" if g else "partial"}')
