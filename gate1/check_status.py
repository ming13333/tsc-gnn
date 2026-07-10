import os, glob, subprocess
from collections import defaultdict

base = 'data/GSE174574'
parts = sorted(glob.glob(base + '/**/*.part*', recursive=True))
gz = [f for f in glob.glob(base + '/**/*.gz', recursive=True) if '.part' not in f]
ok = sum(1 for f in gz if subprocess.run(['gzip', '-t', f], capture_output=True).returncode == 0)
print('part分块数:', len(parts), '| 完整gz OK:', ok, '/ 18')

d = defaultdict(list)
for p in parts:
    d[os.path.dirname(p)].append((os.path.basename(p), os.path.getsize(p)))
for k in sorted(d):
    lst = d[k]
    maxpart = max(lst, key=lambda x: int(x[0].split('part')[-1]))
    print(f'  {os.path.basename(k)}: {len(lst)}块, 最大块={maxpart[0]} 大小={maxpart[1]}')
