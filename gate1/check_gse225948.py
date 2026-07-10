import os, glob, subprocess
b = 'data/GSE225948'
gz = glob.glob(b + '/**/*.gz', recursive=True)
ok = sum(1 for f in gz if subprocess.run(['gzip', '-t', f], capture_output=True).returncode == 0)
print('gz files:', len(gz), '| gzip OK:', ok, '/ 36')
# 显示已完成的样本目录
done = set(os.path.basename(os.path.dirname(f)) for f in gz)
print('样本目录数:', len(done))
for d in sorted(done):
    fs = [os.path.basename(f) for f in glob.glob(os.path.join(b, d, '*.gz'))]
    print(' ', d, '->', len(fs), 'files')
