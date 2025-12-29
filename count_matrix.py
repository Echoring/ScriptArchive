#/usr/bin/env python3
# Usage: ... <countfile1> <countfile2> ...
# Convert FeatureCount result to matrix, each file a sample.
import sys
import collections
countfilelist = sys.argv[1:]

titleline = ''
countdict = collections.defaultdict(list)
for countfile in countfilelist:
    titleline += f'\t{countfile}'
    with open(countfile) as r:
        for line in r:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            gid = line.split()[0]
            count = line.split()[6]
            countdict[gid].append(count)

print(titleline)
for gid, linelist in countdict.items():
    line = '\t'.join(linelist)
    print(f'{gid}\t{line}')
    
with open('sample.txt', 'w') as w:
    w.write('\tcondition\n')
    for countfile in countfilelist:
        countfilec = countfile.replace('-', '.')
        w.write(f'{countfilec}\tg\n')