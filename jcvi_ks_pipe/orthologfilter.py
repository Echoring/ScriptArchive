#!/usr/bin/env python3
import sys
import collections
inanchorAB = sys.argv[1]
inanchorBA = sys.argv[2]
mode = sys.argv[3]

anchordict = collections.defaultdict(list)
used = set()

with open(inanchorAB, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        sidA, sidB, score = line.strip().split()
        if sidA != sidB:
            anchordict[sidA].append((sidB, int(score.strip('L'))))
        
with open(inanchorBA, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        sidB, sidA, score = line.strip().split()
        if sidA != sidB:
            anchordict[sidB].append((sidA, int(score.strip('L'))))
        
for sid, anchorlist in anchordict.items():
    if mode == 'single':
        if len(anchorlist) > 1:
            continue
        target = anchorlist[0][0]
        if target not in anchordict:
            continue
        if len(anchordict[target]) == 1 and sid == anchordict[target][0][0]:
            if f'{sid}\t{target}' in used:
                continue
            print(f'{sid}\t{target}')
            used.add(f'{target}\t{sid}')
    elif mode == 'RBH':
        target = sorted(anchorlist, key=lambda x: x[1], reverse=True)[0][0]
        if target not in anchordict:
            continue
        targetstarget = sorted(anchordict[target], key=lambda x: x[1], reverse=True)[0][0]
        if sid == targetstarget:
            if f'{sid}\t{target}' in used:
                continue
            print(f'{sid}\t{target}')
            used.add(f'{target}\t{sid}')
    elif mode == 'all':
        for anchor in anchorlist:
            target = anchor[0]
            if f'{sid}\t{target}' in used:
                continue
            print(f'{sid}\t{target}')
            used.add(f'{target}\t{sid}')
            used.add(f'{sid}\t{target}')
    else:
        print(f'Unknown mode: {mode}, need to be "single", "RBH" or "all"')
        sys.exit(1)