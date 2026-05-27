#/usr/bin/env python3
# Usage: python3 rearrange_assembly.py hic_assembly quartet_assembly
# Rearrange each chr in hic_assembly to match the quartet_assembly.
import sys

hicassembly = sys.argv[1]
qamassembly = sys.argv[2]

with open(hicassembly, 'r') as f:
    hicseqdict = {}
    hicreturndict = {}
    hicseqorder = []
    for line in f:
        if line.startswith('>'):
            sid = line.split()[0][1:]
            sidx = int(line.split()[1])
            hicseqdict[sidx] = sid
            hicreturndict[sid] = sidx
            print(line.strip())
        else:
            hicseqorder.append(['-' + hicseqdict[abs(int(i))] if i[0] == '-' else hicseqdict[abs(int(i))] for i in line.strip().split()])
            
with open(qamassembly, 'r') as f:
    qamseqdict = {} 
    qamseqorder = []
    for line in f:
        if line.startswith('>'):
            sid = line.split()[0][1:]
            sidx = int(line.split()[1])
            qamseqdict[sidx] = sid
        else:
            qamseqorder.append(['-' + qamseqdict[abs(int(i))] if i[0] == '-' else qamseqdict[abs(int(i))] for i in line.strip().split()])

def rearrange_sublist(sub1, sub2):
    sub1_dict = {elem.strip('-'): elem for elem in sub1}
    sub2_dict = {elem.strip('-'): elem for elem in sub2}
    
    rearranged = []
    remaining = []
    
    def flipminus(elem):
        if elem[0] == '-':
            return elem.strip('-')
        else:
            return '-' + elem
    
    for elem in sub2:
        key = elem.strip('-')
        if key in sub1_dict:
            rearranged.append(sub2_dict[key])
            
    for elem in sub1:
        key = elem.strip('-')
        if key not in sub2_dict:
            remaining.append(sub1_dict[key])
    
    while len(remaining) > 0:
        for elem in remaining:
            rm = False
            previousleft = sub1[sub1.index(elem)-1]
            if previousleft in remaining:
                continue
            elif previousleft in rearranged:
                leftnewidx = rearranged.index(previousleft)
                rearranged.insert(leftnewidx+1, elem)
                rm = True
            elif flipminus(previousleft) in rearranged:
                leftnewidx = rearranged.index(flipminus(previousleft))
                rearranged.insert(leftnewidx, flipminus(elem))
                rm = True
            break
        if rm == True:
            remaining.remove(elem)

    
    return rearranged

def find_best_match(sub1, list2):
    best_match = None
    max_common = 0
    
    for sub2 in list2:
        if len(sub2) > 1:
            sub1_set = {elem.strip('-') for elem in sub1}
            sub2_set = {elem.strip('-') for elem in sub2}
            common = len(sub1_set.intersection(sub2_set))
            
            if common > max_common:
                max_common = common
                best_match = sub2
    
    return best_match

outseqorder = []
for sub1 in hicseqorder:
    best_sub2 = find_best_match(sub1, qamseqorder)
    if best_sub2:
        rearranged_sub1 = rearrange_sublist(sub1, best_sub2)
        outseqorder.append(rearranged_sub1)
    else:
        outseqorder.append(sub1)

for sub in outseqorder:
    print(' '.join(['-' + str(hicreturndict[elem.strip('-')]) if elem[0] == '-' else str(hicreturndict[elem.strip('-')]) for elem in sub]))
