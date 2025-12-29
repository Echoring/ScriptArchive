#!/usr/bin/env python3
# Usage: python3 globle_rescue_alt.py <input Allele.ctg.table> <all contig vs ref paf> <output Allele.ctg.new.txt>
# alternative script of ALLHiC's globle rescue, resolve unitigs that do not contain genes annotated on reference
import sys
intable = sys.argv[1]
inalignment = sys.argv[2]
outtable = sys.argv[3]


# Read init Allele.ctg.table, get contig-chromosome assignment
oldtargetdict = {}
linelist = []
with open(intable, 'r') as allele:
    for line in allele:
        linelist.append(line)
        refid = line.split()[0]
        tiglist = line.split()[2:]
        for tigid in tiglist:
            oldtargetdict[tigid] = refid

# Read alignment, add unassigned contigs            
uselist = []
addlinelist = []
with open(inalignment, 'r') as aln:
    for line in aln:
        tid = line.split()[0]
        if tid not in uselist:
            rid = line.split()[5]
            start = line.split()[7]
            if tid not in oldtargetdict:
                addlinelist.append(f'{rid}\t{start}\t{tid}\n')
            uselist.append(tid)
linelist.extend(addlinelist)
linelist.sort(key=lambda x:(x.split()[0], int(x.split()[1])))

with open(outtable, 'w') as w:
    for line in linelist:
        w.write(line)