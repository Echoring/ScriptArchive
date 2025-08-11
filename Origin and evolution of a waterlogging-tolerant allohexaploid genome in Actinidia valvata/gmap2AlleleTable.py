#!/usr/bin/env python3
# Usage: python3 gmap2AlleleTable.py <reference gff3> <gmap gff3> <output Allele.ctg.table>
# Alternative script of ALLHiC gmap2AlleleTable.pl, deal with unstandard reference gff3 format (identifier field name do not match).
from collections import defaultdict
import sys
refgfffile = sys.argv[1]
gmapgff3 = sys.argv[2]
outfile = sys.argv[3]

infordb = defaultdict(str)
with open(gmapgff3, 'r') as qrygff:
    for line in qrygff:
        if line.startswith('#'):
            continue
        tid, annor, field, start, end, score, strand, phase, attributes = line.split()[:9]
        if field == 'gene':
            for attribute in attributes.split(';'):
                if attribute and attribute.split('=')[0] == 'Name':
                    name = attribute.split('=')[1].split('.')[0]
                    infordb[name] += tid + '\t'
                    
with open(refgfffile, 'r') as refgff, open(outfile, 'w') as w:
    for line in refgff:
        if line.startswith('#'):
            continue
        tid, annor, field, start, end, score, strand, phase, attributes = line.split()[:9]
        if field == 'gene':
            for attribute in attributes.split(';'):
                if attribute and attribute.split('=')[0] == 'ID':
                    name = attribute.split('=')[1]
            if name not in infordb:
                continue
            w.write(f'{tid}\t{start}\t{infordb[name]}\n')
            