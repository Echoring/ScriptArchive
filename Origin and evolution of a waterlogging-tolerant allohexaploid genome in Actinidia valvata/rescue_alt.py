#!/usr/bin/env python3
# Alternative script for ALLHiC_rescue
# args: 
# $1=infasta (seq.fasta/ChrXX.fa, from partition_gmap.pl/partition_gmap.py)
# $2=inbam (prunning.sub.bam/ChrXX.bam, from partition_gmap.pl/partition_gmap.py)
# $3=mintiglen (discard contig < this bp)
# $4=maxoverlapalignlen (consider as overlap if alignment length >= this bp, tips 50000)
# $5=thread (max thread for minimap2 and samtools)
# $6=allgroup (prunning.counts_GATC.txt, from allhic partition)
# $7~=ingroups (prunning.counts_GATC.*g*.txt, from allhic partition)
# Generate: (inter) tmp.avsa.paf tmp.mrq.sam tmp.prunedhiclink (result) group*.txt
import collections
import sys
import subprocess
def readFastaAsDict(fastafile):
    fastaDict = {}
    allline = open(fastafile, 'r').read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

infasta = sys.argv[1]
inbam = sys.argv[2]
mintiglen = sys.argv[3]
maxoverlapalignlen = sys.argv[4]
thread = sys.argv[5]
allgroup = sys.argv[6]
ingroup = sys.argv[7:]

print('Reading FASTA...')
fastadict = readFastaAsDict(infasta)

# Overlap alignment, those had large overlap are considered allelic contigs and should not be clustered into the same group
print('Generating all vs all overlap alignment...')
overlapdict = collections.defaultdict(list)
subprocess.run(f'minimap2 -x ava-pb -t {thread} -o tmp.avsa.paf {infasta} {infasta}', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
with open(f'tmp.avsa.paf', 'r') as r:
    for line in r:
        if int(line.split()[10]) >= int(maxoverlapalignlen):
            overlapdict[line.split()[0]].append(line.split()[5])
            overlapdict[line.split()[5]].append(line.split()[0])

# Read Hi-C link density
print('Phasing Hi-C linkage...')
subprocess.run(f'samtools view -o tmp.mrq.sam -@ {thread} {inbam}', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
with open(f'tmp.mrq.sam', 'r') as r:
    signaldict = collections.defaultdict(int)
    for line in r:
        if not line.startswith('@'):
            qryid = line.split()[0]
            refid = line.split()[2]
            ref2id = line.split()[6]
            if ref2id != '=' and ref2id != '*':
                tmplist = [refid, ref2id]
                tmplist.sort()
                key = '-'.join(tmplist)
                signaldict[key] += 1
with open(f'tmp.prunedhiclink', 'w') as w:
    out = sorted(signaldict.items(), key=lambda x:x[1], reverse=True)
    for [link, value] in out:
        leftid, rightid = link.split('-')[0], link.split('-')[1]
        if leftid in overlapdict and rightid in overlapdict and leftid not in overlapdict[rightid] and rightid not in overlapdict[leftid]:
            w.write(f'{link}\t{value}\n')
with open(f'tmp.prunedhiclink', 'r') as hi:
    linkdict = collections.defaultdict(list)
    for line in hi:
        linkdict[line.split()[0].split('-')[0]].append([line.split()[0].split('-')[1], line.split()[1]])
        linkdict[line.split()[0].split('-')[1]].append([line.split()[0].split('-')[0], line.split()[1]])

# Import and filter init grouped contigs
print('Reading group...')
with open(allgroup, 'r') as ag:
    odict = {}
    for line in ag:
        if not line.startswith('#'):
            odict[line.split()[0]] = line.strip()
        
tigdict = {}
i = 1
for groupfile in ingroup:
    with open(groupfile, 'r') as g:
        for line in g:
            if not line.startswith('#'):
                tigid = line.split()[0]
                group = str(i)
                if len(fastadict[tigid]) >= int(mintiglen): # remove short
                    tigdict[tigid] = group
    i += 1
# add unclustered contigs
for tigid in fastadict:
    if tigid not in tigdict:
        tigdict[tigid] = 'X'

tiglistfromlongest = sorted(tigdict.keys(), key=lambda x: len(fastadict[x]), reverse=True)
tiglistfromshortest = tiglistfromlongest[::1]

# If two contigs in a group are allelic, remove the one with less Hi-C link density.
print('Removing overlap contig in group...')           
for tigid in tiglistfromshortest:
    if tigdict[tigid] != 'X':
        grouplist = []
        for otid in tigdict:
            if tigdict[otid] == tigdict[tigid]:
                grouplist.append(otid)
        for otid in grouplist:
            if tigid != otid:
                if otid in overlapdict[tigid]: # Large overlap detected, decide remove which one
                    tlink, olink = 0, 0
                    for [targetid, score] in linkdict[tigid]:
                        if targetid != otid:
                            tlink += int(score)
                    for [targetid, score] in linkdict[otid]:
                        if targetid != tigid:
                            olink += int(score)
                    if tlink >= olink:
                        print(f'{otid} (length {len(fastadict[otid])}, score {olink}) removed from group {tigdict[otid]} for conflict with {tigid} (length {len(fastadict[tigid])}, score {tlink}).')
                        tigdict[otid] = 'X'
                    else:
                        print(f'{tigid} (length {len(fastadict[tigid])}, score {tlink}) removed from group {tigdict[tigid]} for conflict with {otid} (length {len(fastadict[otid])}, score {olink}).')
                        tigdict[tigid] = 'X' 

# Assign unplaced contigs by Hi-C link density, while avoid cluster allelic contigs together.
print('Assigning unplaced contig...')               
for tigid in tiglistfromlongest:
    if tigdict[tigid] == 'X' and tigid in linkdict:
        sortedlinklist = sorted(linkdict[tigid], key=lambda x: int(x[1]), reverse=True)
        for [targetid, score] in sortedlinklist:
            flag = False
            if targetid in tigdict and tigdict[targetid] != 'X':
                grouplist = []
                for otid in tigdict:
                    if tigdict[otid] == tigdict[targetid]:
                        grouplist.append(otid)
                for otid in grouplist:
                    if tigid != otid:
                        if otid in overlapdict[tigid]: # Large overlap detected, try next group.
                            print(f'{tigid} (length {len(fastadict[tigid])}, score {score}) cannot join group {tigdict[targetid]} for conflict with {otid}.')
                            flag = True
                if flag == True:
                    continue
                tigdict[tigid] = tigdict[targetid]
                print(f'{tigid} (length {len(fastadict[tigid])}, score {score}) joins group {tigdict[tigid]} based on link with {targetid}.')
                break

print('Writing result...')
groupdict = collections.defaultdict(list)
for tigid, group in tigdict.items():
    groupdict[group].append(tigid)
for group in groupdict:
    if group != 'X':
        with open(f'group{group}.txt', 'w') as w:
            w.write('#Contig\tRECounts\tLength\n')
            for tigid in groupdict[group]:
                w.write(odict[tigid]+'\n')
                