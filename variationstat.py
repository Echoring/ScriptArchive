import sys
import collections
inalignedfasta = sys.argv[1]
refseqid = sys.argv[2]

def readFastaAsDict(fastafile):
    fastaDict = {}
    with open(fastafile, 'r') as f:
        allline = f.read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

fastadict = readFastaAsDict(inalignedfasta)
refseq = fastadict[refseqid]
matchdict = collections.defaultdict(int)
snpdict = collections.defaultdict(int)
insertdict = collections.defaultdict(int)
deletedict = collections.defaultdict(int)
lendict = collections.defaultdict(int)
for sid, seq in fastadict.items():
    for i in range(len(seq)):
        if refseq[i] == seq[i] and refseq[i] != '-' and seq[i] != '-':
            matchdict[sid] += 1
        elif refseq[i] != seq[i] and refseq[i] != '-' and seq[i] != '-':
            snpdict[sid] += 1
        elif refseq[i] != seq[i] and refseq[i] == '-' and seq[i] != '-':
            insertdict[sid] += 1
        elif refseq[i] != seq[i] and refseq[i] != '-' and seq[i] == '-':
            deletedict[sid] += 1
        if seq[i] != '-':
            lendict[sid] += 1

for sid, count in matchdict.items():
    print(f'Match: {sid} {count}')
for sid, count in snpdict.items():
    print(f'SNP: {sid} {count}')
for sid, count in insertdict.items():
    print(f'Insert: {sid} {count}')
for sid, count in deletedict.items():
    print(f'Delete: {sid} {count}')
for sid, count in lendict.items():
    print(f'Length: {sid} {count}'  )