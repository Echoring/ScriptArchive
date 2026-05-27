#!/usr/bin/env bash
# $1=Hi-C_R1.fq $2=Hi-C_R2.fq $3=ploidy $4=prefix $5=thread $6=telomere $7+=HiFi.fq
# mamba activate hifiasm
set -euxo pipefail

if [ $3 == 2 ]; then
    nhap="--dual-scaf"
else
    nhap="--n-hap $3"
fi
hifiasm --h1 "$1" --h2 "$2" $nhap -o $4 -t $5 --telo-m $6 ${@:7}
python3 ~/script/gfa2fasta.py *tg.gfa
seqkit stat -a -b -N 50,90 -o stat.txt -j 40 --quiet *tg.fasta