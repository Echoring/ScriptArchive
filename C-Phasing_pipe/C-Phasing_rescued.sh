#!/usr/bin/env bash
# $1=tig.fasta $2=HiC.R1.fastq/"-pqs" $3=HiC.R2.fastq/pairs.pqs $4=haplotype chr number $5=ploidy $6=RE site $7=thread $8=prefix $9=reference.fasta $10=noseq.gfa
# mamba activate cphasing
set -euxo pipefail

if [ $2 == "-pqs" ]; then
    cphasing pipeline -f $1 --pairs $3 -n ${4}:${5} -hcr -p $6 -t $7 --chimeric-correct -ss 4,5
else
    cphasing pipeline -f $1 -hic1 $2 -hic2 $3 -n ${4}:${5} -hcr -p $6 -t $7 --chimeric-correct -ss 4,5
fi
cphasing collapse from-gfa ${10}
hg=$(ls cphasing_output/3.hyperpartition/*.corrected.pairs.q1.e5m.hg 2>/dev/null || ls cphasing_output/3.hyperpartition/*.pairs.q1.e5m.hg)
size=$(ls cphasing_output/*.corrected.contigsizes 2>/dev/null || ls cphasing_output/*.contigsizes)
allele=$(ls cphasing_output/3.hyperpartition/*.corrected.allele.table 2>/dev/null || ls cphasing_output/3.hyperpartition/*.allele.table)
cphasing collapse rescue $hg $size cphasing_output/3.hyperpartition/output.clusters.txt contigs.collapsed.contig.list -n $5 -at $allele
re=$(ls cphasing_output/2.prepare/*.corrected.counts_$6.txt 2>/dev/null || ls cphasing_output/2.prepare/*.counts_$6.txt)
clm=$(ls cphasing_output/2.prepare/*.corrected.clm.gz 2>/dev/null || ls cphasing_output/2.prepare/*.clm.gz)
contacts=$(ls cphasing_output/2.prepare/*.corrected.split.contacts 2>/dev/null || ls cphasing_output/2.prepare/*.split.contacts)
cphasing scaffolding collapsed.rescue.clusters.txt $re $clm -at $allele -sc $contacts -f $1 -t $7 -o $8.rescue.agp -m precision --corrected
cphasing agp2fasta $8.rescue.agp $1 --contigs -o $8.dup.tig.fasta
ragp=$(ls $8.rescue.corrected.agp 2>/dev/null || ls $8.rescue.agp)
cphasing collapse agp-dup $ragp -o $8.dup.agp
pairs=$(find cphasing_output -maxdepth 1 -name "*.corrected.pairs.pqs" 2>/dev/null | head -n 1)
if [ -z "$pairs" ]; then
    pairs=$(find cphasing_output -maxdepth 1 -name "*.pairs.pqs" 2>/dev/null | head -n 1)
fi
cphasing collapse pairs-dup $pairs collapsed.rescue.contigs.list -o $8.dup.pairs.pqs
fasta=$8.dup.tig.fasta
agp=$8.dup.agp
cphasing rename -r $9 -f $fasta -a $agp -t $7 --unphased -o $8.rename.agp
cphasing pairs2cool $8.dup.pairs.pqs $8.dup.pairs.pqs/_contigsizes $8.10k.cool --threads $7
cphasing plot -a $8.rename.agp -m $8.10k.cool -o $8.500k.HiC.png
cphasing pairs2mnd $8.dup.pairs.pqs -o $8.mnd.txt
cphasing utils agp2assembly $8.rename.agp > $8.assembly
bash ~/tools/3ddna/3d-dna/visualize/run-assembly-visualizer.sh $8.assembly $8.mnd.txt

# require quartet and rearrange_assembly
quartet.py am -r $9 -q $fasta -c 20000 -l 5000 -t $7 -a unimap --nofilter --keep --noplot -p $8
cphasing utils agp2assembly $8.draftgenome.agp > $8.quarTeT.assembly
python3 ~/script/rearrange_assembly.py $8.assembly $8.quarTeT.assembly > $8.preview.assembly
zip -j $8.juicebox.zip $8.preview.assembly $8.assembly $8.hic $8.500k.HiC.png