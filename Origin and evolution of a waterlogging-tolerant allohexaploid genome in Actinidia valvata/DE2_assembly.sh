#!/usr/bin/env bash
# Required packages: hifiasm bwa samtools ALLHiC gmap matlock python3 minimap2 juicer 3d-dna 
# Custom scripts: gmap2AlleleTable.py globle_rescue_alt.py rescue_alt.py
# Inital files: HiFi and Hi-C data (sample.hifi.fq.gz sample_hic_raw_1.fq.gz sample_hic_raw_2.fq.gz); Reference genome and annotation (reference.genome.fasta reference.gff3 reference.cds.fasta)

# Hifiasm assemble unitigs
hifiasm -o sample -t 64  --h1 sample_hic_raw_1.fq.gz --h2 sample_hic_raw_2.fq.gz sample.hifi.fq.gz
awk '/^S/{print ">"$2;print $3}' sample.hic.p_utg.gfa >sample.hic.p_utg.gfa.fasta

# Hi-C alignment
bwa index -a bwtsw sample.hic.p_utg.gfa.fasta
samtools faidx sample.hic.p_utg.gfa.fasta
bwa aln -t 40 sample.hic.p_utg.gfa.fasta sample_hic_raw_1.fq.gz > sample_R1.sai
bwa aln -t 40 sample.hic.p_utg.gfa.fasta sample_hic_raw_2.fq.gz > sample_R2.sai
bwa sampe sample.hic.p_utg.gfa.fasta sample_R1.sai sample_R2.sai sample_hic_raw_1.fq.gz sample_hic_raw_2.fq.gz > sample.bwa_aln.sam
ALLHiC/scripts/PreprocessSAMs.pl sample.bwa_aln.sam sample.hic.p_utg.gfa.fasta MBOI # ALLHiC scripts
ALLHiC/scripts/filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam # ALLHiC scripts
samtools view -bt sample.hic.p_utg.gfa.fasta.fai sample.clean.sam > sample.clean.bam
samtools sort -m 10G -o sample.clean.sorted.bam -@ 40 sample.clean.bam
matlock bam2 juicer sample.clean.sorted.bam merged_nodups.txt

# Chromosome partition
gmap_build -D . -d DB sample.hic.p_utg.gfa.fasta
gmap -D . -d DB -t 40 -f 2 -n 4 reference.cds.fasta > gmap.gff3
python3 gmap2AlleleTable.py reference.gff3 gmap.gff3 Allele.ctg.table # Alternative script of ALLHiC gmap2AlleleTable.pl, deal with unstandard reference gff3 format.
minimap2 -x asm5 -t 40 -o avsref.utg.paf reference.genome.fasta sample.hic.p_utg.gfa.fasta
python3 globle_rescue_alt.py Allele.ctg.table avsref.utg.paf Allele.ctg.new.txt # Alternative script of ALLHiC's globle rescue, resolve unitigs that do not contain genes annotated on reference
ALLHiC/scripts/partition_gmap.py -r sample.hic.p_utg.gfa.fasta -g Allele.ctg.new.txt -b sample.clean.sorted.bam -d chrpart -t 40 # ALLHiC scripts

# ALLHiC phasing
for i in {01..29}; do
    cd chrpart/Chr$i
    ALLHiC_prune -i ../../Allele.ctg.new.txt -b Chr$i.bam -r Chr$i.fa
    allhic extract prunning.bam Chr$i.fa --RE GATC
    allhic partition prunning.counts_GATC.txt prunning.pairs.txt 6 --minREs $minREs --maxLinkDensity 100 # $minREs: Chromosome 05:1000; 03,04,11,25:2500; 06:3000; remaining:1500
    python3 rescue_alt.py Chr$i.fa Chr$i.bam 50000 50000 40 prunning.counts_GATC.txt prunning.counts_GATC.*g*.txt # Alternative script for ALLHiC_rescue, resolve ALLHiC_rescue ignore allelic contig information
    allhic extract Chr$i.bam Chr$i.fa --RE GATC
    for j in `ls group*.txt`; do 
        allhic optimize $j Chr${i%%bam}clm
        done
    ALLHiC_build Chr$i.fa

    juicer/juicebox_scripts/juicebox_scripts/agp2assembly.py groups.agp groups.assembly # juicer script
    bash 3d-dna/visualize/run-assembly-visualizer.sh -q 1 -p true groups.assembly ../../merged_nodups.txt # 3d-dna script
    cd ../..
done

# Manual curation using juicebox software with groups.hic and groups.assembly file, generate groups.review.assembly

# Generate final assembly
for i in {01..29}; do
    python3 juicer/juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py -a chrpart/Chr$i/groups.review.assembly -f chrpart/Chr$i/Chr$i.fa -p chrpart/Chr$i/Chr$i.review -s # juicer script    
    awk '/^>/ { if (seq) {print seq} print; seq = ""; next; } { seq = seq $0 } END { if (seq) {print seq} }' chrpart/Chr$i/Chr$i.review.fasta > chrpart/Chr$i/Chr$i.nowrap.fasta
    grep -A 1 --no-group-separator '^>Chromosome[0-9]$' chrpart/Chr$i/Chr$i.nowrap.fasta | sed "s/^>Chromosome\([0-9]\)/>Avh\1Chr$i/" >> genome.fasta
    grep -A 1 --no-group-separator '^>utg' chrpart/Chr$i/Chr$i.nowrap.fasta >> unassigned_utg.fasta
    grep "Chromosome[0-9]" chrpart/Chr$i/Chr$i.review.agp | sed "s/^Chromosome\([0-9]\)/Avh\1Chr$i/" >> genome.agp
done

