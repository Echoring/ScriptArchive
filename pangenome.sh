##1.Genome Assembly, Annotation, and Evaluation

# Hifiasm assemble contigs
hifiasm -o sample -t 64 sample.fq 1>sample.hifiasm.log 2>sample.hifiasm.err

#contigs to chromosome
python quartet_assemblymapper.py -r Hongyang.hap1.fasta -q Acsample.fasta
python quartet_assemblymapper.py -r Midao.hap1.fasta -q Aesample.fasta

##BUSCO evaluation
chrfile=$1;
busco -i $chrfile -f -c 32 -o $2 -m geno -l embryophyta_odb10/

##LTR Assembly Index
ltr_finder  $1 >$1.scn
LTR_retriever -threads 32 -genome $1 -infinder $1.scn

#genome structure annotation
genome=$1
sample=${genome##*/}
mkdir $sample
braker	--useexisting	--genome=$genome	--prot_seq=Actinida.pep	--rnaseq_sets_ids=Acstorge	--rnaseq_sets_dir=rnaseq_dir	--species=Hongyang	--gff3	--workingdir=$sample

#gene function annotation
python emapper.py -i sample.pep -o output --itype proteins -m diamond --cpu 32
pfam_scan.pl -fasta sample.pep -dir pfamdatabase/ -outfile sample.pfam -as


##2. Pangenome family analysis

#orthogroups  identification
orthofinder -f data -S diamond -og -t 64

#saturation curve
Rscript CorePan_geneAnalysis.R Orthogroups_PAV.tsv 500 sim500
#Classification of gene families
Rscript Core.R Orthogroups_PAV.tsv


##3.NLRpan analysis

#NLRs identification
pfam_scan.pl -fasta sample.pep -dir ~/soft/pfamdatabase/ -outfile sample.pfam -as
grep "NB-ARC" $sample.pfam |cut -d " " -f1 |sort -u >$sample.nlrlist
java -Xmx80000M -jar NLR-Annotator-v2.1b.jar -i $sample.cds -g $sample.gff -x mot.txt  -y store.txt -t 32

#NLRpan construction
orthofinder -f 01.data -S diamond -og -t 32


##4.RNA-seq analysis
ref=$1
gtf=$2
r1=$3
r2=$4
s=${r1##*/}
sample=${s%%_*}
hisat2 --new-summary -p 64 -x $ref -1 $r1 -2 $r2 |samtools sort -@ 32 --output-fmt BAM -o $sample.sorted.bam
samtools flagstat $sample.sorted.bam >$sample.flagstat
Rscript run-featurecounts.R -b $sample.sorted.bam  -g   $gtf -o $sample


## 5. construction of reference-unbiased pangenome graph 

#run pggb
partition-before-pggb -i all44.chr.rename -o 02.partition -n 44 -t 64 -p 90 -s 50k  -V 'MD31#1:1000'
pggb -i 02.partition/all_4_haplotypes.pansn.fa.dac1d73.community.0.fa \
     -o 02.partition/all_4_haplotypes.pansn.fa.dac1d73.community.0.fa.out \
     -s 50000 -l 250000 -p 90 -c 1 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10M \
     -n 4 -j 0 -e 0 -G 700,1100 -P 1,4,6,2,26,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V MD31#1:1000 --threads 32 --poa-threads 64

#SVs stat and annotation
perl vcf_freq1.pl	all.vcf >all.vcf.freq
awk '{if($17=="SV" && $16<0.03) print $_}' all.vcf.freq >all.freq.filter.vcf
perl sv_gt_filter.pl all.freq.filter.vcf >all.freq.gtfiter.vcf
java -jar snpEff.jar  -c snpEff.config  -ud 2000 -o vcf Midaohap1 all.freq.gtfiter.vcf >all.freq.gtfiter.vcf.anno


