#$1=fastq1 $2=fastq2 $3=outbam $4=outcount
reference=$1
fastq1=$2
fastq2=$3
outprefix=$4
if [ ! -f "$reference.1.ht2" ]; then
    hisat2-build $reference $reference -p 40
fi
hisat2 -p 40 -x $reference -1 $fastq1 -2 $fastq2 | samtools sort -o $outprefix.bam --threads 40 -O BAM
featureCounts -a reference.gtf -o $outprefix.count -M --fraction -T 40 -p $outprefix.bam
