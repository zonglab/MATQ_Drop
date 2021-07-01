#!/bin/bash

# set max wallclock time
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
# set name of job
#PBS -N AXX

sample=AXX
seqdir=BXX
datadir=CXX/${sample}
scriptdir=DXX


echo ${sample}
sample_name=${sample//_/-}
mkdir -p $datadir && cd $datadir

cat $seqdir/${sample_name}_S*_R1_001.fastq.gz > ${sample}_R1.fastq.gz
cat $seqdir/${sample_name}_S*_R2_001.fastq.gz > ${sample}_R2.fastq.gz

##Poly A trimming
cutadapt -j 4 -a AAAAAAAAAAAAAAAAAAAA --minimum-length=30 --pair-filter=any -o ${sample}_R1_trimmed.fastq.gz -p ${sample}_R2_trimmed.fastq.gz ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz > trimming_report.txt
mv trimming_report.txt ${sample}_trimming_report.txt

##Barcode assignment
zcat ${sample}_R2_trimmed.fastq.gz | awk '(NR%4==2)' | cut -c 1-36 > barcode36.txt
python ${scriptdir}/MATQ_Drop_Barcode_count.py barcode36.txt ${scriptdir}/barcode144x144.dat > barcode.count

sort -k2 -rn barcode.count |awk -v OFS='\t' '{print $1, $2}' > ${sample}_barcode_read_count.tsv
sort barcode.count | awk -v OFS='\t' '{print $1, $3, $2}' > ${sample}_whitelist_all.txt

rm barcode36.txt
rm barcode.count


umi_tools extract --stdin ${sample}_R2_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{10})(?P<discard_1>TGGTAGGTGGTAGAGA){s<=16}(?P<cell_2>.{10})(?P<umi_1>.{6}).*"  --stdout ${sample}_R2_trimmed_extracted.fastq.gz --read2-in ${sample}_R1_trimmed.fastq.gz --read2-out=${sample}_R1_trimmed_extracted.fastq.gz --whitelist=${sample}_whitelist_all.txt --filter-cell-barcode --error-correct-cell

##STAR Mapping
STAR --runThreadN 4  --runMode alignReads --genomeDir ~/hg19/starindex --readFilesCommand zcat --readFilesIn ${sample}_R1_trimmed_extracted.fastq.gz  --outSAMtype BAM SortedByCoordinate

mv Aligned.sortedByCoord.out.bam  Aligned.r1.bam
samtools view -b -q 250 Aligned.r1.bam | samtools sort - > r1.bam
mv Log.final.out Log.r1.final.out


##Transcript-based gene expression matrix
featureCounts -t transcript --extraAttributes gene_name -a ~/hg19/gencode/gencode.v19.annotation.gtf  -o gene_assigned.hg -R BAM r1.bam -T 4 -s 2
samtools sort r1.bam.featureCounts.bam -o assigned.sorted.bam
rm r1.bam.featureCounts.bam
samtools index assigned.sorted.bam
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I assigned.sorted.bam -S counts.tsv
mv counts.tsv ${sample}.counts.tsv

python3 ${scriptdir}/UMI_Gene_count.py ${sample}.counts.tsv ${sample}_barcode_read_count.tsv ${sample}_UMI_Gene_count.tsv


sed -i '1d' ${sample}_UMI_Gene_count.tsv
sort -k3 -nr ${sample}_UMI_Gene_count.tsv > ${sample}_UMI_Gene_count_sorted.tsv


##Exon based gene expression matrix
samtools view -H assigned.sorted.bam > SAM_header
samtools view assigned.sorted.bam | grep XT > assigned.sam.body
cat SAM_header assigned.sam.body | sed 's/\tXT\:Z\:[^\t]*//' > assigned.sam
samtools view -b assigned.sam > transcript.assigned.sorted.bam
rm assigned.sam
rm assigned.sam.body
rm SAM_header

featureCounts -t exon --extraAttributes gene_name -a ~/hg19/gencode/gencode.v19.annotation.gtf  -o gene_assigned.hg -R BAM transcript.assigned.sorted.bam -T 4 -s 2
samtools sort transcript.assigned.sorted.bam.featureCounts.bam -o exon_assigned.sorted.bam
rm transcript.assigned.sorted.bam.featureCounts.bam
samtools index exon_assigned.sorted.bam
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I exon_assigned.sorted.bam -S counts.tsv
mv counts.tsv ${sample}_exon.counts.tsv


