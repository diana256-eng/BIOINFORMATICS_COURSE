#!/bin/bash
# Prealignment QC
#Rename the fastq files
mv ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.qz ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.gz
mv ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.qz ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.gz
#1.0 Perform fastqc on raw reads
fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/*.fastq.gz
mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_untrimmed_reads
mv ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_untrimmed_reads/
#2.0 Trim raw reads
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.gz ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50
mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P.fastq
mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P.fastq
#3.0 Perform fastqc on trimmed data
fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P.fastq \
        /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P.fastq
mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads 
mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/
#Alignment
#1.0 Create aligned_data directory, perform alignment and transfer aligned data to the directory
mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:Nextera-NGS0001-blood\tDT:2017-02-23\tPU:11V6WR1' -I 250,50  ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P.fastq ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P.fastq > ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.sam  
#2.0 Convert aligned_data to BAM, sort it and index it
#2.1 Convert sam to bam
samtools view -h -b ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.sam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.bam
#2.2 Sort bam
samtools sort ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001.bam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam
#2.3 Index bam
samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam
#3.0 Post Alignment QC and Filtering
#3.1 Mark duplicates on sorted bam
cd ~/ngs_course/dnaseq_pipeline/data/aligned_data/
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
#3.2 Index the sorted bam
samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam
#3.3 Filter the sorted dupicate marked bam
samtools view -F 1796  -q 20 -o ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam \
 ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam 
#3.4 Index the filtered sorted duplicate marked bam 
samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam
#3.5 Collect insert size metrics  with PICARD
cd ~/ngs_course/dnaseq_pipeline/data/aligned_data/
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=NGS0001_insert_size_metrics.txt H=NGS0001_insert_size_histogram.pdf
mv ~/ngs_course/dnaseq_pipeline/data/aligned_data/*.txt ~/ngs_course/dnaseq_pipeline/data/logs/NGS0001_insert_size_metrics.txt
mv ~/ngs_course/dnaseq_pipeline/data/aligned_data/*.pdf ~/ngs_course/dnaseq_pipeline/data/logs/NGS0001_insert_size_histogram.pdf
#Variant Calling
#1.0 Unpack and index the reference genome
zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa   
samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
#2.0 Call variants using freebayes
freebayes --bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam  --fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa \
 --vcf ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf
#3.0 Compress the vcf file 
bgzip ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf
#4.0 Index the vcf file
tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf.gz
#5.0 Filter the vcf
#5.1 Hard filter
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
        ~/ngs_course/dnaseq_pipeline/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered.vcf
#5.2 Soft filter the variants
bedtools intersect -header -wa -a ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered.vcf -b ~/ngs_course/dnaseq_pipeline/data/annotation.bed \
        > ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf

#5.3 Compress the  filtered vcf file
bgzip ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf
#5.4 Index the compressed vcf file
tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf.gz
#Annotate Variants
cd ~/annovar
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered_annotation.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered.avinput
#1.0 Run a table function on the annovar output
cd ~/annovar
~/annovar/table_annovar.pl ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered.avinput humandb/ -buildver hg19 -out \
 ~/ngs_course/dnaseq_pipeline/results/NGS0001_filtered -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro, -operation g,g,f,f,f -otherinfo -nastring . -csvout



