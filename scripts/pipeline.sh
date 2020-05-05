#!/bin/bash
# Prealignment QC
#1.0 Perform fastqc on raw reads
fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/*.fastq.gz 
mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_untrimmed_reads
mv ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_untrimmed_reads/
#2.0 Trim raw reads
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  /home/ubuntu/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/*.fastq.gz \
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50
#3.0 Perform fastqc on trimmed data
fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
	/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P
mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads
mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/
#Alignment
#1.0 Create aligned_data directory, perform alignment and transfer aligned data to the directory
mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data
bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 \
 ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
 ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P > ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.sam
#2.0 Convert aligned_data to BAM, sort it and index it
#2.1 Convert sam to bam
samtools view -h -b ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.sam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.bam
#2.2 Sort bam
samtools sort ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.bam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted.bam
#2.3 Index bam
samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted.bam
#3.0 Post Alignment QC and Filtering
#3.1 Mark duplicates on sorted bam
cd ~/ngs_course/dnaseq_pipeline/data/aligned_data/
picard MarkDuplicates I=WES01_chr22m_sorted.bam O=WES01_chr22m_sorted_marked.bam M=marked_dup_metrics.txt
#3.2 Index the sorted bam
samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_marked.bam
#3.3 Filter the sorted dupicate marked bam
samtools view -F 1796  -q 20 -o ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam \
 ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_marked.bam
#3.4 Index the filtered sorted duplicate marked bam 
samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam
#Variant Calling
#1.0 Unzip the reference genome
zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
#2.0 Index the reference genome  
samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa
#3.0 Call variants using freebayes
freebayes --bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa \
 --vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf
#4.0 Compress the vcf file 
bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf
#5.0 Index the vcf file
tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz
#6.0 Filter the vcf
#6.1 Hard filter
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf
#6.2 Soft filter 
bedtools intersect -header -wa -a ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf -b ~/ngs_course/dnaseq_pipeline/data/chr22.genes.hg19.bed \
	> ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf
#6.3 Compress the  filtered vcf file
bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf
#6.4 Perform a table function of the  filtered vcf
tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz
#Annotate Variants
cd ~/annovar
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput
#1.0 Run a table function on the annovar output
cd ~/annovar
~/annovar/table_annovar.pl ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput humandb/ -buildver hg19 -out \
 ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro, -operation g,g,f,f,f -otherinfo -nastring . -csvout
