#! /bin/bash  -v
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=40G

date +%Y%m%d%H%M

source $HOME/.bash_profile


##step1 filtering
#cutadapti

############/Share/home/tiangeng/software/cap-miRNA/bin/cutadapt -m 20 -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGAATACATCTCGTATGCCGTCTTCTGCTTG -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R1.clean.fastq.gz  -p /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R2.clean.fastq.gz    /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I15_R1.fastq.gz   /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I15_R2.fastq.gz

#trimommatic QC
#java   -jar /Share/home/tiangeng/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE  /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R1.clean.fastq.gz /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R2.clean.fastq.gz /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R1.clean.paired.fastq.gz   /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R1.clean.unpaired.fastq.gz /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R2.clean.paired.fastq.gz  /Share/home/zhangqw/SELECT/20161220/clean_data/LHT_ZM_I15_R2.clean.unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

#step2
#bwa
#/Share/home/zhangqw/local/bwa-0.7.12/bwa mem -t 48 -M /Share/home/zhangqw/SELECT/ref/hg19.fa /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.fastq.gz   /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.fastq.gz   >/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.sam

#bam to bed
#/Share/home/zhangqw/local/samtools-1.3.1/samtools  view -b -h  -S /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.sam >/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.bam

#version1
#/Share/home/zhangqw//local/bedtools2/bin/bedtools bamtobed  -i /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.bam > /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.bed

#/Share/home/zhangqw//local/bedtools2/bin/bedtools genomecov -i /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.bed -g /Share/home/zhangqw/SELECT/ref/hg19.genomesize.txt >/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.genomecov.txt

#version2
#sort -k1,1 -k2,2n  /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.bed  > /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.sorted.bed

#/Share/home/zhangqw//local/bedtools2/bin/bedtools genomecov -i /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.sorted.bed -g /Share/home/zhangqw/SELECT/ref/hg19.genomesize.txt -bga  >/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.genomecov.txt

#reads percentage rank
#python  /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_percent_rank.py

#reads
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT-ZM-I25.rmdup.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_fq.barplot.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_percent_rank.v2.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_fq.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/fbline.py
#R CMD BATCH --args /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/histo.ggplot2.qsub.R
#Rscript  /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/histo.ggplot2.qsub.R

###new round with fq2 the main file
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT-ZM-I25.mappping_ratio.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT-ZM-I25.filtering_sam.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.fq2_rank50.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/all_samples.chimeric_rank.step1.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step2_fq2.merge_reads.pos.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step3_base_depth_for_each_region.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test-test.py

###std1
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT-ZM-I25.filtering_sam.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.fq2_rank50.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step1_contig_region.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step2_fq2.merge_reads.pos.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step3_base_depth_for_each_region.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step3_base_depth_for_each_351_region.py


###std2 
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT-ZM-I25.filtering_sam.std2.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.fq2_rank50.std2.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step1_contig_region.std2.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step2_fq2.merge_reads.pos.std2.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step3_base_depth_for_each_region.std2.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step3_base_depth_for_each_351_region.std2.py

#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test-big_read_freq.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test-fbline.py

##merge reads for biochemical test
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/1.sort_fq_strand.py

#/Share/home/zhangqw/local/bwa-0.7.12/bwa mem -t 48 -M /Share/home/zhangqw/SELECT/ref/hg19.fa /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.new.fastq   /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.new.fastq   >/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.new.sam
#echo "filtering clip with updated new.sam"
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/2.filtering_sam.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/3.merge_reads.pos_and_strand.py
#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/4.final_merge.py

#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/3.1.rank50.py
echo "python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/get_seq_for_bio_test.py"
python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/get_seq_for_bio_test.py

#python /Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test-RNAME_in_sam_and_fq.py




date +%Y%m%d%H%M
