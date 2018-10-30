#!usr/bin/python


##locate reads in genome
from collections import Counter

#readin fastq and sam
#fastq1
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.fastq","r")
file1=f.readlines()
f.close()
reads1=file1[1::4]

#fastq2
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.fastq","r")
file2=f.readlines()
f.close()
reads2=file2[1::4]

#pre-process fq 

for i in range(len(reads1)):
    reads1[i]=reads1[i].strip("\n")
    reads2[i]=reads2[i].strip("\n")

#frequency ; top 50 reads
reads=reads1+reads2 #all reads array
sorted_reads=Counter(reads).most_common()   #ranked reads list(mixture)


#ouput
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_fq.txt","w")
for i in sorted_reads:
    seq=i[0].strip("\n")
    freq=str(i[1])
    f.write("\t".join([seq,freq]))   #sequence, frequency, chr, pos
    f.write("\n")
f.close()
