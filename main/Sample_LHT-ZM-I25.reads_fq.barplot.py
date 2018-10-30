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

#reads frequency 
reads=reads1+reads2 #all reads array
sorted_reads=Counter(reads).most_common()   #ranked reads list(mixture)
readsf=[]
for i in sorted_reads:
    freq=i[1]    #int
    readsf.append(freq)

##count
nsub=1000  #n number of subsections for barplot
subsection=[0,]*nsub
for i in readsf:
    for j in range(nsub):
        if i//2500==j:     #length of subsections
            subsection[j]=subsection[j]+1
xlab=range(2500,2502500,2500)  #len(xlab)==nsub

#ouput
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_fq.barplot.txt","w")
for i in range(nsub):
    col1=str(xlab[i])
    col2=str(subsection[i])
    f.write("\t".join([col1,col2]))   #sequence, frequency, chr, pos
    f.write("\n")
f.close()
