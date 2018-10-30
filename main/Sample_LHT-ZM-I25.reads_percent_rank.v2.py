#!usr/bin/python
###adjusted frequency of sample25 with sample16 as background
#adjuxt every seq in sample 25 then rank 50
##locate reads in genome
from collections import Counter

#sample 16 as background
#fastq1
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I16/LHT_ZM_I16_R1.clean_pe.fastq","r")
file1=f.readlines()
f.close()
reads1=file1[1::4]

#fastq2
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I16/LHT_ZM_I16_R2.clean_pe.fastq","r")
file2=f.readlines()
f.close()
reads2=file2[1::4]
for i in range(len(reads1)):
    reads1[i]=reads1[i].strip("\n")
    reads2[i]=reads2[i].strip("\n")

#frequency ; top 50 reads
reads=reads1+reads2 #all reads array
sorted_reads=Counter(reads).most_common()   #ranked reads list(mixture)
bg_fq_dict={}   #background frequency dict
for i in sorted_reads:  #for the first 50th
    seq=i[0].strip("\n")
    freq=i[1]   #type:int
    bg_fq_dict[seq]=freq

####sample25

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

#sam
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.sam","r")
tsam=f.readlines()
f.close()

##pre-process fq and sam

for i in range(len(reads1)):
    reads1[i]=reads1[i].strip("\n")
    reads2[i]=reads2[i].strip("\n")

#frequency ; top 50 reads
reads=reads1+reads2 #all reads array
nreads=len(reads)#total reads number
sorted_reads=Counter(reads).most_common()   #ranked reads list(mixture)
sorted_reads=sorted_reads[:1000]
##adjust frquecny according to sample 25
m=0
for i in sorted_reads:  #for the first 50th
    n=sorted_reads.index(i)
    seq=i[0].strip("\n")
    if seq in bg_fq_dict:
        freq=i[1]-bg_fq_dict[seq]   #adjust
        sorted_reads[n]=[seq,freq]
    else:
        print seq
        print "\n"
        m=m+1
print m
sorted_reads=sorted(sorted_reads,reverse=True,key=lambda x:x[1])# high to low according to the 2nd element
sorted_reads=sorted_reads[:50]  #list

#build pos dict according to sam
pos_dict={}
for i in tsam:
    if not( "@" in i):
        i=i.strip("\n")
        temp=i.split()
        key=temp[9]    #sequence
        if temp[6]=="=":
            temp[6]=temp[2]
       # if (not(key in pos_dict)) & (len(key)==150):
        value1=",".join([temp[2],temp[3]]) #chr,pos of current reads
        value2=",".join([temp[6],temp[7]]) #chr,pos of next reads
        value=";".join([value1,value2])
        pos_dict[key]=value
#if "ACAGAGCCTCTCAGGCTCTAATAAGAGGAGTGGGCAAGGGTGGAATGTTTGAGACTGAGACCTCCAAGGACTGGCCAGCTTTGCAGCCCCACACCTCTGTGCTGCTCAGCCCCTGCACCAAGGAGGCTTCCTGGGTGACCAGGGCAGGAG" in pos_dict:
#    print pos_dict["ACAGAGCCTCTCAGGCTCTAATAAGAGGAGTGGGCAAGGGTGGAATGTTTGAGACTGAGACCTCCAAGGACTGGCCAGCTTTGCAGCCCCACACCTCTGTGCTGCTCAGCCCCTGCACCAAGGAGGCTTCCTGGGTGACCAGGGCAGGAG"]

##locate reads in genome
from string import atoi
#make annotation dict
f=open(r"/Share/home/zhangqw/SELECT/ref/HG19-RefSeq.txt","r")
csv=f.readlines()
f.close()
for i in csv:
    n=csv.index(i)
    temp=i.strip("\n")
    templ=temp.split("\t")
    csv[n]=";".join([templ[12],templ[1],templ[2],templ[6],templ[7],templ[4],templ[5],templ[9],templ[10]])
               # gene ;gene id; chr;cds1;cds2;tx1;tx2;exon1;exon2
#function #input chr_pos and csv
def annote(chr_pos,csv):

    chr_pos=chr_pos.split(";")[0]
    cpl=chr_pos.split(",")
    chr=cpl[0]
    pos=atoi(cpl[1])

    for i in csv:
        temp=i.split(";")
        if chr==temp[2] and pos>=atoi(temp[3]) and pos<=atoi(temp[4]): #chr &pos
            return(temp[0])
           # return("coding region")
    return("non-coding region")

#add keys to pos_dict
for i in sorted_reads:  #for the first 50th
    seq=i[0].strip("\n")
    freq=str(i[1])
    if (not(seq in pos_dict)):
        if seq in reads1:                        #sequence in reads1    
            n=reads1.index(seq)                  #index of the current sequence
            j=reads2[n].strip("\n")                       # j is the paired sequence of current sequence
        else: #seq in reads2:
            n=reads2.index(seq)
            j=reads1[n].strip("\n")
        if j in pos_dict:
            temp=pos_dict[j].split(";")
            temp.reverse()
            pos_dict[seq]=";".join(temp)
        else:
            pos_dict[seq]="NA"

#ouput
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.reads_percent_rank.v2.txt","w")
for i in sorted_reads:
    seq=i[0].strip("\n")
    freq=str(i[1])
    if (seq in pos_dict) and (pos_dict[seq]!="NA"):
        chr_pos=pos_dict[seq]
        annotation=annote(chr_pos,csv)
    else:
        chr_pos="NA"
        annotation="NA"
#    f.write("\t".join([seq,freq,chr_pos]))   #sequence, frequency, chr, pos
    f.write("\t".join([seq,freq,chr_pos,annotation]))   #sequence, frequency, chr, pos
    f.write("\n")
f.write("\t".join(["total reads number:",str(nreads)]))
f.write("\n")
f.close()
