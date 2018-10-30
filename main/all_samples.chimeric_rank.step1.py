#!/sur/bin/python

#usage: "for number in []" gives the sample numbers; seperste files will be generated in seperate dir according to samples

""" 
 this script generate the top 1000 chimeric reads and fq of 6 samples 
 output will be represented in seperated files      
"""
from collections import Counter
from string import atoi
from string import atof


for number in [16,20,21,22,23,25]:
    fq1="".join(["/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I",str(number),"/LHT_ZM_I",str(number),"_R1.clean_pe.fastq"])
    fq2="".join(["/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I",str(number),"/LHT_ZM_I",str(number),"_R2.clean_pe.fastq"])
    sam="".join(["/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I",str(number),"/LHT_ZM_I",str(number),".sam"])
    output="".join(["/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I",str(number),"/Sample_LHT-ZM-I",str(number),".chimeric_rank.txt"])

#readin fastq and sam
#fastq1
    f=open(fq1,"r")
    file1=f.readlines()
    f.close()
    reads1=file1[1::4]

#fastq2
    f=open(fq2,"r")
    file2=f.readlines()
    f.close()
    reads2=file2[1::4]

#sam
    f=open(sam,"r")
    tsam=f.readlines()
    f.close()

##pre-process fq and sam

    for i in range(len(reads1)):
        reads1[i]=reads1[i].strip("\n")
        reads2[i]=reads2[i].strip("\n")

#frequency ; top 50 reads
    reads=reads2 #all reads array
    nreads=len(reads)#total reads number
    sorted_reads=Counter(reads).most_common()   #ranked reads list(mixture)
    sorted_reads=sorted_reads[1:1000]

#build pos dict according to sam
    chimeric=[]  #sequence
    for i in tsam:
        if not( "@" in i):
            i=i.strip("\n")
            temp=i.split()
            key=temp[9]    #sequence
            if temp[6]!="=":
                chimeric.append(key)
            
        
#ouput
    f=open(output,"w")

    s_chimeric=0
    s_total=0
    for i in sorted_reads:
        seq=i[0].strip("\n")
        freq=str(i[1])
        s_total=s_total+i[1]
        if (seq in chimeric):
            s_chimeric=s_chimeric+i[1]
            f.write("\t".join([seq,freq]))   #sequence, frequency, chr, pos, strand
            f.write("\n")
    f.write("\t".join(["total chimeric reads number:",str(s_chimeric)]))
    f.write("\n")
    f.write("\t".join(["total reads number of top 1000:",str(s_total)]))
    f.write("\n")
    f.write("\t".join(["percentage of chimeric in total:",str(atof(s_chimeric)/atof(s_total))]))
    f.write("\n")
    f.close()
