#!usr/bin/python

##find reads in genome
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.fastq","r")
txt1=f.readlines()
f.close()

f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.fastq","r")
txt2=f.readlines()
f.close()

f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/fbline.txt","w")
for i in txt1:
    temp=i.strip("\n")
    if temp=="TGCTACTGGGAGGGCTAAGCACAGTGCTCACAGGCCCTGTGTGGTGAGGGCCGGCTGCCAGGGAACAGGAGGGGCTTCCCTTTGTCCTCTCCTGCCCTGGTCACCCAGCAAGCCTCCTGCCCTGGTCACCCAGGAAGCCTCCTTGGTGCA":
        n=txt1.index(i)
        f.write("fq1:\n")
        f.write(txt1[n-1])
        f.write(txt1[n])
        f.write("fq2:\n")
        f.write(txt2[n-1])
        f.write(txt2[n])

for i in txt2:
    temp=i.strip("\n")
    if temp=="TGCTACTGGGAGGGCTAAGCACAGTGCTCACAGGCCCTGTGTGGTGAGGGCCGGCTGCCAGGGAACAGGAGGGGCTTCCCTTTGTCCTCTCCTGCCCTGGTCACCCAGCAAGCCTCCTGCCCTGGTCACCCAGGAAGCCTCCTTGGTGCA":
        n=txt2.index(i)
        f.write("fq2:\n")
        f.write(txt2[n-1])
        f.write(txt2[n])
        f.write("fq1:\n")
        f.write(txt1[n-1])
        f.write(txt1[n])
f.close()





