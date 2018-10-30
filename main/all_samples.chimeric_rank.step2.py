#!/sur/bin/python

#usage: "for number in []" gives the sample numbers; seperste files will be generated in seperate dir according to samples
"""
generate vector for R plot distribution according to 6 files generateed in step1
"""
from string import atof


nsam=[12433632,9754361,6717384,6771421,8311477,8118433]
fq=[]

n=0
for number in [16,20,21,22,23,25]:
    input="".join(["/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I",str(number),"/Sample_LHT-ZM-I",str(number),".chimeric_rank.txt"])

    f=open(input,"r")
    txt=f.readlines()
    f.close()
   
    txt.pop(len(txt)-1)
    txt.pop(len(txt)-1)
    txt.pop(len(txt)-1)
    print len(txt)
    rank=1
    for i in txt:
        i=i.strip("\n")
        fq.append("\t".join(["".join(["sample",str(number)]),str(rank),str(atof(i.split()[1])/atof(nsam[n]))]))
        rank=rank+1         
    n=n+1



f=open("/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/all_samples.chimeric_distribution_for_Rplot.txt","w")
for i in fq:
    f.write(i)
    f.write("\n")
f.close()
