#!usr/bin/python
#coding=utf-8
"annotate according to chr_pos;def function as annote with input (chr_pos,csv)"

f=open(r"/Share/home/zhangqw/SELECT/ref/HG19-RefSeq.test.txt","r")
csv=f.readlines()
f.close()
for i in csv:
    n=csv.index(i)
    temp=i.strip("\n")
    templ=temp.split("\t")
    csv[n]=";".join([templ[12],templ[1],templ[2],templ[6],templ[7],templ[4],templ[5],templ[9],templ[10]])

##find reads in genome
from string import atoi
#make annotation dict
def annote(chr_pos,csv):

    chr_pos=chr_pos.split(";")[0]
    cpl=chr_pos.split(",")
    chr=cpl[0]
    pos=atoi(cpl[1])

    for i in csv:
        temp=i.split(";")
        if chr==temp[2] and pos>=atoi(temp[3]) and pos<=atoi(temp[4]): #chr &pos
            return("coding region")
    return("non-coding region")

g="global variant"


if __name=="__main__":
    chr_pos="chr1,67100000;=,111798798"
    print annote(chr_pos,csv)

