#!usr/bin/python


##find reads in genome
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.sam","r")
txt=f.readlines()
f.close()

def unique(txt): #reads name only cccurred once
    reads=[]
    for i in txt:
        if not ("@" in i):
            reads.append(i.split()[0])

    reads=list(set(list(reads)))
    return "\t".join(["unique reads number",str(len(reads))])

def rmh(txt):   #remove header
    reads=[]
    for i in txt:
        if not ("@" in i):
            reads.append(i.split()[0])

    return "\t".join(["reads number",str(len(reads))])


if __name__=="__main__":
    print unique(txt)
    print rmh(txt)
