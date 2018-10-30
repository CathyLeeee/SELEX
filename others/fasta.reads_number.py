#!usr/bin/python


##find reads in genome




import gzip  
  

def reads_number(fasta_path):
    f=gzip.open(fasta_path)
    txt=f.readlines()
    f.close()
    
    count=0
    for i in txt:
        if "@" in i:
            count=count+1
    return(count)



f=open(r"/Share/home/zhangqw/SELECT/20160803/LHT-ZM_HVNKNCCXX_L4/fasta.reads_number.txt","w")
f.write("\t".join(["reads in fasta1",str(reads_number("/Share/home/zhangqw/SELECT/20160803/LHT-ZM_HVNKNCCXX_L4/LHT-ZM_HVNKNCCXX_L4_1.clean.fq.gz")),"\n"]))
f.write("\t".join(["reads in fasta1",str(reads_number("/Share/home/zhangqw/SELECT/20160803/LHT-ZM_HVNKNCCXX_L4/LHT-ZM_HVNKNCCXX_L4_2.clean.fq.gz")),"\n"]))
f.close()




            



