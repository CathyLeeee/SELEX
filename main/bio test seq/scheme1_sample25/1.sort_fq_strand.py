#!usr/bin/python
"""function: pick up and sort out reads in fq file on strandness(+/-)
   input: fq1 and fq2, test sequecne with strandness
   output: fq1(+) and fq2(-)
"""

from string import atoi
def make_seq_strandness_dict(inpath):
    """
    function: make seq_strandness_dict according to sequecne with strandness file
    input: file path(QNAME in file coordinate with sam file)
    output: seq_strandness_dict
    #eg.seq_strandness_dict={"ATCGGGG":"-"}   #query seq : strandness
    attention: seq in this dict is for fq rather than sam!
    """
    global seq_strandness_dict
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    txt.pop(0)

    seq_strandness_dict_plus={x.split(",")[1].strip():x.split(",")[2].strip("\r\n")  for x in txt if x.split(",")[2].strip("\r\n")=="+"}
    seq_strandness_dict_minus={x.split(",")[1].strip():x.split(",")[2].strip("\r\n")  for x in txt if x.split(",")[2].strip("\r\n")=="-"}
    seq_strandness_dict=dict(seq_strandness_dict_plus,**seq_strandness_dict_minus)
    return seq_strandness_dict


def make_l_fq(inpath):
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    l_fq=[x.strip("\n") for x in txt]
    return l_fq

def write_new_fq(outpath,new_l_fq):
    f=open(outpath,"w")
    for i in new_l_fq:
        f.write(i)
        f.write("\n")
    f.close()
    

#main
seq_strandness_dict=make_seq_strandness_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test_sequence_of_sample25.csv")

l_fq1=make_l_fq(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.fastq")
l_fq2=make_l_fq(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.fastq")

new_l_fq1=[]
new_l_fq2=[]
for seq in seq_strandness_dict:
    print seq
    n0=len(new_l_fq1)
    if seq_strandness_dict[seq]=="+":    #strandness
        if  (seq in l_fq1):
            l_index=[i for i, x in enumerate(l_fq1) if x==seq]
            for index in l_index: 
                new_l_fq1.append("\n".join([l_fq1[index-1],l_fq1[index],l_fq1[index+1],l_fq1[index+2]]))
                new_l_fq2.append("\n".join([l_fq2[index-1],l_fq2[index],l_fq2[index+1],l_fq2[index+2]]))
        if  (seq in l_fq2):
            l_index=[i for i, x in enumerate(l_fq2) if x==seq]
            for index in l_index:
                new_l_fq1.append("\n".join([l_fq2[index-1],l_fq2[index],l_fq2[index+1],l_fq2[index+2]]))
                new_l_fq2.append("\n".join([l_fq1[index-1],l_fq1[index],l_fq1[index+1],l_fq1[index+2]]))
    elif seq_strandness_dict[seq]=="-":    #strandness
        if  (seq in l_fq2):
            l_index=[i for i, x in enumerate(l_fq2) if x==seq]
            for index in l_index:
                new_l_fq1.append("\n".join([l_fq1[index-1],l_fq1[index],l_fq1[index+1],l_fq1[index+2]]))
                new_l_fq2.append("\n".join([l_fq2[index-1],l_fq2[index],l_fq2[index+1],l_fq2[index+2]]))
        if  (seq in l_fq1):
            l_index=[i for i, x in enumerate(l_fq1) if x==seq]
            for index in l_index:
                new_l_fq1.append("\n".join([l_fq2[index-1],l_fq2[index],l_fq2[index+1],l_fq2[index+2]]))
                new_l_fq2.append("\n".join([l_fq1[index-1],l_fq1[index],l_fq1[index+1],l_fq1[index+2]]))
    n1=len(new_l_fq1)
    print n1-n0

print "".join(["len of new new_l_fq1:",str(len(new_l_fq1))])
print "".join(["len of new new_l_fq2:",str(len(new_l_fq2))])
            

write_new_fq(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.new.fastq",new_l_fq1)
write_new_fq(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.new.fastq",new_l_fq2)





    
