#!usr/bin/python

from string import atoi
def make_seq_strandness_dict(inpath):
    """
    function: make seq_strandness_dict according to sequecne with strandness file
    input: file path(QNAME in file coordinate with sam file)
    output: seq_strandness_dict
#eg.seq_strandness_dict={"ATCGGGG":"-"}   #for "+", the seq; for "-",the rc seq
    """
    global seq_strandness_dict
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    txt.pop(0)

    seq_strandness_dict={x.split(",")[1].strip():x.split(",")[2].strip("\r\n")  for x in txt}
    return seq_strandness_dict

def big_read(pe):
    """
    function: merge reads1+reads2 to big_reads
    input: "chr2,240026885;chr2,240026798"  (without order)
    ouput: "chr2,240026798,240027034"
    """
    global pos
    if pe!="NA":
        pos=pe.split(";")
    chr=pos[0].split(",")[0]
    pos1=min(atoi(pos[0].split(",")[1]),atoi(pos[1].split(",")[1]))
    pos2=max(atoi(pos[0].split(",")[1]),atoi(pos[1].split(",")[1]))
    return ",".join([chr,str(pos1),str(pos2+149)])


def rc(seq):
    """
    function: give the reverse complementary seq

    """
    rc_seq=[]
    for i in seq:
        if i=="A" or i=="a":
            rc_seq.append("T")
        elif i=="T" or i=="t":
            rc_seq.append("A")
        elif i=="C" or i=="c":
            rc_seq.append("G")
        elif i=="G" or i=="g":
            rc_seq.append("C")
        elif i=="N":
            rc_seq.append("N")
        else:
            rc_seq.append("*")
    rc_seq.reverse()
    return "".join(rc_seq)

def make_dict_hg19(inpath):
    """
    function: return dict_hg19
    input: path of ref
    output: dict_hg19={"chr1":"GGGGCCCC"}
    """
    f=open(inpath,"r")
    hg19=f.readlines()
    f.close()
    index=[hg19.index(i) for i in hg19 if ">" in i]
    index.append(len(hg19))

    dict_hg19={hg19[index[i]].strip("\n").lstrip(">"):"".join([x.strip("\n") for x in hg19[(index[i]+1):(index[i+1])]]) for i in range(len(index)-1)}
    return dict_hg19

def get_seq(dict_hg19,big_read,strandness):
    """
    function: get base sequece in ref according to big_read pos and strandness
    input: (dict_hg19,"chr2,240026798,240027034","+")
    output: "GGGGTTACCC"
    """
    chr=big_read.split(",")[0]
    start=atoi(big_read.split(",")[1])
    end=atoi(big_read.split(",")[2])
    if chr in dict_hg19:
        if strandness=="+":
            return dict_hg19[chr][start-1:end]
        elif strandness=="-":
            return rc(dict_hg19[chr][start-1:end])    
    else:
        return "NA"
    
def output(strand,inpath):
    """
    function: give output accoding to strandness
    input:"+"/"-"
    output:["seq,pos,strand"...]
    """
    global l_seq
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    
    txt.pop(len(txt)-1)
    
    txt=[line.strip("\n")  for line in txt  if  len(line.strip("\n").split("\t")[0])<500 and line.strip("\n").split("\t")[3]==strand]
    output=[]
    l_seq=[]
    for temp  in txt:
        seq=temp.split("\t")[0]
        l_seq.append(seq)
        pos=big_read(temp.split("\t")[2])
        if seq in seq_strandness_dict:
            output.append("\t".join([get_seq(dict_hg19,pos,strand),pos,strand]))
    return output

###main
seq_strandness_dict=make_seq_strandness_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test_sequence_of_sample25.csv")
dict_hg19=make_dict_hg19(r"/Share/home/zhangqw/SELECT/ref/hg19.fa")
output_plus=output("+",r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.fq2_rank50.std2.txt")
output_minus=output("-",r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.fq2_rank50.std2.txt")
output=output_plus+output_minus

f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/0.test_seq_on_rank50.txt","w")
f.write("seq\tpos\tstrand\n")
for line in output:
    f.write(line)
    f.write("\n")
f.close()
