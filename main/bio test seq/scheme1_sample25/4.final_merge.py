#!usr/bin/python

from string import atoi


def  merge_read(l_big_read):
    """
    function: merge read in l_big_read
    input:["chr2,240026798,240027034",...]
    output:{"chr2":"240026798,240027034",...}
    """
    global chr,dict_chr_contig,temp_chr,temp_contig,contig,extend
    chr=list(set([i.split(",")[0]  for i in l_big_read]))
    chr=sorted(chr)

    func = lambda x,y:x if y in x else x + [y]
    dict_chr_contig={}
    for i in chr:
        temp_chr=[]
        for j in l_big_read:
            if j.split(",")[0]==i:
                temp_chr.append(j)
        temp_contig=[]
        for j in temp_chr:
            contig=[atoi(j.split(",")[1]),atoi(j.split(",")[2])]
            for m in temp_chr:
                extend=[atoi(m.split(",")[1]),atoi(m.split(",")[2])]
                if extend[0]<=contig[0] and extend[1]<=contig[1] and extend[1]>=contig[0]:
                    contig[0]=extend[0]
                elif extend[0]>=contig[0] and extend[0]<=contig[1] and extend[1]>=contig[1]:
                    contig[1]=extend[1]
                elif extend[0]<=contig[0] and extend[1]>=contig[1]:
                    contig[0]=extend[0]
                    contig[1]=extend[1]
            temp_contig.append(contig)
        temp_contig=reduce(func, [[], ] + temp_contig)   #rmdup
        dict_chr_contig[i]=temp_contig
    return chr , dict_chr_contig

def T2U(dna):
    """
    function: substitue T with U to transform from DNA to RNA

    """
    rna=[]
    for i in dna:
        if i=="A" or i=="a":
            rna.append("A")
        elif i=="T" or i=="t":
            rna.append("U")
        elif i=="C" or i=="c":
            rna.append("C")
        elif i=="G" or i=="g":
            rna.append("G")
        elif i=="N":
            rna.append("N")
        else:
            rna.append("*")
    return "".join(rna)

def rc(seq):
    """
    function: give the reverse complementary seq  & T is substitued with U

    """
    rc_seq=[]
    for i in seq:
        if i=="A" or i=="a":
            rc_seq.append("U")
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

def get_seq(dict_hg19,big_read):
    """
    function: grep sequece in ref according to big_read pos
    input: (dict_hg19,"chr2,240026798,240027034")
    output: "GGGGTTACCC"
    """
    chr=big_read.split(",")[0]
    start=atoi(big_read.split(",")[1])
    end=atoi(big_read.split(",")[2])
    if chr in dict_hg19:
        return dict_hg19[chr][start-1:end]
    else:
        return "NA"

def return_chr_len(pos):
    """
    function: return the fragement length
    input: chr2,240026904,240027072
    output: 169  #240027072-240026904+1=169
    """
    return atoi(pos.split(",")[2])-atoi(pos.split(",")[1])+1

    
def output(strand,inpath):
    """
    function: give output accoding to strandness
    input:"+"/"-"
    output:["seq,pos,strand"...]
    """
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    
    txt.pop(0)
   
    txt=[line.strip("\n")  for line in txt  if  return_chr_len(line.strip("\n").split("\t")[0])<500 and line.strip("\n").split("\t")[1]==strand]
    l_big_read=[i.split("\t")[0]   for i in txt ]
 
    l_chr , dict_chr_contig= merge_read(l_big_read)
    output=[]
    for chr in l_chr:
        for start_and_end in dict_chr_contig[chr]:
            pos=",".join([chr]+[str(i)  for i in start_and_end])
            if strand=="+":
                output.append("\t".join([T2U(get_seq(dict_hg19,pos)),pos,strand]))
            if strand=="-":
                output.append("\t".join([rc(get_seq(dict_hg19,pos)),pos,strand]))
    return output

###main
dict_hg19=make_dict_hg19(r"/Share/home/zhangqw/SELECT/ref/hg19.fa")
output_plus=output("+",r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/3.merge_reads.pos_and_strand.txt")
output_minus=output("-",r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/3.merge_reads.pos_and_strand.txt")
output=output_plus+output_minus

f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/4.final_merge.txt","w")
f.write("seq of RNA\tpos\tstrand\n")
for line in output:
    f.write(line)
    f.write("\n")
f.close()
    
    

