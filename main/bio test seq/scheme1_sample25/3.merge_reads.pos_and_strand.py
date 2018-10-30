#!usr/bin/python

from string import atoi
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

####sample25
def make_seq_strandness_dict(inpath):
    """
    function: make seq_strandness_dict according to sequecne with strandness file
    input: file path(QNAME in file coordinate with sam file)
    output: seq_strandness_dict
    #eg.seq_strandness_dict={"ATCGGGG":"-"}   #query seq : strandness
    attention: "-" seq takes rc operation to match sam file
    """
    global seq_strandness_dict
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    txt.pop(0)

    seq_strandness_dict_plus={x.split(",")[1].strip():x.split(",")[2].strip("\r\n")  for x in txt if x.split(",")[2].strip("\r\n")=="+"}
    seq_strandness_dict_minus={rc(x.split(",")[1].strip()):x.split(",")[2].strip("\r\n")  for x in txt  if x.split(",")[2].strip("\r\n")=="-"}
    seq_strandness_dict=dict(seq_strandness_dict_plus,**seq_strandness_dict_minus)
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

def return_pos_and_freq_dict(sam_path,seq_strandness_dict):
    """
    function: return pos_dict according to tsam
    input: path of sam file(QNAME in sam coordinate with seq_strandness_dict)
           seq_strandness_dict  #for filtering
    output: pos_dict  ;freq_dict
            pos_dict["ST-E00251:265:HFHNCALXX:6:1101:5832:1784"]=["chr2,240026885","chr2,240026798"...]
    """

    global pos_dict,freq_dict
    f=open(sam_path,"r")
    tsam=f.readlines()
    f.close()

    pos_dict={}
    freq_dict={}    

    for i in  tsam:
        i=i.strip("\n")
        temp=i.split()
        read=temp[9]
        if temp[6]=="=":
            temp[6]=temp[2]
        pos=";".join([",".join([temp[2],temp[3]]),",".join([temp[6],temp[7]])])
        if read in pos_dict:
            freq_dict[read]=freq_dict[read]+1
            if pos not in pos_dict[read]:
                pos_dict[read].append(big_read(pos))
        else:
            freq_dict[read]=0
            pos_dict[read]=[big_read(pos)]

    pos_dict={k:list(set(v)) for k,v in pos_dict.iteritems() if k in seq_strandness_dict }  #only seq exactly in fq
    return pos_dict,freq_dict

def rc(seq):
    """
    function: return the reverse complementary seq
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



def get_seq(big_read,strandness):
    """
    function: get base sequece in ref according to big_read pos and strandness based on Bio 
    input: ("chr2,240026798,240027034","+")
    output: "GGGGTTACCC"  (get seq from hg19 acooding to pos)
    """
    records = SeqIO.to_dict(SeqIO.parse(open('/Share/home/zhangqw/SELECT/ref/hg19.fa'), 'fasta'))
    chr,start,end=big_read.split(",")
    start,end=atoi(start),atoi(end)
   
    long_seq_record = records[chr]
    long_seq = long_seq_record.seq
    short_seq = str(long_seq)[start-1:end]
    return short_seq     

    #if strandness=="+":
    #    return dict_hg19[chr][start-1:end]
    #elif strandness=="-":
    #    return rc(dict_hg19[chr][start-1:end])

def  merge_read(l_big_read):
    """
    function: merge read in l_big_read
    input:["chr2,240026798,240027034",...]
    output:["chr2,240026798,240027034",...]
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

def write_output(outpath,chr_plus,dict_chr_contig_plus,chr_minus,dict_chr_contig_minus):
    f=open(outpath,"w")
    f.write("pos\tseq\tstrandness\n")
    for i in chr_plus:
        for j in dict_chr_contig_plus[i]:
            pos=",".join([i,str(j[0]),str(j[1])])
            seq=get_seq(pos,"+")
            f.write("\t".join([pos,seq,"+"]))
            f.write("\n")
    for i in chr_minus:
        for j in dict_chr_contig_minus[i]:
            pos=",".join([i,str(j[0]),str(j[1])])
            seq=get_seq(pos,"-")
            f.write("\t".join([pos,seq,"-"]))
            f.write("\n")
    f.close()
    
    
    

###main
seq_strandness_dict=make_seq_strandness_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test_sequence_of_sample25.csv")
pos_dict,freq_dict=return_pos_and_freq_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.filtered.new.sam",seq_strandness_dict)


f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/3.merge_reads.pos_and_strand.txt","w")
f.write("pos\tstrand\n")
for key in pos_dict:
    for i in pos_dict[key]:
        big_read=i
        strandness=seq_strandness_dict[key]
        f.write("\t".join([big_read,strandness]))
        f.write("\n")
f.close()

#chr_plus,dict_chr_contig_plus=merge_read([i.split("\t")[0]  for i in temp if i.split("\t")[1]=="+"])
#chr_minus,dict_chr_contig_minus=merge_read([i.split("\t")[0]  for i in temp if i.split("\t")[1]=="-"])

#write_output(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/3.merge_reads.pos_and_strand.txt",chr_plus,dict_chr_contig_plus,chr_minus,dict_chr_contig_minus)


#seq_strandness_dict=make_seq_strandness_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test_sequence_of_sample25.csv")
#f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/step2_fq2.merge_reads.pos.csv","w")
#f.write("chr,start,end,strand\n")
#for seq in pos_dict:
#    if seq 
#    value=pos_dict[key]
#    if len(value)==2:
#        pe=";".join(value)
#        f.write(",".join([big_read(pe),seq_strandness_dict[key],get_seq(big_read(pe),seq_strandness_dict[key])]))
#        f.write("\n")
#f.close()
    
    
    
    
    

