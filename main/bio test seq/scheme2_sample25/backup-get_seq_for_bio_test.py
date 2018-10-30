#!usr/bin/python
"""function:get seq for bio test 
   data precess: 1.resort sam according to strandness by blat  
                 (input file:fq1, fq2, blat;  
                  ouput file: resorted  fq1 and fq2)
                 2.filter sam
                 (input & output :sam)
                 3.merge(collect RNAME for each seq, merge pos for each seq, get seq according to merged pos)
   **get seq from hg19 accodidng to pos is time-costing step, so do as few times as possible**
"""

from string import atoi
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from compiler.ast import flatten


def resort_fq(path_of_raw_fq1,path_of_raw_fq2,path_of_resorted_fq1,path_of_resorted_fq2,seq_strandness_dict):
    """
    function: resort sam according to strandness by blat
    input: path_of_raw_fq1,path_of_raw_fq2,seq_strandness_dict
    output: path_of_resorted_fq1,path_of_resorted_fq2
    """
    global l_fq1,l_fq2,new_l_fq1,new_l_fq2
    l_fq1=make_l_fq(path_of_raw_fq1)
    l_fq2=make_l_fq(path_of_raw_fq2)
    new_l_fq1=[]
    new_l_fq2=[]
    for seq in seq_strandness_dict:
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
    print "".join(["len of new new_l_fq1:",str(len(new_l_fq1))])
    print "".join(["len of new new_l_fq2:",str(len(new_l_fq2))])
    write_new_fq(path_of_resorted_fq1,new_l_fq1)
    write_new_fq(path_of_resorted_fq2,new_l_fq2)

def make_seq_RNAME_dict(path_of_resorted_fq1 , path_of_resorted_fq2, seq_strandness_dict):
    """
    function: give the relationship between seq and RNAME in fq
    input: path_of_resorted_fq1 , path_of_resorted_fq2, seq_strandness_dict
    output: seq_RNAME_dict={"ATTTGCCC":[RNAME1,RNAME2...]}
    attention: seq in seq_RNAME_dict comes from seq_strandness_dict
    """
    global l_temp_QNAME_for_each_seq,seq_RNAME_dict
    l_fq1=make_l_fq(path_of_resorted_fq1)
    l_fq2=make_l_fq(path_of_resorted_fq2)
    l_fq2[1::4]=[rc(x) for x in l_fq2[1::4]]   #take rc to be insistent with sam
    l_fq=l_fq1+l_fq2
    seq_RNAME_dict={}
    for seq in seq_strandness_dict :
        l_temp_QNAME_for_each_seq=[]
        l_index=[i for i, x in enumerate(l_fq) if x==seq]
        for index in l_index:
            l_temp_QNAME_for_each_seq.append(l_fq[index-1].split()[0].strip("@"))
        l_temp_QNAME_for_each_seq=list(set(l_temp_QNAME_for_each_seq))
        seq_RNAME_dict[seq]=l_temp_QNAME_for_each_seq
    return seq_RNAME_dict
        
        
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

def make_seq_QNAME_dict(inpath):
    """
    function: make seq_QNAME_dict according to sequecne with strandness file
    input: file path(QNAME in file coordinate with sam file)
    output: seq_QNAME_dict
    #eg.seq_strandness_dict={"ATCGGGG":"seq1"}   #query seq : QNAME
    attention: "-" seq takes rc operation to match sam file
    """
    global seq_QNAME_dict
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()
    txt.pop(0)

    seq_QNAME_dict_plus={x.split(",")[1].strip():x.split(",")[0]  for x in txt if x.split(",")[2].strip("\r\n")=="+"}
    seq_QNAME_dict_minus={rc(x.split(",")[1].strip()):x.split(",")[0]  for x in txt  if x.split(",")[2].strip("\r\n")=="-"}
    seq_QNAME_dict=dict(seq_QNAME_dict_plus,**seq_QNAME_dict_minus)
    return seq_QNAME_dict



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
    
def filter_sam(inpath,outpath):
    """
    function: filter sam on 4 cateria (non-chimeric by "SA"; mapped; non-clip; pe are mapped to the same chr)
    input: inpath of sam to be filtered
    outpath
    """
    f=open(inpath,"r")
    txt=f.readlines()
    f.close()

    n_total=0
    n_chimeric=0
    n_unmap=0
    n_clip=0
    n_pe_mapped_to_diff_chr=0
    f=open(outpath,"w")
    for i in txt:
        n_total+=1

        temp=i.strip("\n").split()
        rname=i.split()[0]   #read name
        if not ("@" in i):  #skip header
            i=i.strip("\n")
            line=i.split()
            mrnm=line[6]    #MRNM; "=" indicate mate1 and mate2 are mapped on the same chr
            flag=atoi(line[1])    #FLAG; bitwise

            ###decision1: 2 mates are mapped to the same chr
            if not("SA" in i):        ###decision1: "SA" tag for filtering chimeric alignment
                l=list(bin(flag))   #decimal to binary
                l.reverse()
                if l[2]=="0": ###decison2: current read is  mapped
                    if not("S" in temp[5]) and not("H" in temp[5]): ###decison3: neither hard clip nor soft clip
                            if  mrnm=="=":   ###decison4: pe are mapped to the same chr
                                f.write(i)
                                f.write("\n")
                            else:
                                n_pe_mapped_to_diff_chr+=1
                    else:
                        n_clip+=1
                else:
                    n_unmap+=1
            else:
                n_chimeric+=1


    f.close()
    print "".join(["n_total:\t",str(n_total)])
    print "".join(["n_clip:\t",str(n_clip)])
    print "".join(["n_unmap:\t",str(n_unmap)])
    print "".join(["n_chimeric:\t",str(n_chimeric)])
    print "".join(["n_pe_mapped_to_diff_chr:\t",str(n_pe_mapped_to_diff_chr)])

def return_pos_dict(sam_path):
    """
    function: return pos_dict according to tsam
    input: path of sam file(QNAME in sam coordinate with seq_strandness_dict)
    output: pos_dict  
            pos_dict["ST-E00251:265:HFHNCALXX:6:1101:5832:1784"]=['chr2,240026875,240027060', 'chr2,240026869,240027060'...]
    """

    global pos_dict
    f=open(sam_path,"r")
    tsam=f.readlines()
    f.close()

    pos_dict={}

    for i in  tsam:
        i=i.strip("\n")
        temp=i.split()
        RNAME=temp[0]
        if temp[6]=="=":
            temp[6]=temp[2]
        pos=big_read(";".join([",".join([temp[2],temp[3]]),",".join([temp[6],temp[7]])]))
        if RNAME in pos_dict:
            if pos not in pos_dict[RNAME]:
                pos_dict[RNAME].append(pos)
        else:
            pos_dict[RNAME]=[pos]

    #pos_dict={k:list(set(v)) for k,v in pos_dict.iteritems() if k in seq_strandness_dict }  #only seq exactly in fq
    return pos_dict

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

def get_seq(big_read):
    """
    function: grep sequece in ref according to big_read pos dependent on Bio
    input: ("chr2,240026798,240027034")
    output: "GGGGTTACCC"  (get seq from hg19 acooding to pos)
    """
    records = SeqIO.to_dict(SeqIO.parse(open('/Share/home/zhangqw/SELECT/ref/hg19.fa'), 'fasta'))
    chr,start,end=big_read.split(",")
    start,end=atoi(start),atoi(end)

    long_seq_record = records[chr]
    long_seq = long_seq_record.seq
    short_seq = str(long_seq)[start-1:end]
    return short_seq

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

def  merge_read(l_big_read):
    """
    function: merge read in l_big_read
    input:["chr2,240026798,240027034",...]
    output:["chr2,240026798,240027034",...]
    """
    global chr,dict_chr_contig,temp_chr,temp_contig,contig,extendi,l_merge_read
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
        ##dict_chr_contig={'chr7': [[151084201, 151084500]]}
    l_merge_read=[]
    for i in chr:
        for j in dict_chr_contig[i]:
            l_merge_read.append(",".join([i,",".join(map(str,j))]))
    return l_merge_read

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


def merge_and_write_output(pos_dict,seq_strandness_dict,seq_QNAME_dict,seq_RNAME_dict,outpath):
    """
    function :write output
    input: four dict (pos_dict,seq_strandness_dict,seq_QNAME_dict,seq_RNAME_dict)
           outpath
    """
    global merge_read
    f=open(outpath,"w")
    for seq in seq_strandness_dict:
        strandness=seq_strandness_dict[seq]
        f.write("\t".join([seq_QNAME_dict[seq],strandness]))
        f.write("\n")
        f.write("position in hg19\tRNA sequence")
        f.write("\n")
        l_big_read=flatten([pos_dict[RNAME]  for RNAME in seq_RNAME_dict[seq] if RNAME in pos_dict])
        print l_big_read
        l_merge_read=merge_read(l_big_read)
        sequence_for_test="NA"
        for pos in l_merge_read:
            if strandness=="+":
                sequence_for_test=T2U(get_seq(pos))
            elif strandness=="-":
                sequence_for_test=T2U(rc(get_seq(pos)))
            f.write("\t".join([pos,sequence_for_test]))
            f.write("\n")
    f.close()



if __name__=="__main__":
    ###step1 resort fq
    seq_strandness_dict=make_seq_strandness_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test_sequence_of_sample25.csv")

    path_of_raw_fq1=r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.fastq"
    path_of_raw_fq2=r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.fastq"
    path_of_resorted_fq1=r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R1.clean_pe.new.fastq"
    path_of_resorted_fq2=r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25_R2.clean_pe.new.fastq"
    #resort_fq(path_of_raw_fq1,path_of_raw_fq2,path_of_resorted_fq1,path_of_resorted_fq2,seq_strandness_dict)

    ###step2 filter sam
    #filter_sam(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.new.sam",r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.filtered.new.sam")

    ###step3 merge and write output
    ##step3.1 relate data in different file by biuld dict
    seq_QNAME_dict=make_seq_QNAME_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/test_sequence_of_sample25.csv")
    seq_RNAME_dict=make_seq_RNAME_dict(path_of_resorted_fq1 , path_of_resorted_fq2, seq_strandness_dict)
    pos_dict=return_pos_dict(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.filtered.new.sam")
    ##step3.1 merge and write output
    #merge_and_write_output(pos_dict,seq_strandness_dict,seq_QNAME_dict,seq_RNAME_dict,r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/get_seq_for_bio_test.txt") 
    merge_and_write_output(pos_dict,seq_strandness_dict,seq_QNAME_dict,seq_RNAME_dict,r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/get_seq_for_bio_test.without_seq.txt") 







