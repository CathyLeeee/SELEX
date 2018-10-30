#!usr/bin/python


from string import atoi

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
                        print i
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

if __name__=="__main__":
    filter_sam(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.new.sam",r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.filtered.new.sam")
