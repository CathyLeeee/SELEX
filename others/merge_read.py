from string import atoi

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
    l_merge_read=[]
    for i in chr:
        for j in dict_chr_contig[i]:
            l_merge_read.append(",".join([i,",".join(map(str,j))]))
    return l_merge_read
        

if __name__=="__main__":
    minus1=["chr2,240026798,240027072","chr2,240026776,240027076","chr2,240026784,240027059","chr2,240026776,240027072","chr2,240026786,240027069","chr2,240026798,240027076","chr2,240026784,240027076"]
    print merge_read(minus1)
