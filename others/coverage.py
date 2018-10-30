#!usr/bin/python


import string 
##find reads in genome
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.genomecov.txt","r")
txt=f.readlines()
f.close()




chr=['chrY', 'chrX', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17', 'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chrM', 'chr22', 'chr20', 'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1', 'chr9', 'chr8']



depth_gt_1_dict={}
for i in chr:
    depth_gt_1_dict[i]=0

for i in txt:
    temp=i.strip("\n")
    temp=temp.split()
    l=string.atoi(temp[2])-string.atoi(temp[1])+1  #length
    c=temp[0]   #chromosome
    for j in chr:
        if len(temp[3])<=5:
            if (c==j) & (string.atoi(temp[3])>=1):
                depth_gt_1_dict[c]=l+depth_gt_1_dict[c]    
            
f=open(r"/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/LHT_ZM_I25.coverage.txt","w")
for i in chr:
    f.write("\t".join([i,str(depth_gt_1_dict[i])]))
    f.write("\n")
f.close()



