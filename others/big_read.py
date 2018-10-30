from string import atoi

def big_read(pe):
    """
    function: merge reads1+reads2 to big_reads
    input: "chr2,240026885;chr2,240026798"
    ouput: "chr2,240026798,240027034"
    """
    global pos
    if pe!="NA":
        pos=pe.split(";")
    chr=pos[0].split(",")[0]
    pos1=min(atoi(pos[0].split(",")[1]),atoi(pos[1].split(",")[1]))
    pos2=max(atoi(pos[0].split(",")[1]),atoi(pos[1].split(",")[1]))
    return ",".join([chr,str(pos1),str(pos2+149)])
        

if __name__=="__main__":
    pe="chr18,38954034;chr18,38954106"
    print big_read(pe)
