def strand(flag):
    """
    function: give +/- strand
    
    """
    l=list(bin(flag))
    l.reverse()
    if l[4]==0:
        strand="+"
    else:
        strand="-"


    return(strand)

if __name__=="__main__":
    flag=[99,147]
    print map(strand,flag)
