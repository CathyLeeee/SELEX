from string import atoi

def whetehr_in(small,big):
    """
    function: test whetehr small region is covered by big region
    input: "chr2,240026798,240027072","chr2,45968738,45969123"
    ouput:  "samll region is covered by big region"
            or "samll region is not covered by big region" 
    """
    samll_start=atoi(small.split(",")[1])
    samll_end=atoi(small.split(",")[2])
    big_start=atoi(big.split(",")[1])
    big_end=atoi(big.split(",")[2])

    if    big_start<small_strat and big_end>small_end:
        print "samll region is covered by big region"
    else:
        print "samll region is not covered by big region"

if __name__=="__main__":
    small=["chr2,240026798,240027072","chr2,240026776,240027076","chr2,240026784,240027059","chr2,240026776,240027072","chr2,240026786,240027069","chr2,240026798,240027076",""]
