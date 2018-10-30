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
        

if __name__=="__main__":
    seq="AAAGAGAAGGTTTCTCAGATCAGAGAGCCTTCATGGGAAATTTTCTTACTTGAAAGCACTCTTAGATGTTCTCTTCCCTTTTACTCAACAAAGAGCCAGAGTGCTTGTGATCTGGCAAAGGAAGATGAATGGGGCAGAATCACTGCCTGGAAAGATCTTTGCAGAGAGGGACGTAACCCTCCGTCCTGAAAAAGGCCAGACCTAACCTGCTCGGCCTCGCCTCCCTCCCTGCCCTCACAGAGGCCCTCGGCTGGGTGAGGGCCTGGCCTTGTTGCTCAGAATTTCCAGACCAGCCTCTTTTTTGCCCTCCTGGTTACATTCCCCAGGT"
    print rc(seq)
