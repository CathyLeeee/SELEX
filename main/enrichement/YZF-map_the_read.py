import sys
from pyensembl import EnsemblRelease
from matplotlib import pyplot as plt
from matplotlib import patches as patches
import math
from matplotlib import gridspec

chr_relative_size ={'chrY': 0.2382082971821362, 'chrX': 0.6229495412169906, 'chr13': 0.4620645579053542,
                    'chr12': 0.5370172979428605, 'chr11': 0.5416496675448604, 'chr10': 0.5437689441102737,
                    'chr17': 0.3257573027270412, 'chr16': 0.36250562842128287, 'chr15': 0.41135862205133683,
                    'chr14': 0.4306891576410556, 'chr19': 0.2372270237994713, 'chr18': 0.313247957765369,
                    'chr22': 0.20583525848065992, 'chr20': 0.2528600319916555, 'chr21': 0.19309839553017602,
                    'chr7': 0.638468471458693, 'chr6': 0.6865181170401177, 'chr5': 0.7258367472633097,
                    'chr4': 0.7669159468212519, 'chr3': 0.7944711600136796, 'chr2': 0.9757222350110012,
                    'chr1': 1.0, 'chr9': 0.5665519726027082, 'chr8': 0.587216278189333, 'chr1_length':249250621}


def map_the_read(infile,outfile):
    f = open(infile,"r")
    l = f.readlines()
    l = [x[:-1] for x in l]
    l = [x.split("\t") for x in l]
    total_reads = int(l[-1][1])
    l = l[:-1]
    # import the ensembl gene annotation
    genome = EnsemblRelease(75)
    # parse the file
    seq_lib = {}
    num = 0
    num_chimeric = 0 #chimeric alignment #
    num_no_map = 0 # seqs without mapping info
    num_chrM = 0 # read that mapped to the mitochondria
    l_out = []
    for seq in l:
        num += 1
        frequency = seq[1]
        try:
            lsub = seq[2]
            lsub = lsub.split(";")
            lsub = [x.split(",") for x in lsub]
            if lsub[0][0] == " * chrM":
                num_chrM +=1
                l_out.append(' '.join(seq)+"chrM\n")
            elif lsub[0][0] == "*":
                num_no_map+=1
                l_out.append(' '.join(seq)+" * no match\n")
            else:
                if lsub[0][0] == lsub[1][0]:
                    chr = lsub[0][0][3:]
                    seq_start = lsub[0][1]
                    seq_end = lsub[1][1]
                    seq_position = (int(seq_start)+int(seq_end))/2   # use the average position for annotation
                    gene_name = genome.gene_names_at_locus(contig=chr,position=seq_position)
                    if gene_name == []:
                        gene = "non-coding" # in case the sequence can't be mapped to a gene but none-coding part
                        l_out.append(' '.join(seq) + " non-coding\n")
                    else:
                        gene = gene_name[0] # seq that can be mapped to the gene
                        l_out.append(' '.join(seq) + "%s\n"%gene)
                    if (chr,gene) in seq_lib:
                        seq_lib[(chr,gene)].append((num,frequency,seq_position))
                    else:
                        seq_lib[(chr,gene)] = [(num,frequency,seq_position),]
                else:
                    num_chimeric += 1
                    l_out.append(' '.join(seq) + " * chimeric\n")
        except IndexError:
            num_no_map += 1
            l_out.append(' '.join(seq) + " * no match\n")
    #write the file out
    f = open(outfile,"w")
    for i in l_out:
        f.write(i)
    f.write("total reads number: %d"%total_reads)
    f.close()

    return seq_lib,num_chimeric,num_no_map,num_chrM,total_reads

def plot_the_chrs(lib,ax,length=50.0,x_default=10.0,y_default=10.0):
    l_chr = []
    for i in lib:
        l_chr.append("chr"+i[0])
    l_chr = list(set(l_chr))
    l_chr = sorted(l_chr, key=lambda x: int(x[3:]) if x[3:].isdigit() else x)
    # set the axes
    width = length+20
    height = float(y_default*(len(l_chr)+1))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    # draw the chrs
    lib_chr_y = {}
    for i in range(len(l_chr)):
        chr_length = length*chr_relative_size[l_chr[i]]
        chr_y = (i+1)*y_default
        lib_chr_y[l_chr[i]] = chr_y
        ax.text(
            x_default-2,
            chr_y,
            l_chr[i],
            horizontalalignment='right',
            verticalalignment='center',
            fontsize = 25,
            fontweight="bold"
        )
        ax.plot(
            (x_default, chr_length+x_default),
            (chr_y,chr_y),
            color="black",
            linewidth=8,
            alpha=0.5,
        )
    return lib_chr_y, height

def plot_the_genes_mapping(lib,lib_char_y,ax,total_read,length=50.0,x_default=10.0,frequency_parameter=0.05):
    #parse the data
    l_gene = []
    for i in lib:
        chr = "chr"+i[0]
        gene_name = i[1]
        if gene_name == "non-coding":
            for j in lib[i]:
                frequency = [float(j[1])*100/total_read,]
                position = j[2]
                l_gene.append((chr,None,frequency,position))
        else:
            frequency = [float(x[1])*100/total_read for x in lib[i]]
            position = [int(x[2]) for x in lib[i]]
            l_gene.append((chr,gene_name,frequency,position))
    # draw the genes
    for i in l_gene:
        gene_name = i[1]
        gene_y = lib_char_y[i[0]]
        frequency = sum(i[2])
        chr_length = chr_relative_size[i[0]] * chr_relative_size["chr1_length"]
        #find the x coordinate of the gene in the corresponding chrs
        if gene_name != None:
            genome = EnsemblRelease(75)
            gene = genome.genes_by_name(gene_name)
            gene = gene[0]
            gene_start = gene.start
            gene_end = gene.end
            gene_position_in_chr = (gene_start+gene_end)/2
            gene_x = length*chr_relative_size[i[0]]*gene_position_in_chr/chr_length
        else:
            gene_x = length*chr_relative_size[i[0]]*i[3]/chr_length
        # draw the genes (size of circle stand for the freqency
        if frequency < 0.05:
            radius = 0.5
        else:
            radius = abs(math.log(frequency / frequency_parameter))
        if gene_name == None:
            color = "deeppink"
        else:
            color = "dodgerblue"
        # gene with frequency related radius
        ax.add_patch(
            patches.Circle(
                (gene_x+x_default,gene_y),
                radius,
                color=color,
                alpha=0.5,
                linewidth=0,
            )
        )
        # gene position
        ax.add_patch(
            patches.Circle(
                (gene_x+x_default,gene_y),
                0.5,
                color=color,
                linewidth=0.5,
                edgecolor="black",
            )
        )

        if gene_name != None:
            ax.text(
                gene_x + x_default,
                gene_y-radius-0.8,
                gene_name,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=18,
                fontweight="bold",
            )
        if frequency < 0.0001:
            number = "%.2e"%float(frequency)
        else:
            number="%.4f"%float(frequency)
        ax.text(
            gene_x + x_default,
            gene_y + radius +0.5,
            number,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=12,
        )
    return l_gene

def get_exon_location():
    # import the exon location info
    f = open("ensembl_exon_boundary.tsv",'r')
    l = f.readlines()
    l = l[1:]
    l = [x.split("\t") for x in l]
    exon_lib = {}
    for i in l:
        exon_lib[i[3][:-1]] = (i[1],i[2])
    return exon_lib

def make_the_gene_with_exon_and_map(chr,gene_name,ax,exon_lib,x,y,locus,frequency,width=30.0,height=5.0,scale=1.0):
    genome = EnsemblRelease(75)
    genes = genome.genes_by_name(gene_name)
    gene = genes[0]
    gene_start = gene.start
    gene_end = gene.end
    gene_length = float(gene_end-gene_start)
    # draw the gene
    ax.add_patch(
        patches.Rectangle(
            (x, y-height/2),
            width,
            height,
            alpha = 0.5,
            color="dodgerblue",
            linewidth=0,
        )
    )
    # find the exon ids
    exon_id = genome.exon_ids_of_gene_name(gene_name)
    # plot the gene with exons
    for i in exon_id:
        position = exon_lib[i]
        exon_start = int(position[0])
        exon_end = int(position[1])
        exon_length = exon_end - exon_start
        plot_length = exon_length/gene_length*width
        plot_x = (exon_start-gene_start)/gene_length*width
        ax.add_patch(
            patches.Rectangle(
                (x + plot_x, y - height / 2-0.5*scale),
                plot_length,
                height+1*scale,
                color ="black",
                alpha=0.4,
                linewidth=0,
            )
        )
    # map the read to the plot
    for i in range(len(locus)):
        plot_x = (locus[0]-gene_start)/gene_length*width
        if frequency[i] < 0.0001:
            bar_length = 1
        else:
            bar_length = abs(math.log(frequency[i] / 0.01) * 2)
        ax.plot(
            (x+plot_x,x+plot_x),
            (y+height/2+scale*0.5,y+height/2+scale*0.5+bar_length),
            color = "darkblue",
            alpha=0.5,
            linewidth=2,
        )
    # gene name at the left
    ax.text(
        x-0.5,
        y,
        gene_name+'\n'+chr,
        horizontalalignment='right',
        verticalalignment='center',
        fontsize=25,
        fontweight="bold"
    )
    # gene start
    ax.text(
        x,
        y-height/2-0.8*scale,
        gene_start,
        horizontalalignment='center',
        verticalalignment='top',
        fontsize=12,
    )
    # gene end
    ax.text(
        x+width,
        y-height/2-0.8*scale,
        gene_end,
        horizontalalignment='center',
        verticalalignment='top',
        fontsize=12,
    )

def plot_the_genes(l_gene,height,ax,exon_lib,x_default=7.0,y_default=15.0):
    l_gene = [x for x in l_gene if x[1] != None]
    l_gene = sorted(l_gene,key=lambda x: sum(x[2]))
    if height < (len(l_gene)+1)*y_default:
        scale = float(height)/((len(l_gene)+1)*y_default)
        y_default = float(height)/(len(l_gene)+1)
    else:
        scale=1
    width = 40
    ax.set_xlim(0,width)
    ax.set_ylim(0,height)
    num = 1
    for i in l_gene:
        chr = i[0]
        gene_name = i[1]
        frequency = i[2]
        locus = i[3]
        x = x_default
        y = num*y_default
        make_the_gene_with_exon_and_map(chr,gene_name,ax,exon_lib,x,y,locus,frequency,height=y_default/3,scale=scale)
        num+=1

if __name__ == "__main__":
    # the input and output files
    infile = sys.argv[1]
    if len(sys.argv) >= 3:
        outfile = sys.argv[2]
    else:
        outfile = infile[:-4]+'.png'
    outfile2 = infile[:-4]+'_annotation.txt'
    # plot the figure
    fig = plt.figure()
    # two columns of the figure 1 for mapping to the chr, 2 for specific gene
    gs = gridspec.GridSpec(1,2,width_ratios=[1,3])
    ax = fig.add_subplot(gs[0],frameon=False) # hide the frame
    ax2 = fig.add_subplot(gs[1],frameon=False) # hide the frame
    seq_lib, num_chimeric, num_no_map, num_chrM, total_reads = map_the_read(infile,outfile2)
    # write the mapping info
    ax.text(
        40,
        0,
        "No.seq: 50, No.map: %d, No.chrM: %d, No.chimera: %d, No.no map: %d"%(50-num_no_map-num_chimeric-num_chrM,num_chrM,num_chimeric,num_no_map),
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=25,
    )

    # plot the chr
    lib_chr_y, height = plot_the_chrs(seq_lib,ax)
    # map the gene to the chr by chr id and position and frequency
    l_gene = plot_the_genes_mapping(seq_lib,lib_chr_y,ax,total_reads)
    # find the exomn boundary
    exon_lib = get_exon_location()
    # plot the genes
    plot_the_genes(l_gene,height,ax2,exon_lib)
    # hide the ticks
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax2.get_xaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])
    plt.tight_layout()
    #plt.show()
    # save the picture
    height=2*len(lib_chr_y)
    width=50
    fig.set_size_inches(width,height)
    fig.savefig(outfile,dpi=200,bbox_inches="tight")
