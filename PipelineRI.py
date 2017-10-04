from CGAT import IOTools, GTF
from CGATPipelines.Pipeline import cluster_runnable
from CGATPipelines.Pipeline import run 
import itertools

def exportGTF(infile, outfile):
    '''This function will take an (optionally gziped) GTF
    and output is as a bgzipped, sorted, indexed GTF that
    can be easily visualised with IGV or UCSC'''


    if infile.endswith(".gz"):
        cat = "zcat"
    else:
        cat = "cat"

    statement = '''%(cat)s %(infile)s
                 | sort -k1,1 -k4,4n
                 | bgzip > %(outfile)s;
 
                 checkpoint;

                 tabix -p gff %(outfile)s'''

    run()


@cluster_runnable
def findRetainedIntrons(infile, outfile):

    outf = IOTools.openFile(outfile, "w")

    for gene in GTF.gene_iterator(
            GTF.iterator(IOTools.openFile(infile))):
        
        gene_out = []
        introns_out = []

        # now find if any of the transcripts are retained intron
        # versions of any of the others
        for first, second in itertools.product(gene, gene):
            
            first = sorted([entry for entry in first
                            if entry.feature == "exon"],
                           key=lambda x: x.start)
            second = sorted([entry for entry in second
                             if entry.feature == "exon"],
                            key=lambda x: x.start)

            first_introns = set(GTF.toIntronIntervals(first))
            second_introns = set(GTF.toIntronIntervals(second))
            
            if len(first_introns-second_introns) > 0 and \
               len(second_introns-first_introns) == 0:
                novel_introns = list(first_introns-second_introns)

                def _filterIntron(intron):
                    return intron[0] > second[0].start and \
                        intron[1] < second[-1].end

                novel_introns = filter(_filterIntron, novel_introns)

                if len(novel_introns) > 0:
                    gene_out.extend(first)

                for intron in novel_introns:
                    introns_out.append(intron)

        introns_out = Intervals.combine(introns_out)
        template = gene[0][0]
        template.feature = "exon"
        for gff in introns_out:
            entry = GTF.Entry().copy(template)
            entry.start = gff[0]
            entry.end = gff[1]
            outf.write("%s\n" % str(entry))
