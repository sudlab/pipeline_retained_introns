##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline retained introns
===========================

:Author: Ian Sudbery
:Release: 0.0.1
:Date: |today|
:Tags: Python

This pipeline finds significant retension of introns between two
condiditons. It looks uses DEXSeq to find differential usage of
regions that are one of 1) Consituative introns 2) Retained introns
already annotated in the annotations 3) Sequences that are sometimes
eonx, sometimes intron.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

   python <srcdir>/pipeline_retained_introns.py config

Input files
-----------

* BAM files for each RNA-seq sample. Should be named -
  CELL-CONDITION-REPLICATE.bam (e.g. HEK293-Control-R1.bam)

* (optional) A design file with three columns. The first is the name
  of a comparison. The next two specify regex patterns for determining
  test and control samples respectively. Ig no design file present,
  some of the files are expected to have Control in the Condition
  slot.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:


Requirements:

* samtools >= 1.1
* CGAT
* CGATPipelines
* bedtools 
* R >= 3.1
    * DEXSeq
    * rtracklayer
    * GenomicRanges
    * experimentr (included)
* featureCounts
* bedGraphToBigWig

Pipeline output
===============

The primary output of the pipeline is a database, the name of which
can be specified by the "name" setting in the "database" section of
the report, or otherwise is "csvdb" in the working directory. The key
table in the database is "dexseq_results", which as the name suggests
contains the results of the DEXSeq run. Other tables reffer to
annotation of chunuks of transcripts as introns, exons or retained
introns.

The export directory will contain GTFs for the transcripts broken up
into chunks, the chunks annotated as retained introns and introns that
had a significant change between conditions, as well as bigwig files
for the read coverage.


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
from CGATPipelines import PipelineRnaseq
from CGAT import IOTools
import PipelineRI

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
   
    return dbh

REFERENCE_GENESET = PARAMS["geneset"]

# ---------------------------------------------------
# Specific pipeline tasks
@follows(mkdir("transcript_chunks.dir"))
@transform(REFERENCE_GENESET,
           regex("(.+)"),
                 r"transcript_chunks.dir/reference_chunks.gtf.gz")
def get_transcript_chunks(infile, outfile):
    '''Take reference transcriptome and transform so that each part of an
    exon or intron is an "exon" '''


    statement=''' cgat gtf2gtf --method=genes-to-unique-chunks
                  -I %(infile)s
                  -L %(outfile)s.log
                  -S %(outfile)s'''

    P.run()

###################################################################
@transform(get_transcript_chunks,
           formatter(),
           r"transcript_chunks.dir/filtered_chunks.gtf.gz")
def filter_overlapping_genes(infile, outfile):
    '''Truncate chunks where there is an overlap with more than one gene'''

    if PARAMS["stranded"] == 0:
        stranded = ""
    else:
        stranded = "-s"
        
    statement = '''bedtools intersect -a %(infile)s -b %(infile)s -c %(stranded)s
                 | sed -E 's/;([0-9]+)$/\\t\\1/g'
                 | awk -F'\\t' '$10==1'
                 | sed -E 's/(.+;)(\\t[0-9]+)$/\\1/g'
                 | gzip > %(outfile)s'''

    P.run()

    
###################################################################
@transform(REFERENCE_GENESET,
           formatter(),
           "transcript_chunks.dir/reference_retained_introns.gtf.gz")
def find_retained_introns(infile, outfile):

    PipelineRI.findRetainedIntrons(infile, outfile)


###################################################################
@transform(filter_overlapping_genes,
           suffix(".gtf.gz"),
           add_inputs(REFERENCE_GENESET),
           ".exons.tsv.gz")
def annotate_exon_chunks(infiles, outfile):
    '''Find those transcript chunks that represent  exons in some transcripts'''

    chunks, transcripts = infiles

    job_memory="6G"
    
    statement = ''' zcat %(transcripts)s
                 | awk '$3=="exon"'
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(\w+)\\".+exon_id \\"(\w+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\texon'
                 | gzip > %(outfile)s'''
    
    P.run()

    
###################################################################
@transform(filter_overlapping_genes,
           suffix(".gtf.gz"),
           add_inputs(REFERENCE_GENESET),
           ".introns.tsv.gz")
def annotated_intron_chunks(infiles, outfile):
    ''' Find those chunks that are intronic sequence in at least one transcript'''

    chunks, transcripts = infiles

    statement = '''
                   cgat gtf2gtf
                        -I %(transcripts)s
                        -L %(outfile)s.log
                        --method=set-gene-to-transcript
                 | cgat gtf2gtf
                        -L %(outfile)s.log
                        --method=exons2introns
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(\w+)\\".+exon_id \\"(\w+)\\".+([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\texon'
                 | gzip > %(outfile)s'''
    
    P.run()


###################################################################
@transform(filter_overlapping_genes,
           suffix(".gtf.gz"),
           add_inputs(find_retained_introns),
           ".retained_introns.tsv.gz")
def annotate_retained_introns(infiles, outfile):
    '''find chunks that overlap with annotated retained introns'''

    chunks, introns = infiles

    statement = ''' bedtools intersect -a %(chunks)s -b %(introns)s -s -c
                  | sed -E 's/.+gene_id \\"(\w+)\\".+exon_id \\"(\w+)\\".+([0-9]+)$/\\1\\t\\2\\t\\3/' 
                  | sed '1i gene_id\\texon_id\\texon'
                  | gzip > %(outfile)s'''

    P.run()

    
###################################################################
@transform([annotate_retained_introns,
            annotate_exon_chunks,
            annotated_intron_chunks],
           suffix(".tsv.gz"),
           ".load")
def load_chunk_annotations(infile, outfile):

    P.load(infile, outfile, "-i gene_id -i exon_id")

    tablename = P.toTable(outfile)
    connect().executescript('''DROP INDEX IF EXISTS %(tablename)s_joint;
                               CREATE INDEX %(tablename)s_joint ON
                                   %(tablename)s(gene_id,exon_id)'''
                            % locals())

###################################################################    
@follows(load_chunk_annotations)
def prepare_chunks():
    pass


###################################################################
###################################################################
# Differential expression
if os.path.exists("design.tsv"):
    @split("design.tsv", "*.design.tsv")
    def generate_dexseq_design_files(infile, outfiles):
        '''take the design specification for the pipeline and convert 
        into dexseq design matricies'''

        bamfiles = glob.glob("*.bam")
        bamfiles = [P.snip(os.path.basename(f), ".bam") for f in bamfiles]
        comparisons = [line.split() for line in IOTools.openFile("design.tsv")
                       if  not line.startswith("#")]

        
        for name, pat1, pat2 in comparisons:
            condition1_files = [(f, "test") for f in bamfiles
                                if re.match(f, pat1)]
            condition2_files = [(f, "control") for f in bamfiles
                                if re.match(f, pat2)]
            IOTools.writeLines("alt_utr_anlysis.dir/%s.design.tsv" % name,
                               condition1_files + condition2_files,
                                header=["track","condition"])

else:

    @collate("*.bam",
             regex("(.+)-((?!Control).+)-(.+).bam"),
             add_inputs(r"\1-Control-\3.bam"),
             r"\1-\2.design.txt")
    def generate_dexseq_design_files(infiles, outfile):

        track = os.path.basename(infiles[0][0]).split("-")[0]
        files = [P.snip(os.path.basename(f), ".bam")
                 for p in infiles for f in p]
        files = [(f, f.split("-")[1]) for f in files]
        IOTools.writeLines(outfile, files, header=["track","condition"])

@follows(mkdir("counts.dir"))        
@transform("*.bam", formatter(),
           add_inputs(filter_overlapping_genes),
           "counts.dir/{basename[0]}.tsv.gz")
def count_chunks(infiles, outfile):

    gtffile = infiles[1]
    bamfile = infiles[0]

    PipelineRnaseq.runFeatureCounts(
        gtffile,
        bamfile,
        outfile,
        job_threads=PARAMS["featurecounts_threads"],
        strand=PARAMS["stranded"],
        options="-f " + PARAMS["featurecounts_options"])


###################################################################
@merge(count_chunks,
       "counts.dir/chunk_counts.tsv.gz")
def merge_chunk_counts(infiles, outfile):

    infiles = " ".join(infiles)
    job_memory = "10G"
    statement=''' python %(scriptsdir)s/combine_tables.py
                         -c 1,2,3,4,5,6
                         -k 7
                         --regex-filename='(.+).tsv'
                         --use-file-prefix
                         --merge-overlapping
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


###################################################################
@follows(mkdir("dexseq.dir"))
@transform(generate_dexseq_design_files,
           regex("(.+).design.txt"),
           add_inputs(merge_chunk_counts,
                      filter_overlapping_genes,
                      annotated_intron_chunks),
           r"dexseq.dir/\1.dexseq.tsv")
def run_dexseq(infiles, outfile):
    '''run dexseq on the chunks'''

    design, counts, models, introns = infiles

    infiles = ",".join([models, counts, design, introns])
    outfile = P.snip(outfile, ".tsv")

    job_threads = PARAMS["dexseq_threads"]
    job_memory=PARAMS["dexseq_memory"]

    pipeline_src = os.path.dirname(__file__)
    script = os.path.join(pipeline_src, "run_dexseq_all.R")
    statement = ''' Rscript %(script)s
                            --padj=%(dexseq_padj)s
                            --lfc=%(dexseq_lfc)s 
                            --infiles %(infiles)s
                            --outfiles %(outfile)s.tsv,%(outfile)s.gtf.gz,%(outfile)s.RData
                             -p %(dexseq_threads)s
                    &> %(outfile)s.log '''

    P.run()


# -----------------------------------------------------------------
@merge(run_dexseq,
       "dexseq.dir/dexseq_results.load")
def load_dexseq(infiles, outfile):

    statement = " checkpoint;".join(
        [" sed 's/log2fold_\S+/log2fold/' %s > %s.tmp;" % (f, f)
         for f in infiles])

    P.run()

    infiles = ["%s.tmp" % f for f in infiles]
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).dexseq.tsv.tmp",
                         options="-i groupID -i featureID -i track -i padj",
                         job_memory="6G")

    for f in infiles:
        os.unlink(f)

        
@transform(load_dexseq, suffix(".load"), ".index")
def joint_index_dexseq(infile, outfile):

    db = connect()
    db.executescript('''
             DROP INDEX IF EXISTS dexseq_results_joint;
             CREATE INDEX dexseq_results_joint
                    ON dexseq_results(groupID,featureID);''')
    P.touch(outfile)


# -----------------------------------------------------------------
@follows(joint_index_dexseq,
         load_dexseq)
def dexseq():
    pass


###################################################################
###################################################################
# export

@follows(mkdir("export"))
@transform(filter_overlapping_genes,
           formatter(),
           "export/{basename[0]}.{ext[0]}")
def export_chunks(infile, outfile):
    '''Export transcript chunks in indexed fashion'''

    PipelineRI.exportGTF(infile, outfile)

    
###################################################################
@follows(mkdir("export"))
@transform(find_retained_introns,
           formatter(),
           "export/{basename[0]}.{ext[0]}")
def export_retained_introns(infile, outfile):
    '''Export annotated retained introns'''

    PipelineRI.exportGTF(infile, outfile)


###################################################################
@follows(mkdir("export"))
@transform(run_dexseq,
           formatter(),
           inputs("{path[0]}/{basename[0]}.gtf"),
           "export/{basename[0]}.sig_introns.gtf.gz")
def export_sig_introns(infile, outfile):

    PipelineRI.exportGTF(infile, outfile)


@follows(mkdir("export"))
@transform("*.bam",
           formatter(),
           "export/{basename[0]}.bw")
def export_bigwigs(infile, outfile):
    '''Export coverage plots of bam files to bigWigs for browser
    visualisation'''
 
    genome_file = PARAMS["contigs_tsv"]

    tmp = P.getTempFilename()
    statement = ''' genomeCoverageBed -split -bg -ibam %(infile)s
                                      -g %(genome_file)s  2> %(outfile)s.log
                   | sort -k1,1 -k2,2n > %(tmp)s ;
                    
                    checkpoint;

                    bedGraphToBigWig %(tmp)s
                                     %(genome_file)s
                                     %(outfile)s 2>>%(outfile)s.log;

                    checkpoint;

                    rm %(tmp)s'''

    P.run()
    
###################################################################
@follows( export_chunks,
           export_retained_introns,
           export_sig_introns,
           export_bigwigs)
def export():
    pass

# ---------------------------------------------------
# Generic pipeline tasks
@follows(prepare_chunks,
         dexseq,
         export)
def full():
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
