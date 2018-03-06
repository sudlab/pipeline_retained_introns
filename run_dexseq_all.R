library(rtracklayer)
library(DEXSeq)
library(experimentr)
library(GenomicRanges)
library(optparse)

opts = list(make_option(c("-p","--processes"),
                        type="integer",
                        dest="processes",
                        default=1,
                        help="Number of procsses to use. Using multiple processes requires the BiocParallel library"),
	   make_option(c("--padj"),
	               type="double",
		       dest="padj",
		       default=0.05,
		       help="Adjusted pvalue threshold to use when outputting retained intron GTF"),
	   make_option(c("--lfc"),
	               type="double",
		       dest="lfc",
		       default=1,
		       help="Adjusted log fold change threshold to use when outputting retained intron GTF"))

args = Experiment.start(opts)
gffs <- import(args$infiles[1])

col_data = read.delim(args$infile[3], header=T, row.names = 1)
col_data$condition <- relevel(col_data$condition, "Control")
rownames(col_data) <- gsub("-","_",rownames(col_data))

counts <- read.delim(args$infiles[2])
names(counts) <- gsub(".","_",names(counts), fixed=TRUE)
counts <- counts[, rownames(col_data)]
counts <- as.matrix(counts)

cat("# design is:\n")
print(col_data)
cat('# got ', dim(counts)[2], ' samples\n')
cat('# building data set\n')

introns <- read.delim(args$infiles[4])

dxd <- DEXSeqDataSet(counts, col_data, design= ~ sample + exon + condition:exon,
                           featureID = mcols(gffs)$exon_id, 
                           groupID = mcols(gffs)$gene_id,
                           featureRanges=gffs)

cat('# computing results\n')
if (args$processes > 1) {
  BPPARAM = MulticoreParam(workers=args$processes)
  dxd_results <- DEXSeq(dxd, BPPARAM=BPPARAM, quiet=F)
} else {
  dxd_results <- DEXSeq(dxd, quiet=F)
}

save.image(args$outfiles[3])
cat('#Subset results to only introns')
keep = paste(introns$gene_id[introns$exon>0], introns$exon_id[introns$exon>0], sep=":")
dxd_results <- dxd_results[keep,]
dxd_results$padj <- p.adjust(dxd_results$pvalue, method="BH")
geneQs <- perGeneQValue(dxd_results)
dxd_results$geneQ = geneQs[match(dxd_results$groupID, names(geneQs))]

sig_exons <- rownames(subset(dxd_results, padj < args$padj & abs(dxd_results[,10]) > args$lfc ))
cat('# got', length(sig_exons), ' significant exons\n')
sig_granges <- rowRanges(dxd[sig_exons,])
dxd_df <-as.data.frame(dxd_results)
cat('# outputting results')
write.table(dxd_df[,1:15],args$outfiles[1],
            quote=F, sep = "\t", row.names=FALSE)
export(sig_granges, args$outfiles[2])
save.image(args$outfiles[3])
Experiment.stop()
