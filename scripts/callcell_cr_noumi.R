library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Calculate the sum of a set of numbers")
p <- add_argument(p, arg = "--input",short = '-i', help="sparse barcode-gene matrix file to be input", type="character")
p <- add_argument(p, arg = "--output",short = '-o', help="list of barcodes called as cells", type="character")
p <- add_argument(p, arg = "--threshold",short = '-t', help="threshold of gene counts for barcodes to be specified as background", type="numeric",default = 20)
p <- add_argument(p, arg = "--fdr",short = '-f', help="FDR for testing whether a barcode is a empty cell", type="numeric",default = 0.01)
p <- add_argument(p, arg = "--expectedcell",short = '-e', help="Expected number of cells", type="numeric")
p <- add_argument(p, arg = "--upperquant",short = '-u', help="Upper quantile value for determine the cell number", type="numeric",default = 0.99)
p <- add_argument(p, arg = "--lowerprop",short = '-l', help="lower proprotion value for determine the cell number", type="numeric",default = 0.1)
p <- add_argument(p, arg = "--genome", help="Genome",type="character")

# Parse the command line arguments
argv <- parse_args(p)

# Calling cells
#f = read.table(argv$input,header = F,sep = '\t')
f = read.table(argv$input,header = F,sep = '\t', skip = 1)
library(Matrix,quietly=TRUE)
gene = as.factor(as.character(f[['V1']]))
A <- sparseMatrix(as.numeric(gene),as.numeric(as.factor(f[['V2']])),  x = f[['V3']])

library(DropletUtils,quietly=TRUE)

bcrank = barcodeRanks(A,lower = argv$threshold)

uniq = !duplicated(bcrank$rank)
pdf(paste(argv$output,"/Barcode_rank_plot.pdf",sep = ""))
    par(mar=c(5,4,2,1), bty="n")
    plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", 
              xlab="Rank", ylab="Total UMI count", cex=0.5, cex.lab=1.2)

    abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
    abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

    legend("left", legend=c("Inflection", "Knee"), bty="n", 
                  col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

dev.off()

O <- emptyDrops(A,lower = argv$threshold)
noempty <- (O$FDR <= argv$fdr)
noempty[is.na(noempty)] <- FALSE

trad <- defaultDrops(A,expected = argv$expectedcell,upper.quant=argv$upperquant, lower.prop=argv$lowerprop)
trad[is.na(trad)] <- FALSE

cell = trad | noempty

write.table(levels(factor(f[['V2']]))[cell], file = paste(argv$output,"/cell.list",sep = ""), sep = ',',row.names = F, col.names = F,quote=F)



colnames(A) <- levels(factor(f[['V2']]))
rownames(A) <- levels(factor(f[['V1']]))



A_cell <- A[,cell]
A_cell <- A_cell[apply(A_cell,1,sum) > 0,]
C_cell <- A_cell
C_cell[A_cell>0] <- 1


cellumi = apply(A_cell,2,sum)

cellgene = apply(C_cell,2,sum)

out <- data.frame("Items" = c("Estimated Number of Cells",'Total Genes Detected','Mean Genes per Cell','Median Genes per Cell','Mean UMI Counts per Cell','Median UMI Counts per Cell','Fraction Reads in Cells'),"Value" = c(dim(A_cell)[2],dim(A_cell)[1],mean(cellgene),median(cellgene),mean(cellumi),median(cellumi),sum(cellumi)/sum(f[['V4']])))
write.table(out,file = paste(argv$output,"/stat.tsv",sep = ""),row.names = FALSE,,sep = "\t")
write10xCounts(A_cell, path = paste(argv$output,"/mtx.h5",sep = ""), type = "HDF5", version = "3",genome = argv$genome,overwrite = TRUE)

dir_name=paste("mkdir -p ",argv$output,"/10x",sep = "")
system(dir_name)
write10xCounts(A_cell, path = paste(argv$output,"/10x",sep = ""),genome = argv$genome,overwrite = TRUE)
