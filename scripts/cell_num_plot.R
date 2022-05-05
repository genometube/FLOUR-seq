## Rscript to determine cell numbers. The result should be reviewed by the plot which clearly shows the steep decline. The algorithm is very simple. If the cell numbers should be determined more accurately, please refer the algorithm of 10x cellranger.
# usage:Rscript cell_num_plot.R path win slide
path <- commandArgs(T)[1]
win <- as.numeric(commandArgs(T)[2])
slide <- as.numeric(commandArgs(T)[3])
barcode_xls <- as.matrix(read.table(path,row.names = 1))
barcode_xls <- barcode_xls[-nrow(barcode_xls),]
for(i in seq(win,1,slide))
{
  front=floor(10^(max(log10(barcode_xls[,2]))*(i-win)))
  back=floor(10^(max(log10(barcode_xls[,2]))*i))
  if((log10(barcode_xls[,1])[front]-log10(barcode_xls[,1])[back])>(max(log10(barcode_xls[,1]))-min(log10(barcode_xls[,1])))*0.1) 
  {
    cat("inflection barcode rank:",back,"\n")
    break
  }
}
png(paste(path,".png",sep = ""))
plot(log10(barcode_xls[,2]),log10(barcode_xls[,1]),pch=19)
abline(v=log10(back), col='red')
abline(v=log10(front), col='blue')
legend("topright",legend = c("pre-inflection","inflection"),fill = c("blue","red"))
dev.off()
write(back,file=paste(path,".inflection_barcode_rank.txt",sep = ""))
