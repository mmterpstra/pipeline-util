#!/usr/bin/env Rscript
#options(echo=TRUE) # if you want see commands in output file

warning("args=datatable.txt Sample Sample");



args <- commandArgs(trailingOnly = TRUE)
#print(args)
tsv  <- args[1]
sample1<-args[2]
sample2<-args[3]

#test like Rscript pipeline-util/bin/PlotSegsSampleVsSample.R /scratch/umcg-mterpstra/projects/WES_novogene_path_P8//varscan.T11_20692//multi/merged.txt T11_20692 T11_23706
#tsv    <-"/scratch/umcg-mterpstra/projects/WES_novogene_path_P8//varscan.T11_20692//multi/merged.txt"
#sample1<-"T11_20692"
#sample2<-"T11_23706"

seg.data.frame<-read.table(tsv, sep="\t", header=TRUE)

base <- sub(pattern=".txt", replacement="", tsv, fixed=TRUE)

chromtags.data.frame<-read.table(paste(base,"chrlen.txt", sep="."), sep="\t", header=TRUE)


#print(paste(base,sample1,sample2,"pdf", sep="."))

row.names(chromtags.data.frame)<-chromtags.data.frame$chrom
chromtags.data.frame<-chromtags.data.frame[with(chromtags.data.frame, order(offset)), ]
#headchromtags.data.frame<-head(chromtags.data.frame)
#print(headchromtags.data.frame)


pdf(file=paste(base,sample1,sample2,"pdf", sep="."), height = 7, width = 28)

pairwise.table.tsv <-  data.frame(seg.data.frame[ ,c(paste(sample1,'ID', sep="."),
                                                     'chrom',
                                                     'start',
                                                     'end',
                                                     paste(sample1,"num.mark", sep="."),
                                                     paste(sample1,"seg.mean", sep="."),
                                                     paste(sample1,"num.mark", sep="."),
                                                     paste(sample1,"seg.mean", sep="."),
                                                     paste(sample2,"num.mark", sep="."),
                                                     paste(sample2,"seg.mean", sep=".")) ],
                                   (seg.data.frame[ , paste(sample1,"seg.mean", sep=".")] - seg.data.frame[ , paste(sample2,"seg.mean", sep=".")]) )


colnames(pairwise.table.tsv)<-c('ID',
                                'chrom',
                                'start',
                                'end',
                                'num.mark',
                                'seg.mean',
                                paste(sample1,"num.mark", sep="."),
       	       	       	       	paste(sample1,"seg.mean", sep="."),
       	       	       	       	paste(sample2,"num.mark", sep="."),
       	       	       	        paste(sample2,"seg.mean", sep="."),
                                paste('log2( ',sample1,' / normal ) - log2( ',sample2,' / normal )', sep=""))


write.table(pairwise.table.tsv, sep="\t", file=paste(base,sample1,sample2,"tsv", sep="."), row.names=FALSE, quote=FALSE)
write.table(pairwise.table.tsv[,c('ID', 'chrom', 'start', 'end', 'num.mark', 'seg.mean')], sep="\t", file=paste(base,sample1,"seg", sep="."), row.names=FALSE, quote=FALSE)


plot(seg.data.frame$start.offset, (seg.data.frame[ , paste(sample1,"seg.mean", sep=".")] - seg.data.frame[ , paste(sample2,"seg.mean", sep=".")]), pch=".", ylim=c(-2,2), xaxt = "n", xlab = "Genome position (Chromosomes and ticks at 25Mbp)", ylab = paste("log2(",sample1,") - log2(",sample2,")", sep = ' '), main = "Tumor vs tumor CNV plot")
abline(a=log(x=2,base=2),b=0, col="red")
abline(a=log(x=1.5,base=2),b=0, col="orange")
abline(a=log(x=1,base=2),b=0, col="green")
abline(a=log(x=0.5,base=2),b=0, col="blue")
abline(a=log(x=0.25,base=2),b=0, col="black")

legend("bottomright", title="copynumber relative to normal",c("4n","3n","2n","1n","0n","detected"), fill=c("red","orange","green","blue","black","grey"), horiz=TRUE)

segments(x0 = seg.data.frame$start.offset, x1 = seg.data.frame$end.offset, y0 = seg.data.frame[ , paste(sample1,"seg.mean", sep=".")] - seg.data.frame[ , paste(sample2,"seg.mean", sep=".")], y1 = seg.data.frame[ , paste(sample1,"seg.mean", sep=".")] - seg.data.frame[ , paste(sample2,"seg.mean", sep=".")])

#create minor ticks
#note list specify without the c() construct
axis(1, at=chromtags.data.frame$offset, labels=chromtags.data.frame$chrom, las=2)
for (chrI in 1:length(chromtags.data.frame$chrom)) {
  tickPosIncrement=25000000
  tickPos=tickPosIncrement
  while (chromtags.data.frame$length[chrI] > tickPos) {
    axis(1, at=c(chromtags.data.frame$offset[chrI] + tickPos), labels=c(""), las=2 ,tck=-0.008)
    tickPos = tickPos + tickPosIncrement
  }
}

devcurrentmsg<-dev.off()
