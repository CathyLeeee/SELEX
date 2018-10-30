#!/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/hist.R
sample25<-read.table("/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/Sample_LHT-ZM-I25.reads_fq.barplot.txt",sep="\t")
#pdf("/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/sample25.pdf")
#barplot(log(sample25[,2]+1),xlab="frequency of reads",ylab="log(reads number+1)",main="Distribution of reads frequecny")
#dev.off()

sample25<-sample25[1:100,]
pdf("/Share/home/zhangqw/SELECT/20161220/Sample_LHT-ZM-I25/sample25.zoomin.pdf")
barplot(log(sample25[,2]+1),xlab="frequency of reads",ylab="log(reads number+1)",main="Distribution of reads frequecny")
dev.off()
