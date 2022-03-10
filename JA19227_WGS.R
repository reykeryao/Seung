library(stringr)
library(tidyr)
library(RColorBrewer)
dat<-read.delim("../JA19227_WGS/WGS_coverage.summary")
#dat<-separate(dat,Chr_contig,into = c("Contig","Number"),sep="_",remove = F)
#dat[is.na(dat)]<-0
#dat$Number<-as.integer(dat$Number)
#dat<-dat[order(dat$Number,dat$Start),]

#write.table(dat,"../JA19227_WGS/WGS_coverage.summary",quote=F,sep="\t",row.names=F)
pdf("../JA19227_WGS/coverage.pdf",height=8,width=11)
par(bty="n",mfrow=c(5,1),mar=c(1,2,1.5,0.5))
Bcol<-rgb(100,100,100,50, maxColorValue=255)
plot(dat$WT_depth,xlim=c(0,14000),ylim=c(0,250),type="l",col="green4",axes=F,xlab=NA,ylab=NA,main="WT")
row110<-rownames(dat)[dat$Chr_contig=="contig_110"]
rect(row110[1],0,tail(row110,n=1),250,col=Bcol,border=F)
axis(1,at=seq(0,14000,2000),labels = F)
axis(2,at=c(0,250),las=2,pos=-100)
lines(x = c(-100,14000), y = c(mean(dat$WT_depth),mean(dat$WT_depth)), lty = 3,col="red")

plot(dat$KO1_depth,xlim=c(0,14000),ylim=c(0,250),type="l",col="salmon",axes=F,xlab=NA,ylab=NA,main="KO1")
row110<-rownames(dat)[dat$Chr_contig=="contig_110"]
rect(row110[1],0,tail(row110,n=1),250,col=Bcol,border=F)
axis(1,at=seq(0,14000,2000),labels = F)
axis(2,at=c(0,250),las=2,pos=-100)
lines(x = c(-100,14000), y = c(mean(dat$KO1_depth),mean(dat$KO1_depth)), lty = 3,col="red")

plot(dat$KO3_depth,xlim=c(0,14000),ylim=c(0,250),type="l",col="darkorchid1",axes=F,xlab=NA,ylab=NA,main="KO3")
row110<-rownames(dat)[dat$Chr_contig=="contig_110"]
rect(row110[1],0,tail(row110,n=1),250,col=Bcol,border=F)
axis(1,at=seq(0,14000,2000),labels = F)
axis(2,at=c(0,250),las=2,pos=-100)
lines(x = c(-100,14000), y = c(mean(dat$KO3_depth),mean(dat$KO3_depth)), lty = 3,col="red")
text(grconvertX(0.01, "ndc", "user"), grconvertY(.75, "ndc", "user"), 
     srt=90,"WGS coverage (mean of 500-bp bin)", xpd=NA)

plot(dat$KO1_depth/dat$WT_depth,xlim=c(0,14000),ylim=c(0,5),type="l",col="salmon",axes=F,xlab=NA,ylab=NA,main="KO1 vs WT")
row110<-rownames(dat)[dat$Chr_contig=="contig_110"]
rect(row110[1],0,tail(row110,n=1),5,col=Bcol,border=F)
axis(1,at=seq(0,14000,2000),labels = F)
axis(2,at=c(0,1,5),las=2,pos=-100)
lines(x = c(-100,14000), y = c(1,1), lty = 3,col="red")
par(mar=c(5,2,1.5,0.5))
plot(dat$KO3_depth/dat$WT_depth,xlim=c(0,14000),ylim=c(0,5),type="l",
     col="darkorchid1",axes=F,xlab = "Genomic location (500-bp bins)",ylab=NA,main="KO3 vs WT")
row110<-rownames(dat)[dat$Chr_contig=="contig_110"]
rect(row110[1],0,tail(row110,n=1),5,col=Bcol,border=F)
axis(1,at=seq(0,14000,2000))
axis(2,at=c(0,1,5),las=2,pos=-100)
lines(x = c(-100,14000), y = c(1,1), lty = 3,col="red")
text(grconvertX(0.01, "ndc", "user"), grconvertY(.25, "ndc", "user"), 
     srt=90,"Fold change (WGS coverage)", xpd=NA)
dev.off()