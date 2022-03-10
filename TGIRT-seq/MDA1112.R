rm(list=ls())
gc()
setwd("~/Documents/NGS/Seung/Seung_git/")
library(goseq)
library(RColorBrewer)
library(tidyr)
library(DESeq2)
library(pheatmap)
library(Biobase)
library(grid)
load("../MDA_1112/MDA1112.image")


PAO1<-read.delim("../MDA_1112/combined_PAO1.counts")
PAO1<-separate(PAO1,col="Info",into = c("ID","Info"),sep = ":",extra="merge")
PAO1$Strand<-as.factor(PAO1$Strand)
PAO1[PAO1$Type=="CDS",6]<-"Protein coding"
PAO1[PAO1$Type=="pseudo",6]<-"Pseudogene"
PAO1<-PAO1[rowSums(PAO1[,9:24])>0,]
PAO1$Type<-factor(PAO1$Type,levels = c("Protein coding","Pseudogene","tmRNA","ncRNA","tRNA","rRNA"))
tRNA<-PAO1[PAO1$Type=="tRNA",]
PAO1<-PAO1[PAO1$Type!="tRNA",]
tRNA<-aggregate(.~Name, data=tRNA[,c(4,9:24)],FUN=sum)
tRNA$Type="tRNA"
tRNA$ID<-tRNA$Name
tRNA<-tRNA[,c(19,18,1:17)]
tRNA$Info<-tRNA$ID
AZPAE<-read.delim("../MDA_1112/combined_AZPAE12409.counts")
AZPAE$Strand<-as.factor(AZPAE$Strand)
AZPAE[AZPAE$Type=="gene",6]<-"CDS"
AZPAE[AZPAE$Type=="CDS",6]<-"Protein coding"
AZPAE[AZPAE$Type=="pseudo",6]<-"Pseudogene"
AZPAE<-AZPAE[rowSums(AZPAE[,7:22])>0,]
AZPAE[AZPAE$Type=="ncRNA",6]<-"RNase P"
AZPAE$Type<-factor(AZPAE$Type,levels = c("Protein coding","Pseudogene","RNase P","tRNA","rRNA"))
AZPAE<-AZPAE[AZPAE$Type!="tRNA",]

Combined<-read.delim("../MDA_1112/combined.counts")
Combined[Combined$Type=="gene",2]<-"CDS"
Combined[Combined$Type=="CDS",2]<-"Protein coding"
Combined[Combined$Type=="pseudo",2]<-"Pseudogene"
Combined[Combined$Name=="NQ18_RS09985",2]<-"RNaseP RNA"
Combined[Combined$Type=="ncRNA",2]<-"Other sncRNA"
Combined<-Combined[rowSums(Combined[,4:19])>0,]
Combined$Type<-factor(Combined$Type,levels = 
                        c("Protein coding","Pseudogene","Other sncRNA","tmRNA","tRNA","RNaseP RNA","rRNA"))
Combined<-Combined[Combined$Type!="tRNA",]
Combined$Info<-""
for (i in 1:dim(Combined)[1]){
  ID<-Combined[i,1]
  if (length(agrep("^NQ18_RS",ID,value=T))>0){
    Combined[i,20]<-AZPAE[AZPAE$Name==ID,23]
  } else if (length(agrep("^PA",ID,value=T))>0){
    Combined[i,20]<-PAO1[PAO1$ID==ID,8]
  }
}
Combined<-rbind(Combined,tRNA)
rownames(Combined)<-Combined$ID


#barchart is based on combined AZPAE and PAO1 counts
bcol<-brewer.pal(7, "Set1")
pdf("barplot.pdf",height=11,width=6)
par(mfrow=c(2,1))
agg<-aggregate(.~Type,Combined[,c(2,4:19)],FUN=sum)
agg1<-prop.table(as.matrix(agg[,2:17]),margin = 2)
row.names(agg1)<-agg$Type
barplot(100*as.matrix(agg1[,c(1:4,9:12,5:8,13:16)]),col=bcol,legend.text = F)
agg<-aggregate(.~Type,Combined[Combined$Type!="Protein coding" & Combined$Type!="Pseudogene" & Combined$Type!="rRNA"
                               ,c(2,4:19)],FUN=sum)
agg1<-prop.table(as.matrix(agg[,2:17]),margin = 2)
row.names(agg1)<-agg$Type
barplot(100*as.matrix(agg1[,c(1:4,9:12,5:8,13:16)]),col=bcol[3:6],legend.text = F)
dev.off()

#scatter plots
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(2^x-1, 2^y-1,method = "spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor.p <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(2^x-1, 2^y-1,method = "pearson"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,col="red",...)
}

tiff("scatter_all.tiff",res = 300,width = 4800,height=4800)
par(pch=".")
CPM<-apply(Combined[,4:19],MARGIN = 2,FUN=function(x){1e6*x/sum(x)})
CPM<-log2(CPM+1)
pairs(CPM,xlim=c(0,20),ylim=c(0,20),upper.panel=panel.cor,lower.panel = panel.line)
dev.off()
tiff("scatter_pro.tiff",res = 300,width = 4800,height=4800)
par(pch=".")
CPM<-apply(Combined[Combined$Type=="Protein coding",4:19],MARGIN = 2,FUN=function(x){1e6*x/sum(x)})
CPM<-log2(CPM+1)
pairs(CPM,xlim=c(0,20),ylim=c(0,20),upper.panel=panel.cor,lower.panel = panel.line)
dev.off()
tiff("scatter_snc.tiff",res = 300,width = 4800,height=4800)
par(pch=".")
CPM<-apply(Combined[Combined$Type!="Protein coding" & Combined$Type!="Pseudogene" & Combined$Type!="rRNA"
                    ,4:19],MARGIN = 2,FUN=function(x){1e6*x/sum(x)})
CPM<-log2(CPM+1)
pairs(CPM,xlim=c(0,20),ylim=c(0,20),upper.panel=panel.cor,lower.panel = panel.line)
dev.off()

#DESeq2 is based on combined AZPAE and PAO1 counts, but only PAO1 genes will have GO annotations
coldata<-data.frame("Strain"=rep(c("WT","KO"),each=8),
                    "Growth"=rep(rep(c("Log","Sta"),each=4),2),
                    "Sample"=colnames(Combined)[4:19])
coldata$Strain<-factor(coldata$Strain,levels=c("WT","KO"))
coldata$Growth<-factor(coldata$Growth,levels=c("Log","Sta"))
#mRNA DE
all_dds <- DESeqDataSetFromMatrix(countData = Combined[,4:19],
                                  colData = coldata,design = ~Strain+Growth+Strain:Growth)
pro_dds <- DESeqDataSetFromMatrix(countData = Combined[Combined$Type=="Protein coding",4:19],
                                   colData = coldata,design = ~Strain+Growth+Strain:Growth)
#snc_dds <- DESeqDataSetFromMatrix(countData = Combined[Combined$Type!="Protein coding" & Combined$Type!="Pseudogene" & Combined$Type!="rRNA",
#                                                       4:19],
#                                  colData = coldata,design = ~Strain+Growth+Strain:Growth)
all_dds <- DESeq(all_dds,parallel = T)
pro_dds <- DESeq(pro_dds,parallel = T)
#snc_dds <- DESeq(snc_dds,parallel = T)

plot_volcanoS<-function(dat_non_sig,dat_sig,x_lim,y_lim,xlab_seq,ylab_seq,title){
  plot(dat_non_sig$log2FoldChange, dat_non_sig$log10, xlim=x_lim,ylim=y_lim,
       main=title, xlab=NA,xaxt="n",bty="n", yaxt="n",
       ylab=NA,col="gray50",
       pch=20, cex=.7)
  axis(side = 1,at = xlab_seq,labels = NA,tck=-0.02,lwd=0.5)
  axis(side = 2,labels = NA,at = ylab_seq,las=1,tck=-0.02,lwd=0.5)
  abline(v=0, col="black", lty=3, lwd=0.5)
#  abline(v=-log2(1.5), col="black", lty=4, lwd=0.5)
#  abline(v=log2(1.5), col="black", lty=4, lwd=0.5)
  abline(h=-log10(0.05), col="black", lty=3, lwd=0.5)
  pro<-dat_sig[dat_sig$Type=="Protein coding",]
  points(pro$log2FoldChange,pro$log10, pch=21, 
         bg=bcol[pro$Type], cex=.75,lwd=0.1)
  non_pro<-dat_sig[dat_sig$Type!="Protein coding",]
  points(non_pro$log2FoldChange,non_pro$log10, pch=21, 
           bg=bcol[non_pro$Type], cex=.75,lwd=0.1)
  points(dat_non_sig$log2FoldChange[dat_non_sig$Type=="tmRNA"],
         dat_non_sig$log10[dat_non_sig$Type=="tmRNA"], pch=21, 
         bg=bcol[4], cex=.75,lwd=0.1)
}

plot_volcanoP<-function(dat_non_sig,dat_sig,x_lim,y_lim,xlab_seq,ylab_seq,title){
  plot(dat_non_sig$log2FoldChange, dat_non_sig$log10, xlim=x_lim,ylim=y_lim,
       main=title, xlab=NA,xaxt="n",bty="n", yaxt="n",
       ylab=NA,col="gray50",
       pch=20, cex=1)
  axis(side = 1,at = xlab_seq,labels = NA,tck=-0.02,lwd=0.5)
  axis(side = 2,labels = NA,at = ylab_seq,las=1,tck=-0.02,lwd=0.5)
  abline(v=0, col="black", lty=3, lwd=0.5)
#  abline(v=-log2(1.5), col="black", lty=4, lwd=0.5)
#  abline(v=log2(1.5), col="black", lty=4, lwd=0.5)
  abline(h=-log10(0.05), col="black", lty=3, lwd=0.5)
  points(dat_sig$log2FoldChange,dat_sig$log10, pch=21, 
         bg=bcol[1], cex=1.2,lwd=0.1)
}
#all genes
bcol<-brewer.pal(7, "Set1")
bcol<-paste0(bcol,"A0")

pdf("DESeq_all_genes.pdf",height=4,width=6)
par(mfrow=c(1,2),mar = c(1, 1, 2, 0.1))
#log phase KOvsWT
res<-data.frame(results(all_dds,contrast = c("Strain","KO","WT"),alpha = 0.05))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05 )
sig<-subset(res, padj<=0.05 )
sig<-sig[order(sig$log10,decreasing = T),]
plot_volcanoS(non_sig,sig,c(-10,10),c(0,30),seq(-10,10,5),seq(0,30,10),"KO vs WT (log phase)")
legend(-10,25,legend = c(levels(res$Type)[c(1,3,4,5)],"non-significant"),
       pt.bg=c(bcol[c(1,3,4,5)],"gray50"),cex=.5,
       pch=21,pt.lwd=0.1,bty="n")
Up<-sig[sig$log2FoldChange>0,]
UP_tRNA<-unique(substr(Up$Row.names[Up$Type=="tRNA"],1,6))
Up<-Up[1:5,]
Down<-sig[sig$log2FoldChange<0,]
Down<-Down[1:5,]
text(Up$log2FoldChange+1,Up$log10+3,labels = paste(Up$Name,Up$Info,sep=": "),cex=0.5,pos = 4)
text(Down$log2FoldChange-1,Down$log10+3,labels = paste(Down$Name,Down$Info,sep=": "),cex=0.5,pos = 2)
name_l<-c("rpoS","katE","dnaK","algU","recN","clpB","asrA","groES")
name_l<-res[res$Name%in%name_l,]
arrows(name_l$log2FoldChange,name_l$log10,name_l$log2FoldChange+1,name_l$log10+3,code=2,length=0.05)
text(name_l$log2FoldChange+1,name_l$log10+3,labels = name_l$Name,cex=0.5,pos = 4)

res<-data.frame(results(all_dds,list(c("Strain_KO_vs_WT","StrainKO.GrowthSta"))))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05)
sig<-sig[order(sig$log10,decreasing = T),]
plot_volcanoS(non_sig,sig,c(-10,10),c(0,90),seq(-10,10,5),seq(0,90,15),"KO vs WT (stationary phase)")
Up<-sig[sig$log2FoldChange>0,]
Up<-Up[1:5,]
Down<-sig[sig$log2FoldChange<0,]
Down<-Down[1:5,]
text(Up$log2FoldChange,Up$log10,labels = paste(Up$Name,Up$Info,sep=": "),cex=0.5,pos = 4)
text(Down$log2FoldChange,Down$log10,labels = paste(Down$Name,Down$Info,sep=": "),cex=0.5,pos = 2)
name_l<-c("rpoS","katA","ibpA","algU","recN","trxB2","gor","recA","htpX")
name_l<-res[res$Name%in%name_l,]
name_up<-name_l[name_l$log2FoldChange>0,]
arrows(name_up$log2FoldChange,name_up$log10,name_up$log2FoldChange+1,name_up$log10+3,code=2,length=0.05)
text(name_up$log2FoldChange+1,name_up$log10+3,labels = name_up$Name,cex=0.5,pos = 4)
name_down<-name_l[name_l$log2FoldChange<0,]
arrows(name_down$log2FoldChange,name_down$log10,name_down$log2FoldChange+1,name_down$log10+3,code=2,length=0.05)
text(name_down$log2FoldChange-1,name_down$log10+3,labels = name_down$Name,cex=0.5,pos = 4)
dev.off()

#sncRNA and tRNA volcano

pdf("DESeq_tRNA.pdf",height=4,width=6)
par(mfrow=c(1,2),mar = c(2, 2, 2, 0.1))
#log phase KOvsWT
res<-data.frame(results(all_dds,contrast = c("Strain","KO","WT"),alpha = 0.05))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
res<-res[complete.cases(res),]
res<-res[res$Type=="tRNA",]
non_sig<-subset(res, padj>0.05 )
sig<-subset(res, padj<=0.05 )
sig<-sig[order(sig$log10,decreasing = T),]
plot_volcanoS(non_sig,sig,c(-2,2),c(0,6),seq(-2,2,1),seq(0,6,2),"KO vs WT (log phase)")
axis(1,at=seq(-2,2,1),labels = seq(-2,2,1))
axis(2,las=2,at=seq(0,6,2),labels = seq(0,6,2))
legend(-2,5,legend = c("Significant tRNA","non-significant"),
       pt.bg=c(bcol[5],"gray50"),cex=0.5,
       pch=21,pt.lwd=0.1,bty="n")
text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Row.names,cex=0.5)

res<-data.frame(results(all_dds,list(c("Strain_KO_vs_WT","StrainKO.GrowthSta"))))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
res<-res[complete.cases(res),]
res<-res[res$Type=="tRNA",]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05)
sig<-sig[order(sig$log10,decreasing = T),]
plot_volcanoS(non_sig,sig,c(-3,3),c(0,4),seq(-3,3,1),seq(0,4,1 ),"KO vs WT (stationary phase)")
axis(1,at=seq(-3,3,1),labels = seq(-3,3,1))
axis(2,las=2,at=seq(0,4,1),labels = seq(0,4,1))
text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Row.names,cex=0.5)
dev.off()

pdf("DESeq_sncRNA.pdf",height=4,width=6)
par(mfrow=c(1,2),mar = c(2, 2, 2, 0.1))
#log phase KOvsWT
res<-data.frame(results(all_dds,contrast = c("Strain","KO","WT"),alpha = 0.05))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
res<-res[complete.cases(res),]
res<-res[res$Type=="Other sncRNA"| res$Type=="tmRNA"|res$Type=="RNaseP RNA",]
non_sig<-subset(res, padj>0.05 )
sig<-subset(res, padj<=0.05 )
sig<-sig[order(sig$log10,decreasing = T),]
plot_volcanoS(non_sig,sig,c(-2,2),c(0,9),seq(-2,2,1),seq(0,9,3),"KO vs WT (log phase)")
axis(1,at=seq(-2,2,1),labels = seq(-2,2,1))
axis(2,las=2,at=seq(0,9,3),labels = seq(0,9,3))
legend(-2,5,legend = c(sig$Info,"tmRNA","non-significant"),
       pt.bg=c(bcol[c(3,4)],"gray50"),cex=0.5,
       pch=21,pt.lwd=0.1,bty="n")

res<-data.frame(results(all_dds,list(c("Strain_KO_vs_WT","StrainKO.GrowthSta"))))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
res<-res[complete.cases(res),]
res<-res[res$Type=="Other sncRNA"| res$Type=="tmRNA"|res$Type=="RNaseP RNA",]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05)
sig<-sig[order(sig$log10,decreasing = T),]
plot_volcanoS(non_sig,sig,c(-5,5),c(0,6),seq(-5,5,1),seq(0,6,2),"KO vs WT (stationary phase)")
legend(-4,5,legend = c(sig$Info,"tmRNA","non-significant"),
       pt.bg=c(bcol[c(1,3,4)],"gray50"),cex=0.5,
       pch=21,pt.lwd=0.1,bty="n")
points(sig$log2FoldChange[1],sig$log10[1],col=bcol[1],pch=19)
axis(1,at=seq(-5,5,1),labels = seq(-5,5,1))
axis(2,las=2,at=seq(0,6,2),labels = seq(0,6,2))
dev.off()

#protein coding genes
#diff between KO and WT, in log phase
res<-data.frame(results(pro_dds,contrast = c("Strain","KO","WT"),alpha = 0.05))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
#res$Type<-Combined$Type
all_res<-list(res01=res)
#diff between KO and WT, in sta phase
res<-data.frame(results(pro_dds,list(c("Strain_KO_vs_WT","StrainKO.GrowthSta"))))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
all_res$res02<-res
#diff between Sta and log phase, in WT
res<-data.frame(results(pro_dds,contrast = c("Growth","Sta","Log")))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
all_res$res03<-res
#diff between Sta and log phase, in KO
res<-data.frame(results(pro_dds,list(c("Growth_Sta_vs_Log","StrainKO.GrowthSta"))))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
all_res$res04<-res
#effect between Sta and log phase that are different inKO and WT
res<-data.frame(results(pro_dds,name="StrainKO.GrowthSta"))
res$padj<-p.adjust(res$pvalue,"BH")
res$log10<- -log10(res$padj)
res<-merge(res,Combined[,c(2:3,20)],by=0)
all_res$res05<-res
#res01: diff between KO and WT, in log phase
#res02: diff between KO and WT, in sta phase
#res03: diff between Sta and log phase, in WT
#res04: diff between Sta and log phase, in KO
#res05: effect of growth condition (Sta vs Log) that are different between KO and WT
#cairo_ps("DE_ana.eps",height=9,width=6)
#par(mfrow=c(3,2),mar = c(1, 1, 2, 0.1))
res<-all_res[[1]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05 )
heatsig<-list(sig01=sig)
#plot_volcanoS(non_sig,sig,c(-10,10),c(0,30),seq(-10,10,5),seq(0,30,10),"KO vs WT (log phase)")
#legend(-10,25,legend = c(levels(res$Type),"non-significant"),
#       pt.bg=c(bcol,"gray50"),pch=21,cex=1.2,pt.lwd=0.1,bty="n")

res<-all_res[[2]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05 )
sig<-subset(res, padj<=0.05 )
#plot_volcanoS(non_sig,sig,c(-10,10),c(0,90),seq(-10,10,5),seq(0,90,15),"KO vs WT (stationary phase)")
heatsig$sig02<-sig

res<-all_res[[3]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05 )
sig<-subset(res, padj<=0.05 )
heatsig$sig03<-sig
sig$log10[sig$log10=="Inf" | sig$log10>90]<-90
#plot_volcanoS(non_sig,sig,c(-10,10),c(0,90),seq(-10,10,5),seq(0,90,15),"Stationary vs Log phase (WT)")

res<-all_res[[4]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05 )
sig<-subset(res, padj<=0.05 )
heatsig$sig04<-sig
sig$log10[sig$log10=="Inf" | sig$log10>90]<-90
#plot_volcanoS(non_sig,sig,c(-10,10),c(0,90),seq(-10,10,5),seq(0,90,15),"Stationary vs Log phase (KO)")

res<-all_res[[5]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05)
heatsig$sig05<-sig
#plot_volcanoS(non_sig,sig,c(-10,10),c(0,45),seq(-10,10,5),seq(0,45,15),"Stationary vs Log phase (diff in KO and WT)")
#dev.off()

pdf("DESeq_protein.pdf",height=4,width=6)
par(mfrow=c(1,2),mar = c(1, 1, 2, 0.1))
res<-all_res[[1]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05)
sig<-sig[order(sig$log10,decreasing = T),]
Up<-sig[sig$log2FoldChange>0,]
Up<-Up[1:5,]
Down<-sig[sig$log2FoldChange<0,]
Down<-Down[1:5,]
plot_volcanoP(non_sig,sig,c(-10,10),c(0,30),seq(-10,10,5),seq(0,30,10),"KO vs WT (log phase)")
legend(-10,25,legend = c("Significant","non-significant (padj > 0.05)"),cex=1,
       pt.bg=c(bcol[1],"gray50"),pch=21,pt.lwd=0.1,bty="n")
#text(sig$log2FoldChange,sig$log10,labels = paste(sig$Name,sig$Info,sep=": "),cex=0.5,pos=4)
text(Up$log2FoldChange,Up$log10,labels = paste(Up$Name,Up$Info,sep=": "),cex=0.5,pos = 4)
text(Down$log2FoldChange,Down$log10,labels = paste(Down$Name,Down$Info,sep=": "),cex=0.5,pos = 2)

res<-all_res[[2]]
res<-res[complete.cases(res),]
non_sig<-subset(res, padj>0.05)
sig<-subset(res, padj<=0.05)
sig<-sig[order(sig$log10,decreasing = T),]
Up<-sig[sig$log2FoldChange>0,]
Up<-Up[1:5,]
Down<-sig[sig$log2FoldChange<0,]
Down<-Down[1:5,]
plot_volcanoP(non_sig,sig,c(-10,10),c(0,90),seq(-10,10,5),seq(0,90,15),"KO vs WT (stationary phase)")
text(Up$log2FoldChange,Up$log10,labels = paste(Up$Name,Up$Info,sep=": "),cex=0.5,pos = 4)
text(Down$log2FoldChange,Down$log10,labels = paste(Down$Name,Down$Info,sep=": "),cex=0.5,pos = 2)
dev.off()



#GO
GO<-read.delim("gene_ontology_tab.txt")
colnames(GO)[1]<-"Gene"
PAO1_from_combined<-merge(Combined[grep("^PA",Combined$ID),],PAO1[,c(7,2:3)],by="ID",all=F)
PAO1_from_combined$length<-PAO1_from_combined$Ed-PAO1_from_combined$St
dat<-PAO1_from_combined[,c(1:3,23,20)]
dat$Info<-gsub("\"","",dat$Info)
rownames(dat)<-dat$ID
gene_cat<-GO[,c(1,5)]

dat$Sig<-0
res<-all_res[[1]]
sig_list<-subset(res,res$padj<=0.05 & res$log2FoldChange>0)
dat$Sig[dat$ID%in%sig_list$Row.names]<-1
pwf_sig<-dat$Sig
names(pwf_sig)<-dat$ID
pwf<-nullp(pwf_sig,bias.data = dat$length)
Go_ana<-data.frame(goseq(pwf,gene2cat = gene_cat,method ="Hypergeometric"))
Go_ana<-subset(Go_ana,Go_ana$ontology=="BP")
All_go<-list(Up1=Go_ana)
Over_list<-Go_ana$category[Go_ana$over_represented_pvalue<=0.05]
sig_list<-subset(res,res$padj<=0.05 & res$log2FoldChange<0)
dat$Sig<-0
dat$Sig[dat$ID%in%sig_list$Row.names]<-1
pwf_sig<-dat$Sig
names(pwf_sig)<-dat$ID
pwf<-nullp(pwf_sig,bias.data = dat$length)
Go_ana<-goseq(pwf,gene2cat = gene_cat,method ="Hypergeometric" )
Go_ana<-subset(Go_ana,Go_ana$ontology=="BP")
All_go$Down1<-Go_ana
Down_list<-Go_ana$category[Go_ana$over_represented_pvalue<=0.05]

dat$Sig<-0
res<-all_res[[2]]
sig_list<-subset(res,res$padj<=0.05 & res$log2FoldChange>0)
dat$Sig[dat$ID%in%sig_list$Row.names]<-1
pwf_sig<-dat$Sig
names(pwf_sig)<-dat$ID
pwf<-nullp(pwf_sig,bias.data = dat$length)
Go_ana<-goseq(pwf,gene2cat = gene_cat,method ="Hypergeometric" )
Go_ana<-subset(Go_ana,Go_ana$ontology=="BP")
All_go$Up2<-Go_ana
Over_list<-union(Over_list,Go_ana$category[Go_ana$over_represented_pvalue<=0.05])
sig_list<-subset(res,res$padj<=0.05 & res$log2FoldChange<0)
dat$Sig<-0
dat$Sig[dat$ID%in%sig_list$Row.names]<-1
pwf_sig<-dat$Sig
names(pwf_sig)<-dat$ID
pwf<-nullp(pwf_sig,bias.data = dat$length)
Go_ana<-goseq(pwf,gene2cat = gene_cat,method ="Hypergeometric" )
Go_ana<-subset(Go_ana,Go_ana$ontology=="BP")
All_go$Down2<-Go_ana
Down_list<-union(Down_list,Go_ana$category[Go_ana$over_represented_pvalue<=0.05])

Go_ana<-All_go$Up1
Up<-Go_ana[Go_ana$category%in%Over_list,c(1,6,2)]
Go_ana<-All_go$Up2
Up<-merge(Up,Go_ana[Go_ana$category%in%Over_list,c(1,6,2)],by=1:2)
colnames(Up)[3:4]<-c("Log","Sta")
Up[,3:4]<- -log10(Up[,3:4])

Go_ana<-All_go$Down1
Down<-Go_ana[Go_ana$category%in%Down_list,c(1,6,2)]
Go_ana<-All_go$Down2
Down<-merge(Down,Go_ana[Go_ana$category%in%Down_list,c(1,6,2)],by=1:2)
colnames(Down)[3:4]<-c("Log","Sta")
Down[,3:4]<-log10(Down[,3:4])
Up[Up$category=="GO:0006548",3]<-Down[Down$category=="GO:0006548",3]
Up[Up$category=="GO:0043093",3]<-Down[Down$category=="GO:0043093",3]
Down<-Down[Down$category!="GO:0006548",]
Down<-Down[Down$category!="GO:0043093",]
All<-rbind(Up,Down)
All[All$Sta<= -5,4]<- -5
breaksList = seq(-5, 5, by = 0.1)
s_col<-colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList))
pdf ("Go_ana.pdf")
pheatmap(scale(as.matrix(All[,c(3:4)])), cluster_rows = T,fontsize_row =5,color=s_col,
         cluster_cols = F,labels_row = All$term, show_colnames = F,breaks = breaksList)
dev.off()

#heatmap
# up or down
tmp<-heatsig$sig02
rownames(tmp)<-tmp$Row.names
counts_norm<-data.frame(counts(pro_dds,normalized=T))
#name_list_down<-tmp[tmp$log2FoldChange<= -log2(1.5),]
name_list_down<-tmp[tmp$log2FoldChange<0,]
name_list_down<-name_list_down[order(name_list_down$Name),]
name_list_down<-rownames(name_list_down)
#name_list_up<-tmp[tmp$log2FoldChange>= log2(1.5),]
name_list_up<-tmp[tmp$log2FoldChange>0,]
#name_list_up<-name_list_up[order(name_list_up$Name),]
name_list_up<-rownames(name_list_up)
gene<-rbind(counts_norm[name_list_down,],
            counts_norm[name_list_up,])
order_name<-rownames(gene)
gene<-gene[,c(1:4,9:12,5:8,13:16)]
gene$WT_15<-rowMedians(as.matrix(gene[,1:4]))
gene$KO_15<-rowMedians(as.matrix(gene[,5:8]))
gene$WT_30<-rowMedians(as.matrix(gene[,9:12]))
gene$KO_30<-rowMedians(as.matrix(gene[,13:16]))
#gene[,1:16]<-t(scale(t(as.matrix(gene[,1:16])),center = T,scale = T))
#gene[,17:20]<-t(scale(t(as.matrix(gene[,17:20])),center = T,scale = T))
gene<-merge(gene,Combined[,c(3,20)],by=0)
rownames(gene)<-gene$Row.names
gene<-gene[order_name,]
gene$color<-"black"
Repair<-GO[agrep("repair",GO$GO.Term),]
Repair<-Repair[Repair$GO.Term!="protein repair",]
#unique(Repair$GO.Term)
# [1] "DNA repair"                                              
# [2] "DNA synthesis involved in DNA repair"                    
# [3] "base-excision repair"                                    
# [4] "DNA dealkylation involved in DNA repair"                 
# [5] "nucleotide-excision repair"                              
# [6] "double-strand break repair via nonhomologous end joining"
# [7] "excinuclease repair complex"                             
# [8] "nucleotide-excision repair, DNA incision"                
# [9] "regulation of DNA repair"                                
#[10] "mismatch repair"                                         
#[11] "double-strand break repair"                              
#[12] "double-strand break repair via homologous recombination" 
Repair<-unique(Repair$Gene)
#unique(GO[agrep("oxidative stress",GO$GO.Term),6])
#[1] "response to oxidative stress"                       
#[2] "positive regulation of response to oxidative stress"
#[3] "cellular response to oxidative stress"              
#[4] "regulation of cellular response to oxidative stress"
Oxidative_stress<-unique(GO[agrep("oxidative stress",GO$GO.Term),1])
gene[gene$Row.names%in%Repair,24]<-"red"
gene[gene$Row.names%in%Oxidative_stress,24]<-"orange"
gene$Info<-gsub("\"","",gene$Info)
#row_anno<-data.frame("Type"=gene[,18])
pdf("heatmap_all.pdf",height=11,width=4)
p<-pheatmap(as.matrix(gene[,c(18:21)]),scale="row",cluster_rows = F,fontsize_row =1,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), show_colnames = T)
cols<-gene[order(match(rownames(gene), p$gtable$grobs[[3]]$label)), ]$color
name_col<-c("black","red","orange")
p$gtable$grobs[[3]]$gp = gpar(col = cols, fontsize = 1)
p
dev.off()

pdf("heatmap_all1.pdf",height=11,width=4)
p<-pheatmap(log10(as.matrix(gene[,c(18:21)])+1),cluster_rows = F,fontsize_row =1,
            cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), show_colnames = T)
cols<-gene[order(match(rownames(gene), p$gtable$grobs[[3]]$label)), ]$color
name_col<-c("black","red","orange")
p$gtable$grobs[[3]]$gp = gpar(col = cols, fontsize = 1)
p
dev.off()
#heatmap_all in log
tmp<-heatsig$sig01
rownames(tmp)<-tmp$Row.names
counts_norm<-data.frame(counts(pro_dds,normalized=T))
#name_list_down<-tmp[tmp$log2FoldChange<= -log2(1.5),]
name_list_down<-tmp[tmp$log2FoldChange<0,]
name_list_down<-name_list_down[order(name_list_down$Name),]
name_list_down<-rownames(name_list_down)
#name_list_up<-tmp[tmp$log2FoldChange>= log2(1.5),]
name_list_up<-tmp[tmp$log2FoldChange>0,]
#name_list_up<-name_list_up[order(name_list_up$Name),]
name_list_up<-rownames(name_list_up)
gene<-rbind(counts_norm[name_list_down,],
            counts_norm[name_list_up,])
order_name<-rownames(gene)
gene<-gene[,c(1:4,9:12,5:8,13:16)]
gene$WT_15<-rowMedians(as.matrix(gene[,1:4]))
gene$KO_15<-rowMedians(as.matrix(gene[,5:8]))
gene$WT_30<-rowMedians(as.matrix(gene[,9:12]))
gene$KO_30<-rowMedians(as.matrix(gene[,13:16]))
#gene[,1:16]<-t(scale(t(as.matrix(gene[,1:16])),center = T,scale = T))
#gene[,17:20]<-t(scale(t(as.matrix(gene[,17:20])),center = T,scale = T))
gene<-merge(gene,Combined[,c(3,20)],by=0)
rownames(gene)<-gene$Row.names
gene<-gene[order_name,]
gene$color<-"black"
Repair<-GO[agrep("repair",GO$GO.Term),]
Repair<-Repair[Repair$GO.Term!="protein repair",]
#unique(Repair$GO.Term)
# [1] "DNA repair"                                              
# [2] "DNA synthesis involved in DNA repair"                    
# [3] "base-excision repair"                                    
# [4] "DNA dealkylation involved in DNA repair"                 
# [5] "nucleotide-excision repair"                              
# [6] "double-strand break repair via nonhomologous end joining"
# [7] "excinuclease repair complex"                             
# [8] "nucleotide-excision repair, DNA incision"                
# [9] "regulation of DNA repair"                                
#[10] "mismatch repair"                                         
#[11] "double-strand break repair"                              
#[12] "double-strand break repair via homologous recombination" 
Repair<-unique(Repair$Gene)
#unique(GO[agrep("oxidative stress",GO$GO.Term),6])
#[1] "response to oxidative stress"                       
#[2] "positive regulation of response to oxidative stress"
#[3] "cellular response to oxidative stress"              
#[4] "regulation of cellular response to oxidative stress"
Oxidative_stress<-unique(GO[agrep("oxidative stress",GO$GO.Term),1])
gene[gene$Row.names%in%Repair,24]<-"red"
gene[gene$Row.names%in%Oxidative_stress,24]<-"orange"
gene$Info<-gsub("\"","",gene$Info)
#row_anno<-data.frame("Type"=gene[,18])
pdf("heatmap_all_log.pdf",height=11,width=4)
p<-pheatmap(as.matrix(gene[,c(18:21)]),scale="row",cluster_rows = F,fontsize_row =1,
            cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), show_colnames = T)
cols<-gene[order(match(rownames(gene), p$gtable$grobs[[3]]$label)), ]$color
name_col<-c("black","red","orange")
p$gtable$grobs[[3]]$gp = gpar(col = cols, fontsize = 1)
p
dev.off()

pdf("heatmap_all1.pdf",height=11,width=4)
p<-pheatmap(log10(as.matrix(gene[,c(18:21)])+1),cluster_rows = F,fontsize_row =1,
            cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), show_colnames = T)
cols<-gene[order(match(rownames(gene), p$gtable$grobs[[3]]$label)), ]$color
name_col<-c("black","red","orange")
p$gtable$grobs[[3]]$gp = gpar(col = cols, fontsize = 1)
p
dev.off()

'''
lag_list<-c("PA3007","PA3622")
Mod_O_list<-intersect(Oxidative_stress,order_name)
Mod_O_list<-c(Mod_O_list,"PA0140","PA0848")
Mod_O_list<-Mod_O_list[c(1:3,5:8)]
gene<-rbind(counts_norm[rownames(counts_norm)=="NQ18_RS30375",],
            counts_norm[rownames(counts_norm)%in%lag_list,],
            counts_norm[rownames(counts_norm)%in%intersect(Repair,order_name),],
            counts_norm[rownames(counts_norm)%in%Mod_O_list,])
gene<-merge(gene,Combined[,c(3,20)],by=0)
gene<-gene[c(1,8,11,4,10,13,15,9,6,2,3,5,12,14,7),]
rownames(gene)<-gene$Row.names
gene$Name[1]<-"G2L4 RT"
gene$BP<-"DNA repair"
gene$BP[1]<-"RT"
gene$BP[9:15]<-"Oxidative stress"
gene$BP[2:3]<-"Stationary marker"

gene$WT_15<-rowMedians(as.matrix(gene[,2:5]))
gene$KO_15<-rowMedians(as.matrix(gene[,10:13]))
gene$WT_30<-rowMedians(as.matrix(gene[,6:9]))
gene$KO_30<-rowMedians(as.matrix(gene[,14:17]))
#gene[,2:17]<-t(scale(t(as.matrix(gene[,2:17])),center = T,scale = T))
#gene[,21:24]<-t(scale(t(as.matrix(gene[,21:24])),center = T,scale = T))
row_anno<-data.frame("BP"=gene[,20])
rownames(row_anno)<-gene$Name

pdf("heatmap.pdf")
pheatmap(log10(as.matrix(gene[,c(21:24)])+1), annotation_col = NA,cluster_rows = F,fontsize_row =5,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), 
         show_colnames = T)
dev.off()

#scaled -heatmap
gene[,2:17]<-t(scale(t(as.matrix(gene[,2:17])),center = T,scale = T))
gene[,21:24]<-t(scale(t(as.matrix(gene[,21:24])),center = T,scale = T))

pdf("heatmap1.pdf")
pheatmap(as.matrix(gene[,c(21:24)]), annotation_col = NA,cluster_rows = F,fontsize_row =5,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), 
         show_colnames = T)
dev.off()

'''
#boxplot
plot_box<-function(dat_frame,name,type){
  pcol=c("black","red","black","red")
  boxplot(dat_frame$Counts~dat_frame$Group,cex=0.5,
          main=paste0(name," (",type,")"),
          frame=F,outline=F,xlab=NA,ylab="Normalizied counts")
  stripchart(dat_frame$Counts~dat_frame$Group,vertical=T,method="jitter",
             jitter=0.15,add=T,col=pcol,pch=16,cex=1.5)
}

'''
name_list<-gene$Name
ID_list<-as.character(gene$Row.names)
comp_dat<-all_res[[2]]
comp_dat1<-all_res[[1]]
pdf("box.pdf",height=9,width=15)
par(mfrow=c(3,5))
pcol=c("black","red","black","red")
for (i in 1:15){
  dat_tmp<-data.frame((t(counts_norm[ID_list[i],1:16])))
  colnames(dat_tmp)<-"Counts"
  dat_tmp$Group<-unlist(paste(coldata$Strain,coldata$Growth,sep="_"))
  dat_tmp$Group<-factor(dat_tmp$Group,levels = c("WT_Log","KO_Log","WT_Sta","KO_Sta"))
  digits<-round(log10(max(dat_tmp$Counts)))
  first_num<-ceiling(max(dat_tmp$Counts)/10^digits)
  ymax<-3
  dat_tmp$Counts<-dat_tmp$Counts/10^digits
  boxplot(dat_tmp$Counts~dat_tmp$Group,cex=0.5,,ylim=c(0,ymax),
          main=name_list[i],frame=F,outline=F,xlab=NA,
          ylab=paste0("Normalizied counts (10^",digits,")"))
  stripchart(dat_tmp$Counts~dat_tmp$Group,vertical=T,method="jitter",
             jitter=0.15,add=T,col=pcol,pch=16,cex=1.5)

  sig_value<-comp_dat$padj[comp_dat$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"*")
  }
  sig_value<-comp_dat1$padj[comp_dat1$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"*")
  }
}
dev.off()
'''
#SOS heat
SOS_list<-c("lexA","recA","recN")
stress_list<-c("algU","rpoS")
Oxi_list<-c("ahpB","ahpC","ahpF","fumC1","gor","katA","katB","trxB2")
SOS_related_list<-c("cyaB","dps","htpX","ibpA","mexR","oprG","pslB","rhlR","tse5","recD")
counts_norm<-data.frame(counts(pro_dds,normalized=T))
counts_norm<-merge(counts_norm,Combined[,c(3,20)],by=0)
gene<-rbind(counts_norm[counts_norm$Name=="NQ18_RS30375",],
            counts_norm[counts_norm$Name%in%Oxi_list,],
            counts_norm[counts_norm$Name%in%SOS_list,],
            counts_norm[counts_norm$Name%in%SOS_related_list,],
            counts_norm[counts_norm$Name%in%stress_list,])
gene<-gene[c(1:2,4,3,8,6,7,9,5,10:12,19,14,17,18,13,21,15,20,16,22:24),]
rownames(gene)<-gene$Row.names
gene$Name[1]<-"G2L4 RT"
gene$BP<-c("",rep("Oxidative stress",8),rep("SOS response",3),
           rep("SOS related global response",10),rep("Stress response sigma factor",2))

gene$WT_15<-rowMedians(as.matrix(gene[,2:5]))
gene$KO_15<-rowMedians(as.matrix(gene[,10:13]))
gene$WT_30<-rowMedians(as.matrix(gene[,6:9]))
gene$KO_30<-rowMedians(as.matrix(gene[,14:17]))
#heatmap1
row_anno<-data.frame("BP"=gene[,20])
rownames(row_anno)<-gene$Name
pdf("heatmap1.pdf")
pheatmap(log10(as.matrix(gene[,c(21:24)])+1), annotation_col = NA,cluster_rows = F,fontsize_row =5,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), 
         show_colnames = T)
dev.off()
#gene[,2:17]<-t(scale(t(as.matrix(gene[,2:17])),center = T,scale = T))
#gene[,21:24]<-t(scale(t(as.matrix(gene[,21:24])),center = T,scale = T))
pdf("heatmap2.pdf")
pheatmap(as.matrix(gene[,c(21:24)]),scale="row", annotation_col = NA,cluster_rows = F,fontsize_row =5,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), 
         show_colnames = T)
dev.off()
#boxplot
name_list<-gene$Name
ID_list<-as.character(gene$Row.names)
comp_dat<-all_res[[2]]
comp_dat1<-all_res[[1]]
pdf("box1.pdf",height=9,width=15)
par(mfrow=c(3,5))
pcol=c("black","red","black","red")
for (i in 1:24){
  dat_tmp<-data.frame((t(gene[ID_list[i],2:17])))
  colnames(dat_tmp)<-"Counts"
  dat_tmp$Group<-unlist(paste(coldata$Strain,coldata$Growth,sep="_"))
  dat_tmp$Group<-factor(dat_tmp$Group,levels = c("WT_Log","KO_Log","WT_Sta","KO_Sta"))
  digits<-round(log10(max(dat_tmp$Counts)))
  first_num<-ceiling(max(dat_tmp$Counts)/10^digits)
  ymax<-3
  dat_tmp$Counts<-dat_tmp$Counts/10^digits
  boxplot(dat_tmp$Counts~dat_tmp$Group,cex=0.5,,ylim=c(0,ymax),
          main=name_list[i],frame=F,outline=F,xlab=NA,
          ylab=paste0("Normalizied counts (10^",digits,")"))
  stripchart(dat_tmp$Counts~dat_tmp$Group,vertical=T,method="jitter",
             jitter=0.15,add=T,col=pcol,pch=16,cex=1.5)
  
  sig_value<-comp_dat$padj[comp_dat$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"*")
  }
  sig_value<-comp_dat1$padj[comp_dat1$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"*")
  }
}
dev.off()
#newheatmap

SOS_list<-c(paste0("PA",c(3413,2288,3414,3616,1044,1045,"0069","0922","0669","0670","0671")),
            "sulA","lexA","recA","recN","muxB","hupB")
counts_norm<-data.frame(counts(pro_dds,normalized=T))
counts_norm<-merge(counts_norm,Combined[,c(3,20)],by=0)
gene<-rbind(counts_norm[counts_norm$Name=="NQ18_RS30375",],
            counts_norm[counts_norm$Name%in%SOS_list,])
gene<-gene[c(1,9,12,11,2:8,10,14:16,13,17,18),]
rownames(gene)<-gene$Row.names
gene$Name[1]<-"G2L4 RT"
gene$WT_15<-rowMedians(as.matrix(gene[,2:5]))
gene$KO_15<-rowMedians(as.matrix(gene[,10:13]))
gene$WT_30<-rowMedians(as.matrix(gene[,6:9]))
gene$KO_30<-rowMedians(as.matrix(gene[,14:17]))
#heatmap1
pdf("heatmap3.pdf")
pheatmap(log10(as.matrix(gene[,c(20:23)])+1), annotation_col = NA,cluster_rows = F,fontsize_row =5,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), 
         show_colnames = T)
dev.off()
#gene[,2:17]<-t(scale(t(as.matrix(gene[,2:17])),center = T,scale = T))
#gene[,21:24]<-t(scale(t(as.matrix(gene[,21:24])),center = T,scale = T))
pdf("heatmap4.pdf")
pheatmap(as.matrix(gene[,c(20:23)]),scale="row", annotation_col = NA,cluster_rows = F,fontsize_row =5,
         cluster_cols = F,labels_row = paste(gene$Name,gene$Info,sep=": "), 
         show_colnames = T)
dev.off()
#bocplot
name_list<-gene$Name
ID_list<-as.character(gene$Row.names)
comp_dat<-all_res[[2]]
comp_dat1<-all_res[[1]]
pdf("box2.pdf",height=9,width=15)
par(mfrow=c(3,5))
pcol=c("black","red","black","red")
for (i in 2:18){
  dat_tmp<-data.frame((t(gene[ID_list[i],2:17])))
  colnames(dat_tmp)<-"Counts"
  dat_tmp$Group<-unlist(paste(coldata$Strain,coldata$Growth,sep="_"))
  dat_tmp$Group<-factor(dat_tmp$Group,levels = c("WT_Log","KO_Log","WT_Sta","KO_Sta"))
  digits<-round(log10(max(dat_tmp$Counts)))
  first_num<-ceiling(max(dat_tmp$Counts)/10^digits)
  ymax<-3
  dat_tmp$Counts<-dat_tmp$Counts/10^digits
  boxplot(dat_tmp$Counts~dat_tmp$Group,cex=0.5,,ylim=c(0,ymax),
          main=name_list[i],frame=F,outline=F,xlab=NA,
          ylab=paste0("Normalizied counts (10^",digits,")"))
  stripchart(dat_tmp$Counts~dat_tmp$Group,vertical=T,method="jitter",
             jitter=0.15,add=T,col=pcol,pch=16,cex=1.5)
  
  sig_value<-comp_dat$padj[comp_dat$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"*")
  }
  sig_value<-comp_dat1$padj[comp_dat1$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"*")
  }
}
dev.off()

#new box of log pagse

log_list<-c("asrA","clpB","groES","dnaK")
counts_norm<-data.frame(counts(pro_dds,normalized=T))
counts_norm<-merge(counts_norm,Combined[,c(3,20)],by=0)
gene<-counts_norm[counts_norm$Name%in%log_list,]
rownames(gene)<-gene$Row.names
gene$WT_15<-rowMedians(as.matrix(gene[,2:5]))
gene$KO_15<-rowMedians(as.matrix(gene[,10:13]))
gene$WT_30<-rowMedians(as.matrix(gene[,6:9]))
gene$KO_30<-rowMedians(as.matrix(gene[,14:17]))

name_list<-gene$Name
ID_list<-as.character(gene$Row.names)
comp_dat<-all_res[[2]]
comp_dat1<-all_res[[1]]
pdf("box3.pdf",height=9,width=15)
par(mfrow=c(3,5))
pcol=c("black","red","black","red")
for (i in 1:4){
  dat_tmp<-data.frame((t(gene[ID_list[i],2:17])))
  colnames(dat_tmp)<-"Counts"
  dat_tmp$Group<-unlist(paste(coldata$Strain,coldata$Growth,sep="_"))
  dat_tmp$Group<-factor(dat_tmp$Group,levels = c("WT_Log","KO_Log","WT_Sta","KO_Sta"))
  digits<-round(log10(max(dat_tmp$Counts)))
  first_num<-ceiling(max(dat_tmp$Counts)/10^digits)
  ymax<-3
  dat_tmp$Counts<-dat_tmp$Counts/10^digits
  boxplot(dat_tmp$Counts~dat_tmp$Group,cex=0.5,,ylim=c(0,ymax),
          main=name_list[i],frame=F,outline=F,xlab=NA,
          ylab=paste0("Normalizied counts (10^",digits,")"))
  stripchart(dat_tmp$Counts~dat_tmp$Group,vertical=T,method="jitter",
             jitter=0.15,add=T,col=pcol,pch=16,cex=1.5)
  
  sig_value<-comp_dat$padj[comp_dat$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=3,x1=4,y0=2.5)
    text(3.5,2.6,"*")
  }
  sig_value<-comp_dat1$padj[comp_dat1$Row.names==ID_list[i]]
  if(sig_value<=0.01){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"**")
  }else if (sig_value<=0.05){
    segments(x0=1,x1=2,y0=2.5)
    text(1.5,2.6,"*")
  }
}
dev.off()
