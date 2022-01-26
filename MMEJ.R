rm(list=ls())
gc(verbose = F)
setwd("~/Documents/NGS/Seung/MDA_1128/")
library(dplyr)
library(RColorBrewer)
name1<-c("10-MgMn-GCA2H","1-Mg-GCNoRT","2-Mg-GCWT","3-Mg-GCA2V","4-Mg-GCA2R",
         "6-MgMn-GCNoRT","7-MgMn-GCWT","8-MgMn-GCA2V","Mg-AA-G2DA","Mg-AA-G2IA",
         "Mg-AA-G2WT","Mg-AA-GsAI","Mg-AA-GsDA","Mg-AA-GsWT","Mg-GG-G2DA",
         "Mg-GG-G2IA","Mg-GG-G2WT","Mg-GG-GsAI","Mg-GG-GsDA","Mg-GG-GsWT",
         "Mg-GG-NoRT","Mn-AA-G2DA","Mn-AA-G2IA","Mn-AA-G2WT","Mn-AA-GsAI",
         "Mn-AA-GsDA","Mn-AA-GsWT","Mn-AA-NoRT","Mn-GG-G2DA","Mn-GG-G2IA",
         "Mn-GG-G2WT","Mn-GG-GsAI","Mn-GG-GsDA","Mn-GG-GsWT","Mn-GG-NoRT")
#function to make RC sequences
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"        
    if(xx[bbb]=="C") y[bbb]<-"G"        
    if(xx[bbb]=="G") y[bbb]<-"C"        
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)    
}

#Core has to be at the 5' end, not in the middle
Core <- "^GTAAGAGCCTACTCATGGATCCTCCTTGTGATGTAAGGT"
#Core RC has to be at the 3' end, not in the middle
CoreRC <-"ACCTTACATCACAAGGAGGATCCATGAGTAGGCTCTTAC$"
#template has to be at the 5' end, not in the middle
Temp<-c(rep("GTAAGAGCCTACTCATGGATCCTCCTTGTGATGTAAGGTCCCGG",8),
        rep("GTAAGAGCCTACTCATGGATCCTCCTTGTGATGTAAGGTCCCAA",6),
        rep("GTAAGAGCCTACTCATGGATCCTCCTTGTGATGTAAGGTCCCGG",7),
        rep("GTAAGAGCCTACTCATGGATCCTCCTTGTGATGTAAGGTCCCAA",7),
        rep("GTAAGAGCCTACTCATGGATCCTCCTTGTGATGTAAGGTCCCGG",7))
TempRC<-unlist(lapply(Temp,rev.comp))
Temp<-paste0("^",Temp)
TempRC<-paste0(TempRC,"$")

for (j in c(1:35)){
  dat<-read.table(paste0("RAW_seq/",name1[j],".seq"),col.names="Seq")
  print (paste(name1[j],dim(dat)[1]))
  dat<-data.frame(table(dat$Seq))
  colnames(dat)<-c("Seq","Frq")
  print (sum(dat$Frq))
  dat$Seq<-as.character(dat$Seq)
  dat$Len<-nchar(dat$Seq)
  dat$Core_match<-0
  dat$CoreRC_match<-0
  dat$Temp_match<-0
  dat$TempRC_match<-0
  for (i in 1:dim(dat)[1]){
    Seq<-dat[i,1]
    Core_dist<-CoreRC_dist<-Temp_dist<-TempRC_dist<-dist1<-dist2<-Inf
    Pos<-aregexec(Core,Seq,max=list(all = 5 ,del = 3, ins = 5, sub = 5),fixed=F)
    if (Pos[[1]][1]> 0){
      dat[i,4]<-1
      Core_dist<-adist(Core,unlist(regmatches(Seq, Pos)),fixed = F)[1,1]
    }
    Pos<-aregexec(CoreRC,Seq,max=list(all = 5 ,del = 3, ins = 5, sub = 5),fixed=F)
    if (Pos[[1]][1]>0){
      dat[i,5]<-1
      CoreRC_dist<-adist(CoreRC,unlist(regmatches(Seq, Pos)),fixed = F)[1,1]
    }
    Pos<-aregexec(Temp[j],Seq,max=list(all = 5 ,del = 3, ins = 5, sub = 5),fixed=F)
    if (Pos[[1]][1]>0){
      Temp_dist<-adist(Temp[j],unlist(regmatches(Seq, Pos)),fixed = F)[1,1]
      dist1<-Temp_dist-Core_dist
      if ((dist1<=3) & (Temp_dist<4)){
        dat[i,6]<-1
      }
    }
    Pos<-aregexec(TempRC[j],Seq,max=list(all = 5 ,del = 3, ins = 5, sub = 5),fixed=F)
    if (Pos[[1]][1]>0){
      TempRC_dist<-adist(TempRC[j],unlist(regmatches(Seq, Pos)),fixed = F)[1,1]
      dist2<-TempRC_dist-CoreRC_dist
      if ((dist2<=3) & (TempRC_dist<4)) {
        dat[i,7]<-1
      }
    }
    #RC flag 1
    Flag_for_RC1<- (dat$TempRC_match[i] == 1) & 
      (((dat$Temp_match[i] == 1) & ((dist2<dist1) | (TempRC_dist < Temp_dist))) | (dat$Temp_match[i]==0))
    Flag_for_RC2<- ((dat$TempRC_match[i] == 0) & (dat$CoreRC_match[i] == 1) &
      ((CoreRC_dist<Core_dist) | (dat$Core_match[i]==0)))
    Flag_for_RC<- (Flag_for_RC1 | Flag_for_RC2)
    
    if (Flag_for_RC) {
        dat$Seq[i]<-rev.comp(Seq)
        tmp<-dat$TempRC_match[i]
        dat$TempRC_match[i]<-dat$Temp_match[i] 
        dat$Temp_match[i]<-tmp
        tmp<-dat$CoreRC_match[i]
        dat$CoreRC_match[i]<-dat$Core_match[i] 
        dat$Core_match[i]<-tmp
    }
  }
  dat$Match<-rowSums(dat[,4:7])
  dat$Type<-case_when(
    #the type will be catogorized as
    #  Core_match CoreRC_match Temp_match TempRC_match Match	Type
    #           1            1          1            0   3	5-6   Terminal trasnfer/snap (2 subtype)
    #           1            0          1            0   2	5-6   Terminal trasnfer/snap (2 subtype)
    #           0            1          0            0   1	1     Other
    #           1            0          0            0   1	1     Other
    #           1            1          0            0   2	1     Other
    #           1            1          1            1   4	2-4   MMEJ (3 subtype)
    #           0            0          0            0   0	1     Other
    #none of the end matches core and coreRC, it's some product should be discarded, 
    dat$Match < 2 ~ 1,
    #Type 2-4 MMEJ, no trim back, depends on the length, suggest no, 1 , or multiple snap/jump
    dat$Match == 4 & dat$Len < 90 ~ 2,
    dat$Match == 4 & dat$Len < 120 ~ 3,
    dat$Match == 4 & dat$Len >= 120 ~ 4,
    #Type 5-6 snap/jump w/o or with terminal transfer
    dat$Match == 3 & dat$Len < 60 ~ 5,
    dat$Match == 3 & dat$Len >= 60 ~ 6,
    #match= 2 
    dat$Match == 2 & dat$Temp_match == 0 ~ 1,
    dat$Match == 2 & dat$Temp_match == 1 & dat$Len < 60 ~ 5,
    dat$Match == 2 & dat$Temp_match == 1 & dat$Len >= 60 ~ 6,
  )
  dat<-data.frame(dat%>%group_by(Seq)%>%summarize(Frq=sum(Frq),Type=unique(Type)))
  dat<-dat[order(dat$Frq,decreasing=T),]
  dat$ID<-paste0("ID",1:dim(dat)[1])
  print (sum(dat$Frq))
  write.table(dat,paste0(name1[j],".info"),quote=F,sep="\t",row.names=F)
  if (j==1){
    Type<-data.frame("Type"=1:6)
  }
  Type<-merge(Type,data.frame(aggregate(Frq~Type,data=dat,FUN=sum)),by=1,all=T)
  colnames(Type)[j+1]<-name1[j]
}
Type[is.na(Type)]<-0
write.table(Type,"MMEJ_by_type.txt",quote=F,sep="\t",row.names=F)
#check reads length distribution by Type, Type1 jump w or w/o terminal transferase, 
#Type2 MMEJ w or w/o terminal transferase
pdf("length_density.pdf",height=15,width=12,onefile = T)
par(mfrow=c(5,4))
for (i in 1:2){
  for (j in c(1:35)){
    dat<-read.delim(paste0(name1[j],".info"))
    if (i==1){
      dat<-dat[dat$Type<=4 & dat$Type>1,]
      
    } else {
      dat<-dat[dat$Type>4,]
    }
    plot(density(rep(nchar(dat$Seq),dat$Frq)),main = paste0(name1[j]," (Type ",i,")"),xlab=NA,ylab=NA,xlim=c(40,160))
  }
}
dev.off()


pdf("length_density1.pdf",height=15,width=12,onefile = T)
par(mfrow=c(5,4))
for (j in c(1:35)){
    dat<-read.delim(paste0(name1[j],".info"))
    dat<-dat[dat$Type>1,]
     
    plot(density(rep(nchar(dat$Seq),dat$Frq)),main = paste0(name1[j]),xlab=NA,ylab=NA,xlim=c(40,160))
}
dev.off()


pdf("Frq_by_type.pdf",width=11,height=8)
par(mar=c(10,3,1,8))
Type<-read.delim("MMEJ_by_type.txt",row.names =1,check.names = FALSE)
Type_frq<-100*prop.table(as.matrix(Type),margin = 2)
bcol=c("black",brewer.pal(9,"Reds")[c(6,4,2)],brewer.pal(9,"Blues")[c(8,6,4)])
barplot(Type_frq,col=bcol,las=2,xlim = c(0, 45),
        legend = TRUE, args.legend = list(bty = "n", x = "right", ncol = 1))
#legend(x=42,y=70,legend = rownames(Type_frq), fill = bcol,bty="n")
dev.off()


pdf("terminal-MMEJ.pdf",height=15,width=12,onefile = T)
for (i in c(1:35)){
  dat<-read.delim(paste0(name1[i],".info"))
  dat<-dat[dat$Type==2,]
  dat$Len<-nchar(dat$Seq)
  if (i==1){
    dat1<-data.frame(table(rep(dat$Len,dat$Frq)))
    colnames(dat1)<-c("Len",name1[i])
  } else {
    tmp<-data.frame(table(rep(dat$Len,dat$Frq)))
    colnames(tmp)<-c("Len",name1[i])
    dat1<-merge(dat1,tmp,by=1,all=T)
  }
}
dat1<-dat1[2:6,]
dat1[is.na(dat1)]<-0
write.table(dat1,"terminal_MMEJ_length.txt",quote=F,sep="\t",row.names=F)
dev.off()