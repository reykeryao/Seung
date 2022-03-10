rm(list=ls())
gc()
library(tidyr)
setwd("~/Documents/NGS/Seung/JA20095/")
name<-c("GsI","GsIAI","RTcas","YADD","YIDD")
type<-paste0(rep(c("Mg","Mn"),each=2),c("DNA","RNA"))
coln<-paste0(rep(name,each=4),type)
i=2
j=1
paste0(name[i],"/",type[j])
dat<-read.table(paste0(name[i],"/",type[j],".txt"),fill = 0)

c(as.integer(row.names(dat)[dat$V1=="Snap"])+2,as.integer(row.names(dat)[dat$V1=="Snap"])-1)
c(as.integer(row.names(dat)[dat$V1=="Multiple"|dat$V1=="multiple"])+2,
  as.integer(row.names(dat)[dat$V1=="Multiple"|dat$V1=="multiple"])-1)
c(as.integer(row.names(dat)[dat$V1=="Terminal"])+1,as.integer(row.names(dat)[dat$V1=="Terminal"])-1)
dim(dat)[1]

dat_snap<-dat[9:16,]
dat_MJ<-dat[3:6,]
dat_TT<-dat[18:29,]
dat_snap<-separate(dat_snap,V1,c("ID","Counts"),sep="_|:",convert=T)
dat_MJ<-separate(dat_MJ,V1,c("ID","Counts"),sep="_|:",convert=T)
dat_TT<-separate(dat_TT,V1,c("ID","Counts"),sep=":|_",convert=T)

combined_dat<-cbind(combined_dat,data.frame(c(sum(dat_snap$Counts),sum(dat_MJ$Counts),sum(dat_TT$Counts)),
                                            row.names = c("Snap","MJ","TT")))

combined_dat<-cbind(combined_dat,data.frame(c(sum(dat_snap$Counts),0,sum(dat_TT$Counts)),
                                            row.names = c("Snap","MJ","TT")))
