











library(ggcorrplot)
library(dplyr)
library(magrittr)



library(reshape2)





setwd("~/Desktop/master_thesis/padloc/defense_count")

alldf<-read.csv("allspe_genome.csv",header = TRUE)

alldf<-alldf[alldf$X!="",]
rownames(alldf)<-alldf$X
alldf<-alldf[,colnames(alldf)!="X"]
########
df1<-alldf[,colnames(alldf)!="spename"]

df2<-data.frame(apply(df1, 2, as.numeric))
row.names(df2)<-row.names(df1)


zero=rep(0,nrow(df2))
df2['systems_count']<-zero
df2['family_count']<-zero
for (i in 1:nrow(df2)){
 df2$systems_count[i]=sum(df2[i,1:58])
 df2$family_count[i]=sum(df2[i,1:58]!=0)
  
}
spename<-unique(alldf$spename)
spename<-droplevels(spename)
options(digits=2)
sys_per<-rep(0,58)
for (i in 1:58){
  
  sys_per[i]<-100*sum(df2[,i]!=0)/nrow(df2)
  
}

df22=rbind(df2,sys_per)
row.names(df22)[nrow(df22)]<-"system_percen"


library(ggplot2)

dff<-df22[nrow(df22),1:58]
dff1 <- as.data.frame(t(dff))



ggplot(data=dff2, aes(x=reorder(rownames(dff2), -system_percen), y=system_percen)) +ylab("% of genome encoding")+xlab("system")+
  geom_bar(stat="identity",fill="steelblue")+ theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1))+geom_text(aes(label = round(system_percen,1), vjust = -0.2),size=2)




df_cor1=df2[df2$systems_count>0,1:58]


df_cor<-cor(df_cor1)


model.matrix(~0+., data=as.data.frame(df_cor)) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2)

for (i in 1:58){
  print(df_cor[df_cor[i,]^2>0.25,i])
  
  
  
}


ggplot(df2, aes(x=family_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.5) 

ggplot(df2, aes(x=systems_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.5) 

df2$spename<-alldf$spename
#############make the heatmap for all


heatdf<-data.frame(matrix(nrow = 0, ncol = length(df2)))
colnames(heatdf)<-colnames(df2)    
heatdf<-heatdf[,1:58]


plot_lst1 <- vector("list", length = 10)
plot_lst2 <- vector("list", length = 10)
plot_lst3 <- vector("list", length = 10)
nn=1
for (k in spename){
  ddf<-df2[df2$spename==k,]
  
  zero=rep(0,nrow(ddf))
  ddf['systems_count']<-zero
  ddf['family_count']<-zero
  for (i in 1:nrow(ddf)){
    ddf$systems_count[i]=sum(ddf[i,1:58])
    ddf$family_count[i]=sum(ddf[i,1:58]!=0)
    
  }
  for (j in 1:58){
    
    sys_per[j]<-100*sum(ddf[,j]!=0)/nrow(ddf)
    
  }
  
  heatdf[nrow(heatdf)+1,]=sys_per
  
  rownames(heatdf)[nrow(heatdf)]=k
  
  
  
  
  p1<-ggplot(ddf, aes(x=family_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(size=2,binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.1) +ggtitle(paste(k," mean=",round(sum(ddf$family_count)/nrow(ddf),1)) )+theme(plot.title = element_text(hjust = 0.5,size=8))
    
  
  p2<-ggplot(ddf, aes(x=systems_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(size=2,binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.2) +ggtitle(paste(k," mean=" ,round(sum(ddf$systems_count)/nrow(ddf),1)))+theme(plot.title = element_text(hjust = 0.5,size=8))
  
 
  dff1 <- as.data.frame(t(heatdf[nrow(heatdf),]))
  colnames(dff1)<-"system_percen"
  
  dff2<-as.data.frame(dff1[dff1$system_percen>1,])
  rownames(dff2)<-rownames(dff1)[dff1$system_percen>1]
  colnames(dff2)<-"system_percen"
  
  p3<-ggplot(data=dff2, aes(x=reorder(rownames(dff2), -system_percen), y=system_percen)) +ylab("% of genome encoding")+xlab("system")+
    geom_bar(stat="identity",fill="steelblue")+ theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1,size = 6))+geom_text(aes(label = round(system_percen,1), vjust = -0.2),size=2)+ggtitle(k)+theme(plot.title = element_text(hjust = 0.5,size=8))
  
  
 
  
  plot_lst1[[nn]]=p1
  plot_lst2[[nn]]=p2
  plot_lst3[[nn]]=p3
  nn=nn+1
  
}


cowplot::plot_grid(plotlist = plot_lst1[1:10], nrow = 2)


cowplot::plot_grid(plotlist = plot_lst2[1:10], nrow = 2)



cowplot::plot_grid(plotlist = plot_lst3[10], nrow = 2)


a1<-as.data.frame(summary(as.factor(merdf6$species_name)))
a1$id<-rownames(a1)
heatdf$id<-rownames(heatdf)
aa1<-merge(heatdf,a1,by="id",all.x=TRUE)


for (i in 1:nrow(heatdf)) {
  aa1$id[i]=paste(aa1$id[i],aa1$`summary(as.factor(merdf6$species_name))`[i],sep = " N=")
  
}
aa1<-aa1[colnames(aa1)!="summary(as.factor(merdf6$species_name))"]


library(reshape2)


aa1<-aa1[,-c(2,6,12,14,16,17,19,22,25,27,40,41,43,44,49,53,56,57)]
mm <- melt(aa1,id="id")
colnames(mm)<-c("species","system","value")
ggplot(mm, aes(x = system, y =species, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  theme(axis.text.x=element_text(angle=90, hjust=1))



################################################################################
#########################################################################


#################################################for plasmid

alldf<-read.csv("defense_count/allspe_plasmids.csv",header = TRUE)

alldf<-alldf[alldf$X!="",]
rownames(alldf)<-alldf$X
alldf<-alldf[,colnames(alldf)!="X"]
########
df1<-alldf[,colnames(alldf)!="spename"]

df2<-data.frame(apply(df1, 2, as.numeric))
row.names(df2)<-row.names(df1)


zero=rep(0,nrow(df2))
df2['systems_count']<-zero
df2['family_count']<-zero
for (i in 1:nrow(df2)){
  df2$systems_count[i]=sum(df2[i,1:58])
  df2$family_count[i]=sum(df2[i,1:58]!=0)
  
}
alldf$spename


df21<-df2[df2$systems_count>0,]
nrow(df21)#1650

spename<-unique(alldf$spename)
spename<-droplevels(spename)


sys_per<-rep(0,58)


for (i in 1:58){
  
  sys_per[i]<-100*sum(df21[,i]!=0)/nrow(df21)
  
}

df22=rbind(df21,sys_per)
row.names(df22)[nrow(df22)]<-"system_percen"


dff<-df22[nrow(df22),1:58]
dff2 <- as.data.frame(t(dff))

dff2$systems<-row.names(dff2)

dff2<-dff2[dff2$system_percen>1,]

selsys<-dff2$systems

ggplot(data=dff2, aes(x=reorder(rownames(dff2), -system_percen), y=system_percen)) +ylab("% of plasmid encoding")+xlab("system")+
  geom_bar(stat="identity",fill="steelblue")+ theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1))+geom_text(aes(label = round(system_percen,1), vjust = -0.2),size=2)



ggplot(df21, aes(x=family_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.5) 

ggplot(df21, aes(x=systems_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.5) 

df21$spename<-alldf$spename[df2$systems_count>0]
#############make the heatmap for all


heatdf<-data.frame(matrix(nrow = 0, ncol = length(df21)))
colnames(heatdf)<-colnames(df21)    
heatdf<-heatdf[,1:58]

sys_per<-rep(0,58)
plot_lst1 <- vector("list", length = 10)
plot_lst2 <- vector("list", length = 10)
plot_lst3 <- vector("list", length = 10)
nn=1
df21<-droplevels(df21)

spename<-levels(df21$spename)
for (k in spename){
  ddf<-df21[df21$spename==k,]
  
  zero=rep(0,nrow(ddf))
  ddf['systems_count']<-zero
  ddf['family_count']<-zero
  for (i in 1:nrow(ddf)){
    ddf$systems_count[i]=sum(ddf[i,1:58])
    ddf$family_count[i]=sum(ddf[i,1:58]!=0)
    
  }
  for (j in 1:58){
    
    sys_per[j]<-100*sum(ddf[,j]!=0)/nrow(ddf)
    
  }
  
  heatdf[nrow(heatdf)+1,]=sys_per
  k1<-strsplit(k,"plas_")[[1]][2]
  rownames(heatdf)[nrow(heatdf)]=paste(k1,"N=",nrow(ddf))
  
  
  
  
  p1<-ggplot(ddf, aes(x=family_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(size=2,binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.1) +ggtitle(paste(k," mean=",round(sum(ddf$family_count)/nrow(ddf),1),"N=",nrow(ddf)) )+theme(plot.title = element_text(hjust = 0.5,size=8))
  
  
  p2<-ggplot(ddf, aes(x=systems_count)) + theme_minimal()+  geom_histogram(binwidth=1, colour="white", fill="steelblue")+ stat_bin(size=2,binwidth= 1,aes(y=..count.., label=..count..), geom="text", vjust=-.2) +ggtitle(paste(k," mean=" ,round(sum(ddf$systems_count)/nrow(ddf),1),"N=",nrow(ddf)))+theme(plot.title = element_text(hjust = 0.5,size=8))
  
  
  dff1 <- as.data.frame(t(heatdf[nrow(heatdf),]))
  colnames(dff1)<-"system_percen"
  dff1$systems<-row.names(dff1)
  dff1<-dff1[dff1$system_percen>1,]

  
  p3<-ggplot(data=dff1, aes(x=reorder(systems, -system_percen), y=system_percen)) +ylab("% of genome encoding")+xlab("system")+
    geom_bar(stat="identity",fill="steelblue")+ theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1,size = 6))+geom_text(aes(label = round(system_percen,1), vjust = -0.2),size=2)+ggtitle(k)+theme(plot.title = element_text(hjust = 0.5,size=8))
  
  
  
  
  plot_lst1[[nn]]=p1
  plot_lst2[[nn]]=p2
  plot_lst3[[nn]]=p3
  nn=nn+1
  
}


cowplot::plot_grid(plotlist = plot_lst1[1:6], nrow = 2)


cowplot::plot_grid(plotlist = plot_lst2[1:6], nrow = 2)



cowplot::plot_grid(plotlist = plot_lst3[1:6], nrow = 2)



heatdf<-heatdf[selsys]
heatdf$id<-row.names(heatdf)


mm <- melt(heatdf,id="id")
colnames(mm)<-c("species","systems_plasmid","value")
ggplot(mm, aes(x = systems_plasmid, y =species, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(colors = hcl.colors(20, "Dark2")) +geom_text(aes(label = round(value, 1))) +
  theme(axis.text.x=element_text(angle=90, hjust=1))+ theme(axis.text.y=element_text(angle=0, hjust=1,size = 8))


ss1<-df21[df21$CRISPR.Cas>0 &df21$systems_count >3,]



ID<-rownames(ss1)

for (i in 1:length(ID)){ID[i]=strsplit(ID[i],"_genomic")[[1]]}
ss1$ID<-ID


subme<-merge(ss1,merdf6[c("n_plasmid","ID")],by="ID",all.x = TRUE)

###################################################################













#######################merge with other to test the genome size




merdf<-read.csv("~/Desktop/master_thesis/results/merdf8.csv")

allsize<-merdf[c('size','ID',"species_name","n_plasmid")]

id<-rownames(df2)
for (i in 1:length(id)){id[i]=strsplit(id[i],"_genomic")[[1]]}
  
df2$ID<-id
bdf<-merge(allsize,df2,by=c("ID"))

l1<-lm(size~,data=bdf)



library("ggpubr")
spename<-levels(merdf$species_name)
bdf<-bdf[,colnames(bdf)!="ID"]
for (spe in spename){
ndf<-bdf[bdf$species_name==spe,]
ndf<-droplevels(ndf)
p1<-ggscatter(ndf, x = "Retron", y = "DRT", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "systems_count", ylab = "size")+ggtitle(spe)
ndf<-bdf[,colnames(ndf)!="species_name"]
print(p1)
}
allsystems<-colnames(bdf)


lindf<-bdf[,c(2,5:62)]
aver<-sys_per

for (i in 2:59){
  sys_per[i]<-sum(lindf[,i]>0)/nrow(lindf)
  ddff<-lindf[lindf[,i]>0,]
  aver[i]<-sum(ddff[,i])/nrow(ddff)
  
}




fldf<-lindf[,c(TRUE,sys_per>0.01)]


lm1<-lm(size~.,data = ddf)
summary(lm1)

kk="Pseudomonas_aeruginosa"  
ddf<-fldf[bdf$species_name==kk,]
for (i in colnames(ddf)){
  
  p1<-ggscatter(ddf, x = "CRISPR.Cas", y = i, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "systems_count", ylab = "size")+ggtitle(i)
  ndf<-bdf[,colnames(ndf)!="species_name"]
  print(p1)
}

cormat<-round(cor(fldf),2)
melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()



d1=melted_cormat[melted_cormat$value>0.5,]
d2=d1[d1$Var1!=d1$Var2,]
ggscatter(fldf, x = "RM", y = "DMS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "systems_count", ylab = "size")+ggtitle(i)
ndf<-bdf[,colnames(ndf)!="species_name"]



##########test the correlation between crispr and other systems.



df1<-as.data.frame(cbind(average=aver,system=colnames(lindf)))

df1$average<-as.numeric(as.character((df1$average)))
df1$average<-round(df1$average,2)
df1<-na.exclude(df1)


ggplot(data=df1, aes(x=reorder(system, -average), y=average)) +ylab("average number of systems")+xlab("system")+
  geom_bar(stat="identity",fill="steelblue")+ theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1))+geom_text(aes(label = round(average,1), vjust = -0.2),size=2)


###################################################3
