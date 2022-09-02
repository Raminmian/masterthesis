

library("ggplot2")
library("ggpubr")

setwd("~/Desktop/master_thesis/results/spacer_detection/")

#all.txt or orphan.txt or nearcas

#alldf<-read.csv("summary/all.txt",sep = ";",header = FALSE)
#alldf<-read.csv("summary/nearcas.txt",sep = ";",header = FALSE)
alldf<-read.csv("summary/orphan.txt",sep = ";",header = FALSE)


colnames(alldf)<-c("ID","spacer_target","spacer_all","species","blast_type")
alldf$ID<-as.character(alldf$ID)
for (i in 1:nrow(alldf)){
alldf$ID[i]<-strsplit(alldf$ID[i],"_genomic")[[1]][1]
}


blast_types<-unique(alldf$blast_type)
spenames=unique(alldf$species)


asdf<-as.data.frame(matrix(ncol = 5,nrow = 0))

colnames(asdf)<-c("ids","plasmid","ices","viral","spacer_n")




##############################
plot_lst<- vector("list", length = 9)
j=1

ll<-rep(NA,9)
lll<-as.data.frame(cbind(ll,ll,ll))
colnames(lll)<-c("plasmid","ices","viral")
########################################################################

for (aa in spenames){
subdf<-alldf[alldf$species==aa,][c(1,2,3,5)]

subdf$blast_type=as.character(subdf$blast_type)
ids=unique(subdf$ID)
plasmid<-rep(0,length(ids))
ices<-rep(0,length(ids))
viral<-rep(0,length(ids))
spacer_n<-rep(0,length(ids))
for (k in 1:length(ids)){
  p1=subdf[subdf$ID==ids[k]&subdf$blast_type=="blast_plasmid",]
  if(nrow(p1)>0){ plasmid[k]=p1$spacer_target
  }
  
  i1=subdf[subdf$ID==ids[k]&subdf$blast_type=="blast_ices",]
  v1=subdf[subdf$ID==ids[k]&subdf$blast_type=="blast_viral",]
  
  if(nrow(i1)>0){
    
    ices[k]=i1$spacer_target
  }
 
  if(nrow(v1)>0){
  viral[k]=v1$spacer_target
  }
  
  spacer_n[k]=subdf[subdf$ID==ids[k],]$spacer_all[1]
}
newdf<-data.frame(cbind(ids,plasmid,ices,viral,spacer_n))



total_spacer<-sum(as.numeric(as.character(newdf$spacer_n)))
per_plas<-sum(as.numeric(as.character(newdf$plasmid)))/total_spacer


per_ice<-sum(as.numeric(as.character(newdf$ices)))/total_spacer

per_viral<-sum(as.numeric(as.character(newdf$viral)))/total_spacer

per_unknown<-1-per_plas-per_ice-per_viral

df_per<-data.frame(t(cbind(per_plas,per_ice,per_viral,per_unknown)))
colnames(df_per)<-("percentage%")

df_per$`percentage%`<-round(df_per$`percentage%`*100,2)

types<-paste(rownames(df_per),":",df_per$`percentage%`,"%",sep = "")
df_per<-cbind(df_per,types)


pieplot<-ggplot(df_per, aes(x="", y=`percentage%`, fill=types)) +
  geom_col()+
  coord_polar(theta = "y")+ggtitle(paste(aa,"n_spacer=",total_spacer))  +theme_void() +theme(plot.title = element_text(hjust = 0.5,size=10)) + theme(legend.position = "right")

plot_lst[[j]]<- pieplot



plasmid<-round(sum(newdf$plasmid!="0")/nrow(newdf),4)
ices<-round(sum(newdf$ices!="0")/nrow(newdf),4)
viral<-round(sum(newdf$viral!="0")/nrow(newdf),4)

lll[j,]<-c(plasmid,ices,viral)
row.names(lll)[j]=paste(aa," N=",nrow(newdf),collapse = "")
j=j+1
asdf<-rbind(asdf,newdf)
}
#cowplot::plot_grid(plotlist =plot_lst[c(1:6,9)], nrow = 3)

colnames(asdf)<-paste("orphan_",colnames(asdf))
colnames(asdf)[1]<-"ID"

##############################################################
###############################################################
###############################################################




alldf<-read.csv("summary/nearcas.txt",sep = ";",header = FALSE)


colnames(alldf)<-c("ID","spacer_target","spacer_all","species","blast_type")
alldf$ID<-as.character(alldf$ID)
for (i in 1:nrow(alldf)){
  alldf$ID[i]<-strsplit(alldf$ID[i],"_genomic")[[1]][1]
}


blast_types<-unique(alldf$blast_type)
spenames=unique(alldf$species)


asdf1<-as.data.frame(matrix(ncol = 5,nrow = 0))

colnames(asdf1)<-c("ids","plasmid","ices","viral","spacer_n")




##############################
plot_lst<- vector("list", length = 9)
j=1

ll<-rep(NA,9)
lll<-as.data.frame(cbind(ll,ll,ll))
colnames(lll)<-c("plasmid","ices","viral")
########################################################################

for (aa in spenames){
  subdf<-alldf[alldf$species==aa,][c(1,2,3,5)]
  
  subdf$blast_type=as.character(subdf$blast_type)
  ids=unique(subdf$ID)
  plasmid<-rep(0,length(ids))
  ices<-rep(0,length(ids))
  viral<-rep(0,length(ids))
  spacer_n<-rep(0,length(ids))
  for (k in 1:length(ids)){
    p1=subdf[subdf$ID==ids[k]&subdf$blast_type=="blast_plasmid",]
    if(nrow(p1)>0){ plasmid[k]=p1$spacer_target
    }
    
    i1=subdf[subdf$ID==ids[k]&subdf$blast_type=="blast_ices",]
    v1=subdf[subdf$ID==ids[k]&subdf$blast_type=="blast_viral",]
    
    if(nrow(i1)>0){
      
      ices[k]=i1$spacer_target
    }
    
    if(nrow(v1)>0){
      viral[k]=v1$spacer_target
    }
    
    spacer_n[k]=subdf[subdf$ID==ids[k],]$spacer_all[1]
  }
  newdf<-data.frame(cbind(ids,plasmid,ices,viral,spacer_n))
  
  
  
  total_spacer<-sum(as.numeric(as.character(newdf$spacer_n)))
  per_plas<-sum(as.numeric(as.character(newdf$plasmid)))/total_spacer
  
  
  per_ice<-sum(as.numeric(as.character(newdf$ices)))/total_spacer
  
  per_viral<-sum(as.numeric(as.character(newdf$viral)))/total_spacer
  
  per_unknown<-1-per_plas-per_ice-per_viral
  
  df_per<-data.frame(t(cbind(per_plas,per_ice,per_viral,per_unknown)))
  colnames(df_per)<-("percentage%")
  
  df_per$`percentage%`<-round(df_per$`percentage%`*100,2)
  
  types<-paste(rownames(df_per),":",df_per$`percentage%`,"%",sep = "")
  df_per<-cbind(df_per,types)
  
  
  pieplot<-ggplot(df_per, aes(x="", y=`percentage%`, fill=types)) +
    geom_col()+
    coord_polar(theta = "y")+ggtitle(paste(aa,"n_spacer=",total_spacer))  +theme_void() +theme(plot.title = element_text(hjust = 0.5,size=10)) + theme(legend.position = "right")
  
  plot_lst[[j]]<- pieplot
  
  
  
  plasmid<-round(sum(newdf$plasmid!="0")/nrow(newdf),4)
  ices<-round(sum(newdf$ices!="0")/nrow(newdf),4)
  viral<-round(sum(newdf$viral!="0")/nrow(newdf),4)
  
  lll[j,]<-c(plasmid,ices,viral)
  row.names(lll)[j]=paste(aa," N=",nrow(newdf),collapse = "")
  j=j+1
  asdf1<-rbind(asdf1,newdf)
}
#cowplot::plot_grid(plotlist =plot_lst[1:9], nrow = 3,labels = c(1:9))

colnames(asdf1)<-paste("near_cas",colnames(asdf))

colnames(asdf1)[1]<-"ID"

allspacerdf<-merge(asdf,asdf1,by="ID",all = TRUE)


fisdf<- as.data.frame(lapply(allspacerdf[,c(2,3,4,5,6,7,8,9)], function(x) as.numeric(as.character(x))))

fisdf$orphan_.spacer_n[is.na(fisdf$orphan_.spacer_n)]=0
fisdf$near_cas.orphan_.spacer_n[is.na(fisdf$near_cas.orphan_.spacer_n)]=0

fffdf<-cbind(ID=allspacerdf$ID,fisdf)

write.csv(fffdf,"spacer_detection_final.csv")
fffdf$spacer_all_true<-fffdf$orphan_.spacer_n+fffdf$near_cas.orphan_.spacer_n



merdf6<-merge(merdf5,fffdf,by="ID",all.x = TRUE)

merdf6$orphan_.spacer_n[is.na(merdf6$orphan_.spacer_n)]=0


ordf<-merdf6[merdf6$acr_inhibit=="FALSE"& merdf6$CRISPRCas=="yes"& merdf6$orphan_.spacer_n!=0,]

###################################################plot for orphan array

alldf<-ordf[,c(1,27,28,29,30)]
j=1

ll<-rep(NA,9)
lll<-as.data.frame(cbind(ll,ll,ll))
colnames(lll)<-c("plasmid","ices","viral")
for (aa in spenames){
 
  newdf<-alldf[ordf$species_name==aa,]
  
  
  total_spacer<-sum(newdf$orphan_.spacer_n)
  per_plas<-sum(newdf$orphan_.plasmid)/total_spacer
  
  
  per_ice<-sum(newdf$orphan_.ices)/total_spacer
  
  per_viral<-sum(newdf$orphan_.viral)/total_spacer
  
  per_unknown<-1-per_plas-per_ice-per_viral
  
  df_per<-data.frame(t(cbind(per_plas,per_ice,per_viral,per_unknown)))
  colnames(df_per)<-("percentage%")
  
  df_per$`percentage%`<-round(df_per$`percentage%`*100,2)
  
  types<-paste(rownames(df_per),":",df_per$`percentage%`,"%",sep = "")
  df_per<-cbind(df_per,types)
  
  
  pieplot<-ggplot(df_per, aes(x="", y=`percentage%`, fill=types)) +
    geom_col()+
    coord_polar(theta = "y")+ggtitle(paste(aa,"n_spacer=",total_spacer))  +theme_void() +theme(plot.title = element_text(hjust = 0.5,size=10)) + theme(legend.position = "right")
  
  plot_lst[[j]]<- pieplot
  
  
  
  plasmid<-round(sum(newdf$plasmid!="0")/nrow(newdf),4)
  ices<-round(sum(newdf$ices!="0")/nrow(newdf),4)
  viral<-round(sum(newdf$viral!="0")/nrow(newdf),4)
  
  #lll[j,]<-c(plasmid,ices,viral)
  #row.names(lll)[j]=paste(aa," N=",nrow(newdf),collapse = "")
  
  
  
  plasmid<-round(sum(newdf$plasmid!="0")/nrow(newdf),4)
  ices<-round(sum(newdf$ices!="0")/nrow(newdf),4)
  viral<-round(sum(newdf$viral!="0")/nrow(newdf),4)
  
  lll[j,]<-c(plasmid,ices,viral)
  row.names(lll)[j]=paste(aa," N=",nrow(newdf),collapse = "")
  j=j+1
  j=j+1
  
}
cowplot::plot_grid(plotlist =plot_lst[1:9], nrow = 3,labels = c(1:9))




###############################################################


ordf<-merdf6[merdf6$acr_inhibit=="FALSE"& merdf6$CRISPRCas=="yes"& merdf6$near_cas.orphan_.spacer_n!=0,]

alldf<-ordf[,c(1,31,32,33,34)]
j=1

ll<-rep(NA,9)
lll<-as.data.frame(cbind(ll,ll,ll))
colnames(lll)<-c("plasmid","ices","viral")

colnames(alldf)<-c("ID","plasmid","ices","viral","spacers_n")
for (aa in spenames){
  
  newdf<-alldf[ordf$species_name==aa,]
  
  
  
  plasmid<-round(sum(newdf$plasmid!="0")/nrow(newdf),4)
  ices<-round(sum(newdf$ices!="0")/nrow(newdf),4)
  viral<-round(sum(newdf$viral!="0")/nrow(newdf),4)
  
  #lll[j,]<-c(plasmid,ices,viral)
  #row.names(lll)[j]=paste(aa," N=",nrow(newdf),collapse = "")
  
  
  
  plasmid<-round(sum(newdf$plasmid!="0")/nrow(newdf),4)
  ices<-round(sum(newdf$ices!="0")/nrow(newdf),4)
  viral<-round(sum(newdf$viral!="0")/nrow(newdf),4)
  
  lll[j,]<-c(plasmid,ices,viral)
  row.names(lll)[j]=paste(aa," N=",nrow(newdf),collapse = "")
  j=j+1
 
  
}




merdf6$orphan_.spacer_n[is.na(merdf6$orphan_.spacer_n)]=0

merdf6$spacer_all_true[is.na(merdf6$spacer_all_true)]=0

merdf6$targted_ice<-(merdf6$near_cas.orphan_.ices>0)&(merdf6$acr_inhibit=="FALSE")
subdf<-merdf6[merdf6$CRISPRCas=="yes"&merdf6$acr_inhibit=="FALSE",]


i=3
subdf<-subdf[subdf$species_name==spename1[i] ,]






title=paste(spename1[i],"N_false=",sum(subdf$targted_ice==FALSE))
p2<-ggboxplot(subdf, x = "targted_ice", y = "spacer_all_true",
              color="targted_ice", palette = c("#00AFBB", "#E7B800"),title = title
)+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                              element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))
print(p2)







subdf=subdf[subdf$CRISPRCas=="yes"]

ggscatter(subdf, x = "near_cas.orphan_.ices", y = "size", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "spacers_target_ices", ylab = "genome size")+ggtitle(title)








allids <- as.data.frame(read_excel("~/Desktop/master_thesis/results/allids.xlsx")[2])
allids$ID<-as.character(allids$ID)
for (i in 1:nrow(allids)){
  allids$IDS[i]<-strsplit(allids$ID[i],"_")[[1]][2]
}




prodf<-as.data.frame(prophage[,2])

pids<-unique(prodf[,1])

ccdf<-as.data.frame(summary(prodf,maxsum = length(pids)))

ccdf$Freq<-as.character(ccdf$Freq)

for (i in 1:nrow(ccdf)){
  ccdf$ID[i]<-strsplit(ccdf$Freq[i],":")[[1]][1]
  ccdf$ID[i]<-strsplit( ccdf$ID[i],"_")[[1]][2]
  ccdf$Prophage_count[i]<-strsplit(ccdf$Freq[i],":")[[1]][2]
  
}
fpdf<-ccdf[,c(4,5)]

colnames(fpdf)[1]="IDS"
merdf6$ID<-as.character(merdf6$ID)
for (i in 1:nrow(merdf6)){
  merdf6$IDS[i]=strsplit(merdf6$ID[i],"_")[[1]][2]
  
  
  
  
}


merdf7<-merge(merdf6,allids,by="IDS",all= FALSE)
merdf8<-merge(merdf7,fpdf,all.x = TRUE)
merdf8$Prophage_count[is.na(merdf8$Prophage_count)]=0
for (i in 1:9){

subdf<-merdf8[merdf8$species_name==spename1[i],]

subdf$Prophage_count=as.numeric(as.character(subdf$Prophage_count))

title=paste(spename1[i],"N=",nrow(subdf))

subdf$acr_crispr[subdf$acr_crispr!="CRISPRCas+Acr+"]="no"

p2<-ggboxplot(subdf, x = "acr_crispr", y = "size",
              color="acr_crispr", palette = c("#00AFBB", "#E7B800"),title = title
)+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                              element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))


g3<-ggscatter(subdf, x = "Prophage_count", y = "size", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "spacers_target_ices", ylab = "genome size")+ggtitle(title)




print(p2)


}
