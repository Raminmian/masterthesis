
library("ggplot2")
library("ggpubr")
library("plyr")
library(sjmisc)
library(stringr)
library(scales)
library("plyr")
library(gt)
library(rstatix)
library(ggprism)
setwd("~/Desktop/master_thesis/results/cctyper/final_out")




df1=read.csv("all_allspe.csv",header = TRUE)
fdf=df1[df1$ID!="ID",]
rownames(fdf)=fdf$ID
fdf$near_Spacer=as.numeric(as.character(fdf$near_Spacer))
fdf$all_Spacer=as.numeric(as.character(fdf$all_Spacer))
orphan_spacer=fdf$all_Spacer-fdf$near_Spacer
findf=cbind(fdf,orphan_spacer)

findf<-merdf6
spename=list(unique(findf$species_name))[[1]]


plot_lst1 <- vector("list", length = 10)

findf$size=as.numeric(as.character(findf$size))


colnamess=c("species","percentage %")
perc=data.frame(matrix(ncol = 2, nrow = 10))
colnames(perc)=colnamess
##########################################percentage of with crisprcas
findf
spename=list(unique(findf$species_name))[[1]]




for (i in 1:length(spename)){
  subdf=findf[findf$species_name==spename[i],]
  presence=sum(subdf$CRISPRCas=="yes")
  absence=sum(subdf$CRISPRCas=="no")
  # title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)
  title=paste(spename[i]," \n", "N=",nrow(subdf))
  cdat <- ddply(subdf, "CRISPRCas", summarise, size.mean=mean(size))
  #subdf<-merdf6[merdf6$species_name==spename[i],]
  # subdf=subdf[subdf$near_cas.orphan_.spacer_n>0,]
  
  subdf<-droplevels(subdf)
  #p1=ggplot(subdf, aes(x=near_cas.orphan_.spacer_n,fill=type))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8)) +  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=CRISPRCas),
  #                                                                                                                                         linetype="dashed", size=1)
  #p1=ggplot(subdf, aes(x=size))  + geom_density(alpha=.8, fill="#69b3a2", color="#e9ecef")  +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8)) 
  p1=ggplot(subdf, aes(x=size))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=12)) + geom_density(alpha=.8, fill="#00AFBB", color="#00AFBB")
  plot_lst1[[i]]=p1
  
}

cowplot::plot_grid(plotlist = plot_lst1, nrow = 3)






for (i in 1:length(spename)){
  subdf=findf[findf$species_name==spename[i],]
  presence=sum(subdf$CRISPRCas=="yes")
  absence=sum(subdf$CRISPRCas=="no")
 # title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)
  title=paste(spename[i]," \n", "N=",nrow(subdf))
  cdat <- ddply(subdf, "CRISPRCas", summarise, size.mean=mean(size))
  #subdf<-merdf6[merdf6$species_name==spename[i],]
 # subdf=subdf[subdf$near_cas.orphan_.spacer_n>0,]
  
  subdf<-droplevels(subdf)
  #p1=ggplot(subdf, aes(x=near_cas.orphan_.spacer_n,fill=type))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8)) +  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=CRISPRCas),
                                   #                                                                                                                                         linetype="dashed", size=1)
  #p1=ggplot(subdf, aes(x=size))  + geom_density(alpha=.8, fill="#69b3a2", color="#e9ecef")  +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8)) 
  p1=ggplot(subdf, aes(x=size,fill=CRISPRCas))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=12)) +  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=CRISPRCas),
                                                                                                                                                                                                                                                                                                                              linetype="dashed", size=1)
  plot_lst1[[i]]=p1
  print(p1)
  title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)}

cowplot::plot_grid(plotlist = plot_lst1, nrow = 3)




###########################boxplot


spename1=spename[spename!='Bordetella_pertussis'&spename!='Mycobacterium_tuberculosis']

boxplot_lst1 <- vector("list", length = 8)
findf$size=as.numeric(as.character(findf$size))
for (i in 1:length(spename1)){
  subdf=findf[findf$species_name==spename1[i],]
  presence=sum(subdf$CRISPRCas=="yes")
  absence=sum(subdf$CRISPRCas=="no")
  title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)
  
  cdat <- ddply(subdf, "CRISPRCas", summarise, size.mean=mean(size))
 
  title=paste(spename1[i])
  subdf$size[subdf$CRISPRCas=="yes"]<-subdf$size[subdf$CRISPRCas=="yes"]-0.01
  
  subdf<-merdf6[merdf6$species_name=="Klebsiella_pneumoniae",]
  subdf<-droplevels(subdf)
  subdf$n_plasmid<-as.numeric(subdf$n_plasmid)
  subdf$iva_ie<-subdf$chromo_CRISPRCas=="yes" & subdf$plas_CRISPRCas =="yes"
  
  bb1<-ggboxplot(subdf, x = "iva_ie", y = "n_plasmid",
                 color = "iva_ie", palette = "jco")+ggtitle("Klebsiella_pneumoniae")+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                           element_text(hjust = 0.5))

boxplot_lst1[[i]]=b1

  
}
cowplot::plot_grid(plotlist = boxplot_lst1[1:4], nrow = 1,labels = c(1:4))

cowplot::plot_grid(plotlist = boxplot_lst1[5:8], nrow = 1,labels = c(5:8))




#plot the distribution of spacer
for (i in 1:length(spename)){
  subdf=findf[findf$species_name==spename[i],]
  levels(droplevels(subdf$type))
  presence=sum(subdf$CRISPRCas=="yes")
  absence=sum(subdf$CRISPRCas=="no")
  title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)
  
  cdat <- ddply(subdf, "CRISPRCas", summarise, size.mean=mean(size))
  
  p1=ggplot(subdf, aes(x=size,fill=CRISPRCas))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=12)) +  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=CRISPRCas),
                                                                                                                                                                            linetype="dashed", size=1)
  plot_lst1[[i]]=p1
}



cowplot::plot_grid(plotlist = plot_lst1, nrow = 5,labels = c(1:10))



spename=spename[spename!='Bordetella_pertussis']
plot_lst2 <- vector("list", length = 9)

#plot the spacer number
for (i in 1:length(spename)){
  subdf=merdf6[merdf6$species_name==spename[i],]
  subdf<-droplevels(subdf)
  subdf=subdf[subdf$spacer_all_true>0,]
 
  absence=nrow(subdf)
  title=paste(spename[i]," \n", " N=",absence, " mean=", round(mean(subdf$spacer_all_true),0))
  print(spename[i])
  print(sum(subdf$orphan_.spacer_n>0)/nrow(subdf))
  p1=ggplot(subdf, aes(x=spacer_all_true))  + geom_density(alpha=.8, fill="#69b3a2", color="#e9ecef") +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8))+labs( x = "Spacers_number")
  plot_lst1[[i]]=p1
}
cowplot::plot_grid(plotlist = plot_lst1, nrow = 3,labels = c(1:9))




#plot the orphan number

#plot the spacer number
for (i in 1:length(spename)){
  subdf=merdf6[merdf6$species_name==spename[i],]
  subdf<-droplevels(subdf)
  subdf$near_cas.orphan_.spacer_n<-as.numeric(as.character( subdf$near_cas.orphan_.spacer_n))
  subdf$near_cas.orphan_.spacer_n[is.na(subdf$near_cas.orphan_.spacer_n)]=0
  subdf=subdf[subdf$near_cas.orphan_.spacer_n>0,]
  
  absence=nrow(subdf)
  title=paste(spename[i]," \n", " N=",absence, " mean=", round(mean(subdf$near_cas.orphan_.spacer_n),0))
  print(spename[i])

  p1=ggplot(subdf, aes(x=near_cas.orphan_.spacer_n))  + geom_density(alpha=.8, fill="#00AFBB", color="#e9ecef") +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8))+labs( x = "Spacers_number")
  plot_lst1[[i]]=p1
}
cowplot::plot_grid(plotlist = plot_lst1, nrow = 3,labels = c(1:9))


######################################plot the types of crisprcas of spacer number near cas spacers

for (i in 1:length(spename)){
 
  subdf<-merdf6[merdf6$species_name==spename[i],]
  subdf=subdf[subdf$near_cas.orphan_.spacer_n>0,]
  
  subdf<-droplevels(subdf)
  p1=ggplot(subdf, aes(x=near_cas.orphan_.spacer_n,fill=type))  + geom_density(alpha=.5) +theme_classic()+ggtitle(spename[i])+theme(plot.title = element_text(hjust = 0.5,size=8))
  plot_lst1[[i]]=p1
  
  }

cowplot::plot_grid(plotlist = plot_lst1[6:9], nrow = 1,labels = c(1:5))

















##########################################################


spename1=spename[spename!='Bordetella_pertussis']
bbplot_lst<- vector("list", length = 9)
plot_lst2<- vector("list", length = 9)
pieplot_lst <- vector("list", length = 9)
findf$size=as.numeric(as.character(findf$size))
findf$ID=as.character(findf$ID)
for (i in 1:length(spename1)){
  subdf=findf[findf$species_name==spename1[i],]
  presence=sum(subdf$CRISPRCas=="yes")
  absence=sum(subdf$CRISPRCas=="no")
  k=i
  title1<-paste(spename1[i], " N=", nrow(subdf))
  tt<-levels(droplevels(subdf$type))

  tt=tt[tt!=" NA "]
  newt=paste(tt,collapse = '')
  newfactor=strsplit(newt,"\\|")
  
  alltypes=unique(newfactor[[1]])
  print(alltypes)
  
  
  
  subdf$crisprcas_type_n<-rep(0,nrow(subdf))
  
  nn=lapply(subdf$type, f<-function(i)str_count(i,"\\|"))
  for (i in 1:nrow(subdf)){
    subdf$crisprcas_type_n[i]=nn[[i]][1]
  }
  
  for (i in alltypes){
    tyy=rep("no",nrow(subdf))
    for (j in 1:nrow(subdf)){
      
      if(str_contains(subdf$type[j],i)){
        tyy[j]="yes"
      }
    }
    subdf=cbind(subdf,tyy)
    colnames(subdf)[length(colnames(subdf))]=i
    
  }
  
  type_count<-data.frame(matrix(ncol = 2, nrow = 2))
  colnames(type_count)<-c("types","count")
  
  end=dim(subdf)[2]
  allty=rep(NA,nrow(subdf))
  for (i in 1:nrow(subdf)){
  types<-colnames(subdf)[13:end][subdf[i,13:end]=="yes"]
  ttype<-paste(types,collapse = "&")
  allty[i]=ttype
  }
  
  ##count the pieplot
  allty=sub("^$", "no", allty)
  subdf=cbind(subdf,allty)
  CRISPRCas_type=list(unique(sub("^$", "no", allty)))[[1]]

  counts<-rep(NA,length(CRISPRCas_type))
  for (i in 1:length(CRISPRCas_type)){
    counts[i]=sum(allty==CRISPRCas_type[i])
  }
  
  
  counts<-round(counts/sum(counts)*100,2)
  
  
  ddf<-data.frame(cbind(CRISPRCas_type,counts))
  ddf$counts=as.numeric(as.character(ddf$counts))
  ddf$CRISPRCas_types<-paste(ddf$CRISPRCas_type,":",counts,"%",sep = "")
  
  pieplot<-ggplot(ddf, aes(x="", y=counts, fill=CRISPRCas_types)) +
    geom_col() +
    coord_polar(theta = "y")+ggtitle(title1)  +theme_void() +theme(plot.title = element_text(hjust = 0.5,size=12)) + theme(legend.position = "right")
  
  print(k)
  pieplot_lst[[k]]=pieplot
  ttt<-ddf$CRISPRCas_type[ddf$counts>3]
  
  
  ##############################type and genome size, smaller than 10% will not be considered
  
  subdf=subdf[subdf$allty %in% ttt,]
  subdf=droplevels(subdf)
  
  cdat <- ddply(subdf, "allty", summarise, size.mean=mean(size))
  
  p1=ggplot(subdf, aes(x=size,fill=allty,color=allty))  + geom_density(alpha=.4) +theme_classic()+ggtitle(spename1[k])+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12)) +  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=allty),
                                                                                                                                                                             linetype="dashed", size=1)
  
  print(p1)
  plot_lst2[[k]]=p1
  
  
  
  #######################boxplot for it
  bxp <- ggboxplot(
    subdf, x = "allty", y = "size", 
    color = "allty", palette = c("#00AFBB", "#E7B800")
  )
  colnames(subdf)[colnames(subdf)=="allty"]<-"CRISPRCas_type"
  
  
if (k %in% c(2,3,5,6,7)){
  df_p_val <- subdf %>% 
    rstatix::wilcox_test(size ~ CRISPRCas_type) %>% 
    rstatix::add_xy_position()
  
  
  subdf<-subdf[subdf$acr_crispr!="CRISPRCas+Acr+",]
  bb1<-ggboxplot(subdf, x = "CRISPRCas_type", y = "size",
            color = "CRISPRCas_type", palette = "jco")+ggtitle(spename1[k])+             
    add_pvalue(df_p_val,label.size = 3,bracket.size = 0.3)+theme(plot.title = element_text(hjust = 0.5,size=12))
                                                                                                        
    
  bbplot_lst[[k]]=bb1
  
    
  
}}
cowplot::plot_grid(plotlist = plot_lst2[c(2,3,5,6,7)], nrow = 3,labels = c(1:5))
cowplot::plot_grid(plotlist = bbplot_lst[c(2,3,5,6,7)], nrow = 2,labels = c(1:5))
cowplot::plot_grid(plotlist = pieplot_lst[c(1,2,3,4,5,6,7,8,9)], nrow = 3,labels = c(1:9))




####################################################plasmid information





pldf=read.csv("../../plasmid_extract/allcount.txt",header = FALSE,sep = ";")
colnames(pldf)<-c("ID","n_plasmid","len_chro","len_all","species_name")
spename<-list(unique(pldf$species_name))[[1]]




percentage<-c()
for (i in spename){
  per<-pldf[pldf$species==i,]
  perpl<-round(sum(per$n_plasmid>0)/nrow(per)*100,2)
  percentage<-c(percentage,perpl)
}

species<-c("A.baumannii","E.coli","K.pneumoniae","L.monocytogenes"," M.tuberculosis","P.aeruginosa","S.enterica","S.aureus","S.pyogenes")
perc<-data.frame(cbind(species,percentage))


perc$percentage=as.numeric(as.character(perc$percentage))

ggbarplot(perc, "species", "percentage",
          fill = "steelblue", color = "steelblue",
          label = TRUE, lab.pos = "out", lab.col = "black")+theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=0.5,size=10))+scale_x_discrete("species")
################################
ssname<-c("Acinetobacter_baumannii","Escherichia_coli","Klebsiella_pneumoniae","Staphylococcus_aureus","Salmonella_enterica")

merdf <- join_all(list(findf,pldf[,1:(ncol(pldf)-1)]), by = 'ID')

plot_lst2<- vector("list", length = 5)
k=0
for (i in ssname){
  df_sub<-merdf[merdf$species_name==i,]
  
  print(nrow(df_sub))
 
  df_sub$n_plasmid[df_sub$n_plasmid!="0"]=">0"
  
  p1<-ggboxplot(df_sub, x = "n_plasmid", y = "size",
            color = "n_plasmid", palette = "jco")+ggtitle(i)+ stat_compare_means()+theme(plot.title = element_text(hjust = 0.5,size=12))
  k=k+1
  plot_lst2[[k]]=p1
  
}
cowplot::plot_grid(plotlist = plot_lst2[c(1:5)], nrow = 1,labels = c(1:5))


#########################################only plasmid

ppdf<-read.csv("onlyplasmid/plas_allspe.csv",header = TRUE)
ppdf<-ppdf[ppdf$ID!="ID",]
spename<-list(unique(ppdf$species_name))[[1]]

merdf$len_plas=(merdf$len_all-merdf$len_chro)/1000
merdf1 <- join_all(list(ppdf,merdf[,c("ID",'n_plasmid','len_plas','len_chro')]), by = 'ID')



parpd<-ppdf[,c(2,3)]
colnames(parpd)[2]<-"plas_CRISPRCas"
merdf2 <- join_all(list(merdf,parpd), by = 'ID')

########################chromo
chdf<-read.csv("chromo/chromo_allspe.csv",header = TRUE)
chdf<-chdf[chdf$ID!="ID",]


ccpd<-chdf[,c(2,3)]
colnames(ccpd)[2]<-"chromo_CRISPRCas"

merdf3 <- join_all(list(merdf2,ccpd), by = 'ID')

#####check chromosone and plasmid crisprcas on kle and e.coli



kk=spename[3]
qqdf<-merdf3[merdf3$species_name==kk,]


kk=spename[4]
qqdf<-merdf3[merdf3$species_name==kk,]


a1=qqdf[(qqdf$chromo_CRISPRCas=="no"&qqdf$plas_CRISPRCas=="yes"),]
a2<-qqdf[(qqdf$chromo_CRISPRCas=="yes"&qqdf$plas_CRISPRCas=="yes"),]



c1=rep("t1",nrow(a1))
c2=rep("t2",nrow(a2))
combb<-c(c1,c2)
yesno<-cbind(rbind(a1,a2),combb)



p1=ggplot(yesno, aes(x=len_chro,fill=combb))  + geom_density(alpha=.5) +theme_classic()+ggtitle(kk)+theme(plot.title = element_text(hjust = 0.5,size=12)) 
print(p1)

ggboxplot(yesno, x = "combb", y = "len_chro",
          color = "combb", paette = "jco")+ggtitle(i)+ stat_compare_means()+theme(plot.title = element_text(hjust = 0.5,size=12))


p1=ggplot(yesno, aes(x=size,fill=combb))  + geom_density(alpha=.5) +theme_classic()+ggtitle(kk)+theme(plot.title = element_text(hjust = 0.5,size=12)) 
print(p1)

p1=ggplot(yesno, aes(x=len_plas,fill=combb))  + geom_density(alpha=.5) +theme_classic()+ggtitle(kk)+theme(plot.title = element_text(hjust = 0.5,size=12)) 
print(p1)



withcris<-c()

pieplot_lst<-vector("list", length = 4)
kk=0
for (i in spename){
  subdf=merdf1[as.character(merdf1$species_name)==i,]
  
  subdf=subdf[!is.na(subdf$CRISPRCas),]
  subdf$len_chro=as.numeric(as.character(subdf$len_chro))
  #######################
  p1=ggplot(subdf, aes(x=len_chro,fill=CRISPRCas))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=12)) 
  print(p1)
  

  ##############################################################################################################3
    subdf=subdf[subdf$n_plasmid!="0",]
    withcris<-c(withcris,round(100*sum(subdf$CRISPRCas=="yes")/nrow(subdf),2))
  
  kk=kk+1
  title11=paste(i,"\n"," N=",nrow(subdf))
  
  
  tt<-levels(droplevels(subdf$type))
  tt=tt[tt!=" NA "]
  newt=paste(tt,collapse = '')
  newfactor=strsplit(newt,"\\|")
  
  alltypes=unique(newfactor[[1]])
  print(alltypes)
  
  
  
  subdf$crisprcas_type_n<-rep(0,nrow(subdf))
  
  nn=lapply(subdf$type, f<-function(i)str_count(i,"\\|"))
  for (i in 1:nrow(subdf)){
    subdf$crisprcas_type_n[i]=nn[[i]][1]
  }
  
  for (i in alltypes){
    tyy=rep("no",nrow(subdf))
    for (j in 1:nrow(subdf)){
      
      if(str_contains(subdf$type[j],i)){
        tyy[j]="yes"
      }
    }
    subdf=cbind(subdf,tyy)
    colnames(subdf)[length(colnames(subdf))]=i
    
  }
  
  type_count<-data.frame(matrix(ncol = 2, nrow = 2))
  colnames(type_count)<-c("types","count")
  
  end=dim(subdf)[2]
  allty=rep(NA,nrow(subdf))
  for (i in 1:nrow(subdf)){
    types<-colnames(subdf)[13:end][subdf[i,13:end]=="yes"]
    ttype<-paste(types,collapse = "&")
    allty[i]=ttype
  }
  
  ##count the pieplot
  allty=sub("^$", "no", allty)
  subdf=cbind(subdf,allty)
  CRISPRCas_type=list(unique(sub("^$", "no", allty)))[[1]]
  
  counts<-rep(NA,length(CRISPRCas_type))
  for (i in 1:length(CRISPRCas_type)){
    counts[i]=sum(allty==CRISPRCas_type[i])
  }
  
  

  
  
  ddf<-data.frame(cbind(CRISPRCas_type,counts))
  ddf$counts=as.numeric(as.character(ddf$counts))
  
##########################################crisprcas typer on plasmid
  pieplot<-ggplot(ddf, aes(x="", y=counts, fill=CRISPRCas_type)) +
    geom_col() +
    geom_text(aes(label = paste(counts,"%")),
              position = position_stack(vjust = 0.5),size=3) +
    coord_polar(theta = "y")+ggtitle(title11)  +theme_void() +theme(plot.title = element_text(hjust = 0.5,size=12)) + theme(legend.position = "bottom")
  
  print(pieplot)
  pieplot_lst[[kk]]=pieplot

  subdf$size=as.numeric(as.character(subdf$size))
 
 # ggb<-ggboxplot(subdf, x ="CRISPRCas", y = "n_plasmid",
           # color = "CRISPRCas", palette = "jco")+ggtitle(title11)+ stat_compare_means()+theme(plot.title = element_text(hjust = 0.5,size=12))
  
# print(ggb)

}



ddf<-data.frame(cbind(as.character(spename),withcris))


############################################chromosal
spename<-list(unique(merdf3$species_name))[[1]]
spename<-spename[spename!="Bordetella_pertussis"&spename!="Mycobacterium_tuberculosis"]


plot_lst1 <- vector("list", length = 8)

boxplot_lst1 <- vector("list", length = 8)

for (i in 1:length(spename)){
  subdf=merdf3[merdf3$species_name==spename[i],]
  presence=sum(subdf$CRISPRCas=="yes")
  absence=sum(subdf$CRISPRCas=="no")
  title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)
  
  cdat <- ddply(subdf, "CRISPRCas", summarise, size.mean=mean(len_chro))
  subdf$len_chro<-as.numeric(subdf$len_chro)
  subdf=subdf[subdf$len_chro>0,]
  
  p1=ggplot(subdf, aes(x=len_chro,fill=chromo_CRISPRCas))  + geom_density(alpha=.5) +theme_classic()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=8)) +  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=CRISPRCas),
                                                                                                                                                                            linetype="dashed", size=1)
  plot_lst1[[i]]=p1
  
  title=paste(spename[i]," \n", "Yes=", presence, " No=",absence)
  
  
  b1=ggboxplot(subdf, x = "CRISPRCas", y = "len_chro",
               color="CRISPRCas", palette = c("#00AFBB", "#E7B800"),title = title
  )+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                element_text(hjust = 0.5))
  
  boxplot_lst1[[i]]=b1
  
  }

cowplot::plot_grid(plotlist = plot_lst1, nrow = 4,labels = c(1:8))


cowplot::plot_grid(plotlist = boxplot_lst1[1:8], nrow = 2,labels = c(1:8))



######################################################chromosal size







#################################################################acr
acrdf<-read.csv("../../acrfinder/final_out/allacr.csv")
acrdf<-acrdf[,colnames(acrdf)!="X"]

mergdf4<-merge(merdf3,acrdf,by="ID",all.x = TRUE)
###################



mergdf4[is.na(mergdf4)]="no"

homo_count<-rep(NA,9)
high<-rep(NA,9)
medium<-rep(NA,9)
low<-rep(NA,9)
only_gba<-rep(NA,9)
gba_homology<-rep(NA,9)

i=1
for (k in spename1){
  ndf<-mergdf4[mergdf4$species_name==k,]
  ndf<-droplevels(ndf)
  high[i]=sum(ndf$high_condence_acr_GBA!="no")/nrow(ndf)*100
  medium[i]=sum(ndf$medium_condence_acr_GBA!="no")/nrow(ndf)*100
  low[i]=sum(ndf$low_condence_acr_GBA!="no")/nrow(ndf)*100
  gba<-(ndf$medium_condence_acr_GBA=="no")&(ndf$high_condence_acr_GBA=="no")&(ndf$low_condence_acr_GBA=="no")
  homolgy_count[i]=sum((ndf$acr_homology!="no"))/nrow(ndf)*100
  gba_homology[i]=sum((ndf$acr_homology!="no")&gba==FALSE)/nrow(ndf)*100
  
  i=i+1
}
d2f<-as.data.frame(cbind(spename=as.character(spename1),homo_count,gba_low_confidence_acr=low,gba_medium_confidence_acr=medium,gba_high_confidence_acr=high,gba_homology))

df3<-melt(d2f,id="spename")
colnames(df3)<-c("species","acr_count","percentage")
df3_so<-arrange(df3, species,desc(acr_count) )

df3$percentage<-round(as.numeric(df3$percentage),2)


ggplot(data=df3, aes(x=species, y=percentage,fill=acr_count)) +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ coord_flip()



spename1
stadf<-mergdf4[mergdf4$species_name=="Staphylococcus_aureus",]

sum(stadf$CRISPR !="no")/nrow(stadf)

k="Staphylococcus_aureus"


###########################################################show if it target on the homolgy
hhdf<-mergdf4[(mergdf4$acr_homology!="no")&(mergdf4$CRISPRCas!="no"),]
qdf<-hhdf[c("type","acr_homology")]
qdf$type<-as.character(qdf$type)
qdf$acr_homology<-as.character(qdf$acr_homology)
inhibit<-rep(NA,nrow(hhdf))

for (i in 1:nrow(qdf)){
  ssp<-strsplit(qdf[i,1],"\\|")[[1]]
  inh<-TRUE
  for (k in ssp){
    ie<-strsplit(k,"-")[[1]]
    tar<-paste(c("Acr",ie),collapse = "")
    inh<-inh& !str_contains(qdf[i,2],tar)
  }
  inhibit[i]=!inh
  
}
homo<-as.data.frame(cbind(ID=as.character(hhdf$ID),homo_acr_inhibit=as.character(inhibit)))


###################
hhdf<-mergdf4[(mergdf4$high_condence_acr_GBA!="no")|(mergdf4$medium_condence_acr_GBA!="no"),]

inhibit_GBA<-rep(TRUE,nrow(hhdf))

gba_inh<-as.data.frame(cbind(ID=as.character(hhdf$ID),gba_acr_inhibit=as.character(inhibit_GBA)))

both<-merge(gba_inh,homo,by="ID",all = TRUE)
merdf5<-merge(mergdf4,both,by="ID",all.x = TRUE)
merdf5$homo_acr_inhibit=as.character(merdf5$homo_acr_inhibit)
merdf5$gba_acr_inhibit=as.character(merdf5$gba_acr_inhibit)
merdf5[is.na(merdf5)]<-"FALSE"
merdf5$acr_inhibit<-(merdf5$gba_acr_inhibit=="TRUE")|(merdf5$homo_acr_inhibit=="TRUE")





##############make the plot

k=1
merdf5$species_name<-as.character(merdf5$species_name)
merdf5$len_chro<-as.numeric(as.character(merdf5$len_chro))/1000000
for (i in 1:length(spename1)){
  subdf=merdf5[merdf5$species_name==spename1[i],]
  subdf=subdf[subdf$CRISPRCas=="yes",]
  yess=paste("TRUE: N=",sum(subdf$acr_inhibit=="TRUE"),collapse = "")
  noo=paste("FALSE: N=",sum(subdf$acr_inhibit=="FALSE"),collapse = "")
  bb<-paste(yess,noo,collapse = "")
  if (sum(subdf$acr_inhibit==TRUE)>0){
    p1=ggplot(subdf, aes(x=size,fill=acr_inhibit))  + geom_density(alpha=.5) +theme_classic()+ggtitle(spename1[i])+theme(plot.title = element_text(hjust = 0.5,size=8))+scale_fill_discrete(labels=c(noo,yess))
    p2<-ggboxplot(subdf, x = "acr_inhibit", y = "size",
                  color="acr_inhibit", palette = c("#00AFBB", "#E7B800"),title = paste(spename1[i],"\n",bb,collapse = "")
    )+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                        element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))
      
    
    
    plot_lst1[[k]]=p1
    plot_lst2[[k]]=p2
    k=k+1
  }
}

cowplot::plot_grid(plotlist = plot_lst1[1:7], nrow = 2,labels = c(1:7))
cowplot::plot_grid(plotlist = plot_lst2[1:4], nrow = 1)




############################################################crispr + without acr and 

##############make the plot

k=1
merdf5$species_name<-as.character(merdf5$species_name)
merdf5$acr_crispr=rep("no",nrow(merdf5))
merdf5$acr_crispr[(merdf5$CRISPRCas=="yes")&(merdf5$acr_inhibit=="TRUE")]="CRISPRCas+Acr+"
merdf5$acr_crispr[(merdf5$CRISPRCas=="yes")&(merdf5$acr_inhibit=="FALSE")]="CRISPRCas+Acr-"
merdf51<-merdf5[merdf5$acr_crispr!="no",]
merdf51$GC<-as.numeric(as.character(merdf51$GC))

for (i in 1:length(spename1)){
 
  subdf=merdf51[merdf51$species_name==spename1[i],]
  subdf$acr_crispr <- factor(subdf$acr_crispr , levels=c("CRISPRCas+Acr+", "CRISPRCas+Acr-"))
  subdf<-subdf[subdf$GC>10,]
  subdf=subdf[subdf$type!="IV-A3|",]
  yess=paste("CRISPRCas+ Acr+: N=",sum(subdf$acr_crispr=="CRISPRCas+Acr+"),collapse = "")
  noo=paste("CRISPRCas+ Acr-: N=",sum(subdf$acr_crispr=="CRISPRCas+Acr-"),collapse = "")
  
  bb<-paste(yess,noo,collapse = "")
    p2<-ggboxplot(subdf, x = "acr_crispr", y = "GC",
                  color="acr_crispr", palette = c("#00AFBB", "#E7B800"),title = paste(spename1[i],"\n",bb,collapse = "")
    )+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                  element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))
    
    
    print(p2)
    plot_lst1[[k]]=p1
    plot_lst2[[k]]=p2
    k=k+1
  
}

cowplot::plot_grid(plotlist = plot_lst1[1:7], nrow = 2,labels = c(1:6))
cowplot::plot_grid(plotlist = plot_lst2[c(2,3,5,6)], nrow = 1)
cowplot::plot_grid(plotlist = plot_lst2[c(1,7,9)], nrow = 1)

cowplot::plot_grid(plotlist = plot_lst2[c(5:7,9)], nrow = 1)


#################################### with Acr tends to high larger size?


merdf5$species_name<-as.character(merdf5$species_name)
merdf5$acr_homology<-as.character(merdf5$acr_homology)
k=1
for (i in 1:length(spename1)){
  subdf=merdf5[merdf5$species_name==spename1[i],]

  subdf=subdf[subdf$CRISPRCas=="no",]
  yess=paste("TRUE: N=",sum(subdf$acr_homology!="no"),collapse = "")
  
  noo=paste("FALSE: N=",sum(subdf$acr_homology=="no"),collapse = "")
  bb<-paste(yess,noo,collapse = "")
  subdf$acr_homology[subdf$acr_homology!="no"]="yes"
  
  p1=ggplot(subdf, aes(x=size,fill=acr_homology))  + geom_density(alpha=.5) +theme_classic()+ggtitle(spename1[i])+theme(plot.title = element_text(hjust = 0.5,size=8))+scale_fill_discrete(labels=c(noo,yess))
  p2<-ggboxplot(subdf, x = "acr_homology", y = "size",
                color="acr_homology", palette = c("#00AFBB", "#E7B800"),title = paste(spename1[i],"\n",bb,collapse = "")
  )+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))
  
  
  
  plot_lst1[[k]]=p1
  plot_lst2[[k]]=p2
  k=k+1
  
}

cowplot::plot_grid(plotlist = plot_lst1[1:7], nrow = 2,labels = c(1:7))
cowplot::plot_grid(plotlist = plot_lst2[5:9], nrow = 1)


cons<-rep(FALSE,nrow(merdf5))
for (i in 1:nrow(merdf5)){
  
  
  con<-str_contains(merdf5$acr_homology[i],"AcrIIA2")
  cons[i]<-con
  
  
}
merdf5[cons,c(6,19)]

##############################################################################################genome size and 


k=1
merdf5$species_name<-as.character(merdf5$species_name)
merdf5$acr_crispr=rep("no",nrow(merdf5))
merdf5$acr_crispr[(merdf5$CRISPRCas=="yes")&(merdf5$acr_inhibit=="FALSE")&(merdf5$acr_homology!="ni")]="CRISPRCas+Acr-"
merdf5$acr_crispr[(merdf5$CRISPRCas=="no")&(merdf5$acr_homology=="no")]= "CRISPRCas-Acr-"
merdf51<-merdf5[merdf5$acr_crispr!="no",]
for (i in 1:length(spename1)){
  
  subdf=merdf51[merdf51$species_name==spename1[i],]
  subdf$acr_crispr <- factor(subdf$acr_crispr , levels=c("CRISPRCas-Acr-", "CRISPRCas+Acr-"))
  
  subdf=subdf[subdf$type!="IV=A3|",]
  yess=paste("CRISPRCas- Acr-: N=",sum(subdf$acr_crispr=="CRISPRCas-Acr-"),collapse = "")
  noo=paste("CRISPRCas+ Acr-: N=",sum(subdf$acr_crispr=="CRISPRCas+Acr-"),collapse = "")
  
  bb<-paste(yess,noo,collapse = "")
  p2<-ggboxplot(subdf, x = "acr_crispr", y = "size",
                color="acr_crispr", palette = c("#00AFBB", "#E7B800"),title = paste(spename1[i],"\n",bb,collapse = "")
  )+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))
  
  
  print(p2)
  plot_lst1[[k]]=p1
  plot_lst2[[k]]=p2
  k=k+1
  
}

cowplot::plot_grid(plotlist = plot_lst1[1:7], nrow = 2,labels = c(1:6))
cowplot::plot_grid(plotlist = plot_lst2[1:4], nrow = 1)

cowplot::plot_grid(plotlist = plot_lst2[c(5:7,9)], nrow = 1)



a1<-merdf5[merdf5$species_name!="Staphylococcus_aureus",]

sum(a1$CRISPRCas=="no"&a1$CRISPR=="yes")/nrow(a1)
a1<-a1[a1$acr_homology!="no",]
a2<-a1[a1$CRISPRCas=="no",]
a3<-a2[a2$CRISPR=="yes",]

nrow(a3)/nrow(a1)
type_count<-rep(NA,nrow(merdf5))
for (i in 1:nrow(merdf5)){
bb=merdf5$type[i]
bb<-as.character(bb)
type_count[i]<-length(grep("-" , strsplit(bb,"|")[[1]]))
}

ddf<-merdf5[type_count==0,]

ddf<-ddf[ddf$species_name!="Staphylococcus_aureus",]

ddf<-ddf[ddf$plas_CRISPRCas=="no",]
dim(ddf)
sum(ddf$acr_inhibit=='TRUE')/nrow(ddf)
ddf<-ddf[ddf$acr_inhibit=='TRUE',]
ddf[,c(6,19,20,22)]


kll<-merdf5[merdf5$species_name=="Klebsiella_pneumoniae",]
kll<-kll[kll$CRISPR=="yes",]
sum(kll$acr_inhibit=="TRUE")/nrow(kll)



ddf[(ddf$high_condence_acr_GBA!="no")|(ddf$medium_condence_acr_GBA!="no")|ddf$acr_homology!="no",c(6,19,20,22)]

sum(merdf5$acr_inhibit=='TRUE')/sum(merdf5$CRISPRCas=="yes")



########################################################GC contente



k=1
merdf5$species_name<-as.character(merdf5$species_name)
merdf5$GC<-as.numeric(as.character(merdf5$GC))

merdf5$acr_crispr=rep("no",nrow(merdf5))
merdf5$acr_crispr[(merdf5$CRISPRCas=="yes")&(merdf5$acr_inhibit=="TRUE")]="CRISPRCas+Acr+"
merdf5$acr_crispr[(merdf5$CRISPRCas=="yes")&(merdf5$acr_inhibit=="FALSE")]="CRISPRCas+Acr-"
merdf51<-merdf5[merdf5$acr_crispr!="no",]
for (i in 1:length(spename1)){
  
  subdf=merdf51[merdf51$species_name==spename1[i],]
  subdf$acr_crispr <- factor(subdf$acr_crispr , levels=c("CRISPRCas+Acr+", "CRISPRCas+Acr-"))
  
  subdf=subdf[subdf$type!="IV-A3|",]
  yess=paste("CRISPRCas+ Acr+: N=",sum(subdf$acr_crispr=="CRISPRCas+Acr+"),collapse = "")
  noo=paste("CRISPRCas+ Acr-: N=",sum(subdf$acr_crispr=="CRISPRCas+Acr-"),collapse = "")
  
  bb<-paste(yess,noo,collapse = "")
  p2<-ggboxplot(subdf, x = "acr_crispr", y = "GC",
                color="acr_crispr", palette = c("#00AFBB", "#E7B800"),title = paste(spename1[i],"\n",bb,collapse = "")
  )+ theme(plot.title = element_text(size=12))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                element_text(hjust = 0.5,size = 10))+theme(legend.position="bottom",text = element_text(size=10))
  
  
  print(p2)
  plot_lst1[[k]]=p1
  plot_lst2[[k]]=p2
  k=k+1
  
}

cowplot::plot_grid(plotlist = plot_lst1[1:7], nrow = 2,labels = c(1:6))
cowplot::plot_grid(plotlist = plot_lst2[1:4], nrow = 1)

cowplot::plot_grid(plotlist = plot_lst2[c(5:7,9)], nrow = 1)



############################################################################################################################
##########################################################################
############################################################################
###############################################################################
#for spacers


merdf6$type<-as.character(merdf6$type)
sum(strsplit(merdf6$type[1],"|")[[1]]=="-")
ccount=rep(0,nrow(merdf6))
for (i in 1:nrow(merdf6)){
  ccount[i]<-sum(strsplit(merdf6$type[i],"|")[[1]]=="-")
}
ccount[ccount>2]=2
merdf6$CRISPRCas_count<-ccount
merdf6$CRISPRCas_count<-as.factor(merdf6$CRISPRCas_count)
merdf6$near_cas.orphan_.spacer_n[is.na(merdf6$near_cas.orphan_.spacer_n)]=0

for (i in 1:length(spename)){
  
  subdf<-merdf6[merdf6$species_name==spename[i],]
  subdf<-subdf[subdf$near_cas.orphan_.spacer_n>0,]
  subdf<-subdf[subdf$CRISPRCas_count==1,]
  subdf<-subdf[subdf$type!="IV-A1|" & subdf$type!="IV-A2|" & subdf$type!="I-C|" ,]
  subdf<-droplevels(subdf)
  b1<-ggboxplot(subdf, x = "type", y = "near_cas.orphan_.spacer_n",
                 color = "type", palette = "jco")+ggtitle(spename[i])+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                                                                                                 element_text(hjust = 0.5))
  
  boxplot_lst1[[i]]=b1
  
}
cowplot::plot_grid(plotlist = boxplot_lst1[c(2,3,5,6,7)], nrow = 1,labels = c(1:5))







########################plasmid number and chromosal CRISPR-cas

merdf6$n_plasmid<-as.numeric(merdf6$n_plasmid)
merdf6$n_plasmid[is.na(merdf6$n_plasmid)]=0
for (i in 1:length(spename)){
 
  subdf<-merdf6[merdf6$species_name==spename[i],]
  if(i!=7){
    subdf$chromo_CRISPRCas<-subdf$CRISPRCas
  }
  subdf<-droplevels(subdf)
  b1<-ggboxplot(subdf, x = "chromo_CRISPRCas", y = "n_plasmid",
                color = "chromo_CRISPRCas", palette = "jco")+ggtitle(spename[i])+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = 
                                                                                                                                                                                                                 element_text(hjust = 0.5))+ theme(legend.position="right")
    scale_fill_discrete(labels=c("no","CRISPRCas+Acr-","CRISPRCas+Acr+"))
  
  boxplot_lst1[[i]]=b1
  
}
cowplot::plot_grid(plotlist = boxplot_lst1[c(1,2,3,4,5,6,7,9)], nrow = 2,labels = c(1:5))



for (i in 1:length(spename)){
  
  subdf<-merdf6[merdf6$species_name==spename[i],]
 
  subdf<-droplevels(subdf)
  
  df_p_val <- subdf %>% 
    rstatix::wilcox_test(size ~ acr_crispr) %>% 
    rstatix::add_xy_position()
  
  
  
  if(i %in% c(1,2,3,4,5,6,7,9)){
  level_order <-c("CRISPRCas+Acr-","no","CRISPRCas+Acr+")
  subdf$acr_crispr <- factor(subdf$acr_crispr,levels = level_order)
  b1<-ggboxplot(subdf, x = "acr_crispr", y = "len_chro",
                color = "acr_crispr", palette = "jco")+ggtitle(spename[i])+ add_pvalue(df_p_val,label.size = 3,bracket.size = 0.3)+ theme(plot.title = element_text(size=14)) +theme(legend.text = element_text(size=10),axis.text.x = element_text(size = 8))+ theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position="right")
  
  boxplot_lst1[[i]]=b1
  }
  
}
cowplot::plot_grid(plotlist = boxplot_lst1[c(1,2,3,4,5,6,7,9)], nrow = 2,labels = c(1:8))




merdf6$types<-merdf6$type
merdf6$types<-as.character(merdf6$types)
merdf6$types[merdf6$types=="I-B|"]="I-B"
merdf6$types[merdf6$types=="I-B|II-A|"]="I-B&II-A"
merdf6$types[merdf6$types=="II-A|I-B|"]="I-B&II-A"
merdf6$types[merdf6$types=="I-C|"]="I-C"
merdf6$types[merdf6$types=="I-C|II-A|"]="I-C&II-A"
merdf6$types[merdf6$types=="II-A|I-C|"]="I-C&II-A"
merdf6$types[merdf6$types=="I-C|I-F|"]="I-C&I-F"
merdf6$types[merdf6$types=="I-F|I-C|"]="I-C&I-F"
merdf6$types[merdf6$types=="I-E|I-E|"]="I-E"
merdf6$types[merdf6$types=="I-E|"]="I-E"
merdf6$types[merdf6$types=="I-F|"]="I-F"
merdf6$types[merdf6$types=="I-F|I-E|"]="I-F&I-E"
merdf6$types[merdf6$types=="I-E|I-F|"]="I-F&I-E"
merdf6$types[merdf6$types=="I-F|I-F|"]="I-F"
merdf6$types[merdf6$types=="I-E|IV-A3|"]="I-E&IV-A3"
merdf6$types[merdf6$types=="IV-A3|I-E|"]="I-E&IV-A3"
merdf6$types[merdf6$types=="I-E|IV-A3|IV-A3|"]="I-E&IV-A3"
merdf6$types[merdf6$types=="III-A|"]="III-A"
merdf6$types[merdf6$types=="II-A|"]="II-A"
merdf6$types[merdf6$types=="IV-A1|"]="IV-A1"
merdf6$types[merdf6$types=="IV-A2|"]="IV-A2"
merdf6$types[merdf6$types=="IV-A3|"]="IV-A3"
merdf6$types[merdf6$types==" NA " ]="no"





merdf6$types<-as.factor(merdf6$types)




merdf6$gba_acr_inhibit1<-merdf6$high_condence_acr_GBA!="no"
merdf6$acr_inhibit1<-merdf6$gba_acr_inhibit1==TRUE | merdf6$homo_acr_inhibit=="TRUE"
merdf6$acr_crispr1[merdf6$ acr_inhibit1==TRUE & merdf6$CRISPRCas=="yes"]="CRISPRCas+Acr+"
merdf6$acr_crispr1[merdf6$ acr_inhibit1==FALSE & merdf6$CRISPRCas=="yes"]="CRISPRCas+Acr-"
merdf6$acr_crispr1[merdf6$CRISPRCas=="no"]="no"


spename<-spename[spename!="Bordetella_pertussis"]

for (k in c(1,2,3,5,6,7,9)){
  
  subdf<-merdf6[merdf6$species_name==spename[k],]
  subdf<-droplevels(subdf)
 
  a1<-as.data.frame(summary(subdf$types))
 # subdf<-subdf[!is.na(subdf$len_chro),]
  #subdf<-subdf[subdf$len_chro>1,]
 ##### typer<-row.names(a1)[a1$`summary(subdf$types)`>20]
  subdf<-subdf[subdf$types %in% typer,]
  #subdf<-subdf[subdf$acr_crispr!="CRISPRCas+Acr+" & subdf$CRISPRCas_count<3,]
  #subdf<-subdf[subdf$n_plasmid==0,]
  subdf<-droplevels(subdf)
  #subdf<-subdf[abs(subdf$size-as.numeric(subdf$len_all)/1000000)<0.1,]
  
  df_p_val <- subdf %>% 
    rstatix::wilcox_test(size ~ types) %>% 
    rstatix::add_xy_position()
  
  
  bb1<-ggboxplot(subdf, x = "types", y = "size",
                 color = "types", palette = "jco")+ggtitle(spename1[k])+             
    add_pvalue(df_p_val,label.size = 3,bracket.size = 0.3)+theme(plot.title = element_text(hjust = 0.5,size=12))+xlab("CRISPRCas_type")+ theme(legend.position = "right")
  
  
  bbplot_lst[[k]]=bb1
  
  
  
  a1<-as.data.frame(summary(subdf$types))
  typer<-row.names(a1)[a1$`summary(subdf$types)`>0.05*nrow(subdf)]
  subdf<-subdf[subdf$types %in% typer,]
  subdf<-droplevels(subdf)
  typer<-typer[typer!="no"]
  subdf<-subdf[subdf$acr_crispr!="no" & subdf$CRISPRCas_count<3,]
 
  for (tt in typer){
    subdf<-subdf[subdf$types==tt,]
    
    bb1<-ggboxplot(subdf, x = "acr_crispr", y = "size",
                   color = "acr_crispr", palette = "jco")+ggtitle(spename1[k])+stat_compare_means(method = "wilcox.test") +theme(plot.title = element_text(hjust = 0.5,size=12))+xlab("CRISPRCas_type")+ theme(legend.position = "right")
    
    print(bb1)
  }
  
  
  
  
  
}
cowplot::plot_grid(plotlist = bbplot_lst[c(1,2,3,5,6,7,9)], nrow = 2)




CRISPRCas_count<-rep(0,nrow(merdf6))
for (i in 1:nrow(merdf6)) {
  CRISPRCas_count[i]<-length(strsplit(as.character(merdf6$type[i]),"-")[[1]])-1
  
}
merdf6$CRISPRCas_count<-CRISPRCas_count


j=1
for (k in c(1,2,3,5,6,7,9)){
  
  subdf<-merdf6[merdf6$species_name==spename[k],]
  
  subdf<-subdf[subdf$acr_crispr!="no" & subdf$CRISPRCas_count==1,]
  a1<-as.data.frame(summary(subdf$types))
  typer<-row.names(a1)[a1$`summary(subdf$types)`>20]
  subdf<-subdf[subdf$types %in% typer,]
  subdf<-droplevels(subdf)
  
  typer<-typer[typer!="no"]
  
  for (tt in typer){
    af<-subdf[subdf$types==tt,]
    
    bb1<-ggboxplot(af, x = "acr_crispr", y = "size",
                   color = "acr_crispr", palette = "jco")+ggtitle(paste(spename1[k],tt))+stat_compare_means(method = "wilcox.test") +theme(plot.title = element_text(hjust = 0.5,size=12))+xlab("CRISPRCas_type")+ theme(legend.position = "right")
    
   # bbplot_lst[[j]]=bb1
    
    
    ##print(bb1)
    p1=ggplot(af, aes(x=size,fill=acr_crispr,color=acr_crispr))  + geom_density(alpha=.4) +theme_classic()+ggtitle(paste(spename1[k],tt))+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12)) 
   print(p1)
   
   bbplot_lst[[j]]=p1
   j=j+1
    }

  
}
cowplot::plot_grid(plotlist = bbplot_lst[1:6], nrow = 2)


cowplot::plot_grid(plotlist = bbplot_lst[7:12], nrow = 2)


