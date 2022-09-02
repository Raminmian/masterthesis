



setwd("~/Desktop/master_thesis/results/spacer_detection")
merdf6<-read.csv("merdf6.csv")
merdf6<-merdf6[colnames(merdf6!="X.1")]

##################merge

ID<-rownames(df2)
merdf6$species_name<-as.character(merdf6$species_name)
spename1<-unique(merdf6$species_name)

for (i in 1:length(ID)){ID[i]=strsplit(ID[i],"_genomic")[[1]]}
systems<-colnames(df2)
df2$ID<-ID


amerd<-merge(merdf6,df2,by="ID",all.x = TRUE)
amerd<-amerd[!is.na(amerd$CRISPR.Cas),]




##################genome size/accout
dege<-read.csv("~/Desktop/master_thesis/padloc/per_size/all_defense_size.csv",header = FALSE)

colnames(dege)<-c("ID","defense_size")
amerd1<-merge(amerd,dege,by="ID",all.x = TRUE)
amerd1<-amerd1[!is.na(amerd1$defense_size),]
amerd1$per_densen_size<-100*(amerd1$defense_size/(amerd1$size*1000000))
amerd1$active_CRISPRCas<-amerd1$acr_crispr=="CRISPRCas+Acr-"

ppd1<-amerd1[amerd1$species_name!="Bordetella_pertussis" &amerd1$species_name!="Mycobacterium_tuberculosis", ]
cdat <- ddply(amerd1, "species_name", summarise, size.mean=mean(per_densen_size))
ggplot(ppd1,aes(x=per_densen_size,fill=species_name))  + geom_density(alpha=.5) +xlab("defense genes size/genome size (%)")+
  geom_vline(data=cdat, aes(xintercept=size.mean,  colour=species_name),
            linetype="dashed", size=1)+
  facet_wrap(~species_name)+
  theme_classic()+ggtitle("Distribution of percentage of defense genes size")+theme(plot.title = element_text(hjust = 0.5,size=12))
                                                                                                                                                                                       
j=1
for (i in spename){
  subdf<-amerd1[amerd1$species_name==i,]

  cdat <- ddply(subdf, "active_CRISPRCas", summarise, size.mean=mean(per_densen_size))
  
  p1=ggplot(subdf,aes(x=per_densen_size,fill=active_CRISPRCas))  + geom_density(alpha=.5) +
    theme_classic()+ggtitle(i)+theme(plot.title = element_text(hjust = 0.5,size=10))+geom_vline(data=cdat, aes(xintercept=size.mean,  colour=active_CRISPRCas),
                                                                                               linetype="dashed", size=1)
  
  
  b1<-ggboxplot(subdf, x = "active_CRISPRCas", y = "per_densen_size",
                color = "active_CRISPRCas", palette = "jco")+ggtitle(i)+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position = "right")
  
  
  
  plot_lst1[[j]]=b1
  

  j=j+1
}



cowplot::plot_grid(plotlist = plot_lst1[c(1,2,3,4,5,7,8,9)], nrow = 2,labels = c(1:10))





sum(amerd$CRISPRCas_count==amerd$CRISPR.Cas)/nrow(amerd)



systems<-systems[!systems %in% c("ID","systems_count","family_count","spename")]


###################correlation between CRISPRCas and other defene systems

cor_cris<-as.data.frame(matrix(ncol=58))
colnames(cor_cris)<-systems
j=1
for (i in 1:nrow(amerd)){
  if(amerd$acr_crispr[i]=="CRISPRCas+Acr+"){
    amerd$CRISPR.Cas[i]=0}
  
  
}
for(i in c(1,2,3,4,5,7,8)){
  
  subdf<-amerd[amerd$species_name==spename[i],]
  sysdf<-subdf[systems]
  subdf=subdf[subdf$CRISPR.Cas!=0,]
  print(spename[i])
  print(cor.test(subdf$qatABCD,subdf$gabija,method = "kendall"))
  
  
}


selsys<-colnames(cor_cris)[colSums(is.na(cor_cris))<3]
heatdf<-cor_cris[selsys]
heatdf$id<-row.names(cor_cris)


mm <- melt(heatdf,id="id")
colnames(mm)<-c("species","systems","value")
ggplot(mm, aes(x = systems, y =species, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +geom_text(aes(label = round(value, 2)),size=3) +
  theme(axis.text.x=element_text(angle=90, hjust=1))+ theme(axis.text.y=element_text(angle=0, hjust=1,size = 8))












if(amerd$acr_crispr[i]=="CRISPRCas+Acr+"){
  amerd$CRISPR.Cas[i]=0}




#######################CRISPR and RM
for(i in 1:10){
  
  subdf<-amerd[amerd$species_name==spename1[i],]
  subdf$active_CRISPR<-subdf$acr_crispr=="CRISPRCas+Acr-"
  subdf<-droplevels(subdf)
  
  bb1<-ggboxplot(subdf, x = "active_CRISPR", y = "kiwa",
                 color = "active_CRISPR", palette = "jco")+ggtitle(spename1[i])+stat_compare_means(method = "wilcox.test") +theme(plot.title = element_text(hjust = 0.5,size=12))+xlab("CRISPRCas")+ theme(legend.position = "right")
  
  
  
  bbplot_lst[[i]]=bb1
  
}
cowplot::plot_grid(plotlist = bbplot_lst[c(1,2,3,4,5,6,7,9)], nrow = 2)



#################################CRISPRCas and diversity and abudnace of defense systems

j=1
for (i in spename){
  subdf<-amerd1[amerd1$species_name==i,]
  

  b1<-ggboxplot(subdf, x = "active_CRISPRCas", y = "systems_count",
                color = "active_CRISPRCas", palette = "jco")+ggtitle(i)+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position = "right")
  
  
  
  plot_lst1[[j]]=b1
  
  
  j=j+1
}

cowplot::plot_grid(plotlist = plot_lst1[c(1,2,3,4,5,7,8,9)], nrow = 2)


j=1
for (i in spename){
  subdf<-amerd1[amerd1$species_name==i,]
  
  
  b1<-ggboxplot(subdf, x = "active_CRISPRCas", y = "family_count",
                color = "active_CRISPRCas", palette = "jco")+ggtitle(i)+ theme(plot.title = element_text(size=14))+stat_compare_means(method = "wilcox.test")+  theme(legend.text = element_text(size=10))+ theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position = "right")
  
  
  
  plot_lst1[[j]]=b1
  
  
  j=j+1
}

cowplot::plot_grid(plotlist = plot_lst1[c(1,2,3,4,5,7,8,9)], nrow = 2)

####################################################################################



amerd$CRISPR.Cas<-amerd$CRISPR.Cas>0

amerd$CRISPR.Cas[amerd$CRISPR.Cas==TRUE]="yes"
amerd$CRISPR.Cas[amerd$CRISPR.Cas==FALSE]="no"

sum(amerd$CRISPRCas==amerd$CRISPR.Cas)/nrow(amerd)


aright<-amerd[amerd$CRISPRCas!=amerd$CRISPR.Cas,]

summary(as.factor(aright$species_name))


###########do crisprcas influence the system number 

spename=spename1[spename1!="Mycobacterium_tuberculosis" & spename1!="Bordetella_pertussis"]
for(i in 1:8){
 
  subdf<-amerd[amerd$species_name==spename[i],]
  subdf$active_CRISPR<-subdf$acr_crispr=="CRISPRCas+Acr-"
  subdf<-droplevels(subdf)
  
  bb1<-ggboxplot(subdf, x = "acr_crispr", y = "family_count",
                color = "acr_crispr", palette = "jco")+ggtitle(spename[i])+stat_compare_means(method = "wilcox.test") +theme(plot.title = element_text(hjust = 0.5,size=12))+xlab("CRISPRCas")+ theme(legend.position = "right")
  
  
  
 # corr <- round(cor(mtcars), 1) 
#  ggcorrplot(corr)
  #cor_pmat(mtcars)
  bbplot_lst[[i]]=bb1
  
}
cowplot::plot_grid(plotlist = bbplot_lst[1:8], nrow = 2)


for(i in 1:10){
  
  subdf<-amerd[amerd$species_name==spename1[i],]
  subdf$active_CRISPR<-subdf$acr_crispr=="CRISPRCas+Acr-"
  subdf<-droplevels(subdf)
  
  bb1<-ggscatter(amerd, x = "size", y = "systems_count", size=0.2,
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                color = "species_name", palette = "jco",
                facet.by = "species_name", #scales = "free_x",
                xlab = "size", ylab = "systems_count")+ggtitle(spename1[i])
 
  # corr <- round(cor(mtcars), 1) 
  #  ggcorrplot(corr)
  #cor_pmat(mtcars)
  bbplot_lst[[i]]=bb1
  
}
cowplot::plot_grid(plotlist = bbplot_lst[1:10], nrow = 2)

if(amerd$acr_crispr[i]=="CRISPRCas+Acr+"){
  amerd$CRISPR.Cas[i]=0}


for(i in 1:10){
  
  subdf<-amerd[amerd$species_name==spename[i],]
  subdf$active_CRISPR<-subdf$acr_crispr=="CRISPRCas+Acr-"
  subdf<-droplevels(subdf)
  
  bb1<-ggscatter(subdf, x = "size", y = "systems_count", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "genome size", ylab = "systems_count")+ggtitle(spenamea[i])
  
 
  # corr <- round(cor(mtcars), 1) 
  #  ggcorrplot(corr)
  #cor_pmat(mtcars)
  bbplot_lst[[i]]=bb1
  
}
cowplot::plot_grid(plotlist = bbplot_lst[1:10], nrow = 2)


systems<-systems[!systems %in% c("ID","systems_count","family_count","spename")]


cor_cris<-as.data.frame(matrix(nco=58))
colnames(cor_cris)<-systems





bb1<-ggscatter(amerd, x = "size", y = "family_count", size=0.2,
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               color = "species_name", palette = "jco",
               facet.by = "species_name", #scales = "free_x",
               xlab = "size", ylab = "systems_count")+ggtitle("Genome_size ~ total_systems")+theme(plot.title = element_text(hjust = 0.5,size=12))+ scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))




bb1<-ggscatter(amerd, x = "size", y = "family_count", size=0.2,
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "size", ylab = "family_count")+ggtitle("Genome_size ~ total_families")+theme(plot.title = element_text(hjust = 0.5,size=12))



bb1<-ggscatter(amerd, x = "size", y = "systems_count", size=0.2,
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "size", ylab = "systems_count")+ggtitle("Genome_size ~ total_systems")+theme(plot.title = element_text(hjust = 0.5,size=12))









