
for (i in 1:3){
l1=snakemake@input[[i]]

setwd("~/Desktop/master_thesis/results/spacer_detection/")
#l1="results/all_spacers/blast_plasmid/Staphylococcus_aureus"
#l2<-"spacers/spacer_nearcas_list/Staphylococcus_aureus.txt"
l2<-snakemake@input[[4]]

alls<-read.csv(l2,header=FALSE,sep = ";")
blasttype<-strsplit(l1,"/")[[1]][3]

allfiles<-list.files(l1)
fileout<-snakemake@output[[1]]
#fileout<-"test.txt"

spit<-strsplit(l1,"/")
spename<-spit[[1]][length(spit[[1]])]

  
for (name in allfiles){
fileloc<-paste(l1,name,sep = "/")
df1<-read.csv(fileloc,header = FALSE,sep ="\t")
colnames(df1)<-c("qseqid","sseqid","pident","length","mismatch","qstart","qend","evalue","qcovs")


df2<-df1[df1$pident>95&df1$qcovs>95,]


############to know if it is spacers from orphan crispr array or not

id<-strsplit(name,".txt")[[1]]
alls$V2<-as.character(alls$V2)
newalls<-alls[alls$V1==id,]
if (nrow(newalls)!=0){
  newalls$V2<-as.character(newalls$V2)
  isin<-rep(FALSE,nrow(df2))
  df2$qseqid<-as.character(df2$qseqid)
  for (k in 1:nrow(df2)){
    isin[k]=df2$qseqid[k] %in% strsplit(newalls$V2[1]," ")[[1]]
    
  }
  
  df2<-df2[isin,]
  ###set threshold
  
  nu<-length(unique(df2$qseqid))
  all<-length(strsplit(newalls$V2[1]," ")[[1]])
  
  
  results<-paste(id,nu,all,spename,blasttype,sep = ";")
  #write(results,file="test.txt",append = TRUE)
  write(results, file = fileout,append = TRUE)
}


}
}