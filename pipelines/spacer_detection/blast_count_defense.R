
for (i in 1:2){
l1=snakemake@input[[i]]

setwd("~/Desktop/master_thesis/results/spacer_detection/")
#l1="results/defense1/Escherichia_coli"




blasttype<-strsplit(l1,"/")[[1]][2]

allfiles<-list.files(l1)
fileout<-snakemake@output[[1]]
#fileout<-"test.txt"

spit<-strsplit(l1,"/")
spename<-spit[[1]][length(spit[[1]])]

  
for (name in allfiles){

fileloc<-paste(l1,name,sep = "/")
dd<-read.csv(fileloc)

df1<-read.csv(fileloc,header = FALSE,sep ="\t")
colnames(df1)<-c("qseqid","sseqid","pident","length","mismatch","qstart","qend","evalue","qcovs")
###set threshold
sdf<-df1[df1$pident>95&df1$qcovs>95,]
sdf$sseqid<-as.character(sdf$sseqid)

qq<-unique(sdf$qseqid)
nu<-length(unique(sdf$qseqid))
all<-length(unique(df1$qseqid))
sseqid<-paste(sdf$sseqid,collapse  = "|")
tars<-"no"
if (nu>0){
ssids<-rep("no",nu)
tar<-rep("no",nu)
for (ss in 1:nu){
  
  adf<-sdf[sdf$qseqid==qq[ss],]
  a1df<-adf[adf$evalue==min(adf$evalue),][1,]
  ssids[ss]<-a1df$sseqid
  tar[ss]<-strsplit(a1df$sseqid,",")[[1]][2]
}
sseqid<-paste(ssids,collapse  = "|")
tars<-paste(tar,collapse  = "|")
}

id<-strsplit(name,".txt")[[1]]
results<-paste(id,nu,all,spename,blasttype,sseqid,tars,sep = ";")

write(results, file = fileout,append = TRUE)
}
}