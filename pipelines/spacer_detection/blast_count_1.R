
for (i in 1:3){
l1=snakemake@input[[i]]

setwd("~/Desktop/master_thesis/results/spacer_detection/")
#l1="results/all_spacers/blast_plasmid/Acinetobacter_baumannii"

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
###set threshold
nu<-length(unique(df1[df1$pident>95&df1$qcovs>95,]$qseqid))
all<-length(unique(df1$qseqid))

id<-strsplit(name,".txt")[[1]]
results<-paste(id,nu,all,spename,blasttype,sep = ";")
write(results, file = fileout,append = TRUE)
}
}