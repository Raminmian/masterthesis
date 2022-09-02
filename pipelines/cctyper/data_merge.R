

setwd("~/Desktop/master_thesis/results/cctyper")


#file1="cctyper_result/Escherichia_coli/summary.txt"

sum="summary.txt"
file1=paste(snakemake@input[[1]],sum,sep = "/")

file2<-snakemake@input[[2]]
file_out<-snakemake@output[[1]]
df1<-read.csv(file1,header = TRUE,sep = ";")
df1$CRISPRCas[df1$CRISPRCas!=0]="yes"
df1$CRISPRCas[df1$CRISPRCas==0]="no"
df1$Cas[df1$Cas!=0]="yes"
df1$Cas[df1$Cas==0]="no"
df1$CRISPR [df1$CRISPR!=0]="yes"
df1$CRISPR [df1$CRISPR==0]="no"





#read file containing genome size and gc contents
read_size<-function(file2){
df2<-read.csv(file2)
df2=df2[c("RefSeq.FTP","Size.Mb.","GC.","X.Organism.Name")]
colnames(df2)<-c("ID","size","GC","species_name")
df2$ID=as.character(df2$ID)
for (i in 1:nrow(df2)){
a1=df2[i,1]
a2=strsplit(a1,"/")[[1]]
df2[i,1]=a2[length(a2)]
}
return(df2)
}

dff<-read_size(file2)
final_df<-merge(df1, dff, by="ID", all.x =FALSE)
write.csv(final_df,file_out)



