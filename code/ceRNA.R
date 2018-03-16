setwd("C:/projects/ceRNA")
cancergene=read.csv("Census_allSun Nov 13 17-37-13 2016.tsv",stringsAsFactors=F,sep="\t",skip=0,header = T)
cancergene=unique(cancergene$Gene.Symbol)
net=read.csv("HI3_Y2H_102416.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
CGC_neig=c()
for (i in 1:dim(net)[1]) {
  x=which(cancergene==net$Symbol.for.A[i])
  y=which(cancergene==net$Symbol.for.B[i])
  if(length(x)>0|length(y)>0){
    CGC_neig=rbind(CGC_neig,net[i,])
  }
}
CGC_neig_gene=union(CGC_neig$Symbol.for.A,CGC_neig$Symbol.for.B)
Finalgene=union(CGC_neig_gene,cancergene)
write.table(Finalgene,file = "Finalgene.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

###################################
miRNA_gene=read.csv("miRNA_gene_ensem.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
miRNA_lncRNA_MCF=read.csv("miRNA_lncRNA_Hela.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
miRNA_reg=rbind(miRNA_gene,miRNA_lncRNA_MCF)
gene=unique(miRNA_reg$V2)
N=386
ceRNA_candidate=c()
for (i in 1:(length(gene)-1)) {
  for (j in (i+1):length(gene)) {
    print(c(i,j))
    x1=which(miRNA_reg$V2==gene[i])
    x2=which(miRNA_reg$V2==gene[j])
    miRNA1=unique(miRNA_reg$V1[x1])
    miRNA2=unique(miRNA_reg$V1[x2])
    n1=length(miRNA1)
    n2=length(miRNA2)
    overmiRNA=intersect(miRNA1,miRNA2)
    n=length(overmiRNA)
    if(n>=3){
      P=phyper(n-1,n1,N-n1,n2,lower.tail=FALSE, log.p = FALSE)
      KK=cbind(gene[i],gene[j],n1,n2,n,P)
      ceRNA_candidate=rbind(ceRNA_candidate,KK)
    }
  }
}
FDR=p.adjust(ceRNA_candidate[,6],method = "fdr")
ceRNA_candidate=cbind(ceRNA_candidate,FDR)
sigce=which(FDR<0.05)
MCF_ceRNA=ceRNA_candidate[sigce,]
write.table(MCF_ceRNA,file = "Hela_ceRNA.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)


#####################Download the gene expression from TCGA
library(TCGAbiolinks)
# Downloading and prepare 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
BRCAdata <- GDCprepare(query, save = TRUE, 
                       save.filename = "BRCAExpression.rda",
                       remove.files.prepared = TRUE)
library(SummarizedExperiment)
BRCAmatrix<-assay(data,1,"FPKM")

gene=rownames(BRCAmatrix)

############3
BRCAclinical=read.csv("BRCA_clinical.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
x=which(BRCAclinical$ER.Status=="Positive"&BRCAclinical$PR.Status=="Positive"&BRCAclinical$HER2.Final.Status=="Negative")
mcf_sample=BRCAclinical[x,]
x=which(BRCAclinical$ER.Status=="Negative"&BRCAclinical$PR.Status=="Negative"&BRCAclinical$HER2.Final.Status=="Negative")
mda_sample=BRCAclinical[x,]
###########load 
load("C:/projects/ceRNA/BRCAExpression.rda")
library(SummarizedExperiment)
BRCAmatrix<-assay(data,1,"FPKM")
allBRCAsample=colnames(BRCAmatrix)
C_N=data@colData@listData$definition#############sample information tumor vs normal
Normal_s=which(C_N=="Solid Tissue Normal")
Cancer_s=which(C_N=="Primary solid Tumor"|C_N=="Metastatic")
BRCAthree=c()
for (i in 1:length(allBRCAsample)) {
  aa=strsplit(allBRCAsample[i],"-")
  bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
  BRCAthree=rbind(BRCAthree,bb)
}
Cancer_s=cbind(BRCAthree[Cancer_s],Cancer_s)
Normal_s=cbind(BRCAthree[Normal_s],Normal_s)
inter_mcf=merge(mcf_sample,Cancer_s,by.x="Complete.TCGA.ID",by.y="V1")
inter_mda=merge(mda_sample,Cancer_s,by.x="Complete.TCGA.ID",by.y="V1")

###########################################co-expression analysis
MCF_ceRNA=read.csv("MCF_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)

##########MCF-7
library(Hmisc)
allgene=rownames(BRCAmatrix)
R_ceRNA_MCF=c()
for (i in 1:dim(MCF_ceRNA)[1]) {
  print(i)
  xxa=which(allgene==MCF_ceRNA$V1[i])
  xxb=which(allgene==MCF_ceRNA$V2[i])
  if(length(xxa)>0&length(xxb)>0){
    expa=BRCAmatrix[xxa,inter_mcf$Cancer_s]
    expb=BRCAmatrix[xxb,inter_mcf$Cancer_s]
    expa=log2(as.numeric(as.character(expa)))
    expb=log2(as.numeric(as.character(expb)))
    R=cor(expa, expb, method = "spearman")
    pp=cor.test(expa, expb, method = "spearman")$p.value
    R_ceRNA_MCF=rbind(R_ceRNA_MCF,cbind(MCF_ceRNA[i,],R,pp))
  }
}
FDR=p.adjust(R_ceRNA_MCF$pp,method = "BH")
R_ceRNA_MCF=cbind(R_ceRNA_MCF,FDR)
Sig_index=which(R_ceRNA_MCF$R>0&R_ceRNA_MCF$FDR<0.01)
write.table(R_ceRNA_MCF[Sig_index,],file = "MCF_ceRNA_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
####################MDA
MDA_ceRNA=read.csv("MDA_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
R_ceRNA_MDA=c()
for (i in 1:dim(MDA_ceRNA)[1]) {
  print(i)
  xxa=which(allgene==MDA_ceRNA$V1[i])
  xxb=which(allgene==MDA_ceRNA$V2[i])
  if(length(xxa)>0&length(xxb)>0){
    expa=BRCAmatrix[xxa,inter_mda$Cancer_s]
    expb=BRCAmatrix[xxb,inter_mda$Cancer_s]
    expa=log2(as.numeric(as.character(expa)))
    expb=log2(as.numeric(as.character(expb)))
    R=cor(expa, expb, method = "spearman")
    pp=cor.test(expa, expb, method = "spearman")$p.value
    R_ceRNA_MDA=rbind(R_ceRNA_MDA,cbind(MDA_ceRNA[i,],R,pp))
  }
}
FDR=p.adjust(R_ceRNA_MDA$pp,method = "BH")
R_ceRNA_MDA=cbind(R_ceRNA_MDA,FDR)
Sig_index=which(R_ceRNA_MDA$R>0&R_ceRNA_MDA$FDR<0.01)
write.table(R_ceRNA_MDA[Sig_index,],file = "MDA_ceRNA_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

######################################################## K562 analysis
library(TCGAbiolinks)
# Downloading and prepare 
query <- GDCquery(project = "TCGA-LAML", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
LAMLdata <- GDCprepare(query, save = TRUE, 
                       save.filename = "LAMLExpression.rda",
                       remove.files.prepared = TRUE)
library(SummarizedExperiment)
LAMLmatrix<-assay(LAMLdata,1,"FPKM")
LAMLclinal<-LAMLdata@colData@listData$definition#############sample information tumor vs normal
LAMLclinal=t(t(LAMLclinal)) #############NO normal samples
allgene=rownames(LAMLmatrix)
K562_ceRNA=read.csv("K562_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
R_ceRNA_K562=c()
for (i in 1:dim(K562_ceRNA)[1]) {
  print(i)
  xxa=which(allgene==K562_ceRNA$V1[i])
  xxb=which(allgene==K562_ceRNA$V2[i])
  if(length(xxa)>0&length(xxb)>0){
    expa=LAMLmatrix[xxa,]
    expb=LAMLmatrix[xxb,]
    expa=log2(as.numeric(as.character(expa)))
    expb=log2(as.numeric(as.character(expb)))
    R=cor(expa, expb, method = "spearman")
    pp=cor.test(expa, expb, method = "spearman")$p.value
    R_ceRNA_K562=rbind(R_ceRNA_K562,cbind(K562_ceRNA[i,],R,pp))
  }
}
FDR=p.adjust(R_ceRNA_K562$pp,method = "BH")
R_ceRNA_K562=cbind(R_ceRNA_K562,FDR)
Sig_index=which(R_ceRNA_K562$R>0&R_ceRNA_K562$FDR<0.01)
write.table(R_ceRNA_K562[Sig_index,],file = "K562_ceRNA_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

######################################################## Hela analysis

library(TCGAbiolinks)
# Downloading and prepare 
query <- GDCquery(project = "TCGA-UCEC", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
UCECdata <- GDCprepare(query, save = TRUE, 
                       save.filename = "UCECExpression.rda",
                       remove.files.prepared = TRUE)
library(SummarizedExperiment)
UCECmatrix<-assay(UCECdata,1,"FPKM")
UCECclinal<-UCECdata@colData@listData$definition#############sample information tumor vs normal
UCECclinal=t(t(UCECclinal)) #############NO normal samples
Cancer_s=which(UCECclinal=="Primary solid Tumor"|UCECclinal=="Recurrent Solid Tumor")
allgene=rownames(UCECmatrix)
Hela_ceRNA=read.csv("Hela_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
R_ceRNA_Hela=c()
for (i in 1:dim(Hela_ceRNA)[1]) {
  print(i)
  xxa=which(allgene==Hela_ceRNA$V1[i])
  xxb=which(allgene==Hela_ceRNA$V2[i])
  if(length(xxa)>0&length(xxb)>0){
    expa=UCECmatrix[xxa,Cancer_s]
    expb=UCECmatrix[xxb,Cancer_s]
    expa=log2(as.numeric(as.character(expa)))
    expb=log2(as.numeric(as.character(expb)))
    R=cor(expa, expb, method = "spearman")
    pp=cor.test(expa, expb, method = "spearman")$p.value
    R_ceRNA_Hela=rbind(R_ceRNA_Hela,cbind(Hela_ceRNA[i,],R,pp))
  }
}
FDR=p.adjust(R_ceRNA_Hela$pp,method = "BH")
R_ceRNA_Hela=cbind(R_ceRNA_Hela,FDR)
Sig_index=which(R_ceRNA_Hela$R>0&R_ceRNA_Hela$FDR<0.01)
write.table(R_ceRNA_Hela[Sig_index,],file = "Hela_ceRNA_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

######################################################## U87 analysis

library(TCGAbiolinks)
# Downloading and prepare 
query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
GBMdata <- GDCprepare(query, save = TRUE, 
                       save.filename = "GBMExpression.rda",
                       remove.files.prepared = TRUE)
library(SummarizedExperiment)
GBMmatrix<-assay(GBMdata,1,"FPKM")
GBMclinal<-GBMdata@colData@listData$definition#############sample information tumor vs normal
GBMclinal=t(t(GBMclinal)) #############NO normal samples
Cancer_s=which(GBMclinal=="Primary solid Tumor"|GBMclinal=="Recurrent Solid Tumor")
allgene=rownames(GBMmatrix)
U87_ceRNA=read.csv("U87_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
R_ceRNA_U87=c()
for (i in 1:dim(U87_ceRNA)[1]) {
  print(i)
  xxa=which(allgene==U87_ceRNA$V1[i])
  xxb=which(allgene==U87_ceRNA$V2[i])
  if(length(xxa)>0&length(xxb)>0){
    expa=GBMmatrix[xxa,Cancer_s]
    expb=GBMmatrix[xxb,Cancer_s]
    expa=log2(as.numeric(as.character(expa)))
    expb=log2(as.numeric(as.character(expb)))
    R=cor(expa, expb, method = "spearman")
    pp=cor.test(expa, expb, method = "spearman")$p.value
    R_ceRNA_U87=rbind(R_ceRNA_U87,cbind(U87_ceRNA[i,],R,pp))
  }
}
FDR=p.adjust(R_ceRNA_U87$pp,method = "BH")
R_ceRNA_U87=cbind(R_ceRNA_U87,FDR)
Sig_index=which(R_ceRNA_U87$R>0&R_ceRNA_U87$FDR<0.01)
write.table(R_ceRNA_U87[Sig_index,],file = "U87_ceRNA_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

########################################### candidate ceRNA overlap analysis
MCF7_ceRNA=read.csv("MCF_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
MDA_ceRNA=read.csv("MDA_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
K562_ceRNA=read.csv("K562_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
Hela_ceRNA=read.csv("Hela_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
U87_ceRNA=read.csv("U87_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)

MCF7_ceRNA_sorted=c()
for (i in 1:dim(MCF7_ceRNA)[1]) {
  print(i)
  if((MCF7_ceRNA$V1[i]>MCF7_ceRNA$V2[i])==TRUE){
    MCF7_ceRNA_sorted=rbind(MCF7_ceRNA_sorted,MCF7_ceRNA[i,c(1,2)])
  } else{
    MCF7_ceRNA_sorted=rbind(MCF7_ceRNA_sorted,MCF7_ceRNA[i,c(2,1)])
  }
}

MDA_ceRNA_sorted=c()
for (i in 1:dim(MDA_ceRNA)[1]) {
  print(i)
  if((MDA_ceRNA$V1[i]>MDA_ceRNA$V2[i])==TRUE){
    MDA_ceRNA_sorted=rbind(MDA_ceRNA_sorted,MDA_ceRNA[i,c(1,2)])
  } else{
    MDA_ceRNA_sorted=rbind(MDA_ceRNA_sorted,MDA_ceRNA[i,c(2,1)])
  }
}

K562_ceRNA_sorted=c()
for (i in 1:dim(K562_ceRNA)[1]) {
  print(i)
  if((K562_ceRNA$V1[i]>K562_ceRNA$V2[i])==TRUE){
    K562_ceRNA_sorted=rbind(K562_ceRNA_sorted,K562_ceRNA[i,c(1,2)])
  } else{
    K562_ceRNA_sorted=rbind(K562_ceRNA_sorted,K562_ceRNA[i,c(2,1)])
  }
}

Hela_ceRNA_sorted=c()
for (i in 1:dim(Hela_ceRNA)[1]) {
  print(i)
  if((Hela_ceRNA$V1[i]>Hela_ceRNA$V2[i])==TRUE){
    Hela_ceRNA_sorted=rbind(Hela_ceRNA_sorted,Hela_ceRNA[i,c(1,2)])
  } else{
    Hela_ceRNA_sorted=rbind(Hela_ceRNA_sorted,Hela_ceRNA[i,c(2,1)])
  }
}

U87_ceRNA_sorted=c()
for (i in 1:dim(U87_ceRNA)[1]) {
  print(i)
  if((U87_ceRNA$V1[i]>U87_ceRNA$V2[i])==TRUE){
    U87_ceRNA_sorted=rbind(U87_ceRNA_sorted,U87_ceRNA[i,c(1,2)])
  } else{
    U87_ceRNA_sorted=rbind(U87_ceRNA_sorted,U87_ceRNA[i,c(2,1)])
  }
}

write.table(MCF7_ceRNA_sorted,file = "MCF7_ceRNA_sorted_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(MDA_ceRNA_sorted,file = "MDA_ceRNA_sorted_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(K562_ceRNA_sorted,file = "K562_ceRNA_sorted_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(Hela_ceRNA_sorted,file = "Hela_ceRNA_sorted_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(U87_ceRNA_sorted,file = "U87_ceRNA_sorted_final.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

x1=unique(MCF7_ceRNA_sorted$V1)
x2=unique(MCF7_ceRNA_sorted$V2)
MCFgene=union(x1,x2)

x1=unique(MDA_ceRNA_sorted$V1)
x2=unique(MDA_ceRNA_sorted$V2)
MDAgene=union(x1,x2)

x1=unique(K562_ceRNA_sorted$V1)
x2=unique(K562_ceRNA_sorted$V2)
K562gene=union(x1,x2)

x1=unique(Hela_ceRNA_sorted$V1)
x2=unique(Hela_ceRNA_sorted$V2)
Helagene=union(x1,x2)

x1=unique(U87_ceRNA_sorted$V1)
x2=unique(U87_ceRNA_sorted$V2)
U87gene=union(x1,x2)

MCF7_ceRNA=read.csv("MCF_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
MDA_ceRNA=read.csv("MDA_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
K562_ceRNA=read.csv("K562_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
Hela_ceRNA=read.csv("Hela_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
U87_ceRNA=read.csv("U87_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
lncRNA=read.csv("lncRNAwithmiRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)

x1=unique(MCF7_ceRNA$V1)
x2=unique(MCF7_ceRNA$V2)
MCFgene=union(x1,x2)

x1=unique(MDA_ceRNA$V1)
x2=unique(MDA_ceRNA$V2)
MDAgene=union(x1,x2)

x1=unique(K562_ceRNA$V1)
x2=unique(K562_ceRNA$V2)
K562gene=union(x1,x2)

x1=unique(Hela_ceRNA$V1)
x2=unique(Hela_ceRNA$V2)
Helagene=union(x1,x2)

x1=unique(U87_ceRNA$V1)
x2=unique(U87_ceRNA$V2)
U87gene=union(x1,x2)

Conserved=read.csv("conservedceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
x1=unique(Conserved$V1)
x2=unique(Conserved$V2)
Conserved_gene=union(x1,x2)
write.table(Conserved_gene,file = "Conserved_gene.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

#############################################################network analysis
MCF7_ceRNA=read.csv("MCF_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
MDA_ceRNA=read.csv("MDA_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
K562_ceRNA=read.csv("K562_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
Hela_ceRNA=read.csv("Hela_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
U87_ceRNA=read.csv("U87_ceRNA_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
x1=unique(MCF7_ceRNA$V1)
x2=unique(MCF7_ceRNA$V2)
MCFgene=union(x1,x2)

x1=unique(MDA_ceRNA$V1)
x2=unique(MDA_ceRNA$V2)
MDAgene=union(x1,x2)

x1=unique(K562_ceRNA$V1)
x2=unique(K562_ceRNA$V2)
K562gene=union(x1,x2)

x1=unique(Hela_ceRNA$V1)
x2=unique(Hela_ceRNA$V2)
Helagene=union(x1,x2)

x1=unique(U87_ceRNA$V1)
x2=unique(U87_ceRNA$V2)
U87gene=union(x1,x2)

###############degree
degree=c()
for (i in 1:length(MCFgene)) {
  x1=which(MCF7_ceRNA$V1==MCFgene[i])
  x2=which(MCF7_ceRNA$V2==MCFgene[i])
  degree=rbind(degree,length(x1)+length(x2))
}
MCF_degree=cbind(MCFgene,degree)

degree=c()
for (i in 1:length(MDAgene)) {
  x1=which(MDA_ceRNA$V1==MDAgene[i])
  x2=which(MDA_ceRNA$V2==MDAgene[i])
  degree=rbind(degree,length(x1)+length(x2))
}
MDA_degree=cbind(MDAgene,degree)

degree=c()
for (i in 1:length(K562gene)) {
  x1=which(K562_ceRNA$V1==K562gene[i])
  x2=which(K562_ceRNA$V2==K562gene[i])
  degree=rbind(degree,length(x1)+length(x2))
}
K562_degree=cbind(K562gene,degree)

degree=c()
for (i in 1:length(Helagene)) {
  x1=which(Hela_ceRNA$V1==Helagene[i])
  x2=which(Hela_ceRNA$V2==Helagene[i])
  degree=rbind(degree,length(x1)+length(x2))
}
Hela_degree=cbind(Helagene,degree)

degree=c()
for (i in 1:length(U87gene)) {
  x1=which(U87_ceRNA$V1==U87gene[i])
  x2=which(U87_ceRNA$V2==U87gene[i])
  degree=rbind(degree,length(x1)+length(x2))
}
U87_degree=cbind(U87gene,degree)


lncRNA=read.csv("lncRNAwithmiRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
CGC=read.csv("CGC_ens.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)

lgene<-function(lncRNA,CGC,MCF_degree){
  LncM=c()
  for (i in 1:dim(MCF_degree)[1]) {
  x=which(lncRNA$V1==MCF_degree[i,1])
x1=which(CGC$V1==MCF_degree[i,1])
 if(length(x)>0){
  LncM=rbind(LncM,2)
} else if (length(x1)>0){
     LncM=rbind(LncM,1)
           } else{
               LncM=rbind(LncM,0)
             }
     }
  MCF_degree=cbind(MCF_degree,LncM)
  return(MCF_degree)
}


MCF_degree=lgene(lncRNA,CGC,MCF_degree)
MDA_degree=lgene(lncRNA,CGC,MDA_degree)
K562_degree=lgene(lncRNA,CGC,K562_degree)
Hela_degree=lgene(lncRNA,CGC,Hela_degree)
U87_degree=lgene(lncRNA,CGC,U87_degree)

library(vioplot)
par(mfrow=c(3,3))

x1=which(MCF_degree[,3]==0)
Y1=as.numeric(as.character(MCF_degree[x1,2]))

x2=which(MCF_degree[,3]==1)
Y2=as.numeric(as.character(MCF_degree[x2,2]))

x3=which(MCF_degree[,3]==2)
Y3=as.numeric(as.character(MCF_degree[x3,2]))


vioplot(Y1,Y2,Y3)

x1=which(MDA_degree[,3]==0)
Y1=as.numeric(as.character(MDA_degree[x1,2]))

x2=which(MDA_degree[,3]==1)
Y2=as.numeric(as.character(MDA_degree[x2,2]))

x3=which(MDA_degree[,3]==2)
Y3=as.numeric(as.character(MDA_degree[x3,2]))
vioplot(Y1,Y2,Y3)

wilcox.test(Y1,Y2)
wilcox.test(Y1,Y3)
wilcox.test(Y2,Y3)

x1=which(K562_degree[,3]==0)
Y1=as.numeric(as.character(K562_degree[x1,2]))

x2=which(K562_degree[,3]==1)
Y2=as.numeric(as.character(K562_degree[x2,2]))

x3=which(K562_degree[,3]==2)
Y3=as.numeric(as.character(K562_degree[x3,2]))
vioplot(Y1,Y2,Y3)

wilcox.test(Y1,Y2)
wilcox.test(Y1,Y3)
wilcox.test(Y2,Y3)

x1=which(Hela_degree[,3]==0)
Y1=as.numeric(as.character(Hela_degree[x1,2]))

x2=which(Hela_degree[,3]==1)
Y2=as.numeric(as.character(Hela_degree[x2,2]))

x3=which(Hela_degree[,3]==2)
Y3=as.numeric(as.character(Hela_degree[x3,2]))
vioplot(Y1,Y2,Y3)

wilcox.test(Y1,Y2)
wilcox.test(Y1,Y3)
wilcox.test(Y2,Y3)


x1=which(U87_degree[,3]==0)
Y1=as.numeric(as.character(U87_degree[x1,2]))

x2=which(U87_degree[,3]==1)
Y2=as.numeric(as.character(U87_degree[x2,2]))

x3=which(U87_degree[,3]==2)
Y3=as.numeric(as.character(U87_degree[x3,2]))
vioplot(Y1,Y2,Y3)

wilcox.test(Y1,Y2)
wilcox.test(Y1,Y3)
wilcox.test(Y2,Y3)

write.table(MCF_degree,file = "MCF7_ceRNA_degree.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(MDA_degree,file = "MDA_ceRNA_degree.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(K562_degree,file = "K562_ceRNA_degree.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(Hela_degree,file = "Hela_ceRNA_degree.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
write.table(U87_degree,file = "U87_ceRNA_degree.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

#########################hub analysis
conserved_hub=read.csv("conserved_hub.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
cancer=c("MCF","MDA","K562","Hela","U87")
Result=c()
for (i in 1:dim(conserved_hub)[1]) {
  print(i)
  XX=c()
  hub=conserved_hub$V1[i]
  for (j in 1:4) {
    for (k in (j+1):5) {
      net1=read.csv(paste(cancer[j],"_ceRNA_final.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
      net2=read.csv(paste(cancer[k],"_ceRNA_final.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
      xx1=which(net1$V1==hub)
      xx2=which(net1$V2==hub)
      n1=unique(net1$V2[xx1])
      n2=unique(net1$V1[xx2])
      n=union(n1,n2)
      yy1=which(net2$V1==hub)
      yy2=which(net2$V2==hub)
      m1=unique(net2$V2[yy1])
      m2=unique(net2$V1[yy2])
      m=union(m1,m2)
      Over=intersect(m,n)
      XX=cbind(XX,length(Over)/min(length(m),length(n)))
    }
  }
  Result=rbind(Result,XX)
}
library(pheatmap)
rownames(Result)=conserved_hub$V1
colnames(Result)=c(1:10)
pheatmap(Result,cluster_rows = F,show_colnames = T)

image(t(Result),col =terrain.colors(60))
axis(1, 1:10)
axis(2, 1:13)
for (x in 1:13) {
  for (y in 1:10) {
    text(x,y,sprintf("%0.2f",Result[x,y]))
  }
}
  
  

#####################################3mutation analysis
cancer=c("MCF","MDA","K562","Hela","U87")
cancername=c("BRCA","BRCA","LAML","UCEC","GBM")
##########Whether the mutation freq of ceRNA are higher than other genes
P=c()
R=c()
P_cor=c()
par(mfrow=c(2,1))
plot(0,0,type="n",xlim=c(0.5,10.5), ylim=c(0,0.3))
library(vioplot)
m=1
AllX1=c()
AllX2=c()
for (i in 1:5) {
  setwd("C:/projects/ceRNA")
  ceRNA=read.csv(paste(cancer[i],"_degree_symbol.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  setwd("C:/projects/TCGAmutation")
  mutation_freq<-read.csv(paste(cancername[i],"_freq.txt",sep=""),stringsAsFactors=F,sep="\t",skip=1,header = FALSE)
  ceRNA_freq=merge(ceRNA,mutation_freq,by.x = "V5",by.y = "V1")
  allgene=unique(mutation_freq$V1)
  otherg=setdiff(allgene,ceRNA$V5)
  otherg=t(t(otherg))
  other_freq=merge(otherg,mutation_freq,by.x="V1",by.y="V1")
  X1=as.numeric(as.character(ceRNA_freq$V2.y))
  X2=as.numeric(as.character(other_freq$V2))
  X3=as.numeric(as.character(ceRNA_freq$V3))
  p_v=wilcox.test(X1,X2,alternative = "greater")$p.value
  rr=cor(X1,X3)
  R=rbind(R,rr)
  pp=cor.test(X1,X3)$p.value
  P_cor=rbind(P_cor,p_v)
  P=rbind(P,p_v)
  vioplot(X1,at = m, add = T,col ="Darkblue" ,colMed="paleturquoise",wex=1)
  vioplot(X2,at = m+1, add = T,col ="Gray" ,colMed="paleturquoise",wex=1)
  AllX1=rbind(AllX1,X1)
  AllX2=rbind(AllX2,X2)
  m=m+2
}
AllX1=matrix(AllX1,ncol = 1)
AllX2=matrix(AllX2,ncol = 1)
vioplot(AllX1,AllX2)
wilcox.test(AllX1,AllX2,alternative = "greater")
##########################whether these mutations in ceRNAs are functional importance----CADD and consvertion
cancer=c("MCF","MDA","K562","Hela","U87")
cancername=c("BRCA","BRCA","LAML","UCEC","GBM")
##########Whether the mutation freq of ceRNA are higher than other genes
P=c()
#par(mfrow=c(2,1))
#plot(0,0,type="n",xlim=c(0.5,10.5), ylim=c(0,1))
library(vioplot)
m=1
AllCADD=c()
Allrand=c()
for (i in 1:5) {
  setwd("C:/projects/ceRNA")
  ceRNA=read.csv(paste(cancer[i],"_degree_symbol.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  setwd("C:/projects/TCGAmutation")
  mutation_freq<-read.csv(paste(cancername[i],".maf",sep=""),stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
  setwd("C:/projects/Map")
  Mutfun<-read.csv(paste("TCGA.",cancername[i],".maf.hg38_multianno.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
  ceRNA_mut=merge(ceRNA,mutation_freq,by.x = "V5",by.y = "Hugo_Symbol")
  mutscore=merge(ceRNA_mut,Mutfun,by.x = c("Chromosome","Start_Position","End_Position","Tumor_Seq_Allele1","Tumor_Seq_Allele2"),
                 by.y = c("Chr","Start","End","Ref","Alt"))
  cadd_score<-mutscore$integrated_fitCons_score ###real score--------------------------
  cadd_score<-as.numeric(as.character(cadd_score))  
  ks<-which(is.na(cadd_score)==FALSE)
  cadd_score<-cadd_score[ks]
  myn<-length(cadd_score)
  allmutn<-dim(Mutfun)[1]
  RR_score<-c()
  for (j in 1:100) {
    randn<-sample(c(1:allmutn),myn)
    rand_score<-Mutfun$integrated_fitCons_score[randn]
    rand_score<-as.numeric(as.character(rand_score))
    sss<-which(is.na(rand_score)==FALSE)
    rand_score<-rand_score[sss]
    #xp<-wilcox.test(cadd_score,rand_score,alternative = "greater")$p.value
    #p1<-rbind(p1,xp)
    RR_score<-rbind(RR_score,t(t(rand_score)))
  }
  xp<-wilcox.test(cadd_score,RR_score,alternative = "greater")$p.value
  P=rbind(P,xp)
  AllCADD=rbind(AllCADD,cadd_score)
  Allrand=rbind(Allrand,RR_score)
 # vioplot(cadd_score,at = m, add = T,col ="Skyblue" ,colMed="paleturquoise",wex=1)
  #vioplot(RR_score,at = m+1, add = T,col ="Gray" ,colMed="paleturquoise",wex=1)
  m=m+2
}
AllCADD=matrix(AllCADD,ncol = 1)
Allrand=matrix(Allrand,ncol = 1)
#####################################################3
#######################################################
########################functional analysis of lncRNA ceRNA
cancer=c("MCF","MDA","K562","Hela","U87")
cancername=c("BRCA","BRCA","LAML","UCEC","GBM")
All_result=c()
for (i in 1:5) {
  setwd("C:/projects/ceRNA")
  ceRNA=read.csv(paste(cancer[i],"_degree_symbol.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  x=which(ceRNA$V4==2)
  ceRNAnet=read.csv(paste(cancer[i],"_ceRNA_sorted_final.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  ensg2symbol=read.csv("ensmble2symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  lncRNA=ceRNA[x,c(2,5)]
  lncRNA_Result=c()
  for (j in 1:dim(lncRNA)[1]) {
    xx1=which(ceRNAnet$V1==lncRNA$V2[j])
    xx2=which(ceRNAnet$V2==lncRNA$V2[j])
    gene1=unique(ceRNAnet$V2[xx1])
    gene2=unique(ceRNAnet$V1[xx2])
    gene=union(gene2,gene1)
    gene=t(t(gene))
    genename=merge(gene,ensg2symbol,by.x="V1",by.y="V1")
    gene=genename$V2
    Result=c()
    Biofunction <- file("go2hallmark1.txt", "r")
    i<-0
    N <-18999
    line=readLines(Biofunction,n=1)
    while( length(line) != 0 ) {
      i <- i+1
      msig.go <- strsplit(line,"\t")
      msig.go <- msig.go[[1]]
      m <- length(msig.go)
      goterm <- msig.go[1]
      msig.go <- msig.go[3:m]
      n1 <- length(gene)
      intergene <- intersect(msig.go,gene)
      k1 <- length(intergene)
      print(k1)
      if (k1>0){
        Result <- cbind(Result,1)
      } else {
        Result <- cbind(Result,0)
      }
      line=readLines(Biofunction,n=1);
    }
    close(Biofunction)
    WW=cbind(lncRNA[j,],Result)
    names(WW)=c(1:36)
    lncRNA_Result=rbind(lncRNA_Result,WW)
  }
  All_result=rbind(All_result,lncRNA_Result)
}
write.table(All_result,file = "lncRNA_hallmark.txt",sep = "\t",quote = F,row.names = F,col.names = F)

#######################Extract the lncRNA-gene pair in each cancer
cancer=c("MCF","MDA","K562","Hela","U87")
cancername=c("BRCA","BRCA","LAML","UCEC","GBM")
All_result=c()
for (i in 1:5) {
  setwd("C:/projects/ceRNA")
  ceRNA=read.csv(paste(cancer[i],"_degree_symbol.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  x=which(ceRNA$V4==2)
  ceRNAnet=read.csv(paste(cancer[i],"_ceRNA_sorted_final.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  ensg2symbol=read.csv("ensmble2symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  lncRNA=ceRNA[x,c(2,5)]
  Result=c()
  for (j in 1:dim(lncRNA)[1]) {
    xx1=which(ceRNAnet$V1==lncRNA$V2[j])
    xx2=which(ceRNAnet$V2==lncRNA$V2[j])
    gene1=unique(ceRNAnet$V2[xx1])
    gene2=unique(ceRNAnet$V1[xx2])
    gene=union(gene2,gene1)
    gene=t(t(gene))
    genename=merge(gene,ensg2symbol,by.x="V1",by.y="V1")
    gene=genename$V2
    Result=rbind(Result,cbind(lncRNA[j,],gene))
  }
  All_result=rbind(All_result,cbind(cancer[i],Result))
}
##############extract gene-go-hallmark pair
gene_hallmark=c()

Biofunction <- file("go2hallmark1.txt", "r")
line=readLines(Biofunction,n=1)
while( length(line) != 0 ) {
  msig.go <- strsplit(line,"\t")
  msig.go <- msig.go[[1]]
  m <- length(msig.go)
  goterm <- msig.go[1]
  hallmark=msig.go[2]
  msig.go <- msig.go[3:m]
  aaa=c()
  for (i in 1:length(msig.go)) {
    aaa=rbind(aaa,cbind(msig.go[i],goterm,hallmark))
  }
  gene_hallmark=rbind(gene_hallmark,aaa)
  line=readLines(Biofunction,n=1);
}
close(Biofunction)

lncRNA_hallmark=merge(All_result,gene_hallmark,by.x = "gene",by.y = "V1")

write.table(lncRNA_hallmark,file = "lncRNA_hallmark_net.txt",sep = "\t",quote = F,row.names = F,col.names = F)


###################################lncRNA survival analysis
#######################Breast cancer analysis
###########################MCF
MCF_ceRNA=read.csv("MCF_degree_symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
XX=which(MCF_ceRNA$V4==2)
MCF_lncRNA=MCF_ceRNA[XX,]
Result_clinical=c()
for (i in 1:dim(MCF_lncRNA)[1]) {
  print(i)
  lnc=MCF_lncRNA$V2[i]
  xx=which(allgene==lnc)
  lnc_EXP=BRCAmatrix[xx,inter_mcf$Cancer_s]
  lnc_EXP=as.numeric(as.character(lnc_EXP))
  low=which(lnc_EXP<=mean(lnc_EXP))
  high=which(lnc_EXP>mean(lnc_EXP))
  low_clinical=inter_mcf[low,]
  high_clinical=inter_mcf[high,]
  lnc_clinical=c()
  p_age=wilcox.test(high_clinical$Age.at.Initial.Pathologic.Diagnosis,low_clinical$Age.at.Initial.Pathologic.Diagnosis)$p.value
  ##########stage
  a1=which(high_clinical$AJCC.Stage=="Stage I"|high_clinical$AJCC.Stage=="Stage IA"|high_clinical$AJCC.Stage=="Stage IB")
  a2=which(high_clinical$AJCC.Stage=="Stage II"|high_clinical$AJCC.Stage=="Stage IIA"|high_clinical$AJCC.Stage=="Stage IIB")
  a3=which(high_clinical$AJCC.Stage=="Stage III"|high_clinical$AJCC.Stage=="Stage IIIA"|high_clinical$AJCC.Stage=="Stage IIIB"|high_clinical$AJCC.Stage=="Stage IIIC")
  a4=which(high_clinical$AJCC.Stage=="Stage IV"|high_clinical$AJCC.Stage=="Stage X")
  
  b1=which(low_clinical$AJCC.Stage=="Stage I"|low_clinical$AJCC.Stage=="Stage IA"|low_clinical$AJCC.Stage=="Stage IB")
  b2=which(low_clinical$AJCC.Stage=="Stage II"|low_clinical$AJCC.Stage=="Stage IIA"|low_clinical$AJCC.Stage=="Stage IIB")
  b3=which(low_clinical$AJCC.Stage=="Stage III"|low_clinical$AJCC.Stage=="Stage IIIA"|low_clinical$AJCC.Stage=="Stage IIIB"|low_clinical$AJCC.Stage=="Stage IIIC")
  b4=which(low_clinical$AJCC.Stage=="Stage IV"|low_clinical$AJCC.Stage=="Stage X")
  
  Stage=matrix(c(length(a1),length(a2),length(a3),length(a4),length(b1),length(b2),length(b3),length(b4)),nrow = 2)
  p_stage=fisher.test(Stage,simulate.p.value = TRUE,B=1.0e-5)$p.value
  KK=c()
  KK=c(1:(dim(high_clinical)[1]+dim(low_clinical)[1]))
  KK[low]=0
  KK[high]=1
  TimeSur=c()
  TimeSur=rbind(t(t(low_clinical$OS.Time)),t(t(high_clinical$OS.Time)))
  LD=c()
  LD=rbind(t(t(low_clinical$OS.event)),t(t(high_clinical$OS.event)))
  library(survival)
  fit<-survdiff(Surv(TimeSur,LD==1)~KK)
  p_surv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
  lnc_clinical=cbind(lnc,p_age,p_stage,p_surv)
  Result_clinical=rbind(Result_clinical,lnc_clinical)
}
############################MDA
MDA_ceRNA=read.csv("MDA_degree_symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
XX=which(MDA_ceRNA$V4==2)
MDA_lncRNA=MDA_ceRNA[XX,]
Result_clinical=c()
for (i in 1:dim(MDA_lncRNA)[1]) {
  print(i)
  lnc=MDA_lncRNA$V2[i]
  xx=which(allgene==lnc)
  lnc_EXP=BRCAmatrix[xx,inter_mda$Cancer_s]
  lnc_EXP=as.numeric(as.character(lnc_EXP))
  low=which(lnc_EXP<=median(lnc_EXP))
  high=which(lnc_EXP>median(lnc_EXP))
  low_clinical=inter_mda[low,]
  high_clinical=inter_mda[high,]
  lnc_clinical=c()
  p_age=wilcox.test(high_clinical$Age.at.Initial.Pathologic.Diagnosis,low_clinical$Age.at.Initial.Pathologic.Diagnosis)$p.value
  ##########stage
  a1=which(high_clinical$AJCC.Stage=="Stage I"|high_clinical$AJCC.Stage=="Stage IA"|high_clinical$AJCC.Stage=="Stage IB")
  a2=which(high_clinical$AJCC.Stage=="Stage II"|high_clinical$AJCC.Stage=="Stage IIA"|high_clinical$AJCC.Stage=="Stage IIB")
  a3=which(high_clinical$AJCC.Stage=="Stage III"|high_clinical$AJCC.Stage=="Stage IIIA"|high_clinical$AJCC.Stage=="Stage IIIB"|high_clinical$AJCC.Stage=="Stage IIIC")
  a4=which(high_clinical$AJCC.Stage=="Stage IV"|high_clinical$AJCC.Stage=="Stage X")
  
  b1=which(low_clinical$AJCC.Stage=="Stage I"|low_clinical$AJCC.Stage=="Stage IA"|low_clinical$AJCC.Stage=="Stage IB")
  b2=which(low_clinical$AJCC.Stage=="Stage II"|low_clinical$AJCC.Stage=="Stage IIA"|low_clinical$AJCC.Stage=="Stage IIB")
  b3=which(low_clinical$AJCC.Stage=="Stage III"|low_clinical$AJCC.Stage=="Stage IIIA"|low_clinical$AJCC.Stage=="Stage IIIB"|low_clinical$AJCC.Stage=="Stage IIIC")
  b4=which(low_clinical$AJCC.Stage=="Stage IV"|low_clinical$AJCC.Stage=="Stage X")
  
  Stage=matrix(c(length(a1),length(a2),length(a3),length(a4),length(b1),length(b2),length(b3),length(b4)),nrow = 2)
  p_stage=fisher.test(Stage,simulate.p.value = TRUE,B=1.0e-5)$p.value
  KK=c()
  KK=c(1:(dim(high_clinical)[1]+dim(low_clinical)[1]))
  KK[low]=0
  KK[high]=1
  TimeSur=c()
  TimeSur=rbind(t(t(low_clinical$OS.Time)),t(t(high_clinical$OS.Time)))
  LD=c()
  LD=rbind(t(t(low_clinical$OS.event)),t(t(high_clinical$OS.event)))
  library(survival)
  fit<-survdiff(Surv(TimeSur,LD==1)~KK)
  p_surv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
  lnc_clinical=cbind(lnc,p_age,p_stage,p_surv)
  Result_clinical=rbind(Result_clinical,lnc_clinical)
}

########################K562
K562_ceRNA=read.csv("K562_degree_symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
XX=which(K562_ceRNA$V4==2)
K562_lncRNA=K562_ceRNA[XX,]
Result_clinical=c()
for (i in 1:dim(K562_lncRNA)[1]) {
  print(i)
  lnc=K562_lncRNA$V2[i]
  xx=which(allgene==lnc)
  lnc_EXP=LAMLmatrix[xx,]
  lnc_EXP=as.numeric(as.character(lnc_EXP))
  low=which(lnc_EXP<=median(lnc_EXP))
  high=which(lnc_EXP>median(lnc_EXP))
  
  low_clinical=inter_mda[low,]
  high_clinical=inter_mda[high,]
  lnc_clinical=c()
  p_age=wilcox.test(high_clinical$Age.at.Initial.Pathologic.Diagnosis,low_clinical$Age.at.Initial.Pathologic.Diagnosis)$p.value
  Result_clinical=rbind(Result_clinical,lnc_clinical)
}

####################lncRNA RBP analysis
HepG2=read.csv("RBP_gene_inter_HepG2.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
K562=read.csv("RBP_gene_inter_K562.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
x1=which(HepG2$V3=="lincRNA")
Hep_linc=HepG2[x1,]
x1=which(K562$V3=="lincRNA")
K562_linc=K562[x1,]
Essential=read.csv("essential lncRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
##########hepg2 lincRNA degree
lincRNA_degree_hep=c()
lin_hep=unique(Hep_linc$V2)
for (i in 1:length(lin_hep)) {
  aa=which(Hep_linc$V2==lin_hep[i])
  bb=unique(Hep_linc$V1[aa])
  lincRNA_degree_hep=rbind(lincRNA_degree_hep,cbind(lin_hep[i],length(bb)))
}
###########k562 lincRNa degree
lincRNA_degree_k562=c()
lin_k562=unique(K562_linc$V2)
for (i in 1:length(lin_k562)) {
  aa=which(K562_linc$V2==lin_k562[i])
  bb=unique(K562_linc$V1[aa])
  lincRNA_degree_k562=rbind(lincRNA_degree_k562,cbind(lin_k562[i],length(bb)))
}

ess_hep=intersect(Essential$V1,lincRNA_degree_hep[,1])
ess_hep=t(t(ess_hep))
ess_other=setdiff(lincRNA_degree_hep[,1],Essential$V1)
ess_other=t(t(ess_other))
ess_hep_degree=merge(lincRNA_degree_hep,ess_hep,by.x="V1",by.y="V1")
ess_other_degree=merge(lincRNA_degree_hep,ess_other,by.x="V1",by.y="V1")
ess_hep_degree=cbind(ess_hep_degree,1)
ess_other_degree=cbind(ess_other_degree,0)

X1=as.numeric(as.character(ess_hep_degree$V2))
X2=as.numeric(as.character(ess_other_degree$V2))
RR_h=c()
n=length(X1)
N=length(X2)
KK=c()
for (i in 1:10000) {
  xx=sample(c(1:N),n)
  RR_h=rbind(RR_h,X2[xx])
  KK=rbind(KK,mean(X2[xx]))
}
RR_h=matrix(RR_h,ncol = 1)
wilcox.test(X1,RR_h,alternative = "greater") ###p=0.07


ess_kk=intersect(Essential$V1,lincRNA_degree_k562[,1])
ess_kk=t(t(ess_kk))
ess_other1=setdiff(lincRNA_degree_k562[,1],Essential$V1)
ess_other1=t(t(ess_other1))
ess_kk_degree=merge(lincRNA_degree_k562,ess_kk,by.x="V1",by.y="V1")
ess_other_degree1=merge(lincRNA_degree_k562,ess_other1,by.x="V1",by.y="V1")
ess_kk_degree=cbind(ess_kk_degree,1)
ess_other_degree1=cbind(ess_other_degree1,0)

X3=as.numeric(as.character(ess_kk_degree$V2))
X4=as.numeric(as.character(ess_other_degree1$V2))
RR_h_k=c()
n=length(X3)
N=length(X4)
for (i in 1:10000) {
  xx=sample(c(1:N),n)
  RR_h_k=rbind(RR_h_k,X4[xx])
}
RR_h_k=matrix(RR_h_k,ncol = 1)
wilcox.test(as.numeric(X3),as.numeric(RR_h_k),alternative = "greater")#######p=0.18


par(mfrow=c(3,3))
boxplot(X1,RR_h)
boxplot(X3,RR_h_k)

EE_inter_HEP=merge(Essential,Hep_linc,by.x = "V1",by.y = "V2")
EE_inter_KK=merge(Essential,K562_linc,by.x = "V1",by.y = "V2")
write.table(EE_inter_KK,file ="RBP_linc_K562",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
write.table(EE_inter_HEP,file ="RBP_linc_HEP",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

write.table(K562_linc,file ="K562_RBP_linc.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
write.table(Hep_linc,file ="Hep_RBP_linc.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
###############combined RBP-lincRNA
RLnet=rbind(Hep_linc[,1:2],K562_linc[,1:2])
RLnet=RLnet[!duplicated(RLnet),]
lincRNA_degree_combine=c()
lin_combine=unique(RLnet$V2)
for (i in 1:length(lin_combine)) {
  aa=which(RLnet$V2==lin_combine[i])
  bb=unique(RLnet$V1[aa])
  lincRNA_degree_combine=rbind(lincRNA_degree_combine,cbind(lin_combine[i],length(bb)))
}

ess_kk=intersect(Essential$V1,lincRNA_degree_combine[,1])
ess_kk=t(t(ess_kk))
ess_other1=setdiff(lincRNA_degree_combine[,1],Essential$V1)
ess_other1=t(t(ess_other1))
ess_kk_degree=merge(lincRNA_degree_combine,ess_kk,by.x="V1",by.y="V1")
ess_other_degree1=merge(lincRNA_degree_combine,ess_other1,by.x="V1",by.y="V1")

X3=as.numeric(as.character(ess_kk_degree$V2))
X4=as.numeric(as.character(ess_other_degree1$V2))
RR_h_k=c()
n=length(X3)
N=length(X4)
for (i in 1:10000) {
  xx=sample(c(1:N),n)
  RR_h_k=rbind(RR_h_k,X4[xx])
}
RR_h_k=matrix(RR_h_k,ncol = 1)
wilcox.test(as.numeric(X3),as.numeric(RR_h_k),alternative = "greater")#######p=0.18

####################rbp degree analysis
HepG2=read.csv("RBP_gene_inter_HepG2.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
K562=read.csv("RBP_gene_inter_K562.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
x1=which(HepG2$V3=="lincRNA")
Hep_linc=HepG2[x1,]
x1=which(K562$V3=="lincRNA")
K562_linc=K562[x1,]
Essential=read.csv("essential lncRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)

RBP_hep=unique(Hep_linc$V1)
RBP_hep_Degree=c()
for (i in 1:length(RBP_hep)) {
  x1=which(Hep_linc$V1==RBP_hep[i])
  lc=unique(Hep_linc$V2[x1])
  aa=length(intersect(lc,Essential$V1))
  RBP_hep_Degree=rbind(RBP_hep_Degree,cbind(RBP_hep[i],length(lc),aa))
}

RBP_kk=unique(K562_linc$V1)
RBP_kk_Degree=c()
for (i in 1:length(RBP_kk)) {
  x1=which(K562_linc$V1==RBP_kk[i])
  lc=unique(K562_linc$V2[x1])
  aa=length(intersect(lc,Essential$V1))
  RBP_kk_Degree=rbind(RBP_kk_Degree,cbind(RBP_kk[i],length(lc),aa))
}

setwd("C:/projects/PPI/mapnet")
APMS=read.csv("APMS.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
HI3=read.csv("HI3_Y2H_102416.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
DD=c()
for (i in 1:dim(RBP_hep_Degree)[1]) {
  x1=which(APMS$Symbol.A==RBP_hep_Degree[i,1])
  g1=unique(APMS$Symbol.B[x1])
  x2=which(APMS$Symbol.B==RBP_hep_Degree[i,1])
  g2=unique(APMS$Symbol.A[x2])
  g=union(g1,g2)
  y1=which(HI3$Symbol.for.A==RBP_hep_Degree[i,1])
  gg1=unique(HI3$Symbol.for.B[y1])
  y2=which(HI3$Symbol.for.B==RBP_hep_Degree[i,1])
  gg2=unique(HI3$Symbol.for.A[y1])
  gg=union(gg1,gg2)
  DD=rbind(DD,cbind(length(g),length(gg)))
}
RBP_hep_f=cbind(RBP_hep_Degree,DD)

DD1=c()
for (i in 1:dim(RBP_kk_Degree)[1]) {
  x1=which(APMS$Symbol.A==RBP_kk_Degree[i,1])
  g1=unique(APMS$Symbol.B[x1])
  x2=which(APMS$Symbol.B==RBP_kk_Degree[i,1])
  g2=unique(APMS$Symbol.A[x2])
  g=union(g1,g2)
  y1=which(HI3$Symbol.for.A==RBP_kk_Degree[i,1])
  gg1=unique(HI3$Symbol.for.B[y1])
  y2=which(HI3$Symbol.for.B==RBP_kk_Degree[i,1])
  gg2=unique(HI3$Symbol.for.A[y1])
  gg=union(gg1,gg2)
  DD1=rbind(DD1,cbind(length(g),length(gg)))
}
RBP_kk_f=cbind(RBP_kk_Degree,DD1)
##################drug analysis
setwd("C:/projects/ceRNA")
cancer=c("MCF","MDA","K562","Hela","U87")
ensembe2symbol=read.csv("Ens2sym.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essential=read.csv("essential lncRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
setwd("C:/projects/ceRNA/drug")
drug=read.csv("drug_target.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
drug_t=merge(drug,ensembe2symbol,by.x = "V2",by.y = "V2")
colnames(drug_t)=c("t_symbol","drug","t_ensg")
drug=unique(drug_t$V1.x)
N=length(unique(ensembe2symbol$V2))
for (i in 1:5) {
  setwd("C:/projects/ceRNA")
  CeRNAnet=read.csv(paste(cancer[i],"_ceRNA_sorted_final.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = F)
  x=which(essential$V2==cancer[i])
  ess_linc=essential$V1[x]
  ess_linc=t(t(ess_linc))
  colnames(ess_linc)=c("lincRNA")
  colnames(CeRNAnet)=c("ceRNA1","ceRNA2")
  xx1=merge(CeRNAnet,ess_linc,by.x = "ceRNA1",by.y = "lincRNA")
  Dru1=merge(xx1,drug_t,by.x = "ceRNA2",by.y = "t_ensg")

  xx2=merge(CeRNAnet,ess_linc,by.x = "ceRNA2",by.y = "lincRNA")
  Dru2=merge(xx2,drug_t,by.x = "ceRNA1",by.y = "t_ensg")
  
  Dru=rbind(Dru1,Dru2)
  # drug_lncRNA=c()
  # for (k1 in 1:length(drug)) {
  #   for (k2 in 1:length(ess_linc)) {
  #     xx1=which(drug_t$V1.x==drug[k1])
  #     tt1=unique(drug_t$V1.y[xx1])
  #     
  #     xx2=which(CeRNAnet$V1==ess_linc[k2])
  #     t21=unique(CeRNAnet$V2[xx2])
  #     xx3=which(CeRNAnet$V2==ess_linc[k2])
  #     t22=unique(CeRNAnet$V1[xx3])
  #     tt2=union(t21,t22)
  #     
  #     aaa=intersect(tt1,tt2)
  #     print(length(tt2))
  #     if(length(aaa)>2){
  #       n=length(tt1)
  #       m=length(tt2)
  #       o=length(aaa)
  #       
  #       p.value <- phyper(o-1,m,N-m,n,lower.tail=FALSE, log.p = FALSE)
  #       drug_lncRNA=rbind(drug_lncRNA,cbind(drug[k1],ess_linc[k2],p.value))
  #     }
  #   }
  # }
}

################SURVIVAL ANALYSIS
setwd("C:/projects/ceRNA")
load("C:/projects/ceRNA/K562finalceRNA.RData")
xx=which(allgene=="ENSG00000179818")
lncS=LAMLmatrix[xx,]
Labb=rep(0,151)
highs=which(lncS>median(lncS))
lows=which(lncS<=median(lncS))
Labb[highs]=1
Sampleall=colnames(LAMLmatrix)
Survitime1=LAMLdata@colData@listData$days_to_last_follow_up
Survitime2=LAMLdata@colData@listData$days_to_death
Surtime=c()
for (i in 1:length(Survitime2)) {
  # if((Survitime1[i]=="NA")=TRUE){
  #   Surtime=rbind(Surtime,Survitime2[i])
  # } else if ((Survitime2[i]=="NA")=TRUE){
  #   Surtime=rbind(Surtime,Survitime1[i])
  # } else {
     Surtime=rbind(Surtime,max(Survitime1[i],Survitime2[i],na.rm = TRUE))
  # }
}
XX=which((Surtime=="-Inf")==F)
Virtu=LAMLdata@colData@listData$vital_status
library(survival)
fit<-survdiff(Surv(Surtime[XX],Virtu[XX]=="dead")~Labb[XX])
fit

###########################funseq2
setwd("C:/projects/ceRNA/mutation")
mutscore=read.csv("TCGA_funseq2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
XX=which(mutscore$V9!=0)
mutscore=mutscore[XX,]
Esslinc=read.csv("Esse_mut.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Esscore=merge(Esslinc,mutscore,by.x = c("V1","V2","V3","V6","V7"),
              by.y = c("V1","V2","V5","V3","V4"))
Randscore=c()
n=dim(Esscore)[1]
N=dim(mutscore)[1]
for (i in 1:100) {
  XX=sample(c(1:N),n)
  Randscore=rbind(Randscore,mutscore$V9[XX])
}
Randscore=matrix(Randscore,ncol = 1)
realscore=as.numeric(as.character(Esscore$V9))
Randscore=as.numeric(as.character(Randscore))
library(vioplot)
par(mfrow=c(3,3))
vioplot(Randscore,realscore)
wilcox.test(realscore,Randscore,alternative = "greater")

#################################################### miRNA regulation comparsion
miRNA_lnc_ensg=read.csv("miRNA_lncRNA_ensg.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lncRNA=unique(miRNA_lnc_ensg$V2)
lnc_degree=c()
for (i in 1:length(lncRNA)) {
  xx=which(miRNA_lnc_ensg$V2==lncRNA[i])
  lnc_degree=rbind(lnc_degree,length(xx))
}
lnc_degree=cbind(lncRNA,lnc_degree)
Esslinc=read.csv("essential lncRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lnc_label=c()
for (i in 1:length(lncRNA)) {
  xx=which(Esslinc$V1==lncRNA[i])
  if(length(xx)>0){
    lnc_label=rbind(lnc_label,1)
  } else {
    lnc_label=rbind(lnc_label,0)
  }
}
par(mfrow=c(3,3))
boxplot(as.numeric(lnc_degree[,2])~lnc_label,outline=FALSE,boxwex=c(0.5,0.5))
aa=which(lnc_label==0)
bb=which(lnc_label==1)
wilcox.test(as.numeric(lnc_degree[bb,2],as.numeric(lnc_degree[aa,2])),alternative = "greater")
cancer=c("MCF","MDA","K562","Hela","U87")
lnc_degree_label=cbind(lnc_degree,lnc_label)
for (i in 1:5) {
  xx=which(Esslinc$V2==cancer[i])
  cancerlnc=Esslinc[xx,]
  cancerlab=c()
  for (j in 1:length(lncRNA)) {
    a=which(cancerlnc$V1==lncRNA[j])
    if(length(a)>0){
      cancerlab=rbind(cancerlab,1)
    } else {
      cancerlab=rbind(cancerlab,0)
    }
  }
  boxplot(as.numeric(lnc_degree[,2])~cancerlab,outline=FALSE,boxwex=c(0.5,0.5))
  aa=which(cancerlab==0)
  bb=which(cancerlab==1)
  pp=wilcox.test(as.numeric(lnc_degree[bb,2],as.numeric(lnc_degree[aa,2])),alternative = "greater")$p.value
  print(pp)
}
########################CGC and neighbors miRNA regulation
miRNA_gene=read.csv("miRNA_gene3.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene=read.csv("Finalgene.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
U_gene=unique(miRNA_gene$V2)
U_miRNAs=c()
for (i in 1:length(U_gene)) {
  xx=which(miRNA_gene$V2==U_gene[i])
  Umir=unique(miRNA_gene$V1[xx])
  U_miRNAs=rbind(U_miRNAs,cbind(U_gene[i],length(Umir)))
}
gene_mir=merge(gene,U_miRNAs,by.x="V1",by.y="V1")
othergene=setdiff(U_gene,gene$V1)
othergene=t(t(othergene))
other_mir=merge(othergene,U_miRNAs,by.x="V1",by.y="V1")
par(mfrow=c(3,3))
AA=c()
for (i in 1:100) {
  xx=sample(c(1:5037),696)
  AA=rbind(AA,other_mir[xx,])
}
boxplot(as.numeric(AA$V2),as.numeric(gene_mir$V2),outline=FALSE,boxwex=c(0.5,0.5))
wilcox.test(as.numeric(gene_mir$V2),as.numeric(AA$V2),alternative = "greater")

################################subtype ananlysis
U87_hub=read.csv("U87_elnc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
hubid=c()
for (i in 1:dim(U87_hub)[1]) {
  xx=which(allgene==U87_hub$V1[i])
  hubid=rbind(hubid,xx)
}
can=which(GBMclinal!="Solid Tissue Normal")
hubexp=GBMmatrix[hubid,]
aa=rep(0,length(GBMclinal))
aa[can]=1
Colcol=rep("green",length(GBMclinal))
Colcol[can]="red"
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
heatmap(log2(hubexp+0.05),labCol = aa,col=my_palette,ColSideColors=Colcol)
hubcancer=GBMmatrix[hubid,can]
Canname=colnames(GBMmatrix)[can]
Threesampe=c()
for (i in 1:length(Canname)) {
  aa=strsplit(Canname[i],"-")
  bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
  Threesampe=rbind(Threesampe,cbind(i,bb))
}
hubcancer=log2(hubcancer+0.05)
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(hubcancer,maxK=10,reps=100,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf",title="GBM")



Hela_hub=read.csv("Hela_hub.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
hubid=c()
for (i in 1:dim(Hela_hub)[1]) {
  xx=which(allgene==Hela_hub$V1[i])
  hubid=rbind(hubid,xx)
}
can=which(UCECclinal!="Solid Tissue Normal")
hubexp=UCECmatrix[hubid,]
aa=rep(0,length(UCECclinal))
aa[can]=1
Colcol=rep("green",length(UCECclinal))
Colcol[can]="red"
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
heatmap(log2(hubexp+0.05),labCol = aa,col=my_palette,ColSideColors=Colcol)

K562_hub=read.csv("K562_hub.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
hubid=c()
for (i in 1:dim(K562_hub)[1]) {
  xx=which(allgene==K562_hub$V1[i])
  hubid=rbind(hubid,xx)
}
can=which(LAMLclinal!="Solid Tissue Normal")
hubexp=LAMLmatrix[hubid,]
aa=rep(0,length(LAMLclinal))
aa[can]=1
Colcol=rep("green",length(LAMLclinal))
Colcol[can]="red"
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
heatmap(log2(hubexp+0.05),labCol = aa,col=my_palette,ColSideColors=Colcol)
hubcancer=LAMLmatrix[hubid,can]
Canname=colnames(LAMLmatrix)[can]
Threesampe=c()
for (i in 1:length(Canname)) {
  aa=strsplit(Canname[i],"-")
  bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
  Threesampe=rbind(Threesampe,cbind(i,bb))
}
hubcancer=log10(hubcancer+0.05)
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(hubcancer,maxK=10,reps=100,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="spearman",seed=1262118388.71279,plot="pdf",title="LAML")
##############################essential lncRNA and immune score
setwd("C:/projects/ceRNA")
essen_lnc=read.csv("all_ensial.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essen_lnc=essen_lnc[!duplicated(essen_lnc),]
immune_score=read.csv("immunescore.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
#cancer=unique(immune_score$TumorType)
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
cancer<-c("BRCA","UCEC","LAML","GBM")
setwd("C:/projects/AS/geneexpression")
genecode=read.csv("gencode_gene_pro.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)

Allexp=c()
Allgeneexp=c()
for (i in 1:length(cancer)) {
  print(i)
  load(paste(cancer[i],"Expression.rda",sep = ""))
  library(SummarizedExperiment)
  Expmatrix<-assay(data,1,"FPKM")
  C_N=data@colData@listData$definition#############sample information tumor vs normal
  Cancer_s=which(C_N!="Solid Tissue Normal")
  Expmatrix=Expmatrix[,Cancer_s]
  allgene=rownames(Expmatrix)
  SSam=colnames(Expmatrix)
  Canlnc=c()
  for (j in 1:dim(essen_lnc)[1]) {
    xx=which(allgene==essen_lnc$V1[j])
    if(length(xx)>0){
      Canlnc=rbind(Canlnc,xx)
    }
  }
  
  Canexp=Expmatrix[Canlnc,]
  Allexp=cbind(Allexp,Canexp)
  genename=c()
  Cangene=c()
  for (j in 1:dim(genecode)[1]) {
    xx=which(allgene==genecode$V1[j])
    if(length(xx)>0){
      Cangene=rbind(Cangene,xx)
      genename=rbind(genename,genecode$V2[j])
    }
  }
  Canexp_g=Expmatrix[Cangene,]
  Allgeneexp=cbind(Allgeneexp,Canexp_g)
}
Allexp_1=log2(Allexp+0.05)
WW=kmeans(t(Allexp_1),2)
cclabel=WW$cluster
Low=which(cclabel==1)
High=which(cclabel==2)
N=dim(Allgeneexp)[1]
generank=c()
Allgeneexp_1=log2(Allgeneexp+0.05)
for (k in 1:N) {
  print(k)
  M1=mean(Allgeneexp_1[k,Low])
  M2=mean(Allgeneexp_1[k,High])
  FC=M1/M2
  p=wilcox.test(Allgeneexp_1[k,Low],Allgeneexp_1[k,High])$p.value
  Sco=-log10(p)*sign(log2(FC))
  generank=rbind(generank,Sco)
}
KWW=cbind(genename,generank)
setwd("C:/projects/ceRNA")
write.table(KWW,file = "gene_rank314.rnk",sep = "\t",row.names = F, col.names = F,quote = FALSE)

plot(t(Allexp_1),col=WW$cluster)
LncClu=Allexp_1[,cbind(Low,High)]

par(mfrow=c(3,3))
AA=apply(LncClu, 2, mean)
boxplot(AA[1:1652],AA[1653:3304],boxwex=c(0.5,0.5))
LncClu_Norm=c()
for (i in 1:dim(LncClu)[1]) {
  xx=LncClu[i,]-min(LncClu[i,])/(max(LncClu[i,])-min(LncClu[i,]))
  LncClu_Norm=rbind(LncClu_Norm,xx)
}
aa=which(LncClu_Norm[,1]=="Inf")
LncClu_Norm1=LncClu_Norm[-aa,]
library(pheatmap)
pheatmap(LncClu_Norm1,cluster_cols = F)

Allexp=t(Allexp)
KK=apply(Allexp, 2, sum)
x=which(KK==0)
Allexp=Allexp[,-x]
Ss=rownames(Allexp)
Threesampe=c()
for (i in 1:length(Ss)) {
  aa=strsplit(Ss[i],"-")
  bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
  Threesampe=rbind(Threesampe,bb)
}
Allexp=cbind(Allexp,Threesampe)
write.table(Allexp,file = "LncRNA_pan_cancer.txt",sep = "\t",row.names = T, col.names = T,quote = FALSE)
Allexp=read.csv("LncRNA_pan_cancer.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)

AA=merge(Allexp,immune_score,by.x="X",by.y="TumorSample")
ccc=unique(AA$TumorType)
R_mat=c()
P_mat=c()
Lnc=colnames(AA)[3:118]
for (i in 1:length(ccc)) {
  for (j in 3:118) {
    xxx=which(AA$TumorType==ccc[i])
    lncRNAC=AA[xxx,j]
    lncRNAC=log2(as.numeric(as.character(lncRNAC))+0.05)
    Imm=AA$Immune.Signature.Score[xxx]
    Imm=as.numeric(as.character(Imm))
    R=cor(lncRNAC,Imm,method = "spearman",use = "complete")
    P=cor.test(lncRNAC,Imm,method = "spearman")$p.value
    R_mat=rbind(R_mat,R)
    P_mat=rbind(P_mat,P)
  }
}
R_mat1=matrix(R_mat,nrow = 12)
P_mat1=matrix(P_mat,nrow = 12)
library(pheatmap)
P_mat2=P_mat1
xx=which(P_mat2>=0.05)
P_mat2[xx]=1
xx=which(P_mat2<0.05)
P_mat2[xx]=0

pheatmap(R_mat1,cluster_rows = F,cluster_cols = F)
pheatmap(P_mat2,cluster_rows = F,cluster_cols = F)
KK=apply(P_mat2, 2, sum)
KK1=12-KK
hist(KK1,12)
xx=which(KK1==10)
###############################DNA repair gene
immune=read.csv("DNArepari.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
ensg2sym=read.csv("ensmble2symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
immune_eng=merge(immune,ensg2sym,by.x = "V1",by.y = "HGNC.symbol")
allens=read.csv("all_ensial.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
cancer=c("MCF","MDA","K562","Hela","U87")
Allimu=c()
allens=allens[!duplicated(allens),]
for (i in 1:length(cancer)) {
  immnet=c()
  cenet=read.csv(paste(cancer[i],"_ceRNA_sorted_final.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = F)
  for (j in 1:dim(cenet)[1]) {
    x1=which(immune_eng$Ensembl.Gene.ID==cenet$V1[j])
    x2=which(immune_eng$Ensembl.Gene.ID==cenet$V2[j])
    x3=which(allens$V1==cenet$V1[j])
    x4=which(allens$V1==cenet$V2[j])
    x=length(x1)+length(x2)+length(x3)+length(x4)
    if(x==2){
      immnet=rbind(immnet,cenet[j,])
    }
    
  }
  print(dim(immnet)[1])
  Allimu=rbind(Allimu,immnet)
}
aa=merge(Allimu,ensg2sym,by.x = "V1",by.y = "Ensembl.Gene.ID")
bb=merge(aa,ensg2sym,by.x = "V2",by.y = "Ensembl.Gene.ID")
write.table(bb,file = "repair_ceRNA.txt",sep = "\t",row.names = T, col.names = T,quote = FALSE)


##################new figures
install.packages("UpSetR")
library(UpSetR)
MCF_lncRNA=read.csv("MCF7_lnc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
MDA_lncRNA=read.csv("MDA_lnc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
K562_lncRNA=read.csv("K562_lnc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Hela_lncRNA=read.csv("Hela_lnc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
U87_lncRNA=read.csv("U87_lnc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Ulnc=rbind(MCF_lncRNA,MDA_lncRNA,K562_lncRNA,Hela_lncRNA,U87_lncRNA)
Ulnc=unique(Ulnc)
Celline=c("MCF7","MDA","K562","Hela","U87")
lncMAT=c()
for (i in 1:5) {
  cc=c()
  lnc=read.csv(paste(Celline[i],"_lnc.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = F)
  for (j in 1:dim(Ulnc)[1]) {
    xx=which(lnc==Ulnc$V1[j])
    if(length(xx)>0){
      cc=rbind(cc,1)
    } else{
      cc=rbind(cc,0)
    }
  }
  lncMAT=cbind(lncMAT,cc)
}
LL=cbind(Ulnc$V1,lncMAT)
write.table(LL,file = "overlap_ceRNA.txt",sep = "\t",row.names = F, col.names = F,quote = FALSE)
LL=read.csv("overlap_ceRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
upset(LL, sets = c("V6", "V5", "V4", "V3", "V2"), keep.order=T,sets.bar.color = "#56B4E9", empty.intersections = "on")


MCF_ceRNA=read.csv("MCF_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
MDA_ceRNA=read.csv("MDA_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
K562_ceRNA=read.csv("K562_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Hela_ceRNA=read.csv("Hela_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
U87_ceRNA=read.csv("U87_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Ulnc=rbind(MCF_ceRNA,MDA_ceRNA,K562_ceRNA,Hela_ceRNA,U87_ceRNA)
Ulnc=Ulnc[!duplicated(Ulnc),]
Celline=c("MCF","MDA","K562","Hela","U87")
lncMAT=c()
for (i in 1:5) {
  cc=c()
  lnc=read.csv(paste(Celline[i],"_ceRNA_sorted_final.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = F)
  for (j in 1:dim(Ulnc)[1]) {
    xx=which(lnc$V1==Ulnc$V1[j]&lnc$V2==Ulnc$V2[j])
    if(length(xx)>0){
      cc=rbind(cc,1)
    } else{
      cc=rbind(cc,0)
    }
  }
  lncMAT=cbind(lncMAT,cc)
}
LL=cbind(Ulnc,lncMAT)
write.table(LL,file = "overlap_ceRNAnet_final.txt",sep = "\t",row.names = F, col.names = F,quote = FALSE)
LL=read.csv("overlap_ceRNAnet_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
upset(LL, sets = c("V7", "V6", "V5", "V4", "V3"), keep.order=T,sets.bar.color = "#56B4E9", empty.intersections = "on")

Celline=c("MCF","MDA","K562","Hela","U87")
ZZ=YY
for (i in 4:5) {
  aa=read.csv(paste(Celline[i],"_degree_symbol.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = F)
  xx=sort.int(as.numeric(aa$V3),index.return = T,decreasing = T)
  YY=cbind(aa$V5[xx$ix],xx$x/max(xx$x))
  ZZ=merge(ZZ,YY,by.x="V1",by.y="V1")
}

write.table(ZZ,file = "ceRNA_degree_rank.txt",sep = "\t",row.names = F, col.names = F,quote = FALSE)
######################################################crispr-cas9 enrichment analysis
#################protein coding gene 33 cell lines
setwd("C:/projects/ceRNA/new")
essentgene=read.csv("core_fitness_all.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
Cell_line=colnames(essentgene)[2:6]
for (i in 1:5) {
  xx=cbind(essentgene$Gene,essentgene[,(i+1)])
  write.table(xx,file = paste(Cell_line[i],"_gene_rank.rnk",sep = ""),sep = "\t",row.names = F, col.names = F,quote = FALSE)
}
############average rank
Average_rank=c()
for (i in 1:5) {
  xx=cbind(essentgene$Gene,essentgene[,(i+1)])
  aa=rank(as.numeric(xx[,2]))
  aa=aa/max(aa)
  Average_rank=cbind(Average_rank,aa)
}
Average_rank_1=apply(Average_rank, 1, mean)
Average_rank_1=cbind(essentgene$Gene,Average_rank_1)
write.table(Average_rank_1,file = "average_gene_rank_cell.rnk",sep = "\t",row.names = F, col.names = F,quote = FALSE)
##############corfitness enrichment
corefit=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Biofunction <- file("c2.cp.v5.2.symbols.gmt", "r")
i<-0
N <-19022
line=readLines(Biofunction,n=1)
Result1 <- c()
while( length(line) != 0 ) {
  i <- i+1
  msig.go <- strsplit(line,"\t")
  msig.go <- msig.go[[1]]
  m <- length(msig.go)
  goterm <- msig.go[1]
  msig.go <- msig.go[3:m]
  n1 <- length(corefit$V1)
  intergene1 <- intersect(msig.go,corefit$V1)
  k1 <- length(intergene1)
  print(k1)
  if (k1>2){
    p.value1 <- phyper(k1-1,m,N-m,n1,lower.tail=FALSE, log.p = FALSE)
    go <- c(goterm,k1,m,n1,p.value1)
    Result1 <- rbind(Result1,go)
  } else {
    go <- c(goterm,k1,m,n1,1)
    Result1 <- rbind(Result1,go)
  }
  line=readLines(Biofunction,n=1);
}
close(Biofunction)
FD=p.adjust(as.numeric(Result1[,5]))
Result1=cbind(Result1,FD)
write.table(Result1,"All_core_pathway.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
###################gene and lncRNA
setwd("C:/projects/ceRNA/new/KEGG")
gene=read.csv("gene.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lncRNA=read.csv("lncRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene_OP=read.csv("gene_OP.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lnc_OP=read.csv("lnc_OP.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Over_gene=intersect(gene_OP$V2,lnc_OP$V2)
Over_gene=cbind(Over_gene,"red")
L_gene=setdiff(lnc_OP$V2,gene_OP$V2)
L_gene=cbind(L_gene,"yellow")
G_gene=setdiff(gene_OP$V2,lnc_OP$V2)
G_gene=cbind(G_gene,"blue")
ZZ=rbind(Over_gene,L_gene,G_gene)
write.table(ZZ,"OP_KEGG.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
aa=intersect(gene$V1,lncRNA$V1)
##################################################identify the ceRNAs in gene expression
MCF=read.csv("MCF_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
MDA=read.csv("MDA_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
K562=read.csv("K562_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Hela=read.csv("Hela_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
U87=read.csv("U87_ceRNA_sorted_final.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
ens2sym=read.csv("Ens2sym.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
XX=rbind(MCF,MDA)
XX=rbind(XX,K562)
XX=rbind(XX,Hela)
XX=rbind(XX,U87)
XX=XX[!duplicated(XX),]
Un_net=c()
for (i in 1:dim(XX)[1]) {
  print(i)
  x1=which(ens2sym$V1==XX$V1[i])
  x2=which(ens2sym$V1==XX$V2[i])
  if(length(x1)>0&length(x2)>0){
    Un_net=rbind(Un_net,cbind(ens2sym$V2[x1],ens2sym$V2[x2]))
  }
}
x1=unique(Un_net[,1])
x2=unique(Un_net[,2])
gene=union(x1,x2)
rgene=read.csv("gene_rank314.rnk",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene=t(t(gene))
new_rank=merge(rgene,gene,by.x = "V1",by.y = "V1")
write.table(new_rank,"ceRNA_gene_rank.rnk",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
Finalgene=read.csv("Finalgene.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene_r=merge(rgene,Finalgene,by.x = "V1",by.y = "V1")
write.table(gene_r,"primary_gene_rank.rnk",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
#######################essential gene analysis
################hub analysis
setwd("C:/projects/ceRNA/new")
corefit=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header =F )
#PPI=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
HI3=read.csv("SIN.csv",header = TRUE,sep = ",")
othergene=setdiff(HI3$name,corefit)
othergene=t(t(othergene))
core_topo=merge(HI3,corefit,by.x = "name",by.y = "V1")
other_topo=merge(HI3,othergene,by.x = "name",by.y = "V1")
X1=as.numeric(as.character(core_topo$Degree))
X2=as.numeric(as.character(other_topo$Degree))
X3=as.numeric(as.character(core_topo$BetweennessCentrality))
X4=as.numeric(as.character(other_topo$BetweennessCentrality))
par(mfrow=c(3,3))
boxplot(X1,X2,boxwex=c(0.5,0.5),ylab="Degree",outline = FALSE,col=c("MediumSeaGreen","Gray"))
wilcox.test(X1,X2,alternative = "greater")##########
boxplot(X3,X4,boxwex=c(0.5,0.5),ylab="Betweenness",outline = FALSE,col=c("MediumSeaGreen","Gray"))
wilcox.test(X3,X4,alternative = "greater")##########
####################mutation freq
corefit=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header =F )
hg38=read.csv("hg38freq_all.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
ORF=read.csv("human_ORF.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
ORF_f=merge(ORF,hg38,by.x = "ORF_ID",by.y = "ORF_ID")
cor_freq=merge(ORF_f,corefit,by.x = "Gene_symbol",by.y = "V1")
othergene=setdiff(ORF_f$Gene_symbol,corefit$V1)
othergene=t(t(othergene))
other_freq=merge(ORF_f,othergene,by.x = "Gene_symbol",by.y = "V1")
par(mfrow=c(3,3))
X1=as.numeric(as.character(cor_freq$CADD_raw))
X2=as.numeric(as.character(other_freq$CADD_raw))
boxplot(X1,X2,boxwex=c(0.5,0.5),ylab="Conservation",outline = FALSE,col=c("MediumSeaGreen","Gray"))
wilcox.test(X1,X2,alternative = "greater")##########
###########expression 
corefit=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header =F )
setwd("C:/projects/genefusion/data/HPA")
Rnaexp=read.csv("rna_tissue.csv",stringsAsFactors=F,sep=",",skip=0,header = T)
Gene=unique(Rnaexp$Gene.name)
Avgexp=c()
for (i in 1:length(Gene)) {
  print(i)
  xx=which(Rnaexp$Gene.name==Gene[i])
  GGexp=Rnaexp$Value[xx]
  GGexp=as.numeric(as.character(GGexp))
  Avgexp=rbind(Avgexp,cbind(Gene[i],mean(GGexp)))
}
core_exp=merge(corefit,Avgexp,by.x = "V1",by.y = "V1")
othergene=setdiff(Avgexp[,1],corefit)
othergene=t(t(othergene))
otherexp=merge(othergene,Avgexp,by.x="V1",by.y="V1")
par(mfrow=c(3,3))
X1=as.numeric(as.character(core_exp[,2]))
X2=as.numeric(as.character(otherexp[,2]))
boxplot(X1,X2,boxwex=c(0.5,0.5),ylab="Expression",outline = FALSE,col=c("MediumSeaGreen","Gray"))
wilcox.test(X1,X2,alternative = "greater")##########
setwd("C:/projects/ceRNA/new")
###############science2
SC2=read.csv("aac7041_SM_Table_S3.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
xx1=which(SC2$KBM7.CS<0&SC2$KBM7.adjusted.p.value<0.05)
gene1=unique(SC2$Gene[xx1])
xx2=which(SC2$K562.CS<0&SC2$K562.adjusted.p.value<0.05)
gene2=unique(SC2$Gene[xx2])
xx3=which(SC2$Jiyoye.CS<0&SC2$Jiyoye.adjusted.p.value<0.05)
gene3=unique(SC2$Gene[xx3])
xx4=which(SC2$Raji.CS<0&SC2$Raji.adjusted.p.value<0.05)
gene4=unique(SC2$Gene[xx4])
Gene=union(gene1,gene2)
Gene=union(Gene,gene3)
Gene=union(Gene,gene4)
Num=c()
for (i in 1:length(Gene)) {
  x1=which(gene1==Gene[i])
  x2=which(gene2==Gene[i])
  x3=which(gene3==Gene[i])
  x4=which(gene4==Gene[i])
  Num=rbind(Num,length(x1)+length(x2)+length(x3)+length(x4))
}
xx=which(Num>=2)
Core2=Gene[xx]
write.table(Core2,"core_2.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
###########################essential gene analysis
####
setwd("C:/projects/ceRNA/new")
cell_ess=read.csv("core_fitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
Core_num=c()
for (i in 1:5) {
  print(i)
  xx=which(cell_ess$numTKOHits==i)
  n=length(xx)
  Core_num=rbind(Core_num,n)
}
par(mfrow=c(3,3))
barplot(t(Core_num))
SC1=c()
for (i in 1:max(Num)) {
  print(i)
  aa=which(Num==i)
  SC1=rbind(SC1,length(aa))
}
barplot(t(SC1))

###############essential genes in tissue specific PPI
setwd("C:/projects/ceRNA/new")
fusion_gene=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
###pan-cancer analysis
Ens2sym=read.csv("Ens2sym.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
fusion_gene=merge(fusion_gene,Ens2sym,by.x="V1",by.y="V2")
colnames(fusion_gene)=c("gene","ensg")
setwd("C:/projects/genefusion/data/TIN/HPAR")
tissue=list.files(pattern=".tsv _degree",recursive=T)
fus_degree=c()
other_degree=c()
All_d=c()
for (i in 1:length(tissue)) {
  print(i)
  degree=read.csv(tissue[i],stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
  fus_d=merge(degree,fusion_gene,by.x="V1",by.y="ensg")
  others=setdiff(degree$V1,fusion_gene$ensg)
  others=t(t(others))
  other_d=merge(degree,others,by.x="V1",by.y="V1")
  n1=mean(as.numeric(as.character(fus_d$V2)))
  n2=mean(as.numeric(as.character(other_d$V2)))
  fus_degree=rbind(fus_degree,n1)
  other_degree=rbind(other_degree,n2)
  aa=cbind(fus_d$V2,"Fus",tissue[i])
  bb=cbind(other_d$V2,"Other",tissue[i])
  cc=rbind(aa,bb)
  All_d=rbind(All_d,cc)
}
write.table(All_d,"All_degree.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
All_d=read.csv("All_degree.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
colnames(All_d)=c("Leg","type","tiss")
cdata <- ddply(All_d, c("type", "tiss"), summarise,
               N    = length(Leg),
               mean = mean(Leg),
               sd   = sd(Leg),
               se   = sd / sqrt(N)
)

setwd("C:/projects/ceRNA/new")
library(ggplot2)
pd <- position_dodge(0.1)
ggplot(cdata, aes(x=tiss, y=mean, colour=type)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)

KK=cbind(tissue,fus_degree,other_degree)
aa=sort.int(as.numeric(as.character(KK[,2])),index.return=T)
KK_sort=KK[aa$ix,]
x=c(1:29)
par(mfrow=c(3,1))
plot(x,as.numeric(KK_sort[,2])-10,lty=1,lwd=1,col="Red",pch=16,ylim=c(20,50))
grid (29,15, lty = 6, col = "cornsilk2")
axis(1, at=seq(1,29,by=1),labels=KK_sort[,1], las = 2)
lines(x,as.numeric(KK_sort[,2])-10,lty=1,lwd=1)
points(x,as.numeric(KK_sort[,3]),lty=1,lwd=1,col="gray",pch=16)
lines(x,as.numeric(KK_sort[,3]),lty=1,lwd=1)
points(x,as.numeric(KK_sort[,2])-10,lty=1,lwd=1,col="Red",pch=16,ylim=c(20,50))
##################################
#########KEGG and gene ontology analysis
setwd("C:/projects/ceRNA/new/PathDB")
lncBP=read.csv("lncRNA_BP.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
cellBP=read.csv("cell_BP.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
S1BP=read.csv("S1_BP.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
S2BP=read.csv("S2_BP.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
x1=intersect(lncBP$term_name,cellBP$term_name)
x2=intersect(x1,S1BP$term_name)
xterm=intersect(x2,S2BP$term_name)


lncPA=read.csv("lncRNA_KEGG.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
cellPA=read.csv("cell_KEGG.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
S1PA=read.csv("S1_KEGG.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
S2PA=read.csv("S2_KEGG.tab",stringsAsFactors=F,sep="\t",skip=0,header = T)
x1=intersect(lncPA$pathway,cellPA$pathway)
x2=intersect(x1,S1PA$pathway)
xpathway=intersect(x2,S2PA$pathway)
P_pathway=c()
KEGG <- function(pathway, enrich ){
  P=c()
  for (i in 1:length(pathway)) {
    xx1=which(enrich$pathway==pathway[i])
    P=rbind(P,enrich$q.value[xx1])
  }
  return(P)
}
P1=KEGG(xpathway,lncPA)
P2=KEGG(xpathway,cellPA)
P3=KEGG(xpathway,S1PA)
P4=KEGG(xpathway,S2PA)
P_kegg=cbind(P1,P2,P3,P4)
pp=-log10(P_kegg)

#################convert to function-gene pair
xx=which(lncBP$size>15&lncBP$size<500)
lncBP_1=lncBP[xx,]

##############kegg color
setwd("C:/projects/ceRNA/new")
esslnc=read.csv("neighbor.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
celgene=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
s1=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
s2=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
ss1=union(celgene$V1,s1$V1)
ss1=union(ss1,s2$V1)
allgene=union(ss1,esslnc$V1)
genecol=c()
for (i in 1:length(allgene)) {
  x1=which(esslnc$V1==allgene[i])
  x2=which(ss1==allgene[i])
  if(length(x1)>0&length(x2)>0){
    genecol=rbind(genecol,"red")
  } else if (length(x1)>0){
    genecol=rbind(genecol,"blue")
  } else {
    genecol=rbind(genecol,"yellow")
  }
}
Colgene=cbind(allgene,genecol)
write.table(Colgene,"genecol_KEGG.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

#########################essential gene and neighbor network
setwd("C:/projects/ceRNA/new")
neg=read.csv("neighbor.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
cell=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci1=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci2=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essentgene=union(cell$V1,sci1$V1)
essentgene=union(essentgene,sci2$V1)
Allg=union(neg$V1,essentgene)
SCIN=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
E_net=c()
for (i in 1:dim(SCIN)[1]) {
  print(i)
  x1=which(Allg==SCIN$V1[i])
  x2=which(Allg==SCIN$V2[i])
  if(length(x1)>0&length(x2)>0){
    E_net=rbind(E_net,SCIN[i,])
  }
}
xx=which(E_net$V1==E_net$V2)
E_net_N=E_net[-xx,]
write.table(E_net_N,"E_net_N.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
ID2sym=read.csv("ID_Sym_Uniprot.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
Ensg_net=c()
for (i in 1:dim(E_net_N)[1]) {
  print(i)
  xx1=which(ID2sym$HGNC.symbol==E_net_N$V1[i])
  xx2=which(ID2sym$HGNC.symbol==E_net_N$V2[i])
  if(length(xx1)>0&length(xx2)>0){
    ENSG1=unique(ID2sym$UniProt.SwissProt.Accession[xx1])
    ENSG2=unique(ID2sym$UniProt.SwissProt.Accession[xx2])
    if(length(ENSG1)>0&length(ENSG2)>0){
      Ensg_net=rbind(Ensg_net,cbind(ENSG1,ENSG2))
    }
  }
}
write.table(Ensg_net,"Uniprot_net_N.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
Allg=t(t(Allg))
EntrezID=merge(Allg,ID2sym,by.x="V1",by.y="HGNC.symbol")
write.table(EntrezID,"Essential_ID.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
############################functional enrichment analysis by script
setwd("C:/projects/ceRNA/new")
neg=read.csv("neighbor.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
cell=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci1=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci2=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)

funenrich=function(gene){
  Biofunction <- file("c2.cp.v5.2.symbols.gmt", "r")
  i<-0
  N <-19022
  line=readLines(Biofunction,n=1)
  Result <- c()
  while( length(line) != 0 ) {
    i <- i+1
    msig.go <- strsplit(line,"\t")
    msig.go <- msig.go[[1]]
    m <- length(msig.go)
    goterm <- msig.go[1]
    msig.go <- msig.go[3:m]
    n1 <- length(gene)
    intergene <- intersect(msig.go,gene)
    k1 <- length(intergene)
    print(k1)
    if (k1>2){
      p.value1 <- phyper(k1-1,m,N-m,n1,lower.tail=FALSE, log.p = FALSE)
      go <- c(goterm,k1,m,n1,p.value1)
      Result <- rbind(Result,go)
    } else {
      go <- c(goterm,k1,m,n1,1)
      Result <- rbind(Result,go)
    }
    line=readLines(Biofunction,n=1);
  }
  close(Biofunction)
  FDR=p.adjust(as.numeric(Result[,5]))
  Result=cbind(Result,FDR)
  return(Result)
}

Result_neg=funenrich(neg$V1)
Result_cell=funenrich(cell$V1)
Result_sci1=funenrich(sci1$V1)
Result_sci2=funenrich(sci2$V1)
x=which(as.numeric(Result_neg[,6])<0.05)
y1=Result_neg[x,]
x=which(as.numeric(Result_cell[,6])<0.05)
y2=Result_cell[x,]
x=which(as.numeric(Result_sci1[,6])<0.05)
y3=Result_sci1[x,]
x=which(as.numeric(Result_sci2[,6])<0.05)
y4=Result_sci2[x,]
II=intersect(y1[,1],y2[,1])
II=intersect(II,y3[,1])
II=intersect(II,y4[,1])
II=t(t(II))
write.table(Result_neg,"KEGG_lncRNA.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
write.table(Result_cell,"KEGG_cell.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
write.table(Result_sci1,"KEGG_sci1.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
write.table(Result_sci2,"KEGG_sci2.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)




####################randomly selected genes
Ens2sym=read.csv("Ens2sym.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Result_rand=c()
N=dim(Ens2sym)[1]
for (i in 1:100) {
  print(i)
  xx=sample(c(1:N),length(neg$V1))
  gene=Ens2sym$V2[xx]
  RR=funenrich(gene)
  Result_rand=cbind(Result_rand,RR[,6])
}
###############from hub gene
SCIN=read.csv("SIN.csv",stringsAsFactors=F,sep=",",skip=0,header = T)
hub=SCIN$name[1:1310]
Result_hub=funenrich(hub)
Nebhub=setdiff(neg$V1,hub)
Result_neghub=funenrich(Nebhub)
cellhub=setdiff(cell$V1,hub)
Result_cellhub=funenrich(cellhub)


####################Enrichment by TCGAbiolinks
setwd("C:/projects/ceRNA/new")
neg=read.csv("neighbor.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
cell=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci1=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci2=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
######
library("TCGAbiolinks")
system.time(ansEA <- TCGAanalyze_EAcomplete(TFname = "neg", cell$V1))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),GOBPTab = ansEA$ResBP, GOCCTab = ansEA$ResCC,GOMFTab = ansEA$ResMF,PathTab = ansEA$ResPat, nRGTab = cell$V1, nBar = 10)

########################
#############identify all the prognosistic genes in each cancer type
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")

Three.sam=function(fusion_sample){
  fusionthree=c()
  for (i in 1:length(fusion_sample)) {
    aa=strsplit(fusion_sample[i],"-")
    bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
    fusionthree=rbind(fusionthree,bb)
  }
  return(fusionthree)
}
library(survival)
for(kc in 29:33){
  print(kc)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep=""))
  library(SummarizedExperiment)
  LIHCmatrix<-assay(data,1,"FPKM")############gene exprssion
  LIHCclinal<-data@colData@listData$definition#############sample information tumor vs normal
  LIHCclinal=t(t(LIHCclinal))
  
  cancer.index=which(LIHCclinal!="Solid Tissue Normal")
  All.sample=colnames(LIHCmatrix)
  LIHCmatrix=LIHCmatrix[,cancer.index]
  All.sample=All.sample[cancer.index]
  All.sample=Three.sam(All.sample)
  N=dim(LIHCmatrix)[2]
  filter_line<-c()
  for (i in 1:dim(LIHCmatrix)[1]) {
    n=which(LIHCmatrix[i,]<0.1)
    if(length(n)>N*0.3){
      filter_line<-rbind(filter_line,i)
    }
  }
  Allgene=rownames(LIHCmatrix)
  setwd("C:/projects/AS/tcgaCLI")
  Clinical=read.csv(paste(cancer[kc],"_Mysample.txt",sep=""),stringsAsFactors = FALSE, header = T,sep = "\t")
  All.sample=cbind(All.sample,c(1:dim(All.sample)[1]))
  colnames(All.sample)=c("sample_name","SID")
  Sur.sam=merge(All.sample,Clinical,by.x="sample_name",by.y="ID1")
  LIHCmatrix=LIHCmatrix[,as.numeric(Sur.sam$SID)]
  Sur.sam$SurvObj <- with(Sur.sam, Surv(FinalTime, vital_status == "dead"))
  Gene.pro=c()
  for(i in 1:length(Allgene)){
    print(i)
    gex=log2(as.numeric(as.character(LIHCmatrix[i,]))+0.05)
    N=length(gex)
    xx=length(unique(gex))
    if(xx>N*0.3){
      train_input1=data.frame(score_train=as.numeric(as.character(gex)),
                              SurvObj=Sur.sam$SurvObj,stringsAsFactors =F)
      cox_model = coxph(SurvObj ~ score_train, data=train_input1)
      summcph <- summary(cox_model)
      SS=summcph$coefficients
      hr=as.numeric(SS[1])
      p.v=as.numeric(SS[5])
      Gene.pro=rbind(Gene.pro,cbind(Allgene[i],hr,p.v))
    }
  }
  FDR=p.adjust(as.numeric(Gene.pro[,3]),method = "BH")
  Gene.pro=cbind(Gene.pro,FDR)
  setwd("C:/projects/ceRNA/Sur")
  write.table(Gene.pro,paste(cancer[kc],"_Survival_HR.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
}
 

###############find the protective and risk genes and lncRNAs
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
setwd("C:/projects/ceRNA/Sur")
Num.sur=c()
Allpro=c()
Allris=c()
for(kc in 1:33){
  print(kc)
  Coxg=read.csv(paste(cancer[kc],"_Survival_HR.txt",sep = ""),stringsAsFactors = FALSE, header = F,sep = "\t")
  x1=which(Coxg$V3<0.01&Coxg$V2>0)
  Proc=Coxg[x1,]
  x2=which(Coxg$V3<0.01&Coxg$V2<0)
  Risk=Coxg[x2,]
  kk=cbind(cancer[kc],length(x1),length(x2))
  Num.sur=rbind(Num.sur,kk)
  Proc=cbind(Proc,cancer[kc])
  Risk=cbind(Risk,cancer[kc])
  Allpro=rbind(Allpro,Proc)
  Allris=rbind(Allris,Risk)
}

setwd("C:/projects/ceRNA/new")
esslnc=read.csv("essential lncRNA.txt",stringsAsFactors = FALSE, header = F,sep = "\t")
esslnc=unique(esslnc$V1)
core.cell=read.csv("corefitness.txt",stringsAsFactors = FALSE, header = F,sep = "\t")
core.s1=read.csv("core_1.txt",stringsAsFactors = FALSE, header = F,sep = "\t")
core.s2=read.csv("core_2.txt",stringsAsFactors = FALSE, header = F,sep = "\t")
ens2sym=read.csv("Ens2sym.txt",stringsAsFactors = FALSE, header = F,sep = "\t")
core.cell.id=merge(core.cell,ens2sym,by.x="V1",by.y="V2")
core.s1.id=merge(core.s1,ens2sym,by.x="V1",by.y="V2")
core.s2.id=merge(core.s2,ens2sym,by.x="V1",by.y="V2")

a1=intersect(esslnc,Allpro$V1)
a2=intersect(esslnc,Allris$V1)

b1=intersect(core.cell.id[,2],Allpro$V1)
b2=intersect(core.cell.id[,2],Allris$V1)

c1=intersect(core.s1.id[,2],Allpro$V1)
c2=intersect(core.s1.id[,2],Allris$V1)

d1=intersect(core.s2.id[,2],Allpro$V1)
d2=intersect(core.s2.id[,2],Allris$V1)

N.pro=unique(Allpro$V1)
N.risk=unique(Allris$V1)

s1=intersect(b1,c1)
s2=intersect(s1,d1)

s2=t(t(s2))
s2.name=merge(s2,ens2sym,by.x="V1",by.y="V1")
a1=t(t(a1))
a1.name=merge(a1,ens2sym,by.x="V1",by.y="V1")

#################################ppi Network analysis
setwd("C:/projects/ceRNA/new")

neg=read.csv("neighbor.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
cell=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci1=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
sci2=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essentgene=union(cell$V1,sci1$V1)
essentgene=union(essentgene,sci2$V1)
Allg=union(neg$V1,essentgene)
SCIN=read.csv("HI3_Y2H_102416.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
E_net=c()
for (i in 1:dim(SCIN)[1]) {
  print(i)
  x1=which(Allg==SCIN$Symbol.for.A[i])
  x2=which(Allg==SCIN$Symbol.for.B[i])
  if(length(x1)>0|length(x2)>0){
    E_net=rbind(E_net,SCIN[i,])
  }
}
xx=which(E_net$Symbol.for.A==E_net$Symbol.for.B)
E_net_N=E_net[-xx,]
write.table(E_net_N,"E_net_N_HI3_or.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

G1=unique(E_net_N$Symbol.for.A)
G2=unique(E_net_N$Symbol.for.B)
G=union(G1,G2)
g.Lable=c()
for(i in 1:length(G)){
  x1=which(essentgene==G[i])
  x2=which(neg$V1==G[i])
  if(length(x1)>0&length(x2)>0){
    g.Lable=rbind(g.Lable,3)
  } else if(length(x1)>0&length(x2)==0){
    g.Lable=rbind(g.Lable,2)
  } else if (length(x1)==0&length(x2)>0){
    g.Lable=rbind(g.Lable,1)
  } else {
    g.Lable=rbind(g.Lable,0)
  }
}
GG=cbind(G,g.Lable)
write.table(GG,"E_net_N_HI3_or_node.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

E_net=c()
for (i in 1:dim(SCIN)[1]) {
  print(i)
  x1=which(Allg==SCIN$Symbol.for.A[i])
  x2=which(Allg==SCIN$Symbol.for.B[i])
  if(length(x1)>0&length(x2)>0){
    E_net=rbind(E_net,SCIN[i,])
  }
}
xx=which(E_net$Symbol.for.A==E_net$Symbol.for.B)
E_net_N=E_net[-xx,]
N.real=dim(E_net_N)[1]
ORF=read.csv("human_ORF.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
N.rand=c()
for(r in 1:1000){
  print(r)
  XX=sample(c(1:14840),2887)
  Allg=ORF$Gene_symbol[XX]
  E_net=c()
  for (i in 1:dim(SCIN)[1]) {
    print(i)
    x1=which(Allg==SCIN$Symbol.for.A[i])
    x2=which(Allg==SCIN$Symbol.for.B[i])
    if(length(x1)>0&length(x2)>0){
      E_net=rbind(E_net,SCIN[i,])
    }
  }
  E_net_N=E_net[-xx,]
  N.rand=rbind(N.rand,dim(E_net_N)[1])
}
d <- density(N.rand)
par(mfrow=c(3,1))
plot(d,xlim=c(280,1300))
polygon(d, col="gray", border="gray",xlim=c(280,1300)) 
points(1237,0,type = "o")

############Essential lncRNA and gene survival analysis
setwd("C:/projects/ceRNA/new")
esslnc=read.csv("essential lncRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
library(SummarizedExperiment)
threesam=function(allBRCAsample){
  BRCAthree=c()
  for (i in 1:length(allBRCAsample)) {
    aa=strsplit(allBRCAsample[i],"-")
    bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
    BRCAthree=rbind(BRCAthree,bb)
  }
  return(BRCAthree)
}
library(survival)
Pan.P=c()
for(kc in 1:33){
  print(kc)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep = ""))
  BRCAmatrix<-assay(data,1,"FPKM")
  allBRCAsample=colnames(BRCAmatrix)
  C_N=data@colData@listData$definition#############sample information tumor vs normal
  Cancer_s=which(C_N!="Solid Tissue Normal")
  Cancer.mat=BRCAmatrix[,Cancer_s]
  Cancer.sam=colnames(Cancer.mat)
  Ca.three=threesam(Cancer.sam)
  allgene=rownames(Cancer.mat)
  P.lnc=c()
  for(i in 1:dim(esslnc)[1]){
    x1=which(allgene==esslnc$V1[i])
    if(length(x1)>0){
      a=which(as.numeric(Cancer.mat[x1,])>median(as.numeric(Cancer.mat[x1,])))
      b=which(as.numeric(Cancer.mat[x1,])<=median(as.numeric(Cancer.mat[x1,])))
      if(length(a)>0&length(b)>0) {
        Hi.sam=Ca.three[a]
        Lo.sam=Ca.three[b]
        Hi.sam=t(t(Hi.sam))
        Lo.sam=t(t(Lo.sam))
        SS=read.csv(paste("C:/projects/AS/revised/",cancer[kc],"_Mysample_revised.txt",sep =""),stringsAsFactors=F,sep="\t",skip=0,header = T)
        Hi.sur=merge(Hi.sam,SS,by.x="V1",by.y="sample_re")
        Lo.sur=merge(Lo.sam,SS,by.x="V1",by.y="sample_re")
        AA=cbind(Hi.sur$FinalTime,Hi.sur$vital_status,1)
        BB=cbind(Lo.sur$FinalTime,Lo.sur$vital_status,0)
        CC=rbind(AA,BB)
        fit<-survdiff(Surv(as.numeric(CC[,1]),CC[,2]=="dead")~CC[,3])
        p_surv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
        P.lnc=rbind(P.lnc,p_surv)
      } else {
        P.lnc=rbind(P.lnc,1)
      }
    } else {
      P.lnc=rbind(P.lnc,1)
    }
  }
  Pan.P=cbind(Pan.P,P.lnc)
}
Pan.P.n=matrix(as.numeric(Pan.P)<0.05,nrow = 161)*1
library(pheatmap)
pheatmap(Pan.P.n,cluster_cols = T,labels_col =cancer)
##############cell gene
esslnc=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
ens2sym=read.csv("Ens2sym.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
esslnc=merge(esslnc,ens2sym,by.x="V1",by.y="V2")
colnames(esslnc)=c("V2","V1")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
library(SummarizedExperiment)
threesam=function(allBRCAsample){
  BRCAthree=c()
  for (i in 1:length(allBRCAsample)) {
    aa=strsplit(allBRCAsample[i],"-")
    bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
    BRCAthree=rbind(BRCAthree,bb)
  }
  return(BRCAthree)
}
library(survival)
Pan.P=c()
for(kc in 1:33){
  print(kc)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep = ""))
  BRCAmatrix<-assay(data,1,"FPKM")
  allBRCAsample=colnames(BRCAmatrix)
  C_N=data@colData@listData$definition#############sample information tumor vs normal
  Cancer_s=which(C_N!="Solid Tissue Normal")
  Cancer.mat=BRCAmatrix[,Cancer_s]
  Cancer.sam=colnames(Cancer.mat)
  Ca.three=threesam(Cancer.sam)
  allgene=rownames(Cancer.mat)
  P.lnc=c()
  for(i in 1:dim(esslnc)[1]){
    x1=which(allgene==esslnc$V1[i])
    if(length(x1)>0){
      a=which(as.numeric(Cancer.mat[x1,])>median(as.numeric(Cancer.mat[x1,])))
      b=which(as.numeric(Cancer.mat[x1,])<=median(as.numeric(Cancer.mat[x1,])))
      if(length(a)>0&length(b)>0) {
        Hi.sam=Ca.three[a]
        Lo.sam=Ca.three[b]
        Hi.sam=t(t(Hi.sam))
        Lo.sam=t(t(Lo.sam))
        SS=read.csv(paste("C:/projects/AS/revised/",cancer[kc],"_Mysample_revised.txt",sep =""),stringsAsFactors=F,sep="\t",skip=0,header = T)
        Hi.sur=merge(Hi.sam,SS,by.x="V1",by.y="sample_re")
        Lo.sur=merge(Lo.sam,SS,by.x="V1",by.y="sample_re")
        if(dim(Hi.sur)[1]>0&dim(Lo.sur)[1]>0){
          AA=cbind(Hi.sur$FinalTime,Hi.sur$vital_status,1)
          BB=cbind(Lo.sur$FinalTime,Lo.sur$vital_status,0)
          CC=rbind(AA,BB)
          fit<-survdiff(Surv(as.numeric(CC[,1]),CC[,2]=="dead")~CC[,3])
          p_surv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
          P.lnc=rbind(P.lnc,p_surv)
        } else {
          P.lnc=rbind(P.lnc,1)
        }
      } else {
        P.lnc=rbind(P.lnc,1)
      }
    } else {
      P.lnc=rbind(P.lnc,1)
    }
  }
  Pan.P=cbind(Pan.P,P.lnc)
}
Pan.P.n=matrix(as.numeric(Pan.P)<0.05,nrow = 1927)*1 ####1724,1885,1927
library(pheatmap)
pheatmap(Pan.P.n,cluster_cols = T,labels_col =cancer)
#######CANCER PROPORTION
load("C:/projects/ceRNA/new/S2_ess.RData")
A3=apply(Pan.P.n,2,sum)
A3=A3/1927
load("C:/projects/ceRNA/new/S1_ess.RData")
A2=apply(Pan.P.n,2,sum)
A2=A2/1885
load("C:/projects/ceRNA/new/cell_ess.RData")
A1=apply(Pan.P.n,2,sum)
A1=A1/1724
load("C:/projects/ceRNA/new/lncRNA_survival_pan.RData")
A0=apply(Pan.P.n,2,sum)
A0=A0/161
AA=rbind(A3,A2,A1,A0)
PPP=read.csv("Pan.P.txt",
             stringsAsFactors=F,sep="\t",skip=0,header = T)
PP1=c()
for(i in 1:33){
  print(i)
  RR=c()
  XX=as.numeric(as.character(PPP[,i]))
  N=length(XX)
  for(j in 1:100){
    ww=sample(c(1:N),1927)
    n=sum(XX[ww]<0.05)
    RR=rbind(RR,n/1927)
  }
  p1=sum(RR>A3[i])/100
  PP1=cbind(PP1,p1)
}
PP2=c()
for(i in 1:33){
  print(i)
  RR=c()
  XX=as.numeric(as.character(PPP[,i]))
  N=length(XX)
  for(j in 1:100){
    ww=sample(c(1:N),1885)
    n=sum(XX[ww]<0.05)
    RR=rbind(RR,n/1885)
  }
  p1=sum(RR>A2[i])/100
  PP2=cbind(PP2,p1)
}
PP3=c()
for(i in 1:33){
  print(i)
  RR=c()
  XX=as.numeric(as.character(PPP[,i]))
  N=length(XX)
  for(j in 1:100){
    ww=sample(c(1:N),1724)
    n=sum(XX[ww]<0.05)
    RR=rbind(RR,n/1724)
  }
  p1=sum(RR>A1[i])/100
  PP3=cbind(PP3,p1)
}
PP4=c()
for(i in 1:33){
  print(i)
  RR=c()
  XX=as.numeric(as.character(PPP[,i]))
  N=length(XX)
  for(j in 1:100){
    ww=sample(c(1:N),161)
    n=sum(XX[ww]<0.05)
    RR=rbind(RR,n/161)
  }
  p1=sum(RR>A0[i])/100
  PP4=cbind(PP4,p1)
}
P.rand=rbind(PP1,PP2,PP3,PP4)

ff.S=c()
for(i in 1:4){
  for(j in 1:33){
    aa=cbind(i,j,AA[i,j])
    ff.S=rbind(ff.S,aa)
  }
}
par(mfrow=c(3,1))
symbols(ff.S[,2],ff.S[,1],ff.S[,3],inches = .25,bg="steelblue2", fg=NULL,
        xlim = c(0,34))
grid(nx = NULL, ny =NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

library(corrplot)
corrplot(AA, ccol = terrain.colors(10), cl.lim = c(0, 1))
corrplot(AA, method = "number",col="black", cl.lim = c(0, 1))

#####TOP lncRNA and genes
load("C:/projects/ceRNA/new/S2_ess.RData")
A3=apply(Pan.P.n,1,sum)
esslnc=cbind(esslnc,A3)
xx1=sort.int(as.numeric(esslnc$A3),decreasing = T,index.return = T)
barplot(xx1$x[10:1],horiz = T)
nn=esslnc$V2[xx1$ix[10:1]]
load("C:/projects/ceRNA/new/S1_ess.RData")
A2=apply(Pan.P.n,1,sum)
esslnc=cbind(esslnc,A2)
xx1=sort.int(as.numeric(esslnc$A2),decreasing = T,index.return = T)
barplot(xx1$x[10:1],horiz = T)
nn=esslnc$V2[xx1$ix[10:1]]
nn
load("C:/projects/ceRNA/new/cell_ess.RData")
A2=apply(Pan.P.n,1,sum)
esslnc=cbind(esslnc,A2)
xx1=sort.int(as.numeric(esslnc$A2),decreasing = T,index.return = T)
barplot(xx1$x[10:1],horiz = T)
nn=esslnc$V2[xx1$ix[1:10]]
load("C:/projects/ceRNA/new/lncRNA_survival_pan.RData")
A0=apply(Pan.P.n,1,sum)
esslnc=cbind(esslnc,A0)

esslnc=esslnc[,c(1,3)]
esslnc=esslnc[!duplicated(esslnc),]
xx1=sort.int(as.numeric(esslnc$A0),decreasing = T,index.return = T)
barplot(xx1$x[10:1],horiz = T)
nn=esslnc$V1[xx1$ix[1:10]]
nn

#####ENSG00000255717
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
esslnc="ENSG00000265206"####"ENSG00000249859"
esslnc=t(t(esslnc))
library(SummarizedExperiment)
threesam=function(allBRCAsample){
  BRCAthree=c()
  for (i in 1:length(allBRCAsample)) {
    aa=strsplit(allBRCAsample[i],"-")
    bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
    BRCAthree=rbind(BRCAthree,bb)
  }
  return(BRCAthree)
}
library(survival)
Pan.P=c()
for(kc in 1:33){
  print(kc)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep = ""))
  BRCAmatrix<-assay(data,1,"FPKM")
  allBRCAsample=colnames(BRCAmatrix)
  C_N=data@colData@listData$definition#############sample information tumor vs normal
  Cancer_s=which(C_N!="Solid Tissue Normal")
  Cancer.mat=BRCAmatrix[,Cancer_s]
  Cancer.sam=colnames(Cancer.mat)
  Ca.three=threesam(Cancer.sam)
  allgene=rownames(Cancer.mat)
  P.lnc=c()
  for(i in 1:length(esslnc)){
    x1=which(allgene==esslnc[i])
    if(length(x1)>0){
      a=which(as.numeric(Cancer.mat[x1,])>median(as.numeric(Cancer.mat[x1,])))
      b=which(as.numeric(Cancer.mat[x1,])<=median(as.numeric(Cancer.mat[x1,])))
      if(length(a)>0&length(b)>0) {
        Hi.sam=Ca.three[a]
        Lo.sam=Ca.three[b]
        Hi.sam=t(t(Hi.sam))
        Lo.sam=t(t(Lo.sam))
        SS=read.csv(paste("C:/projects/AS/revised/",cancer[kc],"_Mysample_revised.txt",sep =""),stringsAsFactors=F,sep="\t",skip=0,header = T)
        Hi.sur=merge(Hi.sam,SS,by.x="V1",by.y="sample_re")
        Lo.sur=merge(Lo.sam,SS,by.x="V1",by.y="sample_re")
        if(dim(Hi.sur)[1]>0&dim(Lo.sur)[1]>0){
          AA=cbind(Hi.sur$FinalTime,Hi.sur$vital_status,1)
          BB=cbind(Lo.sur$FinalTime,Lo.sur$vital_status,0)
          CC=rbind(AA,BB)
          ddd<-survfit(Surv(as.numeric(CC[,1]),CC[,2]=="dead")~CC[,3])
          plot(ddd,col = c("red","blue"))
          title(cancer[kc])
          fit<-survdiff(Surv(as.numeric(CC[,1]),CC[,2]=="dead")~CC[,3])
          p_surv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
          P.lnc=rbind(P.lnc,p_surv)
        } else {
          P.lnc=rbind(P.lnc,1)
        }
      } else {
        P.lnc=rbind(P.lnc,1)
      }
    } else {
      P.lnc=rbind(P.lnc,1)
    }
  }
  Pan.P=cbind(Pan.P,P.lnc)
}

##############predict the essential genes
setwd("C:/projects/ceRNA/new")
essentgene=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
SCIN=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene1=unique(SCIN$V1)
gene2=unique(SCIN$V2)
gene=union(gene1,gene2)
Ess.pro=c()
for(i in 1:length(gene)){
  print(i)
  xx=which(SCIN$V1==gene[i]|SCIN$V2==gene[i])
  sub=SCIN[xx,]
  X1=unique(sub[,1])
  X2=unique(sub[,2])
  XX=union(X1,X2)
  XX=setdiff(XX,gene[i])
  N=length(XX)
  n=length(intersect(XX,essentgene$V1))
  if(N>0){
    kk=cbind(gene[i],n/N,N)
    Ess.pro=rbind(Ess.pro,kk)
  }
}
##########ROC curve
Lab=c()
for(i in 1:dim(Ess.pro)[1]){
  xx=which(essentgene$V1==Ess.pro[i,1])
  if(length(xx)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
Ess.pro=cbind(Ess.pro,Lab)
library(ROCR)
pred_PRO1 <- prediction(as.numeric(Ess.pro[,2]), Ess.pro[,4])
perf_PRO1 <- performance(pred_PRO1,"tpr","fpr")
plot(perf_PRO1,ylim=c(0,1),lwd=2,lty=1)
AUC_PRO1<-performance(pred_PRO1,"auc")

pred_deg1 <- prediction(as.numeric(Ess.pro[,3]), Ess.pro[,4])
perf_deg1 <- performance(pred_deg1,"tpr","fpr")
plot(perf_deg1,ylim=c(0,1),lwd=2,lty=1,add=TRUE, col="MediumSeaGreen")
AUC_deg1<-performance(pred_deg1,"auc")

setwd("C:/projects/ceRNA/new")
essentgene=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
SCIN=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene1=unique(SCIN$V1)
gene2=unique(SCIN$V2)
gene=union(gene1,gene2)
Ess.pro=c()
for(i in 1:length(gene)){
  print(i)
  xx=which(SCIN$V1==gene[i]|SCIN$V2==gene[i])
  sub=SCIN[xx,]
  X1=unique(sub[,1])
  X2=unique(sub[,2])
  XX=union(X1,X2)
  XX=setdiff(XX,gene[i])
  N=length(XX)
  n=length(intersect(XX,essentgene$V1))
  if(N>0){
    kk=cbind(gene[i],n/N,N)
    Ess.pro=rbind(Ess.pro,kk)
  }
}
##########ROC curve
Lab=c()
for(i in 1:dim(Ess.pro)[1]){
  xx=which(essentgene$V1==Ess.pro[i,1])
  if(length(xx)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
Ess.pro=cbind(Ess.pro,Lab)

pred_PRO2 <- prediction(as.numeric(Ess.pro[,2]), Ess.pro[,4])
perf_PRO2 <- performance(pred_PRO2,"tpr","fpr")
plot(perf_PRO2,ylim=c(0,1),lwd=2,lty=1,add=TRUE, col="red")
AUC_PRO2<-performance(pred_PRO2,"auc")

pred_deg2 <- prediction(as.numeric(Ess.pro[,3]), Ess.pro[,4])
perf_deg2 <- performance(pred_deg2,"tpr","fpr")
plot(perf_deg2,ylim=c(0,1),lwd=2,lty=1,add=TRUE, col="gray")
AUC_deg2<-performance(pred_deg2,"auc")

setwd("C:/projects/ceRNA/new")
essentgene=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
SCIN=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene1=unique(SCIN$V1)
gene2=unique(SCIN$V2)
gene=union(gene1,gene2)
Ess.pro=c()
for(i in 1:length(gene)){
  print(i)
  xx=which(SCIN$V1==gene[i]|SCIN$V2==gene[i])
  sub=SCIN[xx,]
  X1=unique(sub[,1])
  X2=unique(sub[,2])
  XX=union(X1,X2)
  XX=setdiff(XX,gene[i])
  N=length(XX)
  n=length(intersect(XX,essentgene$V1))
  if(N>0){
    kk=cbind(gene[i],n/N,N)
    Ess.pro=rbind(Ess.pro,kk)
  }
}
##########ROC curve
Lab=c()
for(i in 1:dim(Ess.pro)[1]){
  xx=which(essentgene$V1==Ess.pro[i,1])
  if(length(xx)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
Ess.pro=cbind(Ess.pro,Lab)

pred_PRO3 <- prediction(as.numeric(Ess.pro[,2]), Ess.pro[,4])
perf_PRO3 <- performance(pred_PRO3,"tpr","fpr")
plot(perf_PRO3,ylim=c(0,1),lwd=2,lty=1,add=TRUE, col="yellow")
AUC_PRO3<-performance(pred_PRO3,"auc")

pred_deg3 <- prediction(as.numeric(Ess.pro[,3]), Ess.pro[,4])
perf_deg3 <- performance(pred_deg3,"tpr","fpr")
plot(perf_deg3,ylim=c(0,1),lwd=2,lty=1,add=TRUE, col="blue")
AUC_deg3<-performance(pred_deg3,"auc")

barplot(c(0.773,0.704,0.792,0.713,0.758,0.709),ylim = c(0.65,0.8))

##############a b c 
a=seq(0,1,0.1)
b=seq(0,1,0.1)
RR=c()
for(i in 1:11){
  for(j in 1:11){
    RR=rbind(RR,cbind(a[i],b[j],as.numeric(1-a[i]-b[j])))
  }
}
B=apply(RR,1,sum)
x=which(B==1&RR[,1]>=0&RR[,1]<=1&RR[,2]>=0&RR[,2]<=1&RR[,3]>=0&RR[,3]<=1)
RR_t=RR[x,]
##################Prediction
setwd("C:/projects/ceRNA/new")
essentgene0=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essentgene1=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essentgene2=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
essentgene=union(essentgene0$V1,essentgene1$V1)
essentgene=union(essentgene,essentgene2$V1)
essentgene=data.frame(V1=essentgene)
SCIN=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene1=unique(SCIN$V1)
gene2=unique(SCIN$V2)
gene=union(gene1,gene2)
Ess.pro=c()
for(i in 1:length(gene)){
  print(i)
  xx=which(SCIN$V1==gene[i]|SCIN$V2==gene[i])
  sub=SCIN[xx,]
  X1=unique(sub[,1])
  X2=unique(sub[,2])
  XX=union(X1,X2)
  XX=setdiff(XX,gene[i])
  N=length(XX)
  n=length(intersect(XX,essentgene$V1))
  if(N>0){
    kk=cbind(gene[i],n,n/N,N)
    Ess.pro=rbind(Ess.pro,kk)
  }
}
##########ROC curve
Lab=c()
for(i in 1:dim(Ess.pro)[1]){
  xx=which(essentgene$V1==Ess.pro[i,1])
  if(length(xx)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
Ess.pro=cbind(Ess.pro,Lab)
###################gene expression
#corefit=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header =F )
corefit=essentgene
setwd("C:/projects/genefusion/data/HPA")
Rnaexp=read.csv("rna_tissue.csv",stringsAsFactors=F,sep=",",skip=0,header = T)
Gene=unique(Rnaexp$Gene.name)
Avgexp=c()
for (i in 1:length(Gene)) {
  print(i)
  xx=which(Rnaexp$Gene.name==Gene[i])
  GGexp=Rnaexp$Value[xx]
  GGexp=as.numeric(as.character(GGexp))
  Avgexp=rbind(Avgexp,cbind(Gene[i],mean(GGexp)))
}

Feat.mat=merge(Avgexp,Ess.pro,by.x="V1",by.y="V1")
Feat.mat$V2=rank(as.numeric(as.character(Feat.mat$V2)))
Feat.mat$V2=Feat.mat$V2/max(Feat.mat$V2)

Feat.mat$V3=rank(as.numeric(as.character(Feat.mat$V3)))
Feat.mat$V3=Feat.mat$V3/max(Feat.mat$V3)

Feat.mat$N=rank(as.numeric(as.character(Feat.mat$N)))
Feat.mat$N=Feat.mat$N/max(Feat.mat$N)
colnames(Feat.mat)=c("gene","exp_rank","numofess","pro_rank","degree_rank","Lab")
F.score=c()
for(i in 1:dim(Feat.mat)[1]){
  a=as.numeric(Feat.mat$exp_rank[i])
  b=as.numeric(Feat.mat$pro_rank[i])
  c=as.numeric(Feat.mat$degree_rank[i])
  FF=0.2*a+0.6*b+0.2*c
  F.score=rbind(F.score,FF)
}
Feat.mat=cbind(Feat.mat,F.score)
library(ROCR)
pred_PRO4 <- prediction(as.numeric(F.score), Feat.mat$Lab)
perf_PRO4 <- performance(pred_PRO4,"tpr","fpr")
plot(perf_PRO4,ylim=c(0,1),lwd=2,lty=1, col="red")
AUC_PRO4<-performance(pred_PRO4,"auc")

#############com---red---0.814
pred_PRO1 <- prediction(as.numeric(Feat.mat$degree_rank), Feat.mat$Lab)
perf_PRO1 <- performance(pred_PRO1,"tpr","fpr")
plot(perf_PRO1,ylim=c(0,1),lwd=2,lty=1, col="yellow",add=TRUE)
AUC_PRO1<-performance(pred_PRO1,"auc")
###########exp--green--0.665
###########pro--blue---0.776
##########degree---yellow---0.705

barplot(c(0.814,0.776,0.705,0.665),ylim = c(0.6,0.83))
setwd("C:/projects/ceRNA/new")
write.table(Feat.mat,"Final_AUC1230.txt",
            sep = "\t",row.names = F, col.names = T,quote = FALSE)

############################validation
Fscore=read.csv("Final_AUC1230.txt",
                stringsAsFactors=F,sep="\t",skip=0,header = T)

A33=read.csv("C:/projects/ceRNA/new/33cell/Achilles_v3.3.8.gct",
             stringsAsFactors=F,sep="\t",skip=2,header = T)
WW=merge(Fscore,A33,by.x = "gene",by.y = "Description")
x=as.numeric(as.character(WW$F.score))
RR=c()
for(i in 9:41){
  y=as.numeric(as.character(WW[,i]))
  r1=cor(x,y,method = "spearman")
  p1=cor.test(x,y)$p.value
  RR=rbind(RR,cbind(r1,p1))
}
RR=data.frame(cancer=colnames(A33)[3:35],
              x1=as.numeric(as.character(RR[,1])),
              x2=as.numeric(as.character(RR[,2])))
RR$x2=0.01
RR$x1=-RR$x1
RR$Pl=-log10(RR$x2)
library(ggplot2)
ggplot(RR, aes(x=reorder(cancer,x1,FUN = abs), y=x1)) + 
  geom_point(aes(size=Pl),col="green") + 
  geom_segment(aes(x=cancer, 
                   xend=cancer, 
                   y=0, 
                   yend=x1))+ylim(0,0.22) 
RK=Fscore[,c(1,7)]
write.table(RK,"Final_Ess.rnk",
            sep = "\t",row.names = F, col.names = F,quote = FALSE)


Pos=which(Feat.mat$V4==1)
Neg=which(Feat.mat$V4==0)
N.pos=length(Pos)
N.neg=length(Neg)
library(ROCR)
AUC_rand=c()
for(j in 1:100){
  xx=sample(c(1:N.neg),N.pos)
  Feat.mat.r=rbind(Feat.mat[Pos,],Feat.mat[Neg[xx],])
  AUC_k=c()
  for(i in 1:dim(RR_t)[1]){
    print(i)
    KW=as.numeric(as.character(RR_t[i,1]))*as.numeric(as.character(Feat.mat.r[,2]))+
      as.numeric(as.character(RR_t[i,2]))*as.numeric(as.character(Feat.mat.r[,3]))+
      as.numeric(as.character(RR_t[i,3]))*as.numeric(as.character(Feat.mat.r[,4]))
    Lab=as.numeric(as.character(Feat.mat.r[,5]))
    pred_deg3 <- prediction(KW, Lab)
    perf_deg3 <- performance(pred_deg3,"tpr","fpr")
    #plot(perf_deg3,ylim=c(0,1),lwd=2,lty=1,add=TRUE, col="blue")
    AUC_deg3<-performance(pred_deg3,"auc")
    AUC=AUC_deg3@y.values[[1]]
    AUC_k=rbind(AUC_k,AUC)
  }
  AUC_rand=cbind(AUC_rand,AUC_k)
}
rownames(AUC_rand)=paste("[",RR_t[,1],",",RR_t[,2],",",RR_t[,3],"]",sep="")
AUC_rand=t(AUC_rand)
Med.AUC=apply(AUC_rand,2,median)
ks=sort.int(Med.AUC,index.return = T)
boxplot(AUC_rand[,ks$ix],las=2)



xx=which(AUC_k==max(AUC_k))
AUC_thr=cbind(RR_t,AUC_k)
library(pheatmap)
pheatmap(AUC_thr,cluster_rows = F,cluster_cols = F)


library(rgl)
ww=sort.int(as.numeric(as.character(AUC_thr[,4])),index.return = T)
ww2=ww$ix
AUC_thr=AUC_thr[ww2,]
library(pheatmap)
pheatmap(AUC_thr,cluster_rows = F,cluster_cols = F)
x=AUC_thr[,1]
y=AUC_thr[,2]
z=AUC_thr[,3]

c = AUC_thr[,4]
c = cut(c, breaks=60)
cols = redblue(60)[as.numeric(c)]
plot3d(x,y,z,col=cols,size = 15)
texts3d(x,y,z,AUC_thr[,4])

diamonds=data.frame(V1=as.numeric(as.character(Feat.mat[,2])),V2=as.numeric(as.character(Feat.mat[,3])),
                    V3=as.numeric(as.character(Feat.mat[,4])),V4=as.numeric(as.character(Feat.mat[,5])))
library(caret)
model <- train(
  V4 ~ ., diamonds,
  method = "lm",
  trControl = trainControl(
    method = "cv", number = 10,
    verboseIter = TRUE
  )
)


###################################Figure-1a
setwd("C:/projects/ceRNA/figure1")
Achi=read.csv("Achilles_v3.3.8_genesoln.gct",
              stringsAsFactors=F,sep="\t",skip=2,header = T)
Ugene=unique(Achi$Description)
Ugene=t(t(Ugene))
Gene.loc=read.csv("gene_loc.txt",
                  stringsAsFactors=F,sep="\t",skip=0,header = T)
Ess.gene=merge(Ugene,Gene.loc,by.x="V1",by.y="HGNC.symbol")
x=(grep("CH",Ess.gene$Chromosome.scaffold.name))
Ess.gene=Ess.gene[-x,]
Leg=c()
for(i in 1:dim(Ess.gene)[1]){
  print(i)
  st=as.numeric(as.character(Ess.gene$Gene.start..bp.[i]))
  ed=as.numeric(as.character(Ess.gene$Gene.end..bp.[i]))
  ll=abs(ed-st)
  Leg=rbind(Leg,ll)
}
colnames(Leg)="length"
Ess.gene=cbind(Ess.gene,Leg)
Ess=Ess.gene[,c(1,5)]
colnames(Ess)=c("gene","length")

##########################TCGA mutations
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
N.cancer=c()
for(kc in 1:33){
  print(kc)
  mut=read.csv(paste("C:/projects/DNA_repair/Mutation/MC3/",cancer[kc],"_mut.txt",sep=""),
               stringsAsFactors=F,sep="\t",skip=0,header =T)
  CC=c()
  for(i in 1:dim(Ess)[1]){
    xa=which(mut$Hugo_Symbol==Ess$gene[i])
    n1=length(xa)/Ess$length[i]
    CC=rbind(CC,n1)
  }
  N.cancer=cbind(N.cancer,CC)
}
colnames(N.cancer)=cancer
Mut.cancer=cbind(Ess,N.cancer)
Avg.m1=apply(Mut.cancer[,3:35],1,sum)
Avg.m1=Avg.m1*1000
WW=cbind(as.character(Mut.cancer[,1]),as.numeric(Avg.m1))
write.table(WW,file = "TCGA_gene.rnk",sep = "\t",row.names = F, col.names = F,quote = FALSE)
#######################essential gmt
fileConn<-file("Ess.gmt","a")
cor1=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor2=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor3=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor1=t(cor1)
cor2=t(cor2)
cor3=t(cor3)
CC1=intersect(cor1,cor2)
CC2=intersect(cor3,CC1)
kw=cbind("Ess1","Ess1",cor1)
#cat(kw)
cat(kw,file=fileConn,append = TRUE,sep = "\t")
cat("\n",file=fileConn,append=T)
kw=cbind("Ess2","Ess2",cor2)
#cat(kw)
cat(kw,file=fileConn,append = TRUE,sep = "\t")
cat("\n",file=fileConn,append=T)
kw=cbind("Ess3","Ess3",cor3)
#cat(kw)
cat(kw,file=fileConn,append = TRUE,sep = "\t")
cat("\n",file=fileConn,append=T)

CC2=t(CC2)
kw=cbind("Ess4","Ess4",CC2)
#cat(kw)
cat(kw,file=fileConn,append = TRUE,sep = "\t")
cat("\n",file=fileConn,append=T)
close(fileConn)

for(i in 1:length(cat)){
  #print(i)
  fileConn<-file("Ess.gmt","a")
  mm=dim(FT.sam)[1]
  xx=which(FT.sam[,i]>0.001)
  ge=unique(rownames(FT.sam)[xx])
  Li=paste("patient",i,sep = "")
  kw=cbind(Li,cat[i],t(ge))
  #cat(kw)
  cat(kw,file=fileConn,append = TRUE,sep = "\t")
  cat("\n",file=fileConn,append=T)
  close(fileConn)
}

#######


Achi=Achi[,-1]
Mut.ess=merge(Mut.cancer,Achi,by.x = "gene",by.y = "Description")
Avg.e=apply(-Mut.ess[,36:68],1,min)
Avg.m=apply(Mut.ess[,3:35],1,sum)
Avg.m=Avg.m*1000
cor(Avg.e,Avg.m)
library(Hmisc)
Wg=as.numeric(cut2(Avg.e, g=10))
i=c(1:10)
XX=tapply(as.numeric(Avg.m), Wg, mean)
scatter.smooth(i,XX)


###################ExAc
setwd("C:/projects/ceRNA/figure1")
Achi=read.csv("Achilles_v3.3.8_genesoln.gct",
              stringsAsFactors=F,sep="\t",skip=2,header = T)
Ugene=unique(Achi$Description)
Ugene=t(t(Ugene))
Gene.loc=read.csv("gene_loc.txt",
                  stringsAsFactors=F,sep="\t",skip=0,header = T)
Ess.gene=merge(Ugene,Gene.loc,by.x="V1",by.y="HGNC.symbol")
x=(grep("CH",Ess.gene$Chromosome.scaffold.name))
Ess.gene=Ess.gene[-x,]
Leg=c()
for(i in 1:dim(Ess.gene)[1]){
  print(i)
  st=as.numeric(as.character(Ess.gene$Gene.start..bp.[i]))
  ed=as.numeric(as.character(Ess.gene$Gene.end..bp.[i]))
  ll=abs(ed-st)
  Leg=rbind(Leg,ll)
}
colnames(Leg)="length"
Ess.gene=cbind(Ess.gene,Leg)
Ess=Ess.gene[,c(1,5)]
colnames(Ess)=c("gene","length")

Exac=read.csv("fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt",
              stringsAsFactors=F,sep="\t",skip=0,header = T)
NN=c()
for(i in 1:dim(Ess)[1]){
  xx=which(Exac$gene==Ess$gene[i])
  if(length(xx)>0){
    k1=Exac$n_syn[xx]/Ess$length[i]
    k2=Exac$n_mis[xx]/Ess$length[i]
  } else {
    k1=0
    k2=0
  }
  K=cbind(k1,k2)
  NN=rbind(NN,K)
}
Ess=cbind(Ess,NN)

Achi=Achi[,-1]
Mut.ess=merge(Ess,Achi,by.x = "gene",by.y = "Description")
Avg.e=apply(-Mut.ess[,5:37],1,min)
Avg.m=Mut.ess$k2
Avg.m=Avg.m*1000
cor(Avg.e,Avg.m)
library(Hmisc)
Wg=as.numeric(cut2(Avg.e, g=100))
i=c(1:100)
XX=tapply(as.numeric(Avg.m), Wg, mean)
scatter.smooth(i,XX,ylim = c(2,12))




Clivar=read.csv("HGVS_summary.txt",
                stringsAsFactors=F,sep="\t",skip=0,header = F)
Pdeg=c()
for(i in 1:dim(Ess)[1]){
  print(i)
  xx=which(Clivar$V2==Ess$gene[i])
  n1=length(unique(Clivar$V6[xx]))
  deg=n1/(as.numeric(as.character(Ess$length[i])))
  Pdeg=rbind(Pdeg,deg)
}
Ess=cbind(Ess,Pdeg)
Ess$Pdeg=Ess$Pdeg*1000
Achi=Achi[,-1]
Achi.ess=merge(Achi,Ess,by.x="Description",by.y="gene")
i=c(1:10)
plot(i,i,type="l")
library(Hmisc)
for(j in 1:33){
  Achi.ess$group=as.numeric(cut2(Achi.ess[,j+1], g=50))
  XX=tapply(Achi.ess$Pdeg, Achi.ess$group, mean)
  qplot(i,XX, geom='smooth', span =0.9)
}

WW=apply(-Achi.ess[,2:34],1,median)
Wg=as.numeric(cut2(WW, g=100))
i=c(1:100)
XX=tapply(as.numeric(Achi.ess$Pdeg), Wg, mean)
scatter.smooth(i,XX,ylim = c(0.2,0.3))

R=c()
for(j in 1:33){
r=cor(as.numeric(as.character(Achi.ess[,j+1])),as.numeric(as.character(Achi.ess$Pdeg)))
R=rbind(R,r)

}

##################gene mutation frequency
setwd("C:/projects/ceRNA/figure1")
gene=read.csv("gencode_gene_pro.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
N.cancer=c()
for(kc in 1:33){
  print(kc)
  mut=read.csv(paste("C:/projects/DNA_repair/Mutation/MC3/",cancer[kc],"_mut.txt",sep=""),
               stringsAsFactors=F,sep="\t",skip=0,header =T)
  C.n=length(mut$Tumor_Sample_Barcode)
  CC=c()
  for(i in 1:dim(gene)[1]){
    xa=which(mut$Hugo_Symbol==gene$V2[i])
    n1=length(unique(mut$Tumor_Sample_Barcode[xa]))
    rr=n1/C.n
    CC=rbind(CC,rr)
  }
  N.cancer=cbind(N.cancer,CC)
}
Avg.e=apply(N.cancer,1,mean)
gene=cbind(gene,Avg.e)
write.table(gene,"TCGA_gene_allfreq.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

cor1=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor2=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor3=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor1=t(cor1)
cor2=t(cor2)
cor3=t(cor3)
CC1=intersect(cor1,cor2)
CC2=intersect(cor3,CC1)

library(Hmisc)
Wg=as.numeric(cut2(Avg.e, g=5))
boxplot(Avg.e~Wg)
#gene=cbind(gene,Wg)
R.r=c()
for(i in 1:5){
  x1=which(Wg==i)
  M=length(x1)
  aa=intersect(gene[x1,2],CC2)
  m=length(aa)
  R.r=rbind(R.r,m/M)
}
R.r_s=1-R.r[,1]
R.r=cbind(R.r,R.r_s)
barplot(t(R.r))
Lab=c()
for(i in 1:dim(gene)[1]){
  xx=which(CC2==gene[i,2])
  if(length(xx)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
x1=which(Lab==1)
Y1=as.numeric(gene[x1,3])
x2=which(Lab==0)
Y2=as.numeric(gene[x2,3])
boxplot(gene[,3]~Lab)
library(vioplot)
boxplot(Y1*1000,Y2*1000,outline = F)
vioplot::vioplot(Y1*1000,Y2*1000)
wilcox.test(Y1,Y2,alternative = "greater")

####################################
setwd("C:/projects/ceRNA/figure1")
gene=read.csv("gencode_gene_pro.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
N.cancer=c()
for(kc in 1:33){
  print(kc)
  mut=read.csv(paste("C:/projects/DNA_repair/Mutation/MC3/",cancer[kc],"_mut.txt",sep=""),
               stringsAsFactors=F,sep="\t",skip=0,header =T)
  mm2=read.csv(paste("C:/projects/DNA_repair/Mutation/Annovar/",cancer[kc],"_annovar_input.txt.hg19_multianno.txt",sep=""),
               stringsAsFactors=F,sep="\t",skip=0,header =T)
  mut=cbind(mut,mm2)
  xx=which(mut$MutationTaster_pred=="D")
  C.n=length(mut$Tumor_Sample_Barcode)
  CC=c()
  for(i in 1:dim(gene)[1]){
    xa=which(mut$Hugo_Symbol==gene$V2[i])
    xa=intersect(xa,xx)
    n1=length(unique(mut$Tumor_Sample_Barcode[xa]))
    rr=n1/C.n
    CC=rbind(CC,rr)
  }
  N.cancer=cbind(N.cancer,CC)
}
Avg.e=apply(N.cancer,1,mean)
gene=cbind(gene,Avg.e)
write.table(gene,"TCGA_gene_Delfreq_MutationTaster.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

cor1=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor2=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor3=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor1=t(cor1)
cor2=t(cor2)
cor3=t(cor3)
CC1=intersect(cor1,cor2)
CC2=intersect(cor3,CC1)

Lab=c()
for(i in 1:dim(gene)[1]){
  xx=which(CC2==gene[i,2])
  if(length(xx)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
x1=which(Lab==1)
Y1=as.numeric(gene[x1,3])
x2=which(Lab==0)
Y2=as.numeric(gene[x2,3])

boxplot(Y1*1000,Y2*1000,outline = F)
wilcox.test(Y1,Y2,alternative = "greater")


######################### Conservation
conser=read.csv("conservatinScore_ExonSore_Pcg.txt",
                stringsAsFactors=F,sep="\t",skip=0,header =F)
genecode=read.csv("gencode_gene_pro.txt",
                  stringsAsFactors=F,sep="\t",skip=0,header =F)
gene.con=merge(conser,genecode,by.x = "V1",by.y = "V1")
TCGA=read.csv("TCGA_gene_allfreq.txt",
              stringsAsFactors=F,sep="\t",skip=0,header =T)
TCGA2=read.csv("TCGA_gene_Delfreq.txt",
              stringsAsFactors=F,sep="\t",skip=0,header =T)
colnames(TCGA)=c("ENSG","Symbol","Allfreq")
colnames(TCGA2)=c("ENSG","Symbol","Delfreq")
TCGA=merge(TCGA,TCGA2,by.x = "ENSG",by.y = "ENSG")
TCGA.con=merge(TCGA,gene.con,by.x = "Symbol.x",by.y = "V2.y")

TCGA.con$Allfreq=as.numeric(as.character(TCGA.con$Allfreq))*1000
TCGA.con$Delfreq=as.numeric(as.character(TCGA.con$Delfreq))*1000
library(Hmisc)
Wg=as.numeric(cut2(TCGA.con$V2.x, g=10))

boxplot(TCGA.con$V2.x~Wg)
AA=boxplot(TCGA.con$Allfreq~Wg,outline=F)
bb=boxplot(TCGA.con$Delfreq~Wg,outline=F)
library(pastecs)
trend.test(as.numeric(as.character(bb$stats[3,])))
trend.test(as.numeric(as.character(AA$stats[3,])))
######################gene ages analyses
setwd("C:/projects/ceRNA/figure1")
TCGA.f=read.csv("TCGA_gene_allfreq.txt",
                stringsAsFactors=F,sep="\t",skip=0,header =T)
TCGA.d=read.csv("TCGA_gene_Delfreq.txt",
                stringsAsFactors=F,sep="\t",skip=0,header =T)
gene.age=read.csv("human-gene-age-class.txt",
                  stringsAsFactors=F,sep="_",skip=0,header =F)
G.a.f=merge(gene.age,TCGA.f,by.x="V2",by.y = "V2")
G.a.f$age=-ceiling(as.numeric(as.character(G.a.f$V1.x))/3)
G.a.f$Avg.e1=log2(G.a.f$Avg.e*1000)
AA=boxplot(G.a.f$Avg.e1~G.a.f$age,outline=F)
Ga.d=merge(gene.age,TCGA.d,by.x = "V2",by.y = "V2")
Ga.d$age1=-ceiling(as.numeric(as.character(Ga.d$V1.x))/3)
Ga.d$Avg.e1=log2(Ga.d$Avg.e*1000)
BB=boxplot(Ga.d$Avg.e1~Ga.d$age1,outline=F)
trend.test(as.numeric(as.character(BB$stats[3,])))
trend.test(as.numeric(as.character(AA$stats[3,])))



p <- ggplot(df2, aes(variable, value,fill=Tool))
p + geom_boxplot() + labs(title = "CMP")


gene.age=read.csv("human-gene-age-class.txt",
                  stringsAsFactors=F,sep="_",skip=0,header =F)
cor1=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor2=read.csv("core_1.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor3=read.csv("core_2.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
cor1=t(cor1)
cor2=t(cor2)
cor3=t(cor3)
CC1=intersect(cor1,cor2)
CC2=intersect(cor3,CC1)

gene.age$V3=-ceiling(as.numeric(as.character(gene.age$V1))/3)
G.ratio=c()
for(i in -9:-1){
  x=which(gene.age$V3==i)
  aa=length(intersect(gene.age$V2[x],CC2))
  kw=cbind(aa/length(x),1-aa/length(x))
  G.ratio=rbind(G.ratio,kw)
}
barplot(t(G.ratio[,1]))
library(pastecs)
trend.test(G.ratio[,1])

