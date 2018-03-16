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
  load(paste(cancer[kc],"Expression.rda",sep = ""))
  BRCAmatrix<-assay(data,1,"FPKM")
  allBRCAsample=colnames(BRCAmatrix)
  C_N=data@colData@listData$definition#############sample information tumor vs normal
  Cancer_s=which(C_N!="Solid Tissue Normal")
  Cancer.mat=BRCAmatrix[,Cancer_s]
  Cancer.sam=colnames(Cancer.mat)
  Ca.three=threesam(Cancer.sam)
  allgene=rownames(Cancer.mat)
  P.lnc=c()
  esslnc=t(t(allgene))
  for(i in 1:dim(esslnc)[1]){
    print(c(kc,i))
    x1=which(allgene==esslnc[i,1])
    if(length(x1)>0){
      a=which(as.numeric(Cancer.mat[x1,])>median(as.numeric(Cancer.mat[x1,])))
      b=which(as.numeric(Cancer.mat[x1,])<=median(as.numeric(Cancer.mat[x1,])))
      if(length(a)>0&length(b)>0) {
        Hi.sam=Ca.three[a]
        Lo.sam=Ca.three[b]
        Hi.sam=t(t(Hi.sam))
        Lo.sam=t(t(Lo.sam))
        SS=read.csv(paste(cancer[kc],"_Mysample.txt",sep =""),stringsAsFactors=F,sep="\t",skip=0,header = T)
        Hi.sur=merge(Hi.sam,SS,by.x="V1",by.y="V3")
        Lo.sur=merge(Lo.sam,SS,by.x="V1",by.y="V3")
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
rownames(Pan.P)=esslnc[,1]
colnames(Pan.P)=cancer
write.table(Pan.P,"Pan.P.txt",row.names = FALSE, col.names = T,quote = FALSE,sep="\t")
