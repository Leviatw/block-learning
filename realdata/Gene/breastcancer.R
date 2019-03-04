#R 3.5.1
#################################################################################################
# Data preprocessing, not need run again, the final data is dat.exp
###############################################################################################
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("affy")
biocLite("hgu95av2.db")
biocLite("hgu133a.db")
biocLite("annotate")
biocLite("genefilter")
biocLite("clusterProfiler")
library(affy)
library("hgu95av2.db")
library("Biobase")
library("annotate")
library(GEOquery)
library(genefilter)
library(hgu133a.db)
library("clusterProfiler")

#After install and load necessary packages, we can start query gene expression data on GEO with accession number, take GSE37418 as an example.

getGEOSuppFiles("GSE11121")
#or directly download the data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse11121
#Then we unzip the raw data files into a directory named "data"
untar("GSE11121/GSE11121_RAW.tar", exdir="gse11121/data")
cels <- list.files("gse11121/data/", pattern = "[gz]")

setwd('gse11121')
sapply(paste("data", cels, sep="/"), gunzip)


#Read in all CEL files and normalize data with RMA method (Highly recommended!!!)
setwd('data/')
AB <- ReadAffy()
ESet <- rma(AB,normalize=T) 
ESet.filtered <- nsFilter(ESet, require.entrez=FALSE, remove.dupEntrez=FALSE)
ESet.filtered$filter.log

#Map probes to gene symbols (convenient for future analysis)
dat <- exprs(ESet.filtered$eset)
prb.list <- row.names(dat)
enID <- select(hgu133a.db, prb.list, c("SYMBOL", "ENTREZID", "GENENAME","UNIPROT"))
enID1 <- enID[!duplicated(enID[,1]), 3]
gs_map <- enID1[!is.na(enID1)]
dat_map <- dat[!is.na(enID1),]


#Remove duplicated observations
dupgene <- names(table(gs_map))[table(gs_map)>1]
indupgene <- names(table(gs_map))[table(gs_map)==1]
expfile.indup <- matrix(0,length(indupgene),ncol(dat_map))
for(i in 1:length(indupgene)){
  expfile.indup[i,] <- as.numeric(dat_map[which(gs_map==indupgene[i],arr.ind=T),])
}
expfile.dup <- matrix(0,length(dupgene),ncol(dat_map))
for(i in 1:length(dupgene)){
  #expfile.dup[i,] <- apply(as.matrix(dat_map[which(gs_map==dupgene[i],arr.ind=T),]),2,mean)
  Mat.dup<-as.matrix(dat_map[which(gs_map==dupgene[i],arr.ind=T),])
  rvs<-apply(Mat.dup,1,sd)
  expfile.dup[i,]<-Mat.dup[which.max(rvs),]
}
exp.output <- rbind(expfile.indup,expfile.dup)
row.names(exp.output) <- c(indupgene,dupgene)


#remove genes not in hgu133plus2ENTREZID

#expo<-exp.output[row.names(exp.output)%in% unique(toTable(hgu133aENTREZID)[,2]),]
expo<-exp.output
dim(expo)
#selected genes with the variance larger than 0.2
rv<-rowVars(expo)
exp.subset<-expo[rv>0.2,]
exp.subset.names<-row.names(exp.subset)
dim(exp.subset)
#4506  200

dat.exp<-t(exp.subset)
colnames(dat.exp)<-exp.subset.names
dim(dat.exp)

saveRDS(dat.exp, "dat.exp")
###############################################################################################

#clustering
dat.exp<-readRDS("dat.exp")
n <- nrow(dat.exp)
d <- ncol(dat.exp)

#Fast estimation of Kendall's tau rank correlation coefficient
#Tau.hat<-cor.fk(dat.exp)
Tau.hat<-cor(dat.exp,method="spearman")
Tau.eg<-eigen(Tau.hat/d)
plot(Tau.eg$values,xlab="Index",ylab="Eigenvalues",pch=20,main="(a) Plot of eigenvalues")   





#using EET to estimate the block size 
k<-sum(abs(Tau.eg$values)>0.1/max(log(n),log(p)))
k
hc<-hclust(dist(Tau.eg$vectors[,1:k]))
km<-cutree(hc,k=k)
table(km)

clust<-lapply(1:k, function(nc) colnames(dat.exp)[km==nc])  
names(clust)<-paste("C",1:k,sep="")


ck2 <- compareCluster(geneCluster = clust, fun = "enrichKEGG")
head(as.data.frame(ck2))

dotplot(ck2)




###############################################################################
##End,2019.1.25
 



