#R 3.5.1
library(rvest)
library(quantmod)
library(clusterCrit)
library(pcaPP)
library(ape)
library(data.table)
library(fGarch)

########################################################################################################################
#Do not need run, download and preprocess the data, notice that the list of NASDAQ 100 stocks  is rebalanced once a year.
#We download the data at December, 2018.

symbols <- read_html("https://www.cnbc.com/nasdaq-100/")

symbols <- symbols %>%
  html_nodes(".text a") %>%
  html_text()
symbols.nasdaq<-symbols
 
# Get stocks value at close from google finance 
options("getSymbols.yahoo.warning"=FALSE)
getSymbols.yahoo(symbols,
                 env = .GlobalEnv,
                 return.class = 'xts',
                 from = "2018-01-02",
                 to = "2018-12-29")

# Create data matrix 

X <- do.call(cbind, sapply(symbols, function(s){
  print(s)
  eval(parse(text = paste("`",s,"`", sep = "")))[,4]
}, simplify = FALSE))
rm(symbols)

colnames(X) <- symbols.nasdaq

lg.X<-diff(log(X))[-1,]

Xr <- sapply(1:(dim(lg.X)[2]),function(i){
  print(symbols.nasdaq[i])
  garchFit(formula = ~garch(1, 1), data = lg.X[,i],
           cond.dist="QMLE", trace = FALSE, delta = 2, skew = 1, shape = 10)@residuals})


NASDAQ <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nasdaq&render=download")

NASDAQ100 <- t(sapply(symbols.nasdaq, function(s){
  NASDAQ[which(NASDAQ[,1] == s),]
},simplify = FALSE))

NASDAQ100 <- rbindlist(NASDAQ100)

#save the stock price data
saveRDS(X, "X") 
#save the residuals for GARCH(1,1) model of log-returns
saveRDS(Xr,"Xr")
#save NASDAQ symbols
saveRDS(symbols.nasdaq,"symbols") 
#save NASDAQ classfication and company detailed information
saveRDS(NASDAQ100,"NASDAQ100")


#######################################################################################################
# Run from here to get consistent results with the paper
######################################################################################################

Xr <- as.matrix(fread("Xr")) #read the stock returns residuals
NASDAQ100 <- fread("NASDAQ100") # read the stock information

#load the stock full names and industry and sector  classification
full.names <- as.vector(NASDAQ100$Name)
sector <- as.vector(NASDAQ100$Sector)
industry <- as.vector(NASDAQ100$industry)
table(sector)
#sector
#        Capital Goods Consumer Non-Durables     Consumer Services           Health Care 
#                    4                     5                    23                    18 
#        Miscellaneous      Public Utilities            Technology        Transportation 
#                    6                     3                    41                     3 
#load the stock names
colnames(Xr) <- readRDS("symbols")
symbols.nasdaq<-colnames(Xr)
n <- nrow(Xr)
d <- ncol(Xr)

 
###show the sector classification
kmt<-as.integer(factor(sector))
lapply(1:max(kmt),function(nc) names.nasdaq[nc==kmt])

#kendall tau correlation matrix
Tau.hat<-cor.fk(Xr)
#show the correlation levels 
range(Tau.hat[lower.tri(Tau.hat)])

Tau.hat['GOOG','GOOGL']
Tau.hat['FOX','FOXA']
Tau.hat['LBTYA','LBTYK']

#compute the eigenvaluese
Tau.eg<-eigen(Tau.hat/d)

plot(Tau.eg$values) 
Tau.eg$values[1:15]

k<-sum(Tau.eg$values>0.08/max(log(n),log(p)))
k


row.names(Tau.eg$vectors)<-symbols.nasdaq
hc<-hclust(dist(Tau.eg$vectors[,1:k]))
km<-cutree(hc,k=k)
#table(cutree(hc,h=sqrt(2/(d-18))))

lapply(1:k,function(nc) symbols.nasdaq[nc==km])



#compare to NASDAQ sector
extCriteria(km,as.integer(factor(sector)),"Precision")$precision
extCriteria(km,as.integer(factor(sector)),"Recall")$recall

#plot the Spearman correlation matrix and compare to the ordered matrix
par(mfrow = c(1,2), mar = c(1,1,1,1))
colfunc <- colorRampPalette(c("darkred","lightyellow","forestgreen"))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(d)) # original Tau.hat

ind<-order(factor(km))
newTau<-Tau.hat[ind,ind] 
image(t(newTau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(d)) # original Tau.hat



#plot the dendrogram
par(mar=c(2,2,1,1))
cols<-c("black","tan1","blue","green3","cyan","magenta", "red","gray")
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,symbols.nasdaq)
     
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=0.8,pch=c(0,1,2,5,6,15,16)[km[ligne]],lab.col=cols[kmt[ligne]],lab.font=1,lab.cex=1))
  }
  return(n)
}

hcd = as.dendrogram(hc)
hcdd<-dendrapply(hcd,colLab)
plot(hcdd)
legend(84,0.8,
legend = names(table(sector)),col = mycols0, 
 bty = "n",  pt.cex = 1.5, cex = 1 , 
text.col = cols, horiz = FALSE, inset = c(0.1, 0.1))



 
