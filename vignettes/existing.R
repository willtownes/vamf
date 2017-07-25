# Functions for running existing methods such as tSNE, PCA, MDS
# slight modification of https://github.com/willtownes/vamf-paper/blob/master/algs/existing.R
library(Matrix)

rm_zero_rowcol<-function(Y){
  #remove all rows and columns containing all zeros
  Y<-Y[rowSums(Y>0)>0,] #remove rows with zeros all the way across
  Y<-Y[,colSums(Y>0)>0]
  Y
}

norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}

pca<-function(Y,L=2,center=TRUE,scale=TRUE){
  Y<-rm_zero_rowcol(Y)
  #if(scale) Y<-scale(Y)
  factors<-prcomp(as.matrix(t(Y)),center=center,scale=scale)$x
  factors<-factors[,1:L]
  colnames(factors)<-paste0("dim",1:L)
  as.data.frame(factors)
}

mds<-function(Y,L=2,metric=TRUE,distance="euclidean",scale=TRUE){
  #Multidimensional Scaling
  #Y is data
  #L is desired latent dimension
  #metric=TRUE means use cmdscale(). metric=FALSE means use isoMDS()
  #see http://www.statmethods.net/advstats/mds.html
  #distance is choice of distance function passed to dist()
  Y<-rm_zero_rowcol(Y)
  Yt<-scale(t(Y),scale=scale)
  d <- dist(as.matrix(Yt),distance) # euclidean distances between the cols
  if(metric){
    fit<-cmdscale(d,k=L)
  } else {
    fit<-MASS::isoMDS(d,k=L)$points
  }
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}

tsne<-function(Y,L=2,center=TRUE,scale=TRUE,...){
  Yt<-t(rm_zero_rowcol(Y))
  if(center || scale){
    Yt<-scale(Yt,center=center,scale=scale)
  }
  Yt<-as.matrix(Yt)
  fit<-Rtsne::Rtsne(Yt,dims=L,...)$Y
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}
