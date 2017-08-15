N<-50; G<-100; Gsignal<-20; Gnoise<-G-Gsignal
theta<-seq(from=0,to=2*pi,length.out=N)
true_cell_positions<-data.frame(dim1=3*cos(theta),dim2=3*sin(theta))
with(true_cell_positions,plot(dim1,dim2))
informative_rows<-as.matrix(true_cell_positions)%*%matrix(rnorm(Gsignal*2),nrow=2)
noise_rows<-matrix(rnorm(Gnoise*N),nrow=Gnoise)
Y<-rbind(t(informative_rows),noise_rows)+rnorm(G)+10
Z<-matrix(rbinom(G*N,1,.8),nrow=G)
Y<-Y*Z
pca_factors<-prcomp(t(Y),center=TRUE,scale=TRUE)$x
plot(pca_factors[,1:2])
vamf_factors<-vamf(Y,5,nrestarts=2,log2trans=FALSE)$factors
with(vamf_factors,plot(dim1,dim2))
