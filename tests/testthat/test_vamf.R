set.seed(100); N<-20; G<-60; Gsignal<-20; Gnoise<-G-Gsignal
theta<-seq(from=0,to=2*pi,length.out=N)
true_cell_positions<-data.frame(dim1=5*cos(theta),dim2=5*sin(theta))
informative_rows<-as.matrix(true_cell_positions)%*%matrix(2*rnorm(Gsignal*2),nrow=2)
noise_rows<-matrix(.5*rnorm(Gnoise*N),nrow=Gnoise)
Y<-rbind(t(informative_rows),noise_rows)+rnorm(G)+10
Z<-matrix(rbinom(G*N,1,.8),nrow=G)
Y<-Y*Z
L<-5
vamf_res<-vamf(Y,L,nrestarts=2,log2trans=FALSE)
rn<-colNorms(vamf_res$factors)
print(rn)
print(sort(rn,decreasing=TRUE))

test_that("Loadings matrix has orthonormal rows", {
  expect_equal(tcrossprod(vamf_res$loadings), diag(L))
})

#test_that("Factors matrix has decreasing row norms", {
#  expect_equal(sort(rn,decreasing=TRUE),rn)
#})