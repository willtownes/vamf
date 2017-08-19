#' @useDynLib vamf, .registration = TRUE
#' @import Rcpp
#' @import methods
NULL

#' Remove zero rows and columns.
#' 
#' Removes all rows and columns that have only zeros from a sparse Matrix.
#' 
#' @param Y A sparse Matrix
#' @return The same Matrix without zero rows and columns
rm_zero_rowcol<-function(Y){
  #remove all rows and columns containing all zeros
  Y<-Y[Matrix::rowSums(Y>0)>0,] #remove rows with zeros all the way across
  Y<-Y[,Matrix::colSums(Y>0)>0]
  Y
}

STMOD<-stanmodels[["vamf"]]

#' L2 norm of an object.
#' 
#' L2 norm of an object.
#' 
#' @param v An vector or array-like object
#' @return The L2 (Euclidean) norm of v
norm<-function(v){sqrt(sum(v^2))}

#' L2 norms of columns of a matrix.
#' 
#' L2 norms of columns of a matrix.
#' @param x a matrix
#' @return a vector containing the L2 norms of the columns of x.
#' @export
colNorms<-function(x){
  apply(x,2,norm)
}

#' Effective dimension.
#' 
#' Finds the effective dimension of a latent factor matrix by first computing L2 norms of each column.
#' 
#' Columns with L2 norm greater than the maximum norm times 'thresh' are nonzero.
#' Columns with L2 norm less than the cutoff (above) are 'zero'.
#' @param x a matrix whose columns represent coordinates in a latent space
#' @param thresh a threshold close to zero
#' @return a count of the number of nonzero columns
effective_dimension<-function(x,thresh=.05){
  l2norms<-colNorms(x)
  sum( l2norms > max(l2norms)*thresh )
}

#' Initialize sufficient statistics.
#' 
#' Calculate and cache in a list various sufficient statistics of the data used by the VAMF algorithm.
#' 
#' @param Y a data matrix, typically non-negative values such as counts or TPM
#' @param L desired latent dimension (upper bound)
#' @param log2trans should nonzero values of Y be converted to log-scale?
#' @param pseudocount small offset for log transformation. Transformation is g(y)=log2(pseudocount+y) for all nonzero y values
#' @param b1_range rough guess of typical slopes for the dropout mechanism
#' @return a list of sufficient statistics used by other vamf functions
init_ss<-function(Y,L,log2trans=TRUE,pseudocount=0.0,b1_range=c(0.3,.7)){
  Y<-Matrix::Matrix(Y,sparse=TRUE)
  Y<-rm_zero_rowcol(Y) #remove rows with all zeros
  Z<-Y>0
  #convert sparse matrix to stan-friendly long format
  Ys<-Matrix::summary(Y) #converts to triplet format
  if(log2trans) Ys[,3]<-log2(pseudocount+Ys[,3])
  #stan data variables & hyperparameters
  ss<-list(gg=Ys[,1],nn=Ys[,2],y=Ys[,3],nvals=nrow(Ys),L=L,N=ncol(Y),G=nrow(Y))
  ss$Z<-t(matrix(as.integer(as.matrix(Z)),nrow=ss$G)) #dense integer matrix for stan
  #note ss$Z is transpose of Z
  ss$ymn<-stats::median(ss$y) #automatically calculated in stan program
  Yctr<-Y-ss$ymn*Z
  #a<-colSums(Yctr)/colSums(Z)
  #Yctr<-t(t(Yctr)-a*t(Z)) #take advantage of recycling
  w<-Matrix::rowSums(Yctr)/Matrix::rowSums(Z)
  #ss$sa<-stats::mad(a)
  ss$sw<-stats::mad(w)
  ### this block inefficient, improve later!
  Yctr<-as.matrix(Yctr)-w
  #sd_cols<-apply(Yctr,2,function(x){stats::mad(x[x>0])})
  #ss$su<-mean(sd_cols[!is.na(sd_cols)])
  #rough estimate for variation in the row factors
  sd_rows<-apply(Yctr,1,function(x){stats::mad(x[x>0])})
  ss$sv<-mean(sd_rows[!is.na(sd_rows)])
  ### end inefficient block
  ss$Q<-Matrix::colMeans(Z) #detection rates
  ss$b1_mn<-mean(b1_range)
  ss$b1_sd<-diff(b1_range)/4 #ie, range is +/- 2 standard devs from mean
  ss
}

#' Variational Bayes wrapper function.
#' 
#' Convenience function for parallelizing calls to VB.
#' 
#' @param svmult Scalar multiplier to increase or decrease the sigma_v scale hyperparameter.
#' @param stmod An object of class \code{\linkS4class{stanmodel}}.
#' @param ss A list of sufficient statistics from \code{\link{init_ss}}.
#' @param resnames A string vector of parameter names to be included in output.
#' @param output_samples How many samples from the approximate posterior to estimate the posterior mean of each parameter.
#' @return A list with components stan_vb_fit (the result from Stan) and logtxt (the log of Stan output).
vb_wrap<-function(svmult,stmod,ss,resnames,output_samples){
  ss$sv_rate<- 1/(svmult*ss$sv) #empirical bayes hyperparam set
  #since Gamma shape is 2, the mode is svmult*ss$sv
  logtxt<-utils::capture.output(stan_vb_fit<-rstan::vb(stmod,data=ss,pars=resnames,output_samples=output_samples),type="output",split="true")
  mget(c("stan_vb_fit","logtxt"))
}

#' Varying-Censoring Aware Matrix Factorization via Stan.
#' 
#' Runs VAMF with many random initializations in parallel. Parallelism will only work if the number of cores is set as a global option. For example \code{options(mc.cores=4)}. Discards any random restarts that resulted in Stan errors, usually due to numeric convergence failures.
#' 
#' @param ss A list of sufficient statistics from \code{\link{init_ss}}.
#' @param svmult Vector of multipliers to increase or decrease the sigma_v scale hyperparameter. Length of svmult determines number of random restarts to run.
#' @param output_samples How many samples from the approximate posterior to estimate the posterior mean of each parameter.
#' @return List of results from each random initialization excluding numerical convergence failures. Each element of the list is an object with structure described in \code{\link{vb_wrap}}.
vamf_stan<-function(ss, svmult=rep.int(1.0,4), output_samples=100){
  resnames<-c("U","w","sy","y0","V","sv","b0","b1")
  #variational bayes
  res<-parallel::mclapply(svmult,vb_wrap,STMOD,ss,resnames,output_samples)
  #get rid of model runs that resulted in errors
  return(Filter(function(x){class(x)!="try-error"},res))
}

#' Extract evidence lower bound (ELBO).
#' 
#' Takes the log text from running variational Bayes in Stan and extracts out the value of the evidence lower bound (ELBO) at the final iteration. Higher ELBO is better.
#' @param logtxt A list of strings. Each element is a single line of the output log text of a Stan run.
#' @return The scalar value of the ELBO.
extract_elbo<-function(logtxt){
  elbo<- -Inf
  x<-grep("ELBO CONVERGED",logtxt,value=TRUE,fixed=TRUE)
  if(length(x)!=1){
    if(length(x)==0){
      warning("Failed to converge!")
    } else {
      warning("Multiple convergence points found!")
    }
    return(elbo)
  }
  x<-unlist(strsplit(x,"\\s+",perl=TRUE))
  elbo<-as.numeric(x[3])
  #parse out scientific notation eg -3e+04
  #x<-grep("\\de\\+\\d+",x,perl=TRUE,value=TRUE) 
  #elbo<-as.numeric(x)
  return(elbo)
}

#' Extract latent factors from a Stan model fit.
#' 
#' Extracts out the posterior means of the latent factors and all other parameter values from a Stan model fit.
#' @param stan_vb_fit Stan variational Bayes object obtained from the output of \code{\link{vb_wrap}}.
#' @param ss List of sufficient statistics from \code{\link{init_ss}}.
#' @return List of posterior means of all model parameters.
extract_factors<-function(stan_vb_fit,ss){
  #stan_vb_fit<-vamf_stan_fit$stan_vb_fit
  resnames<-stan_vb_fit@sim$pars_oi
  resnames<-resnames[resnames != "lp__"]
  vars_tmp<-rstan::summary(stan_vb_fit)$summary[,"mean"]
  varnames<-names(vars_tmp)
  vars<-lapply(resnames,function(x){vars_tmp[gdata::startsWith(varnames,x)]})
  #vars<-lapply(resnames,function(x){get_posterior_mean(stan_vb_fit,x)})
  names(vars)<-resnames
  #stan output indexes first by rows, then cols (inconvenient for recycling)
  vars$U<-t(matrix(vars$U,ncol=ss$L)) #output dim: LxN
  vars$V<-matrix(vars$V,nrow=ss$L) #output dim: LxG
  vars
}

#' Orthogonalize and normalize latent factors.
#' 
#' Rotate and transform the latent factors to an equivalent representation that facilitates interpretability by using an orthonormal basis in the loading matrix.
#' @param vars List of posterior means of model parameters from \code{\link{extract_factors}}.
#' @return The same list but with additional components 'factors' and 'loadings' whose interpretation is analogous to Principal Components Analysis.
ortho_vamf<-function(vars){
  v<-vars$V #LxG
  u<-vars$U #LxN, recycles the sv vector
  svd_v<-svd(v)
  A<-svd_v$u #LxL
  D<-if(length(svd_v$d)>1) diag(svd_v$d) else svd_v$d #LxL
  Q<-svd_v$v #GxL
  loadings<-t(Q) #LxG
  factors<-data.frame(crossprod(u,A%*%D))
  colnames(factors)<-paste0("dim",1:ncol(factors))
  c(vars,mget(c("factors","loadings")))
}

ortho_extract<-function(stfit,ss){
  ortho_vamf(extract_factors(stfit$stan_vb_fit,ss))
}

#' Varying-Censoring Aware Matrix Factorization.
#' 
#' VAMF is a probabilistic dimension reduction method intended for single cell RNA-Seq datasets.
#' @param Y Sparse Matrix of gene expression measurements, with G genetic features (genes) in the rows and N samples (typically, individual cells) in the columns.
#' @param L Upper bound on the dimensionality of the latent space to be learned. Automatic relevance determination is used to shrink away unnecessary dimensions.
#' @param nrestarts Number of independent random initializations of the algorithm. Can be parallelized by setting e.g. \code{options(mc.cores=4)}.
#' @param log2trans Should the data be log transformed prior to analysis? Set to FALSE if the data have already been log transformed.
#' @param pseudocount Optional small offset to be added to data before log transformation.
#' @param output_samples Number of samples from approximate posterior used to estimate the posterior means of all parameters.
#' @param save_restarts If multiple initializations are used, set this to TRUE if you want to return the list of all results. Set to FALSE to choose only the best result based on the highest evidence lower bound (ELBO).
#' @param svmult Scalar or vector of multipliers to increase or decrease the sigma_v scale hyperparameter.
#' @return Named list of posterior means for model parameters. The 'factors' and 'loadings' are analogous to PCA. Cell positions in latent space can be plotted by using the 'factors' matrix. If save_restarts is set to TRUE, returns a list of lists, each from a separate VAMF run.
#' \describe{
#'   \item{factors}{NxL matrix whose columns are analogous to principle components. The L2 norm of each column indicates the significance level of the component. The transposed orientation is to facilitate plotting.}
#'   \item{loadings}{LxG matrix whose rows are analogous to principle component loadings. The rows are orthonormal.}
#'   \item{effdim}{Effective dimensionality of the latent space. Computed by L2 norms of the 'factors' matrix}
#'   \item{elbo}{Evidence lower bound, the objective function for variational inference. See \href{http://mc-stan.org/users/documentation/index.html}{Stan user manual}}
#'   \item{b0}{Censoring mechanism random intercepts for each cell (vector of length N)}
#'   \item{b1}{Censoring mechanism random slopes for each cell (vector of length N)}
#'   \item{U}{Raw version of the factors matrix (without rotations and scaling) dimension LxN. Note the factors matrix has a transposed orientation relative to U}
#'   \item{V}{Raw version of loadings matrix (without rotations and scaling) dimension LxG}
#'   \item{w}{Vector of length G with row-specific random intercepts}
#'   \item{y0}{Global intercept (scalar)}
#'   \item{sy}{Standard deviation of global noise (scalar)}
#'   \item{sv}{Standard deviations of each latent dimension (vector of length L), interpretable only with U and V, not interpretable with 'factors' and 'loadings'}
#'   \item{svmult}{Same as the input parameter(s)}
#' }
#' @export
#' @examples
#' set.seed(100); N<-20; G<-60; Gsignal<-20; Gnoise<-G-Gsignal
#' theta<-seq(from=0,to=2*pi,length.out=N)
#' true_cell_positions<-data.frame(dim1=5*cos(theta),dim2=5*sin(theta))
#' with(true_cell_positions,plot(dim1,dim2))
#' informative_rows<-as.matrix(true_cell_positions)%*%matrix(2*rnorm(Gsignal*2),nrow=2)
#' noise_rows<-matrix(.5*rnorm(Gnoise*N),nrow=Gnoise)
#' Y<-rbind(t(informative_rows),noise_rows)+rnorm(G)+10
#' Z<-matrix(rbinom(G*N,1,.8),nrow=G)
#' Y<-Y*Z
#' pca_factors<-prcomp(t(Y),center=TRUE,scale=TRUE)$x
#' plot(pca_factors[,1:2])
#' vamf_factors<-vamf(Y,5,nrestarts=2,log2trans=FALSE)$factors
#' with(vamf_factors,plot(dim1,dim2))
vamf<-function(Y, L, nrestarts=4, log2trans=TRUE, pseudocount=0.0, output_samples=100, save_restarts=FALSE,svmult=1){
  #convenience wrapper for running instances of selection factorization
  ss<-init_ss(Y,L,log2trans=log2trans,pseudocount=pseudocount)
  svmult_vals<-rep_len(svmult,nrestarts)
  stan_fits<-vamf_stan(ss, svmult_vals, output_samples=output_samples)
  elbos<-vapply(stan_fits, function(x){extract_elbo(x$logtxt)}, 1.0)
  if(save_restarts){
    factor_list<-lapply(stan_fits, ortho_extract, ss)
  } else { #keep only max elbo values
    good<- elbos==max(elbos)
    factor_list<-lapply(stan_fits[good], ortho_extract, ss)
    svmult_vals<-svmult_vals[good]
    elbos<-elbos[good]
  }
  eff_dim_table<-rep.int(0,L)
  eff_dim_idx<-rep.int(NA,length(elbos))
  for(i in seq_along(elbos)){
    factor_list[[i]]$elbo<-elbos[i]
    factor_list[[i]]$svmult<-svmult_vals[i]
    L_eff<-effective_dimension(factor_list[[i]]$factors)
    factor_list[[i]]$effdim<-L_eff
    eff_dim_table[L_eff]<-eff_dim_table[L_eff]+1
    eff_dim_idx[i]<-L_eff
  }
  if(save_restarts){
    return(factor_list)
  } else if(length(factor_list)==1){
    return(factor_list[[1]])
  } else { #if ELBOs are tied, choose result with most common effective dimension
    L_eff_popular<-which.max(eff_dim_table) #most popular effective dimension
    return(factor_list[eff_dim_idx==L_eff_popular][[1]])
  }
}
