# ============================================================
# TSLRF GWAS Model with Nested CV, Variable Importance & Permutation Test
# ============================================================
#Based on
#Sun, J., Wu, Q., Shen, D. et al. TSLRF: Two-Stage Algorithm Based on Least Angle Regression and Random Forest in genome-wide association studies. Sci Rep 9, 18034 (2019). https://doi.org/10.1038/s41598-019-54519-x

# --- Setup ------------------------------------------------------------------
setwd()

# --- Load Required Packages -------------------------------------------------
library(mlr3verse)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3misc)
library(mlr3mbo)
library(mlr3viz)
library(paradox)
library(data.table)
library(readr)
library(dplyr)
library(lars)
# --- Data -------------------------------------------------------------------
set.seed(123)
#options(repos = c(CRAN = "https://cloud.r-project.org/"))


{
  emma.eigen.L <- function(Z,K,complete=TRUE) {
    if ( is.null(Z) ) {
      return(emma.eigen.L.wo.Z(K))
    }
    else {
      return(emma.eigen.L.w.Z(Z,K,complete))
    }
  }
  
  emma.eigen.L.wo.Z <- function(K) {
    eig <- eigen(K,symmetric=TRUE)
    return(list(values=eig$values,vectors=eig$vectors))
  }
  
  emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
    return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
  }
  
  emma.eigen.R <- function(Z,K,X,complete=TRUE) {
    if ( ncol(X) == 0 ) {
      return(emma.eigen.L(Z,K))
    }
    else if ( is.null(Z) ) {
      return(emma.eigen.R.wo.Z(K,X))
    }
    else {
      return(emma.eigen.R.w.Z(Z,K,X,complete))
    }
  }
  
  emma.eigen.R.wo.Z <- function(K, X) {
    n <- nrow(X)
    q <- ncol(X)
    S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
    eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
    stopifnot(!is.complex(eig$values))
    return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
  }
  
  emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
    if ( complete == FALSE ) {
      vids <- colSums(Z) > 0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    n <- nrow(Z)
    t <- ncol(Z)
    q <- ncol(X)
    
    SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
    eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
    if ( is.complex(eig$values) ) {
      eig$values <- Re(eig$values)
      eig$vectors <- Re(eig$vectors) 
    }
    qr.X <- qr.Q(qr(X))
    return(list(values=eig$values[1:(t-q)],
                vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                             complete=TRUE)[,c(1:(t-q),(t+1):n)])) 
  }
  
  emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
    n <- length(xi)
    delta <- exp(logdelta)
    return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) ) 
  }
  
  emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
    delta <- exp(logdelta)
    return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
    
  }
  
  emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
    n <- length(xi)
    delta <- exp(logdelta)
    etasq <- etas*etas
    ldelta <- delta*lambda+1
    return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
  }
  
  emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
    delta <- exp(logdelta)
    etasq <- etas.1*etas.1
    ldelta <- delta*lambda+1
    return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
  }
  
  emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
    nq <- length(etas)
    delta <- exp(logdelta)
    return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
  }
  
  emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-
                   sum(log(delta*lambda+1))) ) 
  }
  
  emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
    nq <- length(etas)
    delta <- exp(logdelta)
    etasq <- etas*etas
    ldelta <- delta*lambda+1
    return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
  }
  
  emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
    t <- t1
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    etasq <- etas.1*etas.1
    ldelta <- delta*lambda+1
    return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
  }
  
  emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL)
  {
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    
    if ( det(crossprod(X,X)) == 0 ) {
      warning("X is singular")
      return (list(ML=0,delta=0,ve=0,vg=0))
    }
    
    if ( is.null(Z) ) {
      if ( is.null(eig.L) ) {
        eig.L <- emma.eigen.L.wo.Z(K)
      }
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.wo.Z(K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,n-q,m) 
      Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
      Xis.1<-matrix(eig.L$values,n,m)
      Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1 
      Etasq <- matrix(etas*etas,n-q,m)
      dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-
                          colSums(Xis.1/Xis))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, 
                       etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
        }
      }
    }
    else {
      if ( is.null(eig.L) ) {
        eig.L <- emma.eigen.L.w.Z(Z,K)
      }
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.w.Z(Z,K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      etas.1 <- etas[1:(t-q)]
      etas.2 <- etas[(t-q+1):(n-q)]
      etas.2.sq <- sum(etas.2*etas.2)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,t-q,m)
      Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
      
      Xis.1<-matrix(eig.L$values,t,m)
      Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1 
      Etasq <- matrix(etas.1*etas.1,t-q,m)
      
      dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-
                          colSums(Xis.1/Xis))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, 
                       etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
    }
    
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    if ( is.null(Z) ) {
      maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n 
    }
    else {
      maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
    }
    maxvg <- maxve*maxdelta
    
    return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  }
  
  emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                         esp=1e-10, eig.L = NULL, eig.R = NULL) {
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    
    if ( det(crossprod(X,X)) == 0 ) {
      warning("X is singular")
      return (list(REML=0,delta=0,ve=0,vg=0))
    }
    
    if ( is.null(Z) ) {
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.wo.Z(K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,n-q,m)
      Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
      Etasq <- matrix(etas*etas,n-q,m)
      
      dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-
                          colSums(Lambdas.1/Lambdas))
      
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, 
                       etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
    }
    else {
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.w.Z(Z,K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      etas.1 <- etas[1:(t-q)]
      etas.2 <- etas[(t-q+1):(n-q)]
      etas.2.sq <- sum(etas.2*etas.2)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1 <- matrix(eig.R$values,t-q,m) 
      Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
      Etasq <- matrix(etas.1*etas.1,t-q,m)
      
      dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-
                          colSums(Lambdas.1/Lambdas))
      
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, 
                       etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
    } 
    
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    
    if ( is.null(Z) ) {
      maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q) 
    }
    else {
      maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
    }
    maxvg <- maxve*maxdelta
    return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  }
  
  emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
    if( is.null(Z) ){
      return(emma.maineffects.B.Zo(K,deltahat.g))
    }
    else{
      return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
    }
  }
  
  emma.maineffects.B.Zo <-function(K,deltahat.g){
    t <- nrow(K)
    stopifnot(ncol(K) == t)
    
    B<-deltahat.g*K+diag(1,t)
    eig<-eigen(B,symmetric=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(mC=C,Q=Q,A=A))
  }
  
  emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    
    n <- nrow(Z) 
    B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
    eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(mC=C,Q=Q,A=A,complete=TRUE))
  }
  
  emma.MLE0.c <- function(Y_c,W_c){
    
    n <- length(Y_c)
    
    stopifnot(nrow(W_c)==n)
    
    M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    etas<-crossprod(M_c,Y_c)
    
    LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
    return(list(ML=LL))
    
  }
  
  emma.REMLE0.c <- function(Y_c,W_c){
    
    n <- length(Y_c)
    
    stopifnot(nrow(W_c)==n)
    
    M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    eig <-eigen(M_c)
    t <-qr(W_c)$rank
    v <-n-t
    U_R <-eig$vector[,1:v]
    etas<-crossprod(U_R,Y_c)
    
    LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
    return(list(REML=LL))
  }
}
# ---------- Paths (adjust) ----------
setwd()

# ---------- Read data ----------
#myGD<-read.table("../imputation_test/new_data/myGD_TKW.hmp.txt",sep = ",",header = T)
myGD<-read.table("../imputation_test/new_data/myGD_PH.hmp.txt",sep = ",",header = T)
#pheno<-read.table("../data/imputation/Qtraits.txt",header = T)
pheno<-read.table("../data/imputation/PH.txt",header = F)
#colnames(pheno)[1]<-"taxa"
pheno<-pheno[,-2]
colnames(pheno)<-c("taxa","PH")

#gwas_data<-merge(pheno[,c(1,6)], myGD, by = "taxa")
gwas_data<-merge(pheno[,c(1,2)], myGD, by = "taxa")
#rownames(gwas_data)<-gwas_data[,1]
#gwas_data<-gwas_data[,-1]
#geno <- as.matrix(gwas_data[,-1])


metadata <- read.csv("../imputation_test/new_data/myGM_PH.csv")
# Ensure taxa column name matches
if (!"taxa" %in% colnames(pheno)) {
  colnames(pheno)[1] <- "taxa"
}

# Save sample order and phenotype vector
sample_ids <- gwas_data$taxa
#y <- as.numeric(gwas_data$TKW)
y <- as.numeric(gwas_data$PH)

# Build genotype matrix: rows = samples, cols = SNPs; preserve SNP names
# genotype columns start at column 3 in merged table (after taxa and trait)
geno_cols_start <- 3L
geno_cols <- colnames(gwas_data)[geno_cols_start:ncol(gwas_data)]

geno_dt <- gwas_data[, geno_cols]
# Recode each column to numeric while preserving column names
rownames(geno_dt) <- sample_ids
colnames(geno_dt) <- geno_cols

# Quick sanity checks
if (length(y) != nrow(geno_dt)) stop("phenotype length doesn't match genotype rows")

# ---------- Kinship calculation ----------
# 
M <- as.matrix(geno_dt) # samples x markers

# Center columns
M_centered <- scale(M, center = TRUE, scale = FALSE)
K_raw <- tcrossprod(M_centered) # samples x samples
cc <- mean(diag(K_raw))
K <- K_raw / cc # genetic relationship matrix

# ---------- Adjust for Population Structure using EMMA  ----------
#W and Z are the corresponding designed matrices for α and γ; u
W <- matrix(1, nrow = nrow(M), ncol = 1)

remle_res <- emma.REMLE(y = y, X = W, K = K, Z = NULL, ngrids = 100, llim = -10, ulim = 10, esp = 1e-10)

# fix call to main effects
remle_B <- emma.maineffects.B(K = K, deltahat.g = remle_res$delta)
C2 <- remle_B$mC

# ---------- Transform phenotype and genotype while keeping names ----------
y_trans <- as.numeric(C2 %*% y)

# Transform genotype: we want the result to be samples x markers with same colnames
geno_trans <- as.matrix(C2 %*% M_centered)  # samples x markers
colnames(geno_trans) <- colnames(M)
rownames(geno_trans) <- rownames(M)

# Build design matrix for LARS: intercept + genotype predictors
x_mat <- cbind(Intercept = as.numeric(C2 %*% rep(1, nrow(M))), geno_trans)
# ensure colnames are present
if (is.null(colnames(x_mat))) stop("x_mat should have column names")

# ---------- LARS: fit and map selected predictors back to SNP names ----------
# lars expects predictors (n x p) and response (n)
LAR <- lars(x = as.matrix(x_mat), y = y_trans, type = "lar", trace = FALSE, normalize = TRUE,
            intercept = TRUE, use.Gram = FALSE, max.steps = ncol(x_mat) - 1)

# Get coefficients at final step and selected variable names
final_beta <- LAR$beta[nrow(LAR$beta), , drop = TRUE]
selected_cols_logical <- final_beta != 0
selected_names_all <- colnames(x_mat)[selected_cols_logical]
# Remove intercept if present
selected_snp_names <- setdiff(selected_names_all, "Intercept")

# rf_input: predictors selected by LARS 
rf_input <- x_mat[, selected_names_all, drop = FALSE]

# ---------- Build mlr3 task ----------
combined_df <- as.data.frame(rf_input)
combined_df$y <- as.numeric(y_trans)

task_gwas_rf <- TaskRegr$new(id = "gwas_rf", backend = combined_df, target = "y")

# ---------- Define learner with tuning space ----------
learner_rf = lrn("regr.ranger", num.threads = 1, importance = "impurity",
                 num.trees = to_tune(300, 1500),#1000 # Number of trees
                 mtry = to_tune(1, 300),#100 # Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables
                 min.node.size = to_tune(1,500),# Minimal node size to split at. Default 1 for classification, 5 for regression
                 max.depth = to_tune(1,100))# Maximal tree depth. A value of NULL or 0 (the default) corresponds to unlimited depth


# Tuner (MBO) and terminator
tuner_rf <- tnr("mbo")
terminator_eval <- trm("evals", n_evals = 10)
outer_rs <- rsmp("cv", folds = 5)

# AutoTuner
at <- AutoTuner$new(
  learner = learner_rf,
  resampling = rsmp("cv", folds = 5),#inner resampling
  measure = msr("regr.rmse"),
  terminator = terminator_eval,
  tuner = tuner_rf
)

# --- Run Nested Cross-Validation --------------------------------------------

rr <- resample(task_gwas_rf, at, outer_rs, store_models = TRUE)

# ============================================================
# Model Performance Evaluation
# ============================================================
# Define performance measures
measures <- msrs(c("regr.rmse", "regr.mae", "regr.rsq", "time_train", "time_predict"))

# Aggregate results
cat("Aggregated performance:\n")
print(rr$aggregate(measures = measures))

# Extract best hyperparameters from best outer fold
score_dt <- rr$score(msr("regr.mse"))
best_fold_index <- which.min(score_dt$regr.mse)
cat("Best outer fold index:", best_fold_index, "\n")

best_resampling_learner <- rr$learners[[best_fold_index]]
best_tune_result <- best_resampling_learner$tuning_result
cat("Best-hyperparameters from best outer fold:\n")
print(best_tune_result$x)

# Use the final trained model from the best fold
rf_model <- rr$learners[[best_fold_index]]$model
#learner_rf$param_set$values <- best_hyperpars[[1]]
learner_rf$param_set$values <- c(best_tune_result$x_domain[[1]], list(importance = "impurity"))

# ============================================================
# Variable Importance (Average over Multiple Runs)
# ============================================================

n_reps <- 100                # how many repeated 5-fold CV runs
n_workers_reps <- 10             # how many parallel workers for the 100 reps 

feature_names <- task_gwas_rf$feature_names
n_features <- length(feature_names)

run_one_rep <- function(rep_idx) {
  cat("Starting repetition", rep_idx, "\n")
  set.seed(1000 + rep_idx)
  
  learner_clone <- learner_rf$clone(deep = TRUE)
  resampling_clone <- outer_rs$clone(deep = TRUE)
  
  rr_local <- resample(
    task = task_gwas_rf,
    learner = learner_clone,
    resampling = resampling_clone,
    store_models = TRUE
  )
  
  per_fold_importances <- matrix(
    0, nrow = n_features, ncol = length(rr_local$learners),
    dimnames = list(feature_names, paste0("fold", seq_along(rr_local$learners)))
  )
  
  for (i in seq_along(rr_local$learners)) {
    learner_i <- rr_local$learners[[i]]
    imp_vec <- learner_i$importance()
    
    tmp <- numeric(n_features)
    names(tmp) <- feature_names
    
    if (!is.null(imp_vec) && length(imp_vec) > 0)
      tmp[names(imp_vec)] <- as.numeric(imp_vec)
    
    per_fold_importances[, i] <- tmp
  }
  
  mean_imp_rep <- rowMeans(per_fold_importances, na.rm = TRUE)
  cat("Finished repetition", rep_idx, "\n")
  
  return(list(mean_importance = mean_imp_rep))
}

plan(multisession, workers = n_workers_reps)  # parallel across reps
rep_results <- future_lapply(seq_len(n_reps), function(r) {
  run_one_rep(r)
}, future.seed = TRUE)
plan(sequential)

imp_all_runs <- sapply(rep_results, function(res) res$mean_importance)


imp_mean <- rowMeans(imp_all_runs, na.rm = TRUE)

imp_summary <- data.table(
  marker = feature_names,
  mean_importance = imp_mean)
                       [order(-mean_importance)]

# ============================================================
# Permutation Test for Empirical Threshold
# ============================================================
B <- 1000 #number permutations
n_workers <- 10

compute_mean_imp_from_rr <- function(rr_local, feature_names) {
  n_iter_local <- length(rr_local$learners)
  imp_mat_local <- matrix(0, nrow = length(feature_names), ncol = n_iter_local,
                          dimnames = list(feature_names, paste0("iter", seq_len(n_iter_local))))
  for (i in seq_len(n_iter_local)) {
    imp_vec <- rr_local$learners[[i]]$importance()
    if (!is.null(imp_vec) && length(imp_vec) > 0) {
      tmp <- numeric(length(feature_names))
      names(tmp) <- feature_names
      tmp[names(imp_vec)] <- as.numeric(imp_vec)
      imp_mat_local[, i] <- tmp
    }
  }
  rowMeans(imp_mat_local, na.rm = TRUE)
}

run_one_perm <- function(seed = NULL, resampling_clone) {
  if (!is.null(seed)) set.seed(seed)
  
  data_perm <- as.data.table(task_gwas_rf$data())
  target_name <- task_gwas_rf$target_names
  data_perm[[target_name]] <- sample(data_perm[[target_name]])
  
  task_perm <- TaskRegr$new(
    id = paste0("perm_task_", sample.int(1e8, 1)),
    backend = data_perm,
    target = target_name
  )
  
  learner_clone <- learner_rf$clone(deep = TRUE)
  rr_perm <- resample(
    task = task_perm,
    learner = learner_clone,
    resampling = resampling_clone,
    store_models = TRUE
  )
  
  mean_imp_perm <- compute_mean_imp_from_rr(rr_perm, feature_names)
  max(mean_imp_perm, na.rm = TRUE)
}

set.seed(123)
seeds <- sample.int(.Machine$integer.max, B)
resampling_clone <- outer_rs$clone(deep = TRUE)

cat("Starting permutation threshold with", B, "permutations on", n_workers, "workers\n")

future::plan("multisession", workers = n_workers)
with_progress({
  perm_max_imp <- future_lapply(seq_len(B), function(b) {
    run_one_perm(seed = seeds[b], resampling_clone = resampling_clone)
  }, future.seed = TRUE)
})
future::plan("sequential")

perm_max_imp <- vapply(perm_max_imp, as.numeric, numeric(1))
emp_threshold_perm <- as.numeric(quantile(perm_max_imp, probs = alpha, na.rm = TRUE))

cat(sprintf("permutation threshold (%.1f%% percentile): %g\n",0.95, emp_threshold_perm))



# ============================================================
# Merge with Marker Metadata and Save Results
# ============================================================
#save.image("TSLRF_PH_mlr3.RData")
myGM<-read.csv("myGM_TKW.csv")
avg_importance2<-merge(avg_importance,myGM, by="SNP")
write.csv(imp_summary, file = "TSLRF_avg_importance_TKW.csv",row.names = F)
write.table(emp_threshold_perm, "emp_threshold_TKW_TSLRF_.txt")
