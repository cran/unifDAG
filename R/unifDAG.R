## construct DAG (matrix Q) recursively ###########################
sampleQ <- function(n, K, p.w=1/2) {

  Q <- matrix(0, n, n)
  I <- length(K)
  j <- K[I]
  if(I >= 2) for(i in I:2) { # i = I, I-1, ..., 2
    Ki_1 <- K[i-1L]
    j. <- j + Ki_1
    for(p in (j-K[i]+1):j) {
      repeat{
        for(m in (j+1):j.) { ## FIXME rbinom(1, nn, *)
          Q[m, p] <- rbinom(1, 1, p.w)
        }
        if(sum(Q[(j+1):j., p])>0)
          break
      }
      if(i > 2) {
        for(m in (j.+1):n) { ## FIXME rbinom(1, nn, *)
          Q[m, p] <- rbinom(1, 1, p.w)
        }
      }
    }
    j <- j.
  }
  Q
}

##################################################
## unifDAG
##################################################
## A, B, a, sum, r, t: bigz

unifDAG <- function(n, weighted=FALSE, wFUN=list(runif, min=0.1, max=1)) {
  stopifnot(n>1)
  if (n > 100) stop("Use unifDAG only for n <= 100; for larger n use unifDAG.approx")

  ## step 1
  ## calculate numbers a_{n, k}, b_{n, k} and a_n up to N ##################
  ## is done offline #####################################


  ## step 2
  ## sample an integer between 1 and a_n ##########################
  r <- sampleZ2(.unifDagPreComp$a[n])

  ## step 3
  ## find vector K=c(k_1, ..., k_I) ##############################
  K <- findK.exact(n, r)

  ## step 4
  ## construct DAG (matrix Q) recursively ###########################
  Q <- sampleQ(n, K)

  if(weighted) {
    nrEdge <- sum(Q)
    if(!is.list(wFUN)) {wFUN <- list(wFUN)}
    Q[Q==1] <- do.call(wFUN[[1]], c(nrEdge, wFUN[-1]))
  }

  ## step 5
  ## permute matrix Q and convert to DAG #############################
  perm <- sample.int(n)
  as(Q[perm, perm], "graphNEL")

}


## find vector K=c(k_1, ..., k_I) ##############################
findK.exact <- function(n, r)
{
  K <- rep(0, n) # vector of k_1, ..., k_I
  k <- 1
  while(r>.unifDagPreComp$A[n, k]) {
    r <- r - .unifDagPreComp$A[n, k]
    k <- k+1
  }
  i <- 1
  K[i] <- k
  r <- as.bigz(as.bigq(r, chooseZ(n, k)))+1      #+1: should round to ceil
  m <- n-k
  while(m>0) {
    s <- 1
    t <- (2^k-1)^s * 2^as.bigz(k*(m-s)) * .unifDagPreComp$A[m, s]
    while(r>t) {
      r <- r-t
      s <- s+1
      if(m>=s) {t <- (2^k-1)^as.bigz(s) * 2^as.bigz(k*(m-s)) * .unifDagPreComp$A[m, s]}
      else {t <- r+1}
    }
    if(m>=s) {
      rn.z <- chooseZ(m, s) * (2^k-1)^as.bigz(s) * 2^as.bigz(k*(m-s))
      r.q <-  as.bigq(r, rn.z)
      r <- as.bigz(r.q)  + 1
      nn <- m
      k <- s
      m <- nn-k
      i <- i+1
      K[i] <- k}
    else {
      nn <- m
      k <- min(s, m)
      m <- nn-k
      i <- i+1
      K[i] <- k

    }
  }
  ## I <- i

  K[K!=0]
}

##' @title Sample Uniformly a Large (bigz) Integer
##' @param n a bigz (large) integer
##' @return a random large integer (class \code{"bigz"}) <= n
sampleZ2 <- function(n) {
### numbits <- as.integer(log2(n))+1
  numbits <- as.integer(log2(n-1))+1L
  repeat {
    r.bit <- rbinom(numbits, 1, prob=1/2) # from {0, 1}
    r <- as.bigz(paste0("0b", paste0(r.bit, collapse="")))
    if (r < n)
        return(r + 1)
  }
}


##################################################
## unifDAG.approx
##################################################
unifDAG.approx <- function(n, n.exact = 20, weighted=FALSE,
                           wFUN=list(runif, min=0.1, max=1)) {
  stopifnot(n>1)
  if (n < n.exact) stop("unifDAG.approx: n needs to be at least as big as n.exact!")

  ## step 1&2
  ## calculate numbers a_{n, k}, b_{n, k} and a_n up to N ##################
  ## calculate numbers A_k, B_{s|k} up to N.inf and accuracy #################
  ## is done offline #####################################

  ## step 3
  ## find approx-vector K=c(k_1, ..., k_I) #########################
  K <- findK.approx(n, n.exact)

  ## step 4
  ## construct DAG (matrix Q) recursively ###########################
  Q <- sampleQ(n, K)

  if(weighted) {
    nrEdge <- sum(Q)
    if(!is.list(wFUN)) {wFUN <- list(wFUN)}
    Q[Q==1] <- do.call(wFUN[[1]], c(nrEdge, wFUN[-1]))
  }

  ## step 5
  ## permute matrix Q and convert to DAG #############################
  perm <- sample.int(n)
  as(Q[perm, perm], "graphNEL")

}

## find vector K=c(k_1, ..., k_I) ##############################
findK.approx <- function(n, n.exact)
{
  M <- n
  K1 <- rep(0, n-n.exact)
  i <- 1
  K1[i]  <- sampleZ.cum.vec(.unifDagPreComp$Ak)
  M <- M-K1[i]
  i <- i+1
  while(M>n.exact) {
    K1[i] <- sampleZ.cum.vec(.unifDagPreComp$Bsk[, K1[i-1]])
    M <- M-K1[i]
    i <- i+1
  }
  if(M<n.exact) {
    M <- M+K1[i-1]
    K1[i-1] <- 0
  }
  K1 <- K1[K1!=0]

  K2 <- if(n.exact>=1) {
    ## direct enumeration method with n.exact
    r <- sampleZ2(.unifDagPreComp$a[M])
    findK.exact(M, r)
  }
  else
    0

  K <- c(K1, K2)
  K[K!=0]
}

sampleZ.cum.vec <- function(c) {
  ## c:: bigz-vector; c[i]=numbers of occurance of item i, returns random index, proportional to the numbers in c
  ind <- which(c!=0)
  s <- cumsum(c[ind])
  n <- length(ind)
  r <- sampleZ2(s[n])-1    ## since we want in [1:s[n]]
  ## linear search, since only small c is expected
  i <- 1
  while(s[i]<=r) {
    i <- i+1
  }
  ind[i]
}

##################################################
## Precompute data:
## List unifDagPreComp with elements
## Notation according to: Uniform random generation of large acyclic digraphs
## - A: a_{n, k}
## - B: b_{n, k}
## - a: a_n
## - Ak: A_k
## - Bsk: B_{s|k}
##################################################
if (FALSE) {
    library(gmp)
    setwd("/u/kalischm/research/packages/pcalg/pkg/R")
    source("genRandDAG.R")

    ## Exact --------------------------------
    resExact <- generate.tables(100)
    ##          ---------------
    ## check :
    c1.file <- "/u/kalischm/research/packages/unifDAGs/tables100.RData"
    if(file.exists(c1.file)) {
        load(c1.file)
        stopifnot(identical(resExact[[1]], A),
                  identical(resExact[[2]], B),
                  identical(resExact[[3]], a))
    }

    ## Approx --------------------------------
    resApprox <- approxK(N.inf=100, accuracy=20, A = resExact[["A"]], a = resExact[["a"]])
    ##           -------
    ## check :
    c2.file <- "/u/kalischm/research/packages/unifDAGs/tables_approx100_20.RData"
    if(file.exists(c2.file)) {
        load(c2.file)
        stopifnot(identical(resApprox[[1]], Ak),
                  identical(resApprox[[2]], Bsk))
    }
    ##---- The "precomputed data base" we use ------------------------------

    .unifDagPreComp <- c(resExact, resApprox)
    ##^^^^^^^^^^^^^
    save(.unifDagPreComp,
         file = "/u/kalischm/research/packages/pcalg/pkg/sysdata.rda")
}

## calculate numbers a_{n, k}, b_{n, k} and a_n up to N ##################
## can be done offline ###################################
generate.tables <- function(N, verbose=TRUE)
{
  z0 <- as.bigz(0)
  A <- matrix(z0, N, N) # a_{n, k}
  B <- matrix(z0, N, N) # b_{n, k}
  a <- rep(z0, N)       # a_n

  A[1, 1] <- B[1, 1] <- a[1] <- 1
  for(nn in 2:N) {
    if(verbose) cat(sprintf(" N=%4d / K :", nn))
    for(k in seq_len(nn-1L)) {
      if(verbose) cat(" ", k)
      s <- seq_len(nn-k)
      sum.s <- sum((2^k-1)^as.bigz(s) * 2^as.bigz(k*(nn-k-s)) * A[nn-k, s])
      B[nn, k] <- sum.s
      A[nn, k] <- chooseZ(nn, k) * B[nn, k]
    }
    if(verbose) cat("\n")
    A[nn, nn] <- B[nn, nn] <- 1
    a[nn] <- sum(A[nn, 1:nn])
  }

  ## save(A, B, a, file=paste0(dir, "/tables", N, ".RData"))
  ## cat("\nTables saved in: ", paste0(dir, "/tables", N, ".RData"))
  list(A=A, B=B, a=a)
}

### Construct A_k and B_{s|k} ===================================================

## using  *rational*  arithmetic ("bigq")  "internally" :


approx.Ak <- function(N.inf=100, accuracy=20, A, a) {
  ## Compute  A_k := lim_{n->oo} A_{n, k} / a_n  replacing oo ('Inf') by 'N.inf'

  ## round( 10^acc * A_N / a_N ) :
  Ak <- as.bigz(10^as.bigz(accuracy) * as.vector(A[N.inf,]) / as.vector(a[N.inf]))
  ## typically reducing from 100 to only 10 non-0 ones :
  Ak[Ak != 0]
}

approx.Bsk <- function(Ak) {
  n.k <- length(Ak)
  Bsk <- matrix(as.bigz(0), n.k, n.k)
  for(kk in 1:n.k) {
    ss <- 1:n.k
    ## bug in 'gmp' package: this does nothing !!
    ## Bsk[, kk] <- as.bigz(as.bigq((1-1/(2^kk))^ss) * as.bigq(Ak))
    ##
    ## workaround:
    Bskk <- as.bigz(as.bigq((1-1/(2^kk))^ss) * as.bigq(Ak))
    for(s in ss) Bsk[s,kk] <- Bskk[s]
  }
  Bsk
}


## Need  (A, a) from the exact tables
approxK <- function(N.inf=100, accuracy=20, A, a) {
  Ak <- approx.Ak(N.inf, accuracy, A=A, a=a)
  list(Ak = Ak,
       Bsk= approx.Bsk(Ak))
}


## augment a_{n, k}, b_{n, k} and a_n up form N0 to N ####################
## can be done offline ###################################
## augment.tables <- function(A0, B0, a0, N0, N, dir=getwd(), verbose=FALSE) {
##   A <- as.bigz(matrix(0, N, N))  # a_{n, k}
##   B <- as.bigz(matrix(0, N, N))  # b_{n, k}
##   a <- as.bigz(rep(0, N))       # a_n
##   A[1:N0, 1:N0] <- A0[1:N0, 1:N0]
##   B[1:N0, 1:N0] <- B0[1:N0, 1:N0]
##   a[1:N0] <- a0[1:N0]
##
##   for(nn in ((N0+1):N)) {
##     if(verbose) cat("\n N: ", nn, " K: ")
##     for(k in 1:(nn-1)) {
##       if(verbose) cat(" ", k)
##       sum <- as.bigz(0)
##       for(s in 1:(nn-k)) {
##         sum <- sum + (2^k-1)^as.bigz(s) * 2^as.bigz(k*(nn-k-s)) * A[nn-k, s]
##       }
##       B[nn, k] <-  sum
##       A[nn, k] <- chooseZ(nn, k)*B[nn, k]
##     }
##     A[nn, nn] <- B[nn, nn] <- 1
##     a[nn] <- sum(A[nn, 1:nn])
##   }
##
##   save(A, B, a, file=paste0(dir, "/tables", N, ".RData"))
##   cat("\nAugmented Tables saved in: ", paste0(dir, "/tables", N, ".RData"))
##   list(A, B, a)
## }
