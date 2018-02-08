## Methods to accompany the loggerhead turtle workbook

createProjectionMatrix <- function (surv, fec, d = NULL) {
  ## Create a stage-based population projection matrix
  
  ## Args:
  ##   surv: Stage-based survivorship estimates.
  ##   fec: Stage-based fecundity estimates.
  ##   d: Length of stages. If NULL the matrix reduces to an age-based projection matrix (Leslie matrix).
  
  ## Returns: Projection matrix A in which diagonal elements indicate the
  ## probability of remaining in a given stage; sub-diagonal elements indicate
  ## the probability of surviving and growing into the next stage; and the
  ## first row indicates the contribution to newborns from each stage.
  
  dimA <- length(surv)
  ## expect fec to have same length
  if (length(fec) != dimA)
    stop("expecting equal length survival and fecundity estimates")
  
  if (is.null(d)) {
    d <- rep(1, dimA)
  } else if (length(d) != dimA) {
    stop("expecting equal length survival and stage length vectors")
  }
  
  if(any(surv==1)) {
    warning("replacing survivorship values of 1 with 0.999 to avoid NaNs")
    surv[surv==1] = 0.999
  }
  
  # Calculate P (survival transition) and G (growth transition) values
  P <- ((1 - surv^(d-1))/(1 - surv^d))*surv  # See Crouse 1987, Eq. (1)
  G <- (surv^d*(1 - surv))/(1 - surv^d)  # See Crouse 1987, Eq. (2)
  
  # Create the projection matrix and set the diagonal to P.
  A <- diag(P)
  
  # Set the sub-diagonal to G (discarding the last value which is spurious because
  # the last stage is absorbing).
  diag(A[-1,-dimA]) <- G[1:(dimA-1)]
  
  # Set the first row to the fecundity values, eliding the first entry which we assume
  # corresponds to a non-reproductive stage.
  A[1,2:dimA] <- fec[2:dimA]
  
  return(A)
}

getMPM <- function (surv, fecund, stage.length = NULL) {
  ## Create a stage-based matrix population model
  
  ## Args:
  ##   surv: Stage-based survivorship estimates.
  ##   fecund: Stage-based fecundity estimates.
  ##   stage.length: Length of stages.
  
  ## Returns: List containing the following MPM components:
  ##   A: projection matrix
  ##   eig: object of class eigen containing eigenvectors and eigenvalues for
  ##   the matrix
  ##   lambda: largest eigenvalue (coerced to real number)
  
  A <- createProjectionMatrix(surv, fecund, stage.length)
  eig <- eigen(A, symmetric = FALSE)
  lambda <- Re(eig$values[1])
  
  return(list(A = A, eig = eig, lambda = lambda))
}