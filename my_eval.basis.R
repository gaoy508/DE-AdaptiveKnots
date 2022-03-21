my_eval.basis = function (evalarg, basisobj, Lfdobj = 0, returnMatrix = FALSE) 
{
  if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
    temp <- basisobj
    basisobj <- evalarg
    evalarg <- temp
  }
  if (is.numeric(evalarg)) {
    if (!is.vector(evalarg)) 
      stop("Argument 'evalarg' is not a vector.")
    Evalarg <- evalarg
  }
  else {
    op <- options(warn = -1)
    Evalarg <- as.numeric(evalarg)
    options(op)
    nNA <- sum(is.na(Evalarg))
    if (nNA > 0) 
      stop("as.numeric(evalarg) contains ", nNA, " NA", 
           c("", "s")[1 + (nNA > 1)], ";  class(evalarg) = ", 
           class(evalarg))
  }
  if (!(inherits(basisobj, "basisfd"))) 
    stop("Second argument is not a basis object.")
  Lfdobj <- int2Lfd(Lfdobj)
  nderiv <- Lfdobj$nderiv
  bwtlist <- Lfdobj$bwtlist
  basismat <- my_getbasismatrix(evalarg, basisobj, nderiv, returnMatrix)
  if (nderiv > 0) {
    nbasis <- dim(basismat)[2]
    oneb <- matrix(1, 1, nbasis)
    nonintwrd <- FALSE
    for (j in 1:nderiv) {
      bfd <- bwtlist[[j]]
      bbasis <- bfd$basis
      if (bbasis$type != "constant" || bfd$coefs != 0) 
        nonintwrd <- TRUE
    }
    if (nonintwrd) {
      for (j in 1:nderiv) {
        bfd <- bwtlist[[j]]
        if (!all(c(bfd$coefs) == 0)) {
          wjarray <- eval.fd(evalarg, bfd, 0, returnMatrix)
          Dbasismat <- my_getbasismatrix(evalarg, basisobj, 
                                      j - 1, returnMatrix)
          basismat <- basismat + (wjarray %*% oneb) * 
            Dbasismat
        }
      }
    }
  }
  if ((!returnMatrix) && (length(dim(basismat)) == 2)) {
    return(as.matrix(basismat))
  }
  else {
    return(basismat)
  }
}
