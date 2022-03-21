my_getbasismatrix = function (evalarg, basisobj, nderiv = 0, returnMatrix = FALSE) 
{
  if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
    temp <- basisobj
    basisobj <- evalarg
    evalarg <- temp
  }
  if (is.null(evalarg)) 
    stop("evalarg required;  is NULL.")
  Evalarg <- evalarg
  op <- options(warn = -1)
  evalarg <- as.numeric(Evalarg)
  options(op)
  nNA <- sum(is.na(evalarg))
  if (nNA > 0) 
    stop("as.numeric(evalarg) contains ", nNA, " NA", c("", 
                                                        "s")[1 + (nNA > 1)], ";  class(evalarg) = ", class(Evalarg))
  if (!(inherits(basisobj, "basisfd"))) 
    stop("Second argument is not a basis object.")
  if (!(length(basisobj$basisvalues) == 0 || is.null(basisobj$basisvalues))) {
    if (!is.vector(basisvalues)) 
      stop("BASISVALUES is not a vector.")
    basisvalues <- basisobj$basisvalues
    nvalues <- length(basisvalues)
    N <- length(evalarg)
    OK <- FALSE
    for (ivalues in 1:nvalues) {
      basisvaluesi <- basisvalues[ivalues]
      if (!is.list(basisvaluesi)) 
        stop("BASISVALUES does not contain lists.")
      argvals <- basisvaluesi[[1]]
      if (!length(basisvaluesi) < nderiv + 2) {
        if (N == length(argvals)) {
          if (all(argvals == evalarg)) {
            basismat <- basisvaluesi[[nderiv + 2]]
            OK <- TRUE
          }
        }
      }
    }
    dimnames(basismat) <- list(NULL, basisobj$names)
    if (OK) {
      if (length(dim(basismat)) == 2) {
        return(as.matrix(basismat))
      }
      return(basismat)
    }
  }
  type <- basisobj$type
  nbasis <- basisobj$nbasis
  params <- basisobj$params
  rangeval <- basisobj$rangeval
  dropind <- basisobj$dropind
  if (type == "bspline") {
    if (length(params) == 0) {
      breaks <- c(rangeval[1], rangeval[2])
    }
    else {
      breaks <- c(rangeval[1], params, rangeval[2])
    }
    norder <- nbasis - length(breaks) + 2
    basismat <- my_bsplineS(evalarg, breaks, norder, nderiv)
  }
  else if (type == "const") {
    basismat <- matrix(1, length(evalarg), 1)
  }
  else if (type == "expon") {
    basismat <- expon(evalarg, params, nderiv)
  }
  else if (type == "fourier") {
    period <- params[1]
    basismat <- fourier(evalarg, nbasis, period, nderiv)
  }
  else if (type == "monom") {
    basismat <- monomial(evalarg, params, nderiv)
  }
  else if (type == "polygonal") {
    basismat <- polyg(evalarg, params)
  }
  else if (type == "power") {
    basismat <- powerbasis(evalarg, params, nderiv)
  }
  else {
    stop("Basis type not recognizable")
  }
  dimnames(basismat) <- list(NULL, basisobj$names)
  if (length(dropind) > 0) 
    basismat <- basismat[, -dropind]
  if (length(evalarg) == 1) {
    basismat = matrix(basismat, 1, length(basismat))
  }
  if (length(dim(basismat)) == 2) {
    return(as.matrix(basismat))
  }
  else {
    return(as.matrix(basismat))
  }
}
