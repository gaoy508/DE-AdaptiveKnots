my_bsplineS = function (x, breaks, norder = 4, nderiv = 0, returnMatrix = FALSE) 
{
  x <- as.vector(x)
  n <- length(x)
  tol <- 1e-14
  nbreaks <- length(breaks)
  if (nbreaks < 2) 
    stop("Number of knots less than 2.")
  if (min(diff(breaks)) < 0) 
    stop("Knots are not increasing")
  # if (max(x) > max(breaks) + tol || min(x) < min(breaks) - 
  #     tol) 
  #   stop("Knots do not span the values of X")
  if (x[n] > breaks[nbreaks]) 
    breaks[nbreaks] <- x[n]
  if (x[1] < breaks[1]) 
    breaks[1] <- x[1]
  if (norder > 20) 
    stop("NORDER exceeds 20.")
  if (norder < 1) 
    stop("NORDER less than 1.")
  if (nderiv > 19) 
    stop("NDERIV exceeds 19.")
  if (nderiv < 0) 
    stop("NDERIV is negative.")
  if (nderiv >= norder) 
    stop("NDERIV cannot be as large as order of B-spline.")
  knots <- c(rep(breaks[1], norder - 1), breaks, rep(breaks[nbreaks], 
                                                     norder - 1))
  derivs <- rep(nderiv, n)
  nbasis <- nbreaks + norder - 2
  if (nbasis >= norder) {
    if (nbasis > 1) {
      basismat <- Matrix(splines::spline.des(knots, x, 
                                             norder, derivs)$design)
    }
    else {
      basismat <- as.matrix(splines::spline.des(knots, 
                                                x, norder, derivs)$design)
    }
    if ((!returnMatrix) && (length(dim(basismat)) == 2)) {
      return(as.matrix(basismat))
    }
    return(basismat)
  }
  else {
    stop("NBASIS is less than NORDER.")
  }
}
