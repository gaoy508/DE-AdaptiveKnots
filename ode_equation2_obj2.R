ode_equation2_obj2 <- function(par,tau,Time,lyt,n,lambda){
  
  knots_int = c(min(Time),(max(Time)-min(Time))*tau+min(Time),max(Time))
  knots = unique(knots_int)
  lambda=lambda
  ntimes = length(Time)
  nknots   = length(knots)
  norder   = 4
  nbasis   = length(knots) + norder - 2
  
  v <- par[1]
  gamma <- par[2]
  p <- par[3]
  w <- par[4:(3+nbasis)]
  
  if(is.nan(v)){
    return (10000)
  }else
    
  {####################first part########################
    bsbasis_first = create.bspline.basis(rangeval = c(min(Time),max(Time)),nbasis,norder,knots)
    phi_first=eval.basis(Time, bsbasis_first)
    
    #xhat
    X_first=phi_first%*%w
    
    J=0
    for(t in 1:length(Time)){
      J <- J + sum((lyt[t] - X_first[t])^2)
    }
    
    ####################second part####################
    # Set up quadrature points and weights for Simpson's rule.
    nquad  = 5;
    bsbasis = create.bspline.basis(rangeval = c(min(Time),max(Time)),nbasis,norder,knots)
    eval_delay = sort(c(knots, gamma))
    bsbasis_delay = create.bspline.basis(rangeval = c(min(eval_delay),max(eval_delay)),2 + length(eval_delay),norder,eval_delay)
    quadre_delay = quadset(nquad, breaks=eval_delay)
    quadpts_delay = quadre_delay[,1]
    quadwts= quadre_delay[,2]
    gamma_index = which(quadpts_delay == gamma)[1]
    Q0phi = my_eval.basis(quadpts_delay, bsbasis, 0)
    Q0phi_delay = my_eval.basis((quadpts_delay-gamma)[-(1:gamma_index)], bsbasis, 0)
    
    Q1phi = my_eval.basis(quadpts_delay, bsbasis, 1)
    
    DX1 = Q1phi%*%w
    yhat = exp(DX1)
    
    DX0 = Q0phi%*%w
    DX0_delay = Q0phi_delay%*%w
    
    f = v*(1-exp(DX0_delay)/(1000*p))
    if (length(f)==length(DX1)-1){
      f=c(DX1[1],f)
    }else{
      f=c(DX1[1:(length(DX1)-length(f))],f)
    }
    
    
    DIFF = DX1-f
    PEN = t(DIFF)%*%(DIFF*quadwts)
    
    return(J+lambda*PEN)}
}
