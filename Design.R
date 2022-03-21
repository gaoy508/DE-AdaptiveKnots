library(plyr)
gcd = function(a,b){
  if (a == 0 | b == 0){
    return (a+b)
  }
  else{
    c = 1
    while(c != 0){
      if (a>b){
        c = a - b
        if(c == 0){
          return(a)
        }
        a = b
        b = c
      }else{
        c = b - a
        if(c == 0){
          return(a)
        }
        b = a
        a = c
      }
    }
  }
}
######################MD2 value######################
MD2value = function(XX){
  # the squared mixture discrepancy
  n = nrow(XX)
  m = ncol(XX)
  q = apply(XX,2,max)
  I1 = matrix(1,n,1)
  tq = I1 %*% q
  if (max(q) >=1){
    X=(XX-0.5)/tq
  }else{
    X=XX
  }
  X0 = abs(X-0.5)
  X1 = 5/3-X0/4-X0^2/4  # the first item of MD
  aa2 = vector()
  if (m!=1){
    ss1 = sum(apply(X1, 1, prod))
    for (i in 1:n){
      X2 = I1 %*% X[i,]
      aa2[i]=sum(apply((15/8-X0/4-abs(X2-0.5)/4-3/4*abs(X-X2)+abs(X-X2)^2/2),1,prod))
    }
  }else{
    ss1=apply(X1,1,sum)
    for (i in 1:n){
      X2=I1 %*% X[i,]
      aa2[i]=sum(15/8-X0/4-abs(X2-0.5)/4-3/4*abs(X-X2)+abs(X-X2)^2/2)
    }
  }
  ss2=sum(aa2)
  y=(19/12)^m-2%*%ss1/n+ss2/(n^2)
  return (y)
}
######################Uniform design######################
pglpmMD1_z = function(N1,s,u){
  # find the uniform design by power good lattice point method
  # where N1 is the number of runs (any integer), s is the number of factors
  a=as.matrix(1:N1,N1,1)
  b=0:(s-1)
  c = matrix(0,N1,s)
  for (i in 1:N1){
    d=a[i]
    for (j in 1:s){
      d=(a[i]*d) %% N1
      c[i,j]=gcd(d,N1)
    }
    if (i %% (u/5) == 0){
      print(paste("*loading...",u,"----",round(i/u,2)*100,"%"))
    }
    
  }
  if(s!=1){
    d = a[apply(c,1,sum)==s]# n is a prime
  }else{
    d=a[apply(c,1,sum)==1]
  }
  I1 = matrix(1,1,s)
  X0=a%*%I1;                # the original design
  MD0=MD2value(X0) # evaluation index
  
  dd = vector()    
  for (i in 1:length(d)){
    dd[1]=1
    dd[2]=d[i]
    for (j in 3:s){
      dd[j]=(dd[j-1]*d[i]) %% N1
    }
    x= (a%*%dd-1)%%N1+1
    MD=MD2value(x)
    if (MD<MD0){
      X=x
      MD0=MD
    }
    if (i %% (length(d)/10) == 0){
      print(paste("loading...",length(d),"----",round(i/length(d),2)*100,"%"))
    }
  }
  y=X
  y = y[1:(N1-1),] #N1-1 run size finally;
  
  qt=rep(u,s)
  q=matrix(rep(qt,N1-1),ncol=length(qt),byrow=TRUE)
  I2 = matrix(1,nrow(y),ncol(y))
  y=floor((y-1)/((N1-1)*I2/q))+1
  
  return(y) 
}  
######################generate knots######################
generateUD1s = function(N,N1,min_knots){
  # the used uniform design in the 1th step in SNTO
  U1=list()
  progress.bar <- create_progress_bar("text")  
  progress.bar$init(N)
  for (n in min_knots:N){
    U1[[n]]=pglpmMD1_z(N1+1,n,N1)
    progress.bar$step() 
  }
  return (U1)
}
generateUD2s = function(N,N2,min_knots){
  # the used uniform design in later steps in SNTO
  U2=list()
  progress.bar <- create_progress_bar("text")  
  progress.bar$init(N)
  for (n in min_knots:N){
    U2[[n]]=pglpmMD1_z(N2+1,n+1,N2)
    progress.bar$step() 
  }
  return (U2)
}
######################Initialization######################
U1 = generateUD1s(N,N1,min_knots)
save.image(file = paste("U1_",N1,"_N",N,"_T",Time,".RData",sep = "", collapse = ""))
U2 = generateUD2s(N,N2,min_knots)
save.image(file = paste("U2_",N2,"_N",N,"_T",Time,".RData",sep = "", collapse = ""))
