SNTO_knots = function(lambda,N_list,N,Time,U1,U2,N1,N2,delta){
  #transform the minimal interval to [0,1]
  ell = 1
  ell_d1 = 1/(max(Time)-min(Time))
  
  #Initialization
  CMINn=vector();
  Tn=list();
  
  Rn=list();
  #CMIN[w,n] min{object function} with fixed n,step w
  total_step = 20
  CMIN=matrix(0,total_step,200)
  T=list();
  
  for (n in N_list){
    T[[n]] = matrix(NA,total_step,n)
    #####Initialization: a,b
    a = matrix(0,n,1)
    b = matrix(1,n,1)
    
    c = (b-a)/2
    
    #Initialization Uniform Design: U = (U_{ki})
    U=U1[[n]];
    
    #N1: nrow-U
    qt=matrix(N1,1,n);
    l=matrix(rep(qt,N1),ncol=ncol(qt),byrow=TRUE);  # number of levels
    
    #ncol(C)=n,where the first n columns have u levels.
    C=matrix(0,N1,n); 
    C[,1:n]=U[,1:n]/l;  # transform to [0,1]^(n+1) #u_{ki}'
    
    #ell_d1:l
    B=1-(n-1)*ell_d1;
    
    #####Proposition 1: transform between x_k(A) and c_k(C)#####
    A=matrix(0,N1,n);
    A1=matrix(0,N1,n);
    
    for (i in 1:n){
      ones = matrix(1,N1,1)
      A1[,i]=B*ones;
      for (j in i:n){
        A1[,i]=A1[,i]*(C[,j]^(1/j));
      }
      A[,i]=(i-1)*ell_d1*ones+A1[,i];
    } 
    
    #####Calculate the object function value
    C_value=vector();
    
    for (i in 1:nrow(A)){
      knots_int = c(min(Time),(max(Time)-min(Time))*A[i,]+min(Time),max(Time))
      knots = unique(knots_int)
      nknots = length(knots)
      norder = 4
      nbasis = length(knots) + norder - 2
      bsbasis = create.bspline.basis(rangeval = c(min(Time),max(Time)),nbasis,norder,knots)
      phi=my_eval.basis(Time, bsbasis)
      Mmat = ginv(t(phi)%*%phi)%*%t(phi)
      chat = Mmat%*%lyt
      ini=c(v=v0,gamma=gamma0,p=p0,w=chat)
      par_hat=nlminb(ini, objective = ode_equation2_obj2,tau=A[i,],Time=Time, lyt = lyt,n=n,control = list(rel.tol = 1e-3),lambda=lambda)
      
      vh = par_hat$par[1]
      gammah = par_hat$par[2]
      ph = par_hat$par[3]
      wh = par_hat$par[4:(3+nbasis)]
      par_hat=c(vh,gammah,ph,wh)
      C_value[i]=ode_equation2_obj2(par=par_hat,tau=A[i,1:n],Time=Time,lyt=lyt,n=n,lambda=lambda);
    }
    
    for (l in 1:length(C_value)){
      if (is.nan(C_value[l])|is.infinite(C_value[l])){
        C_value[l] = 1000000
      }
    }
    Cmin=min(C_value);
    
    t_opt=A[which(C_value==Cmin),] 
    
    #step 1
    CMIN[1,n]=Cmin;
    I_ones = matrix(1,1,n)
    T[[n]][1,] = (max(Time)-min(Time))*t_opt[1:n]+min(Time);
    
    #transform from A_(n) back to [a,b]^(n)
    c_opt=vector();
    c_opt[n]=((t_opt[n]-(n-1)*ell_d1)/B)^n;
    for (i in 1:(n-1)){
      c_opt[i]=((t_opt[i]-(i-1)*ell_d1)/(t_opt[i+1]-i*ell_d1))^i;
    }
    
    
    ########################  MSNTO second and later for run size = N2 #######################################
    U=U2[[n]];
    qt=matrix(N2,1,n);
    l=matrix(rep(qt,N2),ncol=ncol(qt),byrow=TRUE);  # number of levels
    
    C=matrix(0,N2,n);
    C[,1:n]=U[,1:n]/l;  # transform to [0,1]^(n)
    
    step=1;
    for (w in 1:total_step-1) {  # total steps to stop
      if(max(b-a)>=delta){
        #update a,b,c
        a = ifelse(as.matrix(c_opt)-sr*c > 0, as.matrix(c_opt)-sr*c,0)
        b = ifelse(as.matrix(c_opt)+sr*c < 1, as.matrix(c_opt)+sr*c,1)
        c=(b-a)/2;
        
        # transform from [0,1]^(n+1) to [a,b]
        # use the same C
        Cab=matrix(rep(a,N2), ncol = ncol(t(a)), byrow= TRUE) + matrix(rep(b-a,N2), ncol = ncol(t(b-a)), byrow= TRUE)*C
        # transform from [a,b] to A_(n+1)
        A=matrix(0,N2,n);
        A1=matrix(0,N2,n);
        for (i in 1:n){
          ones = matrix(1,N2,1)
          B=1-(n-1)*ell_d1
          A1[,i]=B*ones;
          for (j in i:n){
            A1[,i]=A1[,i]*(ifelse(Cab[,j]>0,1,-1)*abs(Cab[,j])^(1/j));
          }
          A[,i]=(i-1)*ell_d1*ones+A1[,i];
        }
        A=rbind(A,t_opt);  # compare with the former t_opt
        
        C_value=vector();
        
        for (i in 1:nrow(A)){
          knots_int = c(min(Time),(max(Time)-min(Time))*A[i,]+min(Time),max(Time))
          knots = unique(knots_int)
          nknots = length(knots)
          norder = 4
          nbasis = length(knots) + norder - 2
          bsbasis = create.bspline.basis(rangeval = c(min(Time),max(Time)),nbasis,norder,knots)
          phi=my_eval.basis(Time, bsbasis)
          Mmat = ginv(t(phi)%*%phi)%*%t(phi)
          chat = Mmat%*%lyt
          ini=c(v=v0,gamma=gamma0,p=p0,w=chat)
          par_hat=nlminb(ini, objective = ode_equation2_obj2,tau=A[i,],Time=Time, lyt = lyt,n=n,control = list(rel.tol = 1e-3),lambda=lambda)
          
          vh = par_hat$par[1]
          gammah = par_hat$par[2]
          ph = par_hat$par[3]
          wh = par_hat$par[4:(3+nbasis)]
          par_hat=c(vh,gammah,ph,wh)
          C_value[i]=ode_equation2_obj2(par=par_hat,tau=A[i,],Time=Time,lyt=lyt,n=n,lambda=lambda);
        }
        print(paste("......loading......",w,"-----",round(w/total_step,2)*100,"%"))
        for (l in 1:length(C_value)){
          if (is.nan(C_value[l])|is.infinite(C_value[l])){
            C_value[l] = 1000000
          }
        }
        Cmin=min(C_value);
        t_opt=A[which(C_value==Cmin),];   # optimal on A
        
        CMIN[w+1,n]=Cmin;  # record min C for every n
        
        
        T[[n]][w+1,]=(max(Time)-min(Time))*t_opt[1:n]+min(Time)
        
        # back to [a,b]^(n+1)
        c_opt=vector();
        c_opt[n]=((t_opt[n]-(n-1)*ell_d1)/B)^n;
        for (i in 1:(n-1)){
          c_opt[i]=((t_opt[i]-(i-1)*ell_d1)/(t_opt[i+1]-i*ell_d1))^i;
        }
        
        step=step+1;  # record steps for each n
      }
      else{
        break
      }
    }
    
    T_opt=vector();
    T_opt[1:n]=((max(Time)-min(Time))*t_opt[1:n]+min(Time));
    
    Tn[[n]]=T_opt;   # record change of opt design for each n
    CMINn[n]=min(abs(CMIN[,n]));   #record change of min C for each n
    if (n %% 10 == 0){
      print(paste("loading......",n,"-----",round(n/N,2)*100,"%"))
    }
  }
  
  Rn=list();
  for (i in 1:length(N_list)){
    Rn[[i]]=c(Tn[[N_list[i]]],CMINn[N_list[i]])
  }
  return(Rn)
  #Rn knots=1:N;CMINn=N+1
}