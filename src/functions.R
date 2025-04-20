neg.log.epsilon = 52 * log(2) # ~ 36

rmvnorm <- function(mean, sigma){
    p = nrow(sigma)
    R = chol(sigma, pivot = TRUE)
    R = R[, order(attr(R, "pivot"))]
    res = crossprod(R,rnorm(p))
    res+mean
}

ldmvnorm <- function(x, mean, sigma) {
  p = length(x)
  chol_sigma = chol(sigma)
  log_det_sigma = 2 * sum(log(diag(chol_sigma)))
  diff = x - mean
  z = backsolve(chol_sigma, diff, transpose = TRUE)
  qf = crossprod(z)[1,1]
  ld = -0.5 * (qf + p * log(2 * pi) + log_det_sigma)
  
  return(ld)
}


vec <- function(x){t(t(as.vector(x)))}

ldet<-function(M){
  determinant(M,logarithm=T)$modulus[1]
}

catn <- function(s){
  cat(deparse(substitute(s)),s,'\n')
}

lse2 = function (x, y) {
  if (is.infinite(x)) return(y)
  m = max(x, y); d = m - min(x, y)
  if(!is.na(d) & !is.nan(d)){if (d > neg.log.epsilon) return(m)} else return(y)
  m + log(1 + exp(-d))
}

lse <- function (l) Reduce(lse2, l)

chol_solve <- function (C, v){backsolve(C, backsolve(C, v, transpose = TRUE))}

dt.scaled <- function(x,nu,mu,s,log=T){
  dt((x - mu)/s, nu, log = T) - log(s)
}

dinvgamma <- function(x, shape, scale, log = FALSE) {
  log_den <- shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
  den <- if (log) log_den else exp(log_den)
  return(den)
}

dinvchisq <- function(x, df, scale, log = FALSE) {
  shape <- df / 2
  return(dinvgamma(x, shape, scale, log = log))
}

rinvchisq <- function(n, df, scale) {
  x <- rchisq(n, df)
  x <- (df * scale) / x
  return(x)
}

create_prior_cov <- function(psi,lam,h,ny){
  #h is num harmonics
  #Y is num years
  phi = psi*exp(lam*(1-(1:h)));
  D = diag(phi)
  #
  Dc = D - tcrossprod(phi)/sum(phi) #condition on sum(param)=0
  Dk = kronecker(diag(2*(ny-1)),Dc[-h,-h])
  return(Dk)
}

compute_Ezlogq <- function(y,X,th,s2,nu){
  #use this to Elogq when k>1
  k = dim(th)[2]
  at = (nu+1)/2
  sapply(1:k, function(j){bt = nu/2 + (y - X%*%th[,j])^2/(2*s2[j]); digamma(at)-log(bt)})
}

compute_Elogq <- function(y,X,th,s2,nu){
  #for k=1 only
  at = (nu+1)/2
  bt = nu/2 + (y - X%*%th)^2/(2*s2)
  digamma(at)-log(bt)
}

bernstein<-function(j,k,t){
 choose(k-1,j-1)*t^(j-1)*(1-t)^(k-j)
}

compute_p_z <- function(t,k){
  n = length(t)
  p_z = array(0,dim=c(n,k))
  p_z[,1] = bernstein(1,k,t)
  for(j in 2:k){
    p_z[,j] = bernstein(j,k,t)
  }
  return(p_z)
}

# FIX compute_p_z calculation of n
compute_p_z_disc <- function(t,k){
  #This n needs to be changed to include the total number of possible observations 365*20/16
  n = length(t) 
  is = 1:n #floor(t*n)
  p_z = sapply(1:k, function(j) choose(is,j-1)*choose(n-is,k-j)/choose(n,k-1))
  return(p_z)
}

compute_p_z_zm1 <- function(t,n,k){
  # The columns of p_z_zm1 are the conditions on zm1.
  # Since for each condition there can only be the current and next state,
  # we only need to store the probability for the current state for each condition.
  # Then the probability of the next state is p(z=c+1 | zm1=c) = 1 - p(z=c | zm1=c)  
  p_z_zm1     = array(dim=c(n-1,k)) 
  for(j in 1:k){
    p_z_zm1[,j] = ((1-t[2:n])/(1-t[1:(n-1)]))^(k-j) / (((1-t[2:n])/(1-t[1:(n-1)]))^(k-j) + (k-j)*(t[2:n]-t[1:(n-1)])/(1-t[1:(n-1)]) * ((1-t[2:n])/(1-t[1:(n-1)]))^(k-j-1));
  }
  p_z_zm1[(n-1),] = c(rep(0,k-1),1)
  return(p_z_zm1)
}

compute_p_z_zm1_nonID <- function(t,n,k){
  # The columns of p_z_zm1 are the conditions on zm1.
  # Since for each condition there can only be the current and next state,
  # we only need to store the probability for the current state for each condition.
  # Then the probability of the next state is p(z=c+1 | zm1=c) = 1 - p(z=c | zm1=c)  
  p_z_zm1     = array(0,dim=c(n-1,k,k))
  for(j in 1:k){
    for(h in j:k)
    p_z_zm1[,j,h] = choose(k-j,h-j)*(1 - (1-t[2:n])/(1-t[1:(n-1)]))^(h-j) * ((1-t[2:n])/(1-t[1:(n-1)]))^(k-h);
  }
  p_z_zm1[,k,1:(k-1)]=0
  p_z_zm1[,k,k]=1
  return(p_z_zm1)
}

compute_p_z_zm1_geometric <- function(t,n,k,ps){
  p_z_zm1 = array(dim=c(n-1,k)) 
  for(j in 1:(k-1)){
    p_z_zm1[,j] = ps[j]
  }
  p_z_zm1[(n-1),] = c(rep(0,k-1),1)
  p_z_zm1[,k] = 1
  return(p_z_zm1)
}

compute_p_z_zm1_disc <- function(n,k,p_z){
  p_z_zm1     = array(dim=c(n-1,k)) 
  p_z_zm1[,1] = p_z[2:n,1]/p_z[1:(n-1),1]; p_z_zm1[n-1,1]=0;
  if(k>2){
  for(j in 2:(k-1)){
    p_z_zm1[,j] = (apply(p_z[2:n,1:j],1,"sum") - apply(p_z[1:(n-1),1:j-1,drop=F],1,"sum"))/p_z[1:(n-1),j]  
    }
  }
  p_z_zm1[,k] = 1
  return(p_z_zm1)
}

compute_logp_z_zm1 <- function(p_z_zm1,n,k){
  #v3 changing this to be nxkxk
  #There are n-1 transitions, k-1 changes, 2 probabilities
  logp_z_zm1 = array(-Inf, dim=c(n-1,k,k))
  for(i in 1:(k-1)){
    logp_z_zm1[,i,i] = log(p_z_zm1[,i])
    logp_z_zm1[,i,i+1] = log(1-p_z_zm1[,i])
  }
  logp_z_zm1[is.na(logp_z_zm1)] = -Inf
  logp_z_zm1[,k,k] = 0
  return(logp_z_zm1)
}

compute_logp_z_zm1_nonID <- function(p_z_zm1,n,k){
  #v3 changing this to be nxkxk
  #There are n-1 transitions, k-1 changes, 2 probabilities
  logp_z_zm1 = log(p_z_zm1)
  logp_z_zm1[is.na(logp_z_zm1)] = -Inf
  logp_z_zm1[,k,k] = 0
  return(logp_z_zm1)
}

compute_logp_y_z <- function(y,X,nu,theta,s2,intercept=T,normal=F){
  n = dim(X)[1]; k = length(s2);
  if(intercept){
    if(normal){
      logp_y_z = sapply(1:k, function(j) dnorm(x = y,mean = X*theta[j],sd = s2[j]^0.5,log=T))
    }else{
      logp_y_z = sapply(1:k, function(j) dt.scaled(y,nu,X*theta[j],s2[j]^0.5,log=T))
    }
  }else{
    if(normal){
      logp_y_z = sapply(1:k, function(j) dnorm(x = y,mean = X%*%theta[,j,drop=F],sd = s2[j]^0.5,log=T))
    }else{
      logp_y_z = sapply(1:k, function(j) dt.scaled(y,nu,X%*%theta[,j,drop=F],s2[j]^0.5,log=T))
    }
  }
  return(logp_y_z)
}

forward <-function(logp_y_z,logp_z_zm1){
  n = dim(logp_y_z)[1]; k = dim(logp_y_z)[2]
  logalpha = array(dim=c(n,k))
  logalpha[1,] = -Inf
  logalpha[1,1] = logp_y_z[1,1]
  #
  for(i in 2:n){
    w = logalpha[i-1,]+logp_z_zm1[i-1,,]
    logalpha[i,] = logp_y_z[i,] + apply(w,2,lse)
  }
  return(logalpha)
}

backward <- function(logp_y_z,logp_z_zm1){
  n = dim(logp_y_z)[1]; k = dim(logp_y_z)[2]
  logbeta = array(dim=c(n,k))
  logbeta[n,] = 0
  #
  for(i in (n-1):1){
    w = logbeta[i+1,] + logp_y_z[i+1,] + t(logp_z_zm1[i,,])
    logbeta[i,] = apply(w,2,lse)
  }
  return(logbeta)
}

compute_Ez <-function(logalpha,logbeta,logpy_k){
  n = dim(logalpha)[1]
  Ez = exp(logalpha + logbeta - logpy_k)
  return(Ez) 
}

compute_Ezzm1 <- function(logalpha,logbeta,logp_y_z,logp_z_zm1){
  n = dim(logalpha)[1]; k = dim(logp_y_z)[2];
  logp_y = lse(logalpha[n,]) #logp_y = apply(logalpha+logbeta,1,lse)
  Ezzm1 = array(dim=c(n-1,k,k)) #n-1,zm1,z
  #for j in z_{i-1}
  for(j in 1:k){
    Ezzm1[,j,] = exp(sweep((logp_y_z[-1,] + logp_z_zm1[,j,] + logbeta[-1,]),1,(logalpha[-n,j] - logp_y),"+"))
  }
  return(Ezzm1)  
}

compute_prior <-function(k,theta,s2,Phi_inv){
  logpr = 0
  #FIXME: Phi_inv was set to 0 when intercept=T, but there's a better solution
  if(sum(Phi_inv)==0){Phi_inv = array(1e-6,dim=c(1,1))}
  ch = chol(Phi_inv)
  Phi = chol_solve(ch,diag(dim(Phi_inv)[1]))
  for(j in 1:k){
    logpr = logpr + ldmvnorm(x=theta[,j],sigma = s2[j]*Phi) - log(sqrt(s2[j]))
    #dnorm(theta[j],0,1/s2[j],log=T)
  }
  return(logpr)
}

am_theta <- function(Eqz,Phi_inv,X,y,th,intercept=F){
  if(intercept){
    k = dim(Eqz)[2]; n = length(y); p = dim(X)[2]; th = array(0,dim=c(p,k))
    for(j in 1:k){
      Xw = X*Eqz[,j]
      XtWy = crossprod(Xw, y)
      ch = try(chol(crossprod(X,Xw) + Phi_inv),silent=T)
      if(class(ch)[1]!='try-error'){th[,j] = chol_solve(ch,XtWy)} else{th[,j] = rep(0,p)}
    }
    return(list(th=th))
  }else{
    k = dim(Eqz)[2]; n = length(y); p = dim(X)[2];
    dw = 6; up = 7;#IAH start at 6,7; int slp 2,3
    Xu = X[,up:p]
    y_tild = y - Xu%*%th[up:p,1]
    th = array(0,dim=c(p,k))
    #
    y_star = rep(y,k); 
    X_star = do.call(rbind, replicate(k, Xu, simplify=FALSE)); 
    w_star = vec(Eqz)[,1]
    #
    Xd = X[,1:dw]
    Phi_inv_d = Phi_inv[1:dw,1:dw]
    for(j in 1:k){
      w = Eqz[,j]
      Xdw = Xd*w
      XdtWy = crossprod(Xdw, y_tild)
      ch = try(chol(crossprod(Xd,Xdw) + Phi_inv_d),silent=T)
      th[1:dw,j] = chol_solve(ch,XdtWy)
      y_star[(1:n)+n*(j-1)] = y - Xd%*%th[1:dw,j]
    }
    Xw_star = X_star*w_star
    XtW_star_y_star = crossprod(Xw_star, y_star)
    ch = try(chol(crossprod(X_star,Xw_star) + Phi_inv[up:p,up:p]),silent=T)
    if(class(ch)[1]!='try-error'){th_star = chol_solve(ch,XtW_star_y_star)}
    th[up:p,] = th_star
    return(list(th=th))
  }
}

am_s2 <- function(Eqz,Ez,y,X,th,lam,Phi_inv){
  k = dim(Eqz)[2]; n = length(y); p = dim(X)[2];
  s2 = array(0,dim=c(k))
  for(j in 1:k){
    s2[j] = crossprod(Eqz[,j],(y-X%*%th[,j])^2) + crossprod(th[,j],Phi_inv)%*%th[,j]
  }
  s2 = rep(sum(s2)/(n + k*p + 2),k)
  return(s2)
}

am_ps<-function(Ezzm1){
  n = dim(Ezzm1)[1]; k = dim(Ezzm1)[2]+1;
  alph = (n/k)/10
  beta = .1
  k = dim(Ezzm1)[2]
  ps = sapply(1:(k-1), function(j){
    num = (sum(Ezzm1[,j,j])+(alph-1))
    den = (sum(Ezzm1[,j,j+1]) + (beta - 1))
    return(num/(num+den))
  })
  #ps[ps<1e-6] = 0.5
  #ps[ps>=1]=0.5; ps[ps<=0]=0.5; ps = sort(ps);
  return(ps)
}

compute_Ey_x <-function(t,t_pred,Ez,th,X_pred,rmj){
  if(length(rmj)>0){Ez = Ez[,-rmj]}
  k = dim(Ez)[2]
  n = dim(X_pred)[1]
  Ey_x = array(0,dim=n)
  if(k==1){
    dim(th) = c(length(th),1)
    Ey_x = Ey_x + (X_pred%*%th[,1])[,1]
  }else{
    for(j in 1:k){
        Ez_pred = approx(t,Ez[,j,drop=F],n = length(t_pred))$y
        Ey_x = Ey_x + Ez_pred*(X_pred%*%th[,j,drop=F])[,1]
    }
  }
  return(Ey_x)
}

taus_from_Ez<-function(Ezzm1,t,rmj=c()){
  #catn(rmj)
  #if(length(rmj)>0){Ezzm1 = Ezzm1[,,-rmj]; Ezzm1 = Ezzm1[,-rmj,]}
  n = dim(Ezzm1)[1]; k = dim(Ezzm1)[2]
  tau=c()
  #Ezzm1 => t,zm1,z
  #cat('taus_from_Ez')
  for(j in 1:(k-1)){
    i=1;p=0
    while(p<0.5){
      #catn(c(i,j,j+1)); catn(p)
      p = p + Ezzm1[i,j,j+1]
      i=i+1
    }
    tau = c(tau,i)
    #catn(tau)
  }
  return(tau)
}

fit <- function(X,y,k,h,psi,lam,discrete=F,geometric=F,intercept=F,normal=F,id=0,tol=1e-3,maxit=10,chtimes=NULL,nu=3){
  #Assume X df has "datetime" as a single column, ordered by date with y
  geometric=F; if(normal){nu=NULL}#else{nu=3}#nu=0.1}
  a = init_data_priors(y,X,psi,lam,h,intercept)
  y=a$y; t=a$t; X=a$X; Phi_inv=a$Phi_inv; X_pred=a$X_pred
  n = length(y); p = dim(X)[2];
  #
  if(k==1){
    ch = try(chol(crossprod(X)+ Phi_inv))
    theta = chol_solve(ch,crossprod(X,y))
    s2 = sum((y-X%*%theta)^2)/n
    logpy_k_old=0
    Eq = rep(1,n)
    for(it in 1:maxit){
      if(!normal){
        Eq = (nu + 1)/(nu + (y-X%*%theta)[,1]^2/s2)
        Elogq = compute_Elogq(y,X,theta,s2,nu)
      }
      ch = try(chol(crossprod(X,diag(Eq))%*%X + Phi_inv))
      theta = chol_solve(ch,crossprod(X,diag(Eq)%*%y))
      s2 = (sum(Eq*(y-X%*%theta)^2)+crossprod(theta,Phi_inv)%*%theta)[1,1]/(n + p + 1)
      if(!normal){ 
        #nu = am_nu(y,theta,s2,nu,Eq,Elogq)
        logpy_k = sum(dt.scaled(y,nu,X%*%theta,s2^0.5,log=T))
        #cat('nu = ',nu,'\n')
      }else{
        logpy_k = sum(dnorm(x = y,mean = X%*%theta,sd = s2^0.5,log=T))
      }
      if(abs(logpy_k_old - logpy_k)/abs(logpy_k_old) < tol) break else logpy_k_old = logpy_k
    }
    return(list(t=t,logpy_k = logpy_k, Ez = array(1,dim=c(n,1)), taus = NULL,theta=theta,s2=s2, X=X,Eq=Eq,ps=NULL,Phi_inv = Phi_inv,nu=nu,X_pred=X_pred))
  }else{
    ps=NULL; theta = array(0,dim=c(p,k))
    if(normal){Eq = ones(n,k)}
    if(discrete){
      p_z = compute_p_z_disc(t,k)
      p_z_zm1 = compute_p_z_zm1_disc(n,k,p_z)
      logp_z_zm1 = compute_logp_z_zm1(p_z_zm1,n,k)
      theta = am_theta(p_z,Phi_inv,X,y,theta,intercept)$th
    }else {
      p_z = compute_p_z(t,k)
      p_z_zm1 = compute_p_z_zm1_nonID(t,n,k)
      theta = am_theta(p_z,Phi_inv,X,y,theta,intercept)$th
      logp_z_zm1 = compute_logp_z_zm1_nonID(p_z_zm1,n,k)#
    }
    if(class(theta)[1]=='try-error')return(theta)
    s2 = am_s2(p_z,p_z,y,X,theta,lam,Phi_inv)
    rmj = which((s2<1e-100)|is.nan(s2)); 
    logpy_k_old = 0
    for(it in 1:maxit){
      # E step
      #
      logp_y_z = compute_logp_y_z(y,X,nu,theta,s2,intercept,normal)
      logalpha = forward(logp_y_z,logp_z_zm1)
      logbeta = backward(logp_y_z,logp_z_zm1)
      logpy_k = lse(logalpha[n,])
      Ez = compute_Ez(logalpha,logbeta,logpy_k)
      #Ezzm1 = compute_Ezzm1(logalpha,logbeta,logp_y_z,logp_z_zm1)
      Eqz = array(0,dim=c(n,k))
      if(!normal){
        for(j in 1:k){
          Eqz[,j] = Ez[,j] * (nu + 1)/(nu + (y-X%*%theta[,j])^2/s2[j]) 
        }
        if(length(rmj)>0){Eq = rowSums(Eqz[,-rmj])}else{Eq = rowSums(Eqz)} #Eqz already multiplies Ez
      } else{Eqz = Ez}
      #
      # M step
      at = am_theta(Eqz,Phi_inv,X,y,theta,intercept)
      theta = at$th
      ch_list = at$ch_list
      if(class(theta)[1]=='try-error')return(theta)
      s2 = am_s2(Eqz,Ez,y,X,theta,lam,Phi_inv)
      rmj = which((s2<1e-100)|is.nan(s2)); 
      if(!normal){
        Ezlogq = compute_Ezlogq(y,X,theta,s2,nu)
        if(length(rmj)>0){Elogq = rowSums(Ezlogq[,-rmj]*Ez[,-rmj])} else{Elogq = rowSums(Ezlogq*Ez)}
        #nu = am_nu(y,theta,s2,nu,Eq,Elogq)
        #cat('nu = ',nu,'\n')
      }
      #
      # Check convergence
      #cat("Marginal log(p(y)): ", logpy_k,"\n\n")
      if(abs(logpy_k_old - logpy_k)/abs(logpy_k_old) < tol) break else logpy_k_old = logpy_k
    }
    logp_y_z = compute_logp_y_z(y,X,nu,theta,s2,intercept,normal)
    logalpha = forward(logp_y_z,logp_z_zm1)
    logbeta = backward(logp_y_z,logp_z_zm1)
    logpy_k = lse(logalpha[n,])
    Ez = compute_Ez(logalpha,logbeta,logpy_k)
    #Ezzm1 = compute_Ezzm1(logalpha,logbeta,logp_y_z,logp_z_zm1)
    #
    return(list(t=t,logpy_k=logpy_k,  Ez=Ez, theta=theta, s2=s2, X=X, Eq=Eq, ps=ps, p_z_zm1=p_z_zm1, Phi_inv = Phi_inv, nu=nu,X_pred=X_pred,rmj=rmj))
  }
}

BPP <- function(X,y,K,nu=3,intercept=F){
  n = length(y); quants=c(0.5); h=2; psi=0.1; lam=1;
  Ez_list = list(); Ey_x_list=list(); theta_list=list();rmj_list=list(); Phi_inv_list=list();sigma_list=list()
  logpy_k = c()
  for(k in 1:K){
    cat("Running EM for k =",k,"regimes\n")
    b = fit(X,y,k,h,psi,lam,nu=nu,discrete=F,geometric=F,intercept=intercept,normal=F)
    #
    Ez_list[[k]] = b$Ez; theta_list[[k]] = b$theta; Phi_inv_list[[k]] = b$Phi_inv; sigma_list[[k]] = b$s2
    logpy_k = c(logpy_k,b$logpy_k)
  }
  logpk_y = compute_pk_y(n,logpy_k,b$t,Phi_inv_list,sigma_list,intercept=intercept,nonInf_pk=F,discrete=F)
  res = lapply(quants, function(qu) bayes_est(qu,n,K,Ez_list,exp(logpk_y),samples=F))
  theta = theta_list[[which.max(logpk_y)]]
  return(list(taus=res[[1]]$taus,logpk_y=logpk_y,theta=theta))
}


fit_EM <- function(X,y,K,h=2,psi=0.1,lam=1,nu=3,quants=c(0.025,0.5,0.975),discrete=F,geometric=F,intercept=F,normal=F,nonInf_pk=F){
  n = length(y)
  Ez_list = list(); Ey_x_list=list(); Eq_list=list(); theta_list=list();rmj_list=list(); Phi_inv_list=list();sigma_list=list()
  logpy_k = c()
  K_sub = K
  for(k in 1:K){
    cat("Running EM for k =",k,"regimes\n")
    b = fit(X,y,k,h,psi,lam,nu=nu,discrete=discrete,geometric=geometric,intercept=intercept,normal=normal)
    #
    Ez_list[[k]] = b$Ez; Eq_list[[k]] = b$Eq; theta_list[[k]] = b$theta; Phi_inv_list[[k]] = b$Phi_inv; sigma_list[[k]] = b$s2
    logpy_k = c(logpy_k,b$logpy_k)
    if(!intercept) Ey_x_list[[k]] = compute_Ey_x(b$t,b$X_pred[,2],b$Ez,b$theta,b$X_pred,b$rmj)
    #if(length(b$rmj)>0){K_sub=k-1; break}
  }
  logpk_y = compute_pk_y(n,logpy_k,b$t,Phi_inv_list,sigma_list,intercept=intercept,nonInf_pk=nonInf_pk,discrete=discrete)
  #catn(logpk_y)
  res = lapply(quants, function(qu) bayes_est(qu,n,K,Ez_list,exp(logpk_y),samples=F))
  Ey_x = array(0,dim=5000)
  q = array(0,dim=n)
  for(k in 1:K){
    if(!intercept) Ey_x = Ey_x + Ey_x_list[[k]]*exp(logpk_y)[k]
    if(!normal) q = q + Eq_list[[k]]*exp(logpk_y)[k]
  }
  theta = theta_list[[which.max(logpk_y)]]
  return(list(res=res,logpk_y=logpk_y,t=b$t,t_pred = b$X_pred[,2],Ey_x=Ey_x,q=q,theta=theta))
}

fit_EM2 <- function(X,y,K,h,psi,lam,nu=3,quants=c(0.025,0.5,0.975),discrete=F,geometric=F,intercept=F,normal=F){
  n = length(y)
  Ez_list = list(); Phi_inv_list=list();sigma_list=list()
  logpy_k = c()
  for(k in 1:K){
    cat("Running EM for k =",k,"regimes\n")
    b = fit(X,y,k,h,psi,lam,nu=nu,discrete=discrete,geometric=geometric,intercept=intercept,normal=normal)
    #
    Ez_list[[k]] = b$Ez; Phi_inv_list[[k]] = b$Phi_inv; sigma_list[[k]] = b$s2
    logpy_k = c(logpy_k,b$logpy_k)
  }
  logpk_y_T = compute_pk_y(n,logpy_k,b$t,Phi_inv_list,sigma_list,intercept=intercept,nonInf_pk=T)
  logpk_y_F = compute_pk_y(n,logpy_k,b$t,Phi_inv_list,sigma_list,intercept=intercept,nonInf_pk=F)
  #
  res_T = lapply(quants, function(qu) bayes_est(qu,n,K,Ez_list,exp(logpk_y_T),samples=F))
  res_F = lapply(quants, function(qu) bayes_est(qu,n,K,Ez_list,exp(logpk_y_F),samples=F))
  return(list(res_T=res_T,res_F = res_F))
}

bayes_est <- function(qu,n,K,Z,pk_y,samples=T){
  sp <- matrix(0, n, K)
  for(k in 1:K){
    Ez <- Z[[k]]
    Ez[is.nan(Ez)] <- 0
    if (samples) {
      Pk <- t(apply(Ez, 1, function(x) tabulate(x, nbins = K) / length(x)))
    } else {
      Pk <- Ez
    }
    js <- seq_len(k)
    sp[, js] <- sp[, js] + pk_y[k] * Pk[, js]
  }
  idx <- sp >= qu
  zm  <- max.col(idx, ties.method = "first")
  taus <- which(c(FALSE, diff(zm) > 0))
  list(taus = taus, zm = zm)
}

full_setup <- function(X,h,psi,lam){
  #Assume X df has "datetime" as a single column, ordered by date with y
  #y = y[order(X[,1])]
  #X = X[order(X[,1]),,drop=F]
  dt = X
  X = design_setup(X,h)
  t = X[,2]
  #
  #ny number of years
  ny = length(unique(X[,"year"]))
  X = X[,colnames(X) != "year"]
  Phi = create_prior_cov(psi,lam,h,ny)
  Phi_inv = solve(Phi) 
  p = dim(X)[2]
  n = dim(X)[1]
  ch = try(chol(crossprod(X)+ Phi_inv))
  #theta = chol_solve(ch,crossprod(X,y))
  return(list(X=X,Phi_inv=Phi_inv,n=n,ch=ch))
}

logprior <- function(K,n,t,Phi_inv_list,sigma_list,nonInf_pk=F,intercept=F,discrete=F){
  nc = 0
  logpr = c()
  for(k in 1:K) {
    if(!intercept){
      Phi_inv = Phi_inv_list[[k]];
      logpr[k] = (-1)^nonInf_pk * ( k*sum(log((1 - t[2:n]+1e-8)/(1 - t[1:(n-1)]))) - k*(dim(Phi_inv)[1]/2)*log(2*pi) + k*(1/2)*ldet(Phi_inv))
      }else if(!discrete){
        logpr[k] = (-1)^nonInf_pk * ( k*sum(log((1 - t[2:n]+1e-8)/(1 - t[1:(n-1)]))) - (k/2)*log(2*pi))
      }else{
        logpr[k] = 1/choose(n-1,k-1) #n counts the first observation; in paper it does not
      }
  }
  logpr = logpr - lse(logpr)
  #
  return(logpr)
}

compute_pk_y <- function(n,logpy_k,t,Phi_inv_list,sigma_list,intercept=F,nonInf_pk=F,discrete=F,gibbs=F){
  K = length(logpy_k);
  if(!gibbs){
    if(!intercept) pen = sapply(1:K,function(k) -(6*k+44 + 1)*log(n)/2 ) else pen = -((1:K)+1)*log(n)/2 #+1 for sigma
  }else{
    pen = 0
  }
  logpr = logprior(K,n,t,Phi_inv_list,sigma_list,nonInf_pk,intercept,discrete)
  logZ = lse(logpy_k + logpr + pen)
  logpk_y = (logpy_k + logpr + pen) - logZ 
  return(logpk_y)
}

design_setup <- function(dt,h){
  #dt is a dataframe with one column "datetime"
  period=365;
  dt$time = as.numeric(dt$datetime)/86400#(24h/d*60m/h*60s/m)
  dt$t = (dt$time - min(dt$time))/(max(dt$time) - min(dt$time) + 0.0001)
  tot_days = as.numeric(max(dt$datetime)-min(dt$datetime))
  dt$intercept = 1
  dt = dt[,c("intercept","t","datetime")]
  X = create_ia_design(dt,h,tot_days,period)
  return(X)
}

create_design <- function(orig_df,dt,h,tot_days,period){
  #dt is datetime
  #assume orig_df has column t for time
  #Changed to be sin,sin2,sin3,cos,cos2,cos3 to accomodate order of prior.
  #orig_df$year = floor(orig_df$t*tot_days/period) - min(floor(orig_df$t*tot_days/period))
  mY = min(as.numeric(format(dt,"%Y")))
  bd = as.POSIXct(paste0(mY,"-01-01"), format="%Y-%m-%d")
  pt = as.numeric(min(dt) - bd)/tot_days
  for (i in 1:h){
    orig_df = cbind(orig_df, sin(2*pi*i*(orig_df$t+pt)*tot_days/period))
    names(orig_df)[3+i] = paste("sin",i,sep = "")
  }
  for (i in 1:h){
    orig_df = cbind(orig_df, cos(2*pi*i*(orig_df$t+pt)*tot_days/period))
    names(orig_df)[3+h+i] = paste("cos",i,sep = "")
  }
  return(orig_df)
}

bdiag<-function(...){
  lmat = list(...)
  nrs = sapply(lmat, nrow)
  ncs = sapply(lmat, ncol)
  nr = sum(nrs)
  nc = sum(ncs)
  M = matrix(0,nr,nc)
  prevr=1; prevc = 1
  for(j in 1:length(lmat)){
    nxtr = (prevr+nrs[j]-1)
    nxtc = (prevc+ncs[j]-1)
    M[prevr:nxtr,prevc:nxtc] = lmat[[j]]
    prevr = nxtr+1
    prevc = nxtc+1
  }
  return(M)
}

create_ia_design <- function(orig_df,h,tot_days,period){
  #assume orig_df = cbind(intercept, t, datetime)
  dt = orig_df$datetime; orig_df$datetime = NULL
  orig_df$year = as.numeric(format(dt,"%Y")); orig_df$year = orig_df$year - min(orig_df$year);
  orig_df = create_design(orig_df,dt,h,tot_days,period)
  nc = dim(orig_df)[2]
  df_years = list() 
  for(i in 1:length(unique(orig_df$year))){
    xi = unique(orig_df$year)[i]
    df_years[[i]]<-as.matrix(orig_df[orig_df$year==xi,(nc-2*h+1):nc,drop=F])
  }
  harmonic_ma = do.call(bdiag,df_years)
  m = length(unique(orig_df$year))
  for(xi in 1:(m-1)){
    #sin_0,sin1_0,sin2_0,cos_0,cos2_0, cos3_0 | sin_1,sin1_1,sin2_1,cos_1,cos2_1,cos3_1
    #Subtract last sin from the remaining sines for each year; repeat for cos
    harmonic_ma[,(2*h*xi+1):(2*h*xi+h-1)] = sweep(harmonic_ma[,(2*h*xi+1):(2*h*xi+h-1),drop=F],1,harmonic_ma[,(2*h*xi+h)],"-") 
    harmonic_ma[,((2*xi+1)*h+1):(2*h*(xi+1)-1)] = sweep(harmonic_ma[,((2*xi+1)*h+1):(2*h*(xi+1)-1),drop=F],1,harmonic_ma[,(2*h*(xi+1))],"-") 
  }
  #Remove last harmonic for each year except the first
  rmve = c()
  for(xi in 1:(m-1)){
    rmve = c(rmve,(2*h*xi+h),(2*h*(xi+1)))
  }
  model_ma = as.matrix(cbind(orig_df,harmonic_ma[,-rmve]))
  model_ma = model_ma[,-((4+2*h):(3+4*h))]
  return(model_ma)
}

am_nu <- function(y,th,s2,nu,Eq,Elogq){
  n = length(y); gr = 0.02; hess=1; i=1
  while(abs(gr/hess) > 0.01){
    gr = sum( 1/2 * log(nu/2) + 1/2 - 1/2 * digamma(nu/2) + 1/2 * Elogq - 1/2 * Eq )
    hess = n * ( 1/2 * 1/nu - 1/4 * trigamma(nu/2) )
    nu = nu - gr/hess
    i = i + 1
  }
  return(nu)
}

total_loss <- function(tau,tau_s,K){
 #tau in (0,1)
 tau = c(tau,rep(1,K-length(tau)))
 tau_s = c(tau_s,rep(1,K-length(tau_s)))
 sum(abs(tau-tau_s))
}

predict_update <- function(k,X,y,q,Phi_inv,theta,s2,p_z_zm1){
  n = dim(X)[1]; pzi_Yim1 = array(0,dim=c(n-1,k)); pzi_Yi = array(0,dim=c(n,k))
  pzi_Yi[1,] = c(1,rep(0,k-1))
  py_zqt = sapply(1:k,function(j) dnorm(x=y,mean = X%*%theta[,j,drop=F],sd=sqrt(s2[j]/q)))
  #
  for(i in 1:(n-1)){
    #predict
    for(j in 1:k){
      if(j==1){
        pzi_Yim1[i,1] = p_z_zm1[i,1,1]*pzi_Yi[i,1]
      }else{
        #pzi_Yim1[i,j] = p_z_zm1[i,j]*pzi_Yi[i,j]+(1-p_z_zm1[i,j-1])*pzi_Yi[i,j-1] #replace this with a sum from 1 to j
        pzi_Yim1[i,j] = sum(sapply(1:j, function(l) p_z_zm1[i,l,j]*pzi_Yi[i,l]))
      }
    }
    #update
    Zi = crossprod(pzi_Yim1[i,],py_zqt[i+1,])    
    pzi_Yi[i+1,] = pzi_Yim1[i,]*py_zqt[i+1,]/Zi
  }
  return(list(pzi_Yi=pzi_Yi,pzi_Yim1=pzi_Yim1,py_zqt=py_zqt))
}

sample_z <- function(k,X,y,q,Phi_inv,theta,s2,p_z_zm1){
  n = dim(X)[1]; z = array(0,dim=n)
  z[n] = k
  pz_ZY = array(0,dim=c(n,k))
  pu = predict_update(k,X,y,q,Phi_inv,theta,s2,p_z_zm1)
  pzi_Yi = pu$pzi_Yi
  pzi_Yim1 = pu$pzi_Yim1
  #
  for(i in (n-1):1){
    r = (1:z[i+1]) #replace this with 1:z[i+1]
    #pzzm1 = c(p_z_zm1[i,z[i+1]],(1-p_z_zm1[i,z[i+1]-1])) #change this to 1:j to j transitions
    pzzm1 = p_z_zm1[i,1:z[i+1],z[i+1]] #change this to 1:j to j transitions
    Zi = crossprod(pzi_Yi[i,r],pzzm1)
    pz_ZY[i,r] = pzi_Yi[i,r]*pzzm1 / Zi #FIXME:
    #z[i] = rcat(n=1,p=sapply(1:k,function(j) max(0,pz_ZY[i,j])))
    z[i] = sample(1:k,1,sapply(1:k,function(j) max(0,pz_ZY[i,j])),replace=T)
    if(z[i]==1){z[1:i] = 1; break}
  }
  return(list(z=z,pzi_Yim1=pzi_Yim1))
}

sample_theta <- function(k,X,y,z,q,s2,Phi_inv,theta_old,init=NULL,extra=F){
  p = dim(Phi_inv)[1]
  theta = array(0,dim=c(p,k))
  logpth_hat_y = 0
  for(j in 1:k){
    zj <- (z==j)
    if(sum(zj)>0){
      Lam = crossprod(X[zj,,drop=F],diag(q[zj],nrow=length(q[zj])))%*%X[zj,,drop=F] + Phi_inv
      ch = try(chol(Lam))
      Lam_inv = s2[j]*chol_solve(ch,diag(dim(Lam)[1]))
      Etheta = chol_solve(ch,crossprod(X[zj,,drop=F],diag(q[zj],nrow=length(q[zj]))%*%y[zj]))
      theta[,j] = rmvnorm(n=1,mean=Etheta,sigma=Lam_inv)
    }else{theta[,j] = rep(0,p)}#theta_old[,j]}
    #
    if(extra) logpth_hat_y = logpth_hat_y + ldmvnorm(init$theta[,j],mean=Etheta,sigma=Lam_inv)
  }
  return(list(theta=theta,logpth_hat_y=logpth_hat_y))
}

sample_sigma <- function(k,X,y,z,q,theta,Phi_inv,init=NULL,extra=F){
  s2 = array(0,dim=k); n = length(y)
  logps2_hat_y = 0
  df = n + k*dim(theta)[2]
  s = sum(sapply(1:k, function(j){zj <- (z==j);(crossprod(q[zj],(y[zj]-X[zj,]%*%theta[,j,drop=F])^2) + crossprod(theta[,j],Phi_inv)%*%theta[,j])} ))/df
  s2 = rinvchisq(1, df, scale=sqrt(s))
  if(extra) logps2_hat_y = logps2_hat_y + dinvchisq(init$s2[1], df, scale=sqrt(s),log=T)
  return(list(s2=rep(s2,k),logps2_hat_y=logps2_hat_y))
}

sample_q <-function(k,X,y,z,nu,theta,s2){
  n = length(y)
  alph = (nu + 1)/2 
  bet = sapply(1:n, function(i) (nu + (y[i]-crossprod(X[i,],theta[,z[i]]))^2/s2[z[i]])/2)
  #bet = 0.5*((nu + (y-(X%*%theta)[,z])^2) / s2[z])
  q = rgamma(n=n,shape = alph, rate = bet)
  return(q)
}

sample_ps <-function(z,k,init=NULL,extra=F){
  logpp_hat_y=NULL; n = length(z)
  sz = as.vector(table(z))[-k]
  a = (n/k)/10; b = 0.1 #page 232 Chib 1996, w/ noninformative n/k regime duration
  ps = sapply(1:(k-1),function(j) rbeta(1, sz[j]+a, b+1))
  if(extra) logpp_hat_y = sum(sapply(1:(k-1),function(j) dbeta(init$ps[j], sz[j]+a, b+1,log=T)))
  return(list(ps=ps,logpp_hat_y=logpp_hat_y))
}

init_data_priors <- function(y,X,psi,lam,h,intercept){
  if(!intercept){
    y = y[order(X[,1])]
    X = X[order(X[,1]),,drop=F]
    dt = X
    X = design_setup(X,h)#[,1:7]
    t = X[,2]
    #
    nn = 5000
    dt_pred = data.frame(seq(min(dt$datetime),max(dt$datetime),length.out=nn))
    colnames(dt_pred) = "datetime"
    X_pred = design_setup(dt_pred,h)#[,1:7]
    X_pred = X_pred[,colnames(X_pred) != "year"]
    #
    #ny number of years
    ny = length(unique(X[,"year"]))
    X = X[,colnames(X) != "year"]
    Phi = create_prior_cov(psi,lam,h,ny)
    Phi_inv = solve(Phi) 
    flat_pr_row = array(0,dim=c(dim(X)[2]-dim(Phi)[2],dim(Phi)[2])) #avg harm has flatpr; contr has norm pr
    flat_pr_col = array(0,dim=c(dim(X)[2],dim(X)[2]-dim(Phi)[2]))
    nip = 1e-6 #noninformative precision for intercept slope and baseline harmonics
    diag(flat_pr_col) = c(nip,5,rep(nip,2*h))#rep(.1/psi,2*h))#exp((1:(2*h))))#0
    Phi_inv = cbind(flat_pr_col, rbind(flat_pr_row,Phi_inv))
    #Phi_inv = diag(dim(X)[2])#array(0,dim=dim(crossprod(X)))
  }else{
    n = length(y)
    t = X[1:n]
    X = as.matrix(rep(1,n),dim=c(n,1))
    Phi_inv = as.matrix(0,dim=c(1))
    X_pred=NULL
  }
  return(list(y=y, t=t,X=X,Phi_inv=Phi_inv,X_pred=X_pred))
}

gibbs_sampler <- function(X,y,k,h,psi,lam,init=NULL,discrete=F,geometric=F,intercept=F,normal=F,nsamp=500){
  # if(!intercept) Assume X df has "datetime" as a single column, ordered by date with y
  # if(intercept) Assume X is time from 0 to 1
  #
  idp = init_data_priors(y,X,psi,lam,h,intercept)
  y=idp$y; t=idp$t; X=idp$X; Phi_inv=idp$Phi_inv; nu = init$nu
  #X = X[,1:6]; Phi_inv = Phi_inv[1:6,1:6]
  n = length(y); p = dim(X)[2]
  #
  if(k==1){
    z = rep(1,n); theta = init$theta; s2 = init$s2; q = init$q
    #
    # Need extra samples of [q|s2_hat] for the posterior density of theta
    modsel_th = c()
    for(it in 1:nsamp){
      q = sample_q(k,X,y,z,nu,theta,s2)
      #
      t_list = sample_theta(k,X,y,z,q,s2,Phi_inv,theta,init,extra=T)
      modsel_th = c(modsel_th, t_list$logpth_hat_y)
      theta = t_list$theta
    }
    # Gibbs sample
    modsel_s2 = c()
    for(it in 1:nsamp){
      theta = sample_theta(1,X,y,z,q,s2,Phi_inv,theta,init)$theta
      s_list = sample_sigma(1,X,y,z,q,theta,Phi_inv,init,extra=T)
      modsel_s2 = c(modsel_s2, s_list$logps2_hat_y)
      s2 = s_list$s2
      if(normal) q = rep(1,n) else q = sample_q(1,X,y,z,nu,theta,s2)
    }
    return(list(Z = array(1,dim=c(n,nsamp)), modsel_th=lse(modsel_th)-log(nsamp), modsel_s2=lse(modsel_s2)-log(nsamp),modsel_ps=0))
  }else{#k>1
    p_z = compute_p_z(t,k)
    p_z_zm1 = compute_p_z_zm1_nonID(t,n,k) #dim = (n-1)xk
    # Storage
    Z = array(0,dim=c(n,nsamp));
    #Initialize from EM
    theta = init$theta; s2 = init$s2; q = init$q; z = init$z
    #
    # Need extra samples of [z|theta_hat,s2_hat] for the posterior density of ps
    modsel_ps=0
    # Need extra samples of [z,q,ps|s2_hat] for the posterior density of theta
    modsel_th = c()
    theta = init$theta
    s2 = init$s2
    q = init$q
    z = init$z
    for(it in 1:nsamp){
      q = sample_q(k,X,y,z,nu,theta,s2)
      z = sample_z(k,X,y,q,Phi_inv,theta,s2,p_z_zm1)$z
      st = sample_theta(k,X,y,z,q,s2,Phi_inv,theta,init,extra=T)
      t_list = st; 
      modsel_th = c(modsel_th, t_list$logpth_hat_y)
      theta = t_list$theta
    }
    #
    # Gibbs sample
    modsel_s2 = c()
    for(it in 1:nsamp){
      st = sample_theta(k,X,y,z,q,s2,Phi_inv,theta,init)
      theta = st$theta
      #
      s_list = sample_sigma(k,X,y,z,q,theta,Phi_inv,init,extra=T)
      modsel_s2 = c(modsel_s2, s_list$logps2_hat_y)
      s2 = s_list$s2
      #
      if(normal) q = rep(1,n) else q = sample_q(k,X,y,z,nu,theta,s2)
      z = sample_z(k,X,y,q,Phi_inv,theta,s2,p_z_zm1)$z
      Z[,it] = z
  }
  return(list(Z=Z,modsel_th=lse(modsel_th)-log(nsamp),modsel_s2=lse(modsel_s2)-log(nsamp),modsel_ps=lse(modsel_ps)-log(nsamp)))
  }
}

fit_gibbs <- function(X,y,K,h,psi,lam,discrete=F,geometric=F,intercept=T,normal=F,nsamp=100){
  # Init:
  n = length(y); Z_list=list(); em_list = list(); modsel_list = list()
  #
  print(K)
  for(k in 1:K){
    cat("Running EM initialization for k =",k,"regimes\n")
    b = fit(X,y,k,h,psi,lam,discrete=discrete,geometric=geometric,intercept=intercept,normal=normal)
    b$z = z_from_taus(n,b$taus); t = b$t
    em_list[[k]] = b
    #
    cat("Running Gibbs Sampler for k =",k,"regimes\n")
    gs = gibbs_sampler(X,y,k,h,psi,lam,init=b,discrete=discrete,geometric=geometric,intercept=intercept,normal=normal,nsamp=nsamp)
    #
    Z_list[[k]]=gs$Z; 
    modsel_list[[k]] = list(modsel_th = gs$modsel_th,
                            modsel_s2 = gs$modsel_s2,
                            modsel_ps = gs$modsel_ps)
  }
  logpk_y = pk_y_gibbs(n,t,K,em_list,modsel_list)
  #catn(logpk_y)
  res = bayes_est(0.5,n,K,Z_list,exp(logpk_y))
  return(list(res=res, logpk_y = logpk_y,K=K,t=t))
}

pk_y_gibbs <- function(n,t,K,em_list,modsel_list,intercept=T,nonInf_pk=F){
  logpy_k = c()
  for(k in 1:K){
    init = em_list[[k]]
    init2 = modsel_list[[k]]
    logpr = compute_prior(k,init$theta,init$s2,init$Phi_inv)
    logpost = init2$modsel_th + init2$modsel_s2 + init2$modsel_ps
    #catn(logpr);catn(logpost);catn(init$logpy_k)
    #
    logpy_k = c(logpy_k, init$logpy_k + logpr - logpost)
  }
  logpk_y = compute_pk_y(n,logpy_k,t,init$Phi_inv,init$s2,intercept=intercept,nonInf_pk=nonInf_pk,gibbs=T) #n,logpy_k,t,Phi_inv_list,sigma_list,intercept=F,nonInf_pk=F,discrete=F
  return(logpk_y)
}

z_from_taus <- function(n,taus){
  k = length(taus)+1
  taus=c(0,taus,n)
  z = unlist(sapply(1:k, function(j) rep(j,taus[j+1]-taus[j]) ))
  return(z)
}

whamming_loss <- function(n,t,z,z_true){
  dt = c(t[1],t[2:n] - t[1:(n-1)])
  sum(abs(z-z_true)*dt)
}

calibrate_uncertainty <- function(taus_med, taus_025,taus_975){
  #length(taus_q025) <= length(taus_med) wp 1
  #length(taus_med) <= length(taus_q975) wp 1
  if(length(taus_025)!=length(taus_med)){
    temp_taus = taus_med
    rm = c()
    for(l in 1:length(taus_025)){
      pos = which.min(abs(taus_med-taus_025[l]))
      rm = c(rm,pos)
    }
    taus_025 = sort(c(taus_025,temp_taus[-rm]))
  }
  #
  if(length(taus_975)!=length(taus_med)){
    temp_taus = c()
    for(l in 1:length(taus_med)){
      pos = which.min(abs(taus_med-taus_975[l]))
      temp_taus = c(temp_taus,taus_975[pos])
    }
    taus_975 = temp_taus
  }
  return(list(taus_025=taus_025,taus_975=taus_975))
}

calc_omit_commit<-function(taus,taus_true){
  commission = 0
  omission = 0
  if(length(taus_true)==2){
    commission = commission + length(taus) - length(taus_true)
  }else{
    for(j in 2:(length(taus_true)-1)){
      ct = sum((taus_true[j]<=taus)&(taus_true[j+1]>taus))
      if(ct == 0){
        omission=omission+1
      }else if(ct > 1){
        commission = commission + ct - 1
      }
    }
  }
  return(list(commission=commission,omission=omission))
}

# calc_omit_commit_window<-function(taus,taus_true,window=0.04353531){
#   commission = 0
#   omission = 0
#   if(length(taus_true)==2){
#     commission = commission + length(taus) - length(taus_true)
#   }else{
#     for(j in 2:(length(taus_true)-1)){
#       ct = sum((abs(taus_true[j]-taus)<=window))
#       catn(sum((abs(taus_true[j]-taus)<=window)))
#       if(ct == 0){
#         omission=omission+1
#       }else if(ct > 1){
#         commission = commission + ct - 1
#       }
#     }
#   }
#   return(list(commission=commission,omission=omission))
# }

calc_omit_commit_window<-function(ta,tt,window=0.04353531){
  #assume ta and tt do NOT include end points
  tau_true = tt
  tau = ta
  commission = 0
  omission = 0
  if(length(tt)==0){
    commission = commission + length(ta) - length(tt)
  }else if(length(ta)>0){
    for(j in 1:length(ta)){
      ct = as.numeric(all((abs(ta[j]-tt)>=window)))
      #catn(all((abs(ta[j]-tt)>=window)))
      commission = commission + ct
      if(ct==0) tt = tt[-which.min(abs(ta[j]-tt))] #remove closest tt so that it is not reused
    }
  }
  tt = tau_true
  if(length(ta)==0){
    omission = omission + length(tt) - length(ta)
  }else if(length(tt)>0){
    for(j in 1:length(tt)){
      ct = as.numeric(all((abs(tt[j]-ta)>=window)))
      #catn(all((abs(tt[j]-ta)>=window)))
      omission = omission + ct
      if(ct==0) ta = ta[-which.min(abs(tt[j]-ta))] #remove closest tt so that it is not reused
    }
  }
  #if(omission>0) omission = omission/length(tau_true)
  #if(commission>0) commission = commission/length(tau)
  return(c(omission,commission))
}

sample_BPP <- function(ti,n,k){
  p_z_zm1 = compute_p_z_zm1_nonID(ti,n,k) #n,j,h
  z = 1; chpts=c()
  for(i in 1:(n-1)){  #i=1 is time=2
    zp1 = sample(z[i]:k,1,prob=p_z_zm1[i,z[i],z[i]:k])
    if(zp1>z[length(z)]){chpts = c(chpts,i)}
    z = c(z, zp1)
    if(zp1==k){break}
  }
  return(chpts)
}
