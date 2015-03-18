### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014
### version 0.6  17Jun2014
### version 0.61 03Sep2014
### version 0.7  14Oct2014
### version 0.8  04Feb2015

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  chol2inv(chol(crossprod(x)))
}

qrreg = function(x,y,w,s2=0,var.comp=TRUE, ...) {
  M.X = sqrt(w)*x
  X.M.X_inv = qrXXinv(M.X) 
  X.M.Y = crossprod(M.X,sqrt(w)*y)
  beta.hat = X.M.X_inv%*%X.M.Y
  Psi.hat=Sigma.hat=0
  if (var.comp==TRUE) {
    Psi.hat = crossprod((w*s2*w)*x,x)
    Sigma.hat = crossprod(Psi.hat%*%X.M.X_inv,X.M.X_inv)
  }
  output = list(X.M.X_inv=X.M.X_inv, X.M.Y=X.M.Y, beta.hat=beta.hat, Psi.hat=Psi.hat, Sigma.hat=Sigma.hat)
  return(output)
}

kweight = function(X, c,  h,  kernel){
  u = (X-c)/h
  if (kernel=="epanechnikov" | kernel=="epa") {
    w = (0.75*(1-u^2)*(abs(u)<=1))/h
  }
  else if (kernel=="uniform" | kernel=="uni") {
    w = (0.5*(abs(u)<=1))/h
  }
  else {
    w = ((1-abs(u))*(abs(u)<=1))/h
  }
  return(w)	
}

bwconst = function(p,v,kernel){
  if (kernel=="epanechnikov" | kernel=="epa" | kernel==3) {
    K.fun = function(u) {(0.75*(1-u^2)*(abs(u)<=1))}
  }
  else if (kernel=="uniform" | kernel=="uni" | kernel==2) {
    K.fun = function(u) {(0.5*(abs(u)<=1))}
  }
  else  {
    K.fun = function(u) {((1-abs(u))*(abs(u)<=1))}
  }
  p1 = p+1  
  Gamma_p = Phi_p = matrix(NA,p1,p1)
  Omega_pq = matrix(NA,p1,1)
  for (i in 1:p1) {
    Omega.fun = function(u) {K.fun(u)*(u^(p1))*(u^(i-1))}
    Omega_pq[i] = integrate(Omega.fun,lower=0,upper=1)$value
    for (j in 1:p1) {
      Gamma.fun = function(u) {K.fun(u)*(u^(i-1))*(u^(j-1))}
      Phi.fun   = function(u) {(K.fun(u)^2)*(u^(i-1))*(u^(j-1))}
      Gamma_p[i,j] = integrate(Gamma.fun,lower=0,upper=1)$value
      Phi_p[i,j] = integrate(Phi.fun,lower=0,upper=1)$value
    }
  }
  B_const = solve(Gamma_p)%*%Omega_pq
  V_const = solve(Gamma_p)%*%Phi_p%*%solve(Gamma_p)
  C1 = B_const[v+1,1]
  C2 = V_const[v+1,v+1]
  return(c(C1,C2))
}

rdvce= function(X,y,frd=NULL,p,h,matches,vce,kernel){
m = matches+1
n = length(X)
p1 = p+1
sigma = matrix(0,n,1)
if (vce=="resid") {
	for (k in 1:n) {
		cutoff = matrix(X[k],n,1)
		cutoff1 = X[k]
		W = kweight(X,cutoff1,h,"kernel")
    ind=W>0
    if (sum(ind)>5) {
      w.p=W[ind]; X.p=X[ind]; y.p=y[ind]
      XX.p = matrix(c((X.p-cutoff1)^0, poly(X.p-cutoff1,degree=p,raw=T)),length(X.p),p+1)
      mu0_phat_y = qr.coef(qr(XX.p*sqrt(w.p), tol = 1e-10), sqrt(w.p)*y.p)[1]
      if (is.null(frd)) {
          sigma[k] = (y[k] - mu0_phat_y)^2
      }
      else if (!is.null(frd)) {
          z.p=frd[ind]
          out=qrreg(XX.p, z.p, w.p, var.comp=FALSE) 
          mu0_phat_z = out$beta.hat[1]
          sigma[k] = (y[k] - mu0_phat_y)*(frd[k] - mu0_phat_z)
          }
      }
	}
}
else  {
#y_match_avg = z_match_avg = matrix(NA,n,1)
for (k in 1:n) {
	diffx = abs(X - X[k])
	m.group = sort(unique(diffx))[2:m]
	ind = which(diffx %in% m.group)
	y_match_avg = z_match_avg = mean(y[ind])
  Ji = length(ind)
	if (is.null(frd)) {
	  sigma[k] = (Ji/(Ji+1))*(y[k] - y_match_avg)^2
	} 
  else if (!is.null(frd)) {
      z_match_avg = mean(frd[ind])
      sigma[k] = (Ji/(Ji+1))*(y[k] - y_match_avg)*(frd[k] - z_match_avg)
      }
	}
}
return(sigma)
}

regconst = function(d,h){
  d2 = 2*d+1
  d1 = d+1
  mu = matrix(0,d2, 1)
  mu[1] = 1
  XX = matrix(0,d1,d1)
  for (j in 2:d2) {
    i = j-1
	    if (j%%2==1) {
		      mu[j] = (1/(i+1))*(h/2)^i
	    }
  }
  for (j in 1:d1) {
	  XX[j,] = t(mu[j:(j+d)])
  }
  invXX =solve(XX)
  return(invXX)
}


