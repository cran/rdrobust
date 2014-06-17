### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014
### version 0.6  17Jun2014

rdbinselect = function(y, x, data, subset = NULL, c=0, p=4, numbinl=NULL, numbinr=NULL, 
                          binselect="es", lowerend=NULL, upperend=NULL, scale=1,
                          hide=FALSE, title=NULL, x.label=NULL, y.label=NULL, 
                          x.lim=NULL, y.lim=NULL, model = FALSE, frame = FALSE) {

  call <- match.call()
  if (missing(data)) 
  data <- environment(formula)
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  na.ok <- complete.cases(x) & complete.cases(y)
  x <- x[na.ok]
  y <- y[na.ok]
  
  if (frame) {
    dat.out <- data.frame(x, y)
  }

  if (is.null(lowerend)) {
	    lowerend = min(x)
	}
  if (is.null(upperend)) {
	  upperend = max(x)
	}
	x_low = lowerend
	x_upp = upperend

	size=sum(x>=x_low & x<=x_upp)
  x=x[x>=x_low & x<=x_upp]
  y=y[x>=x_low & x<=x_upp]

	x_l = x[x<c]; x_r = x[x>=c]	
  y_l = y[x<c];	y_r = y[x>=c]
	x.min = min(x);	x.max = max(x)
	range_l = max(x_l) - min(x_l)
	n_l = length(x_l)
	range_r = max(x_r) - min(x_r)
	n_r = length(x_r)
	n = n_l + n_r
  
  #####********************* ERRORS
  exit=0
	if (c<=x.min | c>=x.max){
		print("c should be set within the range of x")
		exit = 1
	}

	if (p<=0 ){
		print("p should be a positive number")
		exit = 1
	}

	if (scale<=0 ){
		print("scale should be a positive number")
		exit = 1
	}

	p_floor = floor(p)/p

	if (p_floor!=1) {
		print("p should be an integer number")
		exit = 1
	}

	if (exit>0) {
		stop()
	}

	p1 = p+1
  compute =0

  rp_l = matrix(NA,n_l,p1);  rp_r = matrix(NA,n_r,p1)
  for (j in 1:p1) {
    rp_l[,j] = x_l^(j-1)
    rp_r[,j] = x_r^(j-1)
  }
  gamma_p1_l = lm(y_l~rp_l-1)$coeff   
  gamma_p1_r = lm(y_r~rp_r-1)$coeff
  
  mu0_p1_l = rp_l%*%gamma_p1_l;	mu0_p1_r = rp_r%*%gamma_p1_r
  
  J_star_l_orig=numbinl
  J_star_r_orig=numbinr
  
  if (is.null(J_star_l_orig) & is.null(J_star_l_orig)) {
  compute = 1  
  
  y_l.sq = y_l^2
  y_r.sq = y_r^2
  gamma_p2_l = lm(y_l.sq~rp_l-1)$coeff
  gamma_p2_r = lm(y_r.sq~rp_r-1)$coeff
  
	### Bias w/sample

	drp_l = matrix(NA,n_l,p);	drp_r = matrix(NA,n_r,p)
	for (j in 1:p) {
		drp_l[,j] = j*x_l^(j-1)
		drp_r[,j] = j*x_r^(j-1)
	}
  mu1_hat_l = drp_l%*%(gamma_p1_l[2:p1])
	mu1_hat_r = drp_r%*%(gamma_p1_r[2:p1])

  ######################### ES
  ind.l = order(x_l); ind.r = order(x_r)
  x.i.l = x_l[ind.l]; x.i.r = x_r[ind.r]
  y.i.l = y_l[ind.l]; y.i.r = y_r[ind.r]

  dxi.l=(x.i.l[2:length(x.i.l)]-x.i.l[1:(length(x.i.l)-1)])
  dxi.r=(x.i.r[2:length(x.i.r)]-x.i.r[1:(length(x.i.r)-1)])
  dyi.l=(y.i.l[2:length(y.i.l)]-y.i.l[1:(length(y.i.l)-1)])
  dyi.r=(y.i.r[2:length(y.i.r)]-y.i.r[1:(length(y.i.r)-1)])

  x.bar.i.l = (x.i.l[2:length(x.i.l)]+x.i.l[1:(length(x.i.l)-1)])/2
  x.bar.i.r = (x.i.r[2:length(x.i.r)]+x.i.r[1:(length(x.i.r)-1)])/2
  rp.i_l  = matrix(NA,n_l-1,p+1); rp.i_r= matrix(NA,n_r-1,p+1)
  drp.i_l = matrix(NA,n_l-1,p); drp.i_r = matrix(NA,n_r-1,p)

  for (j in 1:p1) {
    rp.i_l[,j] = x.bar.i.l^(j-1)
    rp.i_r[,j] = x.bar.i.r^(j-1)
  }

  for (j in 1:p) {
    drp.i_l[,j] = j*x.bar.i.l^(j-1)
    drp.i_r[,j] = j*x.bar.i.r^(j-1)
  }
  mu0.i_hat_l = rp.i_l%*%gamma_p1_l;  mu0.i_hat_r = rp.i_r%*%gamma_p1_r
  mu1.i_hat_l = drp.i_l%*%(gamma_p1_l[2:p1]);  mu1.i_hat_r = drp.i_r%*%(gamma_p1_r[2:p1])
  mu2.i_hat_l = rp.i_l%*%gamma_p2_l;  mu2.i_hat_r = rp.i_r%*%gamma_p2_r

  sigma2_hat_ul = mu2.i_hat_l - mu0.i_hat_l^2
  sigma2_hat_ur = mu2.i_hat_r - mu0.i_hat_r^2

  B.es.hat.l = ((c-x.min)^2/12)*sum(dxi.l*mu1.i_hat_l^2)
  B.es.hat.r = ((x.max-c)^2/12)*sum(dxi.r*mu1.i_hat_r^2)
  V.es.hat.l = (n/(4*(c-x.min)))*sum(dxi.l^2*dyi.l^2)
  V.es.hat.r = (n/(4*(x.max-c)))*sum(dxi.r^2*dyi.r^2)
  C.es.hat.l = (2*B.es.hat.l)/V.es.hat.l
  C.es.hat.r = (2*B.es.hat.r)/V.es.hat.r
  J.es.hat.l = floor((C.es.hat.l*(n_l+n_r))^(1/3))
  J.es.hat.r = floor((C.es.hat.r*(n_l+n_r))^(1/3))

  V.es.chk.l = (n/(2*(c-x.min)))*sum(dxi.l^2*sigma2_hat_ul)
  V.es.chk.r = (n/(2*(x.max-c)))*sum(dxi.r^2*sigma2_hat_ur)
  C.es.chk.l = (2*B.es.hat.l)/V.es.chk.l
  C.es.chk.r = (2*B.es.hat.r)/V.es.chk.r
  J.es.chk.l = floor((C.es.chk.l*(n_l+n_r))^(1/3))
  J.es.chk.r = floor((C.es.chk.r*(n_l+n_r))^(1/3))

  B.es.hat.dw.l = ((c-x.min)^2/(12*n))*sum(mu1.i_hat_l^2)
  B.es.hat.dw.r = ((x.max-c)^2/(12*n))*sum(mu1.i_hat_r^2)
  V.es.hat.dw.l = (0.5/(c-x.min))*sum(dxi.l*dyi.l^2)
  V.es.hat.dw.r = (0.5/(x.max-c))*sum(dxi.r*dyi.r^2)
  C.es.hat.dw.l = (2*B.es.hat.dw.l)/V.es.hat.dw.l
  C.es.hat.dw.r = (2*B.es.hat.dw.r)/V.es.hat.dw.r
  J.es.hat.dw.l = floor((C.es.hat.dw.l*(n_l+n_r))^(1/3))
  J.es.hat.dw.r = floor((C.es.hat.dw.r*(n_l+n_r))^(1/3))

  ######################### QS
  V.qs.hat.l = (n/(2*n_l))*sum(dxi.l*dyi.l^2)
  V.qs.hat.r = (n/(2*n_r))*sum(dxi.r*dyi.r^2)
  B.qs.hat.l = (n_l^2/72)*sum(dxi.l^3*mu1.i_hat_l^2)
  B.qs.hat.r = (n_r^2/72)*sum(dxi.r^3*mu1.i_hat_r^2)
  C.qs.hat.l = (2*B.qs.hat.l)/V.es.hat.l
  C.qs.hat.r = (2*B.qs.hat.r)/V.es.hat.r
  J.qs.hat.l = floor((C.qs.hat.l*(n_l+n_r))^(1/3))
  J.qs.hat.r = floor((C.qs.hat.r*(n_l+n_r))^(1/3))
	
  V.qs.chk.l = (n/n_l)*sum(dxi.l*sigma2_hat_ul)
  V.qs.chk.r = (n/n_r)*sum(dxi.r*sigma2_hat_ur)
  C.qs.chk.l = (2*B.qs.hat.l)/V.qs.chk.l
  C.qs.chk.r = (2*B.qs.hat.r)/V.qs.chk.r
  J.qs.chk.l = floor((C.qs.chk.l*(n_l+n_r))^(1/3))
  J.qs.chk.r = floor((C.qs.chk.r*(n_l+n_r))^(1/3))

  V.qs.hat.dw.l = (1/(2*n_l))*sum(dxi.l^0*dyi.l^2)
  V.qs.hat.dw.r = (1/(2*n_r))*sum(dxi.r^0*dyi.r^2)
  B.qs.hat.dw.l = (n_l^1/48)*sum(dxi.l^2*mu1.i_hat_l^2)
  B.qs.hat.dw.r = (n_r^1/48)*sum(dxi.r^2*mu1.i_hat_r^2)
  C.qs.hat.dw.l = (2*B.qs.hat.dw.l)/V.es.hat.dw.l
  C.qs.hat.dw.r = (2*B.qs.hat.dw.r)/V.es.hat.dw.r
  J.qs.hat.dw.l = floor((C.qs.hat.dw.l*(n_l+n_r))^(1/3))
  J.qs.hat.dw.r = floor((C.qs.hat.dw.r*(n_l+n_r))^(1/3))
  
  if (binselect=="es" ) {
    J_star_l_orig = J.es.hat.l
    J_star_r_orig = J.es.hat.r
  }
  
  if (binselect=="espr") {
    J_star_l_orig = J.es.chk.l
    J_star_r_orig = J.es.chk.r
  }
  
  if (binselect=="esdw") {
    J_star_l_orig = J.es.hat.dw.l
    J_star_r_orig = J.es.hat.dw.r
  }
  
  if (binselect=="qs") {
    J_star_l_orig = J.qs.hat.l
    J_star_r_orig = J.qs.hat.r
  }
  
  if (binselect=="qspr" ) {
    J_star_l_orig = J.qs.chk.l
    J_star_r_orig = J.qs.chk.r
  }
  
  if (binselect=="qsdw" ) {
    J_star_l_orig = J.qs.hat.dw.l
    J_star_r_orig = J.qs.hat.dw.r
  }
  }

  J_star_l = scale*J_star_l_orig
  J_star_r = scale*J_star_r_orig

  bin_x_l = rep(0,length(x_l)); bin_x_r = rep(0,length(x_r))
  jump_l = range_l/J_star_l;jump_r = range_r/J_star_r;
  
  if (binselect=="es" |binselect=="espr" |binselect=="esdw" ) {
    jumps_l=seq(min(x_l),max(x_l),jump_l)
    jumps_r=seq(min(x_r),max(x_r),jump_r)
    binselect_type="Evenly-Spaced"
  }
  else if (binselect=="qs" |binselect=="qspr" |binselect=="qsdw" ) {
    jumps_l=quantile(x_l,probs=seq(0,1,1/J_star_l))
    jumps_r=quantile(x_r,probs=seq(0,1,1/J_star_r))
    binselect_type="Quantile-Spaced"
  }
  
  for (k in 1:(J_star_l-1)) {
	  bin_x_l[x_l>=jumps_l[k] & x_l<jumps_l[k+1]] = -J_star_l+k-1 
	}
    bin_x_l[x_l>=jumps_l[(J_star_l)]] = -1
    
	for (k in 1:(J_star_r-1)) {
	  bin_x_r[x_r>=jumps_r[k] & x_r<jumps_r[k+1]] = k 
	}
	bin_x_r[x_r>=jumps_r[(J_star_r)]] = J_star_r
  
  bin_xlmean=bin_ylmean=rep(0,J_star_l)
	bin_xrmean=bin_yrmean=rep(0,J_star_r)
  	for (k in 1:(J_star_l)) {
	  bin_xlmean[k]=mean(c(jumps_l[k],jumps_l[k+1]))
	  #bin_xlmean[k]=mean(x_l[bin_x_l==-k])
	  bin_ylmean[J_star_l-k+1]=mean(y_l[bin_x_l==-k])
  	}
	for (k in 1:(J_star_r)) {
	  bin_xrmean[k]=mean(c(jumps_r[k],jumps_r[k+1]))
	  #bin_xrmean[k]=mean(x_r[bin_x_r==k])
	  bin_yrmean[k]=mean(y_r[bin_x_r==k]) 
	}
  
  bin_x=c(bin_x_l,bin_x_r)
  bin_xmean=c(bin_xlmean,bin_xrmean)
	bin_ymean=c(bin_ylmean,bin_yrmean)
	x_sup = c(x_l, x_r)
	#y_hat = c(mu0_p1_l, mu0_p1_r)
  
  if (hide=="FALSE") {
  
  if (is.null(title)){
    title="RD Bin Select"
  }
  
  if (is.null(x.label)){
    x.label="X axis"
  }
  
  if (is.null(y.label)){
    y.label="Y axis"
  }
  
  if (is.null(x.lim)){
    x.lim=c(min(x_l),max(x_r))
  }
  
  if (is.null(y.lim)){
    y.lim=c(min(c(y_l,y_r)),max(c(y_l,y_r)))
  }
    plot(bin_xmean[order(bin_xmean)],bin_ymean[order(bin_xmean)], main=title, xlab=x.label, ylab=y.label, ylim=y.lim, xlim=x.lim, pch=20)
	  points(x_l[order(x_l)],mu0_p1_l[order(x_l)],type="l") 
	  points(x_r[order(x_r)],mu0_p1_r[order(x_r)],type="l")  
    abline(v=c)
  }

  if (compute==1) {
    
    
    tabl1.str=matrix(NA,5,2)
    tabl1.str[1,] = formatC(c(n_l,n_r),digits=0, format="f")
    tabl1.str[2,] = formatC(c(p,p),digits=0, format="f")
    tabl1.str[3,] = formatC(c(J_star_l,J_star_r),digits=0, format="f")  
    tabl1.str[4,] = formatC(c(scale,scale),digits=0, format="f") 
    tabl1.str[5,] = formatC(c(jump_l,jump_r),digits=4, format="f") 
    rownames(tabl1.str)=c("Number of Obs.","Poly. Order","Number of Bins","Scale","Bin Length")
    colnames(tabl1.str)=c("Left","Right")
    
    results=matrix(NA,5,2)
    results[1,] = c(n_l,n_r)
    results[2,] = c(p,p)
    results[3,] = c(J_star_l,J_star_r)
    results[4,] = c(scale,scale)
    results[5,] = c(jump_l,jump_r)
    rownames(results)=c("Number of Obs.","Poly. Order","Number of Bins","Scale","Bin Length")
    colnames(results)=c("Left","Right")
    
    coef=matrix(NA,p+1,2)
    coef[,1] = c(gamma_p1_l)
    coef[,2] = c(gamma_p1_r)
    colnames(coef)=c("Left","Right")
    out=list(method=binselect_type,results=results,coef=coef,tabl1.str=tabl1.str)
    
    out$call <- match.call()
    class(out) <- "rdbinselect"
    return(out)
  }
}

#rdbinselect <- function(y,x, ...) UseMethod("rdbinselect")

#rdbinselect.default <- function(y,x, ...){
#  est <- rdbinselectEst(y,x, ... )
#  est$call <- match.call()
#  class(est) <- "rdbinselect"
#  est
#}

print.rdbinselect <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(paste("Method: ",x$method),quote=F)
  cat("\n\n")
  print(x$tabl1.str,quote=F)
}

summary.rdbinselect <- function(object,...) {
  TAB <- object$results
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.rdbinselect"
  res
}

print.summary.rdbinselect <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
}