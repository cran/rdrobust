### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014
### version 0.6  17Jun2014
### version 0.61 03Sep2014
### version 0.7  14Oct2014
### version 0.8  04Feb2015

rdrobust = function(y, x,  subset = NULL, c=0, p=1, q=2, deriv=0, fuzzy=NULL, h=NULL, b=NULL, rho=NULL, scalepar=1, kernel="tri", bwselect="CCT", scaleregul=1, delta=0.5, cvgrid_min=NULL, cvgrid_max=NULL, cvgrid_length=NULL, cvplot=FALSE, vce="nn", matches=3, level=95, all=FALSE) {
  
  call <- match.call()
  #if (missing(data)) 
  #data <- environment(formula)

  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
    
  if (!is.null(fuzzy)){
    type <- "fuzzy"
    z=fuzzy
    if (!is.null(subset)) 
    z <- z[subset]
    na.ok <- na.ok & complete.cases(z)
  } 
  else {
    type <- "sharp"
  }

  x <- x[na.ok]
  y <- y[na.ok]
  
  if (type == "fuzzy") 
    z <- as.double(z[na.ok])
  #if (frame) {
  #  if (type == "sharp") {
  #    dat.out <- data.frame(x, y)
  #  }
  #  else {
  #    dat.out <- data.frame(x, y, z)
  #  }
  #}
  
  kernel = tolower(kernel)
  bwselect = toupper(bwselect)
  vce = tolower(vce)
  
  X_l=x[x<c];  X_r=x[x>=c]
  Y_l=y[x<c];  Y_r=y[x>=c]
  x_min = min(x);  x_max = max(x)
  N_l = length(X_l);   N_r = length(X_r)
  range_l = abs(max(X_l)-min(X_l));   range_r = abs(max(X_r)-min(X_r))
  N = N_r + N_l
  m = matches + 1
  quant = -qnorm(abs((1-(level/100))/2))
  
  #if (deriv==0 & p==0){
  #  p = 1
  #}
  
  #if (deriv>0 & p==0){
  #  bwselect = "CCT"
  #  p = deriv+1
  #}
  
  if (q==0) {
    q = p+1
  }
  
  p1 = p+1;  q1 = q+1;  exit = 0
  
  #####################################################   CHECK ERRORS

    if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
      print("kernel incorrectly specified")
      exit = 1
    }
    
    if (bwselect!="CCT" & bwselect!="IK" & bwselect!="CV" & bwselect!=""){
      print("bwselect incorrectly specified")
      exit = 1
    }
    
    if (vce!="resid" & vce!="nn" & vce!="" & vce!="s2.pob"){ 
      print("vce incorrectly specified")
      exit = 1
    }
    
    if (c<=x_min | c>=x_max){
      print("c should be set within the range of x")
      exit = 1
    }
    
    if (p<=0 | q<=0 | deriv<0 | matches<=0 ){
      print("p,q,deriv and matches should be positive integers")
      exit = 1
    }
    
    if (p>=q & q>0){
      print("p should be set higher than q")
      exit = 1
    }
    
    if (deriv>p & deriv>0 ){
      print("deriv should be set higher than p")
      exit = 1
    }
    
    p_round = round(p)/p;  q_round = round(q)/q;  d_round = round(deriv+1)/(deriv+1);  m_round = round(matches)/matches
    
    if (p_round!=1 | q_round!=1 | d_round!=1 | m_round!=1 ){
      print("p,q,deriv and matches should be integer numbers")
      exit = 1
    }
    
    if (delta>1 | delta<=0){
      print("delta should be set between 0 and 1")
      exit = 1
    }
    
    if (level>100 | level<=0){
      print("level should be set between 0 and 100")
      exit = 1
    }
    
    if (!is.null(rho)){  
       if (rho<0){
          print("rho should be greater than 0")
          exit = 1
        }
    }
  
    #if (cvgrid_min<0 | cvgrid_max<0 | cvgrid_length<0 ){
    #  print("cvgrid_min, cvgrid_max and cvgrid_length should be positive numbers")
    #  exit = 1
    #}
    
    #if (cvgrid_min>cvgrid_max){
    #  print("cvgrid_min should be lower than cvgrid_max")
    #  exit = 1
    #}
    
    if (exit>0) {
      stop()
    }
      
    if (!is.null(h)) {
      bwselect = "Manual"
    }
  
    if (!is.null(h) & is.null(rho) & is.null(b)) {
      rho = 1
      b = h
    }
  
    if (!is.null(h) & !is.null(rho) ) {
      b = h/rho
    }
    
    ############################################################################################
    #print("Preparing data.") 
    
    if (is.null(h) & bwselect=="IK") {
      rdbws=rdbwselect(y=y, x=x, c=c, rho=rho, p=p, q=q, bwselect="IK", kernel=kernel, precalc=FALSE, scaleregul=scaleregul)
      h = rdbws$bws[1,1];      b = rdbws$bws[1,2]
      bwselect = "IK"
    }
    else if (is.null(h) & bwselect=="CV") {
      rdbws=rdbwselect(y=y, x=x, c=c, rho=rho, p=p, q=q, bwselect="CV", kernel=kernel, precalc=FALSE, delta=delta, cvgrid_min=cvgrid_min, cvgrid_max=cvgrid_max, cvgrid_length=cvgrid_length, scaleregul=scaleregul)
      h = rdbws$bws[1,1];      b = rdbws$bws[1,2]
      bwselect = "CV"
    }
    else if (is.null(h)) {
      rdbws=rdbwselect(y=y, x=x, c=c, rho=rho, p=p, q=q, bwselect="CCT", kernel=kernel, precalc=FALSE, matches=matches, vce=vce, scaleregul=scaleregul, deriv=deriv)
      h = rdbws$bws[1,1];      b = rdbws$bws[1,2]
      bwselect = "CCT"
    }
    
  
    if (kernel=="epanechnikov" | kernel=="epa") {
      kernel_type = "Epanechnikov"
    }
    else if (kernel=="uniform" | kernel=="uni") {
      kernel_type = "Uniform"
    }
    else  {
      kernel_type = "Triangular"
    }
    
  N_l = length(X_l);  N_r = length(X_r)
  wh_l = kweight(X_l,c,h,kernel);  wh_r = kweight(X_r,c,h,kernel);  
  wb_l = kweight(X_l,c,b,kernel);  wb_r = kweight(X_r,c,b,kernel)
  uh_l = (X_l-c)/h;  uh_r = (X_r-c)/h;  uhh_l=uh_l[wh_l>0];  uhh_r=uh_r[wh_r>0]
  Yh_l  = Y_l[wh_l>0];  Yh_r  = Y_r[wh_r>0];  Yb_l  = Y_l[wb_l>0];  Yb_r  = Y_r[wb_r>0]
  Xh_l  = X_l[wh_l>0];  Xh_r  = X_r[wh_r>0];  Xb_l  = X_l[wb_l>0];  Xb_r  = X_r[wb_r>0]
  whh_l = wh_l[wh_l>0]; whh_r = wh_r[wh_r>0]; wbb_l = wb_l[wb_l>0]; wbb_r = wb_r[wb_r>0]
  
  if (type == "fuzzy") {
    T_l  = matrix(z[x<c]);    T_r  = matrix(z[x>=c])
    Th_l = matrix(T_l[wh_l>0]);    Th_r = matrix(T_r[wh_r>0])
    Tb_l = matrix(T_l[wb_l>0]);    Tb_r = matrix(T_r[wb_r>0])
  }
  
  Nh_l = length(Xh_l);  Nb_l = length(Xb_l);  Nh_r = length(Xh_r);  Nb_r = length(Xb_r)
  
  if (Nh_l<5 | Nh_r<5 | Nb_l<5 | Nb_r<5){ 
    stop("Too few observations to compute RD estimates")
  }
      
  X_lp   = matrix(c((X_l-c)^0, poly(X_l-c,degree=p,raw=T)),length(X_l),p+1)
  X_rp   = matrix(c((X_r-c)^0, poly(X_r-c,degree=p,raw=T)),length(X_r),p+1)
  Xh_lp  = matrix(c((Xh_l-c)^0,poly(Xh_l-c,degree=p,raw=T)),length(Xh_l),p+1)
  Xh_rp  = matrix(c((Xh_r-c)^0,poly(Xh_r-c,degree=p,raw=T)),length(Xh_r),p+1)
  
  X_lq   = matrix(c((X_l-c)^0, poly(X_l-c,degree=q,raw=T)),length(X_l),q+1)
  X_rq   = matrix(c((X_r-c)^0, poly(X_r-c,degree=q,raw=T)),length(X_r),q+1)
  Xb_lq  = matrix(c((Xb_l-c)^0,poly(Xb_l-c,degree=q,raw=T)),length(Xb_l),q+1)
  Xb_rq  = matrix(c((Xb_r-c)^0,poly(Xb_r-c,degree=q,raw=T)),length(Xb_r),q+1)
  
  #print("Computing Variance-Covariance Matrix.")
  
    sigmah_l = c(rdvce(X=Xh_l, y=Yh_l, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmah_r = c(rdvce(X=Xh_r, y=Yh_r, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmab_l = c(rdvce(X=Xb_l, y=Yb_l, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmab_r = c(rdvce(X=Xb_r, y=Yb_r, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
  
  #print("Computing RD Estimates.")
  
  factor_p = factorial(seq(0,p,1))
  factor_q = factorial(seq(0,q,1))
  
  out.lp=qrreg(Xh_lp,Yh_l,whh_l,sigmah_l) 
  out.rp=qrreg(Xh_rp,Yh_r,whh_r,sigmah_r)
  out.lq=qrreg(Xb_lq,Yb_l,wbb_l,sigmab_l)
  out.rq=qrreg(Xb_rq,Yb_r,wbb_r,sigmab_r)
  
  tau_lp = factor_p*out.lp$beta.hat
  tau_rp = factor_p*out.rp$beta.hat
  tau_lq = factor_q*out.lq$beta.hat
  tau_rq = factor_q*out.rq$beta.hat
  
  V_lp   = out.lp$Sigma.hat
  V_rp   = out.rp$Sigma.hat
  V_lq   = out.lq$Sigma.hat
  V_rq   = out.rq$Sigma.hat
  
  invGamma_lp = out.lp$X.M.X_inv
  invGamma_rp = out.rp$X.M.X_inv
  invGamma_lq = out.lq$X.M.X_inv
  invGamma_rq = out.rq$X.M.X_inv

  if (b>=h){
    whb_l = wh_l[wb_l>0];   Xb_lp = X_lp[wb_l>0,]
    whb_r = wh_r[wb_r>0];   Xb_rp = X_rp[wb_r>0,]
    Psi_lpq = crossprod(sqrt(whb_l*sigmab_l*wbb_l)*Xb_lp,sqrt(whb_l*sigmab_l*wbb_l)*Xb_lq)
    Psi_rpq = crossprod(sqrt(whb_r*sigmab_r*wbb_r)*Xb_rp,sqrt(whb_r*sigmab_r*wbb_r)*Xb_rq)
  }
  else {
    wbh_l = wb_l[wh_l>0];    Xh_lq = X_lq[wh_l>0,]
    wbh_r = wb_r[wh_r>0];    Xh_rq = X_rq[wh_r>0,]
    Psi_lpq = crossprod(sqrt(whh_l*sigmah_l*wbh_l)*Xh_lp,sqrt(whh_l*sigmah_l*wbh_l)*Xh_lq)
    Psi_rpq = crossprod(sqrt(whh_r*sigmah_r*wbh_r)*Xh_rp,sqrt(whh_r*sigmah_r*wbh_r)*Xh_rq)
  }
  
  Hp = diag(c(1,poly(h,degree=p,raw=T)))

  Cov_l = invGamma_lp%*%Psi_lpq%*%invGamma_lq;  Cov_r = invGamma_rp%*%Psi_rpq%*%invGamma_rq
  v_lp = t(Xh_lp*whh_l)%*%(uhh_l^(p+1));  v_rp = t(Xh_rp*whh_r)%*%(uhh_r^(p+1)) 
  BiasConst_lp = factorial(deriv)*Hp%*%invGamma_lp%*%v_lp;  BiasConst_rp = factorial(deriv)*Hp%*%invGamma_rp%*%v_rp
  Bias_tau = (tau_rq[p+2,1]*BiasConst_rp[deriv+1,1] - tau_lq[p+2,1]*BiasConst_lp[deriv+1,1])*(h^(p+1-deriv)/factorial(p+1))
  
  V_l_cl = factorial(deriv)^2*V_lp[deriv+1,deriv+1] 
  V_r_cl = factorial(deriv)^2*V_rp[deriv+1,deriv+1]  
  V_l_rb = factorial(deriv)^2*V_lp[deriv+1,deriv+1] + factorial(p+1)^2*V_lq[p+2,p+2]*(BiasConst_lp[deriv+1]*h^(p+1-deriv)/factorial(p+1))^2 - 2*factorial(deriv)*factorial(p+1)*Cov_l[deriv+1,p+2]*(BiasConst_lp[deriv+1]*h^(p+1-deriv)/factorial(p+1))
  V_r_rb = factorial(deriv)^2*V_rp[deriv+1,deriv+1] + factorial(p+1)^2*V_rq[p+2,p+2]*(BiasConst_rp[deriv+1]*h^(p+1-deriv)/factorial(p+1))^2 - 2*factorial(deriv)*factorial(p+1)*Cov_r[deriv+1,p+2]*(BiasConst_rp[deriv+1]*h^(p+1-deriv)/factorial(p+1))
  V_cl   = scalepar^2*(V_l_cl + V_r_cl)
  V_rb   = scalepar^2*(V_l_rb + V_r_rb)
    
  tau_Y = c(scalepar*(tau_rp[deriv+1,1] - tau_lp[deriv+1,1]), scalepar*(tau_rp[deriv+1,1] - tau_lp[deriv+1,1] - Bias_tau), scalepar*(tau_rp[deriv+1,1] - tau_lp[deriv+1,1] - Bias_tau))
  se_Y  = c(sqrt(V_cl),sqrt(V_cl),sqrt(V_rb))
  t_Y  =  tau_Y/se_Y
  pv_Y = 2*pnorm(-abs(t_Y))
  ci_Y = matrix(NA,nrow=3,ncol=2)
  rownames(ci_Y)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci_Y)=c("Lower","Upper")
  ci_Y[1,] = c(tau_Y[1] - quant*se_Y[1], tau_Y[1] + quant*se_Y[1])
  ci_Y[2,] = c(tau_Y[2] - quant*se_Y[2], tau_Y[2] + quant*se_Y[2])
  ci_Y[3,] = c(tau_Y[3] - quant*se_Y[3], tau_Y[3] + quant*se_Y[3])
    
  if (type == "fuzzy") {
    # First Stage 
    sigmah_T_l = c(rdvce(X=Xh_l, y=Th_l, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmah_T_r = c(rdvce(X=Xh_r, y=Th_r, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmab_T_l = c(rdvce(X=Xb_l, y=Tb_l, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmab_T_r = c(rdvce(X=Xb_r, y=Tb_r, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    out.lp = qrreg(Xh_lp, Th_l, whh_l, sigmah_T_l) 
    out.rp = qrreg(Xh_rp, Th_r, whh_r, sigmah_T_r)
    out.lq = qrreg(Xb_lq, Tb_l, wbb_l, sigmab_T_l)
    out.rq = qrreg(Xb_rq, Tb_r, wbb_r, sigmab_T_r)
    tau_T_lp = factor_p*out.lp$beta.hat
    tau_T_rp = factor_p*out.rp$beta.hat
    tau_T_lq = factor_q*out.lq$beta.hat
    tau_T_rq = factor_q*out.rq$beta.hat
    V_T_lp   = out.lp$Sigma.hat;    V_T_rp   = out.rp$Sigma.hat
    V_T_lq   = out.lq$Sigma.hat;    V_T_rq   = out.rq$Sigma.hat
    invGamma_T_lp = out.lp$X.M.X_inv;    invGamma_T_rp = out.rp$X.M.X_inv
    invGamma_T_lq = out.lq$X.M.X_inv;    invGamma_T_rq = out.rq$X.M.X_inv
    
    if (b>=h){
      Psi_T_lpq = t(Xb_lp*c(whb_l*sigmab_T_l*wbb_l))%*%Xb_lq
      Psi_T_rpq = t(Xb_rp*c(whb_r*sigmab_T_r*wbb_r))%*%Xb_rq
    }
    else {
      Psi_T_lpq = t(Xh_lp*c(whh_l*sigmah_T_l*wbh_l))%*%Xh_lq
      Psi_T_rpq = t(Xh_rp*c(whh_r*sigmah_T_r*wbh_r))%*%Xh_rq
    } 
    
    Cov_T_l = invGamma_T_lp%*%Psi_T_lpq%*%invGamma_T_lq
    Cov_T_r = invGamma_T_rp%*%Psi_T_rpq%*%invGamma_T_rq
    Bias_T_tau = (tau_T_rq[p+2,1]*BiasConst_rp[deriv+1,1] - tau_T_lq[p+2,1]*BiasConst_lp[deriv+1,1])*(h^(p+1-deriv)/factorial(p+1))
    V_T_l_cl = factorial(deriv)^2*V_T_lp[deriv+1,deriv+1]
    V_T_r_cl = factorial(deriv)^2*V_T_rp[deriv+1,deriv+1]  
    V_T_cl   = V_T_l_cl + V_T_r_cl 
    V_T_l_rb = V_T_l_cl + factorial(p+1)^2*V_T_lq[p+2,p+2]*(BiasConst_lp[deriv+1]*h^(p+1-deriv)/factorial(p+1))^2 - 2*factorial(deriv)*factorial(p+1)*Cov_T_l[deriv+1,p+2]*(BiasConst_lp[deriv+1]*h^(p+1-deriv)/factorial(p+1))
    V_T_r_rb = V_T_r_cl + factorial(p+1)^2*V_T_rq[p+2,p+2]*(BiasConst_rp[deriv+1]*h^(p+1-deriv)/factorial(p+1))^2 - 2*factorial(deriv)*factorial(p+1)*Cov_T_r[deriv+1,p+2]*(BiasConst_rp[deriv+1]*h^(p+1-deriv)/factorial(p+1))
    V_T_rb = V_T_l_rb + V_T_r_rb
    
    tau_T = c(tau_T_rp[deriv+1,1] - tau_T_lp[deriv+1,1], tau_T_rp[deriv+1,1] - tau_T_lp[deriv+1,1] - Bias_T_tau,tau_T_rp[deriv+1,1] - tau_T_lp[deriv+1,1] - Bias_T_tau)
    se_T = c(sqrt(V_T_cl), sqrt(V_T_rb), sqrt(V_T_rb))
    t_T = tau_T/se_T
    pv_T = 2*pnorm(-abs(t_T))
    ci_T=matrix(NA,nrow=3,ncol=2)
    ci_T[1,] = c(tau_T[1] - quant*se_T[1], tau_T[1] + quant*se_T[1])
    ci_T[2,] = c(tau_T[2] - quant*se_T[2], tau_T[2] + quant*se_T[2])
    ci_T[3,] = c(tau_T[3] - quant*se_T[3], tau_T[3] + quant*se_T[3])
        
    ############  Second Stage
    Bias_F_tau = (1/tau_T[1])*Bias_tau - (tau_Y[1]/tau_T[1]^2)*Bias_T_tau 
    tau_F = c(tau_Y[1]/tau_T[1], tau_Y[1]/tau_T[1] - Bias_F_tau, tau_Y[1]/tau_T[1] - Bias_F_tau)
    sigmah_TY_l = c(rdvce(X=Xh_l, y=Yh_l, frd=Th_l, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmah_TY_r = c(rdvce(X=Xh_r, y=Yh_r, frd=Th_r, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmab_TY_l = c(rdvce(X=Xb_l, y=Yb_l, frd=Tb_l, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    sigmab_TY_r = c(rdvce(X=Xb_r, y=Yb_r, frd=Tb_r, p=p, h=h, matches=matches, vce=vce, kernel=kernel))
    out.lp = qrreg(Xh_lp, Th_l, whh_l, sigmah_TY_l); out.rp = qrreg(Xh_rp, Th_r, whh_r, sigmah_TY_r)
    out.lq = qrreg(Xb_lq, Tb_l, wbb_l, sigmab_TY_l); out.rq = qrreg(Xb_rq, Tb_r, wbb_r, sigmab_TY_r)
    V_TY_lp   = out.lp$Sigma.hat; V_TY_rp   = out.rp$Sigma.hat
    V_TY_lq   = out.lq$Sigma.hat; V_TY_rq   = out.rq$Sigma.hat
    invGamma_TY_lp = out.lp$X.M.X_inv; invGamma_TY_rp = out.rp$X.M.X_inv
    invGamma_TY_lq = out.lq$X.M.X_inv; invGamma_TY_rq = out.rq$X.M.X_inv
    
    if (b>=h){
      Psi_TY_lpq = t(Xb_lp*c(whb_l*sigmab_TY_l*wbb_l))%*%Xb_lq
      Psi_TY_rpq = t(Xb_rp*c(whb_r*sigmab_TY_r*wbb_r))%*%Xb_rq
    }
    else {
      Psi_TY_lpq = t(Xh_lp*c(whh_l*sigmah_TY_l*wbh_l))%*%Xh_lq
      Psi_TY_rpq = t(Xh_rp*c(whh_r*sigmah_TY_r*wbh_r))%*%Xh_rq
    } 
    
    Cov_TY_l = invGamma_TY_lp%*%Psi_TY_lpq%*%invGamma_TY_lq
    Cov_TY_r = invGamma_TY_rp%*%Psi_TY_rpq%*%invGamma_TY_rq
    V_TY_cl   = factorial(deriv)^2*(V_TY_lp[deriv+1,deriv+1] + V_TY_rp[deriv+1,deriv+1])
    V_F_cl = (1/tau_T[1]^2)*V_cl + (tau_Y[1]^2/tau_T[1]^4)*V_T_cl -(2*tau_Y[1]/tau_T[1]^3)*V_TY_cl
    C_F_l_pq = factorial(deriv)*factorial(p+1)*((1/tau_T[1]^2)*Cov_l - (2*tau_Y[1]/tau_T[1]^3)*Cov_TY_l + (tau_Y[1]^2/tau_T[1]^4)*Cov_T_l)
    C_F_r_pq = factorial(deriv)*factorial(p+1)*((1/tau_T[1]^2)*Cov_r - (2*tau_Y[1]/tau_T[1]^3)*Cov_TY_r + (tau_Y[1]^2/tau_T[1]^4)*Cov_T_r)
    V_F_rb_t2 = -2*(C_F_l_pq[deriv+1,p+2]*BiasConst_lp[deriv+1] + C_F_r_pq[deriv+1,p+2]*BiasConst_rp[deriv+1])*h^(p+1-deriv)/factorial(p+1)
    C_YY_lq = factorial(p+1)^2*V_lq[p+2,p+2]
    C_YY_rq = factorial(p+1)^2*V_rq[p+2,p+2]
    C_TT_lq = factorial(p+1)^2*V_T_lq[p+2,p+2]
    C_TT_rq = factorial(p+1)^2*V_T_rq[p+2,p+2]
    C_TY_lq = factorial(p+1)^2*V_TY_lq[p+2,p+2]
    C_TY_rq = factorial(p+1)^2*V_TY_rq[p+2,p+2]
    V_F_cl_lq = (1/tau_T[1]^2)*C_YY_lq + (tau_Y[1]^2/tau_T[1]^4)*C_TT_lq -(2*tau_Y[1]/tau_T[1]^3)*C_TY_lq
    V_F_cl_rq = (1/tau_T[1]^2)*C_YY_rq + (tau_Y[1]^2/tau_T[1]^4)*C_TT_rq -(2*tau_Y[1]/tau_T[1]^3)*C_TY_rq
    V_F_rb_t3 = (V_F_cl_lq*BiasConst_lp[deriv+1]^2 + V_F_cl_rq*BiasConst_rp[deriv+1]^2)*(h^(p+1-deriv)/factorial(p+1))^2
    V_F_rb = V_F_cl + V_F_rb_t2 + V_F_rb_t3

    se_F = c(sqrt(V_F_cl),sqrt(V_F_cl), sqrt(V_F_rb))
    t_F =  tau_F/se_F
    pv_F = 2*pnorm(-abs(t_F))
    ci_F=matrix(NA,nrow=3,ncol=2)
    ci_F[1,] = c(tau_F[1] - quant*se_F[1], tau_F[1] + quant*se_F[1])
    ci_F[2,] = c(tau_F[2] - quant*se_F[2], tau_F[2] + quant*se_F[2])
    ci_F[3,] = c(tau_F[3] - quant*se_F[3], tau_F[3] + quant*se_F[3])
  }

  #print("Estimation Completed.") 

  if (type == "sharp") {  
    coef=matrix(tau_Y,3,1)
    se  =matrix(se_Y,3,1)
    z   =matrix(t_Y,3,1)
    pv  =matrix(pv_Y,3,1)
    ci=ci_Y
  }
  else if (type == "fuzzy") {
    coef=matrix(tau_F,3,1)
    se  =matrix(se_F,3,1)
    z   =matrix(t_F,3,1)
    pv  =matrix(pv_F,3,1)
    ci=ci_F
  }

  bws=matrix(c(h,b),1,2)
  rownames(coef)=rownames(se)=rownames(se)=rownames(z)=rownames(pv)=c("Conventional","Bias-Corrected","Robust")
  colnames(coef)="Coeff"
  colnames(se)="Std. Err."
  colnames(z)="z"
  colnames(pv)="P>|z|"
  colnames(bws)=c("h","b")
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("CI Lower","CI Upper")
    
  tabl1.str=matrix(NA,4,1)
  rownames(tabl1.str)=c("Number of Obs", "NN Matches", "BW Type", "Kernel Type")
  dimnames(tabl1.str) <-list(c("Number of Obs", "NN Matches", "BW Type", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=N
  tabl1.str[2,1]=matches
  tabl1.str[3,1]=bwselect
  tabl1.str[4,1]=kernel_type
  
  tabl2.str=matrix(NA,6,2)
  colnames(tabl2.str)=c("Left","Right")
  rownames(tabl2.str)=c("Number of Obs","Order Loc Poly (p)","Order Bias (q)","BW Loc Poly (h)","BW Bias (b)","rho (h/b)")
  tabl2.str[1,]=formatC(c(Nh_l,Nh_r),digits=0, format="f")
  tabl2.str[2,]=formatC(c(p,p),digits=0, format="f")
  tabl2.str[3,]=formatC(c(q,q),digits=0, format="f")
  tabl2.str[4,]=formatC(c(h,h),digits=4, format="f")
  tabl2.str[5,]=formatC(c(b,b),digits=4, format="f")
  tabl2.str[6,]=formatC(c(h/b,h/b),digits=4, format="f")
  
  tabl3.str=matrix("",2,6)
  colnames(tabl3.str)=c("Coef","Std. Err.","z","P>|z|","CI Lower","CI Upper")
  rownames(tabl3.str)=c("Conventional", "Robust")
  tabl3.str[1,1]  =formatC(coef[1],digits=4, format="f")
  tabl3.str[1,2]  =formatC(se[1],  digits=4, format="f")
  tabl3.str[1,3]  =formatC(z[1],   digits=4, format="f")
  tabl3.str[1,4]  =formatC(pv[1],  digits=4, format="f")
  tabl3.str[1,5:6]=formatC(ci[1,], digits=4, format="f")
  tabl3.str[2,4]  =formatC(pv[3],  digits=4, format="f")
  tabl3.str[2,5:6]=formatC(ci[3,] ,digits=4, format="f")
  
  if (all==TRUE){
    tabl3.str=formatC(cbind(coef,se,z,pv,ci),digits=4, format="f")                   
    colnames(tabl3.str)=c("Coef","Std. Err.","z","P>|z|","CI Lower","CI Upper")
  }

  out=list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,tabl3.str=tabl3.str,coef=coef,bws=bws,se=se,z=z,pv=pv,ci=ci,p=p,q=q,h=h,b=b,rho=rho,N=N,N_l=Nh_l,N_r=Nh_r)
  out$call <- match.call()
  class(out) <- "rdrobust"
  return(out)
}

#rdrobust <- function(y,x, ...) UseMethod("rdrobust")

#rdrobust.default <- function(y,x,  ...){
#  est <- rdrobustEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "rdrobust"
#  est
#}

print.rdrobust <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nSummary:\n")
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F)
  cat("\nEstimates:\n")
  print(x$tabl3.str,quote=F)
}

summary.rdrobust <- function(object,...) {
  TAB <- cbind(Estimate    =object$coef,
               "Std. Error"=object$se,
               "z"         =object$z,
               "Pr(>|z|)"  =object$pv,
               "95% CI"    =object$ci)
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.rdrobust"
  res
}

#print.summary.rdrobust <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coef)
#  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
#}