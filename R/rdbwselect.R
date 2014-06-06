### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014

rdbwselect = function(y, x, data, subset = NULL, c=0, p=1, q=2, deriv=0, rho=NULL, kernel="tri", bwselect="CCT", scaleregul=1, delta=0.5, cvgrid_min=NULL, cvgrid_max=NULL, cvgrid_length=NULL, cvplot=FALSE, vce="nn", matches=3, all=FALSE, precalc=TRUE, model = FALSE, frame = FALSE){
  
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
  
  b_calc = 0
  
  if (is.null(rho)){
    b_calc = 1
    rho = 1
  }
  
  X_l = x[x<c];    X_r = x[x>=c]
  Y_l = y[x<c];    Y_r = y[x>=c]
  N_l = length(X_l);   N_r = length(X_r)
  x_min=min(x);  x_max=max(x)
  N = N_r + N_l
  m = matches + 1
  
  if (precalc=="TRUE"){
    
    #if (deriv==0 & p==0){
    #  p = 1
    #}
    
    if (deriv>0 & p==0){
      bwselect = "CCT"
      p = deriv+1
    }
    
    if (q==0) {
      q = p+1
    }
    
    precct=0
    exit=0
    #################  ERRORS
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
      print("p, q, deriv and matches should be positive integers")
      exit = 1
    }
    
    if (p>=q & q>0){
      print("p should be set higher than q")
      exit = 1
    }
    
    if (deriv>=p & deriv>0 ){
      print("deriv should be set higher than p")
      exit = 1
    }
    
    p_round = round(p)/p;    q_round = round(q)/q;    d_round = round(deriv+1)/(deriv+1);    m_round = round(matches)/matches
        
    if (p_round!=1 | q_round!=1 | d_round!=1 | m_round!=1 ){
      print("p,q,deriv and matches should be integer numbers")
      exit = 1
    }
    
    if (delta>1 | delta<=0){
      print("delta should be set between 0 and 1")
      exit = 1
    }
    
    if (rho>1 | rho<0){
      print("rho should be set between 0 and 1")
      exit = 1
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
  
  p1 = p+1;  p2 = p+2;  q1 = q+1;  q2 = q+2;  q3 = q+3  
  h_CCT=b_CCT=h_IK=b_IK=h_CV=NA
  N = length(x)
  
  ct1 = bwconst(p,deriv,kernel)
  C1_h = ct1[1];  C2_h = ct1[2]
  ct2 = bwconst(q,q,kernel)
  C1_b = ct2[1];  C2_b = ct2[2]
  ct3 = bwconst(q1,q1,kernel)
  C1_q = ct3[1];  C2_q = ct3[2]
  
  #***********************************************************************
  #**************************** CCT Approach
  #***********************************************************************
  if (bwselect=="CCT" | all==TRUE) {
    #print("Computing CCT Bandwidth Selector.")
    #### Step 1: q_CCT
    h_pilot_cct = 2.576*min(sd(x),IQR(x))*N^(-1/5)
    X_lq2 = matrix(c((X_l-c)^0, poly(X_l-c,degree=(q3-1),raw=T)),length(X_l),q3)
    X_rq2 = matrix(c((X_r-c)^0, poly(X_r-c,degree=(q3-1),raw=T)),length(X_r),q3)
    m4_l_pilot_cct = qr.coef(qr(X_lq2, tol = 1e-10), Y_l)[q3]
    m4_r_pilot_cct = qr.coef(qr(X_rq2, tol = 1e-10), Y_r)[q3]
    w_hpilot_l = kweight(X_l, c, h_pilot_cct, kernel)
    w_hpilot_r = kweight(X_r, c, h_pilot_cct, kernel)
    Yh_l = Y_l[w_hpilot_l>0]; Xh_l = X_l[w_hpilot_l>0]; w_hpilot_l = w_hpilot_l[w_hpilot_l>0]
    Yh_r = Y_r[w_hpilot_r>0]; Xh_r = X_r[w_hpilot_r>0]; w_hpilot_r = w_hpilot_r[w_hpilot_r>0]
    sigma_l_pilot = c(rdvce(X=Xh_l, y=Yh_l, p=p, h=h_pilot_cct, matches=matches, vce=vce, kernel=kernel))
    sigma_r_pilot = c(rdvce(X=Xh_r, y=Yh_r, p=p, h=h_pilot_cct, matches=matches, vce=vce, kernel=kernel))
    X_lq1  = matrix(c((Xh_l-c)^0,poly(Xh_l-c,degree=q+1,raw=T)),length(Xh_l),q+2)
    X_rq1  = matrix(c((Xh_r-c)^0,poly(Xh_r-c,degree=q+1,raw=T)),length(Xh_r),q+2)
    X_lq  =  matrix(c((Xh_l-c)^0,poly(Xh_l-c,degree=q, raw=T)), length(Xh_l),q+1)
    X_rq  =  matrix(c((Xh_r-c)^0,poly(Xh_r-c,degree=q, raw=T)), length(Xh_r),q+1)
    X_lp  =  matrix(c((Xh_l-c)^0,poly(Xh_l-c,degree=p, raw=T)), length(Xh_l),p+1)
    X_rp  =  matrix(c((Xh_r-c)^0,poly(Xh_r-c,degree=p, raw=T)), length(Xh_r),p+1)
    out.lq1=qr.reg(x=X_lq1, y=Yh_l, w=w_hpilot_l, s2=sigma_l_pilot) 
    out.rq1=qr.reg(x=X_rq1, y=Yh_r, w=w_hpilot_r, s2=sigma_r_pilot) 
    V_m3_hpilot_cct = out.lq1$Sigma.hat[q+2,q+2]+ out.rq1$Sigma.hat[q+2,q+2]
    out.lq=qr.reg(x=X_lq, y=Yh_l, w=w_hpilot_l, s2=sigma_l_pilot) 
    out.rq=qr.reg(x=X_rq, y=Yh_r, w=w_hpilot_r, s2=sigma_r_pilot) 
    V_m2_hpilot_cct   = out.lq$Sigma.hat[q+1,q+1] + out.rq$Sigma.hat[q+1,q+1] 
    out.lp=qr.reg(x=X_lp, y=Yh_l, w=w_hpilot_l, s2=sigma_l_pilot) 
    out.rp=qr.reg(x=X_rp, y=Yh_r, w=w_hpilot_r, s2=sigma_r_pilot) 
    V_m0_hpilot_cct   = out.lp$Sigma.hat[deriv+1,deriv+1]+ out.rp$Sigma.hat[deriv+1,deriv+1]
    N_q_cct=(2*q+3)*N*h_pilot_cct^(2*q+3)*V_m3_hpilot_cct
    D_q_cct=2*(C1_q*(m4_r_pilot_cct-(-1)^(deriv+q)*m4_l_pilot_cct))^2
    q_CCT=(N_q_cct/(N*D_q_cct))^(1/(2*q+5))
    
    ### Step 2: b_CCT
    w_q_l=kweight(X_l,c,q_CCT,kernel);w_q_r=kweight(X_r,c,q_CCT,kernel)
    Yh_l  = Y_l[w_q_l>0];  Yh_r  = Y_r[w_q_r>0];    Xh_l  = X_l[w_q_l>0];  Xh_r  = X_r[w_q_r>0]
    w_q_l = w_q_l[w_q_l>0]; w_q_r = w_q_r[w_q_r>0]
    sigma_l_pilot = c(rdvce(X=Xh_l, y=Yh_l, p=p, h=h_pilot_cct, matches=matches, vce=vce, kernel=kernel))
    sigma_r_pilot = c(rdvce(X=Xh_r, y=Yh_r, p=p, h=h_pilot_cct, matches=matches, vce=vce, kernel=kernel))
    X_lq1  = matrix(c((Xh_l-c)^0,poly(Xh_l-c,degree=q+1,raw=T)),length(Xh_l),q+2)
    X_rq1  = matrix(c((Xh_r-c)^0,poly(Xh_r-c,degree=q+1,raw=T)),length(Xh_r),q+2)
    out.lq1=qr.reg(x=X_lq1, y=Yh_l, w=w_q_l, s2=sigma_l_pilot) 
    out.rq1=qr.reg(x=X_rq1, y=Yh_r, w=w_q_r, s2=sigma_r_pilot) 
    V_m3_q_cct   = out.lq1$Sigma.hat[q+2,q+2]+ out.rq1$Sigma.hat[q+2,q+2]
    m3_l_cct = out.lq1$beta[q2]; m3_r_cct = out.rq1$beta[q2] 
    N_b_cct=  (2*p+3)*N*h_pilot_cct^(2*p+3)*V_m2_hpilot_cct
    D_b_cct=  2*(q-p)*(C1_b*(m3_r_cct - (-1)^(deriv+q+1)*m3_l_cct))^2
    R_b_cct=  scaleregul*2*(q-p)*(C1_b)^2*3*V_m3_q_cct
    b_CCT= (N_b_cct / (N*(D_b_cct+R_b_cct)))^(1/(2*q+3))
  
    
    ### Step 3: h_CCT
    w_b_l = kweight(X_l,c,b_CCT,kernel);w_b_r = kweight(X_r,c,b_CCT,kernel)
    Yh_l=Y_l[w_b_l>0];  Yh_r=Y_r[w_b_r>0];  Xh_l=X_l[w_b_l>0];  Xh_r=X_r[w_b_r>0];w_b_l=w_b_l[w_b_l>0]; w_b_r=w_b_r[w_b_r>0]
    sigma_l_pilot = c(rdvce(X=Xh_l, y=Yh_l, p=p, h=h_pilot_cct, matches=matches, vce=vce, kernel=kernel))
    sigma_r_pilot = c(rdvce(X=Xh_r, y=Yh_r, p=p, h=h_pilot_cct, matches=matches, vce=vce, kernel=kernel))
    X_lq  = matrix(c((Xh_l-c)^0,poly(Xh_l-c,degree=q,raw=T)),length(Xh_l),q+1)
    X_rq  = matrix(c((Xh_r-c)^0,poly(Xh_r-c,degree=q,raw=T)),length(Xh_r),q+1)
    out.lq=qr.reg(x=X_lq, y=Yh_l, w=w_b_l, s2=sigma_l_pilot) 
    out.rq=qr.reg(x=X_rq, y=Yh_r, w=w_b_r, s2=sigma_r_pilot) 
    V_m2_b_cct=out.lq$Sigma.hat[p2,p2] + out.rq$Sigma.hat[p2,p2] 
    m2_l_cct=out.lq$beta[p2];  m2_r_cct=out.rq$beta[p2] 
    N_h_cct = (2*deriv+1)*N*h_pilot_cct^(2*deriv+1)*V_m0_hpilot_cct
    D_h_cct = 2*(p+1-deriv)*(C1_h*(m2_r_cct - (-1)^(deriv+p+1)*m2_l_cct))^2
    R_h_cct = scaleregul*2*(p+1-deriv)*(C1_h)^2*3*V_m2_b_cct
    h_CCT = (N_h_cct / (N*(D_h_cct+R_h_cct)))^(1/(2*p+3))
    
    if (b_calc==0) {
      b_CCT = h_CCT/rho
    }
    results = matrix(NA,1,2)
    colnames(results)=c("h","b")
    rownames(results)=""
    results[1,]=c(h_CCT,b_CCT)
  }
  
  #***************************************************************************************************
  #******************** IK
  #**************************************************************************************************
  if (bwselect=="IK" | all==TRUE) {

    X_lq2 = matrix(c((X_l-c)^0, poly(X_l-c,degree=(q3-1),raw=T)),length(X_l),q3)
    X_rq2 = matrix(c((X_r-c)^0, poly(X_r-c,degree=(q3-1),raw=T)),length(X_r),q3)
    X_lq1 = X_lq2[,1:q2];  X_rq1 = X_rq2[,1:q2]
    X_lq  = X_lq2[,1:q1];  X_rq  = X_rq2[,1:q1]
    X_lp  = X_lq2[,1:p1];  X_rp  = X_rq2[,1:p1]
    
    #print("Computing IK Bandwidth Selector.")
    h_pilot_IK = 1.84*sd(x)*N^(-1/5)
    temp=X_l[X_l>=c-h_pilot_IK]
    n_l_h1 = length(temp)
    temp=X_r[X_r<=c+h_pilot_IK]
    n_r_h1 = length(temp)
    f0_pilot=(n_r_h1+n_l_h1)/(2*N*h_pilot_IK)
    temp=Y_l[X_l>=c-h_pilot_IK]
    s2_l_pilot = var(temp)
    if (s2_l_pilot==0){
      s2_l_pilot=var(Y_l[X_l>=c-2*h_pilot_IK])
    }
    temp=Y_r[X_r<=c+h_pilot_IK]
    s2_r_pilot = var(temp)
    if (s2_r_pilot==0){
      s2_r_pilot=var(Y_r[X_r<=c+2*h_pilot_IK])
    }
    
    V_IK_pilot = (s2_r_pilot+s2_l_pilot)/f0_pilot
    Vm0_pilot_IK = C2_h*V_IK_pilot
    Vm2_pilot_IK = C2_b*V_IK_pilot
    Vm3_pilot_IK = C2_q*V_IK_pilot
    
    x_IK_sort_l = X_l[X_l>=median(X_l)]; y_IK_sort_l = Y_l[X_l>=median(X_l)]
    x_IK_sort_r = X_r[X_r<=median(X_r)]; y_IK_sort_r = Y_r[X_r<=median(X_r)]
    x_IK_sort = c(x_IK_sort_r,x_IK_sort_l); y_IK_sort = c(y_IK_sort_r,y_IK_sort_l)
    sample_IK = length(x_IK_sort)
    X_IK_q2 = matrix(c((x_IK_sort-c)^0, poly(x_IK_sort-c,degree=(q3-1),raw=T)),sample_IK,q3)
    X_IK_q1 = X_IK_q2[,1:q2]
    
    ### First Stage
    m4_l_pilot_IK = m4_r_pilot_IK = qr.coef(qr(X_IK_q2, tol = 1e-10), y_IK_sort)[q+3]
    N_q_r_pilot_IK = (2*q1+1)*C2_q*(s2_r_pilot/f0_pilot)
    N_q_l_pilot_IK = (2*q1+1)*C2_q*(s2_l_pilot/f0_pilot)
    D_q_r_pilot_IK = 2*(q1+1-q1)*(C1_q*m4_r_pilot_IK)^2
    D_q_l_pilot_IK = 2*(q1+1-q1)*(C1_q*m4_l_pilot_IK)^2
    h3_r_pilot_IK = (N_q_r_pilot_IK / (N_r*D_q_r_pilot_IK))^(1/(2*q+5))
    h3_l_pilot_IK = (N_q_l_pilot_IK / (N_l*D_q_l_pilot_IK))^(1/(2*q+5))
    temp=X_l[X_l>=c-h3_l_pilot_IK];    n_l_h3 = length(temp)
    temp=X_r[X_r<=c+h3_r_pilot_IK];    n_r_h3 = length(temp)
    w_h3_l = kweight(X_l,c,h3_l_pilot_IK,"uni")
    w_h3_r = kweight(X_r,c,h3_r_pilot_IK,"uni")
    m3_l_IK=qr.reg(x=X_lq1, y=Y_l, w=w_h3_l,var.comp=FALSE)$beta[q2]
    m3_r_IK=qr.reg(x=X_rq1, y=Y_r, w=w_h3_r,var.comp=FALSE)$beta[q2]
    N_b_IK = (2*p1+1)*Vm2_pilot_IK
    D_b_IK = 2*(q-p)*(C1_b*(m3_r_IK - (-1)^(deriv+q+1)*m3_l_IK))^2
    temp = regconst(q1,1);    con = temp[q2,q2]
    r_l_b = (con*s2_l_pilot)/(n_l_h3*h3_l_pilot_IK^(2*q1))
    r_r_b = (con*s2_r_pilot)/(n_r_h3*h3_r_pilot_IK^(2*q1))
    Vm3_pilot_IK1 = (r_l_b + r_r_b)
    R_b_IK = scaleregul*2*(q-p)*(C1_b)^2*3*Vm3_pilot_IK1
    b_IK   = (N_b_IK / (N*(D_b_IK+R_b_IK)))^(1/(2*q+3))
    
    ### Second Stage
    m3_l_pilot_IK = m3_r_pilot_IK = qr.coef(qr(X_IK_q1, tol = 1e-10), y_IK_sort)[q2]
    N_b_r_pilot_IK = (2*p1+1)*C2_b*(s2_r_pilot/f0_pilot)
    N_b_l_pilot_IK = (2*p1+1)*C2_b*(s2_l_pilot/f0_pilot)
    D_b_r_pilot_IK = 2*(q-p)*(C1_b*m3_r_pilot_IK)^2
    D_b_l_pilot_IK = 2*(q-p)*(C1_b*m3_l_pilot_IK)^2
    h2_l_pilot_IK  = (N_b_l_pilot_IK / (N_l*D_b_l_pilot_IK))^(1/(2*q+3))
    h2_r_pilot_IK  = (N_b_r_pilot_IK / (N_r*D_b_r_pilot_IK))^(1/(2*q+3))
    w_h2_l = kweight(X_l,c, h2_l_pilot_IK, "uni")
    w_h2_r = kweight(X_r,c, h2_r_pilot_IK, "uni")
    m2_l_IK=qr.reg(x=X_lq, y=Y_l, w=w_h2_l,var.comp=FALSE)$beta[p2]
    m2_r_IK=qr.reg(x=X_rq, y=Y_r, w=w_h2_r,var.comp=FALSE)$beta[p2]
    temp=X_l[X_l>=c-h2_l_pilot_IK]; n_l_h2 = length(temp)
    temp=X_r[X_r<=c+h2_r_pilot_IK]; n_r_h2 = length(temp)
    temp = regconst(p1,1);  con = temp[p2,p2]
    r_l_h = (con*s2_l_pilot)/(n_l_h2*h2_l_pilot_IK^(2*p1))
    r_r_h = (con*s2_r_pilot)/(n_r_h2*h2_r_pilot_IK^(2*p1))
    N_h_IK = (2*deriv+1)*Vm0_pilot_IK
    D_h_IK = 2*(p+1-deriv)*(C1_h*(m2_r_IK - (-1)^(deriv+p+1)*m2_l_IK))^2
    R_h_IK = scaleregul*2*(p+1-deriv)*(C1_h)^2*3*(r_l_h + r_r_h)
    h_IK  = (N_h_IK / (N*(D_h_IK+R_h_IK)))^(1/(2*p+3))
    
    #*** DJMC
    D_b_DM  = sqrt(m3_r_IK^2 + m3_l_IK^2); b_DM = (N_b_IK / (N*(D_b_DM)^2))^(1/(2*q+3))
    D_h_DM  = sqrt(m2_r_IK^2 + m2_l_IK^2); h_DM = (N_h_IK / (N*(D_h_DM)^2))^(1/(2*p+3))
    
    if (b_calc==0) {
      b_IK = h_IK/rho
    }
    
    results = matrix(NA,1,2)
    colnames(results)=c("h","b")
    #rownames(results)=c("IK")
    results[1,]=c(h_IK,b_IK)
    
  }
  
  #*********************************************************************
  #********************************** C-V  *****************************
  #*********************************************************************
  if (bwselect=="CV" | all==TRUE) {
    
    #print("Computing CV Bandwidth Selector.")
    norep_l=unique(data.frame(X_l,Y_l))
    norep_r=unique(data.frame(X_r,Y_r))
    X_nr_l = norep_l[,1]; Y_nr_l = norep_l[,2]
    X_nr_r = norep_r[,1]; Y_nr_r = norep_r[,2]
    N_nr_l = length(X_nr_l);N_nr_r = length(X_nr_r)
    
    v_CV_l = order(X_nr_l,decreasing=FALSE)
    x_sort_l = X_nr_l[v_CV_l]
    y_sort_l = Y_nr_l[v_CV_l]
    v_CV_r = order(X_nr_r,decreasing=TRUE)
    x_sort_r = X_nr_r[v_CV_r]
    y_sort_r = Y_nr_r[v_CV_r]
    h_CV_min = 0
    
    if (N_nr_r>20 & N_nr_l>20){
      h_CV_min  = min(c(abs(x_sort_r[N_nr_r]-x_sort_r[N_nr_r-20]),abs(x_sort_l[N_nr_l]-x_sort_l[N_nr_l-20])))
    }
    
    h_CV_max  = min(c(abs(x_sort_r[1]-x_sort_r[N_nr_r]),abs(x_sort_l[1]-x_sort_l[N_nr_l])))
    h_CV_jump = min(c(abs(x_sort_r[1]-x_sort_r[N_nr_r])/10,abs(x_sort_l[1]-x_sort_l[N_nr_l]))/10)
    
    if (is.null(cvgrid_min)) {
        cvgrid_min = h_CV_min
    }

    if (is.null(cvgrid_max)) {
        cvgrid_max = h_CV_max
    }
    
    if (is.null(cvgrid_length)) {
      cvgrid_length = abs(cvgrid_max-cvgrid_min)/20
    }
    
    if (cvgrid_min>=cvgrid_max){
      cvgrid_min = 0
    }
    
    h_CV_seq  = seq(cvgrid_min, cvgrid_max, cvgrid_length)
    s_CV = length(h_CV_seq)
    CV_l = CV_r = matrix(0,1,s_CV)    

    n_CV_l = round(delta*N_nr_l)-3
    n_CV_r = round(delta*N_nr_r)-3
    
    # Set quantile sample
    for (v in 1:s_CV) {
      for (k in 0:n_CV_l) {
        ind_l = N_nr_l-k-1
        x_CV_sort_l = x_sort_l[1:ind_l] 
        y_CV_sort_l = y_sort_l[1:ind_l] 
        w_CV_sort_l = kweight(x_CV_sort_l,x_sort_l[ind_l+1],h_CV_seq[v],kernel)
        x_CV_l = x_CV_sort_l[w_CV_sort_l>0]
        y_CV_l = y_CV_sort_l[w_CV_sort_l>0]
        w_CV_l = w_CV_sort_l[w_CV_sort_l>0]
        XX_CV_l = matrix(c((x_CV_l-x_sort_l[ind_l+1])^0, poly(x_CV_l-x_sort_l[ind_l+1],degree=p,raw=T)),length(w_CV_l),p+1)
        y_CV_hat_l = qr.coef(qr(XX_CV_l*sqrt(w_CV_l), tol = 1e-10), sqrt(w_CV_l)*y_CV_l)[1]
        mse_CV_l = (y_sort_l[ind_l+1] - y_CV_hat_l)^2
        CV_l[v] = CV_l[v] + mse_CV_l
      }
      for (k in 0:n_CV_r) {
        ind_r = N_nr_r-k-1
        x_CV_sort_r = x_sort_r[1:ind_r] 
        y_CV_sort_r = y_sort_r[1:ind_r] 
        w_CV_sort_r = kweight(x_CV_sort_r,x_sort_r[ind_r+1],h_CV_seq[v],kernel)
        x_CV_r = x_CV_sort_r[w_CV_sort_r>0]
        y_CV_r = y_CV_sort_r[w_CV_sort_r>0]
        w_CV_r = w_CV_sort_r[w_CV_sort_r>0]
        XX_CV_r = matrix(c((x_CV_r - x_sort_r[ind_r+1])^0, poly(x_CV_r-x_sort_r[ind_r+1],degree=p,raw=T)),length(w_CV_r),p+1)
        y_CV_hat_r = qr.coef(qr(XX_CV_r*sqrt(w_CV_r), tol = 1e-10), sqrt(w_CV_r)*y_CV_r)[1]
        mse_CV_r = (y_sort_r[ind_r+1] - y_CV_hat_r)^2
        CV_r[v] = CV_r[v] + mse_CV_r
      }
    }
    
    CV_sum = CV_l + CV_r
    CV_sum_order = order(abs(t(CV_sum)))[1]
    h_CV = h_CV_seq[CV_sum_order] 
    h_CV = h_CV[1]
    
    if (cvplot==TRUE){
      plot(h_CV_seq,CV_sum, type="l",main="Cross-Validation Objective Function",xlab="Grid of Bandwidth (h)",ylab="Cross-Validation Objective Function")
      abline(v=h_CV)
    }
    
    #if (b_calc==0) {
    #  b_CV = h_CV/rho
    #}
    results = matrix(NA,1,2)
    colnames(results)=c("h","b")
    #rownames(results)=c("CV")
    results[1,]=c(h_CV,h_CV/rho)
  }
  #results=list(CCT=c(h_CCT,b_CCT),IK=c(h_IK,b_IK),CV=h_CV)
  
  if (all=="TRUE"){
  bwselect="All"
  results = matrix(NA,3,2)
  colnames(results)=c("h","b")
  rownames(results)=c("CCT","IK","CV")
  results[1,]=c(h_CCT,b_CCT)
  results[2,]=c(h_IK,b_IK)
  results[3,1]=h_CV
  #results[4,]=c(h_DM,b_DM)
  }
  
  tabl1.str=matrix(NA,4,1)
  dimnames(tabl1.str) <-list(c("BW Selector", "Number of Obs", "NN Matches", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=bwselect
  tabl1.str[2,1]=N
  tabl1.str[3,1]=matches
  tabl1.str[4,1]=kernel_type
  
  tabl2.str=matrix(NA,3,2)
  colnames(tabl2.str)=c("Left","Right")
  rownames(tabl2.str)=c("Number of Obs","Order Loc Poly (p)","Order Bias (q)")
  tabl2.str[1,]=formatC(c(N_l,N_r),digits=0, format="f")
  tabl2.str[2,]=formatC(c(p,p),digits=0, format="f")
  tabl2.str[3,]=formatC(c(q,q),digits=0, format="f")

  bws=results
  out = list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,bws=bws,bws,bwselect=bwselect,kernel=kernel_type,p=p,q=q)
  out$call <- match.call()
  class(out) <- "rdbwselect"
  return(out)
}


#rdbwselect <- function(y,x, ...) UseMethod("rdbwselect")

#rdbwselect.default <- function(y,x, ...){
#  est <- rdbwselectEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "rdbwselect"
#  est
#}

print.rdbwselect <- function(x,...){
  cat("Call:\n")
  print(x$call)
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F) 
  cat("\n")
  print(x$bws)  
}

summary.rdbwselect <- function(object,...) {
  TAB <- object$bws
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.rdbwselect"
  res
}

print.summary.rdbwselect <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
}