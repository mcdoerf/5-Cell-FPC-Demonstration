####-----------------
###Helper Functions
###-----------------

library(MCMCpack)  #For rdirichlet() function

data_decomp<-function(data){  ###Decomposes data.  Assumes the data comes in the five cell format,
  #as (n15, n2, n4, n6, n37).
  
  n15=data[1]
  n2=data[2]
  n4=data[3]
  n6=data[4]
  n37=data[5]
  
  return(list(n15=n15, n2=n2, n4=n4, n6=n6, n37=n37, n11=n2, n10=n4, n01=n6, nc=n2+n4+n6, n=n15+n2+n6,  Ntot=sum(data)))
}

########Estimators

# RS estimator & variance
rs_estimator <- function(data, type="Unadjusted") {
  
  counts<-data_decomp(data)  ##Decompose the data.
  
  n_rs<- counts$n2+counts$n6   ##Number of individuals that test positive in the random sample.
  
  n<-counts$n15+counts$n2+counts$n6 ##Number of individuals in the random sample.
  
  p_rs<-n_rs/n #Point Estimate of Prevalence
  
  est=(p_rs)*(counts$Ntot) #Point Estimate of Number of Diseased Individuals
  
  p_rs2<-max(n_rs, 0.5)/n #Add 0.5 if the number of diseased in the random sample is 0.
  
  var_rs<-(counts$Ntot)^2*p_rs2*(1-p_rs2)/(n)   ###Non FPC adjusted variance
  
  FPC<-min((counts$n*(counts$Ntot-counts$n))/(counts$Ntot*(counts$n-1)),1) ##Cochran's FPC
  
  var_rs_FPC<-FPC*var_rs  ##FPC adjusted variance of random sample estimator.
  
  if (type=="Unadjusted"){
    return(list(est = est, se = sqrt(var_rs)))
  }
  
  if (type=="FPC"){
    return(list(est = est, se = sqrt(var_rs_FPC)))
  }
  
}


# Chapman estimator & variance
chapman_estimator <- function(data) {
  counts<-data_decomp(data)
  
  n1.<-counts$n11+counts$n10
  n.1<-counts$n11+counts$n01
  
  nchapman<-((n1.+1)*(n.1+1))/(counts$n11+1)-1  ##Chapman's point estimator
  
  var_chap<-(counts$n11+counts$n10+1)*(counts$n11+counts$n01+1)*(counts$n10)*(counts$n01)/((counts$n11+1)^2*(counts$n11+2))  ##Variance of Chapman's estimator.
  
  list(est=nchapman, se=sqrt(var_chap))
}



# 5 cell MLE estimator & variance


five_cell_estimator<-function(data, psi, type="Unadjusted"){
  
  
  
  counts<-data_decomp(data=data)  ##Decompose into counts
  
  N_tot_star=counts$n15+counts$n6+counts$n37
  
  w=(counts$n2+counts$n4)/counts$Ntot
  
  n_rs_star=counts$n15+counts$n6
  
  
  
  
  if (counts$n6==0){
    N_MLE=counts$n2+counts$n4
  } 
  
  else{
    
    eta=(counts$n6)/(counts$n15+counts$n6)
    
    
    pi_all = w+eta*(1-w)
    
    
    ###5 cell MLE estimate
    
    N_MLE = counts$Ntot*pi_all
    
  }
  ##Non FPC Corrected standard error
  
  if (counts$n6==0 || n_rs_star==0){
    var=counts$Ntot^2*w*(1-w)/counts$Ntot
  }
  
  else{
    var_eta=(eta)*(1-eta)/(counts$Ntot*(1-w)*psi)
    var=counts$Ntot^2*((1-w)^2*var_eta+(1-eta)^2*w*(1-w)/counts$Ntot)
  }
 
  
  #FPC1 Corrected Variance
  
  if (n_rs_star==0 || n_rs_star==1){
    var_FPC1=var
    var_FPC2=var
    
  } else if(counts$n6==0){
    
    eta_new=(counts$n6+0.5)/(counts$n6+1+counts$n15) #if n6=0, get the estimate that would be obtained if we
    #used a Jeffreys prior to (p15, p2, p4, p6, p37)^{T}.
    
    var_eta=min(n_rs_star*(N_tot_star-n_rs_star)/(N_tot_star*(n_rs_star-1)) ,1)*(eta_new*(1-eta_new)/n_rs_star)
    
    var_pi_FPC1=(1-w)^2*var_eta
    
    var_FPC1=counts$Ntot^2*var_pi_FPC1  ##FPC Adjustment 1
    
    #FPC corrected standard error 2 (acknowledging variability in w)
    
    var_FPC2<-counts$Ntot^2*(var_pi_FPC1+((1-eta_new)^2)*w*(1-w)/counts$Ntot)
    
  }
  
  else{
    eta_new=(counts$n6)/(counts$n15+counts$n6)
    
    var_eta=min(n_rs_star*(N_tot_star-n_rs_star)/(N_tot_star*(n_rs_star-1)) ,1)*(eta_new*(1-eta_new)/n_rs_star)
    
    var_pi_FPC1=(1-w)^2*var_eta
    
    var_FPC1=counts$Ntot^2*var_pi_FPC1  ##FPC Adjustment 1
    
    #FPC corrected standard error 2 (acknowledging variability in w)
    
    var_FPC2<-counts$Ntot^2*(var_pi_FPC1+((1-eta_new)^2)*w*(1-w)/counts$Ntot)
  }
  
  
  
  
  if (type=="Unadjusted"){
    return(list(est=N_MLE, se=sqrt(var)))
  }
  
  if (type=="FPC1"){
    return(list(est=N_MLE, se=sqrt(var_FPC1)))
  }
  
  if (type == "FPC2"){
    return(list(est=N_MLE, se=sqrt(var_FPC2)))
  }
  
} 


#################################Confidence Intervals


####To make sure no pathological cases occur, we define safe_interval.

safe_interval <- function(lwr, upr, context = "unknown") {
  if (!is.finite(lwr) || !is.finite(upr) || lwr > upr) {
    warning("Degenerate CI in ", context, 
            ": lwr=", lwr, ", upr=", upr, call. = FALSE)
    return(c(NA_real_, NA_real_))
  } else {
    return(c(lwr, upr))
  }
}


###Generic Wald CI function

wald_ci <- function(est, se, alpha=0.05, Ntot) {
  z<-qnorm((1-alpha/2), mean=0, sd=1, lower.tail=TRUE)
  lwr <-  max(est - z * se,0)
  upr <- min(est + z * se, Ntot)
  safe_interval(lwr, upr)
}


###Generic Wald CI function, for cases where we know n15 and Ntot-n15. Truncates the interval from above at
###Ntot-n15 (since we know this is the most that can be diseased) and from below at nc (since we know at least
#this many are diseased).

wald_ci2 <- function(est, se, alpha=0.05, Ntot, n15, nc) {
  z<-qnorm((1-alpha/2), mean=0, sd=1, lower.tail=TRUE)
  lwr <- max(nc, est - z * se)         #Number of diseased has to be greater than nc.
  upr <- min(Ntot-n15, est + z * se)   #Number of diseased cannot exceed Ntot-n15
  if (lwr<upr){
    return(safe_interval(lwr, upr))
  }
  else if (lwr>=upr){
    return(c(nc, Ntot-n15))
  }
}



#Sadinle CI to accompany Chapman's estimate
sadinle_ci<-function(data, alpha=0.05){
  
  counts<-data_decomp(data)
  
  h_n<-(counts$n10+0.5)*(counts$n01+0.5)/(counts$n11+0.5)
  
  sigma_tl<-sqrt(1/(counts$n11+0.5)+1/(counts$n10+0.5)+1/(counts$n01+0.5)+(counts$n11+0.5)*(1/(counts$n10+0.5))*(1/(counts$n01+0.5)))
  
  z<-qnorm(1-(alpha/2), mean=0, sd=1, lower.tail=TRUE)
  
  lwr<-max(counts$nc, (counts$nc-0.5)+h_n*exp(-z*sigma_tl))
  
  upr<-min(counts$Ntot-counts$n15, (counts$nc-0.5)+h_n*exp(z*sigma_tl))
  
  if (upr<=lwr){   ##Sometimes (extremely rarely) the Sadinle CI gives (counts$nc-0.5)+h_n*exp(-z*sigma_tl) and
    #(counts$nc-0.5)+h_n*exp(z*sigma_tl) both less than nc.  In which case, we just set the interval to (nc, Ntot-n15).
    lwr<-counts$nc
    upr<-counts$Ntot-counts$n15
    
  }
  
  safe_interval(lwr, upr)
}



#################R Code for Rivest CI to accompany Chapman's estimator
##The code for pwar and qwar are adapted from the supplementary materials of Rivest and Yuack 2025.
##Rivest, L.P., Yauck, M., 2025. Small sample inference for two-way capture-recapture experiments. International Statistical Review 93, 62–72. doi:https://doi.org/10.1111/insr.12574.

#R code for the generalized Waring cumulative distribution function

pwar <- function(q, a, b, c, lower.tail = TRUE) {
  
  f <- function(x, q) {
    pnbinom(q, size = a, prob = x, lower.tail = lower.tail) *
      dbeta(x, shape1 = c - b - a, shape2 = b)
  }
  
  if (length(q) == 1) {
    integrate(f, lower = 0, upper = 1,
              q = q,
              rel.tol = 1e-10,
              abs.tol = 0)$value
  } else {
    sapply(q, function(qi)
      integrate(f, lower = 0, upper = 1,
                q = qi,
                rel.tol = 1e-10,
                abs.tol = 0)$value)
  }
}




#R code for the generalized Waring quantile function    
# p is a probability in (0,1)
# a b and c are positive real parameters satisfying c>a+b


qwar <- function(p, a, b, c) {
  
  if (p <= 0) return(0)
  
  # Step 1: Find an upper bound
  upper <- 1
  while (pwar(upper, a, b, c) < p) {
    upper <- upper * 2
  }
  
  lower <- 0
  
  # Step 2: Bisection search
  while (upper - lower > 1) {
    mid <- floor((lower + upper) / 2)
    if (pwar(mid, a, b, c) >= p) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  
  return(upper)
}


##Code for calculating rivest CI.  Assumes that data comes in the form
#c(n15, n2, n4, b6, b37).


rivest<-function(data, alpha=0.05){
  counts<-data_decomp(data)
  nc<-counts$nc
  n10<-counts$n10
  n01<-counts$n01
  n15<-counts$n15
  Ntot<-counts$Ntot
  
  rivest_l<-max(nc, nc+qwar(alpha/2, n10+1, n01+1, nc+3))
  rivest_u<-min(Ntot-n15,  nc+qwar(1-alpha/2, n10+1, n01+1, nc+3))
  
  if (rivest_u<=rivest_l){
    rivest_l<-counts$nc
    rivest_u<-counts$Ntot-counts$n15
    
  }
  
  return(c(rivest_l, rivest_u))
}  






#Bayes CI Generator for 5 cell MLE

bayes_ci_5<-function(data, psi, alpha=0.05, type="Unadjusted", postdraws=10000){
  
  counts<-data_decomp(data=data)
  
  n15<-counts$n15
  n2<-counts$n2
  n4<-counts$n4
  n6<-counts$n6
  n37<-counts$n37
  nc<-counts$nc
  Ntot<-counts$Ntot
  
  pstar_post<-rdirichlet(postdraws, c(n15, n2, n4, n6, n37)+0.5)
  
  nstar<-Ntot*pstar_post
  
  Nhat_star<-nstar[,2]+nstar[,3]+nstar[,4]*((nstar[,1]+nstar[,4]+nstar[,5]))/(nstar[,1]+nstar[,4])
  
  if (type=="Unadjusted"){
    
    lwr<-max(as.numeric(quantile(Nhat_star, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star, 1-alpha/2)), Ntot-n15)
    
    if (lwr>=upr){
      lwr<-nc
      upr<-Ntot-n15
    }
    
    safe_interval(lwr, upr)
    
  }
  
  else if (type=="FPC1"){
    
    se_unadj <- five_cell_estimator(data=data, psi=psi, type="Unadjusted")$se
    
    se_FPC1  <- five_cell_estimator(data=data, psi=psi, type="FPC1")$se
    
    a1 <- if (se_unadj > 0) min(se_FPC1 / se_unadj,1) else 1
    
    b1=five_cell_estimator(data=data, psi=psi, type="FPC1")$est *(1-a1)
    
    Nhat_star1<-as.vector(a1)*Nhat_star+as.vector(b1)
    
    lwr<-max(as.numeric(quantile(Nhat_star1, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star1, 1-alpha/2)), Ntot-n15)
    
    if (lwr>=upr){
      lwr<-nc
      upr<-Ntot-n15
    }
    
    safe_interval(lwr, upr)
    
  }
  
  else if (type=="FPC2"){
    
    se_unadj <- five_cell_estimator(data=data, psi=psi, type="Unadjusted")$se
    
    se_FPC2  <- five_cell_estimator(data=data, psi=psi, type="FPC2")$se
    
    a1 <- if (se_unadj > 0) min(se_FPC2 / se_unadj, 1) else 1
    
    b1=five_cell_estimator(data=data, psi=psi, type="FPC2")$est *(1-a1)
    
    Nhat_star1<-as.vector(a1)*Nhat_star+as.vector(b1)
    
    lwr<-max(as.numeric(quantile(Nhat_star1, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star1, 1-alpha/2)), Ntot-n15)
    
    if (lwr>=upr){
      lwr<-nc
      upr<-Ntot-n15
    }
    
    safe_interval(lwr, upr) 
  }
  
}


###Function for Implementing Bayesian Credible CI for Estimator Based on Random Sample Only
bayes_ci_rs<-function(data, alpha=0.05, type="FPC"){
  
  counts<-data_decomp(data)
  
  n15<-counts$n15
  
  n2<-counts$n2
  
  n4<-counts$n4
  
  n6<-counts$n6
  
  n37<-counts$n37
  
  nc<-counts$nc
  
  Ntot<-counts$Ntot
  
  n11<-counts$n11
  
  n10<-counts$n10
  
  n01<-counts$n01
  
  n<-counts$n
  
  Ntot<-counts$Ntot
  
  nrs_pos=n11+n01
  
  prs=rs_estimator(data, type="Unadjusted")$est/Ntot
  
  alphpost = nrs_pos + 0.5
  betpost = n - nrs_pos + 0.5
  
  LL_JeffreysForP = qbeta(alpha/2, alphpost, betpost)
  UL_JeffreysForP = qbeta(1-alpha/2, alphpost, betpost)
  
  if(nrs_pos == 0){ LL_JeffreysForP = 0 }
  if(nrs_pos == n){ UL_JeffreysForP = 1 }
  LL_Jeffreys = max(Ntot*LL_JeffreysForP, nc)
  UL_Jeffreys = min(Ntot*UL_JeffreysForP, Ntot-n15)
  
  if (LL_Jeffreys>=UL_Jeffreys){
    LL_Jeffreys<-nc
    UL_Jeffreys<-Ntot-n15
  }
  
  ####Adjusting Jeffry's CI for FPC
  
  a=sqrt(min((Ntot-n)*n/((Ntot)*(n-1)),1))
  b=prs*(1-a)
  LL_JeffreysForPFPC = a*LL_JeffreysForP + b
  UL_JeffreysForPFPC = a*UL_JeffreysForP + b
  if(nrs_pos == 0){ LL_JeffreysForPFPC = 0 }
  if(nrs_pos == n){ UL_JeffreysForPFPC = 1 }
  LL_JeffreysFPC = max(Ntot*LL_JeffreysForPFPC, nc)
  UL_JeffreysFPC = min(Ntot*UL_JeffreysForPFPC, Ntot-n15)
  
  if (LL_JeffreysFPC>=UL_JeffreysFPC){  ###In some instances, LL_JeffreysFPC>UL_JeffreysFPC
    LL_JeffreysFPC<-nc
    UL_JeffreysFPC<-Ntot-n15
  }
  
  if (type=="Unadjusted"){
    safe_interval(LL_Jeffreys, UL_Jeffreys)
  }
  
  if (type=="FPC"){
    safe_interval(LL_JeffreysFPC, UL_JeffreysFPC)
  }
  
  
}





##########Inference function



inference<-function(data, psi, alpha, postdraws){
  
  counts<-data_decomp(data)
  
  n2<-counts$n2
  n4<-counts$n4
  n6<-counts$n6
  n15<-counts$n15
  n37<-counts$n37
  
  n11<-counts$n11
  n10<-counts$n10
  n01<-counts$n01
  
  Ntot<-counts$Ntot
  
  nc<-counts$nc
  n<-counts$n
  
  
  ##Point Estimators
  
  five_cell_MLE<-five_cell_estimator(data=data, psi=psi, type="Unadjusted")$est
  
  rs_MLE<-rs_estimator(data=data, type="Unadjusted")$est
  
  chapman<-chapman_estimator(data=data)$est
  
  ##Standard Errors
  
  five_cell_MLE_se_unadj<-five_cell_estimator(data=data, psi=psi, type="Unadjusted")$se
  
  five_cell_MLE_se_FPC1<-five_cell_estimator(data=data, psi=psi, type="FPC1")$se
  
  five_cell_MLE_se_FPC2<-five_cell_estimator(data=data, psi=psi, type="FPC2")$se
  
  rs_MLE_se_unadj<-rs_estimator(data=data, type="Unadjusted")$se
  
  rs_MLE_se_FPC<-rs_estimator(data=data, type="FPC")$se
  
  chapman_se<-chapman_estimator(data=data)$se
  
  #Confidence Intervals
  
  #5 cell MLE
  wald_five_cell_MLE_unadj<-wald_ci2(five_cell_MLE, five_cell_MLE_se_unadj, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_five_cell_MLE_FPC1<-wald_ci2(five_cell_MLE, five_cell_MLE_se_FPC1, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_five_cell_MLE_FPC2<-wald_ci2(five_cell_MLE, five_cell_MLE_se_FPC2, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  bayes_five_cell_MLE_unadj<-bayes_ci_5(data=data, psi=psi, alpha=alpha, type="Unadjusted", postdraws=postdraws)
  
  bayes_five_cell_MLE_FPC1<-bayes_ci_5(data=data, psi=psi, alpha=alpha, type="FPC1", postdraws=postdraws)
  
  bayes_five_cell_MLE_FPC2<-bayes_ci_5(data=data, psi=psi, alpha=alpha, type="FPC2", postdraws=postdraws)
  
  #RS Estimator
  wald_rs_unadj<-wald_ci2(rs_MLE, rs_MLE_se_unadj, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_rs_FPC<-wald_ci2(rs_MLE, rs_MLE_se_FPC, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  bayes_rs<-bayes_ci_rs(data=data, alpha=alpha, type="FPC")
  
  #Chapman Estimator
  chapman_ci<-sadinle_ci(data=data, alpha=alpha)
  
  rivest_ci<-rivest(data=data, alpha=alpha)
  
  
  list(
    estimates = list(five_cell_MLE = five_cell_MLE, rs_MLE=rs_MLE,
                     chapman=chapman),
    se = list(
      
      five_cell_MLE=list(five_cell_MLE_se_unadj = five_cell_MLE_se_unadj, five_cell_MLE_se_FPC1=five_cell_MLE_se_FPC1,
                         five_cell_MLE_se_FPC2=five_cell_MLE_se_FPC2),
      
      rs_MLE=list(rs_MLE_se_unadj=rs_MLE_se_unadj, rs_MLE_se_FPC=rs_MLE_se_FPC),
      
      chapman=list(chapman_se=chapman_se)),
    
    ci = list(
      five_cell_MLE=list(wald_five_cell_MLE_unadj=wald_five_cell_MLE_unadj,
                         wald_five_cell_MLE_FPC1=wald_five_cell_MLE_FPC1,
                         wald_five_cell_MLE_FPC2=wald_five_cell_MLE_FPC2,
                         bayes_five_cell_MLE_unadj=bayes_five_cell_MLE_unadj,
                         bayes_five_cell_MLE_FPC1=bayes_five_cell_MLE_FPC1,
                         bayes_five_cell_MLE_FPC2=bayes_five_cell_MLE_FPC2),
      
      rs_MLE=list(wald_rs_unadj=wald_rs_unadj, wald_rs_FPC=wald_rs_FPC, bayes_rs=bayes_rs),
      
      chapman=list(chapman_ci=chapman_ci,
                   rivest_ci=rivest_ci)
    )
    
  )
  
  
}



