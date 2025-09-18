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
  
  data1<-data_decomp(data)  ##Decompose the data.
  
  n_rs<-data1$n2+data1$n6   ##Number of individuals that test positive in the random sample.
  
  n<-data1$n15+data1$n2+data1$n6 ##Number of individuals in the random sample.
  
  p_rs<-n_rs/n #Point Estimate of Prevalence
  
  est=(p_rs)*(data1$Ntot) #Point Estimate of Number of Diseased Individuals
  
  p_rs2<-max(n_rs, 0.5)/n #Add 0.5 if the number of disesaed in the random sample is 0.
  
  var_rs<-(data1$Ntot)^2*p_rs2*(1-p_rs2)/(n)   ###Non FPC adjusted variance
  
  FPC<-min((data1$n*(data1$Ntot-data1$n))/(data1$Ntot*(data1$n-1)),1) ##Cochran's FPC
  
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
  data<-data_decomp(data)
  
  n1.<-data$n11+data$n10
  n.1<-data$n11+data$n01
  
  nchapman<-((n1.+1)*(n.1+1))/(data$n11+1)-1  ##Chapman's point estimator
  
  var_chap<-(data$n11+data$n10+1)*(data$n11+data$n01+1)*(data$n10)*(data$n01)/((data$n11+1)^2*(data$n11+2))  ##Variance of Chapman's estimator.
  
  list(est=nchapman, se=sqrt(var_chap))
}



# 5 cell MLE estimator & variance


five_cell_estimator<-function(data, type="Unadjusted"){
  
  
  data1<-data_decomp(data=data)
  
  w= (data1$n2+data1$n4)/data1$Ntot
  
  eta=(data1$n6)/(data1$n15+data1$n6)
  
  N_tot_star=data1$n15+data1$n6+data1$n37
  
  n_rs_star=data1$n15+data1$n6
  
  pi_all = w+eta*(1-w)
  
  
  ###5 cell MLE estimate
  
  N_MLE = data1$Ntot*pi_all
  
  
  ##Non FPC Corrected standard error (based on Delta Method)
  
  p15hat<-data1$n15/data1$Ntot
  
  p2hat<-data1$n2/data1$Ntot
  
  p4hat<-data1$n4/data1$Ntot
  
  p6hat<-data1$n6/data1$Ntot
  
  p37hat<-data1$n37/data1$Ntot
  
  phat <- c(p15hat, p2hat, p4hat, p6hat, p37hat)
  
  
  gp <- c(-p6hat*p37hat/(p15hat+p6hat)^2, 1, 1, 1+p37hat*p15hat/(p15hat+p6hat)^2, p6hat/(p15hat+p6hat))

  
  Sigma <- -outer(phat, phat) / data1$Ntot
  
  diag(Sigma) <- phat * (1 - phat) / data1$Ntot
  
  
  var <- t(gp)%*%Sigma%*%gp*(data1$Ntot^2)  #Unadjusted Variance
  
  #FPC1 Corrected Variance
  
  eta_new=max(data1$n6, 0.5)/(max(data1$n15, 0.5)+max(data1$n6, 0.5))
  
  var_eta=min(n_rs_star*(N_tot_star-n_rs_star)/(N_tot_star*(n_rs_star-1)) ,1)*(eta_new*(1-eta_new)/n_rs_star)
  
  var_pi_FPC1=(1-w)^2*var_eta
  
  var_FPC1=data1$Ntot^2*var_pi_FPC1  ##FPC Adjustment 1
  
  #FPC corrected standard error 2 (acknowledging variability in w)
  
  var_FPC2<-data1$Ntot^2*(var_pi_FPC1+((1-eta_new)^2)*w*(1-w)/data1$Ntot)
  
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

###Generic Wald CI function

wald_ci <- function(est, se, alpha=0.05) {
  z<-qnorm((1-alpha/2), mean=0, sd=1, lower.tail=TRUE)
  lwr <-  est - z * se
  upr <- est + z * se
  c(lwr, upr)
}


###Generic Wald CI function, for cases where we know n15 and Ntot-n15

wald_ci2 <- function(est, se, alpha=0.05, Ntot, n15, nc) {
  z<-qnorm((1-alpha/2), mean=0, sd=1, lower.tail=TRUE)
  lwr <- max(nc, est - z * se)         #Number of diseased has to be greater than nc.
  upr <- min(Ntot-n15, est + z * se)   #Number of diseased cannot exceed Ntot-n15
  c(lwr, upr)
}



#Sadinle CI to accompany Chapman's estimate
sadinle_ci<-function(data, alpha=0.05){
  data<-data_decomp(data)
  
  h_n<-(data$n10+0.5)*(data$n01+0.5)/(data$n11+0.5)
  
  sigma_tl<-sqrt(1/(data$n11+0.5)+1/(data$n10+0.5)+1/(data$n01+0.5)+(data$n11+0.5)*(1/(data$n10+0.5))*(1/(data$n01+0.5)))
  
  z<-qnorm(1-(alpha/2), mean=0, sd=1, lower.tail=TRUE)
  
  lwr<-max(c(data$nc, (data$nc-0.5)+h_n*exp(-z*sigma_tl)))
  
  upr<-min(data$Ntot, (data$nc-0.5)+h_n*exp(z*sigma_tl))
  
  c(lwr, upr)}


#Bayes CI Generator for 5 cell MLE

bayes_ci_5<-function(data, alpha=0.05, type="Unadjusted", postdraws=10000){
  
  data1<-data_decomp(data=data)
  
  n15<-data1$n15
  n2<-data1$n2
  n4<-data1$n4
  n6<-data1$n6
  n37<-data1$n37
  nc<-data1$nc
  Ntot<-data1$Ntot
  
  pstar_post<-rdirichlet(postdraws, c(n15, n2, n4, n6, n37)+0.5)
  
  nstar<-Ntot*pstar_post
  
  Nhat_star<-nstar[,2]+nstar[,3]+nstar[,4]*((nstar[,1]+nstar[,4]+nstar[,5]))/(nstar[,1]+nstar[,4])
  
  if (type=="Unadjusted"){
    
    lwr<-max(as.numeric(quantile(Nhat_star, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star, 1-alpha/2)), Ntot-n15)
    
    c(lwr, upr)
    
  }
  
  else if (type=="FPC1"){
    
    
    a1=min(five_cell_estimator(data=data, type="FPC1")$se/five_cell_estimator(data=data, type="Unadjusted")$se,1)
    
    b1=five_cell_estimator(data=data, type="FPC1")$est *(1-a1)
    
    Nhat_star1<-as.vector(a1)*Nhat_star+as.vector(b1)
    
    lwr<-max(as.numeric(quantile(Nhat_star1, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star1, 1-alpha/2)), Ntot-n15)
    
    c(lwr, upr)
    
  }
  
  else if (type=="FPC2"){
    
    a1=min(five_cell_estimator(data=data, type="FPC2")$se/five_cell_estimator(data=data, type="Unadjusted")$se,1)
    
    b1=five_cell_estimator(data=data, type="FPC2")$est *(1-a1)
    
    Nhat_star1<-as.vector(a1)*Nhat_star+as.vector(b1)
    
    lwr<-max(as.numeric(quantile(Nhat_star1, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star1, 1-alpha/2)), Ntot-n15)
    
    c(lwr, upr) 
  }
  
}


###Function for Implementing Bayesian Credible CI for Estimator Based on Random Sample Only
bayes_ci_rs<-function(data, alpha=0.05, type="FPC"){
  
  data1<-data_decomp(data)
  
  n11<-data1$n11
  
  n10<-data1$n10
  
  n01<-data1$n01
  
  n<-data1$n
  
  Ntot<-data1$Ntot
  
  nrs_pos=n11+n01
  
  prs=rs_estimator(data, type="Unadjusted")$est/Ntot
  
  alphpost = nrs_pos + 0.5
  betpost = n - nrs_pos + 0.5
  
  LL_JeffreysForP = qbeta(alpha/2, alphpost, betpost)
  UL_JeffreysForP = qbeta(1-alpha/2, alphpost, betpost)
  
  if(nrs_pos == 0){ LL_JeffreysForP = 0 }
  if(nrs_pos == n){ UL_JeffreysForP = 1 }
  LL_Jeffreys = Ntot*LL_JeffreysForP
  UL_Jeffreys = Ntot*UL_JeffreysForP
  
  ####Adjusting Jeffry's CI for FPC
  
  a=sqrt(min((Ntot-n)*n/((Ntot)*(n-1)),1))
  b=prs*(1-a)
  LL_JeffreysForPFPC = a*LL_JeffreysForP + b
  UL_JeffreysForPFPC = a*UL_JeffreysForP + b
  if(nrs_pos == 0){ LL_JeffreysForPFPC = 0 }
  if(nrs_pos == n){ UL_JeffreysForPFPC = 1 }
  LL_JeffreysFPC = Ntot*LL_JeffreysForPFPC
  UL_JeffreysFPC = Ntot*UL_JeffreysForPFPC
  
  if (type=="Unadjusted"){
    c(LL_Jeffreys, UL_Jeffreys)
  }
  
  if (type=="FPC"){
    c(LL_JeffreysFPC, UL_JeffreysFPC)
  }
  
  
}


##########Inference function



inference<-function(data, alpha, postdraws){
  
  data1<-data_decomp(data)
  
  n1<-data1$n1
  n2<-data1$n2
  n3<-data1$n3
  n4<-data1$n4
  n5<-data1$n5
  n6<-data1$n6
  n7<-data1$n7
  
  n15<-data1$n15
  n37<-data1$n37
  
  n11<-data1$n11
  n10<-data1$n10
  n01<-data1$n01
  
  Ntot<-data1$Ntot
  
  nc<-data1$nc
  n<-data1$n
  
  
  ##Point Estimators
  
  five_cell_MLE<-five_cell_estimator(data=data, type="Unadjusted")$est
  
  rs_MLE<-rs_estimator(data=data, type="Unadjusted")$est
  
  chapman<-chapman_estimator(data=data)$est
  
  ##Standard Errors
  
  five_cell_MLE_se_unadj<-five_cell_estimator(data=data, type="Unadjusted")$se
  
  five_cell_MLE_se_FPC1<-five_cell_estimator(data=data, type="FPC1")$se
  
  five_cell_MLE_se_FPC2<-five_cell_estimator(data=data, type="FPC2")$se
  
  rs_MLE_se_unadj<-rs_estimator(data=data, type="Unadjusted")$se
  
  rs_MLE_se_FPC<-rs_estimator(data=data, type="FPC")$se
  
  chapman_se<-chapman_estimator(data=data)$se
  
  #Confidence Intervals
  
  #5 cell MLE
  wald_five_cell_MLE_unadj<-wald_ci2(five_cell_MLE, five_cell_MLE_se_unadj, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_five_cell_MLE_FPC1<-wald_ci2(five_cell_MLE, five_cell_MLE_se_FPC1, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_five_cell_MLE_FPC2<-wald_ci2(five_cell_MLE, five_cell_MLE_se_FPC2, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  bayes_five_cell_MLE_unadj<-bayes_ci_5(data=data, alpha=alpha, type="Unadjusted", postdraws=postdraws)
  
  bayes_five_cell_MLE_FPC1<-bayes_ci_5(data=data, alpha=alpha, type="FPC1", postdraws=postdraws)
  
  bayes_five_cell_MLE_FPC2<-bayes_ci_5(data=data, alpha=alpha, type="FPC2", postdraws=postdraws)
  
  #RS Estimator
  wald_rs_unadj<-wald_ci(rs_MLE, rs_MLE_se_unadj, alpha=alpha)
  
  wald_rs_FPC<-wald_ci(rs_MLE, rs_MLE_se_FPC, alpha=alpha)
  
  bayes_rs<-bayes_ci_rs(data=data, alpha=alpha, type="FPC")
  
  #Chapman Estimator
  chapman_ci<-sadinle_ci(data=data, alpha=alpha)
  
  
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
      
      chapman=list(chapman_ci=chapman_ci)
    )
    
  )
  
  
}














