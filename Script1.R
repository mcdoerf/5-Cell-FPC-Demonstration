####-----------------
###Helper Functions
###-----------------




library(rBeta2009)  #For rdirichlet() function

data_decomp<-function(data, case){  ###Decomposes data.  When data is for 5 cell case, the function assumes input data is
  #of the form (n15, n2, n4, n6, n37).  When data is for 7 cell case, the function assumes input data is of the form
  #(n1, n2, n3, n4, n5, n6, n7). For the 3 cell case, the function assumes input data is of the form (n11, n10, n01, n, Ntot)
  #where n is the number of individuals selected in the anchor stream.
  
  if (case != 3 & case !=5 & case !=7){
    print("Error: Please Enter 3, 5, or 7 for case")
  }
  
  else if (case==5){
    n15=data[1]
    n2=data[2]
    n4=data[3]
    n6=data[4]
    n37=data[5]
    
    return(list(n15=n15, n2=n2, n4=n4, n6=n6, n37=n37, n11=n2, n10=n4, n01=n6, nc=n2+n4+n6, n=n15+n2+n6,  Ntot=sum(data)))}
  
  else if (case==7){
    n1=data[1]
    n2=data[2]
    n3=data[3]
    n4=data[4]
    n5=data[5]
    n6=data[6]
    n7=data[7]
    
    return(list(n1=n1, n2=n2, n3=n3, n4=n4, n5=n5, n6=n6, n7=n7, n11=n2, n10=n4, n01=n6, n15=n1+n5, n37=n3+n7, nc=n2+n4+n6, n=n1+n5+n2+n6, Ntot=sum(data)))
  }
  
  else if (case==3){
    n11=data[1]
    n10=data[2]
    n01=data[3]
    n=data[4]
    Ntot=data[5]
    nc=n11+n10+n01
    n2=n11
    n4=n10
    n6=n01
    
    return(list(n11=n11, n10=n10, n01=n01, n=n, Ntot=Ntot, nc=nc, n2=n2, n4=n4, n6=n6))
  }
  
  
}




########Estimators



# RS estimator & variance
rs_estimator <- function(data, case=5, type="Unadjusted") {
  
  data1<-data_decomp(data, case)
  
  n_rs<-data1$n2+data1$n6   ##Number of individuals that test positive in the random sample.
  
  n<-data1$n15+data1$n2+data1$n6 ##Number of individuals in the random sample.
  
  p_rs<-n_rs/n #Point Estimate of Prevalence
  
  est=(p_rs)*(data1$Ntot) #Point Estimate of Number of Diseased Individuals
  
  p_rs2<-max(n_rs, 0.5)/n
  
  var_rs<-(data1$Ntot)^2*p_rs2*(1-p_rs2)/(n)   ###Non FPC adjusted variance
  
  FPC1<-min((data1$Ntot-n)/(data1$Ntot-1),1)  ##FPC from JSSM Paper
  
  FPC2<-min((data1$n*(data1$Ntot-data1$n))/(data1$Ntot*(data1$n-1)),1)
  
  ##FPC from AJE Paper
  
  var_rs1<-FPC1*var_rs  ##FPC adjusted variance of random sample estimator from JSSM paper
  
  var_rs2<-FPC2*var_rs  ## Cochran FPC adjusted variance of random sample estimator (from AJE paper).. 
  
  
  if (type=="Unadjusted"){
    return(list(est = est, se = sqrt(var_rs)))
  }
  
  if (type=="FPC1"){
    return(list(est = est, se = sqrt(var_rs1)))
  }
  
  if (type=="FPC2"){
    return(list(est = est, se = sqrt(var_rs2)))
  }
  
}


# Chapman estimator & variance
chapman_estimator <- function(data, case=5) {
  data<-data_decomp(data, case)
  
  n1.<-data$n11+data$n10
  n.1<-data$n11+data$n01
  
  nchapman<-((n1.+1)*(n.1+1))/(data$n11+1)-1  ##Chapman's point estimator
  
  var_chap<-(data$n11+data$n10+1)*(data$n11+data$n01+1)*(data$n10)*(data$n01)/((data$n11+1)^2*(data$n11+2))  ##Variance of Chapman's estimator.
  
  list(est=nchapman, se=sqrt(var_chap))
}



# 5 cell MLE estimator & variance



five_cell_estimator<-function(data, case=5, type="Unadjusted"){
  
  
  data1<-data_decomp(data=data, case=case)
  
  p1= (data1$n2+data1$n4)/data1$Ntot
  
  eta=(data1$n6)/(data1$n15+data1$n6)
  
  N_tot_star=data1$n15+data1$n6+data1$n37
  
  n_rs_star=data1$n15+data1$n6
  
  pi_all = p1+eta*(1-p1)
  
  
  ###5 cell MLE estimate
  
  N_MLE = data1$Ntot*pi_all
  
  
  ##Non FPC Corrected standard error (from Wenhao) (based on Delta Method)
  
  p15hat<-data1$n15/data1$Ntot
  
  
  p2hat<-data1$n2/data1$Ntot
  
  
  
  p4hat<-data1$n4/data1$Ntot
  
  
  p6hat<-data1$n6/data1$Ntot
  
  
  
  p37hat<-data1$n37/data1$Ntot
  
  
  gp <- c(-p6hat*p37hat/(p15hat+p6hat)^2, 1, 1, 1+p37hat*p15hat/(p15hat+p6hat)^2, p6hat/(p15hat+p6hat))
  
  
  phat <- c(p15hat, p2hat, p4hat, p6hat, p37hat)
  
  Sigma <- -outer(phat, phat) / data1$Ntot
  
  diag(Sigma) <- phat * (1 - phat) / data1$Ntot
  
  
  var <- t(gp)%*%Sigma%*%gp*(data1$Ntot^2)  #Unadjusted Variance
  
  #FPC1 Corrected Variance
  
  
  eta_new=max(data1$n6, 0.5)/(max(data1$n15, 0.5)+max(data1$n6, 0.5))
  
  var_eta=min(n_rs_star*(N_tot_star-n_rs_star)/(N_tot_star*(n_rs_star-1)) ,1)*(eta_new*(1-eta_new)/n_rs_star)
  
  
  var_pi_FPC1=(1-p1)^2*var_eta
  
  var_FPC1=data1$Ntot^2*var_pi_FPC1  ##FPC Adjustment 1
  
  #FPC corrected standard error 2 (acknowledging variability in p1) (approach outlined in my 7-3-25 Notes)
  
  var_FPC2<-data1$Ntot^2*(var_pi_FPC1+((1-eta_new)^2)*p1*(1-p1)/data1$Ntot)
  
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


###3 Cell Estimator


three_cell_estimator<-function(data, case=3){
  
  data1<-data_decomp(data, case=case)
  n11<-data1$n11
  n10<-data1$n10
  n01<-data1$n01
  Ntot<-data1$Ntot
  n<-data1$n
  nc<-data1$nc
  
  psi_hat<-n/Ntot  ##Estimate of sampling probability into stream 2.
  
  est<-n11+n10+n01/psi_hat  ##MLE based on 3 cell count.
  
  if (n01 != 0){
    
    var<-n01*(1-psi_hat)/(psi_hat^2)  ###Equation (3) JSSM
  }
  
  else if (n01 ==0){
    
    var<-(n01+0.5)*(1-psi_hat)/(psi_hat^2) 
    
  }
  
  
  list(est=est, se=sqrt(var))
  
}


#Confidence Intervals

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
  lwr <- max(nc, est - z * se)
  upr <- min(Ntot-n15, est + z * se)
  c(lwr, upr)
}



#Sadinle CI to accompany Chapman's estimate
sadinle_ci<-function(data, case=5, alpha=0.05){
  data<-data_decomp(data, case)
  
  h_n<-(data$n10+0.5)*(data$n01+0.5)/(data$n11+0.5)
  sigma_tl<-sqrt(1/(data$n11+0.5)+1/(data$n10+0.5)+1/(data$n01+0.5)+(data$n11+0.5)*(1/(data$n10+0.5))*(1/(data$n01+0.5)))
  
  
  z<-qnorm(1-(alpha/2), mean=0, sd=1, lower.tail=TRUE)
  
  lwr<-max(c(data$nc, (data$nc-0.5)+h_n*exp(-z*sigma_tl)))
  upr<-min(data$Ntot, (data$nc-0.5)+h_n*exp(z*sigma_tl))
  
  c(lwr, upr)}


#Bayes CI Generator for 5 cell MLE

bayes_ci_5<-function(data, case=5, alpha=0.05, type="Unadjusted", postdraws=10000){
  
  data1<-data_decomp(data=data, case=case)
  
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
    
    
    a1=five_cell_estimator(data=data, case=case, type="FPC1")$se/five_cell_estimator(data=data, case=case, type="Unadjusted")$se
    
    b1=five_cell_estimator(data=data, case=case, type="FPC1")$est *(1-a1)
    
    Nhat_star1<-as.vector(a1)*Nhat_star+as.vector(b1)
    
    lwr<-max(as.numeric(quantile(Nhat_star1, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star1, 1-alpha/2)), Ntot-n15)
    
    c(lwr, upr)
    
  }
  
  else if (type=="FPC2"){
    
    a1=five_cell_estimator(data=data, case=case, type="FPC2")$se/five_cell_estimator(data=data, case=case, type="Unadjusted")$se
    
    b1=five_cell_estimator(data=data, case=case, type="FPC2")$est *(1-a1)
    
    Nhat_star1<-as.vector(a1)*Nhat_star+as.vector(b1)
    
    lwr<-max(as.numeric(quantile(Nhat_star1, alpha/2)), nc)
    
    upr<-min(as.numeric(quantile(Nhat_star1, 1-alpha/2)), Ntot-n15)
    
    c(lwr, upr) 
  }
  
}


bayes_ci_3<-function(data, case=3, alpha=0.05, type="Unadjusted", postdraws=10000){
  data1<-data_decomp(data, case=case)
  n11<-data1$n11
  n10<-data1$n10
  n01<-data1$n01
  n<-data1$n
  nc<-data1$nc
  Ntot<-data1$Ntot
  
  psi_hat<-n/Ntot 
  pstarcond<-rdirichlet(postdraws, c(n11+1/2, n10+1/2, n01+1/2))
  ptotcheck = rowSums(pstarcond)
  p1starpost = rowSums(pstarcond[,1:2])
  p1post = psi_hat*p1starpost/(psi_hat*p1starpost + pstarcond[,3])
  p2giv1starpost = pstarcond[,1]/(pstarcond[,1] + pstarcond[,2])
  
  p11post = p2giv1starpost*p1post
  p10post = (1-p2giv1starpost)*p1post
  p01post = psi_hat*(1-p1post)
  pcpost = p11post + p10post + p01post 
  
  Nfirst = round((n11+n10+n01)/pcpost)
  Nfirst[Nfirst == 0] <- 1
  ncnew = rbinom(postdraws, size = Nfirst, prob = pcpost)
  
  n11post = ncnew*pstarcond[,1]
  n10post = ncnew*pstarcond[,2]
  n01post = ncnew*pstarcond[,3]
  
  Npsipost = n11post + n10post + n01post/psi_hat
  
  if (type=="Unadjusted"){
    
    lwr<-max(c(nc, as.numeric(quantile(Npsipost, probs=alpha/2)))) ##Lower bound
    
    upr<-min(c(Ntot-(n-n11-n01), as.numeric(quantile(Npsipost, probs=1-alpha/2)))) ##Upper bound
    
    return(c(lwr, upr))
    
  }
  
  else if (type=="COMPA"){
    
    var_chap<-(chapman_estimator(data=data, case=case)$se)^2 ##Variance of Chapman's estimator.
    
    
    var_rs1<-(rs_estimator(data=data, case=case, type="FPC2")$se)^2 ##FPC adjusted variance of random sample estimator from JSSM paper
    
    #var_rs2<-(Ntot)^2*((Ntot-n)*n/((Ntot)*(n-1)))*((p_rs*(1-p_rs))/n)  ## Cochran FPC adjusted variance of random sample estimator (from AJE paper).. 
    
    n_mle3<-three_cell_estimator(data=data, case=case)$est  ##MLE based on 3 cell count.
    
    var_compa<-(var_chap+var_rs1)/4
    
    var1<-(three_cell_estimator(data=data, case=case)$se)^2
    
    a<-sqrt(min(c(var_compa, var1))/var1)
    
    b<-n_mle3*(1-a)
    
    Npsipost<-a*Npsipost+b
    
    lwr<-max(c(nc, as.numeric(quantile(Npsipost, probs=alpha/2)))) ##Lower bound
    
    upr<-min(c(Ntot-(n-n11-n01), as.numeric(quantile(Npsipost, probs=1-alpha/2)))) ##Upper Bound
    
    return(c(lwr, upr))
  }
  
  else if (type=="COMPB"){
    
    n11_new<-max(n11, 0.5)
    n10_new<-max(n10, 0.5)
    n01_new<-max(n01, 0.5)
    
    var_LP<-(n11_new+n10_new)*(n11_new+n01_new)*(n10_new)*(n01_new)/(n11_new^3)
    
    
    var_rs1<-(rs_estimator(data=data, case=case, type="FPC2")$se)^2 ##FPC adjusted variance of random sample estimator from JSSM paper
    
    #var_rs2<-(Ntot)^2*((Ntot-n)*n/((Ntot)*(n-1)))*((p_rs*(1-p_rs))/n)  ## Cochran FPC adjusted variance of random sample estimator (from AJE paper).. 
    
    n_mle3<-three_cell_estimator(data=data, case=case)$est  ##MLE based on 3 cell count.
    
    var_compb<-1/(1/var_rs1 + 1/var_LP)
    
    var1<-(three_cell_estimator(data=data, case=case)$se)^2
    
    a<-sqrt(min(c(var_compb, var1))/var1)
    
    b<-n_mle3*(1-a)
    
    Npsipost<-a*Npsipost+b
    
    lwr<-max(c(nc, as.numeric(quantile(Npsipost, probs=alpha/2)))) ##Lower bound
    
    upr<-min(c(Ntot-(n-n11-n01), as.numeric(quantile(Npsipost, probs=1-alpha/2)))) ##Upper Bound
    
    return(c(lwr, upr))
  }
  
}


###Function for Implementing Bayesian Credible CI for Estimator Based on Random Sample Only as Described in AJE Paper
bayes_ci_rs<-function(data, case, alpha=0.05, type="FPC"){
  
  data1<-data_decomp(data, case=case)
  
  n11<-data1$n11
  
  n10<-data1$n10
  
  n01<-data1$n01
  
  n<-data1$n
  
  Ntot<-data1$Ntot
  
  nrs_pos=n11+n01
  
  prs=rs_estimator(data, case=case, type="Unadjusted")$est/Ntot
  
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



inference<-function(data, case, alpha, postdraws){
  
  data1<-data_decomp(data, case=case)
  
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
  
  five_cell_MLE<-five_cell_estimator(data=data, case=case, type="Unadjusted")$est
  
  three_cell_MLE<-three_cell_estimator(data=data, case=case)$est
  
  rs_MLE<-rs_estimator(data=data, case=case, type="Unadjusted")$est
  
  chapman<-chapman_estimator(data=data, case=case)$est
  
  ##Standard Errors
  
  five_cell_MLE_se_unadj<-five_cell_estimator(data=data, case=case, type="Unadjusted")$se
  
  five_cell_MLE_se_FPC1<-five_cell_estimator(data=data, case=case, type="FPC1")$se
  
  five_cell_MLE_se_FPC2<-five_cell_estimator(data=data, case=case, type="FPC2")$se
  
  three_cell_MLE_se<-three_cell_estimator(data=data, case=case)$se
  
  rs_MLE_se_unadj<-rs_estimator(data=data, case=case, type="Unadjusted")$se
  
  rs_MLE_se_FPC1<-rs_estimator(data=data, case=case, type="FPC1")$se
  
  rs_MLE_se_FPC2<-rs_estimator(data=data, case=case, type="FPC2")$se
  
  chapman_se<-chapman_estimator(data=data, case=case)$se
  
  #Confidence Intervals
  
  #5 cell MLE
  wald_five_cell_MLE_unadj<-wald_ci2(five_cell_MLE, five_cell_MLE_se_unadj, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_five_cell_MLE_FPC1<-wald_ci2(five_cell_MLE, five_cell_MLE_se_FPC1, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  wald_five_cell_MLE_FPC2<-wald_ci2(five_cell_MLE, five_cell_MLE_se_FPC2, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  bayes_five_cell_MLE_unadj<-bayes_ci_5(data=data, case=case, alpha=alpha, type="Unadjusted", postdraws=postdraws)
  
  bayes_five_cell_MLE_FPC1<-bayes_ci_5(data=data, case=case, alpha=alpha, type="FPC1", postdraws=postdraws)
  
  bayes_five_cell_MLE_FPC2<-bayes_ci_5(data=data, case=case, alpha=alpha, type="FPC2", postdraws=postdraws)
  
  #3 Cell MLE
  wald_three_cell_MLE<-wald_ci2(three_cell_MLE, three_cell_MLE_se, alpha=alpha, Ntot=Ntot, n15=n15, nc=nc)
  
  bayes_three_cell_MLE_unadj<-bayes_ci_3(data=data, case=case, alpha=alpha, type="Unadjusted", postdraws=postdraws)
  
  bayes_three_cell_MLE_COMPA<-bayes_ci_3(data=data, case=case, alpha=alpha, type="COMPA", postdraws=postdraws)
  
  bayes_three_cell_MLE_COMPB<-bayes_ci_3(data=data, case=case, alpha=alpha, type="COMPB", postdraws=postdraws)
  
  #RS Estimator
  wald_rs_unadj<-wald_ci(rs_MLE, rs_MLE_se_unadj, alpha=alpha)
  
  wald_rs_FPC1<-wald_ci(rs_MLE, rs_MLE_se_FPC1, alpha=alpha)
  
  wald_rs_FPC2<-wald_ci(rs_MLE, rs_MLE_se_FPC2, alpha=alpha)
  
  bayes_rs<-bayes_ci_rs(data=data, case=case, alpha=alpha, type="FPC")
  
  #Chapman Estimator
  chapman_ci<-sadinle_ci(data=data, case=case, alpha=alpha)
  
  
  list(
    estimates = list(five_cell_MLE = five_cell_MLE, three_cell_MLE = three_cell_MLE, rs_MLE=rs_MLE,
                     chapman=chapman),
    se = list(
      
      five_cell_MLE=list(five_cell_MLE_se_unadj = five_cell_MLE_se_unadj, five_cell_MLE_se_FPC1=five_cell_MLE_se_FPC1,
                         five_cell_MLE_se_FPC2=five_cell_MLE_se_FPC2),
      
      three_cell_MLE=list(three_cell_MLE_se=three_cell_MLE_se),
      
      rs_MLE=list(rs_MLE_se_unadj=rs_MLE_se_unadj, rs_MLE_se_FPC1=rs_MLE_se_FPC1, rs_MLE_se_FPC2=rs_MLE_se_FPC2),
      
      chapman=list(chapman_se=chapman_se)),
    
    ci = list(
      five_cell_MLE=list(wald_five_cell_MLE_unadj=wald_five_cell_MLE_unadj,
                         wald_five_cell_MLE_FPC1=wald_five_cell_MLE_FPC1,
                         wald_five_cell_MLE_FPC2=wald_five_cell_MLE_FPC2,
                         bayes_five_cell_MLE_unadj=bayes_five_cell_MLE_unadj,
                         bayes_five_cell_MLE_FPC1=bayes_five_cell_MLE_FPC1,
                         bayes_five_cell_MLE_FPC2=bayes_five_cell_MLE_FPC2),
      
      three_cell_MLE=list(wald_three_cell_MLE=wald_three_cell_MLE,
                          bayes_three_cell_MLE_unadj=bayes_three_cell_MLE_unadj,
                          bayes_three_cell_MLE_COMPA=bayes_three_cell_MLE_COMPA,
                          bayes_three_cell_MLE_COMPB=bayes_three_cell_MLE_COMPB),
      
      rs_MLE=list(wald_rs_unadj=wald_rs_unadj, wald_rs_FPC1=wald_rs_FPC1, wald_rs_FPC2=wald_rs_FPC2, bayes_rs=bayes_rs),
      
      chapman=list(chapman_ci=chapman_ci)
    )
    
  )
  
  
}
















