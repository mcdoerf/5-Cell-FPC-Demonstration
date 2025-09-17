
###Implementing Simulation Studies.

library(flextable)
library(officer)
library(dplyr)
library(tidyr) 

###This function creates our population
##Ntot is total number of people in the population.  N_flu is the number of flu cases.  psymp_flu is the proportion of
##those that have the flu that are symptotmatic.  psymp_healthy is the proportion of those who do not have the flu that are symptomatic.

population_generator<-function(Ntot, N_flu, psymp_flu, psymp_healthy){
  population<-data.frame(  #Generating population
    id=1:Ntot,
    flu=c(rep(1, N_flu), rep(0, Ntot-N_flu)),
    symptomatic=c(rep(1, ceiling(psymp_flu*N_flu)), rep(0, N_flu-ceiling(psymp_flu*N_flu)), rep(1, ceiling(psymp_healthy*(Ntot-N_flu))), rep(0, Ntot-N_flu-ceiling(psymp_healthy*(Ntot-N_flu))))
  )
  
  
  return(population)
  
}



##This function takes our population, and samples from it (implements streams 1 and 2).  Then, returns a vector of the counts
##n15, n2, n4, n6, and n37.  Note that p1_symp is the probability of being sampled into stream 1 for those with symptoms.  p1_nonysymp is the probability of being sampled into stream 1 for those without symptoms.  In this way, we can make stream 1 a nonrepresentative sample by having it oversample symptomatics, and symptomatics are disproportionately infected.  So, stream 1 would result in infecteds being oversampled.

simulate_data<-function(population, p1_symp, p1_nonsymp, p2){
  
  sim_pop<-population
  
  Ntot<-nrow(population)
  
  for (i in 1:Ntot) {
    sim_pop$stream1[i] <- ifelse(
      population$symptomatic[i] == 1,
      rbinom(1, 1, p1_symp),
      rbinom(1, 1, p1_nonsymp)
    )
  }
  
  n_ones <- round(Ntot * p2)
  sim_pop$stream2 <- sample(c(rep(1, n_ones), rep(0, Ntot - n_ones)))  ##Stream 2 (the anchor stream) is a simple random sample without replacement).
  #sim_pop$stream2<-rbinom(Ntot, 1, p2)
  
  
  n1<-sum(sim_pop$stream1==1 & sim_pop$stream2==1 & sim_pop$flu==0)
  n2<-sum(sim_pop$stream1==1 & sim_pop$stream2==1 & sim_pop$flu==1)
  n3<-sum(sim_pop$stream1==1 & sim_pop$stream2==0 & sim_pop$flu==0)
  n4<-sum(sim_pop$stream1==1 & sim_pop$stream2==0 & sim_pop$flu==1)
  n5<-sum(sim_pop$stream1==0 & sim_pop$stream2==1 & sim_pop$flu==0)
  n6<-sum(sim_pop$stream1==0 & sim_pop$stream2==1 & sim_pop$flu==1)
  n7<-sum(sim_pop$stream1==0 & sim_pop$stream2==0)
  
  return(c(n1+n5, n2, n4, n6, n3+n7))
  
}



########sim_study implements our simulation study.

sim_study <- function(nsim = 1000, Ntot, N_flu,
                      psymp_flu, psymp_healthy, p1_symp, p1_nonsymp, p2, alpha = 0.05, postdraws, seed = 123) {
  
  
  
  set.seed(seed)
  
  population<-population_generator(Ntot=Ntot, N_flu=N_flu, psymp_flu=psymp_flu, psymp_healthy=psymp_healthy) #Generate our fixed,
  #finite population.
  
  # store results for each simulation and method
  results <- vector("list", nsim)
  
  
  for (i in seq_len(nsim)) {
    # 1. simulate data
    dat <- simulate_data(population=population, p1_symp=p1_symp, p1_nonsymp=p1_nonsymp, p2=p2)
    
    # 2. call inference()
    #    tibble(method, est, lwr, upr, se)
    inf <- inference(dat, alpha = alpha, postdraws=postdraws)
    
    
    ## ---- Point estimate table ----
    est_df <- tibble::tibble(
      method = names(inf$estimates),
      est    = unlist(inf$estimates),
      N_flu = N_flu
    )
    
    
    ## ---- SE table ----
    se_df <- purrr::map_dfr(names(inf$se), function(m) {
      purrr::map_dfr(names(inf$se[[m]]), function(se_name) {
        tibble::tibble(
          method = m,
          se_type = se_name,
          se = inf$se[[m]][[se_name]]
        )
      })
    })
    
    
    ## ---- CI table ----
    ci_df <- purrr::map_dfr(names(inf$ci), function(m) {
      purrr::map_dfr(names(inf$ci[[m]]), function(ci_name) {
        ci <- inf$ci[[m]][[ci_name]]
        tibble::tibble(
          method = m,
          ci_type = ci_name,
          lwr = ci[1],
          upr = ci[2],
          coverage = (ci[1] <= N_flu & N_flu <= ci[2]),
          width = ci[2] - ci[1],
          N_flu=N_flu
        )
      })
    })
    
    
    results[[i]] <- list(est_df = est_df, se_df = se_df, ci_df = ci_df)
  }
  
  
  
  
  # ---- Combine results ----
  est_all <- dplyr::bind_rows(lapply(results, `[[`, "est_df"), .id = "sim")
  se_all  <- dplyr::bind_rows(lapply(results, `[[`, "se_df"), .id = "sim")
  ci_all  <- dplyr::bind_rows(lapply(results, `[[`, "ci_df"), .id = "sim")
  
  
  # ---- Summaries ----
  est_summary <- est_all |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      mean_est = mean(est),
      sd_est   = sd(est),
      bias     = mean(est - N_flu),
      .groups = "drop"
    )
  
  
  se_summary <- se_all |>
    dplyr::group_by(method, se_type) |>
    dplyr::summarise(
      mean_se = mean(se),
      sd_se   = sd(se),
      .groups = "drop"
    )
  
  ci_summary <- ci_all |>
    dplyr::group_by(method, ci_type) |>
    dplyr::summarise(
      coverage = mean(coverage),
      avg_width = mean(width),
      .groups = "drop"
    )
  
  
  list(
    est_summary = est_summary,
    se_summary = se_summary,
    ci_summary = ci_summary,
    raw = list(est_all = est_all, se_all = se_all, ci_all = ci_all)
  )
}


##Conducting Simulation Studies

# Step 1: run all scenarios
scenarios <- expand.grid(
  prevalence = c(0.05, 0.1, 0.25, 0.5),
  p2 = c(0.05, 0.1, 0.25, 0.5)   # example
)

all_results1 <- purrr::pmap(scenarios, function(prevalence, p2) {
  sim_study(
    nsim = 10000, 
    Ntot = 10000, 
    N_flu = prevalence * 10000,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000
  )
})

# Step 2: combine summaries with scenario labels
for (i in seq_len(nrow(scenarios))) {
  all_results1[[i]]$scenario <- scenarios[i,]
}

save(all_results1, file = "sim_results11_arXiv.RData")


###Creating Second Set of Tables for arXiv Manuscript with Ntot=1000



##Building Word Document to Store Results

# Step 1: run all scenarios
scenarios <- expand.grid(
  prevalence = c(0.05, 0.1, 0.25, 0.5),
  p2 = c(0.05, 0.1, 0.25, 0.5)   # example
)

all_results2 <- purrr::pmap(scenarios, function(prevalence, p2) {
  sim_study(
    nsim = 10000, 
    Ntot = 1000, 
    N_flu = prevalence * 1000,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=12345
  )
})

# Step 2: combine summaries with scenario labels
for (i in seq_len(nrow(scenarios))) {
  all_results2[[i]]$scenario <- scenarios[i,]
}

save(all_results2, file = "sim_results2_arXiv.RData")


###Third set of scenarios: Ntot even smaller (500).  

# Step 1: run all scenarios
scenarios <- expand.grid(
  prevalence = c(0.05, 0.1, 0.25, 0.5),
  p2 = c(0.05, 0.1, 0.25, 0.5)   # example
)

all_results3 <- purrr::pmap(scenarios, function(prevalence, p2) {
  sim_study(
    nsim = 10000, 
    Ntot = 500, 
    N_flu = prevalence * 500,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=123456
  )
})

# Step 2: combine summaries with scenario labels
for (i in seq_len(nrow(scenarios))) {
  all_results3[[i]]$scenario <- scenarios[i,]
}

save(all_results3, file = "sim_results3_arXiv.RData")





###Fourth set of scenarios: Ntot=1000, but making stream 1 even more nonrepresentative.

# Step 1: run all scenarios
scenarios <- expand.grid(
  prevalence = c(0.05, 0.1, 0.25, 0.5),
  p2 = c(0.05, 0.1, 0.25, 0.5)   # example
)

all_results4 <- purrr::pmap(scenarios, function(prevalence, p2) {
  sim_study(
    nsim = 10000, 
    Ntot = 1000, 
    N_flu = prevalence * 1000,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.9, p1_nonsymp = 0.1,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=1234567
  )
})

# Step 2: combine summaries with scenario labels
for (i in seq_len(nrow(scenarios))) {
  all_results4[[i]]$scenario <- scenarios[i,]
}

save(all_results4, file = "sim_results4_arXiv.RData")


###Fifth set of scenarios: Ntot=1000, make stream 1 even more nonrepresentative.

# Step 1: run all scenarios
scenarios <- expand.grid(
  prevalence = c(0.05, 0.1, 0.25, 0.5),
  p2 = c(0.05, 0.1, 0.25, 0.5)   # example
)

all_results5 <- purrr::pmap(scenarios, function(prevalence, p2) {
  sim_study(
    nsim = 10000, 
    Ntot = 1000, 
    N_flu = prevalence * 1000,
    psymp_flu = 0.9, psymp_healthy = 0.1,
    p1_symp = 0.9, p1_nonsymp = 0.1,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=12345678
  )
})

# Step 2: combine summaries with scenario labels
for (i in seq_len(nrow(scenarios))) {
  all_results5[[i]]$scenario <- scenarios[i,]
}

save(all_results5, file = "sim_results5_arXiv.RData")

