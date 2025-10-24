
###Code for implementing simulation studies.

library(flextable)
library(officer)
library(dplyr)
library(tidyr) 


###Note: the following code refers to the diseased population as those having the flu. As discussed
###in the paper though, the methods developed can be used for monitoring cumulative cases of other diseases,
###not just the flu.

###This function creates our population
##Ntot is total number of people in the population.  N_flu is the number of flu cases.  psymp_flu is the proportion of
##those that have the flu that are symptomatic.  psymp_healthy is the proportion of those who do not have the flu that are symptomatic.

population_generator<-function(Ntot, N_flu, psymp_flu, psymp_healthy){
  population<-data.frame(  #Generating population
    id=1:Ntot,
    flu=c(rep(1, N_flu), rep(0, Ntot-N_flu))
  )
  
  population$symptomatic=rbinom(Ntot, 1, population$flu*psymp_flu+(1-population$flu)*psymp_healthy)
  return(population)
  
}

##This function takes our population and samples from it (implements streams 1 and 2).  Then, returns a vector of the counts
##n15, n2, n4, n6, and n37.  Note that p1_symp is the probability of being sampled into stream 1 for those with symptoms.  
#p1_nonysymp is the probability of being sampled into stream 1 for those without symptoms.  
#In this way, we can make stream 1 a nonrepresentative sample by having it oversample symptomatics, and symptomatics are disproportionately infected.  
#So, stream 1 would result in infecteds being oversampled.

simulate_data <- function(population, p1_symp, p1_nonsymp, p2) {
  Ntot <- nrow(population)
    sim_pop <- population
    
    # stream 1: biased sample depending on symptoms
    for (i in 1:Ntot) {
      sim_pop$stream1[i] <- ifelse(
        population$symptomatic[i] == 1,
        rbinom(1, 1, p1_symp),
        rbinom(1, 1, p1_nonsymp)
      )
    }
    
    # stream 2: SRS without replacement, exactly round(Ntot*p2) individuals
    n_ones <- round(Ntot * p2)
    sim_pop$stream2 <- sample(c(rep(1, n_ones), rep(0, Ntot - n_ones)))
    
    # cell counts
    n1 <- sum(sim_pop$stream1==1 & sim_pop$stream2==1 & sim_pop$flu==0)
    n2 <- sum(sim_pop$stream1==1 & sim_pop$stream2==1 & sim_pop$flu==1)
    n3 <- sum(sim_pop$stream1==1 & sim_pop$stream2==0 & sim_pop$flu==0)
    n4 <- sum(sim_pop$stream1==1 & sim_pop$stream2==0 & sim_pop$flu==1)
    n5 <- sum(sim_pop$stream1==0 & sim_pop$stream2==1 & sim_pop$flu==0)
    n6 <- sum(sim_pop$stream1==0 & sim_pop$stream2==1 & sim_pop$flu==1)
    n7 <- sum(sim_pop$stream1==0 & sim_pop$stream2==0)
    
      return(c(n1+n5, n2, n4, n6, n3 + n7))
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
    
    
    # ---- Cell counts table ----
    counts <- data_decomp(dat)
    counts_df <- tibble::tibble(
      sim  = i,
      n15  = counts$n15,
      n2   = counts$n2,
      n4   = counts$n4,
      n6   = counts$n6,
      n37  = counts$n37,
      Ntot = counts$Ntot
    )
    
    
    
    
    
    results[[i]] <- list(est_df = est_df, se_df = se_df, ci_df = ci_df, counts_df=counts_df)
  }
  
  
  
  
  # ---- Combine results ----
  est_all <- dplyr::bind_rows(lapply(results, `[[`, "est_df"), .id = "sim")
  se_all  <- dplyr::bind_rows(lapply(results, `[[`, "se_df"), .id = "sim")
  ci_all  <- dplyr::bind_rows(lapply(results, `[[`, "ci_df"), .id = "sim")
  counts_all <- dplyr::bind_rows(lapply(results, `[[`, "counts_df"), .id="sim")
  
  
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
    raw = list(est_all = est_all, se_all = se_all, ci_all = ci_all, counts_all=counts_all)
  )
}


##Conducting Simulation Studies
#-----------------------------------------------------------------------------------------------------------
##############---------------------------Table 5---------------------------------------------------------------

# Run all scenarios
scenarios <- expand.grid(
  N_flu = c(500, 1000, 2500, 5000),
  p2 = c(0.05, 0.1, 0.25, 0.5),
  Ntot=c(10000)
)

all_results1 <- purrr::pmap(scenarios, function(N_flu, p2, Ntot) {
  sim_study(
    nsim = 10000, 
    Ntot =Ntot, 
    N_flu = N_flu,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=123
  )
})

#combine summaries with scenario labels
for (i in seq_len(nrow(scenarios))) {
  all_results1[[i]]$scenario <- scenarios[i,]
}

save(all_results1, file = "sim_results1_arXiv.RData")

#-------------------------------------------------------------------------------------------------------------------
###------------------------------------------------------Table 6----------------------


# Step 1: run all scenarios
scenarios2 <- expand.grid(
  N_flu = c(50, 100, 250, 500),
  p2 = c(0.05, 0.1, 0.25, 0.5),
  Ntot=c(1000)
  
)

all_results2 <- purrr::pmap(scenarios2, function(N_flu, p2, Ntot) {
  sim_study(
    nsim = 10000, 
    Ntot = Ntot, 
    N_flu =N_flu,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=1234
  )
})

# combine summaries with scenario labels
for (i in seq_len(nrow(scenarios2))) {
  all_results2[[i]]$scenario <- scenarios2[i,]
}

save(all_results2, file = "sim_results2_arXiv.RData")


#--------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------Table (B1) (Ntot =250)---------


scenariosA3<- expand.grid(
  N_flu = c(13, 25, 63, 125),
  p2 = c(0.05, 0.1, 0.25, 0.5),
  Ntot=c(250)
  
)

all_resultsA3 <- purrr::pmap(scenariosA3, function(N_flu, p2, Ntot) {
  sim_study(
    nsim = 10000, 
    Ntot = Ntot, 
    N_flu =N_flu,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=1234567
  )
})

# combine summaries with scenario labels
for (i in seq_len(nrow(scenariosA3))) {
  all_resultsA3[[i]]$scenario <- scenariosA3[i,]
}

save(all_resultsA3, file = "sim_results_A3_arXiv.RData")


#-----------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------Table (B2) (Ntot =500).----------------------------

scenariosA4<- expand.grid(
  N_flu = c(25, 50, 125, 250),
  p2 = c(0.05, 0.1, 0.25, 0.5),
  Ntot=c(500)
  
)

all_resultsA4 <- purrr::pmap(scenariosA4, function(N_flu, p2, Ntot) {
  sim_study(
    nsim = 10000, 
    Ntot = Ntot, 
    N_flu =N_flu,
    psymp_flu = 0.6, psymp_healthy = 0.1,
    p1_symp = 0.5, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=12345678
  )
})

# combine summaries with scenario labels
for (i in seq_len(nrow(scenariosA4))) {
  all_resultsA4[[i]]$scenario <- scenariosA4[i,]
}

save(all_resultsA4, file = "sim_results_A4_arXiv.RData")



#-----------------------------------------------------------------------------------------------------------------------
#-------------------------------------------Table B3 (Simulating Under Conditions Similar to CRISP) (Ntot =1029)-----------

scenariosCRISP<- expand.grid(
  N_flu = c(156),
  p2 = c(0.194),  ##200/1029
  psymp_flu=c(0.25, 0.5, 0.75, 0.9),
  p1_symp=c(0.5, 0.75, 0.9),
  Ntot=c(1029)
  
)

all_resultsCRISP <- purrr::pmap(scenariosCRISP, function(N_flu, p2, psymp_flu, p1_symp, Ntot) {
  sim_study(
    nsim = 10000, 
    Ntot = Ntot, 
    N_flu =N_flu,
    psymp_flu = psymp_flu, psymp_healthy = 0.1,
    p1_symp = p1_symp, p1_nonsymp = 0.2,
    p2 = p2, alpha = 0.05, postdraws = 10000, seed=888
  )
})

# combine summaries with scenario labels
for (i in seq_len(nrow(scenariosCRISP))) {
  all_resultsCRISP[[i]]$scenario <- scenariosCRISP[i,]
}

save(all_resultsCRISP, file = "sim_results_CRISP_arXiv.RData")

#-----------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------Making LaTeX Tables-----------------------------------------------
# helper to get the first matching row (or NULL)
get_row <- function(df, cond) {
  row <- df[cond, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  row[1, , drop = FALSE]
}

# Build a single cell (left = Mean (SD) \\ {[SEs]}, right = CI block)
make_latex_cell <- function(res, method,
                            ci_order = NULL,   # vector of ci_type in the exact order we want
                            se_types = NULL) { # vector of se_type in the exact order we want
  est <- get_row(res$est_summary, res$est_summary$method == method)
  se  <- res$se_summary
  ci  <- res$ci_summary
  
  # --- left: Mean (SD) \\ { [SE], [SE] } ---
  if (is.null(est)) {
    left <- "--"
  } else {
    if (!is.null(se_types)) {
      se_vals <- sapply(se_types, function(st) {
        v <- se$mean_se[se$method == method & se$se_type == st]
        if (length(v) == 0) NA_real_ else v[1]
      })
      se_str <- paste(sprintf("{[%.1f]}", se_vals), collapse = ", ")
      # Note: "\\\\" in the string literal becomes LaTeX's "\\" (linebreak)
      left <- sprintf("%.1f (%.1f) \\\\ %s", est$mean_est, est$sd_est, se_str)
    } else {
      se_row <- se[se$method == method, ]
      se_val <- if (nrow(se_row)) se_row$mean_se[1] else NA
      left <- sprintf("%.1f (%.1f) \\\\ {[%.1f]}", est$mean_est, est$sd_est, se_val)
    }
  }
  
  # --- right: CI block ---
  if (!is.null(ci_order) && length(ci_order) >= 1) {
    # pull each requested ci_type in order (fill missing with NA)
    ci_entries <- lapply(ci_order, function(t) {
      row <- ci[ci$ci_type == t, ]
      if (nrow(row) == 0) return(list(coverage = NA_real_, avg_width = NA_real_))
      list(coverage = row$coverage[1], avg_width = row$avg_width[1])
    })
    cov <- sapply(ci_entries, function(x) if (is.na(x$coverage)) "" else sprintf("%.3f", x$coverage))
    wid <- sapply(ci_entries, function(x) if (is.na(x$avg_width)) "{[]}" else sprintf("{[%.1f]}", x$avg_width))
    
    # Build the LaTeX mini-table depending on how many CI types we want.
    if (length(ci_order) == 4) {
      # layout: top row = var1 & var2, second row = [wid1] & [wid2],
      #         third row = var3 & var4, fourth row = [wid3] & [wid4]
      right <- sprintf("\\begin{tabular}{cc}\n%s & %s \\\\\n%s & %s \\\\\n%s & %s \\\\\n%s & %s\n\\end{tabular}",
                       cov[1], cov[2], wid[1], wid[2], cov[3], cov[4], wid[3], wid[4])
    } else if (length(ci_order) == 3) {
      # layout: top row var1 & var2; next row wid1 & wid2; next: blank & var3; next: blank & wid3
      right <- sprintf("\\begin{tabular}{cc}\n%s & %s \\\\\n%s & %s \\\\\n & %s \\\\\n & %s\n\\end{tabular}",
                       cov[1], cov[2], wid[1], wid[2], cov[3], wid[3])
    } else { # single CI
      right <- sprintf("%s \\\\ %s", cov[1], wid[1])
    }
  } else {
    right <- ""
  }
  
  list(left = left, right = right)
}

# Build mapping matrix from list-index -> (N_index, p2_index)
build_index_map <- function(n_scenarios, N_len = 4, P_len = 4) {
  map <- matrix(NA_integer_, nrow = N_len, ncol = P_len)
  for (i in seq_len(n_scenarios)) {
    k <- (i-1)%%4+1
    N_idx  <- k
    P_idx  <- ceiling(i/4)
    map[N_idx, P_idx] <- i
  }
  map
}

# ---------------------------
# Main: create full table
# ---------------------------
make_full_table <- function(all_results1, write_to = NULL) {
  N_vals  <- c(500, 1000, 2500, 5000)
  p2_vals <- c(0.05, 0.1, 0.25, 0.5)
  n_scen <- length(all_results1)
  
  # build index map according to the diagonal indexing
  map_idx <- build_index_map(n_scen, N_len = length(N_vals), P_len = length(p2_vals))
  
  out_lines <- character(0)
  out_lines <- c(out_lines, "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}")
  out_lines <- c(out_lines, "\\hline")
  out_lines <- c(out_lines, "& \\multicolumn{2}{|c|}{$p_2=0.05$} & \\multicolumn{2}{|c|}{$p_2=0.1$} & \\multicolumn{2}{|c|}{$p_2=0.25$} & \\multicolumn{2}{|c|}{$p_2=0.5$} \\\\")
  out_lines <- c(out_lines, "\\hline")
  out_lines <- c(out_lines, "Estimator & \\makecell{Mean (SD) \\\\ {[avg. SE]}} & \\makecell{CI coverage \\\\ {[avg. width]}} & \\makecell{Mean (SD) \\\\ {[avg. SE]}} & \\makecell{CI coverage \\\\ {[avg. width]}} & \\makecell{Mean (SD) \\\\ {[avg. SE]}} & \\makecell{CI coverage \\\\ {[avg. width]}} & \\makecell{Mean (SD) \\\\ {[avg. SE]}} & \\makecell{CI coverage \\\\ {[avg. width]}} \\\\")
  out_lines <- c(out_lines, "\\hline")
  
  for (r in seq_along(N_vals)) {
    out_lines <- c(out_lines, sprintf("\\multicolumn{9}{|c|}{$N=%d$}\\\\", N_vals[r]))
    out_lines <- c(out_lines, "\\hline")
    
    # --- five-cell row: we will collect 4 cells (one per p2 in order) ---
    five_cells <- lapply(seq_along(p2_vals), function(c) {
      idx <- map_idx[r, c]
      if (is.na(idx)) {
        list(left="--", right="--")
      } else {
        make_latex_cell(all_results1[[idx]],
                        method = "five_cell_MLE",
                        # ordering for the 2x2: Wald unadj | Wald FPC2
                        #                      Bayes unadj | Bayes FPC2
                        ci_order = c("wald_five_cell_MLE_unadj", "wald_five_cell_MLE_FPC1",
                                     "bayes_five_cell_MLE_unadj", "bayes_five_cell_MLE_FPC1"),
                        se_types = c("five_cell_MLE_se_unadj", "five_cell_MLE_se_FPC1"))
      }
    })
    row_five <- paste(sapply(five_cells, function(x) sprintf("\\makecell{%s} & %s", x$left, x$right)), collapse = " & ")
    out_lines <- c(out_lines, paste0("$\\hat{N}_5$ & ", row_five, " \\\\"))
    
    out_lines<-c(out_lines, "\\hline")
    
    
    # --- RS row ---
    rs_cells <- lapply(seq_along(p2_vals), function(c) {
      idx <- map_idx[r, c]
      if (is.na(idx)) {
        list(left="--", right="--")
      } else {
        make_latex_cell(all_results1[[idx]],
                        method = "rs_MLE",
                        ci_order = c("wald_rs_unadj", "wald_rs_FPC", "bayes_rs"),
                        se_types = c("rs_MLE_se_unadj", "rs_MLE_se_FPC"))
      }
    })
    row_rs <- paste(sapply(rs_cells, function(x) sprintf("\\makecell{%s} & %s", x$left, x$right)), collapse = " & ")
    out_lines <- c(out_lines, paste0("$\\hat{N}_{RS}$ & ", row_rs, " \\\\"))
    
    out_lines<-c(out_lines, "\\hline")
    
    # --- Chapman row ---
    chap_cells <- lapply(seq_along(p2_vals), function(c) {
      idx <- map_idx[r, c]
      if (is.na(idx)) {
        list(left="--", right="--")
      } else {
        make_latex_cell(all_results1[[idx]],
                        method = "chapman",
                        ci_order = c("chapman_ci"),
                        se_types = c("chapman_se"))
      }
    })
    row_chap <- paste(sapply(chap_cells, function(x) sprintf("\\makecell{%s} & \\makecell{%s}", x$left, x$right)), collapse = " & ")
    out_lines <- c(out_lines, paste0("$\\hat{N}_{Chap}$ & ", row_chap, " \\\\"))
    
    out_lines <- c(out_lines, "\\hline")
  }
  
  out_lines <- c(out_lines, "\\end{tabular}")
  
  # print to console
  cat(paste(out_lines, collapse = "\n"), "\n")
  
  # optionally write to a file for \input{}
  if (!is.null(write_to)) writeLines(out_lines, con = write_to)
  
  invisible(out_lines)
}




