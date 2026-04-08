###############################################
## functions_application.R
## Updated: Delta method (Δ = S - B) for npcdens
## FAST VERSION with normal-reference bandwidth
## AUTO-DETECTS when to exclude X1
## X3 REMOVED (constant in filtered data)
##
## USER REQUEST UPDATE:
##  - Whenever estimating something | S, B, ...  ==>  use | delta, ...
##    where delta = S - B
##  - Added "# CHANGED:" comments at modifications
###############################################

library(SuperLearner)
library(np)

expit <- function(x) exp(x)/(1+exp(x))

###############################################
## Smoothed trimming helper functions
###############################################

getSt = function(piHat, epsilon = 1e-2, t = 0.1){
  pnorm(piHat - t, mean = 0, sd = epsilon)
}

getSt.deriv = function(piHat, epsilon = 1e-2, t = 0.1){
  dnorm(piHat - t, mean = 0, sd = epsilon)
}

gaussianKernel = function(S, s, h){
  u = S - s
  (1/sqrt(2*pi*h^2)) * exp(-0.5*(u/h)^2)
}


###############################################################
## Actual Risk in Identified Subpopulation (General Version)
## Target: E[Y(a)] for a in {0, 1}
###############################################################
est.actual.risk.smoothed.function <- 
  function(data, a_target = 1, s = 3, t = 0.1, epsilon = 1e-2, h = 0.8,
           s0.seq = s0.seq, variance = TRUE) {
    
    n <- nrow(data)
    nsplits <- 2
    s0.len <- s0.seq[2] - s0.seq[1]
    
    num_i <- numeric(n)
    den_i <- numeric(n)
    
    # Pre-calculate kernel for the grid (used in the denominator integral)
    kernel.s0 <- gaussianKernel(s0.seq, s, h)
    kernel.s0 <- kernel.s0 / sum(kernel.s0 * s0.len) # Normalize to integrate to 1
    
    sl.lib <- c("SL.glm", "SL.mean", "SL.ranger")
    
    include_X1 <- length(unique(data$X1)) > 1
    X_vars_noA_delta_only <- if (include_X1) c("X1","X2") else c("X2")
    X_vars_dens <- c("A", "S", X_vars_noA_delta_only)
    
    set.seed(123)
    splitIndex <- sample(rep(1:nsplits, length.out = n))
    
    for (split in 1:nsplits) {
      train <- splitIndex != split
      test  <- splitIndex == split
      
      data.train <- data[train, ]
      data.test  <- data[test, ]
      n_test <- sum(test)
      
      data.train$delta_S <- data.train$S - data.train$B
      data.test$delta_S  <- data.test$S  - data.test$B
      
      ## 1) Treatment model (IPW)
      X_ps_train <- data.train[, c(X_vars_noA_delta_only, "S", "delta_S")]
      X_ps_test  <- data.test[,  c(X_vars_noA_delta_only, "S", "delta_S")]
      
      ps.fit <- SuperLearner(
        Y = data.train$A, X = X_ps_train, newX = X_ps_test,
        SL.library = sl.lib, family = binomial()
      )
      pA1_test <- as.numeric(ps.fit$SL.predict)
      pA_test <- if (a_target == 1) pA1_test else (1 - pA1_test)
      pA_test <- pmax(pA_test, 1e-6)
      # wA_test ensures we only look at people who actually took a_target
      wA_test <- as.numeric(data.test$A == a_target) / pA_test
      
      ## 2) Density model (piHat)
      bw.fit.delta <- npcdensbw(
        ydat = data.train$delta_S,
        xdat = data.train[, X_vars_dens],
        bwmethod = "normal-reference"
      )
      
      # A) For Numerator: Pi at observed S_i (given A = a_target)
      data.test_as_target <- data.test
      data.test_as_target$A <- a_target
      piHat.obs <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = data.test_as_target[, X_vars_dens],
        eydat = data.test$delta_S
      ))
      
      # B) For Denominator: Pi over the s0 grid (the integral)
      S_obs_expanded <- rep(data.test$S, each = length(s0.seq))
      delta_obs_expanded <- rep(data.test$delta_S, each = length(s0.seq))
      s0_expanded <- rep(s0.seq, times = n_test)
      delta_s0 <- (s0_expanded - S_obs_expanded) + delta_obs_expanded
      
      valid_idx <- delta_s0 >= 0
      piHat.s0.vec <- rep(1e-10, length(delta_s0))
      
      if (sum(valid_idx) > 0) {
        data_s0_expanded <- data.test[rep(seq_len(n_test), each = length(s0.seq)), ]
        data_s0_expanded$A <- a_target
        piHat.s0.vec[valid_idx] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s0_expanded[valid_idx, X_vars_dens],
          eydat = delta_s0[valid_idx]
        ))
      }
      piHat.s0.test <- matrix(piHat.s0.vec, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      
      ## 3) Calculate components for Equation (4)
      
      # Numerator logic: wA * K_h(S_i - s) * phi(pi_obs) * Y_i
      # We use the same kernel logic for the observed point as the grid
      kernel.obs <- dnorm(data.test$S, mean = s, sd = h) 
      # Note: The denominator integral uses a normalized kernel, 
      # so we ensure kernel.obs is scaled the same way if necessary.
      
      phi_obs <- getSt(piHat.obs, epsilon, t)
      num_i[test] <- wA_test * kernel.obs * phi_obs * data.test$Y
      
      # Denominator logic: Integral of [K_h(s'-s) * phi(pi(s'))] ds'
      # This is averaged over the distribution of B, X (the test sample)
      St.s0.test <- getSt(piHat.s0.test, epsilon, t)
      den_i[test] <- s0.len * rowSums(
        St.s0.test * matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      )
    }
    
    # STWCR is the ratio of the expected numerator and expected denominator
    risk_hat <- mean(num_i) / mean(den_i)
    
    ## Influence-function based variance for the ratio of means
    # phi_ratio = (num - ratio * den) / E[den]
    inf_fn <- (num_i - risk_hat * den_i) / mean(den_i)
    var_hat <- var(inf_fn) / n
    
    list(
      actual_risk.est = risk_hat,
      CI_lower  = risk_hat - qnorm(0.975) * sqrt(var_hat),
      CI_upper  = risk_hat + qnorm(0.975) * sqrt(var_hat)
    )
  }


###############################################################
## MAIN: Smoothed CoR with SuperLearner + 2-fold cross-fitting
## DELTA METHOD: models π(Δ | A,B,X) where Δ = S - B
## AUTO-DETECTS X1 variance and excludes if constant
###############################################################

est.tau.smoothed.function <- 
  function(data, s = 3, t = 0.1, epsilon = 1e-2, h = 0.8,
           s0.seq = s0.seq, variance = TRUE) {
    
    n <- nrow(data)
    nsplits <- 2
    s0.len <- s0.seq[2] - s0.seq[1]
    
    psihat.num.term  <- numeric(n)
    psihat.deno.term <- numeric(n)
    
    kernel.s0 <- gaussianKernel(s0.seq, s, h)
    kernel.s0 <- kernel.s0 / sum(kernel.s0 * s0.len)
    
    sl.lib <- c("SL.glm", "SL.mean", "SL.ranger")
    
    # AUTO-DETECT: Should X1 be included?
    include_X1 <- length(unique(data$X1)) > 1
    
    if (include_X1) {
      X_vars <- c("A","B","X1","X2")
      X_vars_noA <- c("B","X1","X2")
      
      # CHANGED: versions that replace (S,B) by delta_S when both appear in conditioning set
      X_vars_delta_only   <- c("A","X1","X2")
      X_vars_noA_delta_only <- c("X1","X2")
      
      # CHANGED2: density model will condition on (A, S, X) without B
      X_vars_dens <- c("A","S", X_vars_noA_delta_only)
    } else {
      cat("Note: X1 is constant in this dataset, excluding from models\n")
      X_vars <- c("A","B","X2")
      X_vars_noA <- c("B","X2")
      
      # CHANGED: versions that replace (S,B) by delta_S when both appear in conditioning set
      X_vars_delta_only   <- c("A","X2")
      X_vars_noA_delta_only <- c("X2")
      
      # CHANGED2: density model will condition on (A, S, X) without B
      X_vars_dens <- c("A","S", X_vars_noA_delta_only)
    }
    
    set.seed(123)
    splitIndex <- sample(rep(1:nsplits, length.out = n))
    
    for (split in 1:nsplits) {
      
      train <- splitIndex != split
      test  <- splitIndex == split
      
      data.train <- data[train, ]
      data.test  <- data[test, ]
      n_test <- sum(test)
      
      # CHANGED: define delta_S on BOTH train/test (used in models replacing |S,B,... with |delta,...)
      data.train$delta_S <- data.train$S - data.train$B
      data.test$delta_S  <- data.test$S  - data.test$B
      
      ########################################
      ## 1) Treatment model
      ##    was: A | (B, X, S)
      ## CHANGED: A | (X, delta_S)
      ########################################
      X_ps_train <- data.train[, c(X_vars_noA_delta_only, "S", "delta_S")]  # CHANGED2
      X_ps_test  <- data.test[,  c(X_vars_noA_delta_only, "S", "delta_S")]  # CHANGED2
      
      ps.fit <- SuperLearner(
        Y = data.train$A,
        X = X_ps_train,
        newX = X_ps_test,
        SL.library = sl.lib,
        family = binomial()
      )
      est_p_a_test <- as.numeric(ps.fit$SL.predict)
      
      
      ########################################
      ## 2) Conditional mean of S | A,B,X (kept for reference)
      ## (Not of the form |S,B,... so unchanged)
      ########################################
      X_piS_train <- data.train[, X_vars]
      X_piS_test  <- data.test[,  X_vars]
      
      piS.fit <- SuperLearner(
        Y = data.train$S,
        X = X_piS_train,
        newX = X_piS_test,
        SL.library = sl.lib
      )
      
      meanS.test <- as.numeric(piS.fit$SL.predict)
      meanS.train <- as.numeric(predict(piS.fit, newdata = X_piS_train)$pred)
      sdS.est <- sd(data.train$S - meanS.train)
      
      
      ########################################
      ## 3) Outcome model
      ##    was: Y | (A, B, X, S)
      ## CHANGED: Y | (A, X, delta_S)
      ########################################
      X_r_train <- data.train[, c(X_vars_delta_only, "S", "delta_S")]  # CHANGED2
      X_r_test  <- data.test[,  c(X_vars_delta_only, "S", "delta_S")]  # CHANGED2
      
      r.fit <- SuperLearner(
        Y = data.train$Y,
        X = X_r_train,
        newX = X_r_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      rHat.S.test <- as.numeric(r.fit$SL.predict)
      
      
      ########################################
      ## 4) r(a, s0.grid, X)
      ##    was: vary S on grid
      ## CHANGED: vary delta_S on grid using delta_s0 = s0 - B
      ########################################
      data_s0 <- data.test[rep(seq_len(n_test), each = length(s0.seq)), ]
      # CHANGED2: store observed S and delta so we can hold B fixed without explicitly using B
      S_obs_rep <- data_s0$S
      delta_obs_rep <- data_s0$delta_S
      data_s0$S <- rep(s0.seq, times = n_test)  # keep S for kernel usage below
      # CHANGED2: hold B fixed via delta_s0 = s0 - B = (s0 - S_obs) + delta_obs
      data_s0$delta_S <- (data_s0$S - S_obs_rep) + delta_obs_rep
      
      X_r_s0 <- data_s0[, colnames(X_r_train)]  # CHANGED (X_r_train now uses delta_S)
      rHat.s0.test <- matrix(
        as.numeric(predict(r.fit, newdata = X_r_s0)$pred),
        nrow = n_test,
        ncol = length(s0.seq),
        byrow = TRUE
      )
      
      
      ########################################
      ## 5) π(S | X,A) using DELTA METHOD
      ##    Models π(Δ | A,B,X) where Δ = S - B
      ## (Already delta-based; keep conditioning xdat as in your original delta method)
      ########################################
      
      # Create delta = S - B for training data (already created above)
      # data.train$delta_S <- data.train$S - data.train$B
      
      # Fit conditional density on Δ
      bw.fit.delta <- npcdensbw(
        ydat = data.train$delta_S,
        xdat = data.train[, X_vars_dens],  # CHANGED2,            # keep π(Δ | A,B,X)
        bwmethod = "normal-reference"           # Fast bandwidth selection
      )
      
      ## π(S_obs | X) via π(Δ_obs | X)
      delta_obs <- data.test$delta_S
      
      piHat.S.test <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = data.test[, X_vars_dens],  # CHANGED2,
        eydat = delta_obs
      ))
      piHat.S.test <- pmax(piHat.S.test, 1e-10)
      
      
      ## π(s0 | X) via π(delta_s0 | X) where delta_s0 = s0 - B
      S_obs_expanded <- rep(data.test$S, each = length(s0.seq))  # CHANGED2
      delta_obs_expanded <- rep(data.test$delta_S, each = length(s0.seq))  # CHANGED2
      s0_expanded <- rep(s0.seq, times = n_test)
      delta_s0 <- (s0_expanded - S_obs_expanded) + delta_obs_expanded  # CHANGED2
      
      valid_idx <- delta_s0 >= 0
      piHat.s0.vec <- rep(1e-10, length(delta_s0))
      
      if (sum(valid_idx) > 0) {
        data_s0_expanded <- data.test[rep(seq_len(n_test), each = length(s0.seq)), ]
        
        piHat.s0.vec[valid_idx] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s0_expanded[valid_idx, X_vars_dens],  # CHANGED2,
          eydat = delta_s0[valid_idx]
        ))
      }
      
      piHat.s0.test <- matrix(piHat.s0.vec, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      piHat.s0.test <- pmax(piHat.s0.test, 1e-10)
      
      
      ########################################
      ## 6–8 EIF terms (unchanged)
      ########################################
      
      St.S.test        <- getSt(piHat.S.test, epsilon, t)
      St.deriv.S.test  <- getSt.deriv(piHat.S.test, epsilon, t)
      St.s0.test       <- getSt(piHat.s0.test, epsilon, t)
      St.deriv.s0.test <- getSt.deriv(piHat.s0.test, epsilon, t)
      
      kernel.S.test <- gaussianKernel(data.test$S, s, h)
      kernel.S.test <- kernel.S.test / sum(kernel.s0 * s0.len)
      
      ## NUMERATOR
      term1 <- kernel.S.test * (data.test$A / est_p_a_test) *
        St.deriv.S.test * rHat.S.test
      
      term2 <- kernel.S.test * (data.test$A / est_p_a_test) *
        (St.S.test / piHat.S.test) *
        (data.test$Y - rHat.S.test)
      
      term3 <- -s0.len * rowSums(
        (rHat.s0.test * piHat.s0.test * St.deriv.s0.test) *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      ) * (data.test$A / est_p_a_test)
      
      term4 <- s0.len * rowSums(
        (rHat.s0.test * St.s0.test) *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      )
      
      psihat.num.term[test] <- term1 + term2 + term3 + term4
      
      
      ## DENOMINATOR
      den1 <- kernel.S.test * (data.test$A / est_p_a_test) * St.deriv.S.test
      
      den2 <- -s0.len * rowSums(
        (piHat.s0.test * St.deriv.s0.test) *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      ) * (data.test$A / est_p_a_test)
      
      den3 <- s0.len * rowSums(
        St.s0.test *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      )
      
      psihat.deno.term[test] <- den1 + den2 + den3
      
    } # split
    
    tau_num <- mean(psihat.num.term)
    tau_den <- mean(psihat.deno.term)
    tau     <- tau_num / tau_den
    
    var_hat <- var((psihat.num.term - tau * psihat.deno.term) / tau_den) / n
    
    list(
      tau.h.est = tau,
      CI_lower  = tau - qnorm(0.975) * sqrt(var_hat),
      CI_upper  = tau + qnorm(0.975) * sqrt(var_hat)
    )
  }

###############################################################
## CVE Function with DELTA METHOD
## AUTO-DETECTS X1 variance and excludes if constant
###############################################################

est.cve.eif.smooth.function <- function(data, s0, s1,
                                        a0 = 0, a1 = 1,
                                        t = 0.1, h = 0.5,
                                        s0.seq = s0.seq, s1.seq = s1.seq,
                                        epsilon = 1e-2) {
  
  n <- nrow(data)
  nsplits <- 2
  set.seed(123)
  splitIndex <- sample(rep(1:nsplits, length.out = n))
  
  sl.lib <- c("SL.glm", "SL.mean", "SL.ranger")
  
  # AUTO-DETECT: Should X1 be included?
  include_X1 <- length(unique(data$X1)) > 1
  
  if (include_X1) {
    X_vars <- c("A","B","X1","X2")
    X_vars_noA <- c("B","X1","X2")
    
    # CHANGED: versions that replace (S,B) by delta_S when both appear in conditioning set
    X_vars_delta_only     <- c("A","X1","X2")
    X_vars_noA_delta_only <- c("X1","X2")
    
    # CHANGED2: density model will condition on (A, S, X) without B
    X_vars_dens <- c("A","S", X_vars_noA_delta_only)
  } else {
    X_vars <- c("A","B","X2")
    X_vars_noA <- c("B","X2")
    
    # CHANGED: versions that replace (S,B) by delta_S when both appear in conditioning set
    X_vars_delta_only     <- c("A","X2")
    X_vars_noA_delta_only <- c("X2")
    
    # CHANGED2: density model will condition on (A, S, X) without B
    X_vars_dens <- c("A","S", X_vars_noA_delta_only)
  }
  
  ############################################################
  ## Helper: compute f(a0,s0), g(a1,s1), ψ_i via cross-fitting
  ############################################################
  get_den_with_eif_cf <- function(data, s0, s1, a0, a1) {
    
    n <- nrow(data)
    s0.len <- s0.seq[2] - s0.seq[1]
    s1.len <- s1.seq[2] - s1.seq[1]
    
    kernel.s0 <- gaussianKernel(s0.seq, s0, h)
    kernel.s0 <- kernel.s0 / sum(kernel.s0 * s0.len)
    
    kernel.s1 <- gaussianKernel(s1.seq, s1, h)
    kernel.s1 <- kernel.s1 / sum(kernel.s1 * s1.len)
    
    f_vec <- numeric(n)
    g_vec <- numeric(n)
    psi_vec <- numeric(n)
    
    for (split in 1:nsplits) {
      
      train <- splitIndex != split
      test  <- splitIndex == split
      
      dt_train <- data[train, ]
      dt_test  <- data[test, ]
      n_test <- sum(test)
      
      # CHANGED: define delta_S on BOTH train/test (used in models replacing |S,B,... with |delta,...)
      dt_train$delta_S <- dt_train$S - dt_train$B
      dt_test$delta_S  <- dt_test$S  - dt_test$B
      
      #########################
      ## Treatment model
      ## was: A | (B, X, S)
      ## CHANGED: A | (X, delta_S)
      #########################
      X_ps_train <- dt_train[, c(X_vars_noA_delta_only, "S", "delta_S")]  # CHANGED2
      X_ps_test  <- dt_test[,  c(X_vars_noA_delta_only, "S", "delta_S")]  # CHANGED2
      
      psMod <- SuperLearner(
        Y = dt_train$A,
        X = X_ps_train,
        newX = X_ps_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      pA1_test <- as.numeric(psMod$SL.predict)
      pA0_test <- 1 - pA1_test
      
      #########################
      ## S model (kept for reference)
      ## (Not of the form |S,B,... so unchanged)
      #########################
      X_S_train <- dt_train[, X_vars]
      X_S_test  <- dt_test[,  X_vars]
      
      SMod <- SuperLearner(
        Y = dt_train$S,
        X = X_S_train,
        newX = X_S_test,
        SL.library = sl.lib
      )
      
      meanS_test <- as.numeric(SMod$SL.predict)
      meanS_train <- as.numeric(predict(SMod, newdata = X_S_train)$pred)
      sdS.est <- sd(dt_train$S - meanS_train)
      
      # counterfactual conditioning variables
      X_S_test_a0 <- X_S_test; X_S_test_a0$A <- a0
      X_S_test_a1 <- X_S_test; X_S_test_a1$A <- a1
      
      meanS_a0_test <- as.numeric(predict(SMod, newdata = X_S_test_a0)$pred)
      meanS_a1_test <- as.numeric(predict(SMod, newdata = X_S_test_a1)$pred)
      
      #########################
      ## Outcome model
      ## was: Y | (A, B, X, S)
      ## CHANGED: Y | (A, X, delta_S)
      #########################
      X_Y_train <- dt_train[, c(X_vars_delta_only, "S", "delta_S")]  # CHANGED2
      X_Y_test  <- dt_test[,  c(X_vars_delta_only, "S", "delta_S")]  # CHANGED2
      
      YMod <- SuperLearner(
        Y = dt_train$Y,
        X = X_Y_train,
        newX = X_Y_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      # r(a0, S_obs, X)  (implemented via observed delta_S)
      X_Y_test_a0 <- X_Y_test; X_Y_test_a0$A <- a0
      r_S_a0 <- as.numeric(predict(YMod, newdata = X_Y_test_a0)$pred)
      
      # r(a0, s0.grid, X): CHANGED to vary delta_S = s0 - B
      dt_s0 <- dt_test[rep(seq_len(n_test), each = length(s0.seq)), ]
      dt_s0$S <- rep(s0.seq, times = n_test)      # keep S for kernels
      dt_s0$A <- a0
      dt_s0$delta_S <- dt_s0$S - dt_s0$B          # CHANGED
      X_s0 <- dt_s0[, colnames(X_Y_train)]        # CHANGED
      r_s0 <- matrix(
        as.numeric(predict(YMod, newdata = X_s0)$pred),
        nrow = n_test,
        ncol = length(s0.seq),
        byrow = TRUE
      )
      
      # r(a1, s1.grid, X): CHANGED to vary delta_S = s1 - B
      dt_s1 <- dt_test[rep(seq_len(n_test), each = length(s1.seq)), ]
      dt_s1$S <- rep(s1.seq, times = n_test)      # keep S for kernels
      dt_s1$A <- a1
      dt_s1$delta_S <- dt_s1$S - dt_s1$B          # CHANGED
      X_s1 <- dt_s1[, colnames(X_Y_train)]        # CHANGED
      r_s1 <- matrix(
        as.numeric(predict(YMod, newdata = X_s1)$pred),
        nrow = n_test,
        ncol = length(s1.seq),
        byrow = TRUE
      )
      
      #########################
      ## π(s | X,A) with DELTA METHOD
      ## (Already delta-based; keep conditioning xdat as in your original delta method)
      #########################
      
      # Fit conditional density model on delta
      bw.fit.delta <- npcdensbw(
        ydat = dt_train$delta_S,
        xdat = dt_train[, X_vars_dens],  # CHANGED2,
        bwmethod = "normal-reference"
      )
      
      ## π(S_obs | X, a0) via delta
      X_test_a0 <- dt_test[, X_vars_dens]  # CHANGED2
      X_test_a0$A <- a0
      delta_obs_a0 <- dt_test$delta_S
      
      pi_S_a0 <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = X_test_a0,
        eydat = delta_obs_a0
      ))
      pi_S_a0 <- pmax(pi_S_a0, 1e-10)
      
      ## π(S_obs | X, a1) via delta
      X_test_a1 <- dt_test[, X_vars_dens]  # CHANGED2
      X_test_a1$A <- a1
      delta_obs_a1 <- dt_test$delta_S
      
      pi_S_a1 <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = X_test_a1,
        eydat = delta_obs_a1
      ))
      pi_S_a1 <- pmax(pi_S_a1, 1e-10)
      
      St_S_a0 <- getSt(pi_S_a0, epsilon, t)
      St_S_a1 <- getSt(pi_S_a1, epsilon, t)
      
      ## π(s0 | X, a0) on grid via delta
      S_obs_expanded <- rep(dt_test$S, each = length(s0.seq))  # CHANGED2
      delta_obs_expanded <- rep(dt_test$delta_S, each = length(s0.seq))  # CHANGED2
      s0_expanded <- rep(s0.seq, times = n_test)
      delta_s0 <- (s0_expanded - S_obs_expanded) + delta_obs_expanded  # CHANGED2
      valid_s0 <- delta_s0 >= 0
      
      pi_s0_a0_vec <- rep(1e-10, length(delta_s0))
      if (sum(valid_s0) > 0) {
        data_s0_expanded <- dt_test[rep(seq_len(n_test), each = length(s0.seq)), ]
        data_s0_expanded$A <- a0
        data_s0_expanded$S <- s0_expanded  # CHANGED2: needed because X_vars_dens includes S
        
        pi_s0_a0_vec[valid_s0] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s0_expanded[valid_s0, X_vars_dens],  # CHANGED2,
          eydat = delta_s0[valid_s0]
        ))
      }
      pi_s0_a0 <- matrix(pi_s0_a0_vec, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      pi_s0_a0 <- pmax(pi_s0_a0, 1e-10)
      
      ## π(s1 | X, a1) on grid via delta
      S_obs_expanded <- rep(dt_test$S, each = length(s1.seq))  # CHANGED2
      delta_obs_expanded <- rep(dt_test$delta_S, each = length(s1.seq))  # CHANGED2
      s1_expanded <- rep(s1.seq, times = n_test)
      delta_s1 <- (s1_expanded - S_obs_expanded) + delta_obs_expanded  # CHANGED2
      valid_s1 <- delta_s1 >= 0
      
      pi_s1_a1_vec <- rep(1e-10, length(delta_s1))
      if (sum(valid_s1) > 0) {
        data_s1_expanded <- dt_test[rep(seq_len(n_test), each = length(s1.seq)), ]
        data_s1_expanded$A <- a1
        data_s1_expanded$S <- s1_expanded  # CHANGED2: needed because X_vars_dens includes S
        
        pi_s1_a1_vec[valid_s1] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s1_expanded[valid_s1, X_vars_dens],  # CHANGED2,
          eydat = delta_s1[valid_s1]
        ))
      }
      pi_s1_a1 <- matrix(pi_s1_a1_vec, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      pi_s1_a1 <- pmax(pi_s1_a1, 1e-10)
      
      St_s0_a0       <- getSt(pi_s0_a0, epsilon, t)
      St_s1_a1       <- getSt(pi_s1_a1, epsilon, t)
      St_der_s0_a0   <- getSt.deriv(pi_s0_a0, epsilon, t)
      St_der_s1_a1   <- getSt.deriv(pi_s1_a1, epsilon, t)
      
      #########################
      ## Kernels at observed S (estimand still indexed by S, unchanged)
      #########################
      kern_s0_obs <- gaussianKernel(dt_test$S, s0, h)
      kern_s0_obs <- kern_s0_obs / sum(kernel.s0 * s0.len)
      
      kern_s1_obs <- gaussianKernel(dt_test$S, s1, h)
      kern_s1_obs <- kern_s1_obs / sum(kernel.s1 * s1.len)
      
      #########################
      ## Integrals f(a0,s0), g(a1,s1)
      #########################
      f_test <- s0.len * rowSums(
        (St_s0_a0 * r_s0) * matrix(kernel.s0, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      )
      
      g_test <- s1.len * rowSums(
        (St_s1_a1 * r_s1) * matrix(kernel.s1, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      )
      
      #########################
      ## IPW residuals
      #########################
      ipwRes <- (dt_test$Y - r_S_a0) / pi_S_a0
      
      #########################
      ## EIF for f(a0,s0)
      #########################
      Df1 <- kern_s0_obs * ((dt_test$A == a0) / pA0_test) *
        St_S_a0 * r_S_a0
      
      Df2 <- kern_s0_obs * ((dt_test$A == a0) / pA0_test) *
        St_S_a0 * ipwRes
      
      Df3 <- -s0.len * rowSums(
        (r_s0 * pi_s0_a0 * St_der_s0_a0) *
          matrix(kernel.s0, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      ) * ((dt_test$A == a0) / pA0_test)
      
      Df4 <- s0.len * rowSums(
        (r_s0 * St_s0_a0) *
          matrix(kernel.s0, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      )
      
      Df <- Df1 + Df2 + Df3 + Df4
      
      
      #########################
      ## EIF for g(a1,s1)
      #########################
      Dg1 <- kern_s1_obs * ((dt_test$A == a1) / pA1_test) * St_S_a1
      
      Dg2 <- -s1.len * rowSums(
        (pi_s1_a1 * St_der_s1_a1) *
          matrix(kernel.s1, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      ) * ((dt_test$A == a1) / pA1_test)
      
      Dg3 <- s1.len * rowSums(
        (St_s1_a1) *
          matrix(kernel.s1, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      )
      
      Dg <- Dg1 + Dg2 + Dg3
      
      
      #########################
      ## EIF ψ_i
      #########################
      psi <- f_test * Dg + g_test * Df - f_test * g_test
      
      #########################
      ## store
      #########################
      f_vec[test] <- f_test
      g_vec[test] <- g_test
      psi_vec[test] <- psi
      
    } # end split
    
    list(
      f_vec = f_vec,
      g_vec = g_vec,
      psi_vec = psi_vec,
      est = mean(psi_vec),
      true = mean(f_vec * g_vec)
    )
    
  } ## end helper
  
  
  ###############################################################
  ## Numerator (swap s0/s1, a0/a1)
  ###############################################################
  den <- get_den_with_eif_cf(data, s0, s1, a0, a1)
  num <- get_den_with_eif_cf(data, s1, s0, a1, a0)
  
  ratio_hat <- num$est / den$est
  cve_hat <- 1 - ratio_hat

  psi_cve <- -(num$psi_vec - ratio_hat * den$psi_vec) / den$est
  var_cve <- var(psi_cve) / n
  se_cve  <- sqrt(var_cve)

  psi_log_ratio <- num$psi_vec / num$est - den$psi_vec / den$est
  if (is.finite(ratio_hat) && ratio_hat > 0) {
    se_log_ratio <- sqrt(var(psi_log_ratio) / n)
    ci_ratio_log <- ratio_hat * exp(c(-1.96, 1.96) * se_log_ratio)
    ci_ratio_log_to_cve <- c(1 - ci_ratio_log[2], 1 - ci_ratio_log[1])
  } else {
    ci_ratio_log_to_cve <- c(NA_real_, NA_real_)
  }
  
  list(
    cve_est  = cve_hat,
    cve_num = num$est,
    cve_den = den$est,
    se       = se_cve,
    ci       = c(cve_hat - 1.96 * se_cve,
                 cve_hat + 1.96 * se_cve),
    ci_ratio_log_to_cve = ci_ratio_log_to_cve,
    cve_true = 1 - num$true / den$true
  )
  
}


est.tau.smoothed.function.a0 <- 
  function(data, s = 3, t = 0.1, epsilon = 1e-2, h = 0.8,
           s0.seq = s0.seq, variance = TRUE) {
    
    n <- nrow(data)
    nsplits <- 2
    s0.len <- s0.seq[2] - s0.seq[1]
    
    psihat.num.term  <- numeric(n)
    psihat.deno.term <- numeric(n)
    
    kernel.s0 <- gaussianKernel(s0.seq, s, h)
    kernel.s0 <- kernel.s0 / sum(kernel.s0 * s0.len)
    
    sl.lib <- c("SL.glm", "SL.mean", "SL.ranger")
    
    include_X1 <- length(unique(data$X1)) > 1
    
    if (include_X1) {
      X_vars <- c("A","B","X1","X2")
      X_vars_noA <- c("B","X1","X2")
      
      X_vars_delta_only        <- c("A","X1","X2")
      X_vars_noA_delta_only    <- c("X1","X2")
      
      X_vars_dens <- c("A","S", X_vars_noA_delta_only)
    } else {
      cat("Note: X1 is constant in this dataset, excluding from models\n")
      X_vars <- c("A","B","X2")
      X_vars_noA <- c("B","X2")
      
      X_vars_delta_only        <- c("A","X2")
      X_vars_noA_delta_only    <- c("X2")
      
      X_vars_dens <- c("A","S", X_vars_noA_delta_only)
    }
    
    set.seed(123)
    splitIndex <- sample(rep(1:nsplits, length.out = n))
    
    for (split in 1:nsplits) {
      
      train <- splitIndex != split
      test  <- splitIndex == split
      
      data.train <- data[train, ]
      data.test  <- data[test, ]
      n_test <- sum(test)
      
      data.train$delta_S <- data.train$S - data.train$B
      data.test$delta_S  <- data.test$S  - data.test$B
      
      ########################################
      ## 1) Treatment model: P(A=1 | X, S, delta)
      ########################################
      X_ps_train <- data.train[, c(X_vars_noA_delta_only, "S", "delta_S")]
      X_ps_test  <- data.test[,  c(X_vars_noA_delta_only, "S", "delta_S")]
      
      ps.fit <- SuperLearner(
        Y = data.train$A,
        X = X_ps_train,
        newX = X_ps_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      pA1_test <- as.numeric(ps.fit$SL.predict)
      pA1_test <- pmax(pmin(pA1_test, 1 - 1e-6), 1e-6)
      pA0_test <- 1 - pA1_test
      
      # ⭐ a=0 IPW weight
      wA_test <- (1 - data.test$A) / pA0_test
      
      ########################################
      ## 2) Conditional mean of S | A,B,X (unchanged)
      ########################################
      X_piS_train <- data.train[, X_vars]
      X_piS_test  <- data.test[,  X_vars]
      
      piS.fit <- SuperLearner(
        Y = data.train$S,
        X = X_piS_train,
        newX = X_piS_test,
        SL.library = sl.lib
      )
      
      meanS.test <- as.numeric(piS.fit$SL.predict)
      meanS.train <- as.numeric(predict(piS.fit, newdata = X_piS_train)$pred)
      sdS.est <- sd(data.train$S - meanS.train)
      
      ########################################
      ## 3) Outcome model: Y | (A, X, S, delta)
      ########################################
      X_r_train <- data.train[, c(X_vars_delta_only, "S", "delta_S")]
      X_r_test  <- data.test[,  c(X_vars_delta_only, "S", "delta_S")]
      
      r.fit <- SuperLearner(
        Y = data.train$Y,
        X = X_r_train,
        newX = X_r_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      # ⭐ Force outcome prediction at A=0
      X_r_test_a0 <- X_r_test
      X_r_test_a0$A <- 0
      rHat.S.test <- as.numeric(predict(r.fit, newdata = X_r_test_a0)$pred)
      
      ########################################
      ## 4) r(0, s0.grid, X): force A=0 on grid
      ########################################
      data_s0 <- data.test[rep(seq_len(n_test), each = length(s0.seq)), ]
      S_obs_rep <- data_s0$S
      delta_obs_rep <- data_s0$delta_S
      
      data_s0$S <- rep(s0.seq, times = n_test)
      data_s0$delta_S <- (data_s0$S - S_obs_rep) + delta_obs_rep
      data_s0$A <- 0
      
      X_r_s0 <- data_s0[, colnames(X_r_train)]
      rHat.s0.test <- matrix(
        as.numeric(predict(r.fit, newdata = X_r_s0)$pred),
        nrow = n_test,
        ncol = length(s0.seq),
        byrow = TRUE
      )
      
      ########################################
      ## 5) π(Δ | A,S,X) (your delta density model)
      ########################################
      bw.fit.delta <- npcdensbw(
        ydat = data.train$delta_S,
        xdat = data.train[, X_vars_dens],
        bwmethod = "normal-reference"
      )
      
      # ⭐ IMPORTANT: for a=0 estimand, density should be evaluated at A=0
      X_dens_test_a0 <- data.test[, X_vars_dens]
      X_dens_test_a0$A <- 0
      
      delta_obs <- data.test$delta_S
      piHat.S.test <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = X_dens_test_a0,
        eydat = delta_obs
      ))
      piHat.S.test <- pmax(piHat.S.test, 0.05)
      
      # grid density π(s0 | X, A=0)
      S_obs_expanded <- rep(data.test$S, each = length(s0.seq))
      delta_obs_expanded <- rep(data.test$delta_S, each = length(s0.seq))
      s0_expanded <- rep(s0.seq, times = n_test)
      delta_s0 <- (s0_expanded - S_obs_expanded) + delta_obs_expanded
      
      valid_idx <- delta_s0 >= 0
      piHat.s0.vec <- rep(1e-10, length(delta_s0))
      
      if (sum(valid_idx) > 0) {
        data_s0_expanded <- data.test[rep(seq_len(n_test), each = length(s0.seq)), ]
        data_s0_expanded$A <- 0
        data_s0_expanded$S <- s0_expanded  # because X_vars_dens includes S
        
        piHat.s0.vec[valid_idx] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s0_expanded[valid_idx, X_vars_dens],
          eydat = delta_s0[valid_idx]
        ))
      }
      
      piHat.s0.test <- matrix(piHat.s0.vec, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      piHat.s0.test <- pmax(piHat.s0.test, 1e-10)
      
      ########################################
      ## 6–8 EIF terms (same form, but use wA_test)
      ########################################
      St.S.test        <- getSt(piHat.S.test, epsilon, t)
      St.deriv.S.test  <- getSt.deriv(piHat.S.test, epsilon, t)
      St.s0.test       <- getSt(piHat.s0.test, epsilon, t)
      St.deriv.s0.test <- getSt.deriv(piHat.s0.test, epsilon, t)
      
      kernel.S.test <- gaussianKernel(data.test$S, s, h)
      kernel.S.test <- kernel.S.test / sum(kernel.s0 * s0.len)
      
      ## NUMERATOR
      term1 <- kernel.S.test * wA_test * St.deriv.S.test * rHat.S.test
      
      term2 <- kernel.S.test * wA_test *
        (St.S.test / piHat.S.test) *
        (data.test$Y - rHat.S.test)
      
      term3 <- -s0.len * rowSums(
        (rHat.s0.test * piHat.s0.test * St.deriv.s0.test) *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      ) * wA_test
      
      term4 <- s0.len * rowSums(
        (rHat.s0.test * St.s0.test) *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      )
      
      psihat.num.term[test] <- term1 + term2 + term3 + term4
      
      ## DENOMINATOR
      den1 <- kernel.S.test * wA_test * St.deriv.S.test
      
      den2 <- -s0.len * rowSums(
        (piHat.s0.test * St.deriv.s0.test) *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      ) * wA_test
      
      den3 <- s0.len * rowSums(
        St.s0.test *
          matrix(kernel.s0, nrow=n_test, ncol=length(s0.seq), byrow=TRUE)
      )
      
      psihat.deno.term[test] <- den1 + den2 + den3
    }
    
    tau_num <- mean(psihat.num.term)
    tau_den <- mean(psihat.deno.term)
    tau     <- tau_num / tau_den
    
    var_hat <- var((psihat.num.term - tau * psihat.deno.term) / tau_den) / n
    
    list(
      tau.h.est = tau,
      tau_num = tau_num,
      tau_den = tau_den,
      CI_lower  = tau - qnorm(0.975) * sqrt(var_hat),
      CI_upper  = tau + qnorm(0.975) * sqrt(var_hat)
    )
  }


############################################################
## Full drop-in code: STWCRVE / CVE EIF-smoothed estimator
## (Fixes the hard-coded pA0/pA1 bug so you can set a0=a1=1)
############################################################

est.cve.eif.smooth.function.new <- function(data, s0, s1,
                                        a0 = 0, a1 = 1,
                                        t = 0.1, h = 0.5,
                                        s0.seq = s0.seq, s1.seq = s1.seq,
                                        epsilon = 1e-2) {
  
  n <- nrow(data)
  nsplits <- 2
  set.seed(123)
  splitIndex <- sample(rep(1:nsplits, length.out = n))
  
  sl.lib <- c("SL.glm", "SL.mean", "SL.ranger")
  
  # AUTO-DETECT: Should X1 be included?
  include_X1 <- length(unique(data$X1)) > 1
  
  if (include_X1) {
    X_vars <- c("A","B","X1","X2")
    X_vars_noA <- c("B","X1","X2")
    
    # replace (S,B) by delta_S when both appear in conditioning set
    X_vars_delta_only     <- c("A","X1","X2")
    X_vars_noA_delta_only <- c("X1","X2")
    
    # density model conditions on (A, S, X) without B
    X_vars_dens <- c("A","S", X_vars_noA_delta_only)
  } else {
    X_vars <- c("A","B","X2")
    X_vars_noA <- c("B","X2")
    
    X_vars_delta_only     <- c("A","X2")
    X_vars_noA_delta_only <- c("X2")
    
    X_vars_dens <- c("A","S", X_vars_noA_delta_only)
  }
  
  # helper to pick correct propensity for any a in {0,1}
  pA_for <- function(a, pA1, pA0) {
    if (a == 1) pA1 else pA0
  }
  
  ############################################################
  ## Helper: compute f(a0,s0), g(a1,s1), ψ_i via cross-fitting
  ############################################################
  get_den_with_eif_cf <- function(data, s0, s1, a0, a1) {
    
    n <- nrow(data)
    s0.len <- s0.seq[2] - s0.seq[1]
    s1.len <- s1.seq[2] - s1.seq[1]
    
    kernel.s0 <- gaussianKernel(s0.seq, s0, h)
    kernel.s0 <- kernel.s0 / sum(kernel.s0 * s0.len)
    
    kernel.s1 <- gaussianKernel(s1.seq, s1, h)
    kernel.s1 <- kernel.s1 / sum(kernel.s1 * s1.len)
    
    f_vec <- numeric(n)
    g_vec <- numeric(n)
    psi_vec <- numeric(n)
    
    for (split in 1:nsplits) {
      
      train <- splitIndex != split
      test  <- splitIndex == split
      
      dt_train <- data[train, ]
      dt_test  <- data[test, ]
      n_test <- sum(test)
      
      # delta_S on BOTH train/test
      dt_train$delta_S <- dt_train$S - dt_train$B
      dt_test$delta_S  <- dt_test$S  - dt_test$B
      
      #########################
      ## Treatment model: A | (X, S, delta_S)
      #########################
      X_ps_train <- dt_train[, c(X_vars_noA_delta_only, "S", "delta_S")]
      X_ps_test  <- dt_test[,  c(X_vars_noA_delta_only, "S", "delta_S")]
      
      psMod <- SuperLearner(
        Y = dt_train$A,
        X = X_ps_train,
        newX = X_ps_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      pA1_test <- as.numeric(psMod$SL.predict)
      pA1_test <- pmax(pmin(pA1_test, 1 - 1e-6), 1e-6)
      pA0_test <- 1 - pA1_test
      
      # ✅ FIX: use correct propensity for each arm
      pAa0_test <- pA_for(a0, pA1_test, pA0_test)
      pAa1_test <- pA_for(a1, pA1_test, pA0_test)
      
      #########################
      ## S model (kept for reference)
      #########################
      X_S_train <- dt_train[, X_vars]
      X_S_test  <- dt_test[,  X_vars]
      
      SMod <- SuperLearner(
        Y = dt_train$S,
        X = X_S_train,
        newX = X_S_test,
        SL.library = sl.lib
      )
      
      # counterfactual conditioning variables
      X_S_test_a0 <- X_S_test; X_S_test_a0$A <- a0
      X_S_test_a1 <- X_S_test; X_S_test_a1$A <- a1
      
      #########################
      ## Outcome model: Y | (A, X, S, delta_S)
      #########################
      X_Y_train <- dt_train[, c(X_vars_delta_only, "S", "delta_S")]
      X_Y_test  <- dt_test[,  c(X_vars_delta_only, "S", "delta_S")]
      
      YMod <- SuperLearner(
        Y = dt_train$Y,
        X = X_Y_train,
        newX = X_Y_test,
        SL.library = sl.lib,
        family = binomial()
      )
      
      # r(a0, S_obs, X) (implemented via observed delta_S)
      X_Y_test_a0 <- X_Y_test; X_Y_test_a0$A <- a0
      r_S_a0 <- as.numeric(predict(YMod, newdata = X_Y_test_a0)$pred)
      
      # r(a0, s0.grid, X): vary delta_S = s0 - B
      dt_s0 <- dt_test[rep(seq_len(n_test), each = length(s0.seq)), ]
      dt_s0$S <- rep(s0.seq, times = n_test)
      dt_s0$A <- a0
      dt_s0$delta_S <- dt_s0$S - dt_s0$B
      X_s0 <- dt_s0[, colnames(X_Y_train)]
      r_s0 <- matrix(
        as.numeric(predict(YMod, newdata = X_s0)$pred),
        nrow = n_test,
        ncol = length(s0.seq),
        byrow = TRUE
      )
      
      # r(a1, s1.grid, X): vary delta_S = s1 - B
      dt_s1 <- dt_test[rep(seq_len(n_test), each = length(s1.seq)), ]
      dt_s1$S <- rep(s1.seq, times = n_test)
      dt_s1$A <- a1
      dt_s1$delta_S <- dt_s1$S - dt_s1$B
      X_s1 <- dt_s1[, colnames(X_Y_train)]
      r_s1 <- matrix(
        as.numeric(predict(YMod, newdata = X_s1)$pred),
        nrow = n_test,
        ncol = length(s1.seq),
        byrow = TRUE
      )
      
      #########################
      ## π(s | X,A) with DELTA METHOD
      #########################
      bw.fit.delta <- npcdensbw(
        ydat = dt_train$delta_S,
        xdat = dt_train[, X_vars_dens],
        bwmethod = "normal-reference"
      )
      
      ## π(S_obs | X, a0)
      X_test_a0 <- dt_test[, X_vars_dens]
      X_test_a0$A <- a0
      delta_obs_a0 <- dt_test$delta_S
      pi_S_a0 <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = X_test_a0,
        eydat = delta_obs_a0
      ))
      pi_S_a0 <- pmax(pi_S_a0, 1e-10)
      
      ## π(S_obs | X, a1)
      X_test_a1 <- dt_test[, X_vars_dens]
      X_test_a1$A <- a1
      delta_obs_a1 <- dt_test$delta_S
      pi_S_a1 <- fitted(npcdens(
        bws = bw.fit.delta,
        exdat = X_test_a1,
        eydat = delta_obs_a1
      ))
      pi_S_a1 <- pmax(pi_S_a1, 1e-10)
      
      St_S_a0 <- getSt(pi_S_a0, epsilon, t)
      St_S_a1 <- getSt(pi_S_a1, epsilon, t)
      
      ## π(s0 | X, a0) on grid via delta
      S_obs_expanded <- rep(dt_test$S, each = length(s0.seq))
      delta_obs_expanded <- rep(dt_test$delta_S, each = length(s0.seq))
      s0_expanded <- rep(s0.seq, times = n_test)
      delta_s0 <- (s0_expanded - S_obs_expanded) + delta_obs_expanded
      valid_s0 <- delta_s0 >= 0
      
      pi_s0_a0_vec <- rep(1e-10, length(delta_s0))
      if (sum(valid_s0) > 0) {
        data_s0_expanded <- dt_test[rep(seq_len(n_test), each = length(s0.seq)), ]
        data_s0_expanded$A <- a0
        data_s0_expanded$S <- s0_expanded
        pi_s0_a0_vec[valid_s0] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s0_expanded[valid_s0, X_vars_dens],
          eydat = delta_s0[valid_s0]
        ))
      }
      pi_s0_a0 <- matrix(pi_s0_a0_vec, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      pi_s0_a0 <- pmax(pi_s0_a0, 1e-10)
      
      ## π(s1 | X, a1) on grid via delta
      S_obs_expanded <- rep(dt_test$S, each = length(s1.seq))
      delta_obs_expanded <- rep(dt_test$delta_S, each = length(s1.seq))
      s1_expanded <- rep(s1.seq, times = n_test)
      delta_s1 <- (s1_expanded - S_obs_expanded) + delta_obs_expanded
      valid_s1 <- delta_s1 >= 0
      
      pi_s1_a1_vec <- rep(1e-10, length(delta_s1))
      if (sum(valid_s1) > 0) {
        data_s1_expanded <- dt_test[rep(seq_len(n_test), each = length(s1.seq)), ]
        data_s1_expanded$A <- a1
        data_s1_expanded$S <- s1_expanded
        pi_s1_a1_vec[valid_s1] <- fitted(npcdens(
          bws = bw.fit.delta,
          exdat = data_s1_expanded[valid_s1, X_vars_dens],
          eydat = delta_s1[valid_s1]
        ))
      }
      pi_s1_a1 <- matrix(pi_s1_a1_vec, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      pi_s1_a1 <- pmax(pi_s1_a1, 1e-10)
      
      St_s0_a0     <- getSt(pi_s0_a0, epsilon, t)
      St_s1_a1     <- getSt(pi_s1_a1, epsilon, t)
      St_der_s0_a0 <- getSt.deriv(pi_s0_a0, epsilon, t)
      St_der_s1_a1 <- getSt.deriv(pi_s1_a1, epsilon, t)
      
      #########################
      ## Kernels at observed S
      #########################
      kern_s0_obs <- gaussianKernel(dt_test$S, s0, h)
      kern_s0_obs <- kern_s0_obs / sum(kernel.s0 * (s0.seq[2] - s0.seq[1]))
      
      kern_s1_obs <- gaussianKernel(dt_test$S, s1, h)
      kern_s1_obs <- kern_s1_obs / sum(kernel.s1 * (s1.seq[2] - s1.seq[1]))
      
      #########################
      ## Integrals f(a0,s0), g(a1,s1)
      #########################
      f_test <- s0.len * rowSums(
        (St_s0_a0 * r_s0) * matrix(kernel.s0, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      )
      
      g_test <- s1.len * rowSums(
        (St_s1_a1 * r_s1) * matrix(kernel.s1, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      )
      
      #########################
      ## IPW residuals (for a0)
      #########################
      ipwRes <- (dt_test$Y - r_S_a0) / pi_S_a0
      
      #########################
      ## EIF for f(a0,s0)  ✅ uses pAa0_test
      #########################
      Df1 <- kern_s0_obs * ((dt_test$A == a0) / pAa0_test) * St_S_a0 * r_S_a0
      
      Df2 <- kern_s0_obs * ((dt_test$A == a0) / pAa0_test) * St_S_a0 * ipwRes
      
      Df3 <- -s0.len * rowSums(
        (r_s0 * pi_s0_a0 * St_der_s0_a0) *
          matrix(kernel.s0, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      ) * ((dt_test$A == a0) / pAa0_test)
      
      Df4 <- s0.len * rowSums(
        (r_s0 * St_s0_a0) *
          matrix(kernel.s0, nrow = n_test, ncol = length(s0.seq), byrow = TRUE)
      )
      
      Df <- Df1 + Df2 + Df3 + Df4
      
      #########################
      ## EIF for g(a1,s1) ✅ uses pAa1_test
      #########################
      Dg1 <- kern_s1_obs * ((dt_test$A == a1) / pAa1_test) * St_S_a1
      
      Dg2 <- -s1.len * rowSums(
        (pi_s1_a1 * St_der_s1_a1) *
          matrix(kernel.s1, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      ) * ((dt_test$A == a1) / pAa1_test)
      
      Dg3 <- s1.len * rowSums(
        (St_s1_a1) *
          matrix(kernel.s1, nrow = n_test, ncol = length(s1.seq), byrow = TRUE)
      )
      
      Dg <- Dg1 + Dg2 + Dg3
      
      #########################
      ## EIF ψ_i
      #########################
      psi <- f_test * Dg + g_test * Df - f_test * g_test
      
      #########################
      ## store
      #########################
      f_vec[test] <- f_test
      g_vec[test] <- g_test
      psi_vec[test] <- psi
    }
    
    list(
      f_vec = f_vec,
      g_vec = g_vec,
      psi_vec = psi_vec,
      est = mean(psi_vec),
      true = mean(f_vec * g_vec)
    )
  }
  
  ###############################################################
  ## Numerator (swap s0/s1, a0/a1)
  ###############################################################
  den <- get_den_with_eif_cf(data, s0, s1, a0, a1)
  num <- get_den_with_eif_cf(data, s1, s0, a1, a0)
  
  ratio_hat <- num$est / den$est
  cve_hat <- 1 - ratio_hat

  psi_cve <- -(num$psi_vec - ratio_hat * den$psi_vec) / den$est
  var_cve <- var(psi_cve) / n
  se_cve  <- sqrt(var_cve)

  psi_log_ratio <- num$psi_vec / num$est - den$psi_vec / den$est
  if (is.finite(ratio_hat) && ratio_hat > 0) {
    se_log_ratio <- sqrt(var(psi_log_ratio) / n)
    ci_ratio_log <- ratio_hat * exp(c(-1.96, 1.96) * se_log_ratio)
    ci_ratio_log_to_cve <- c(1 - ci_ratio_log[2], 1 - ci_ratio_log[1])
  } else {
    ci_ratio_log_to_cve <- c(NA_real_, NA_real_)
  }

  list(
    cve_est             = cve_hat,
    cve_num             = num$est,
    cve_den             = den$est,
    se                  = se_cve,
    ci                  = c(cve_hat - 1.96 * se_cve,
                             cve_hat + 1.96 * se_cve),
    ci_ratio_log_to_cve = ci_ratio_log_to_cve,
    cve_true            = 1 - num$true / den$true
  )
}

############################################################
## Example usage:
## 1) Original VE-like contrast:
## est.cve.eif.smooth.function(df, s0=..., s1=..., a0=0, a1=1)
##
## 2) Your requested "within A=1" version:
## est.cve.eif.smooth.function(df, s0=..., s1=..., a0=1, a1=1)
############################################################
