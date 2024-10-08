library(dplyr)
library(tidyr)
library(nnet)
library(mgcv)
library(randomForest)

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
  }



tau_fun <- function(data, a, s, t = 0.1){
  # model for estimating r(a, s, B, X) := r(A = a, S = s | B, X) = P(Y | B, X)
  # md_r <- glm(Y ~ X1 + X2 + X3 + B, family = "binomial", data = data, subset = (A == a) & (S == s))
  md_r <- glm(Y ~ S + X2 + X3 + B + S, family = "binomial", data = data, subset = (A == a))
  # I(T0 <= t) = Y binary
  # Calculate estimated event rates r(a, s, B, X) for each individual
  # est_r <- predict(md_r, newdata = data, type = "response")
  
  data_new <- data %>% 
    mutate(S = s)
  est_r <- predict(md_r, newdata = data_new, type = "response")
  
  # models for estimating P(S = s | A = a, B, X) and P(A = a | B, X),
  # md_s <- multinom(S ~ A + B + X1 + X2 + X3, data = data, subset = (A == a), trace = F)
  # md_a <- glm(A ~ B + X1 + X2 + X3, data = data, family = "binomial")

  # est_p_s <- predict(md_s, newdata = data, type = "probs")
  # est_p_a <- predict(md_a, newdata = data, type = "response")
  
  # summary(data[data$S == s & data$A == a, ])
  # table(data$X1[data$S == s & data$A == a], data$B[data$S == s & data$A == a])
  p_s <- function(data, b, x1){
    denominator = sum(data$A == a & data$B == b & data$X1 == x1)
    numerator = sum(data$S == s & data$A == a & data$B == b & data$X1 == x1)
    return(numerator/denominator)
  }
  
  est_p_s <- sapply(1:nrow(data), function(i) {
    row <- data[i, ]
    p_s(data, row$B, row$X1)
  })
  
  md_A <- glm(A ~ X1, data = data, family = 'binomial')
  est_p_a <- predict(md_A, newdata = data, type = "response")
  
  
  weight_ps = est_p_s * est_p_a
  weight_indicator = (weight_ps > t)
  
  tau_weighted_ps <- mean(weight_ps * est_r)/mean(weight_ps)
  tau_weighted_ps_eif <- mean((data$S == s & data$A == a)*data$Y)/mean(data$S == s & data$A == a)
  tau_weighted_indicator <- mean(weight_indicator * est_r)/mean(weight_indicator)
  tau <- mean(est_r)
  
  return(list(tau_weighted_ps = tau_weighted_ps, 
              tau_weighted_ps_eif = tau_weighted_ps_eif,
              tau_weighted_indicator = tau_weighted_indicator, 
              tau = tau))
}


  #---------------------------------#
  #   tau_weighted_indicator_eif.   #
  #---------------------------------#

tau_weighted_indicator_eif <- 
  function(data, a, t = 0.1, epsilon = 10^(-2), s = 3, s.seq = seq(2, 7, by = 0.05), h = 0.3){
  n = nrow(data)
  s.seq.len <- length(s.seq)
  
  md_A <- glm(A ~ X1, data = data, family = 'binomial')
  est_p_a <- predict(md_A, newdata = data, type = "response")
  
  md_pi_S = lm(S ~ A + B + X1 + I(X2^2), data = data)
  meanS.est = md_pi_S$fitted.values
  
  #\hat{\pi}(S | A, B, X)
  piHat.S = dnorm(data$S, mean = meanS.est, sd = 1)
  

  piHat.s0 = sapply(s.seq, FUN = dnorm, mean = meanS.est, sd = 1)
  
  #To implement the trimming estimator,
  #it'll be helpful to have functions for S_t and S'_t
  getSt = function(piHat, epsilon = 10^(-4), t = 0.1){
    St = pnorm(piHat - t, mean = 0, sd = epsilon)
    return(St)
  }
  getSt.deriv = function(piHat, epsilon = 10^(-4), t = 0.1){
    St.deriv = dnorm(piHat - t, mean = 0, sd = epsilon)
    return(St.deriv)
  }
  
  #compute smoothed indicator and its derivative.
  St.S = getSt(piHat = piHat.S, epsilon = epsilon, t = t)
  St.deriv.S = getSt.deriv(piHat = piHat.S, epsilon = epsilon, t = t)
  St.s0 = getSt(piHat.s0, epsilon = epsilon, t = t)
  St.deriv.s0 = getSt.deriv(piHat.s0, epsilon = epsilon, t = t)
  
  
  md_r = glm(Y ~ S + X2 + X3 + B + A, family = "binomial", data = data)
  # hat_r(A, S, B, X) for each individual
  rHat.S = predict(md_r, newdata = data, type = "response")
  
  # hat_r(A, s0, B, X) for each individual
  # 
  data_new = data[rep(seq_len(nrow(data)), each = length(s.seq)), ]
  data_new$S = rep(s.seq, times = nrow(data))
  rHat.s0 = matrix(as.numeric(predict(md_r, newdata = data_new, type = "response")), ncol = length(s.seq), byrow = T)
  
  gaussianKernel = function(S, s, h){
    #the difference is
    u = S - s
    #then, the kernel is
    kernel = (1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
    return(kernel)
  }
  
  
  kernel.S = gaussianKernel(S = data$S, s = s, h = h)
  kernel.s0 = gaussianKernel(S = s.seq, s = s, h = h)
  
  #normalize the kernels
  kernel.S = kernel.S/sum(kernel.s0*s.seq.len)
  kernel.s0 = kernel.s0/sum(kernel.s0*s.seq.len)
  
  ipwRes = (data$Y = rHat.S)/piHat.S
  
  psihat.num.term1 = kernel.S*St.deriv.S*rHat.S
  psihat.num.term2 = kernel.S*ipwRes*St.S
  psihat.num.term3 = -rowSums(t(apply(rHat.s0*piHat.s0*St.deriv.s0*s.seq.len,
                                      MARGIN = 1, FUN = function(x) x*kernel.s0)))
  psihat.num.term4 = rowSums(t(apply(rHat.s0*St.s0*s.seq.len,
                                        MARGIN = 1, FUN = function(x) x*kernel.s0))) 
  
  psihat.num.term = psihat.num.term1 + psihat.num.term2 + psihat.num.term3 + psihat.num.term4
  
  psihat.deno.term1 = kernel.S*St.deriv.S
  psihat.deno.term2 = -rowSums(t(apply(piHat.s0*St.deriv.s0*s.seq.len,
                                      MARGIN = 1, FUN = function(x) x*kernel.s0)))
  psihat.deno.term3 = rowSums(t(apply(St.s0*s.seq.len,
                                     MARGIN = 1, FUN = function(x) x*kernel.s0))) 
  
  psihat.deno.term = psihat.deno.term1 + psihat.deno.term2 + psihat.deno.term3
  
  tau_weighted_indicator_eif = (mean(psihat.num.term)/mean(psihat.deno.term))
  return(tau_weighted_indicator_eif)
  
}
  
  
  

# Simulate some data


generate_data <- function(N){
  
  X1 <- rbinom(N, 1, 0.3) # 1 = non-naive; 0 = naive
  X2 <- runif(N)
  X3 <- runif(N)
  # X <- cbind(X1, X2, X3) # Simulate some covariates
  p_A <- expit(2 * X1 - 1)
  A <- rbinom(N, 1, p_A) # Simulate treatment assignment (depend on binary x1)
  
  p0 <- c(0.2, 0.3, 0.4, 0.05, 0.05)
  p1 <- c(0.1, 0.15, 0.3, 0.3, 0.15)
  B <- rep(0, N)
  B[X1 == 0] <- sample(1:5, sum(1 - X1), replace = TRUE, prob = p0) # Simulate baseline marker levels
  B[X1 == 1] <- sample(1:5, sum(X1), replace = TRUE, prob = p1) # Simulate baseline marker levels
  # hist(B)
  # hist(B[X1 == 0])
  # hist(B[X1 == 1])
  #S <- (B + 2 * A + rnorm(N, 2, 0.2)) # Simulate post-vaccination marker levels
  
  S <- B + 2 * A - 0.5 * X1 + X2^2 + 2 + rnorm(N, 0, 1)
  
  p_Y <- expit(0.5 * X2 - X3 - 0.5 * S + A - 0.3 * B)
  Y <- rbinom(N, 1, p_Y)
  data <- data.frame(A, B, S, Y, X1, X2, X3)
  return(data)
}

# s = 3:9, a = 1
est_ground_truth <- function(s, a = 1, t = 0.1, n = 1e4){
  data1 <- generate_data(n)
  est_p_a = expit(2 * data$X1 - 1)
  p_s <- function(data, b, x1){
    denominator = sum(data$A == a & data$B == b & data$X1 == x1)
    numerator = sum(data$S == s & data$A == a & data$B == b & data$X1 == x1)
    return(numerator/denominator)
  }
  
  est_p_s <- sapply(1:nrow(data1), function(i) {
    row <- data1[i, ]
    p_s(data1, row$B, row$X1)
  })
  
  est_r <- expit(0.5 * data1$X2 - data1$X3 - 0.5 * data1$S  - data1$A - 0.3 * data1$B)
  
  weight_ps = est_p_s * est_p_a
  weight_indicator = (weight_ps > t)
  
  tau_weighted_ps <- mean(weight_ps * est_r)/mean(weight_ps)
  tau_weighted_indicator <- mean(weight_indicator * est_r)/mean(weight_indicator)
  return(list(tau_weighted_ps_true = tau_weighted_ps, 
              tau_weighted_indicator_true = tau_weighted_indicator))
  
}

N <- 1000
set.seed(123) 
data <- generate_data(N)



a <- 1 # Change as needed
s = 3

S_vector <- 1:9
rslt <- as.data.frame(sapply(S_vector,  function(s) tau_fun(data, a, s, 0.05)))
tau_weighted_ps <- unlist(rslt[1, ])
tau_weighted_ps_eif <- unlist(rslt[2, ])

tau_weighted_indicator <- unlist(rslt[3, ])
tau <- unlist(rslt[4, ])

rslt_ground_truth <- as.data.frame(sapply(S_vector,  function(s) est_ground_truth(s, a = 1, t = 0.05, n = 5e4)))
tau_weighted_ps_true <- unlist(rslt_ground_truth[1, ])
tau_weighted_indicator_true <- unlist(rslt_ground_truth[2, ])
# add dgp curve(S - tau)
library(ggplot2)

# dt1 <- data.frame(S = S_vector,
#                   tau_type = c(rep("tau", length(S_vector)), 
#                                rep("tau_weighted_ps", length(S_vector)), 
#                                rep("tau_weighted_indicator", length(S_vector)),
#                                rep("tau_weighted_ps_true", length(S_vector)),
#                                rep("tau_weighted_indicator_true", length(S_vector))),
#                   tau_value = c(tau, tau_weighted_ps, tau_weighted_indicator, 
#                                 tau_weighted_ps_true, tau_weighted_indicator_true))

dt1 <- data.frame(S = S_vector,
                  tau_type = c(rep("tau", length(S_vector)), 
                               rep("tau_weighted_ps", length(S_vector)), 
                               rep("tau_weighted_ps_eif", length(S_vector)), 
                               rep("tau_weighted_indicator", length(S_vector))),
                  tau_value = c(tau, tau_weighted_ps, tau_weighted_ps_eif, tau_weighted_indicator))


library(latex2exp)
dt1$tau_type <- factor(dt1$tau_type, labels = c("tau" = parse(text=TeX('($\\tau$)')),
                                                "tau_weighted_ps" = parse(text=TeX('($\\tau_{\\omega_{1, s, 1}}(1, s)$)')), 
                                                "tau_weighted_ps" = parse(text=TeX('($\\tau_{\\omega_{1, s, 1}}(1, s), eif$)')), 
                                                "tau_weighted_indicator" = parse(text=TeX('($\\tau_{\\omega_{1, s, 2}}(1, s)$)'))))

dt1$estimand <- c(rep("3", length(S_vector)), 
                  rep("4", length(S_vector)), 
                  rep("1", length(S_vector)), 
                  rep("2", length(S_vector)))



curve <- ggplot(dt1, aes(x=S, y=tau_value, group=estimand)) +
  #geom_line(aes(color=tau_type))+
  geom_point(aes(color=estimand), size = 3)+
  scale_color_brewer(palette="Paired")+
  theme_minimal()+
  xlab("The post-vaccination immune marker level S = s")+
  ylab("Controlled risk ") 

curve

# true treatment effect - prob = 4, 5/ = 4, 5
# develop a better estimator - peter's paper (mediation analysis) for literature review (intro) question is not well post, S= 3 subgp, s prob = 3 subgp
# read oliver's paper
# EIF estimators -  functional analysis