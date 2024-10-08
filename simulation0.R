library(dplyr)
library(tidyr)
library(nnet)

tau_fun <- function(data, a, s, epsilon = 0.1){
  # model for estimating r(a, s, B, X) := r(A = a, S = s | B, X) = P(Y | B, X)
  # md_r <- glm(Y ~ X1 + X2 + X3 + B, family = "binomial", data = data, subset = (A == a) & (S == s))
  md_r <- glm(Y ~ X1 + X2 + X3 + B + S, family = "binomial", data = data, subset = (A == a))
  
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
    numerator = sum(data$S == s & A == a & data$B == b & data$X1 == x1)
    return(numerator/denominator)
  }
  
  est_p_s <- sapply(1:nrow(data), function(i) {
    row <- data[i, ]
    p_s(data, row$B, row$X1)
  })
  
  est_p_a = 0.5
  
  weight_ps = est_p_s * est_p_a
  weight_indicator = (weight_ps > epsilon)
  
  tau_weighted_ps <- mean(weight_ps * est_r)/mean(weight_ps)
  tau_weighted_indicator <- mean(weight_indicator * est_r)/mean(weight_indicator)
  tau <- mean(est_r)
  
  return(list(tau_weighted_ps = tau_weighted_ps, 
              tau_weighted_indicator = tau_weighted_indicator, 
              tau = tau))
}
# Simulate some data
set.seed(123) # For reproducibility


N <- 1000 # Number of observations
A <- rbinom(N, 1, 0.5) # Simulate treatment assignment
X1 <- rbinom(N, 1, 0.3) # 1 = non-naive; 0 = naive
X2 <- runif(N)
X3 <- runif(N)
# X <- cbind(X1, X2, X3) # Simulate some covariates
p0 <- c(0.2, 0.3, 0.4, 0.05, 0.05)
p1 <- c(0.1, 0.15, 0.3, 0.3, 0.15)
B <- rep(0, N)
B[X1 == 0] <- sample(1:5, sum(1 - X1), replace = TRUE, prob = p0) # Simulate baseline marker levels
B[X1 == 1] <- sample(1:5, sum(X1), replace = TRUE, prob = p1) # Simulate baseline marker levels
hist(B)
hist(B[X1 == 0])
hist(B[X1 == 1])
S <- (B + 2 * A + rbinom(N, 2, 0.6)) # Simulate post-vaccination marker levels
expit <- function(x) exp(x)/(1+exp(x))
p_Y <- expit(0.5 * X2 - X3 - 0.5 * S * X1 + A - 0.3 * B)
Y <- rbinom(N, 1, p_Y)
data <- data.frame(A, B, S, Y, X1, X2, X3)
summary(data)


a <- 1 # Change as needed
s = 3
table(data$S[data$A == a])

S_vector <- 3:9
rslt <- as.data.frame(sapply(S_vector,  function(s) tau_fun(data, a, s, 0.05)))
tau_weighted_ps <- unlist(rslt[1, ])
tau_weighted_indicator <- unlist(rslt[2, ])
tau <- unlist(rslt[3, ])


# add dgp curve(S - tau)
library(ggplot2)

dt1 <- data.frame(S = S_vector,
                  tau_type = c(rep("tau", length(S_vector)), rep("tau_weighted_ps", length(S_vector)), rep("tau_weighted_indicator", length(S_vector))),
                  tau_value = c(tau, tau_weighted_ps, tau_weighted_indicator))
curve <- ggplot(dt1, aes(x=S, y=tau_value, group=tau_type)) +
  geom_line(aes(color=tau_type))+
  geom_point(aes(color=tau_type))+
  scale_color_brewer(palette="Paired")+
  theme_minimal()

curve

# true treatment effect - prob = 4, 5/ = 4, 5
# develop a better estimator - peter's paper (mediation analysis) for literature review (intro) question is not well post, S= 3 subgp, s prob = 3 subgp
# read oliver's paper
# EIF estimators -  functional analysis
