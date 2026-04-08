#####################################################
# First Part (Superpopulation Estimation)
#   Compute the true dose-response curve and corresponding trimmed estimands using a superpopulation of size N=10^5
#   The results here serve as ground truth benchmarks

# tau.function: ground truth of plug-in of CoR
# tau.smoothed.function: ground truth of EIF of CoR
# cve.function: ground truth of plug-in of CVE
# cve.den.smooth.function: ground truth of EIF of CVE

gaussianKernel = function(S, s, h){
  u = S - s
  kernel = (1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
  return(kernel)
}

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

tau.smoothed.function <- function(n = 1e4, s, t = 0.1, h = 0.2, s0.seq = s0.seq,
                                  epsilon = 10^(-2), additional_info = FALSE,
                                  data = NULL){
  if (is.null(data)) {
    data <- generate_data(n)
  } else {
    n <- nrow(data)
  }
  
  s0.len = s0.seq[2] - s0.seq[1]
  # smoothedNum.itemized = matrix(nrow = length(data[, 1]), ncol = length(s0.seq))
  # smoothedDenom.itemized = matrix(nrow = length(data[, 1]), ncol = length(s0.seq))
  
  #for(s0 in 1:length(s0.seq)){
    kernel.s0 = gaussianKernel(S = s0.seq, s = s, h = h)
    kernel.s0 = kernel.s0/sum(kernel.s0*s0.len)
    
    
    # \hat_r(A, s0, B, X) for each individual
    # md_r = glm(Y ~ S + X2 + X3 + B + A, family = "binomial", data = data)
    # data_new = data[rep(seq_len(nrow(data)), each = length(s0.seq)), ]
    # data_new$S = rep(s0.seq, times = nrow(data))
    # rHat.s0 = matrix(as.numeric(predict(md_r, newdata = data_new, type = "response")), ncol = length(s0.seq), byrow = T)
    
    data_new = data[rep(seq_len(nrow(data)), each = length(s0.seq)), ]
    data_new$s0 = rep(s0.seq, times = nrow(data))
    #rHat.s0 = matrix(expit(0.5 * data_new$X2 + 2 * data_new$X3 - 0.5 * data_new$s0  + data_new$A + 0.3 * data$B), 
    #                 ncol = length(s0.seq), byrow = T)
    rHat.s0 = matrix(expit(0.5 * data_new$X2 + 2 * data_new$X3 - 0.2 * data_new$s0  - data_new$A - 0.3 * data_new$B + 1.5), 
                     ncol = length(s0.seq), byrow = T)

    # S(\pi(s_0), t)
    #piHat.s0 = sapply(s0.seq, FUN = dnorm, mean = data$B + data$A - 0.5 * data$X1 + data$X2^2 + 4, sd = 1)
    piHat.s0 = sapply(s0.seq, FUN = dnorm, mean = data$B + data$A - 0.5 * data$X1 + data$X2^2 + 4, sd = 1)

    St.s0 = getSt(piHat.s0, epsilon = epsilon, t = t)
    smoothedNum.itemized = s0.len * rowSums(t(apply(St.s0*1*rHat.s0,
                                                   MARGIN = 1, FUN = function(x) x*kernel.s0))) 
    smoothedDenom.itemized = s0.len * rowSums(t(apply(St.s0*1,
                                                    MARGIN = 1, FUN = function(x) x*kernel.s0)))
    
  #}
  tau.smoothed = mean(smoothedNum.itemized)/mean(smoothedDenom.itemized)
  
  if(additional_info){
    return(list(tau.smoothed.num = mean(smoothedNum.itemized),
                tau.smoothed.den = mean(smoothedDenom.itemized),
                tau.smoothed = tau.smoothed
                ))
  }else{
    return(list(tau.h = tau.smoothed))
  }
}
