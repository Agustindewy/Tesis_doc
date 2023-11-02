
# Author Federico Cortés (cortes.federico@gmail.com)

#-------------------------------------- Function perry_D ------------------------------------------

# perry_D calculates de distance between f(t) and g(t)

perry_D <- function(AMB, CAP, min_t, max_t, rango_t) {
  
  # Shots frequency per environmental range, or f(t) 
  xn <- seq(from = min_t, to = max_t, by = rango_t)
  frec_lances = numeric(length(xn))
  
  # the sum per range of X for the first interval 
  for (i in 1:length(xn) - 1) {
    idx <- which(AMB >= xn[i] & AMB < xn[i + 1])
    frec_lances[i] <- length(AMB[idx]) / length(AMB)
  }
  
  # and for the last interval
  idx <- which(AMB >= xn[i + 1])
  frec_lances[i + 1] <- length((AMB[idx])) / length(AMB)
  
  Ft <- cumsum(frec_lances) # f(t)

  # Frequency of catch per environmental range, or g(t)
  frec_cap <- numeric(length(xn))
  
  # the sum per range of X for the first interval 
  for (i in 1:length(xn) - 1) {
    idx <- which(AMB >= xn[i] & AMB < xn[i + 1])
    frec_cap[i] <- sum((CAP)[idx]) / sum(CAP)
  }
  
  # and for the last interval
  idx <- which(AMB >= xn[i + 1])
  frec_cap[i + 1] <- sum((CAP)[idx]) / sum(CAP)
  
  Gt<-cumsum(frec_cap) # g(t)

  # Maximum distance between f(t) and g(t)
  Dmax <- max(abs(Gt - Ft))
  
  # Results
  res_perryD <- list(Xambiente = xn, Ft = Ft, Gt = Gt, Dmax = Dmax)
  res_perryD
}


#------------------------------------- Function perry_BOOT ----------------------------------------

# perry_BOOT calculates the p-value for the analysis

perry_BOOT <- function(AMB, CAP, min_t, max_t, rango_t, N_iter) {
  
  # Dmax distribution based on bootstrap and randomization of environment and catch
  DmaxDIST <- numeric(length = N_iter)
  FtRAN <- matrix(NA, ncol = N_iter, nrow = length(seq(from = min_t, to = max_t, by = rango_t)))
  GtRAN <- matrix(NA, ncol = N_iter, nrow = length(seq(from = min_t, to = max_t, by = rango_t)))
  
  for (i in 1:length(DmaxDIST)) {
    
    # bootstrap with and without replacement
    bootAMB <- sample(AMB, replace = T)
    randomCAP <- sample(CAP, replace = F)
    
    # Dmax calculation for bootstrapped data
    perryBOOT <- perry_D(bootAMB, randomCAP, min_t = min_t, max_t = max_t, rango_t = rango_t)
    
    # save iteration results
    DmaxDIST[i] <- perryBOOT$Dmax
    FtRAN[, i] <- perryBOOT$Ft
    GtRAN[, i] <- perryBOOT$Gt
  }
  
  # Observed Dmax (original data, not bootstrapped)
  DmaxOBS <- perry_D(AMB, CAP, min_t = min_t, max_t = max_t, rango_t = rango_t)$Dmax
  
  # if p-value is < than 0.01 we assume association between environment and catch 
  p.value <- length(which(DmaxDIST >= DmaxOBS)) / N_iter
  
  # Results
  res_perryBOOT <- list(Dmax = DmaxOBS, p = p.value, distD = DmaxDIST,
                        Xambiente = perryBOOT$Xambiente, randomFt = FtRAN, randomGt = GtRAN)
  res_perryBOOT
}


#--------------------------------------------- END -------------------------------------
