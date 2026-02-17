# Simulation wrapper for two-sample LSM networks with mesoscale 
# (rectangular) differences

# Mesoscale model specifies at most
# - name (LSM or wLSM)
# - hyp_set (must be a rectangle, specified as list(row_indices,col_indices))
# - directed (default TRUE)
# - self_loops (default TRUE)
# - dispersion (default 1)
# - d 
# - signal (default 0), perturbation variance of pos'ns incident to hyp_set = signal^2 / d
# - similarity ('dist' or 'ip')
# - alpha (intercept if similarity='dist', defaults to 0)

# NOTE: only supports logit link for binary n/ws, identity link/Gaussian edges
# for weighted networks

# main function
Simulate_mesoscale <- function(n,model=list(),m){
  # check hypothesis set
  if(!(length(model$hyp_set)==2)){
    stop('Only supports rectangular hypothesis sets')
  }
  # set default parameters
  # default undirected
  if(is.null(model$directed)){
    model$directed <- TRUE
  }
  # default no self loops
  if(is.null(model$self_loops)){
    model$self_loops <- TRUE
  }
  # default dispersion/sigma^2 = 1
  if(is.null(model$dispersion)){
    model$dispersion <- 1
  }
  # default signal = 0
  if(is.null(model$signal)){
    model$signal <- 0
  }
  # default alpha (dist model intercept) = 0
  if(is.null(model$alpha)){
    model$alpha <- 0
  }
  # default link function (logit for LSM)
  if(model$name=='LSM'){
    model$link <- 'logit'
  }
  # copy models for two-samples
  model1 <- model2 <- model
  # E(A) for each sample
  # recover row/col indices, perturbation std.dev (psd)
  row_indices <- model$hyp_set[[1]]; col_indices <- model$hyp_set[[2]]
  all_indices <- union(row_indices,col_indices)
  psd <- model$signal/sqrt(model$d)
  # generate positions
  if(model$directed){
    # LHS positions
    Z1 <- matrix(stats::rnorm(n*d),n,d); Z2 <- matrix(stats::rnorm(n*d),n,d)
    Z2[row_indices,] <- perturb_mat(Z1[row_indices,],sd=psd)
    # RHS positions
    Y1 <- matrix(stats::rnorm(n*d),n,d); Y2 <- matrix(stats::rnorm(n*d),n,d)
    Y2[col_indices,] <- perturb_mat(Y1[col_indices,],sd=psd)
    # store positions
    model1$Z <- Z1; model2$Z <- Z2
    model1$Y <- Y1; model2$Y <- Y2
  }
  else{
    # LHS positions
    Z1 <- matrix(stats::rnorm(n*d),n,d); Z2 <- matrix(stats::rnorm(n*d),n,d)
    Z2[all_indices,] <- perturb_mat(Z1[all_indices,],sd=psd)
    # store positions
    model1$Z <- Z1; model2$Z <- Z2
  }
  # generate networks
  A1 <- Simulate_netmodel(n,model1,m)$A; A2 <- Simulate_netmodel(n,model2,m)$A
  # store and return
  if(model$directed){
  return(list(A1=A1,A2=A2,
              Z1=Z1,Z2=Z2,
              Y1=Y1,Y2=Y2))
  }
  else{
    return(list(A1=A1,A2=A2,
                Z1=Z1,Z2=Z2))
  }
}

# #### test instances ####
# set.seed(12)
# # all with n=100,m=10,d=3,hyp_set=list(1:20,71:100)
# n <- 100; m <- 10
# H <- list(1:20,71:100)
# 
# # (1) null Gaussian IP
# model1 <- list(name='wLSM',hyp_set=H,d=3,signal=0,dispersion=1,similarity='ip')
# dat1 <- Simulate_mesoscale(n,model1,m)
# image(dat1$A1[[1]] - dat1$A2[[1]])
# 
# # (2) non-null Gaussian IP
# model2 <- list(name='wLSM',hyp_set=H,d=3,signal=1/sqrt(2),dispersion=1,similarity='ip')
# dat2 <- Simulate_mesoscale(n,model2,m)
# image(dat2$A1[[1]] - dat2$A2[[1]])
# 
# # (3) null Gaussian distance, undirected
# model3 <- list(name='wLSM',hyp_set=H,d=3,signal=0,dispersion=1,similarity='dist',directed=FALSE)
# dat3 <- Simulate_mesoscale(n,model3,m)
# image(dat3$Z1 - dat3$Z2)
# 
# # (4) null logistic IP
# model4 <- list(name='LSM',hyp_set=H,d=3,signal=0,similarity='ip')
# dat4 <- Simulate_mesoscale(n,model4,m)
# table(dat4$A1[[1]] - dat4$A2[[1]])
# 
# # (5) non-null logistic IP
# model5 <- list(name='LSM',hyp_set=H,d=3,signal=0.5/sqrt(2),similarity='ip')
# dat5 <- Simulate_mesoscale(n,model5,m)
# table(dat5$A1[[1]] - dat5$A2[[1]])
# 
# # (6) non-null logistic IP, overdispersed
# model6 <- list(name='LSM',hyp_set=H,d=3,signal=0.5/sqrt(2),dispersion=2,similarity='ip')
# dat6 <- Simulate_mesoscale(n,model6,m)
# table(dat6$A1[[1]] - dat6$A2[[1]])