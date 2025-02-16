####################################################################################
####################################################################################
# FUNCTION TO COMPUTE THE GNEDIN URN SCHEME ########################################
####################################################################################
####################################################################################

urn_GN <- function(v_minus,gamma_GN){
  H<-length(v_minus)
  return(c((v_minus+1)*(sum(v_minus)-H+gamma_GN),H^2-H*gamma_GN))
}

####################################################################################
####################################################################################
# GIBBS SAMPLER FOR THE BINARY STOCHASTIC BLOCK MODEL  #############################
####################################################################################
####################################################################################

# Inputs:
# Y_bin = VxV symmetric adjacency matrix
# seed = a seed for reproducibility
# z_init = V-vector of initialization assignment for each node (default = one cluster for each node)
# a,b = parameters of the Beta prior on the pi
# x = V-vector of categorical covariates
# if x is specified, also alpha_xi (a C-vector of parameters for a Dirichlet distribution) must be set

# Output:
# Posterior samples of the community labels for each node v=1,...,V

esbm_binary <- function(Y_bin, seed, N_iter, z_init=c(1:nrow(Y_bin)), a=1, b=1, gamma_GN=NA, x=NULL, alpha_xi=NULL){
  
  # ----------------------------------------------
  # Select the Gnedin urn scheme
  # ----------------------------------------------

    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  
  # ----------------------------------------------
  # Pre-processing of the node attributes
  # ----------------------------------------------
  
  if (!is.null(x)){
    print("Covariates have been provided")
    x <- as.numeric(as.factor(x))
    X <- vec2mat(x)
    if (!is.null(alpha_xi)){
      alpha0 <- sum(alpha_xi)
    } else {
      stop("If covariates x are given, then alpha_xi must be set as well")
    }
  }
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  set.seed(seed)
  
  V <- nrow(Y_bin)
  
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- vec2mat(z_init)
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  
  # create the matrix with block connections
  # connections of each node with those in the different blocks
  temp   <- Y_bin%*%Z
  # connections in each block
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
    for (v in 1:V){
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      # Z_v = cluster indicator matrix without node v
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y_bin[-v,v])      
      # current cluster of node v
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m matrix by difference
      if(length(h_v) == 1){
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
        m            <- m_full - resid1
      } else {m <- m_full} # No need to update m in this case
      
      # number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      # required to compute increase in the number of non-edges in any cluster h when adding v to this cluster
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      # increase in the number of edges in any cluster h when adding v to this cluster
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      # log-likelihoods of old and new clusters
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      log_addit    <- 0
      
      # supervision by attributes
      if(!is.null(x)){
        Vx        <- crossprod(Z_v, X[-v,])
        addit_old <- (Vx[,x[v]] + alpha_xi[x[v]]) / (v_minus+alpha0)
        addit_new <- alpha_xi[x[v]] / alpha0
        log_addit <- log(c(addit_old, addit_new))
      }
      
      # full-conditional probabilities of cluster assignments
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v             <- crossprod(Z[-v,], Y_bin[-v,v])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full
      resid2              <- matrix(0,H,H)
      resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
      m_full              <- m + resid2
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}

####################################################################################
####################################################################################
# GIBBS SAMPLER FOR THE POISSON STOCHASTIC BLOCK MODEL  ############################
####################################################################################
####################################################################################

# Inputs:
# Y_count = VxV symmetric (weighted) adjacency matrix
# seed = a seed for reproducibility
# z_init = V-vector of initialization assignment for each node (default = one cluster for each node)
# a_1,a_2 = parameters of the Gamma prior on the lambdas
# x = V-vector of categorical covariates
# if x is specified, also alpha_xi (a C-vector of parameters for a Dirichlet distribution) must be set

# Output:
# Posterior samples of the community labels for each node v=1,...,V

esbm_count <- function(Y_count, seed, N_iter, z_init=c(1:nrow(Y_count)), a_1=1, a_2=1, gamma_GN=NA, x=NULL, alpha_xi=NULL){
  
  # ----------------------------------------------
  # Select the Gnedin urn scheme
  # ----------------------------------------------

    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  
  # ----------------------------------------------
  # Pre-processing of the node attributes
  # ----------------------------------------------
  
  if (!is.null(x)){
    print("Covariates have been provided")
    x <- as.numeric(as.factor(x))
    X <- vec2mat(x)
    if (!is.null(alpha_xi)){
      alpha0 <- sum(alpha_xi)
    } else {
      stop("If covariates x are given, then alpha_xi must be set as well")
    }
  }
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  set.seed(seed)
  
  V <- nrow(Y_count)
  
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- vec2mat(z_init)
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  
  # create the matrix with block counts
  # sum weighted connections of each node with those in the different blocks
  temp_w   <- Y_count%*%Z
  # sum weighted connections in each block
  m_full_w <- t(Z)%*%temp_w - diag(0.5*colSums(temp_w*Z),ncol(Z))
  # number of node pairs in each block
  n_full_w <- outer(table(z_init),table(z_init))-diag(table(z_init))
  diag(n_full_w) <- diag(n_full_w)/2
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
    for (v in 1:V){
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full_w and n_full_w matrix
        m_full_w <- matrix(m_full_w[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
        n_full_w <- matrix(n_full_w[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      # Z_v = cluster indicator matrix without node v
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v_w = sum of weighted connections from node v to each cluster (no self-loops)
      r_v_w       <- crossprod(Z_v, Y_count[-v,v])      
      # current cluster of node v
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m_w and n_w matrix without node v
      if(length(h_v) == 1){
        resid1_w       <- matrix(0,H,H)
        resid2_w       <- matrix(0,H,H)
        resid1_w[,h_v] <- r_v_w; resid1_w[h_v,] <- r_v_w
        resid2_w[,h_v] <- v_minus; resid2_w[h_v,] <- v_minus
        m_w            <- m_full_w - resid1_w
        n_w            <- n_full_w - resid2_w
      } else {m_w <- m_full_w
      	n_w <- n_full_w} # No need to update m in this case
      
      # increase in the number pairs and weighted connections in any cluster h when adding v to this cluster
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R_w  <- matrix(1,H,1)%*%matrix(r_v_w,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      # log-likelihoods of old and new clusters
      log_lhds_old <- rowSums((a_1+m_w)*log(a_2+n_w)-(a_1+m_w+R_w)*log(a_2+n_w+V_minus)+lgamma(a_1+m_w+R_w)-lgamma(a_1+m_w)) # length H
      log_lhd_new  <- sum((a_1)*log(a_2)-(a_1+r_v_w)*log(a_2+v_minus)+lgamma(a_1+r_v_w)-lgamma(a_1)) # scalar
      log_addit    <- 0
      
      # supervision by attributes
      if(!is.null(x)){
        Vx        <- crossprod(Z_v, X[-v,])
        addit_old <- (Vx[,x[v]] + alpha_xi[x[v]]) / (v_minus+alpha0)
        addit_new <- alpha_xi[x[v]] / alpha0
        log_addit <- log(c(addit_old, addit_new))
      }
      
      # full-conditional probabilities of cluster assignments
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v_w, m_w, n_w and v_minus
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m_w             <- rbind(cbind(m_w,0),0)
        n_w             <- rbind(cbind(n_w,0),0)
        H               <- H + 1
        r_v_w           <- crossprod(Z[-v,], Y_count[-v,v])
        v_minus         <- colSums(Z[-v,])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full_w and n_full_w
      resid3_w              <- matrix(0,H,H)
      resid4_w              <- matrix(0,H,H)
      resid3_w[,new_sample] <- r_v_w; resid3_w[new_sample,] <- r_v_w
      resid4_w[,new_sample] <- v_minus; resid4_w[new_sample,] <- v_minus
      m_full_w              <- m_w + resid3_w
      n_full_w              <- n_w + resid4_w
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}

####################################################################################
####################################################################################
# GIBBS SAMPLER FOR THE ZIP STOCHASTIC BLOCK MODEL  ################################
####################################################################################
####################################################################################

# Inputs:
# Y = VxV symmetric (weighted) adjacency matrix
# seed = a seed for reproducibility
# z_init = V-vector of initialization assignment for each node (default = one cluster for each node)
# a_1,a_2 = parameters of the Gamma prior on the lambdas
# a,b = parameters of the Beta prior on the pi
# x = V-vector of categorical covariates
# if x is specified, also alpha_xi (a C-vector of parameters for a Dirichlet distribution) must be set

# Output:
# Posterior samples of the community labels for each node v=1,...,V

esbm_zip <- function(Y, seed, N_iter, z_init=c(1:nrow(Y)), a=1, b=1,a_1=1, a_2=1, gamma_GN=NA, x=NULL, alpha_xi=NULL){
  
  # ----------------------------------------------
  # Select the Gnedin urn scheme
  # ----------------------------------------------

    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  
  # ----------------------------------------------
  # Pre-processing of the node attributes
  # ----------------------------------------------
  
  if (!is.null(x)){
    print("Covariates have been provided")
    x <- as.numeric(as.factor(x))
    X <- vec2mat(x)
    if (!is.null(alpha_xi)){
      alpha0 <- sum(alpha_xi)
    } else {
      stop("If covariates x are given, then alpha_xi must be set as well")
    }
  }
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  set.seed(seed)
  
  V <- nrow(Y)
  
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- vec2mat(z_init)
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  # initialize Y_bin (X in paper notation)
  Y_bin_star <- (Y==0)*1
  Y_bin_star[Y_bin_star==1] = NA
  upperTriangle(Y_bin_star) <- 0
  
  # initialize Y_count (W in paper notation)
  Y_count_star <- Y
  Y_count_star[Y_count_star==0] <- NA
  upperTriangle(Y_count_star) <- 0

  # select edges that require data augmentation (vectorized form)
  sel_augmented <- which(is.na(lowerTriangle(Y_bin_star)))
  n_augmented   <- length(sel_augmented)
  
  # initalize pi and lambda with mean from the prior
  H   <- ncol(Z)
  pi_block     <- matrix(a/(a+b),H,H)
  lambda_block <- matrix(a_1/a_2,H,H)
  
  # compute probabilities in equation (6)
  pi_miss_block <- pi_block/(pi_block+(1-pi_block)*exp(-lambda_block))
  
  # sample the augmented data
  lambda_matr <- Z%*%lambda_block%*%t(Z)
  pi_miss_matr <- Z%*%pi_miss_block%*%t(Z)
  
  # binary matrix (X)
  Y_bin <- Y_bin_star
  lowerTriangle(Y_bin)[sel_augmented]<-rbinom(n_augmented,1,prob=lowerTriangle(pi_miss_matr)[sel_augmented])
  Y_bin <- Y_bin + t(Y_bin)
  
  # count matrix (W)
  Y_count <- Y_count_star
  lowerTriangle(Y_count)[sel_augmented]<-lowerTriangle(Y_bin)[sel_augmented]*rpois(n_augmented,lambda=lowerTriangle(lambda_matr)[sel_augmented])
  Y_count <- Y_count + t(Y_count)
  
  diag(Y_bin) <- 0
  diag(Y_count) <- 0
  
  # create the matrixes with block connections and block counts
  # sum weighted connections of each node with those in the different blocks
  temp_w   <- Y_count%*%Z
  # connections of each node with those in the different blocks
  temp     <- Y_bin%*%Z

  # sum weighted connections in each block
  m_full_w <- t(Z)%*%temp_w - diag(0.5*colSums(temp_w*Z),ncol(Z))
  # connections in each block
  m_full   <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z)) 
  # number of node pairs in each block
  n_full_w <- outer(table(z_init),table(z_init))-diag(table(z_init))
  diag(n_full_w) <- diag(n_full_w)/2
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
  	 
  	 # ----------------------------------------------
     # Update node assignments
     # ----------------------------------------------
   
    for (v in 1:V){
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full and m_full_w matrix
        m_full_w <- matrix(m_full_w[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
        n_full_w <- matrix(n_full_w[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
        m_full   <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      # Z_v = cluster indicator matrix without node v
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v_w = sum of weighted connections from node v to each cluster (no self-loops)
      r_v_w       <- crossprod(Z_v, Y_count[-v,v])      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y_bin[-v,v]) 
      # current cluster of node v
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m, m_w and n_w matrix without node v
      if(length(h_v) == 1){
        resid1_w       <- matrix(0,H,H)
        resid2_w       <- matrix(0,H,H)
        resid1         <- matrix(0,H,H)
        resid1_w[,h_v] <- r_v_w; resid1_w[h_v,] <- r_v_w
        resid2_w[,h_v] <- v_minus; resid2_w[h_v,] <- v_minus
        resid1[,h_v]   <- r_v; resid1[h_v,] <- r_v
        m_w            <- m_full_w - resid1_w
        n_w            <- n_full_w - resid2_w
        m              <- m_full - resid1
      } else {m_w <- m_full_w
      	n_w <- n_full_w
      	m <- m_full} # No need to update in this case
      
      # number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      # increase in the number pairs and weighted connections in any cluster h when adding v to this cluster
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R_w     <- matrix(1,H,1)%*%matrix(r_v_w,1,H)
      # increase in the number of edges in any cluster h when adding v to this cluster
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      # log-likelihoods of old and new clusters
      log_lhds_old <- rowSums((a_1+m_w)*log(a_2+n_w)-(a_1+m_w+R_w)*log(a_2+n_w+V_minus)+lgamma(a_1+m_w+R_w)-lgamma(a_1+m_w)) + rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # length H
      log_lhd_new  <- sum((a_1)*log(a_2)-(a_1+r_v_w)*log(a_2+v_minus)+lgamma(a_1+r_v_w)-lgamma(a_1)) + sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      log_addit    <- 0
      
      # supervision by attributes
      if(!is.null(x)){
        Vx        <- crossprod(Z_v, X[-v,])
        addit_old <- (Vx[,x[v]] + alpha_xi[x[v]]) / (v_minus+alpha0)
        addit_new <- alpha_xi[x[v]] / alpha0
        log_addit <- log(c(addit_old, addit_new))
      }
      
      # full-conditional probabilities of cluster assignments
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v_w, m_w, n_w, v_minus, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m_w             <- rbind(cbind(m_w,0),0)
        n_w             <- rbind(cbind(n_w,0),0)
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v_w           <- crossprod(Z[-v,], Y_count[-v,v])
        r_v             <- crossprod(Z[-v,], Y_bin[-v,v])
        v_minus         <- colSums(Z[-v,])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full, m_full_w and n_full_w
      resid3_w              <- matrix(0,H,H)
      resid4_w              <- matrix(0,H,H)
      resid2                <- matrix(0,H,H)
      resid3_w[,new_sample] <- r_v_w; resid3_w[new_sample,] <- r_v_w
      resid4_w[,new_sample] <- v_minus; resid4_w[new_sample,] <- v_minus
      resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
      m_full_w              <- m_w + resid3_w
      n_full_w              <- n_w + resid4_w
      m_full                <- m + resid2
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    # ----------------------------------------------
    # pi block and lambda block
    # ----------------------------------------------
    
    H   <- ncol(Z)
 
    # update block pi from Beta full conditional
    pi_block <- matrix(0,H,H)
    lowerTriangle(pi_block,diag=TRUE)<-rbeta(H*(H-1)/2+H,lowerTriangle(m_full,diag=TRUE)+a,lowerTriangle(n_full_w,diag=TRUE)-lowerTriangle(m_full,diag=TRUE)+b)
    upperTriangle(pi_block,byrow=TRUE)<-lowerTriangle(pi_block)
    
    # update block lambda from Gamma full conditional
    lambda_block <- matrix(0,H,H)
    lowerTriangle(lambda_block,diag=TRUE)<-rgamma(H*(H-1)/2+H,lowerTriangle(m_full_w,diag=TRUE)+a_1,lowerTriangle(n_full_w,diag=TRUE)+a_2)
    upperTriangle(lambda_block,byrow=TRUE)<-lowerTriangle(lambda_block)

  
    # ----------------------------------------------
    # update augmented data
    # ----------------------------------------------
 
     pi_miss_block <- pi_block/(pi_block+(1-pi_block)*exp(-lambda_block))
  
     lambda_matr <- Z%*%lambda_block%*%t(Z)
     pi_miss_matr <- Z%*%pi_miss_block%*%t(Z)
  
     Y_bin <- Y_bin_star
     lowerTriangle(Y_bin)[sel_augmented]<-rbinom(n_augmented,1,prob=lowerTriangle(pi_miss_matr)[sel_augmented])
     Y_bin <- Y_bin + t(Y_bin)
  
     Y_count <- Y_count_star
     lowerTriangle(Y_count)[sel_augmented]<-lowerTriangle(Y_bin)[sel_augmented]*rpois(n_augmented,lambda=lowerTriangle(lambda_matr)[sel_augmented])
     Y_count <- Y_count + t(Y_count)
  
     diag(Y_bin) <- 0
     diag(Y_count) <- 0

    # ----------------------------------------------
    # update count matricess
    # ----------------------------------------------
     temp_w   <- Y_count%*%Z
     temp     <- Y_bin%*%Z

     m_full_w <- t(Z)%*%temp_w - diag(0.5*colSums(temp_w*Z),ncol(Z))
     m_full   <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
     n_full_w <- outer(table(z_post[,t]),table(z_post[,t]))-diag(table(z_post[,t]))
     diag(n_full_w) <- diag(n_full_w)/2

    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}

####################################################################################
####################################################################################
# PUT CLUSTER LABELS IN BINARY MATRIX FORM  ########################################
####################################################################################
####################################################################################

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

####################################################################################
####################################################################################
# COMPUTE POSTERIOR CO-CLUSTERING MATRIX  ##########################################
####################################################################################
####################################################################################

pr_cc <- function(z_post){
  # in: posterior sample of assignments (VxN_iter matrix)
  # out: VxV matrix c with elements c[vu]=fraction of iterations in which v and u are in the same cluster
  V <- nrow(z_post)    
  N_iter <- ncol(z_post)
  c <- matrix(0,V,V)
  for (t in 1:N_iter){
    Z <- vec2mat(z_post[,t])
    c <- c + Z%*%t(Z)
  }
  return(c/N_iter)
}

####################################################################################
####################################################################################
# SAMPLE THE BLOCK PROBABILITIES ###################################################
####################################################################################
####################################################################################

sample_pi <- function(z, Y_bin, seed, N_iter, a, b){
 # in: one vector of node labels z, VxV adjancency matrix (Y_bin), seed, number of samples, hyperparameters beta
 # out: HxHxN_iter matrix of posterior samples for pi given z
  
  set.seed(seed)
  memb <- z
  z <- vec2mat(z)
  H <- ncol(z)
  V <- dim(Y_bin)[1]
  
  M <- t(z)%*%Y_bin%*%z
  diag(M) <- diag(M)/2
  Tot <- t(z)%*%matrix(1,V,V)%*%z
  diag(Tot) <- (diag(Tot)-table(memb))/2
  Mbar <- Tot-M
  a_n <- lowerTriangle(M,diag=TRUE)+a
  b_bar_n <- lowerTriangle(Mbar,diag=TRUE)+b
  
  pi_post <- array(0,dim=c(H,H,N_iter))
  
  for (t in 1:N_iter){
  theta <- rbeta(length(a_n),a_n,b_bar_n)
  Theta <- matrix(0,H,H)
  Theta[lower.tri(Theta,diag=TRUE)] <- theta
  Theta <- Theta+t(Theta)
  diag(Theta) <- diag(Theta)/2
  pi_post[,,t] <- Theta
  }
    
  return(pi_post)
}

####################################################################################
####################################################################################
# SAMPLE THE BLOCK LAMBDAS #########################################################
####################################################################################
####################################################################################

sample_lambda <- function(z, Y_count, seed, N_iter, a_1, a_2){
 # in: one vector of node labels z, VxV adjancency matrix (Y_count), seed, number of samples, hyperparameters gamma
 # out: HxHxN_iter matrix of posterior samples for lambda given z
  
  set.seed(seed)
  memb <- z
  z <- vec2mat(z)
  H <- ncol(z)
  V <- dim(Y_count)[1]
  
  M <- t(z)%*%Y_count%*%z
  diag(M) <- diag(M)/2
  N <- outer(table(memb),table(memb))-diag(table(memb))
  diag(N) <- diag(N)/2

  a_1_update <- lowerTriangle(M,diag=TRUE)+a_1
  a_2_update <- lowerTriangle(N,diag=TRUE)+a_2
  
  lambda_post <- array(0,dim=c(H,H,N_iter))
  
  for (t in 1:N_iter){
  theta <- rgamma(length(a_1_update),a_1_update,a_2_update)
  Theta <- matrix(0,H,H)
  Theta[lower.tri(Theta,diag=TRUE)] <- theta
  Theta <- Theta+t(Theta)
  diag(Theta) <- diag(Theta)/2
  lambda_post[,,t] <- Theta
  }
    
  return(lambda_post)
}

####################################################################################
####################################################################################
# SAMPLE THE BLOCK LAMBDAS AND BLOCK PROBABILITIES #################################
####################################################################################
####################################################################################

sample_lambda_pi <- function(z, Y, seed, N_iter, a, b, a_1, a_2){
  # in: one vector of node labels z, VxV adjancency matrix (Y), seed, number of samples, hyperparameters
  # out: HxHxN_iter matrixes of posterior samples for lambda and pi given z
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  set.seed(seed)
  
  V <- nrow(Y)  
  Z <- vec2mat(z)
  
  # initialize Y_bin (X in paper notation)
  Y_bin_star <- (Y==0)*1
  Y_bin_star[Y_bin_star==1] = NA
  upperTriangle(Y_bin_star) <- 0
  
  # initialize Y_count (W in paper notation)
  Y_count_star <- Y
  Y_count_star[Y_count_star==0] <- NA
  upperTriangle(Y_count_star) <- 0

  # select edges that require data augmentation (vectorized form)
  sel_augmented <- which(is.na(lowerTriangle(Y_bin_star)))
  n_augmented   <- length(sel_augmented)
  
  # initalize pi and lambda with mean from the prior
  H   <- ncol(Z)
  pi_block     <- matrix(a/(a+b),H,H)
  lambda_block <- matrix(a_1/a_2,H,H)
  
  # compute probabilities in equation (6)
  pi_miss_block <- pi_block/(pi_block+(1-pi_block)*exp(-lambda_block))
  
  # sample the augmented data
  lambda_matr <- Z%*%lambda_block%*%t(Z)
  pi_miss_matr <- Z%*%pi_miss_block%*%t(Z)
  
  # binary matrix (X)
  Y_bin <- Y_bin_star
  lowerTriangle(Y_bin)[sel_augmented]<-rbinom(n_augmented,1,prob=lowerTriangle(pi_miss_matr)[sel_augmented])
  Y_bin <- Y_bin + t(Y_bin)
  
  # count matrix (W)
  Y_count <- Y_count_star
  lowerTriangle(Y_count)[sel_augmented]<-lowerTriangle(Y_bin)[sel_augmented]*rpois(n_augmented,lambda=lowerTriangle(lambda_matr)[sel_augmented])
  Y_count <- Y_count + t(Y_count)
  
  diag(Y_bin) <- 0
  diag(Y_count) <- 0
  
  # create the matrixes with block connections and block counts
  # sum weighted connections of each node with those in the different blocks
  temp_w   <- Y_count%*%Z
  # connections of each node with those in the different blocks
  temp     <- Y_bin%*%Z

  # sum weighted connections in each block
  m_full_w <- t(Z)%*%temp_w - diag(0.5*colSums(temp_w*Z),ncol(Z))
  # connections in each block
  m_full   <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z)) 
  # number of node pairs in each block
  n_full_w <- outer(table(z),table(z))-diag(table(z))
  diag(n_full_w) <- diag(n_full_w)/2

  # create pi_post and lambda_post matrices containing samples of pi and lambda
  pi_post <- array(0,dim=c(H,H,N_iter))
  lambda_post <- array(0,dim=c(H,H,N_iter))


  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
 
   for (t in 1:N_iter){
  	 
    # ----------------------------------------------
    # pi_block and lambda_block
    # ----------------------------------------------
     
    # update block pi from Beta full conditional
    pi_block <- matrix(0,H,H)
    lowerTriangle(pi_block,diag=TRUE)<-rbeta(H*(H-1)/2+H,lowerTriangle(m_full,diag=TRUE)+a,lowerTriangle(n_full_w,diag=TRUE)-lowerTriangle(m_full,diag=TRUE)+b)
    upperTriangle(pi_block,byrow=TRUE)<-lowerTriangle(pi_block)
    pi_post[,,t] <- pi_block
    
    # update block lambda from Gamma full conditional
    lambda_block <- matrix(0,H,H)
    lowerTriangle(lambda_block,diag=TRUE)<-rgamma(H*(H-1)/2+H,lowerTriangle(m_full_w,diag=TRUE)+a_1,lowerTriangle(n_full_w,diag=TRUE)+a_2)
    upperTriangle(lambda_block,byrow=TRUE)<-lowerTriangle(lambda_block)
    lambda_post[,,t] <- lambda_block
  
    # ----------------------------------------------
    # update augmented data
    # ----------------------------------------------
 
     pi_miss_block <- pi_block/(pi_block+(1-pi_block)*exp(-lambda_block))
  
     lambda_matr <- Z%*%lambda_block%*%t(Z)
     pi_miss_matr <- Z%*%pi_miss_block%*%t(Z)
  
     Y_bin <- Y_bin_star
     lowerTriangle(Y_bin)[sel_augmented]<-rbinom(n_augmented,1,prob=lowerTriangle(pi_miss_matr)[sel_augmented])
     Y_bin <- Y_bin + t(Y_bin)
  
     Y_count <- Y_count_star
     lowerTriangle(Y_count)[sel_augmented]<-lowerTriangle(Y_bin)[sel_augmented]*rpois(n_augmented,lambda=lowerTriangle(lambda_matr)[sel_augmented])
     Y_count <- Y_count + t(Y_count)
  
     diag(Y_bin) <- 0
     diag(Y_count) <- 0

    # ----------------------------------------------
    # update count matrices
    # ----------------------------------------------
     temp_w   <- Y_count%*%Z
     temp     <- Y_bin%*%Z

     m_full_w <- t(Z)%*%temp_w - diag(0.5*colSums(temp_w*Z),ncol(Z))
     m_full   <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
     n_full_w <- outer(table(z),table(z))-diag(table(z))
     diag(n_full_w) <- diag(n_full_w)/2

}
 
return(list(pi_post,lambda_post))}
 
####################################################################################
####################################################################################
# COMPUTE THE DISTRIBUTION OF H AND THE EXPECTED VALUE UNDER the GNEDIN ############
####################################################################################
####################################################################################

# ----------------------------------------------
# GNEDIN
# ----------------------------------------------

HGnedin <- function(V, h, gamma=0.5){
  exp(lchoose(V, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(V+ gamma - h) - lgamma(V +gamma))
}
