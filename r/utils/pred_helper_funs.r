library(Rcpp)
library(inline)
# library(bayesSurv)

# source('r/utils/simDataFuns.r')
# source('prediction/r/utils/pred_plot_funs.r')

# build c matrix
cppFunction('
  NumericMatrix build_c(double rho, double tau, double eta, 
                        NumericMatrix d_inter, NumericMatrix lag) {
    int T=lag.nrow(), N=d_inter.nrow(), N_knots=d_inter.ncol();
    NumericMatrix c(N*T, N_knots*T);
    for (int s = 0; s<N; s++)
      for (int t = 0; t<T; t++)
        for (int u = 0; u<N_knots; u++)
          for (int v = 0; v<T; v++) {
            //std::cout << s << " " << t << " " << u << " " << v << " " << std::endl;
            c(s*T+t,u*T+v) = eta*eta*exp(-d_inter(s,u)/rho)*exp(-1/tau*lag(t,v));
          }
    return c;
  }
', verbose=TRUE, showOutput=TRUE)

# r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
# additive log_ratio transformation
cppFunction('
  NumericMatrix sum2one_constraint(int K, int N, int T, NumericMatrix g, NumericVector sum_exp_g) {
    //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
    NumericMatrix r(N*T, K);
    for (int k = 0; k<(K-1); k++)
      for (int j = 0; j<N*T; j++)
        r(j,k) = exp(g(j,k))/(1+sum_exp_g(j));

    for (int j = 0; j<N*T; j++)
      r(j,K-1) = 1.0 / (1 + sum_exp_g[j]);
   
return r;
  }
', verbose=TRUE)


# additive log_ratio transformation
cppFunction('
  NumericMatrix sum2one_constraint_nb(int K, int N, int T, NumericMatrix g, NumericVector sum_exp_g) {
    //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
    NumericMatrix r(N*T, K);
    for (int k = 0; k<K; k++)
      for (int j = 0; j<N*T; j++)
        r(j,k) = exp(g(j,k))/(sum_exp_g(j));
   
return r;
  }
', verbose=TRUE)

# # matrix vector mult
# cppFunction('
#   NumericMatrix mat_vec_mult(int N, int N_knots, NumericMatrix H, NumericVector alpha) {
#     //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
#     NumericVector Halpha(N);
#     
#     for (int i = 0; i<N; i++)
#         Halpha(i) = 0;
#       for (int j = 0; j<N_knots; j++)
#         Halpha(i) += H(i,j) * alpha(j);    
# 
# return Halpha;
#   }
# ', verbose=TRUE)

load_stan_output <- function(suff_fit){
  # if (!file.exists(paste0('output/', suff_fit,'.rdata'))){
  fname     = sprintf('output/%s.bin', suff_fit)
  object    = read_stanbin(fname) 

  # find first parameter
  col0 = 1
  for (i in 1:10) {
    if (grepl("__", object$params[i], fixed=TRUE) == FALSE) {
      col0 = i
      break
    }
  }
  
  samples   = data.frame(object$samples[,col0:ncol(object$samples)], object$samples[,1])
  post      = array(0, c(nrow(samples), 1, ncol(samples)))   
  post[,1,] = as.matrix(samples) 
  dimnames(post)[[3]] = colnames(samples)
  
  return(list(post=post, par_names=object$par_names))
}



build_r <- function(post_dat, N, T, K){
  
  if (length(dim(post_dat$post)) == 3) {
    post      = post_dat$post[,1,]
  } else {
    post      = post_dat$post
  }
  
  par_names = post_dat$par_names
  
  W       = K-1
  niter   = nrow(post) 
  
  g      = array(NA, dim=c(N*T, W, niter))
  r      = array(NA, dim=c(N*T, K, niter))
  
  # fix the stan2bin!
  g_cols = which(par_names == 'g')#-6
  
  for (k in 1:W){
    
    print(k)
    
    g_k_cols = seq(k, T*N*W, by=W)
    g[,k,]   = t(post[,g_cols[g_k_cols]])    

  }
  
  print("Log-ratio transforming g")
  
  for (i in 1:niter){
    
    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
    
    sum_exp_g = rowSums(exp(g[,,i]))
    
    #ti <- proc.time()
    
    # additive log-ratio transformation
    # for (k in 1:W)
    #   for (j in 1:(N*T))
    #     r[j,k,i] <- exp(g[j,k,i]) / (1 + sum_exp_g[j])
    
    # for (j in 1:(N*T))
    #   r[j,K,i] <- 1 / (1 + sum_exp_g[j])
    # #     
    # tj <- proc.time()
    # print("OLD Build r:")
    # print(tj-ti)
       
    r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
    # tk <- proc.time()
    # print("NEW Build r:")
    # print(tk-tj)
  }
  
  return(list(r=r, g=g))
}

build_r_nb <- function(post_dat, N, T, K){
  
  if (length(dim(post_dat$post)) == 3) {
    post      = post_dat$post[,1,]
  } else {
    post      = post_dat$post
  }
  
  par_names = post_dat$par_names
  
  niter   = nrow(post) 
  
  g      = array(NA, dim=c(N*T, K, niter))
  r      = array(NA, dim=c(N*T, K, niter))
  
  # fix the stan2bin!
  g_cols = which(par_names == 'g')#-6
  
  for (k in 1:K){
    
    print(k)
    
    g_k_cols = seq(k, T*N*K, by=K)
    g[,k,]   = t(post[,g_cols[g_k_cols]])    
    
  }
  
  print("Log-ratio transforming g")
  
  for (i in 1:niter){
    
    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
    
    sum_exp_g = rowSums(exp(g[,,i]))
    
    r[,,i] <- sum2one_constraint_nb(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
    # tk <- proc.time()
    # print("NEW Build r:")
    # print(tk-tj)
  }
  
  return(list(r=r, g=g))
}

build_r_from_mu_g <- function(mu_g, N, T, K){
  
  niter = dim(mu_g)[3]
  r     = array(NA, dim=c(N*T, K, niter))
  
  print("Log-ratio transforming g")
  
  for (i in 1:niter){
    
    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
    
    sum_exp_g = rowSums(exp(mu_g[,,i]))       
    r[,,i] <- sum2one_constraint(K, N, T, as.matrix(mu_g[,,i]), sum_exp_g) 
  }
  
  return(r)
}

build_r_from_mu_g_nb <- function(mu_g, N, T, K){

  niter = dim(mu_g)[3]
  r     = array(NA, dim=c(N*T, K, niter))
  
  print("Log-ratio transforming g")
  
  for (i in 1:niter){
    
    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
    
    sum_exp_g = rowSums(exp(mu_g[,,i]))       
    r[,,i] <- sum2one_constraint_nb(K, N, T, as.matrix(mu_g[,,i]), sum_exp_g) 
  }
  
  return(r)
}

# build_r <- function(post_dat, N, T, K){
#   
#   if (length(dim(post_dat$post)) == 3) {
#     post      = post_dat$post[,1,]
#   } else {
#     post      = post_dat$post
#   }
#   
#   par_names = post_dat$par_names
#   
#   W       = K-1
#   niter   = nrow(post) 
#   
#   g      = array(NA, dim=c(N*T, W, niter))
#   r      = array(NA, dim=c(N*T, K, niter))
#   
#   # fix the stan2bin!
#   g_cols = which(par_names == 'g')#-6
#   
#   for (k in 1:W){
#     
#     print(k)
#     g_k_cols = seq(k, T*N*W, by=W)
#     g[,k,]   = t(post[,g_cols[g_k_cols]])    
#   }
#   for (i in 1:niter){
#     
#     if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
#     
#     sum_exp_g = rowSums(exp(g[,,i]))
#     r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
#   }
#   
#   return(list(r=r, g=g))
# }



build_mu_g_serial <- function(post_dat, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0){
  
  post      = post_dat$post
  par_names = post_dat$par_names
  
  N       = nrow(d_inter)
  N_knots = ncol(d_inter)
  niter   = dim(post[,1,])[1] 
  W       = K-1
  
  ones     = matrix(1, nrow=N, ncol=1)
  mu_g     = array(NA, dim=c(N*T, W, niter))
  Halpha_s = array(NA, dim=c(N, W, niter))
  
  if (mu0){
    Halpha_t = array(NA, dim=c(N, W*(T-1), niter))
  } else {
    Halpha_t = array(NA, dim=c(N, W*T, niter))
  }
  #sumHalpha = array(NA, dim=c(T, W, niter))
  
  print("Done allocating")
  
#   col_names  = colnames(post[,1,])
#   par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\.")[[1]][1]))  
#   #   col_substr = sapply(strsplit(col_names, "\\["), function(x) x[1])
#   #   par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))  

  ksi    = post[,1, which(par_names == 'ksi')]
  
  alpha_s_start = min(which(par_names == 'alpha_s'))
  alpha_t_start = min(which(par_names == 'alpha_t'))
  
  if (od) {
    temp = qr(ones)
    Q    = qr.Q(temp)
    #     R = qr.R(temp)
      
    P = Q %*% t(Q)
  }
  
  mu   = post[,1,which(par_names == 'mu')]
  mu_t = post[,1,which(par_names == 'mu_t')]

#  for (k in 1:W){
  for (k in 1:2){
    print(k)
    
    #mu     = post[,1,which(par_names == 'mu')[k]]
    sigma  = post[,1,which(par_names == 'sigma')[k]]
    lambda = post[,1,which(par_names == 'lambda')[k]]
    
    if (mu0){
      mut_cols = seq(k, (T-1)*W, by=W) 
    } else {
      mut_cols = seq(k, T*W, by=W)
    }
#     col_names[which(par_names == 'mu_t')[mut_cols]]
    mu_t_k = mu_t[,mut_cols]
    
    alpha_s_cols = seq(alpha_s_start + k - 1, alpha_s_start + N_knots*W - 1, by=W)
#     col_names[alpha_s_cols]
    alpha_s = post[,1,alpha_s_cols]

#     alpha_s_t = t(post[,1,alpha_s_cols])
    
    C_s <- exp(-d_knots/rho[k])
    c_s <- exp(-d_inter/rho[k])
    C_s_inv = chol2inv(chol(C_s))
    
    if (od) {
      cs_Csinv = c_s %*% C_s_inv - P %*% c_s %*% C_s_inv
    } else {
      cs_Csinv = c_s %*% C_s_inv
    }
      
    # t=1
    mu_g_idx = seq(1, N*T, by=T)

    for (i in 1:niter){
     
      if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}	      

      ts <- proc.time()  
     
      
#       t1 <- proc.time()
      Halpha_s[,k,i] = cs_Csinv %*% alpha_s[i,]
#       t2 <- proc.time()
#       Halpha_s[,k,i] = cs_Csinv %*% alpha_s_t[,i]
#       t3 <- proc.time()
#       print("OLD mult : ")
#       print(t2-t1)
#       print("NEW mult : ")
#       print(t3-t2)
      
      
      if (mu0){
        mu_g[mu_g_idx,k,i] = mu[i,k] + Halpha_s[,k,i] 
      } else {
        mu_g[mu_g_idx,k,i] = mu[i,k] + mu_t_k[i,1] + Halpha_s[,k,i] 
      }
      
      Q <- exp(-d_knots/lambda[i])
      q <- exp(-d_inter/lambda[i])
      Q_inv = chol2inv(chol(Q))
      
      q_Qinv = q %*% Q_inv
      if (i == 4) {
            print((cs_Csinv %*% alpha_s[i,])[1:10])
}	    

      for (t in 2:T){
        
        alpha_t_cols = seq(alpha_t_start + (k-1)*(T-1) + t-1 - 1, alpha_t_start + N_knots*W*(T-1) - 1, by=W*(T-1))
#         col_names[alpha_t_cols]
        alpha_t = post[,1,alpha_t_cols]
        
        mu_g_idx = seq(t, N*T, by=T)
        if (mu0){
          Halpha_t[,(k-1)*(T-1) + t-1,i] = q_Qinv %*% alpha_t[i,] 
          mu_g[mu_g_idx,k,i] = mu[i,k] + mu_t_k[i,t-1] + cs_Csinv %*% alpha_s[i,] + Halpha_t[,(k-1)*(T-1) + t-1,i] #q_Qinv %*% alpha_t[i,] 
        } else {
          
#           t1 <- proc.time()
          mu_g[mu_g_idx,k,i] = mu[i,k] + mu_t_k[i,t] + cs_Csinv %*% alpha_s[i,] + q_Qinv %*% alpha_t[i,] 
#           t2 <- proc.time()
#           mu_g[mu_g_idx,k,i] = mu[i,k] + mu_t_k[i,t-1] + mat_vec_mult(N, N_knots, cs_Csinv, alpha_s[i,]) + mat_vec_mult(N, N_knots, q_Qinv, alpha_t[i,] )
#         
#                 t3 <- proc.time()
#                 print("OLD mult : ")
#                 print(t2-t1)
#                 print("NEW mult : ")
#                 print(t3-t2)
        }
        
      }

#       te <- proc.time()
#       print("Single iteration:")
#       print(te-ts)
    
    }
  }
  
  return(list(mu_g=mu_g, mu=mu, mu_t=mu_t, Halpha_s=Halpha_s, Halpha_t=Halpha_t))
}

# write mean vals to file
write_par_vals <- function(post_dat, taxa, subDir, N_pars){
  
  if (length(dim(post_dat$post)) == 3) {
    post      = post_dat$post[,1,]
  } else {
    post      = post_dat$post
  }
  
  par_names = post_dat$par_names
  
  ess = apply(post[,c(1:N_pars,dim(post)[3])], 2, ess_rfun)
  
  # sink(sprintf('%s/%s/summary.txt', wd, path_figs1), type='output')
  sink(sprintf('%s/summary.txt', subDir), type='output')
  print('The taxa modelled are:')
  print(taxa)
  cat('\n')
  print('Summary of posterior parameter vals:')
  print(get_quants(post, N_pars))
  cat('\n')
  print('ESS:')
  print(ess)
  sink()
  
}

build_props_full <- function(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0, res){
  
  N = nrow(d_inter)
  N_knots = ncol(d_inter)
  niter   = dim(post[,1,])[1] 
  
  W = K-1
  
  mu_g = array(NA, dim=c(N*T, W, niter))
  g    = array(NA, dim=c(N*T, W, niter))
  r    = array(NA, dim=c(N*T, K, niter))
  #   Halpha = array(NA, dim=c(N_knots, W, niter))
  ones = matrix(1, nrow=N*T, ncol=1)
  
  Halpha_s = array(NA, dim=c(N, W, niter))
  
  if (mu0){
    Halpha_t = array(NA, dim=c(N, W*(T-1), niter))
  } else {
    Halpha_t = array(NA, dim=c(N, W*T, niter))
  }
  #sumHalpha = array(NA, dim=c(T, W, niter))
  
  print("Done allocating")
  
  col_names  = colnames(post[,1,])
#   col_substr = sapply(strsplit(col_names, "\\["), function(x) x[1])
#   par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\.")[[1]][1]))  

  ksi    = post[,1, which(par_names == 'ksi')]
  
  alpha_s_start = min(which(par_names == 'alpha_s'))
  alpha_t_start = min(which(par_names == 'alpha_t'))
   
#   if (od){
#     x = matrix(1, nrow=(N*T), ncol=1)
#     N_p = N*T
#     
#     temp = qr(x)
#     Q = qr.Q(temp)
#     #     R = qr.R(temp)
#     
#     #     P = Q %*% t(Q)
#   }
  
  for (k in 1:W){
    print(k)
    
    mu     = post[,1,which(par_names == 'mu')[k]]
    sigma  = post[,1,which(par_names == 'sigma')[k]]
    lambda = post[,1,which(par_names == 'lambda')[k]]
    
    if (mu0){
      mut_cols = seq(k, (T-1)*W) 
    } else {
      mut_cols = seq(k, T*W, by=W)
    }
    col_names[which(par_names == 'mu_t')[mut_cols]]
    mu_t = post[,1,which(par_names == 'mu_t')[mut_cols]]
    
    alpha_s_cols = seq(alpha_s_start + k - 1, alpha_s_start + N_knots*W - 1, by=W)
    col_names[alpha_s_cols]
    alpha_s = post[,1,alpha_s_cols]
    
    #     g    = array(NA, dim=c(N*T, W, niter))
    g_cols = seq(k, T*N*W, by=W)
    col_names[which(par_names == 'g')][g_cols]
    g[,k,] = t(post[,1,which(par_names == 'g')[g_cols]])

    C_s <- exp(-d_knots/rho[k])
    c_s <- exp(-d_inter/rho[k])
    C_s_inv = chol2inv(chol(C_s))
    
    cs_Csinv = c_s %*% C_s_inv

    for (i in 1:niter){
      
      # t=1
      mu_g_idx = seq(1, N*T, by=T)
      
      Halpha_s[,k,i] = cs_Csinv %*% alpha_s[i,]
      
      if (mu0){
        mu_g[mu_g_idx,k,i] = mu[i] + Halpha_s[,k,i] #cs_Csinv %*% alpha_s[i,]
      } else {
        mu_g[mu_g_idx,k,i] = mu[i] + mu_t[i,1] + Halpha_s[,k,i] #cs_Csinv %*% alpha_s[i,]
      }
        
      Q <- exp(-d_knots/lambda[i])
      q <- exp(-d_inter/lambda[i])
      Q_inv = chol2inv(chol(Q))
      
      q_Qinv = q %*% Q_inv
    
      for (t in 2:T){
      
        alpha_t_cols = seq(alpha_t_start + (k-1)*(T-1) + t-1 - 1, alpha_t_start + N_knots*W*(T-1) - 1, by=W*(T-1))
        col_names[alpha_t_cols]
        alpha_t = post[,1,alpha_t_cols]
        
        mu_g_idx = seq(t, N*T, by=T)
        if (mu0){
          Halpha_t[,(k-1)*(T-1) + t-1,i] = q_Qinv %*% alpha_t[i,] 
          mu_g[mu_g_idx,k,i] = mu[i] + mu_t[i,t] + cs_Csinv %*% alpha_s[i,] + Halpha_t[,(k-1)*(T-1) + t-1,i] #q_Qinv %*% alpha_t[i,] 
        } else {
          mu_g[mu_g_idx,k,i] = mu[i] + mu_t[i,t-1] + cs_Csinv %*% alpha_s[i,] + q_Qinv %*% alpha_t[i,] 
        }
      
      }
    }
  }
  
  print("Log-ratio transforming g")
  for (i in 1:niter){
    
    sum_exp_g = rowSums(exp(g[,,i]))
    
    #     ti <- proc.time()
    
    # additive log-ratio transformation
    # r    = array(NA, dim=c(N*T, K, niter))
    for (k in 1:W)
      for (j in 1:(N*T))
        r[j,k,i] <- exp(g[j,k,i]) / (1 + sum_exp_g[j])
    
    for (j in 1:(N*T))
      r[j,K,i] <- 1 / (1 + sum_exp_g[j])
    #     
    #     tj <- proc.time()
    #     print("OLD Build r:")
    #     print(tj-ti)
    #     
    #     r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
    #     tk <- proc.time()
    #     print("NEW Build r:")
    #     print(tk-tj)
  }
  
  return(list(mu_g=mu_g, r=r, g=g, Halpha_s=Halpha_s, Halpha_t=Halpha_t))
}


get_corrections <- function(post_dat, rho, eta, T, K, d_inter, d_knots){
  
  post      = post_dat$post
  par_names = post_dat$par_names
  
  N     = nrow(d_inter)
  niter = dim(post[,1,])[1] 
  
  W = K-1

  adj_s = array(NA, dim=c(N, W))
  adj_t = array(NA, dim=c(N, W))
  
  sigma  = colMeans(post[,1,which(par_names == 'sigma')])
  lambda = colMeans(post[,1,which(par_names == 'lambda')])
  
  for (k in 1:W){

    C_s <- exp(-d_knots/rho[k])
    c_s <- exp(-d_inter/rho[k])
    C_s_inv = chol2inv(chol(C_s))
    
    adj_s[,k] = eta[k]^2 * (1 - diag(c_s %*% C_s_inv %*% t(c_s)))

    Q <- exp(-d_knots/lambda[k])
    q <- exp(-d_inter/lambda[k])
    Q_inv = chol2inv(chol(Q))
    
    adj_t[,k] = sigma[k]^2 * (1 - diag(q %*% Q_inv %*% t(q)))
  }
  
  return(list(adj_s=adj_s, adj_t=adj_t))
}

get_corrections_ar_nb <- function(post_dat, rho, eta, T, K, d_inter, d_knots){
  
  post      = post_dat$post
  par_names = post_dat$par_names
  
  N     = nrow(d_inter)
  niter = dim(post)[1] 
  
  adj_s = array(NA, dim=c(N, K))
  adj_t = array(NA, dim=c(N, K))
  
  omega  = colMeans(post[,which(par_names == 'omega')])
  
  for (k in 1:W){
    
    C_s <- exp(-d_knots/rho[k])
    c_s <- exp(-d_inter/rho[k])
    C_s_inv = chol2inv(chol(C_s))
    
    adj_s[,k] = eta[k]^2 * (1 - diag(c_s %*% C_s_inv %*% t(c_s)))
    adj_t[,k] = (1-omega[k]^2) * diag(c_s %*% C_s_inv %*% t(c_s))
  }
  
  return(list(adj_s=adj_s, adj_t=adj_t))
}


get_mut <- function(post, W){
  
  niter   = dim(post[,1,])[1] 
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\.")[[1]][1]))
  
  idx_mut = which(par_names == 'mu_t')
  
  mu_t = array(NA, c(T-1, W, niter))
  
  for (k in 1:W){
    print(k)
    idx_taxon = seq(k, W*(T-1), by=W)
    idx = idx_mut[idx_taxon]
    colnames(post[,1,])[idx]
    mu_t[,k,] = t(post[,1,idx])
  }
  return(mu_t)
}

get_mu <- function(post, W){
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\.")[[1]][1]))
  
  idx_mu = which(par_names == 'mu')
  
  mu = t(post[,1,idx_mu])
  
  return(mu)
}

# build_pollen_counts(tmin=tmin, tmax=tmax, int=int, pollen_ts=pollen_ts3, taxa_all, taxa_sub, bacon)
# build pollen counts
build_pollen_counts <- function(tmin, tmax, int, pollen_ts, taxa_all, taxa_sub, age_model){
  
  if (age_model=='bacon') {
    age_col = 'age_bacon'
  } else if  (age_model=='bchron') {
    age_col = 'age_bchron'
  } else {
    age_col = 'age_default'
  }
  
  taxa.start.col = min(match(taxa_all, colnames(pollen_ts)), na.rm=TRUE)

  if (int > 0){
    #   breaks = seq(0,2500,by=int)
    breaks = seq(tmin,tmax,by=int)
    # breaks[1] = min(pollen_ts$age_bacon)
  
    meta_pol  = pollen_ts[which((pollen_ts[, age_col] >= tmin) &
                                (pollen_ts[, age_col] <= tmax)),1:(taxa.start.col-1)]
    counts = pollen_ts[which((pollen_ts[, age_col] >= tmin) &
                             (pollen_ts[, age_col] <= tmax)),taxa.start.col:ncol(pollen_ts)]
    # meta_pol  = pollen_ts[which(pollen_ts[, age_col] <= tmax),1:(taxa.start.col-1)]
    # counts = pollen_ts[which(pollen_ts[, age_col] <= tmax),taxa.start.col:ncol(pollen_ts)]
    
    meta_pol = data.frame(meta_pol, age=rep(NA, nrow(meta_pol)))
    
    counts_agg = matrix(NA, nrow=0, ncol=ncol(counts))
    colnames(counts_agg) = colnames(counts)
    
    meta_agg = matrix(NA, nrow=0, ncol=6)
    colnames(meta_agg) = colnames(meta_pol)[1:6]
    
    meta_all = matrix(NA, nrow=0, ncol=ncol(meta_pol))
    colnames(meta_all) = colnames(meta_pol)
  
    ids = unique(meta_pol$id)
    ncores = length(ids)
  
    for (i in 1:ncores){
      
      print(i)
      
      #print(i)
      core_rows = which(meta_pol$id == ids[i])
      #     core_counts = counts[core_rows,]
    
      for (j in 1:(length(breaks)-1)){
        
        #print(j)
        age = breaks[j] + int/2
      
        age_rows = core_rows[(meta_pol[core_rows, age_col] >= breaks[j]) & 
                             (meta_pol[core_rows, age_col] < breaks[j+1])]
      
        if (length(age_rows)>1){
        
          counts_agg = rbind(counts_agg, colSums(counts[age_rows, ]))
        
          meta_agg      = rbind(meta_agg, data.frame(meta_pol[age_rows[1],1:6], age=age/100, zero=FALSE))
          
          meta_all_row = data.frame(meta_pol[age_rows,])
          meta_all_row$age = rep(age/100)
          meta_all = rbind(meta_all, meta_all_row)
        
      } else if (length(age_rows) == 1){
        
          counts_agg = rbind(counts_agg, counts[age_rows, ])
        
          meta_agg      = rbind(meta_agg, data.frame(meta_pol[age_rows,1:6], age=rep(age/100), zero=FALSE))
          # meta_agg$age[nrow(meta_agg)] = age/100
          
          meta_all_row = data.frame(meta_pol[age_rows,])
          meta_all_row$age = rep(age/100)
          meta_all = rbind(meta_all, meta_all_row)
          
      } else if (length(age_rows) == 0){
          
          #FIX ME
          counts_agg = rbind(counts_agg, rep(0,ncol(counts_agg)))
        
          meta_row = meta_pol[core_rows[1],1:6]
          meta_agg      = rbind(meta_agg, data.frame(meta_row, age=age/100, zero=TRUE))
          # meta_agg$age[nrow(meta_agg)] = age/100
          
          meta_all_row = data.frame(meta_pol[core_rows[1],])
          meta_all_row$age = NA
          meta_all_row$age_bacon = NA
          meta_all_row$age_default = NA
          
          meta_all = rbind(meta_all, meta_all_row)
         
        }
      
      }
    }
  # } else if (int==0){
  #   breaks = seq(tmin,tmax,by=int)
  #   # breaks[1] = min(pollen_ts$age_bacon)
  #   
  #   meta_pol  = pollen_ts[which((pollen_ts[, age_col] >= tmin) &
  #                                 (pollen_ts[, age_col] <= tmax)),1:(taxa.start.col-1)]
  #   counts = pollen_ts[which((pollen_ts[, age_col] >= tmin) &
  #                              (pollen_ts[, age_col] <= tmax)),taxa.start.col:ncol(pollen_ts)]
  #   # meta_pol  = pollen_ts[which(pollen_ts[, age_col] <= tmax),1:(taxa.start.col-1)]
  #   # counts = pollen_ts[which(pollen_ts[, age_col] <= tmax),taxa.start.col:ncol(pollen_ts)]
  #   
  #   meta_pol = data.frame(meta_pol, age=rep(NA, nrow(meta_pol)))
  #   
  #   counts_agg = matrix(NA, nrow=0, ncol=ncol(counts))
  #   colnames(counts_agg) = colnames(counts)
  #   
  #   meta_agg = matrix(NA, nrow=0, ncol=6)
  #   colnames(meta_agg) = colnames(meta_pol)[1:6]
  #   
  #   meta_all = matrix(NA, nrow=0, ncol=ncol(meta_pol))
  #   colnames(meta_all) = colnames(meta_pol)
  #   
  #   ids = unique(meta_pol$id)
  #   ncores = length(ids)
  #   
  #   for (i in 1:ncores){
  #     
  #     print(i)
  #     
  #     #print(i)
  #     core_rows = which(meta_pol$id == ids[i])
  #     #     core_counts = counts[core_rows,]
  #     
  #     for (j in 1:(length(breaks)-1)){
  #       
  #       #print(j)
  #       age = breaks[j] + int/2
  #       
  #       age_rows = core_rows[(meta_pol[core_rows, age_col] >= breaks[j]) & 
  #                              (meta_pol[core_rows, age_col] < breaks[j+1])]
  #       
  #       if (length(age_rows)>1){
  #         
  #         counts_agg = rbind(counts_agg, colSums(counts[age_rows, ]))
  #         
  #         meta_agg      = rbind(meta_agg, data.frame(meta_pol[age_rows[1],1:6], age=age/100, zero=FALSE))
  #         
  #         meta_all_row = data.frame(meta_pol[age_rows,])
  #         meta_all_row$age = rep(age/100)
  #         meta_all = rbind(meta_all, meta_all_row)
  #         
  #       } else if (length(age_rows) == 1){
  #         
  #         counts_agg = rbind(counts_agg, counts[age_rows, ])
  #         
  #         meta_agg      = rbind(meta_agg, data.frame(meta_pol[age_rows,1:6], age=rep(age/100), zero=FALSE))
  #         # meta_agg$age[nrow(meta_agg)] = age/100
  #         
  #         meta_all_row = data.frame(meta_pol[age_rows,])
  #         meta_all_row$age = rep(age/100)
  #         meta_all = rbind(meta_all, meta_all_row)
  #         
  #       } else if (length(age_rows) == 0){
  #         
  #         #FIX ME
  #         counts_agg = rbind(counts_agg, rep(0,ncol(counts_agg)))
  #         
  #         meta_row = meta_pol[core_rows[1],1:6]
  #         meta_agg      = rbind(meta_agg, data.frame(meta_row, age=age/100, zero=TRUE))
  #         # meta_agg$age[nrow(meta_agg)] = age/100
  #         
  #         meta_all_row = data.frame(meta_pol[core_rows[1],])
  #         meta_all_row$age = NA
  #         meta_all_row$age_bacon = NA
  #         meta_all_row$age_default = NA
  #         
  #         meta_all = rbind(meta_all, meta_all_row)
  #         
  #       }
  #       
  #     }
  #   }
   }

#   counts = counts_agg
#   meta_pol = meta_agg
  return(list(counts_agg, meta_agg, meta_all)) 
}


# build pollen counts
build_pollen_counts_fast_core <- function(tmin, tmax, int, pollen_ts){
  
  idx = which(colnames(pollen_ts) == 'oak')
  
  if (int > 0){
    #   breaks = seq(0,2500,by=int)
    breaks = seq(tmin,tmax,by=int)
    
    meta_pol  = pollen_ts[which((pollen_ts[, 'ages'] >= tmin) & 
                                  (pollen_ts[, 'ages'] <= tmax)),1:(idx-1)]
    counts = pollen_ts[which((pollen_ts[, 'ages'] >= tmin) & 
                               (pollen_ts[, 'ages'] <= tmax)),idx:ncol(pollen_ts)]
    
    meta_agg = matrix(NA, nrow=0, ncol=ncol(meta_pol))
    colnames(meta_agg) = colnames(meta_pol)
    
    counts_agg = matrix(NA, nrow=0, ncol=ncol(counts))
    colnames(counts_agg) = colnames(counts)
    
    ids = unique(meta_pol$id)
    ncores = length(ids)
    
    for (j in 1:(length(breaks)-1)){
    
      for (i in 1:ncores){
      
        #print(i)
        core_rows = which(meta_pol$id == ids[i])
        #     core_counts = counts[core_rows,]
              
      #print(j)
        age = breaks[j] + int/2        
        age_rows = core_rows[(meta_pol[core_rows, 'ages'] >= breaks[j]) & 
                               (meta_pol[core_rows, 'ages'] < breaks[j+1])]
        
        if (length(age_rows)>1){
          
          counts_agg = rbind(counts_agg, colSums(counts[age_rows, ]))
          
          meta_agg      = rbind(meta_agg, meta_pol[core_rows[1],])
          meta_agg$ages[nrow(meta_agg)] = age/100
          
        } else if (length(age_rows) == 1){
          
          counts_agg = rbind(counts_agg, counts[age_rows, ])
          
          meta_agg      = rbind(meta_agg, meta_pol[core_rows[1],])
          meta_agg$ages[nrow(meta_agg)] = age/100
          
        } else if (length(age_rows) == 0){
          
          #FIX ME
          counts_agg = rbind(counts_agg, rep(0,ncol(counts_agg)))
          
          meta_agg      = rbind(meta_agg, meta_pol[core_rows[1],])
          meta_agg$ages[nrow(meta_agg)] = age/100
        }
        
      }
    }
  }
  
  #   counts = counts_agg
  #   meta_pol = meta_agg
  return(list(counts_agg, meta_agg)) 
}


# build idx_cores
build_idx_cores <- function(centers_polU, centers_pls, N_cores){

  idx_cores = vector(length=N_cores)

  for (i in 1:nrow(centers_polU)){
    core_site = centers_polU[i,]
    d1 = rdist(matrix(core_site, ncol=2), as.matrix(centers_pls))
    idx_cores[i] = which.min(d1)
  }
return(idx_cores)
}

# build the weight matrix
build_weight_matrix <- function(post, d, idx_cores, N, N_cores, run){
  
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  
  KW     = FALSE
  kernel = run$kernel
  phi    = post[which(par_names == 'phi')][1:K]
  
  if (kernel=='gaussian'){
    one_psi = run$one_psi
    if (one_psi){
      #       psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
      psi = post[which(par_names == 'psi')]
    } else {
      KW  = TRUE
      psi = post[which(par_names == 'psi')]
    }
  } else if (kernel=='pl'){
    one_a = run$one_a
    if (one_a){
      #       a = rep(mean(post[,1,which(par_names == 'a')]), K)
      a = post[which(par_names == 'a')]
    } else {
      KW = TRUE
      a  = post[which(par_names == 'a')]
    }
    
    one_b = run$one_b
    if (one_b){
      #       b = rep(mean(post[,1,which(par_names == 'b')]), K)
      b = post[which(par_names == 'b')]
      if (KW) b = rep(b, K)
    } else {
      KW = TRUE
      b  = post[which(par_names == 'b')]
    }
  }
  
  if (KW){
    w = array(0, c(K, N_cores, N))
    for (k in 1:K){
      print(paste0('k = ', k))
      for (i in 1:N_cores){
        for (j in 1:N){ 
          if (j != idx_cores[i]){
            if (kernel == 'gaussian'){
              w[k,i,j] = exp(-(d[i,j]*d[i,j])/(psi[k]*psi[k]))
            } else if (kernel == 'pl'){
              w[k,i,j] = (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) * (1 + d[i,j] / a[k]) ^ (-b[k])
            }
          } 
        }
      }
    }
  } else {
    w = array(0, c(N_cores, N))
      for (i in 1:N_cores){
        for (j in 1:N){ 
#           if (j != idx_cores[i]){
            if (kernel == 'gaussian'){
              w[i,j] = exp(-(d[i,j]*d[i,j])/(psi*psi))
            } else if (kernel == 'pl'){
              w[i,j] = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d[i,j] / a) ^ (-b)
            }
#           } 
        }
    }
  }
  
  
  #   for (i in 1:N_cores){
  #     for (j in 1:N){
  #       print(paste0("i = ", i))
  #       print(j)
  #       if ( d[idx_cores[i],j] > 0 ) {
  #         w[i,j] <- exp(-(d[idx_cores[i],j]/psi)^2)
  #       } 
  #     }
  #   }
  
  
  return(w)
}

# 
# # build the weight matrix
# build_weight_matrix <- function(post, d, idx_cores, N, N_cores, run){
#   
#   KW     = FALSE
#   kernel = run$kernel
#   phi    = colMeans(post[,1,which(par_names == 'phi')])[1:K]
#   
#   if (kernel=='gaussian'){
#     one_psi = run$one_psi
#     if (one_psi){
#       #       psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
#       psi = mean(post[,1,which(par_names == 'psi')])
#     } else {
#       KW  = TRUE
#       psi = colMeans(post[,1,which(par_names == 'psi')])
#     }
#   } else if (kernel=='pl'){
#     one_a = run$one_a
#     if (one_a){
#       #       a = rep(mean(post[,1,which(par_names == 'a')]), K)
#       a = mean(post[,1,which(par_names == 'a')])
#     } else {
#       KW = TRUE
#       a  = colMeans(post[,1,which(par_names == 'a')])
#     }
#     
#     one_b = run$one_b
#     if (one_b){
#       #       b = rep(mean(post[,1,which(par_names == 'b')]), K)
#       b = mean(post[,1,which(par_names == 'b')])
#       if (KW) b = rep(b, K)
#     } else {
#       KW = TRUE
#       b  = colMeans(post[,1,which(par_names == 'b')])
#     }
#   }
#   
#   if (KW){
#     w = array(0, c(K, N_cores, N))
#     for (k in 1:K){
#       print(paste0('k = ', k))
#       for (i in 1:N_cores){
#         for (j in 1:N){ 
#           if (j != idx_cores[i]){
#             if (kernel == 'gaussian'){
#               w[k,i,j] = exp(-(d[i,j]*d[i,j])/(psi[k]*psi[k]))
#             } else if (kernel == 'pl'){
#               w[k,i,j] = (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) * (1 + d[i,j] / a[k]) ^ (-b[k])
#             }
#           } 
#         }
#       }
#     }
#   } else {
#     w = array(0, c(N_cores, N))
#     for (i in 1:N_cores){
#       for (j in 1:N){ 
#         #           if (j != idx_cores[i]){
#         if (kernel == 'gaussian'){
#           w[i,j] = exp(-(d[i,j]*d[i,j])/(psi*psi))
#         } else if (kernel == 'pl'){
#           w[i,j] = (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d[i,j] / a) ^ (-b)
#         }
#         #           } 
#       }
#     }
#   }
#   
#   
#   #   for (i in 1:N_cores){
#   #     for (j in 1:N){
#   #       print(paste0("i = ", i))
#   #       print(j)
#   #       if ( d[idx_cores[i],j] > 0 ) {
#   #         w[i,j] <- exp(-(d[idx_cores[i],j]/psi)^2)
#   #       } 
#   #     }
#   #   }
#   
#   
#   return(w)
# }


# build the total potential neighborhood weighting
build_sumw_pot <- function(post, K, N_pot, d_pot, run){
  
  d_pot = d_pot[!(d_pot[,1] < 1e-8), ]
  
  KW = FALSE
  kernel = run$kernel

  par_names  = unlist(lapply(names(post), function(x) strsplit(x, "\\[")[[1]][1]))
  
  if (kernel=='gaussian'){
    one_psi = run$one_psi
    if (one_psi){
#       psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
      psi = post[which(par_names == 'psi')]
    } else {
      KW  = TRUE
      psi = post[which(par_names == 'psi')]
    }
    
    if (KW){
      sum_w = rep(NA, K)
      for (k in 1:K)
        sum_w[k] = sum(d_pot[,2] * exp(-d_pot[,1]^2/psi[k]^2))
    } else {
      sum_w = sum(d_pot[,2] * exp(-d_pot[,1]^2/psi^2))
    }
    
  } else if (kernel=='pl'){
    one_a = run$one_a
    one_b = run$one_b
    if (one_a){
#       a = rep(mean(post[,1,which(par_names == 'a')]), K)
      a = post[which(par_names == 'a')]
    } else {
      KW = TRUE
      a = post[which(par_names == 'a')]
    }
    if (one_b){
#       b = rep(mean(post[,1,which(par_names == 'b')]), K)
      b = post[which(par_names == 'b')]
      if (KW) b = rep(b, K)
    } else {
      KW = TRUE
      b  = post[which(par_names == 'b')]
    }

    if (KW){
      sum_w = rep(NA, K)
      for (k in 1:K)
        sum_w[k] = sum( d_pot[,2] * (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) * (1 + d_pot[,1] / a[k])^(-b[k]) )
    } else {
      sum_w = sum( d_pot[,2] * (b-1) * (b-2) / (2 * pi * a  * a) * (1 + d_pot[,1] / a)^(-b) )
    }
  }
  return(sum_w)
}

# # build the total potential neighborhood weighting
# build_sum_w_pot <- function(psi, rescale){
#   
#   x_pot = seq(-528000, 528000, by=8000)
#   y_pot = seq(-416000, 416000, by=8000)
#   coord_pot = expand.grid(x_pot, y_pot)
#   
#   d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)
#   d_pot = unname(as.matrix(count(d_pot)))
#   
#   N_pot = nrow(d_pot)
#   
#   sum_w_pot = 0
#   for (v in 1:N_pot){
#     sum_w_pot = sum_w_pot +  d_pot[v,2] *  exp(-(d_pot[v,1]/psi)^2)
#   }
#   return(sum_w_pot)
# }

# rescale gamma when the coarse grid is used
recompute_gamma <- function(w_coarse, sum_w_pot, gamma){
  
  core = 21

  gamma_new = gamma + ( 1 - gamma) * w_coarse / sum_w_pot
  
  return(gamma_new)
  
}

# determine which knots are in the domain given a list of cutlines
# cutlines are a list of two twoples: a point and a direction vector
choppy <- function(x, y, cutlines) {
  
  right = rep(TRUE, length(x))
  
  for (cl in cutlines) {
    x0 = cl[[1]][1]
    y0 = cl[[1]][2]
    dx = cl[[2]][1]
    dy = cl[[2]][2]
    r1 = sqrt(dx*dx+dy*dy)
    for (i in 1:length(x)) {
      rx = x[i] - x0
      ry = y[i] - y0
      r2 = sqrt(rx*rx+ry*ry)
      theta = acos((rx*dy - ry*dx)/r1/r2)
      #print(c(x0,y0,dx,dy,rx,ry,theta,r1,r2))
      if (theta > pi/2) {
        # point is to the left of the line
        right[i] = FALSE
      }
    }
  }
  which(right == FALSE)
}

#make a subgrid of cells
regular_subgrid <- function(cells, dx, dy, xoff, yoff){
  
  xlo = min(cells[,1])
  xhi = max(cells[,1]) 
  ylo = min(cells[,2])
  yhi = max(cells[,2])
  
  knots = matrix(nrow=0,ncol=2)
  colnames(knots) = c("x", "y")
  
  Nx = floor((xhi - xlo) / dx)
  Ny = floor((yhi - ylo) / dy)
  for (i in 0:Nx) {
    x = xlo + (i+xoff)*dx
    #print(x)
    
    for (j in 0:Ny) {
      y = ylo + (j+yoff)*dy
      #print(y)
      
      knots = rbind(knots, c(x,y))           
    }
  }
  return(knots)
}

get_knots <-function(knot_vals, cells, cell_width, thresh, count){
   
  cell_width = cell_width#diff(unique(cells[,1]))
  
  d = rdist(as.matrix(knot_vals), as.matrix(cells))
  
  knots_int = matrix(nrow=0, ncol=2)
  
#   thresh = 2.1*cell_width
  #thresh = sqrt(2*8^2)
  
  for (i in 1:nrow(knot_vals)){
    
    close = cells[which(d[i,] < thresh),]
    #print(close)
    
    N_nbrs = length(which(d[i,] < thresh))
    
    #print(N_nbrs)
#     if (N_neighbors>=2){
#     par(ask=TRUE)
#     plot(cells[,1],cells[,2])
#     points(knot_vals[i,1], knot_vals[i,2], col='blue',pch=19)
#     points(close[,1], close[,2], pch=19, col='red')
#     }
#     
    if (N_nbrs >= count) {
        
          #print(N_nbrs)
          #if (N_nbrs == 10){
          knots_int = rbind(knots_int, knot_vals[i,])
    }
  }
#     knot_x = knot_vals[i,1]
#     knot_y = knot_vals[i,2]
    

#     
#     x_close_idx = which.min(abs(knot_x-cells[,1]))
#     y_vals  = cells[which(cells[,1] == cells[x_close_idx,1]),2]
#     
#     y_close_idx = which.min(abs(knot_vals[i,2]-cells[,2]))
#     x_vals  = unique(cells[which(cells[,2] == cells[y_close_idx,2]),2])
#     
#     
#     
#     if ((knot_vals[i,2] <= max(y_vals)) & (knot_vals[i,2] >= min(y_vals))){ #&
#         #(knots[i,1] <= max(x_vals)) & (knots[i,1] >= min(x_vals))) {
#             knots_int = rbind(knots_int, knot_vals[i,])
#           }
    
    
    # WORKS!
    
#     knot_x = knot_vals[i,1]
#     knot_y = knot_vals[i,2]
#     
#     xright = knot_x + cell_width
#     xleft  = knot_x - cell_width
#     yup    = knot_y + cell_width
#     ydown  = knot_y - cell_width
#     
#     right = matrix(cbind(xright, knot_y), nrow=1)
#     left = matrix(cbind(xleft, knot_y), nrow=1)
#     up = matrix(cbind(knot_x, yup), nrow=1)
#     down = matrix(cbind(knot_x, ydown), nrow=1)
#     
#     d = rdist(as.matrix(cells), right)
#     nright = length(which(d < cell_width))
#     
#     d = rdist(as.matrix(cells), left)
#     nleft = length(which(d < cell_width))
#     
#     d = rdist(as.matrix(cells), up)
#     nup = length(which(d < cell_width))
#     
#     d = rdist(as.matrix(cells), down)
#     ndown = length(which(d < cell_width))
#     
#     if ((nup > 0) & (ndown > 0) & (nright > 0) & (nleft > 0)) {
#       knots_int = rbind(knots_int, knot_vals[i,])
#     }
  
  knots_int
}


build_domain_objects <- function(centers_pls, dx, cell_width, nclust){
  
  d = rdist(as.matrix(centers_pls))
  diag(d) <- 0
  
  knot_coords = kmeans(centers_pls, nclust, iter.max = 100, algorithm= "Hartigan-Wong")$centers
  
#   subgrid = regular_subgrid(centers_pls, dx=dx, dy=dx)
#   knot_coords = get_knots(subgrid, centers_pls, cell_width=cell_width)
  N_knots = dim(knot_coords)[1]

  d_knots = rdist(knot_coords, knot_coords)
  diag(d_knots) <- 0

  d_inter = rdist(centers_pls, knot_coords)
  d_inter[which(d_inter<1e-8)]=0
  
  return(list(d=d, d_knots=d_knots, d_inter=d_inter, knot_coords=knot_coords))
}

build_alpha_init <- function(W, N_knots, T, rho, tau, eta, d_knots, lag){
  
  alpha_init <- matrix(0, nrow=W, ncol=N_knots*T)
  
  for (k in 1:W){  
    
    C_s <- exp(-d_knots/rho[k]) # construct spatial covariance matrix
    C_t <- exp(-lag/tau)        # construct temporal covariance matrix
#     
#     C_s_inv <- solve(C_s)
#     C_t_inv <- solve(C_t)
    
    C_star     = eta[k]*kronecker(C_s, C_t)
    alpha_init[k,] = rmvnorm(n=1, mean=rep(0,N_knots*T), sigma=C_star)
    
#     C_star_inv = 1/eta[k]*kronecker(C_s_inv, C_t_inv)
#     alpha_init[k,] = rMVNorm(n=1, mean=rep(0,N_knots*T), Q=C_star_inv)
    
  }
  
  return(alpha_init)
}

pred_build_inits <- function(K, N, N_knots, eta, rho, mu, tau, d_knots, d_inter, lag){

  bt = TRUE
  
  if (length(rho) == K){
    bt = FALSE
    W = K
  } else {
    W = K-1
  }
  
  alpha_init <- matrix(0, nrow=W, ncol=N_knots*T)
  g_init     <- matrix(0, nrow=W, ncol=N*T)
  
  #   if (bt){
  #     alpha_init <- matrix(0, nrow=K, ncol=N_knots)
  #     g_init     <- matrix(0, nrow=K, ncol=N)
  #     nfit = K
  #   } else {
  #     alpha_init <- matrix(0, nrow=K-1, ncol=N_knots)
  #     g_init     <- matrix(0, nrow=K-1, ncol=N)
  #     nfit = W
  #   }
  
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  if (length(rho) < (K-1)){
    rho = rep(rho, W)
  }
  
  C_t <- exp(-lag/tau)        # construct temporal covariance matrix
  
  for (k in 1:W){  
    print(k)
    
    C_s = exp(-d_knots/rho[k]) # construct spatial covariance matrix
    c_s = exp(-d_inter/rho[k])
    
    C_s_inv = chol2inv(chol(C_s))
    C_star  = eta[k]^2*kronecker(C_s, C_t)
    
    alpha_init[k,] = rmvnorm(n=1, mean=rep(0,N_knots*T), sigma=C_star)

    H    = kronecker(c_s %*% C_s_inv, diag(T))
    cvar = c_s %*% C_s_inv %*% t(c_s)

    g_init[k,] = rnorm(mu[k] + H %*% alpha_init[k,], sd=sqrt(eta[k]^2 + eta[k]^2 * diag(cvar)))
    
  }
  
  return(list(alpha_init=alpha_init, g_init=g_init))
}

pred_build_inits_mut2 <- function(K, N, N_knots, eta, rho, mu_t, tau, d_knots, d_inter, lag){
  
  bt = TRUE
  
  if (length(rho) == K){
    bt = FALSE
    W = K
  } else {
    W = K-1
  }
  
  alpha_init <- matrix(0, nrow=W, ncol=N_knots*T)
  g_init     <- matrix(0, nrow=W, ncol=N*T)
  
  #   if (bt){
  #     alpha_init <- matrix(0, nrow=K, ncol=N_knots)
  #     g_init     <- matrix(0, nrow=K, ncol=N)
  #     nfit = K
  #   } else {
  #     alpha_init <- matrix(0, nrow=K-1, ncol=N_knots)
  #     g_init     <- matrix(0, nrow=K-1, ncol=N)
  #     nfit = W
  #   }
  
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  if (length(rho) < (K-1)){
    rho = rep(rho, W)
  }
  
  C_t <- exp(-lag/tau)        # construct temporal covariance matrix
  
  for (k in 1:W){  
    print(k)
    
    C_s = exp(-d_knots/rho[k]) # construct spatial covariance matrix
    c_s = exp(-d_inter/rho[k])
    
    C_s_inv = chol2inv(chol(C_s))
    C_star  = eta[k]^2*kronecker(C_s, C_t)
    
    alpha_init[k,] = rmvnorm(n=1, mean=rep(0,N_knots*T), sigma=C_star)
    
    H       = kronecker(c_s %*% C_s_inv, diag(T))
    cvar    = eta[k]^2 + eta[k]^2 * diag(c_s %*% C_s_inv %*% t(c_s))
    H_alpha = H %*% alpha_init[k,]
    
    for (i in 1:N){
      for (t in 1:T){
        g_init[k,(i-1)*T+t] = rnorm(1, mu_t[k, t] + H_alpha[(i-1)*T+t], sd=sqrt(cvar[i]))
      }
    }
    
  }
  
  return(list(alpha_init=alpha_init, g_init=g_init))
}

pred_build_inits_full <- function(K, N, N_knots, eta, rho, mu, mu_t, tau=0.1, d_knots, d_inter, lag){
  
  bt = TRUE
  
  if (length(rho) == K){
    bt = FALSE
    W = K
  } else {
    W = K-1
  }
  
  alpha_s_init <- matrix(0, nrow=W, ncol=N_knots)
  alpha_t_init <- matrix(0, nrow=W*(T-1), ncol=N_knots)
  g_init       <- matrix(0, nrow=W, ncol=N*T)
  
  Halpha_s <- matrix(0, nrow=N, ncol=1)
  Halpha_t <- array(0, c(K, N, T-1))
  
  #   if (bt){
  #     alpha_init <- matrix(0, nrow=K, ncol=N_knots)
  #     g_init     <- matrix(0, nrow=K, ncol=N)
  #     nfit = K
  #   } else {
  #     alpha_init <- matrix(0, nrow=K-1, ncol=N_knots)
  #     g_init     <- matrix(0, nrow=K-1, ncol=N)
  #     nfit = W
  #   }
  
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  if (length(rho) < (K-1)){
    rho = rep(rho, W)
  }
  
  for (k in 1:W){  
    print(paste0("k=", k))
    
    C_s = exp(-d_knots/rho[k]) # construct spatial covariance matrix
    c_s = exp(-d_inter/rho[k])
    
    C_s_inv = chol2inv(chol(C_s))
    
    Q_s = exp(-d_knots/lambda[k]) # construct spatial covariance matrix
    q_s = exp(-d_inter/lambda[k])
    
    Q_s_inv = chol2inv(chol(Q_s))
    
    alpha_s_init[k,] = rmvnorm(n=1, mean=rep(0,N_knots), sigma=eta[k] * eta[k] * C_s)
    
    alpha_t_init[(k-1)*(T-1)+1,] = rmvnorm(n=1, mean=rep(0,N_knots), sigma=sigma[k] * sigma[k] * Q_s)
    for (t in 2:(T-1))
      alpha_t_init[(k-1)*(T-1)+t,] = rmvnorm(n=1, mean=alpha_t_init[(k-1)*(T-1)+t-1,], sigma=sigma[k] * sigma[k] * Q_s)
     
    Halpha_s = c_s %*% C_s_inv %*% matrix(alpha_s_init[k,])
    
    for (t in 1:(T-1)){
      Halpha_t[k, , t] = q_s %*% Q_s_inv %*% matrix(alpha_t_init[(k-1)*(T-1)+t,])  
    }
    
    cvar    = eta[k]^2 - eta[k]^2 * diag(c_s %*% C_s_inv %*% t(c_s))
    qvar    = sigma[k]^2 - sigma[k]^2 * diag(q_s %*% Q_s_inv %*% t(q_s))
    
    cvar[abs(cvar) < 1e-8] = 0
    qvar[abs(qvar) < 1e-8] = 0
    
    for (i in 1:N){
      print(i)
      sqrt_var = sqrt(cvar[i] + tau)
      print(sqrt_var)
      
      
      if (sqrt_var > 0){
        # g_init[k,(i-1)*T+1] = rnorm(1, mu[k] + mu_t[k, 1] + Halpha_s[i,1], sd=sqrt_var)
        g_init[k,(i-1)*T+1] = rnorm(1, mu[k] + Halpha_s[i,1], sd=sqrt_var)
      } else if (sqrt_var == 0) {
        # g_init[k,(i-1)*T+1] = mu[k] + mu_t[k, 1] + Halpha_s[i,1]
        g_init[k,(i-1)*T+1] = mu[k] + Halpha_s[i,1]
      }
        
      sqrt_var = sqrt(cvar[i] + qvar[i] + tau)
      for (t in 2:T){
        if (sqrt_var > 0){
          g_init[k,(i-1)*T+t] = rnorm(1, mu[k] + mu_t[k, t-1] + Halpha_s[i,1] + Halpha_t[k,i,t-1], sd=sqrt_var)
        } else if (sqrt_var == 0){
          g_init[k,(i-1)*T+t] = mu[k] + mu_t[k, t] + Halpha_s[i,1] + Halpha_t[k,i,t-1]
        }
      }
    }
    
  }
  
  return(list(alpha_s_init=alpha_s_init, alpha_t_init=alpha_t_init, g_init=g_init))
}


pred_build_inits_ar <- function(K, N, N_knots, eta, rho, omega, mu, mu_t, tau=0.1, d_knots, d_inter, lag){
  
  bt = TRUE
  
  if (length(rho) == K){
    bt = FALSE
    W = K
  } else {
    W = K-1
  }
  
  alpha_s_init <- matrix(0, nrow=W, ncol=N_knots)
  alpha_t_init <- matrix(0, nrow=W*T, ncol=N_knots)
  g_init       <- matrix(0, nrow=W, ncol=N*T)
  
  Halpha_s <- matrix(0, nrow=N, ncol=1)
  Halpha_t <- array(0, c(K, N, T))
  
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  if (length(rho) < (K-1)){
    rho = rep(rho, W)
  }
  
  for (k in 1:W){  
    print(paste0("k=", k))
    
    C_s = exp(-d_knots/rho[k]) # construct spatial covariance matrix
    c_s = exp(-d_inter/rho[k])
    
    C_s_inv = chol2inv(chol(C_s))
    
    alpha_s_init[k,] = rmvnorm(n=1, mean=rep(0,N_knots), sigma=eta[k] * eta[k] * C_s)
    
    alpha_t_init[(k-1)*T+1,] = rmvnorm(n=1, mean=rep(0,N_knots), sigma=C_s)
    for (t in 2:T)
      alpha_t_init[(k-1)*T+t,] = rmvnorm(n=1, mean=omega[k] * alpha_t_init[(k-1)*T+t-1,], sigma=(1 - omega[k] * omega[k]) * C_s)
    
    Halpha_s = c_s %*% C_s_inv %*% matrix(alpha_s_init[k,])
    
    for (t in 1:T){
      Halpha_t[k, , t] = c_s %*% C_s_inv %*% matrix(alpha_t_init[(k-1)*T+t,])  
    }
    
    allvar = diag(c_s %*% C_s_inv %*% t(c_s))
    cvar    = eta[k]^2 - eta[k]^2 * allvar
    qvar    = (1 - omega[k] * omega[k]) - (1 - omega[k] * omega[k]) * allvar
    
    cvar[abs(cvar) < 1e-8] = 0
    qvar[abs(qvar) < 1e-8] = 0
    
    for (i in 1:N){
      sqrt_var = sqrt(cvar[i] + tau)
      
      if (sqrt_var > 0){
        # g_init[k,(i-1)*T+1] = rnorm(1, mu[k] + mu_t[k, 1] + Halpha_s[i,1], sd=sqrt_var)
        g_init[k,(i-1)*T+1] = rnorm(1, mu[k] + Halpha_s[i,1] + Halpha_t[k,i,t], sd=sqrt_var)
      } else if (sqrt_var == 0) {
        # g_init[k,(i-1)*T+1] = mu[k] + mu_t[k, 1] + Halpha_s[i,1]
        g_init[k,(i-1)*T+1] = mu[k] + Halpha_s[i,1]
      }
      
      sqrt_var = sqrt(cvar[i] + qvar[i] + tau)
      for (t in 2:T){
        if (sqrt_var > 0){
          g_init[k,(i-1)*T+t] = rnorm(1, mu[k] + mu_t[k, t-1] + Halpha_s[i,1] + Halpha_t[k,i,t], sd=sqrt_var)
        } else if (sqrt_var == 0){
          g_init[k,(i-1)*T+t] = mu[k] + mu_t[k, t] + Halpha_s[i,1] + Halpha_t[k,i,t]
        }
      }
    }
    
  }
  
  return(list(alpha_s_init=alpha_s_init, alpha_t_init=alpha_t_init, g_init=g_init))
}


pred_build_inits_space <- function(K, N, N_knots, eta, rho, mu, mu_t, tau=0.1, d_knots, d_inter, lag, 
                                   alpha_s_init, alpha_t_init){
  
  bt = TRUE
  
  if (length(rho) == K){
    bt = FALSE
    W = K
  } else {
    W = K-1
  }
  
  g_init   <- matrix(0, nrow=W, ncol=N*T)
  
  Halpha_s <- matrix(0, nrow=N, ncol=1)
  Halpha_t <- array(0, c(K, N, T-1))
  
  ones     = matrix(1, nrow=N, ncol=1)
  temp = qr(ones)
  Q    = qr.Q(temp)
  P = Q %*% t(Q)
  
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  if (length(rho) < (K-1)){
    rho = rep(rho, W)
  }
  
  for (k in 1:W){  
    print(paste0("k=", k))
    
    C_s = exp(-d_knots/rho[k]) # construct spatial covariance matrix
    c_s = exp(-d_inter/rho[k])
    
    C_s_inv = chol2inv(chol(C_s))
    
    Q_s = exp(-d_knots/lambda[k]) # construct innovation covariance matrix
    q_s = exp(-d_inter/lambda[k])
    
    Q_s_inv = chol2inv(chol(Q_s))
    
    Halpha_s = (c_s %*% C_s_inv - P %*% c_s %*% C_s_inv) %*% matrix(alpha_s_init[k,])
    
    for (t in 1:(T-1)){
      Halpha_t[k, , t] = q_s %*% Q_s_inv %*% matrix(alpha_t_init[(k-1)*(T-1)+t,])  
    }
    
    cvar    = eta[k]^2 - eta[k]^2 * diag(c_s %*% C_s_inv %*% t(c_s))
    qvar    = sigma[k]^2 - sigma[k]^2 * diag(q_s %*% Q_s_inv %*% t(q_s))
    
    cvar[abs(cvar) < 1e-8] = 0
    qvar[abs(qvar) < 1e-8] = 0
    
    for (i in 1:N){
      print(i)
      
      sqrt_var = sqrt(cvar[i] + tau)
      g_init[k,(i-1)*T+1] = rnorm(1, mu[k] + Halpha_s[i,1], sd=sqrt_var)

      sqrt_var = sqrt(cvar[i] + qvar[i] + tau)
      for (t in 2:T){
        g_init[k,(i-1)*T+t] = rnorm(1, mu[k] + mu_t[k, t-1] + Halpha_s[i,1] + Halpha_t[k,i,t-1], sd=sqrt_var)
      }
    }
    
  }
  
  return(g_init)
}

space_g_init_ar <- function(K, N, N_knots, eta, rho, omega, mu, mu_t, tau=0.1, d_knots, d_inter,
                            alpha_s_init, alpha_t_init){
  
  W = K-1
  
  g_init   <- matrix(0, nrow=W, ncol=N*T)
  Halpha_s <- matrix(0, nrow=N, ncol=1)
  Halpha_t <- array(0, c(K, N, T))
  
  ones     = matrix(1, nrow=N, ncol=1)
  temp = qr(ones)
  Q    = qr.Q(temp)
  P = Q %*% t(Q)
  
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  if (length(rho) < (K-1)){
    rho = rep(rho, W)
  }
  
  for (k in 1:W){  
    print(paste0("k=", k))
    
    C_s = exp(-d_knots/rho[k]) # construct spatial covariance matrix
    c_s = exp(-d_inter/rho[k])
    
    C_s_inv = chol2inv(chol(C_s))
    
    Halpha_s = (c_s %*% C_s_inv - P %*% c_s %*% C_s_inv) %*% matrix(alpha_s_init[k,])
    
    for (t in 1:T){
      Halpha_t[k, , t] = c_s %*% C_s_inv %*% matrix(alpha_t_init[(k-1)*T+t,])  
    }
    
    cvar    = eta[k]^2 - eta[k]^2 * diag(c_s %*% C_s_inv %*% t(c_s))
    qvar    = (1-omega[k]^2)*(1 - diag(c_s %*% C_s_inv %*% t(c_s)))
    
    cvar[abs(cvar) < 1e-8] = 0
    qvar[abs(qvar) < 1e-8] = 0
    
    for (i in 1:N){
      print(i)
      
      sqrt_var = sqrt(cvar[i] + tau)
      g_init[k,(i-1)*T+1] = rnorm(1, mu[k] + Halpha_s[i,1] + Halpha_t[k,i,1], sd=sqrt_var)
      
      sqrt_var = sqrt(cvar[i] + qvar[i] + tau)
      for (t in 2:T){
        g_init[k,(i-1)*T+t] = rnorm(1, mu[k] + mu_t[k, t-1] + Halpha_s[i,1] + Halpha_t[k,i,t], sd=sqrt_var)
      }
    }
    
  }
  
  return(g_init)
}

build_alpha_init_v2 <- function(W, N_knots, T, rho, tau, eta, d_knots, lag){
  
  alpha_init <- vector(length=N_knots*T)
  
  #for (k in 1:W){  
    
    C_s <- exp(-d_knots/rho) # construct spatial covariance matrix
    C_t <- exp(-lag/tau)        # construct temporal covariance matrix
    #     
    #     C_s_inv <- solve(C_s)
    #     C_t_inv <- solve(C_t)
    
    C_star     = eta*kronecker(C_s, C_t)
    alpha_init = as.vector(rmvnorm(n=1, mean=rep(0,N_knots*T), sigma=C_star))
    
    #     C_star_inv = 1/eta[k]*kronecker(C_s_inv, C_t_inv)
    #     alpha_init[k,] = rMVNorm(n=1, mean=rep(0,N_knots*T), Q=C_star_inv)
    
  #}
  
  return(alpha_init)
}


# compute effective sample size form a stanfit object
ess <- function(fit){
  ess = summary(fit)$summary[,"n_eff"]
  return(ess)
}

pollen_to_albers <- function(pollen_ts){

  centers_pol = data.frame(x=pollen_ts$long, y=pollen_ts$lat)

  coordinates(centers_pol) <- ~ x + y
  proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

  centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
  centers_polA <- as.matrix(data.frame(centers_polA))/1000000

  pollen_ts$long = centers_polA[,'x']
  pollen_ts$lat = centers_polA[,'y']
 
  colnames(pollen_ts)[grep("lat", colnames(pollen_ts))] = 'y'
  colnames(pollen_ts)[grep("long", colnames(pollen_ts))] = 'x'
  
  return(pollen_ts)
}


other_build <- function(y, counts, other.idx){
  
  if (sum(other.idx)>1){
    if (!is.null(y$other.hardwood))
      other.vec = y$other.hardwood + rowSums(counts[, other.idx])
    else 
      other.vec = rowSums(counts[, other.idx])
  } else if  (sum(other.idx)==1) {
    if (!is.null(y$other.hardwood))
      other.vec = y$other.hardwood + counts[, other.idx]
    else 
      other.vec = counts[, other.idx]
  } else {
    if (!is.null(y$other.hardwood))
      other.vec = y$other.hardwood 
    else 
      other.vec = rep(0, nrow(y))
  }
  
  return(other.vec)
}

y_build <- function(counts, taxa_sub){ 
  
  three.p = data.frame(read.table(file='data/level3P_v0.2.csv', sep=",", row.names=NULL, header=TRUE, stringsAsFactors=FALSE))
  three.p = rbind(three.p, c('Other conifer', 'Other conifer', 'TRUE'))
  three.p[,1:2] = as.data.frame(apply(three.p[,1:2],2,function(x)gsub('\\s+', '.',x)))
  three.p[,1:2] = as.data.frame(apply(three.p[,1:2],2,function(x)gsub('\\/', '.', x)))
  three.p = as.data.frame(apply(three.p,2,function(x)toupper(x)))
  
  taxa_veg = colnames(counts)
  taxa_use = taxa_sub
  taxa_3p = toupper(three.p[,2])
  
  con = as.logical(three.p[match(taxa_veg, taxa_3p),3])
  
  other.hw.idx = !(taxa_veg %in% taxa_use) & !con
  other.con.idx = !(taxa_veg %in% taxa_use) & con
  
  if (sum(taxa_veg %in% taxa_use) < length(taxa_sub)) print('One or more of the taxa in the provided list is not in 
                                                            the pls data or appears under a different name.')
  
  y       = counts[, taxa_veg %in% taxa_use]
  
  # other is super annoying
  
  y$OTHER.HARDWOOD = other_build(y, counts, other.hw.idx)
  y$OTHER.CONIFER  = other_build(y, counts, other.con.idx)
  
  taxa = colnames(y)
  y       = unname(round(as.matrix(y)))
  
  return(list(y=y, taxa=taxa))
}

knots_in_domain4 <-function(knots, cells, cell_width){
  
  knots_int = matrix(nrow=0, ncol=2)
  
  for (i in 1:nrow(knots)){
    x = knots[i,1]
    y = knots[i,2]
    
    xright = x + cell_width
    xleft  = x - cell_width
    yup    = y + cell_width
    ydown  = y - cell_width
    
    right = matrix(cbind(xright, y), nrow=1)
    left = matrix(cbind(xleft, y), nrow=1)
    up = matrix(cbind(x, yup), nrow=1)
    down = matrix(cbind(x, ydown), nrow=1)
    
    d = rdist(as.matrix(cells), right)
    nright = length(which(d < cell_width))
    
    d = rdist(as.matrix(cells), left)
    nleft = length(which(d < cell_width))
    
    d = rdist(as.matrix(cells), up)
    nup = length(which(d < cell_width))
    
    d = rdist(as.matrix(cells), down)
    ndown = length(which(d < cell_width))
    
    if ((nup > 0) & (ndown > 0) & (nright > 0) & (nleft > 0)) {
      knots_int = rbind(knots_int, knots[i,])
    }
    
  }
  
  knots_int
}

# cores_near_domain(centers_polA/rescale, centers_veg, cell_width = res*8000/rescale)

cores_near_domain <-function(knots, cells, cell_width){
    
  knots_int = matrix(nrow=0, ncol=2)
   
   for (i in 1:nrow(knots)){
#      x = knots[i,1]
#      y = knots[i,2]
     
     d = rdist(as.matrix(cells), matrix(t(knots[i,]), nrow=1))
    
#     xright = x + cell_width
#     xleft  = x - cell_width
#     yup    = y + cell_width
#     ydown  = y - cell_width
#     
#     right = matrix(cbind(xright, y), nrow=1)
#     left = matrix(cbind(xleft, y), nrow=1)
#     up = matrix(cbind(x, yup), nrow=1)
#     down = matrix(cbind(x, ydown), nrow=1)
#     
#     d = rdist(as.matrix(cells), right)
#     nright = length(which(d < cell_width))
#     
#     d = rdist(as.matrix(cells), left)
#     nleft = length(which(d < cell_width))
#     
#     d = rdist(as.matrix(cells), up)
#     nup = length(which(d < cell_width))
#     
#     d = rdist(as.matrix(cells), down)
#     ndown = length(which(d < cell_width))
#     
#     if ((nup > 0) | (ndown > 0) | (nright > 0) | (nleft > 0)) {
#       knots_int = rbind(knots_int, knots[i,])
#     }
    
    if (any(d < cell_width)) {
      knots_int = rbind(knots_int, knots[i,])
    }


   }
  
  knots_int
}

convert_counts <- function(counts, tree_type, taxa_sub){
  
  colnames(counts) = toupper(colnames(counts))
  
  y_veg = data.frame(counts[, colnames(counts) %in% taxa_sub])
  
  left = rownames(tree_type)[!(rownames(tree_type) %in% taxa_sub)]
  taxa_other_hw  = rownames(tree_type)[which(tree_type$type == 'HW')]
  taxa_other_con = rownames(tree_type)[which(tree_type$type == 'CON')]
  
  if (sum(left %in% taxa_other_hw) > 1){
    y_veg$OTHER.HARDWOOD = rowSums(counts[,left[left %in% taxa_other_hw]])
  } else {
    y_veg$OTHER.HARDWOOD = counts[,left[left %in% taxa_other_hw]]
  }
  if (sum(left %in% taxa_other_con) > 1){
    y_veg$OTHER.CONIFER  = rowSums(counts[,left[left %in% taxa_other_con]])
  } else {
    y_veg$OTHER.CONIFER  = counts[,left[left %in% taxa_other_con]]
  }
  
  y_veg = round(as.matrix(y_veg[,sort(colnames(y_veg))]))
  
  return(y_veg)
}



split_mi <- function(meta){
  
  centers = data.frame(x=meta$x, y=meta$y)
  
  coordinates(centers) <- ~ x + y
  proj4string(centers) <- CRS('+init=epsg:3175')
  
  centers_ll <- spTransform(centers, CRS('+proj=longlat +ellps=WGS84'))
  centers_ll <- as.matrix(data.frame(centers_ll))
  
  if (!any(colnames(meta) == 'state')) {
    meta$state = map.where(database='state', centers_ll[,1], centers_ll[,2])
  }
  
  idx.mi = which(meta$state=='michigan_north')
  meta$state2 = as.vector(meta$state)
  meta$state2[idx.mi] = map.where(database="state", centers_ll[idx.mi,1], centers_ll[idx.mi,2])
  idx.na = which(is.na(meta$state2))
  idx.not.na = which(!is.na(meta$state2))
  
  idx.mi.s = which(meta$state=='michigan_south')
  meta$state2[idx.mi.s] = 'michigan:south'#map.where(database="state", centersLL[idx.mi.s,1], centersLL[idx.mi.s,2])
  
  for (i in 1:length(idx.na)){
    idx = idx.na[i]
    centers = centers_ll[idx.not.na,]
    dmat = rdist(matrix(centers_ll[idx,], nrow=1) , matrix(centers, ncol=2))
    min.val = dmat[1,which.min(dmat[which(dmat>1e-10)])]
    idx_close = which(dmat == min.val)
    state  = map.where(database="state", centers[idx_close,1], centers[idx_close,2])
    meta$state2[idx] = state
  }
  
  meta$state2[which(meta$state2[idx.mi]=='minnesota')] = 'michigan:north'
  
  if (any(meta$state2 %in% c('illinois', 'indiana', 'ohio', 'iowa', 'north dakota', 'south dakota'))) {
    meta$state2[which(meta$state2 == 'illinois')] = 'wisconsin'
    meta$state2[which(meta$state2 %in% c('indiana', 'ohio'))] = 'michigan:south'
    meta$state2[which((meta$state2 == 'iowa') & (meta$x > 370000))] = 'wisconsin'
    meta$state2[which((meta$state2 == 'iowa') & (meta$x < 370000))] = 'minnesota'
    meta$state2[which(meta$state2 %in% c('north dakota', 'south dakota'))] = 'minnesota'
  }
  
  return(meta)
  
}

get_quants <- function(post, npars){
  
#   quants <- colMeans(post[,1,1:npars])
  
  quants <- t(apply(post[,1:npars], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
  
#   quants <- cbind(summary(fit)$summary[,'mean'][1:npars],
#                   summary(fit)$summary[,'2.5%'][1:(npars)],
#                   summary(fit)$summary[,'50%'][1:(npars)],
#                   summary(fit)$summary[,'97.5%'][1:(npars)],
#                   summary(fit)$summary[,'n_eff'][1:(npars)],
#                   summary(fit)$summary[,'Rhat'][1:(npars)])
#   colnames(quants) = c('mean', '2.5%', '50%', '97.5%', 'n_eff', 'Rhat')
  
  return(quants)
}


build_props <- function(post, rho, eta, tau, mu, alpha, N_knots, T, K, N_pars, mpp){
  
  #if (length(eta == 1)){eta = rep(eta, W)}
  
  #    niter = nrow(post)
  niter=10
  W     = K-1
  
  g = array(NA, dim=c(N*T, W, niter))
  r    = array(NA, dim=c(N*T, K, niter))
  ones = matrix(1, nrow=N*T, ncol=1)
  
  alpha = post[,1,(N_pars+1):(ncol(post[,1,])-1)]
  
  n = seq(1, N_knots)
  
  for (i in 1:niter){
    
    print(i)
    
    for (k in 1:W){  
      #       ta <- proc.time()
      #     
      #       C_s <- exp(-d_knots/rho[k]) # construct spatial covariance matrix
      #       C_t <- exp(-lag/tau)        # construct temporal covariance matrix
      #     
      # #       tb <- proc.time()
      #       
      #       C_s_inv = chol2inv(chol(C_s))
      #       C_t_inv = chol2inv(chol(C_t))
      #       
      # #       tc <- proc.time()
      #       
      #       C_star_inv = 1/(eta[k]*eta[k])*kronecker(C_s_inv, C_t_inv)
      #       
      # #       td <- proc.time()
      #       
      #       c <- build_c(rho[k], tau, eta[k], d_inter, lag)
      
      #       te <- proc.time()
      #       print('Build c:')
      #       print(te-td)
      
      #alpha_idx <- seq(k,ncol(alpha),by=T) 
      alpha_idx <- seq(k,ncol(alpha),by=W) 
      alpha_k   <- alpha[i,alpha_idx]
      #       knot_idx <- function(w, n, t){
      #         6 + (n-1)*T*W + (t-1)*W + w-1
      #       }
      #       alpha_k <- alpha[((W-1)*N_knots*T + 1):(W*N_knots*T)]
      
      #       tf <- proc.time()
      
      H_alpha <- c%*%C_star_inv%*%alpha_k    
      
      #       tg <- proc.time()
      #       print('Matrix mult:')
      #       print(tg-tf)
      
      g[,k,i] <- mu[k]*ones + H_alpha
      
      #       th <-proc.time()
      
    }
    
    sum_exp_g = rowSums(exp(g[,,i]))
    
    #     ti <- proc.time()
    
    # additive log-ratio transformation
    for (k in 1:W)
      for (j in 1:(N*T))
        r[j,k,i] <- exp(g[j,k,i]) / (1 + sum_exp_g[j])
    
    for (j in 1:(N*T))
      r[j,K,i] <- 1 / (1 + sum_exp_g[j])
    #     
    #     tj <- proc.time()
    #     print("OLD Build r:")
    #     print(tj-ti)
    #     
    #     r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
    #     tk <- proc.time()
    #     print("NEW Build r:")
    #     print(tk-tj)
  }
  
  return(list(r=r, g=g))
}


build_props_new <- function(post, rho, eta, T, K, d, d_inter, d_knots, mpp){
  
  N = nrow(d_inter)
  N_knots = ncol(d_inter)
  niter   = dim(post[,1,])[1] 
  
  W = K-1
  
  #   g = array(NA, c(W, N, niters))
  #   Halpha = array(NA, c(W, N, niters))
  
  #if (length(eta == 1)){eta = rep(eta, W)}
  
  #    niter = nrow(post)
  W     = K-1
  
  #   g = array(NA, c(W, N, niters))
  #   Halpha = array(NA, c(W, N, niters))
  
  g    = array(NA, dim=c(N*T, W, niter))
  r    = array(NA, dim=c(N*T, K, niter))
  Halpha = array(NA, dim=c(N_knots, W, niter))
  ones = matrix(1, nrow=N*T, ncol=1)
  
  #   alpha = post[,1,(N_pars+1):(ncol(post[,1,])-1)]
  
  #   n = seq(1, N_knots)
  
  col_substr = substr(colnames(post[,1,]),1,2)
  
  tau   = post[,1,which(col_substr == 'ta')]
  
  for (k in 1:W){
    
    if (mpp){
      
      mu    = post[,1,which(col_substr == 'mu')[k]]
      
      knot_cols = seq((W + 1 + k), (W + k + W*N_knots*T), by=W) 
      
      alpha     = post[,1, knot_cols]
      
      
      for (i in 1:niter){
        C_s <- exp(-d_knots/rho[k]) 
        C_t <- exp(-lag/tau[i])       
        
        C_s_inv = chol2inv(chol(C_s))
        C_t_inv = chol2inv(chol(C_t))
        
        #       tc <- proc.time()
        
        C_star_inv = 1/eta[k]*kronecker(C_s_inv, C_t_inv)
      }
      
      g_cols    = seq((W + W*N_knots*T + k + 1), (W + W*N_knots*T + k + W*N*T), by=W)
      
      g[,k,]         = t(post[,1,g_cols])
      
      #     for (i in 1:niters){
      #       
      #       
      #       
      #     }
    }
  }
  
  
  
  #   alpha     = post[,1, knot_cols]
  # 
  #   for (i in 1:niter){
  #     
  #     print(i)
  #     
  #     for (k in 1:W){  
  #       
  #       C_s <- exp(-d_knots/rho[k]) # construct spatial covariance matrix
  #       C_t <- exp(-lag/tau)        # construct temporal covariance matrix
  #       
  #       #       tb <- proc.time()
  #       
  #       C_s_inv = chol2inv(chol(C_s))
  #       C_t_inv = chol2inv(chol(C_t))
  #       
  #       #       tc <- proc.time()
  #       
  #       C_star_inv = 1/eta[k]*kronecker(C_s_inv, C_t_inv)
  #       
  #       #       td <- proc.time()
  #       
  #       c <- build_c(rho[k], tau, eta[k], d_inter, lag)
  #       
  #       #       te <- proc.time()
  #       #       print('Build c:')
  #       #       print(te-td)
  #       
  #       #alpha_idx <- seq(k,ncol(alpha),by=T) 
  #       alpha_idx <- seq(k,ncol(alpha),by=W) 
  #       alpha_k   <- alpha[i,alpha_idx]
  #       #       knot_idx <- function(w, n, t){
  #       #         6 + (n-1)*T*W + (t-1)*W + w-1
  #       #       }
  #       #       alpha_k <- alpha[((W-1)*N_knots*T + 1):(W*N_knots*T)]
  #       
  #       #       tf <- proc.time()
  #       
  #       H_alpha <- c%*%C_star_inv%*%alpha_k    
  #       
  #       #       tg <- proc.time()
  #       #       print('Matrix mult:')
  #       #       print(tg-tf)
  #       
  #       g[,k,i] <- mu[k]*ones + H_alpha
  #       
  #       #       th <-proc.time()
  #       
  #     }
  
  for (i in 1:niter){
    
    sum_exp_g = rowSums(exp(g[,,i]))
    
    #     ti <- proc.time()
    
    # additive log-ratio transformation
    for (k in 1:W)
      for (j in 1:(N*T))
        r[j,k,i] <- exp(g[j,k,i]) / (1 + sum_exp_g[j])
    
    for (j in 1:(N*T))
      r[j,K,i] <- 1 / (1 + sum_exp_g[j])
    #     
    #     tj <- proc.time()
    #     print("OLD Build r:")
    #     print(tj-ti)
    #     
    #     r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
    #     tk <- proc.time()
    #     print("NEW Build r:")
    #     print(tk-tj)
  }
  
  return(list(r=r, g=g))
}


build_props_mut <- function(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mut){
  
  N = nrow(d_inter)
  N_knots = ncol(d_inter)
  niter   = dim(post[,1,])[1] 
  
  W = K-1
  
  g    = array(NA, dim=c(N*T, W, niter))
  r    = array(NA, dim=c(N*T, K, niter))
  #   Halpha = array(NA, dim=c(N_knots, W, niter))
  ones = matrix(1, nrow=N*T, ncol=1)
  
  Halpha = array(NA, dim=c(N*T, W, niter))
  sumHalpha = array(NA, dim=c(T, W, niter))
  
  col.names = colnames(post[,1,])
  
  #   alpha = post[,1,(N_pars+1):(ncol(post[,1,])-1)]
  
  #   n = seq(1, N_knots)
  
  col_substr = substr(colnames(post[,1,]),1,2)
  #  
  #   tau   = post[,1,which(col_substr == 'ta')]
  #   ksi   = post[,1,which(col_substr == 'ks')]
  #   omega = post[,1,which(col_substr == 'om')]
  tau = post[,1,1]
  
  if (od){
    x = matrix(1, nrow=(N*T), ncol=1)
    N_p = N*T
    
    temp = qr(x)
    Q = qr.Q(temp)
    #     R = qr.R(temp)
    
    #     P = Q %*% t(Q)
  }
  
  for (k in 1:W){
    print(k)
    mu    = post[,1,which(col_substr == 'mu')[k]]
    #       mut_cols = seq((3+W+k), (3+W+k+W*T-1), by=W)
    #       colnames(post[,1,])[mut_cols]
    #       mu_t      = post[,1,mut_cols]
    if (od & !mpp & !mut){
      knot_cols = seq((W + 1 + k), (W + 1 + k + W*N_knots*T - 1), by=W) 
    } else if  (od & mpp & mut) {
      knot_cols = seq((W + 3 + k + W*T), (W + 3 + k + W*T + W*N_knots*T - 1), by=W) 
    }
    
    print(col.names[knot_cols][1:10])
    print(length(col.names[knot_cols]))
    alpha     = post[,1, knot_cols]
    
    C_s <- exp(-d_knots/rho[k])
    c_s <- exp(-d_inter/rho[k])
    C_s_inv = chol2inv(chol(C_s))
    
    cs_Csinv = c_s %*% C_s_inv
    
    c_Cstarinv = kronecker(cs_Csinv, diag(T))
    
    for (i in 1:niter){
      print(i)
      
      #         C_t <- exp(-lag/tau[i])       
      #         
      #         
      #         C_t_inv = chol2inv(chol(C_t))
      #         
      #         #       tc <- proc.time()
      #         
      #         C_star_inv = 1/eta[k]*kronecker(C_s_inv, C_t_inv)
      #         c <- build_c(rho[k], tau[i], eta[k], d_inter, lag)
      #         
      #         Halpha2 <- c %*% C_star_inv %*% alpha[i,]
      if (od){
        c_Cinv_alpha <- c_Cstarinv %*% alpha[i,]
        Halpha[,k,i] <- c_Cinv_alpha - Q %*% (t(Q) %*% c_Cinv_alpha)
      } else {
        Halpha[,k,i] <- c_Cstarinv %*% alpha[i,]
      }
      
      
      
      for (t in 1:T){
        sumHalpha[t, k, i] = sum(Halpha[seq(t, N*T, by=T),k,i])
      }
      
    }
    
    if (od & !mpp & !mut){
      g[,k,] = mu[k] + Halpha[,k,]
    } else if  (od & mpp & mut) {
      g_cols = seq((3 + W + W*T + W*N_knots*T + k), (3 + W + W*T + W*N_knots*T + k + W*N*T - 1), by=W)
      colnames(post[,1,])[g_cols]
      g[,k,]         = t(post[,1,g_cols])
    }
    
  }
  
  
  
  #   alpha     = post[,1, knot_cols]
  # 
  #   for (i in 1:niter){
  #     
  #     print(i)
  #     
  #     for (k in 1:W){  
  #       
  #       C_s <- exp(-d_knots/rho[k]) # construct spatial covariance matrix
  #       C_t <- exp(-lag/tau)        # construct temporal covariance matrix
  #       
  #       #       tb <- proc.time()
  #       
  #       C_s_inv = chol2inv(chol(C_s))
  #       C_t_inv = chol2inv(chol(C_t))
  #       
  #       #       tc <- proc.time()
  #       
  #       C_star_inv = 1/eta[k]*kronecker(C_s_inv, C_t_inv)
  #       
  #       #       td <- proc.time()
  #       
  #       c <- build_c(rho[k], tau, eta[k], d_inter, lag)
  #       
  #       #       te <- proc.time()
  #       #       print('Build c:')
  #       #       print(te-td)
  #       
  #       #alpha_idx <- seq(k,ncol(alpha),by=T) 
  #       alpha_idx <- seq(k,ncol(alpha),by=W) 
  #       alpha_k   <- alpha[i,alpha_idx]
  #       #       knot_idx <- function(w, n, t){
  #       #         6 + (n-1)*T*W + (t-1)*W + w-1
  #       #       }
  #       #       alpha_k <- alpha[((W-1)*N_knots*T + 1):(W*N_knots*T)]
  #       
  #       #       tf <- proc.time()
  #       
  #       H_alpha <- c%*%C_star_inv%*%alpha_k    
  #       
  #       #       tg <- proc.time()
  #       #       print('Matrix mult:')
  #       #       print(tg-tf)
  #       
  #       g[,k,i] <- mu[k]*ones + H_alpha
  #       
  #       #       th <-proc.time()
  #       
  #     }
  
  for (i in 1:niter){
    
    sum_exp_g = rowSums(exp(g[,,i]))
    
    #     ti <- proc.time()
    
    # additive log-ratio transformation
    for (k in 1:W)
      for (j in 1:(N*T))
        r[j,k,i] <- exp(g[j,k,i]) / (1 + sum_exp_g[j])
    
    for (j in 1:(N*T))
      r[j,K,i] <- 1 / (1 + sum_exp_g[j])
    #     
    #     tj <- proc.time()
    #     print("OLD Build r:")
    #     print(tj-ti)
    #     
    #     r[,,i] <- sum2one_constraint(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
    #     tk <- proc.time()
    #     print("NEW Build r:")
    #     print(tk-tj)
  }
  
  return(list(r=r, g=g, sumHalpha=sumHalpha, Halpha=Halpha))
}


# ess_rfun <- function(sims) {
#   # Compute the effective sample size for samples of several chains 
#   # for one parameter; see the C++ code of function  
#   # effective_sample_size2 in chains.cpp 
#   # 
#   # Args:
#   #   sims: a 2-d array _without_ warmup samples (# iter * # chains) 
#   # 
#   # Note: 
#   #   The implementation in R uses acf in R to compute the autocovariance
#   #   and the results might be a little bit different from that in stan. 
#   #   The R function wrapping the C++ implementation is defined in 
#   #   chains.R with name rstan_ess2_cpp 
#   if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
#   chains <- ncol(sims)
#   n_samples <- nrow(sims)
#   
#   acov <- lapply(1:chains, 
#                  FUN = function(i) {
#                    cov <- acf(sims[,i], lag.max = n_samples - 1, 
#                               plot = FALSE, type = c("covariance"), na.action=na.pass) 
#                    cov$acf[,,1]
#                  }) 
#   #   acov <- lapply(1:chains, 
#   #                  FUN = function(i) stan_prob_autocovariance(sims[, i]))
#   acov <- do.call(cbind, acov)
#   chain_mean <- apply(sims, 2, mean, na.rm=TRUE)
#   mean_var <- mean(acov[1,], na.rm = TRUE) * n_samples / (n_samples - 1) 
#   var_plus <- mean_var * (n_samples - 1) / n_samples
#   if (chains > 1) 
#     var_plus <- var_plus + var(chain_mean)
#   rho_hat_sum <- 0
#   for (t in 2:nrow(acov)) {
#     rho_hat <- 1 - (mean_var - mean(acov[t, ])) / var_plus
#     if (is.nan(rho_hat)) rho_hat <- 0
#     if (rho_hat < 0) break
#     rho_hat_sum <- rho_hat_sum + rho_hat
#   } 
#   ess <- chains * n_samples
#   if (rho_hat_sum > 0) ess <- ess / (1 + 2 * rho_hat_sum)
#   ess 
# } 


# right now this only works for chains that are the same length
split_rhat_rfun <- function(sims) {
  # Compute the split rhat for the diagnostics of converging; 
  # see the C++ code of split_potential_scale_reduction in chains.cpp.  
  # 
  # Args:
  #   sims: a 2-d array _without_ warmup samples (# iter * # chains) 
  # 
  # Note: 
  #   The R function wrapping the C++ implementation is defined 
  #   in chains.R with name rstan_splitrhat2_cpp 
  if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  half_n <- floor(n_samples / 2)
  # cat("n_samples=", n_samples, "\n"); cat("chains=", chains, "\n")
  # cat("half_n=", half_n, "\n")
  idx_2nd <- n_samples - half_n + 1
  
  split_chain_mean <- numeric(chains * 2)
  split_chain_var <- numeric(chains * 2)
  
  for (i in 1:chains) {
    split_chain_mean[i] <- mean(sims[1:half_n, i], na.rm=TRUE)
    split_chain_var[i] <- var(sims[1:half_n, i], na.rm=TRUE)
    split_chain_mean[chains + i] <- mean(sims[idx_2nd:n_samples, i], na.rm=TRUE)
    split_chain_var[chains + i] <- var(sims[idx_2nd:n_samples, i], na.rm=TRUE)
  } 
  var_between <- half_n * var(split_chain_mean)
  var_within <- mean(split_chain_var) 
  sqrt((var_between/var_within + half_n -1)/half_n)
} 

ess_rfun <- function(sims) {
  # Compute the effective sample size for samples of several chains 
  # for one parameter; see the C++ code of function  
  # effective_sample_size2 in chains.cpp 
  # 
  # Args:
  #   sims: a 2-d array _without_ warmup samples (# iter * # chains) 
  # 
  # Note: 
  #   The implementation in R uses acf in R to compute the autocovariance
  #   and the results might be a little bit different from that in stan. 
  #   The R function wrapping the C++ implementation is defined in 
  #   chains.R with name rstan_ess2_cpp 
  if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  
  acov <- lapply(1:chains, 
                 FUN = function(i) {
                   cov <- acf(sims[,i], lag.max = n_samples - 1, 
                              plot = FALSE, type = c("covariance"), na.action=na.pass) 
                   cov$acf[,,1]
                 }) 
  #   acov <- lapply(1:chains, 
  #                  FUN = function(i) stan_prob_autocovariance(sims[, i]))
  acov <- do.call(cbind, acov)
  chain_mean <- apply(sims, 2, mean, na.rm=TRUE)
  mean_var <- mean(acov[1,], na.rm = TRUE) * n_samples / (n_samples - 1) 
  var_plus <- mean_var * (n_samples - 1) / n_samples
  if (chains > 1) 
    var_plus <- var_plus + var(chain_mean)
  rho_hat_sum <- 0
  for (t in 2:nrow(acov)) {
    rho_hat <- 1 - (mean_var - mean(acov[t, ])) / var_plus
    if (is.nan(rho_hat)) rho_hat <- 0
    if (rho_hat < 0) break
    rho_hat_sum <- rho_hat_sum + rho_hat
  } 
  ess <- chains * n_samples
  if (rho_hat_sum > 0) ess <- ess / (1 + 2 * rho_hat_sum)
  ess 
} 


monitor <- function(sims, warmup = 0,# floor(dim(sims)[1] / 2), 
                    probs = c(0.025, 0.25, 0.50, 0.75, 0.975), 
                    digits_summary = 1, 
                    print = TRUE, ...) { 
  # print the summary for a general simulation results 
  # of 3-d array: # iter * # chains * # parameters 
  # Args:
  #   sims: a 3-d array described above 
  #   warmup: the number of iterations used for warmup 
  #   probs: probs of summarizing quantiles 
  #   print: print out the results
  # 
  # Return: 
  #   A summary array  
  dim_sims <- dim(sims)
  dimnames_sims <- dimnames(sims)
  parnames <- dimnames_sims[[3]]
  if (length(dim_sims) != 3) 
    stop("'sims' is not a 3-d array")
  if (warmup > dim_sims[1])
    stop("warmup is larger than the total number of iterations")
  num_par <- dim_sims[3]
  
  if (is.null(parnames)) parnames <- paste0("V", 1:num_par)
  sims_wow <- if (warmup >= 1) apply(sims, c(2, 3), FUN = function(x) x[-(1:warmup)]) else sims 
  m <- apply(sims_wow, 3, mean, na.rm=TRUE)
  sd <- sapply(1:num_par, FUN = function(i) sd(as.vector(sims_wow[,,i]), na.rm=TRUE)) 
  quan <- lapply(1:num_par, FUN = function(i) quantile(sims_wow[,,i], probs = probs, na.rm=TRUE))
  probs_str <- names(quan[[1]])
  quan <- do.call(rbind, quan)
  rhat <- sapply(1:num_par, FUN = function(i) split_rhat_rfun(sims_wow[,,i]))
  ess <- sapply(1:num_par, FUN = function(i) ess_rfun(sims_wow[,,i]))
  sem <- sd / sqrt(ess)
  
  summary <- cbind(m, sem, sd, quan, ess, rhat)
  colnames(summary) <- c("mean", "se_mean", "sd", probs_str, "n_eff", "Rhat")
  rownames(summary) <- parnames 
  if (print) {
    cat("Inference for the input samples (")
    cat(dim_sims[2], " chains: each with iter=", dim_sims[1], "; warmup=", warmup, "):\n\n", sep = "")
    # round n_eff to integers
    summary[, 'n_eff'] <- round(summary[, 'n_eff'], 0)
    print(round(summary, digits_summary), ...)
    
    cat("\nFor each parameter, n_eff is a crude measure of effective sample size,\n", 
        "and Rhat is the potential scale reduction factor on split chains (at \n",
        "convergence, Rhat=1).\n", sep = '')
  } 
  invisible(summary) 
} 


# load_cell <- function(dir_name, id, T, K) {
#   #fname = sprintf('output/%s.bin', run$suff_fit)
#   fname = paste0(dir_name, '/output.bin')
#   bin      = file(fname, "rb")
#   nwarmup  = readBin(bin, "integer")
#   nsamples = readBin(bin, "integer")
#   nparams  = readBin(bin, "integer")
#   
#   # seek to start of parameter names
#   seek(bin, where=(nsamples+nwarmup)*nparams*4, origin="current")
#   params    = readBin(bin, "character", n=nparams)
#   par_names = readBin(bin, "character", n=nparams)
#   
#   # find first parameter
#   pcol0 = 1
#   for (i in 1:10) {
#     if (grepl("__", params[i], fixed=TRUE) == FALSE) {
#       pcol0 = i
#       break
#     }
#   }
#   
#   params = params[pcol0:length(params)]
#   par_names = par_names[pcol0:length(par_names)]
#   
#   #   # get T, W
#   #   load(paste0(dir_name, '/input.rdata'))
#   W = K-1  
#   
#   # figure out which columns we want (cols)
#   g_cols    = which(par_names == 'g')
#   col0      = g_cols[(id-1)*T*W + 1]
#   ncols     = T*W
#   
#   # save start of samples
#   start = 3*4 + nwarmup*nparams*4
#   
#   # read appropriate columns
#   out = array(NA, dim=c(nsamples, ncols))
#   for (row in 1:nsamples) {
#     seek(bin, where=start+nparams*(row-1)*4+(pcol0+col0-2)*4, origin="start")
#     out[row,] = readBin(bin, "numeric", n=ncols, size=4)
#   }
#   colnames(out) <- params[col0:(col0+ncols-1)]
#   close(bin)
#   
#   out
# }

# # functions
# build_r_cell <- function(g_cell, N, T, K){
#   
#   W       = K-1
#   niter   = nrow(g_cell) 
#   
#   gc      = array(NA, dim=c(T, W, niter))
#   rc      = array(NA, dim=c(T, K, niter))
#   
#   for (k in 1:W){
#     #print(k)
#     
#     g_k_cols = seq(k, T*W, by=W)
#     gc[,k,]   = t(g_cell[,g_k_cols])    
#   }
#   
#   for (i in 1:niter){
#     
#     #if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
#     
#     sum_exp_g = rowSums(exp(gc[,,i]))
#     
#     rc[,,i] <- sum2one_constraint_cell(K, N, T, as.matrix(gc[,,i]), sum_exp_g) 
#   }
#   
#   return(list(rc=rc, gc=gc))
# }

# r[,,i] <- sum2one_constraint_cell(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
# additive log_ratio transformation
cppFunction('
  NumericMatrix sum2one_constraint_cell(int K, int N, int T, NumericMatrix g, NumericVector sum_exp_g) {
    //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
    NumericMatrix r(T, K);
    for (int k = 0; k<(K-1); k++)
      for (int j = 0; j<T; j++)
        r(j,k) = exp(g(j,k))/(1+sum_exp_g(j));

    for (int j = 0; j<T; j++)
      r(j,K-1) = 1.0 / (1 + sum_exp_g[j]);
   
return r;
  }
', verbose=TRUE, showOutput=TRUE)

cbind_matrices <- function(l) {
  m = matrix(0, nrow=0, ncol=ncol(l[[1]]))
  for (i in 1:length(l))
    m <- rbind(m, l[[i]])
  m
}


# functions
build_r_cell <- function(g_cell, N, T, K){
  
  W = K
  # W       = K-1
  niter   = nrow(g_cell) 
  
  gc      = array(NA, dim=c(T, W, niter))
  rc      = array(NA, dim=c(T, K, niter))
  
  for (k in 1:W){
    print(k)
    
    g_k_cols = seq(k, T*W, by=W)
    gc[,k,]   = t(g_cell[,g_k_cols])    
  }
  
  for (i in 1:niter){
    
    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
    
    sum_exp_g = rowSums(exp(gc[,,i]))
    
    rc[,,i] <- sum2one_constraint_cell_nb(K, N, T, as.matrix(gc[,,i]), sum_exp_g) 
  }
  
  return(list(rc=rc, gc=gc))
}



# r[,,i] <- sum2one_constraint_cell(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
# additive log_ratio transformation
cppFunction('
            NumericMatrix sum2one_constraint_cell(int K, int N, int T, NumericMatrix g, NumericVector sum_exp_g) {
            //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
            NumericMatrix r(T, K);
            for (int k = 0; k<(K-1); k++)
            for (int j = 0; j<T; j++)
            r(j,k) = exp(g(j,k))/(1+sum_exp_g(j));
            
            for (int j = 0; j<T; j++)
            r(j,K-1) = 1.0 / (1 + sum_exp_g[j]);
            
            return r;
            }
            ', verbose=TRUE)

# r[,,i] <- sum2one_constraint_cell(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
# additive log_ratio transformation
cppFunction('
            NumericMatrix sum2one_constraint_cell_nb(int K, int N, int T, NumericMatrix g, NumericVector sum_exp_g) {
            //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
            NumericMatrix r(T, K);
            for (int k = 0; k<K; k++)
            for (int j = 0; j<T; j++)
            r(j,k) = exp(g(j,k))/(1+sum_exp_g(j));
            return r;
            }
            ', verbose=TRUE)

load_cell <- function(run, cell_id, T, K) {
  #fname = sprintf('output/%s.bin', run$suff_fit)
  fname = sprintf('runs/%s/run1/output.bin', run$suff_run)
  bin      = file(fname, "rb")
  magic    = readBin(bin, "integer")  
  nwarmup  = readBin(bin, "integer")
  nsamples = readBin(bin, "integer")
  nparams  = readBin(bin, "integer")
  ndiag    = readBin(bin, "integer")
  
  if (magic > 0) {
    stop(paste0("Can't read ", fname, ": unexpected magic number (probably an old stan2bin file)."))
  }
  
  # seek to start of parameter names
  seek(bin, where=(nwarmup*nparams+1+ndiag+nsamples*nparams)*4, origin="current")
  params    = readBin(bin, "character", n=nparams)
  par_names = readBin(bin, "character", n=nparams)
  
  # find first parameter
  pcol0 = 1
  for (i in 1:10) {
    if (grepl("__", params[i], fixed=TRUE) == FALSE) {
      pcol0 = i
      break
    }
  }
  
  params = params[pcol0:length(params)]
  par_names = par_names[pcol0:length(par_names)]
  
  # get T, W
  load(sprintf('runs/%s/run1/input.rdata', run$suff_run))
  W = K  
  
  # figure out which columns we want (cols)
  g_cols    = which(par_names == 'g')
  col0      = g_cols[(cell_id-1)*T*W + 1]
  ncols     = T*W
  
  # save start of samples
  start = 5*4 + (nwarmup*nparams+1+ndiag)*4
  
  # read appropriate columns
  out = array(NA, dim=c(nsamples, ncols))
  for (row in 1:nsamples) {
    seek(bin, where=start+nparams*(row-1)*4+(pcol0+col0-2)*4, origin="start")
    out[row,] = readBin(bin, "numeric", n=ncols, size=4)
  }
  colnames(out) <- params[col0:(col0+ncols-1)]
  close(bin)
  
  out
}



############################################################################################################################################
## constrain pollen bacon ages
## only keep samples that are presettlement or older
## if no presettlement sample, keep 150 ybp or older
############################################################################################################################################

remove_post_settlement <- function(pollen, age_ps){
  
  # age_ps
  ids = unique(pollen$id)
  ncores = length(unique(pollen$id))
  
  pollen_ps = data.frame(matrix(NA, nrow=0, ncol=ncol(pollen))) 
  colnames(pollen_ps) = colnames(pollen)
  for (i in 1:ncores){

    pol_site = pollen[pollen$id == ids[i],]
    
    if ((any(age_ps$id == ids[i])) & (all(!is.na(pol_site$depth)))) {
      pol_site = pol_site[which(pol_site$depth >= age_ps$set_depth[which(age_ps$id == ids[i])]),]
    } else {
      pol_site = pol_site[pol_site$age_bacon >= 150, ]
    }
    
    pollen_ps = rbind(pollen_ps, pol_site)
  }
  
  if (any(pollen_ps$age_bacon <= 0)){
    pollen_ps$age_bacon[which(pollen_ps$age_bacon <= 0)] = 1
  }
  
  return(pollen_ps)
}


# drop the pollen samples that are 500 years older than the oldest geochron date
constrain_pollen <- function(pollen, age_con, nbeyond){
  
  drop_samples = vector(length=nrow(pollen))
  for (i in 1:nrow(pollen)){
    
    idx = match(pollen$id[i], age_con$id)
    
    if (is.na(idx)){
      print(paste0(i, ' WTF'))
      print(pollen$id[i])
    }
    
    drop_samples[i] = pollen$age_bacon[i] > (age_con$geo_age_max[idx] + nbeyond)
  }
  
  return(drop_samples)
}


##
## predicted pollen based on weighted neighborhoods using estimated pars
## 
pollen_preds_sp <- function(phi, gamma, w, sum_w, d, idx_cores, r){
  
  K       = dim(w)[1]
  N_cells = dim(w)[3]
  N_cores = dim(w)[2]
  
  T = dim(r)[3]
  
  r_new = array(NA, c(N_cores, K, T))
  preds = array(NA, c(N_cores, K, T))
  
  for (i in 1:N_cores){
    for (t in 1:T){
      print(i)
      
      out_sum = rep(0, K)
      for (k in 1:K){
        print(paste0('k = ', k))
        # pl_p1 = (b[k]-1) * (b[k]-2) / (2 * pi * a[k]  * a[k]) 
        
        for (j in 1:N_cells){ # changed N_hood to N_locs
          if (j != idx_cores[i]){
            # w = pl_p1 * (1 + d[j,i] / a[k]) ^ (-b[k])
            out_sum[k] <- out_sum[k] + w[k,i,j] * r[j,k,t]
          }  
        }
      }
      
      for (k in 1:K){
        r_new[i,k,t]  = gamma[k]*r[idx_cores[i],k,t] + (1-gamma[k])*out_sum[k]/sum_w[k]
        preds[i,k,t]  = phi[k]*r_new[i,k,t]        
      }
    }
  }
  
  # alpha = rowSums(preds)   
  
  #convert to proportions
  preds = apply(preds, c(1,3), function(x) x/sum(x))
  
  return(preds=preds)
}
