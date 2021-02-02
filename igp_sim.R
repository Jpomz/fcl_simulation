# Zhou model with IGP simulation

k_n_upstream <- function(k_base = 150,
                   k_c = 10,
                   k_min_exponent = 1.10,
                   k_max_exponent = 1.35,
                   r_max,
                   n_upstream,
                   n_patch){
  # carrying capacity ####
  # power law with exponent between 1.1-1.35 (Koening et al. 2019)
  # k_base = "minimum" carrying capacity
  # k_c = constant multiplier for power law
  
  # k_exp = exponent for carrying capacity function
  k_exp = runif(n = n_patch,
                min = k_min_exponent,
                max = k_max_exponent)
  k <- round(
    k_base + k_c * n_upstream^k_exp)
  b = (r_max - 1) /k
  b
}

# write function to estimate carrying capacity for three species based on different values

k_sp <- function(k_b = 500,
                 ebc = 2,
                 ebp = 1,
                 ecp = 2,
                 alphabc = 4,
                 alphabp = 2,
                 alphacp = 4,
                 betabc = 20,
                 betabp = 40,
                 betacp = 20){
  r_c = ebc * alphabc
  k_c = k_b*(r_c - 1) / betabc
  r_bp = ebp*alphabp
  k_bp = k_b*(r_bp - 1) / betabp
  r_cp = ecp*alphacp
  k_cp = k_c*(r_cp - 1) / betacp
  k_p_total = k_cp + k_bp
  
  c(r_c = r_c, r_bp = r_bp, r_cp = r_cp,
    k_b = k_b, k_c = k_c, k_bp = k_bp, k_cp = k_cp,
    k_p_total = k_p_total)
}

#k_sp()

max_reproduction <- function(ebc = 2,
                             ebp = 1,
                             ecp = 2,
                             alphabc = 4,
                             alphabp = 2,
                             alphacp = 4){
  max_c = ebc * alphabc
  p_b = ebp * alphabp
  p_c = ecp * alphacp
  max_p = p_c +p_b
  c(max_c = max_c, p_b = p_b, p_c = p_c, max_p = max_p)
}

#max_reproduction()


get_fcl <- function(N){# food chain length state
  patch_fcl = N
  patch_fcl[patch_fcl>0] = 1
  patch_fcl = patch_fcl*c(1, 2, 3)
  fcl = colSums(patch_fcl)
  fcl[fcl == 3] <- 2
  fcl[fcl == 4] <- 2.5
  fcl[fcl == 6] <- 3
  fcl
}

fcl_prop <- function(fcl){
  prop_0 = length(fcl[fcl == 0]) / length(fcl)
  prop_1 = length(fcl[fcl == 1]) / length(fcl)
  prop_2 = length(fcl[fcl == 2]) / length(fcl)
  prop_25 = length(fcl[fcl == 2.5]) / length(fcl)
  prop_3 = length(fcl[fcl == 3]) / length(fcl)
  fcl_prop = (c("p0" = prop_0,
                "p1" = prop_1,
                "p2" = prop_2,
                "p2.5" = prop_25,
                "p3" = prop_3))
  fcl_prop
}
# distance matrix, number of steps between patches
# below is from mcbrnet, make sure that output works in this function

# dist_mat <- riv_net$distance_matrix
# riv_patch <- riv_net$df_patch

# need to add rho_ab rho_wb somewhere
# use branch ID for correlation risk
# figure out what I want to record and make outputs list
#  make simulation to run over multiple time steps

# variables
# need to make this pretty strong to limit dispersal
# i.e. p_dispersal = 0.01, theta = 100
# strength of dispersal decline with distance

igp_sim <- function(n_patch = 20,
                    n_0 = NULL,
                    dist_mat = NULL,
                    k_function = NULL, #c("patches upstream", "random", "auto-correlated")
                    n_upstream = NULL,
                    k_base = 150,
                    k_c = 10,
                    k_min_exponent = 1.10,
                    k_max_exponent = 1.35,
                    k = 500,
                    p_dispersal = 0.1,
                    theta = 1,
                    r_max = 2.5,
                    ebc = 2,
                    ebp = 1,
                    ecp = 2,
                    alphabc = 4,
                    alphabp = 2,
                    alphacp = 4,
                    betabc = 20,
                    betabp = 35,
                    betacp = 20,
                    so = 0.75,
                    p_disturb = 1e-4,
                    mag_disturb = 0.5,
                    t = 1000,
                    plot_disturbance = TRUE,
                    plot_fcl = TRUE, 
                    plot_patch_dynamics = TRUE){
  # need to import functions:
  # dplyr:: bind_rows %>%
  # tidyr:: pivot_longer
  # ggplot2
library(tidyverse)
  
  n_sp = 3
  
  ## Dispersal ##
  # distance matrix
  if(is.null(dist_mat)){
    print("No distance matrix supplied, assuming 10x10 square landscape")
    x_coord = runif(n_patch, 0, 10)
    y_coord = runif(n_patch, 0, 10)
    dist_mat = data.matrix(dist(cbind(x_coord, y_coord),
                                diag = TRUE, upper = TRUE))
    #print(plot(x_coord, y_coord))
  }
  
  # dispersal matrix
  m_dispersal <- data.matrix(exp(-theta * dist_mat))
  diag(m_dispersal) <- 0
  
  # probability of dispersal
  if(length(p_dispersal) == 1){
    v_p_dispersal <- rep(x = p_dispersal, times = n_sp)
  }else{
    if(length(p_dispersal)!= n_sp) stop("length of p_dispersal should be 1  or n_sp")
    v_p_dispersal <- p_dispersal
  }
  
  ## carrying capacity ##
  if(is.null(k_function)){
  if(length(k) == 1){
    print("only one value of k supplied, assuming it is the same in all patches")
    b = (r_max - 1) /k
  }} else{  # for branching rier networks
    if(k_function == "patches upstream"){
      if(is.null(n_upstream)) stop("number of patches upstream must be supplied to n_upstream argument")
      if(length(n_upstream)!= n_patch) stop("length of n_upstream needs to equal n_patch")
      b = k_n_upstream(k_base = k_base,
                       k_c = k_c,
                       k_min_exponent = k_min_exponent,
                       k_max_exponent = k_max_exponent,
                       r_max = r_max,
                       n_upstream = n_upstream,
                       n_patch = n_patch)
    }} #==========# add other functions here, i.e. if == "random"...
  # if == "correlated" ...
  # others?
  
  ## initial community ##
  if(is.null(n_0)){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(k*0.8, k*0.5, k*0.25)),
                nrow = n_sp, ncol = n_patch)
  }
  if(length(n_0) == 1){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(n_0, n_0, n_0)),
                nrow = n_sp, ncol = n_patch)
  }
  if(length(n_0) == 3){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(n_0[1], n_0[2], n_0[3])),
                nrow = n_sp, ncol = n_patch)
  }

  # result output
  output <- list()
  fcl_list <- list()
  # record starting conditions
  fcl = get_fcl(N = N)
  out = cbind(1:n_patch, t(N), i = 1, k, patch_extinction = 0, fcl)
  colnames(out) <- c("patch", "B", "C", "P",
                     "time", "basal_k", "disturbance", "fcl")
  output[[1]] <- as.data.frame(out)
  
  ## simulation ##
  for (i in 2:t+1){
    # immigration / emmigration
    m_e_hat <- N * v_p_dispersal
    v_e_sum <- rowSums(m_e_hat)
    m_i_raw <- m_e_hat %*% m_dispersal
    v_i_sum <- rowSums(m_i_raw)
    v_i_sum[v_i_sum == 0] <- 1
    m_i_prob <- m_i_raw / v_i_sum
    m_i_hat <- m_i_prob * v_e_sum
    m_n_prime <- N + m_i_hat - m_e_hat
    
    m_n_prime[is.nan(m_n_prime)] <- 0
    m_n_prime[is.na(m_n_prime)] <- 0
    m_n_prime[m_n_prime <0 ] <- 0
    
    N <- matrix(rpois(n = n_sp * n_patch, lambda = m_n_prime),
                nrow = n_sp, ncol = n_patch)
    
    # predation
    wbc = alphabc * N[1,]*N[2,] / (betabc*N[2,] + N[1,])
    wbp = alphabp *N[1,]*N[3,]/(betabp*N[3,] +N[1,])
    wcp = alphacp*N[2,]*N[3,]/(betacp*N[3,]+N[2,])
    
    # survival
    B_prime = so*(N[1,] - wbc- wbp)
    C_prime = so*(N[2,] - wcp)
    P_prime = so*N[3,]
    
    # reproduction
    B_t1 = (r_max / (1+b*B_prime))*B_prime
    C_t1 = (ebc*alphabc*N[1,] / (betabc*N[2,] + N[1,]))*C_prime
    P_t1 = ((ebp*alphabp*N[1,] / (betabp*N[3,] + N[1,]))+
              (ecp*alphacp*N[2,] / (betacp*N[3,] + N[2,])))*P_prime
    
    # B_t1[is.nan(B_t1)] <- 0
    # C_t1[is.nan(C_t1)] <- 0
    # P_t1[is.nan(P_t1)] <- 0
    
    N = rbind(B_t1, C_t1, P_t1)
    
    N[is.nan(N)] <- 0
    N[is.na(N)] <- 0
    N[N <0 ] <- 0
    
    # Disturbance - reduce patch abundance by some proportion
    patch_extinction <- rbinom(n = n_patch, size = 1, prob = p_disturb)
    N[,which(patch_extinction==1)] <- 
      N[,which(patch_extinction==1)] * (1 - mag_disturb)
    
    # number of individuals  in patch x at time t
    N = matrix(rpois(n = n_sp * n_patch, lambda = N),
               nrow = n_sp, ncol = n_patch)
    
    # if basal species is extinct in a patch
    # make abundance of C and P = 0
    N[,N[1,] == 0] <- 0
    
    
    # food chain length state
    fcl = fcl_prop(get_fcl(N = N))
    
    # path_dynamics output
    out = cbind(1:n_patch, t(N), i, k, patch_extinction, fcl)
        colnames(out) <- c("patch", "B", "C", "P",
                       "time", "basal_k", "disturbance", "fcl")
    output[[i]] <- as.data.frame(out)
    fcl_list[[i]] <- fcl
  }
  
  # make sure that these functions are properly imported
  dat <- dplyr::bind_rows(output)
  dat <- tidyr::pivot_longer(dat, 2:4, names_to = "species")
  
  fcl_df <- data.frame(bind_rows(fcl_list), time = 2:t)
  fcl_df <- pivot_longer(fcl_df, 1:5, names_to = "fcl_state")
  

# Plots -------------------------------------------------------------------

  # plot of food chain length
  if(plot_fcl == TRUE){
    fcl_plot <- 
      ggplot(fcl_df,
             aes(y = value, x = time, color = fcl_state)) +
      geom_step() +
      facet_wrap(.~fcl_state) +
      theme_bw() +
      labs(y = "proportion of patches",
           title = "Food chain length") +
      NULL
    print(fcl_plot)
  }
    # plot of patch dynamics
    if(plot_patch_dynamics == TRUE){
      
      plot_dat <- dat %>% filter(patch %in% sample(1:n_patch, 5))
      dist_dat <- dplyr::filter(plot_dat, disturbance == 1)
    
    if(nrow(dist_dat) >0){
      patch_plot <- ggplot(plot_dat, 
                           aes(y = value,
                               x = time, 
                               color = species))+
        geom_line() +
        geom_point(inherit.aes = FALSE,
                   data = dist_dat,
                   mapping = aes(x = time,
                                 y = -15),
                   shape = 2,
                   size = 2) +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 5 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
      } else {
        patch_plot <- ggplot(plot_dat, 
                             aes(y = value,
                                 x = time, 
                                 color = species))+
          geom_line() +
          facet_wrap(.~patch, labeller = label_both) +
          labs(y = "Abundance",
               title = "Dynamics for 5 random patches") +
          theme_bw() +
          NULL
        print(patch_plot)
      }
    }
    # plot of patches with disturbances
    if(plot_disturbance == TRUE){
      dist_patch <- dplyr::filter(dat, disturbance == 1) %>%
        select(patch) %>% unique() %>% pull()
      dist_dat <- filter(dat, patch %in% dist_patch)
      dist_point <- dplyr::filter(dat, disturbance == 1)
      dist_plot <- ggplot(dist_dat, 
                           aes(y = value,
                               x = time, 
                               color = species))+
        geom_line() +
        geom_point(inherit.aes = FALSE,
                   data = dist_point,
                   mapping = aes(x = time,
                                 y = -15),
                   shape = 2,
                   size = 2) +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "All patches experiencing disturbances") +
        theme_bw() +
        NULL
      print(dist_plot)
    }
  return(list(sp_dynamics = dat,
              fcl = fcl_df))
  }

# d <- igp_sim(k = 100, n_0 = c(100, 50, 20),
#              p_disturb = 0.0001)
# igp_sim(n_patch = 100, k = 100)
<<<<<<< HEAD


igp_sim2 <- function(n_patch = 20,
                    n_0 = NULL,
                    dist_mat = NULL,
                    k_function = NULL, #c("patches upstream", "random", "auto-correlated")
                    n_upstream = NULL,
                    k_base = 150,
                    k_c = 10,
                    k_min_exponent = 1.10,
                    k_max_exponent = 1.35,
                    k = 500,
                    p_dispersal = 0.1,
                    theta = 1,
                    r_max = 2.5,
                    erc = 2,
                    alpha = 4,
                    beta = 20,
                    P_pref = 0.25, # preference of B over C
                    so = 0.75,
                    p_disturb = 1e-4,
                    mag_disturb = 0.5,
                    t = 1000,
                    plot_disturbance = TRUE,
                    plot_fcl = TRUE, 
                    plot_patch_dynamics = TRUE){
  # need to import functions:
  # dplyr:: bind_rows %>%
  # tidyr:: pivot_longer
  # ggplot2
  library(tidyverse)
  
  n_sp = 3
  
  ## Dispersal ##
  # distance matrix
  if(is.null(dist_mat)){
    print("No distance matrix supplied, assuming 10x10 square landscape")
    x_coord = runif(n_patch, 0, 10)
    y_coord = runif(n_patch, 0, 10)
    dist_mat = data.matrix(dist(cbind(x_coord, y_coord),
                                diag = TRUE, upper = TRUE))
    #print(plot(x_coord, y_coord))
  }
  
  # dispersal matrix
  m_dispersal <- data.matrix(exp(-theta * dist_mat))
  diag(m_dispersal) <- 0
  
  # probability of dispersal
  if(length(p_dispersal) == 1){
    v_p_dispersal <- rep(x = p_dispersal, times = n_sp)
  }else{
    if(length(p_dispersal)!= n_sp) stop("length of p_dispersal should be 1  or n_sp")
    v_p_dispersal <- p_dispersal
  }
  
  ## carrying capacity ##
  if(is.null(k_function)){
    if(length(k) == 1){
      print("only one value of k supplied, assuming it is the same in all patches")
      b = (r_max - 1) /k
    }} else{  # for branching rier networks
      if(k_function == "patches upstream"){
        if(is.null(n_upstream)) stop("number of patches upstream must be supplied to n_upstream argument")
        if(length(n_upstream)!= n_patch) stop("length of n_upstream needs to equal n_patch")
        b = k_n_upstream(k_base = k_base,
                         k_c = k_c,
                         k_min_exponent = k_min_exponent,
                         k_max_exponent = k_max_exponent,
                         r_max = r_max,
                         n_upstream = n_upstream,
                         n_patch = n_patch)
      }} #==========# add other functions here, i.e. if == "random"...
  # if == "correlated" ...
  # others?
  
  ## initial community ##
  if(is.null(n_0)){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(k*0.8, k*0.5, k*0.25)),
                nrow = n_sp, ncol = n_patch)
  }
  if(length(n_0) == 1){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(n_0, n_0, n_0)),
                nrow = n_sp, ncol = n_patch)
  }
  if(length(n_0) == 3){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(n_0[1], n_0[2], n_0[3])),
                nrow = n_sp, ncol = n_patch)
  }
  
  # result output
  output <- list()
  fcl_list <- list()
  # record starting conditions
  fcl = get_fcl(N = N)
  out = cbind(1:n_patch, t(N), i = 1, k, patch_extinction = 0, fcl)
  colnames(out) <- c("patch", "B", "C", "P",
                     "time", "basal_k", "disturbance", "fcl")
  output[[1]] <- as.data.frame(out)
  
  ## simulation ##
  for (i in 2:t+1){
    # immigration / emmigration ####
    m_e_hat <- N * v_p_dispersal
    v_e_sum <- rowSums(m_e_hat)
    m_i_raw <- m_e_hat %*% m_dispersal
    v_i_sum <- rowSums(m_i_raw)
    v_i_sum[v_i_sum == 0] <- 1
    m_i_prob <- m_i_raw / v_i_sum
    m_i_hat <- m_i_prob * v_e_sum
    m_n_prime <- N + m_i_hat - m_e_hat
    
    m_n_prime[is.nan(m_n_prime)] <- 0
    m_n_prime[is.na(m_n_prime)] <- 0
    m_n_prime[m_n_prime <0 ] <- 0
    
    N <- matrix(rpois(n = n_sp * n_patch, lambda = m_n_prime),
                nrow = n_sp, ncol = n_patch)
    
    # predation ####
    wbc = alpha * N[1,]*N[2,] / (beta*N[2,] + N[1,])
    wbp = P_pref * (alpha *N[1,]*N[3,]/(beta*N[3,] +N[1,]))
    wcp = (1 - P_pref) * (alpha*N[2,]*N[3,]/(beta*N[3,]+N[2,]))
    
    # survival ####
    B_prime = so*(N[1,] - wbc - wbp)
    C_prime = so*(N[2,] - wcp)
    P_prime = so*N[3,]
    
    # reproduction ####
    B_t1 = (r_max / (1+b*B_prime))*B_prime
    C_t1 = (erc*alpha*N[1,] / 
              (beta*N[2,] + N[1,])) * # B converted to C
      C_prime # number of C 
    P_t1 = ((P_pref * (erc*alpha*N[1,] / 
               (beta*N[3,] + N[1,]))) + # B converted to P
              ((1 - P_pref) * (erc*alpha*N[2,] / 
                 (beta*N[3,] + N[2,])))) * # C converted to P
      P_prime # number of P
    
    # B_t1[is.nan(B_t1)] <- 0
    # C_t1[is.nan(C_t1)] <- 0
    # P_t1[is.nan(P_t1)] <- 0
    
    N = rbind(B_t1, C_t1, P_t1)
    
    N[is.nan(N)] <- 0
    N[is.na(N)] <- 0
    N[N <0 ] <- 0
    
    # Disturbance - reduce patch abundance by some proportion
    patch_extinction <- rbinom(n = n_patch, size = 1, prob = p_disturb)
    N[,which(patch_extinction==1)] <- 
      N[,which(patch_extinction==1)] * (1 - mag_disturb)
    
    # number of individuals  in patch x at time t
    N = matrix(rpois(n = n_sp * n_patch, lambda = N),
               nrow = n_sp, ncol = n_patch)
    
    # if basal species is extinct in a patch
    # make abundance of C and P = 0
    N[,N[1,] == 0] <- 0
    
    
    # food chain length state
    fcl = fcl_prop(get_fcl(N = N))
    
    # path_dynamics output
    out = cbind(1:n_patch, t(N), i, k, patch_extinction, fcl)
    colnames(out) <- c("patch", "B", "C", "P",
                       "time", "basal_k", "disturbance", "fcl")
    output[[i]] <- as.data.frame(out)
    fcl_list[[i]] <- fcl
  }
  
  # make sure that these functions are properly imported
  dat <- dplyr::bind_rows(output)
  dat <- tidyr::pivot_longer(dat, 2:4, names_to = "species")
  
  fcl_df <- data.frame(bind_rows(fcl_list), time = 2:t)
  fcl_df <- pivot_longer(fcl_df, 1:5, names_to = "fcl_state")
  
  
  # Plots -------------------------------------------------------------------
  
  # plot of food chain length
  if(plot_fcl == TRUE){
    fcl_plot <- 
      ggplot(fcl_df,
             aes(y = value, x = time, color = fcl_state)) +
      geom_step() +
      facet_wrap(.~fcl_state) +
      theme_bw() +
      labs(y = "proportion of patches",
           title = "Food chain length") +
      NULL
    print(fcl_plot)
  }
  # plot of patch dynamics
  if(plot_patch_dynamics == TRUE){
    
    plot_dat <- dat %>% filter(patch %in% sample(1:n_patch, 5))
    dist_dat <- dplyr::filter(plot_dat, disturbance == 1)
    
    if(nrow(dist_dat) >0){
      patch_plot <- ggplot(plot_dat, 
                           aes(y = value,
                               x = time, 
                               color = species))+
        geom_line() +
        geom_point(inherit.aes = FALSE,
                   data = dist_dat,
                   mapping = aes(x = time,
                                 y = -15),
                   shape = 2,
                   size = 2) +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 5 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
    } else {
      patch_plot <- ggplot(plot_dat, 
                           aes(y = value,
                               x = time, 
                               color = species))+
        geom_line() +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 5 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
    }
  }
  # plot of patches with disturbances
  if(plot_disturbance == TRUE){
    dist_patch <- dplyr::filter(dat, disturbance == 1) %>%
      select(patch) %>% unique() %>% pull()
    dist_dat <- filter(dat, patch %in% dist_patch)
    dist_point <- dplyr::filter(dat, disturbance == 1)
    dist_plot <- ggplot(dist_dat, 
                        aes(y = value,
                            x = time, 
                            color = species))+
      geom_line() +
      geom_point(inherit.aes = FALSE,
                 data = dist_point,
                 mapping = aes(x = time,
                               y = -15),
                 shape = 2,
                 size = 2) +
      facet_wrap(.~patch, labeller = label_both) +
      labs(y = "Abundance",
           title = "All patches experiencing disturbances") +
      theme_bw() +
      NULL
    print(dist_plot)
  }
  return(list(sp_dynamics = dat,
              fcl = fcl_df))
}
=======
>>>>>>> 5f92ca836e7cefaec374c026ae85c0d882085d85
