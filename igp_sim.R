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

prey_preference <- function(e1, e2, N1, N2){
  # e1 = ebp, N1 = B
  # e2 = ecp, N2 = C
  # returns P preference of B over C
  e1*N1 / (e1*N1 + e2*N2)
}


dispersal_n <- function(N, v_p_dispersal, theta, dist_mat){
  if(length(theta) != 1 & length(theta) !=3){
    stop("Theta must be length 1 or 3")
  }
  if(length(theta) == 1){
    v_theta = rep(theta, 3)
  }
  if(length(theta) == 3){
    v_theta = theta
  }
  
  # dispersal matrices per species
  m_b_dispersal <- data.matrix(exp(-v_theta[1] * dist_mat))
  diag(m_b_dispersal) <- 0
  m_c_dispersal <- data.matrix(exp(-v_theta[2] * dist_mat))
  diag(m_c_dispersal) <- 0
  m_p_dispersal <- data.matrix(exp(-v_theta[3] * dist_mat))
  diag(m_p_dispersal) <- 0
  
  m_e_hat <- N * v_p_dispersal
  v_e_sum <- rowSums(m_e_hat)
  # repeat for each species
  m_i_b_raw <- m_e_hat[1,] %*% m_b_dispersal
  m_i_c_raw <- m_e_hat[2,] %*% m_c_dispersal
  m_i_p_raw <- m_e_hat[3,] %*% m_p_dispersal
  m_i_raw <- rbind(m_i_b_raw, m_i_c_raw, m_i_p_raw)
  
  
  v_i_sum <- rowSums(m_i_raw)
  v_i_sum[v_i_sum == 0] <- 1
  m_i_prob <- m_i_raw / v_i_sum
  m_i_hat <- m_i_prob * v_e_sum
  m_n_prime <- N + m_i_hat - m_e_hat
  
  m_n_prime[is.nan(m_n_prime)] <- 0
  m_n_prime[is.na(m_n_prime)] <- 0
  m_n_prime[m_n_prime <0 ] <- 0
  
  m_n_prime
}


igp_sim <- function(n_patch = 20,
                    n_0 = NULL,
                    dist_mat = NULL,
                    adjacency_matrix = NULL,
                    k_function = NULL, 
                    #c("patches upstream", "random", "auto-correlated")
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
                    ebp = 2,
                    ecp = 2,
                    alphabc = 4,
                    alphap = 4,
                    betabc = 20,
                    betap = 20,
                    P_pref = 0.25, # preference of B over C
                    s0 = 0.75,
                    disturb_p = 1e-4,
                    disturb_mag = 0.5,
                    disturb_rho = 1, # 2D habitats
                    disturb_decay = 0.75, # downstream
                    t = 1000,
                    plot_disturbance = TRUE,
                    plot_fcl = TRUE, 
                    plot_patch_dynamics = TRUE){
  # need to import functions:
  # dplyr:: bind_rows %>%
  # tidyr:: pivot_longer
  # ggplot2
  library(tidyverse)
  
  param_df <- data.frame(
    alphabc = alphabc, betabc = betabc, ebc = ebc,
    alphap = alphap, betap = betap, ebp = ebp, ecp =ecp,
    P_pref = ifelse(is.null(P_pref), NA, P_pref),
    p_dispersal = p_dispersal, s0 = s0,
    disturb_p = disturb_p, disturb_mag = disturb_mag)
  
  n_sp = 3
  
  ## distance matrix ##
  if(!is.null(dist_mat)){
    if (!is.matrix(dist_mat)) stop("distance matrix should be
                                          provided as matrix")
    if (nrow(dist_mat) != n_patch) stop(
    "invalid dimension: distance matrix must have a dimension of 
    n_patch * n_patch")
    if (any(diag(dist_mat) != 0)) stop(
      "invalid distance matrix: diagonal elements must be zero")
    # df_xy_coord <- NULL
    dist_mat <- dist_mat
    # add / make adj_matrix here
  }
  if(is.null(dist_mat)){
    message("No distance matrix supplied, assuming 10x10 square landscape")
    x_coord = runif(n_patch, 0, 10)
    y_coord = runif(n_patch, 0, 10)
    dist_mat = data.matrix(dist(cbind(x_coord, y_coord),
                                diag = TRUE, upper = TRUE))
    # convert to adjacency matrix for disturbance lower
  }

  # probability of dispersal
  if(length(p_dispersal) == 1){
    v_p_dispersal <- rep(x = p_dispersal, times = n_sp)
    message(
      paste("1 value of p_dispersal supplied:",
            p_dispersal, 
            "for all three species"))
  }else{
    if(length(p_dispersal)!= n_sp) stop("length of p_dispersal should be 1  or n_sp")
    v_p_dispersal <- p_dispersal
    message(
      paste("3 values of theta supplied: B_theta =", p_dispersal[1],
            "C_theta =", p_dispersal[2],
            "and P_theta =", p_dispersal[3]))
  }
  
  # distance decay of dispersal, theta
  if(length(theta) == 1){
    message(
      paste("1 value of theta supplied: B_theta = C_theta = P_theta =",
            theta))
  }
  if(length(theta) == 3){
    message(
      paste("3 values of theta supplied: B_theta =", v_theta[1],
            "C_theta =", v_theta[2],
            "and P_theta =", v_theta[3]))
  }
  
  ## carrying capacity ##
  if(is.null(k_function)){
    if(length(k) == 1){
      message("only one value of k supplied, assuming it is the same in all patches")
      b = (r_max - 1) /k
    }} else{  # for branching river networks
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
  
  # predator preference fixed or variable
  if(is.null(P_pref)){
    fixed_P_pref = FALSE
    message("Predator preference varies with resource abundance")
  }
  if(is.numeric(P_pref) & length(P_pref == 1)){
    fixed_P_pref = TRUE
    message(paste("Predator preference is set at", P_pref))
  }
  
  ## initial community ##
  if(is.null(n_0)){
    N <- matrix(rpois(n = n_sp*n_patch,
                      lambda = c(k*0.8, k*0.5, k*0.25)),
                nrow = n_sp, ncol = n_patch)
    message("initial community = 0.8, 0.5 and 0.25% K")
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
    m_n_prime <- dispersal_n(N = N, v_p_dispersal = v_p_dispersal,
                             theta = theta, dist_mat = dist_mat)
    
    N <- matrix(rpois(n = n_sp * n_patch, lambda = m_n_prime),
                nrow = n_sp, ncol = n_patch)
    
    # predation ####
    if(fixed_P_pref == TRUE){
      P_pref = P_pref
    }
    if(fixed_P_pref == FALSE){
      P_pref = prey_preference(e1 = ebp, e2 = ecp, N1 = N[1,], N2 = N[2,])
    }
    
    wbc = alphabc * N[1,]*N[2,] / (betabc*N[2,] + N[1,])
    wbp = P_pref * (alphap *N[1,]*N[3,]/(betap*N[3,] +N[1,]))
    wcp = (1 - P_pref) * (alphap*N[2,]*N[3,]/(betap*N[3,]+N[2,]))
    
    # survival ####
    B_prime = s0*(N[1,] - wbc - wbp)
    C_prime = s0*(N[2,] - wcp)
    P_prime = s0*N[3,]
    
    # reproduction ####
    B_t1 = (r_max / (1+b*B_prime))*B_prime
    C_t1 = (ebc*alphabc*N[1,] / 
              (betabc*N[2,] + N[1,])) * # B converted to C
      C_prime # number of C 
    P_t1 = ((P_pref * (ebp*alphap*N[1,] / 
                         (betap*N[3,] + N[1,]))) + # B converted to P
              ((1 - P_pref) * (ecp*alphap*N[2,] / 
                                 (betap*N[3,] + N[2,])))) * # C converted to P
      P_prime # number of P
    
    N = rbind(B_t1, C_t1, P_t1)
    
    N[is.nan(N)] <- 0
    N[is.na(N)] <- 0
    N[N <0 ] <- 0
    
    # disturbance ####
    ## add an if statement or 2D habitats
    if(is.null(adjacency_matrix)){
      patch_extinction <- rbinom(n = n_patch,
                                 size = 1, prob = disturb_p)
      patch_disturb_id <- which(patch_extinction==1)
      if(length(patch_disturb_id) == 0){
        N <-  N
      } 
      if(length(patch_disturb_id) != 0){
        disturb_m <- data.matrix(exp(-disturb_rho * dist_mat))
        disturb_m <- as.matrix(disturb_m[,patch_disturb_id])
        N <- N*(1 - rowSums(disturb_m)*disturb_mag)[col(N)]
        }
      }
    
    if(!is.null(adjacency_matrix)){
      #message("adj matrix not null")
      # Disturbance "decay" down stream
      patch_extinction <- rbinom(n = n_patch,
                                 size = 1,
                                 prob = disturb_p)
      patch_disturb_id <- which(patch_extinction==1)
      
      m_adj_up <- adjacency_matrix
      m_adj_up[lower.tri(m_adj_up)] <- 0
      
      # vector to hold disturbance magnitudes
      disturb_dummy <- disturb_v <-  rep(0, n_patch)
      # "source" disturbance magnitude
      disturb_dummy[patch_disturb_id] <- disturb_mag
      
      for(np in n_patch:1){
        disturb_v <- disturb_v + disturb_dummy
        disturb_dummy <- m_adj_up %*% (disturb_decay * disturb_dummy)
      }
      
      disturb_v[disturb_v >=1] <- 1
      disturb_v <- as.vector(disturb_v)
      
      N <- N*(1 - disturb_v)[col(N)]
      #message("end of disturb loop")
    }
    
    # N[,which(patch_extinction==1)] <- 
    #   N[,which(patch_extinction==1)] * (1 - mag_disturb)
    
    # end of cycle ####
    # number of individuals  in patch x at time t
    N <-  matrix(rpois(n = n_sp * n_patch, lambda = N),
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
    #message("end of loop")
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
    if(nrow(dist_dat) >0){
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
    } else {message("No disturbances occured")}
  }
  
  # return ####
  return(list(sp_dynamics = dat,
              fcl = fcl_df,
              sim_params = param_df))
}

