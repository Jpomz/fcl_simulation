# igp_sim2() was the original formulation and I know it works, retaining here for the time being. 
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
                     ebc = 2,
                     ebp = 1,
                     ecp = 2,
                     alphabc = 4,
                     alphabp = 2,
                     alphacp = 4,
                     betabc = 20,
                     betabp = 35,
                     betacp = 20,
                     s0 = 0.75,
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
    B_prime = s0*(N[1,] - wbc- wbp)
    C_prime = s0*(N[2,] - wcp)
    P_prime = s0*N[3,]
    
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
  return(list(sp_dynamics = dat,
              fcl = fcl_df))
}

# d <- igp_sim(k = 100, n_0 = c(100, 50, 20),
#              p_disturb = 0.0001)
# igp_sim(n_patch = 100, k = 100)



# original independent patch disturbance
# # Disturbance - reduce patch abundance by some proportion
# patch_extinction <- rbinom(n = n_patch,
#                            size = 1, prob = disturb_p)
# N[,which(patch_extinction==1)] <- 
#   N[,which(patch_extinction==1)] * (1 - disturb_mag)
# 