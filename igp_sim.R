# Zhou model with IGP simulation

k_n_upstream <- function(k_base = 150,
                   k_c = 10,
                   min_exponent = 1.10,
                   max_exponent = 1.35,
                   n_patch_upstream){
  # carrying capacity ####
  # power law with exponent between 1.1-1.35 (Koening et al. 2019)
  # k_base = "minimum" carrying capacity
  # k_c = constant multiplier for power law
  
  # k_exp = exponent for carrying capacity function
  k_exp = runif(n = n_patch,
                min = min_exponent,
                max = max_exponent)
  k <- round(
    k_base + k_c * n_patch_upstream^k_exp)
  b = (a - 1) /k
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

k_sp()

igp_sim <- function(n_patch = 20,
                    n_sp = 3,
                    dist_mat = NULL,
                    k,
                    p_dispersal = 0.0,
                    theta = 1,
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
                    p_diturb = 1e-4,
                    mag_disturb = 0.25,
                    a = 2.5,
                    t = 10000)


# distance matrix, number of steps between patches
dist_mat <- riv_net$distance_matrix
riv_patch <- riv_net$df_patch

# need to add rho_ab rho_wb somewhere
# use branch ID for correlation risk
# figure out what I want to record and make outputs list
#  make simulation to run over multiple time steps

# variables
# need to make this pretty strong to limit dispersal
# i.e. p_dispersal = 0.01, theta = 100
 # strength of dispersal decline with distance

# dispersal
if(length(p_dispersal == 1)){
  v_p_dispersal <- rep(x = p_dispersal, times = n_sp)
}

# just once
m_dispersal <- data.matrix(exp(-theta * dist_mat))
diag(m_dispersal) <- 0

# initial community
N <- matrix(rpois(n = n_sp*n_patch,
                   lambda = c(50, 50, 50)),
             nrow = n_sp, ncol = n_patch)

# result output
output <- list()
# record starting conditions
out = cbind(1:n_patch, t(N), i = 1, k, patch_extinction = 0)
colnames(out) <- c("patch", "B", "C", "P", "time", "basal_k", "disturbance")
output[[1]] <- as.data.frame(out)

# simulation
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
  B_t1 = (a / (1+b*B_prime))*B_prime
  C_t1 = (ebc*alphabc*N[1,] / (betabc*N[2,] + N[1,]))*C_prime
  P_t1 = ((ebp*alphabp*N[1,] / (betabp*N[3,] + N[1,]))+
            (ecp*alphacp*N[2,] / (betacp*N[3,] + N[2,])))*P_prime
  
  B_t1[is.nan(B_t1)] <- 0
  C_t1[is.nan(C_t1)] <- 0
  P_t1[is.nan(P_t1)] <- 0
  
  N = rbind(B_t1, C_t1, P_t1)
  
  # Disturbance - reduce patch abundance by some proportion
  patch_extinction <- rbinom(n = n_patch, size = 1, prob = p_diturb)
  N[,which(patch_extinction==1)] <- 
    N[,which(patch_extinction==1)] * mag_disturb
  
  # number of individuals  in patch x at time t
  N = matrix(rpois(n = n_sp * n_patch, lambda = N),
             nrow = n_sp, ncol = n_patch)
  
  # if basal species is extinct in a patch
  # make abundance of C and P = 0
  N[,N[1,] == 0] <- 0
  
  # path_dynamics output
  out = cbind(1:n_patch, t(N), i, k, patch_extinction)
  colnames(out) <- c("patch", "B", "C", "P", "time", "basal_k", "disturbance")
  output[[i]] <- as.data.frame(out)
}

dat <- bind_rows(output)
dat <- pivot_longer(dat, 2:4, names_to = "species")
dist_dat <- dat %>% filter(disturbance == 1)

if(nrow(dist_dat) >0){
  ggplot(dat, 
         aes(y = value,
             x = time, 
             color = species))+
    geom_line()+
    geom_point(inherit.aes = FALSE,
               data = dist_dat,
               mapping = aes(x = time,
                             y = -15),
               shape = 2,
               size = 2) +
    facet_wrap(.~patch) +
    #facet_wrap(interaction(basal_K+patch)~.) +
    # labs(y = "Abundance",
    #      title = "Patch Carrying Capacity for B") +
    theme_bw() +
    NULL
} else {
  ggplot(dat, 
         aes(y = value,
             x = time, 
             color = species))+
    geom_line()+
    facet_wrap(.~patch) +
    #facet_wrap(interaction(basal_K+patch)~.) +
    # labs(y = "Abundance",
    #      title = "Patch Carrying Capacity for B") +
    theme_bw() +
    NULL
}