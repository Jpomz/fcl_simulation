# ebc = 2,
# ebp = 1,
# ecp = 2,
# alphabc = 4,
# alphabp = 2,
# alphacp = 4,
# betabc = 20,
# betabp = 35,
# betacp = 20,

set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # 
        alphabc = 4,
        alphabp = 2,
        alphacp = 4,
        betabc = 20,
        betabp = 35,
        betacp = 20,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# one (?) extinctions

# Decreasing beta by half
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # reduce population by 75%
        alphabc = 4,
        alphabp = 2,
        alphacp = 4,
        betabc = 10,
        betabp = 17,
        betacp = 10,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# Consistent decrease in patches with P

# decrease alpha by 1/2
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # reduce population by 75%
        alphabc = 2,
        alphabp = 1,
        alphacp = 2,
        betabc = 20,
        betabp = 35,
        betacp = 20,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# P gone, all patches B + C

# double beta
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # 
        alphabc = 4,
        alphabp = 2,
        alphacp = 4,
        betabc = 40,
        betabp = 70,
        betacp = 40,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# P gone, all patches B+C

# double alpha
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # 
        alphabc = 8,
        alphabp = 4,
        alphacp = 8,
        betabc = 20,
        betabp = 35,
        betacp = 20,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# All P Gone, ~50% C gone, some B gone

# double alpha with dispersal
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0.1,
        mag_disturb = 0, # 
        alphabc = 8,
        alphabp = 4,
        alphacp = 8,
        betabc = 20,
        betabp = 35,
        betacp = 20,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# mostly FCL = 3, some = 2

# half beta, double alpha
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # reduce population by 75%
        alphabc = 8,
        alphabp = 4,
        alphacp = 8,
        betabc = 10,
        betabp = 17,
        betacp = 10,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# all extinct

# half beta, double alpha with dispersal
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0.1,
        mag_disturb = 0, # reduce population by 75%
        alphabc = 8,
        alphabp = 4,
        alphacp = 8,
        betabc = 10,
        betabp = 17,
        betacp = 10,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# all extinct


# ~3/4  beta, double alpha 
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0,
        mag_disturb = 0, # reduce population by 75%
        alphabc = 8,
        alphabp = 4,
        alphacp = 8,
        betabc = 15,
        betabp = 27,
        betacp = 15,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# 50% no species, 50% B


# ~3/4  beta, double alpha with dispersal
set.seed(251)
igp_sim(n_patch = 100, 
        p_dispersal = 0.1,
        mag_disturb = 0, # reduce population by 75%
        alphabc = 8,
        alphabp = 4,
        alphacp = 8,
        betabc = 15,
        betabp = 27,
        betacp = 15,
        plot_disturbance = FALSE,
        plot_patch_dynamics = FALSE,
        plot_fcl = TRUE)
# lots of dynamics