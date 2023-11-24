# Below is an example of how to use the code.
par(mfrow=c(1,2))
set.seed(012920)
HMC_beta_samples = HMC(n_samples = 2000, U_type = "beta", epsilon = 0.01, L = 2, 
                       current_q = 0.3, alpha = 50, beta = 50)
hist(HMC_beta_samples, breaks = 'scott', freq = FALSE,
     border = "#ffffff",
     col = "#81ddff",
     xlab = "Values of HMC Samples",
     main = "Histogram of HMC Samples (Beta Dist)")
curve(dbeta(x, 50, 50), col = "#FF6666", add = TRUE, lwd = 2)
set.seed(012920)
HMC_norm_samples = HMC(n_samples = 2000, U_type = "norm", epsilon = 0.2, L = 2, 
                       current_q = 0.3, mu = 0, sigma = 1)
hist(HMC_norm_samples, breaks = 'scott', freq = FALSE,
     border = "#ffffff",
     col = "#f79ee1",
     xlab = "Values of HMC Samples",
     main = "Histogram of HMC Samples (Normal Dist)")
curve(dnorm(x,0, 1), col = "#6699FF", add = TRUE, lwd = 2)