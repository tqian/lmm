# 2018.10.22 Tianchen Qian

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it

dgm_with_treatment <- function(sample_size, total_T, dgm_type) {
    
    # dgm_type is in c(1,2,3,4)
    stopifnot(dgm_type %in% c(1,2,3,4))
    
    # dgm_type = 1 or 3
    alpha_0 <- - 1
    alpha_1 <- - 0.3
    beta_0 <- 0.5
    beta_1 <- 0.1
    sigma_b0 <- 2
    sigma_b1 <- 0
    sigma_b2 <- 1
    sigma_b3 <- 0
    sigma_eps <- 1
    
    if (dgm_type == 2) {
        sigma_b1 <- sigma_b3 <- 0.5
    }
    if (dgm_type == 4) {
        sigma_b2 <- 0
    }
    
    
    prob_a <- 0.5
    
    df_names <- c("userid", "day", "X", "prob_A", "A", "Y", "b0", "b1", "b2", "b3", "eps", "delta")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$day <- rep(1:total_T, times = sample_size)
    
    # uncorrelated random effects
    b_0i <- rnorm(sample_size, mean = 0, sd = sigma_b0)
    b_1i <- rnorm(sample_size, mean = 0, sd = sigma_b1)
    b_2i <- rnorm(sample_size, mean = 0, sd = sigma_b2)
    b_3i <- rnorm(sample_size, mean = 0, sd = sigma_b3)
    
    # b_1i[b_1i > 2] <- 2
    # b_1i[b_1i < -2] <- -2
    # 
    # b_3i[b_3i > 2] <- 2
    # b_3i[b_3i < -2] <- -2
    
    for (t in 1:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        
        if (dgm_type == 1) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- rep(prob_a, sample_size)
        } else if (dgm_type == 2) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- ifelse(dta$X[row_index] > - 1.27, 0.7, 0.3)
        } else if (dgm_type %in% c(3,4)) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size, mean = b_0i) # X involves b_i!!
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size, mean = b_0i) # X involves b_i!!
            }
            dta$prob_A[row_index] <- rep(prob_a, sample_size)
        }
        
        dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
        
        dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = sigma_eps)
        
        dta$delta[row_index] <- beta_0 + beta_1 * dta$X[row_index] + b_2i + b_3i * dta$X[row_index]
        
        dta$Y[row_index] <- alpha_0 + alpha_1 * dta$X[row_index] + b_0i + b_1i * dta$X[row_index] + 
            dta$A[row_index] * dta$delta[row_index] + dta$eps[row_index]
        
        dta$b0[row_index] <- b_0i
        dta$b1[row_index] <- b_1i
        dta$b2[row_index] <- b_2i
        dta$b3[row_index] <- b_3i
        
        row_index_lag1 <- row_index
    }
    
    return(dta)
}


##### example: use of lmer() #####
if( 0 ){
    sample_size <- 1000
    total_T <- 20
    
    dta <- dgm_with_treatment(sample_size, total_T, dgm_type = 1)
    summary(dta)
    # dta$A <- dta$A - dta$prob_A # action centering doesn't matter when prob_A is constant
    
    fit <- lmer(Y ~ X * A + (1 + A | userid), data = dta)
    
    fit <- lmer(Y ~ X * A + (X * A | userid), data = dta)
    fit
    
    summary(fit)$coefficients
    
    attr(summary(fit)$varcor$userid, "stddev") # estimated standard deviation of random effect
}
