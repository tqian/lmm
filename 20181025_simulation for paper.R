
# 2018.10.22 Tianchen Qian

# simulation study for the random effects model paper

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it

rm(list = ls())
set.seed(123)

library(rootSolve)
library(geepack)
library(lme4)
library(mvtnorm)
library(foreach)
library(doMC)
library(doRNG)

dgm_with_treatment <- function(sample_size, total_T, dgm_type) {
    
    if (dgm_type == 1 | dgm_type == 3) {
        alpha_0 <- - 1
        alpha_1 <- - 0.3
        beta_0 <- 0.5
        beta_1 <- 0.1
        sigma_b0 <- 2
        sigma_b1 <- 0
        sigma_b2 <- 1
        sigma_b3 <- 0
        sigma_eps <- 1
    } else if (dgm_type == 2) {
        alpha_0 <- - 1
        alpha_1 <- - 0.3
        beta_0 <- 0.5
        beta_1 <- 0.1
        sigma_b0 <- 2
        sigma_b1 <- 0.5
        sigma_b2 <- 1
        sigma_b3 <- 0.5
        sigma_eps <- 1
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
                dta$X[row_index] <- dta$Y[row_index_lag1]
            }
            dta$prob_A[row_index] <- rep(prob_a, sample_size)
        } else if (dgm_type == 2) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- ifelse(dta$X[row_index] > -0.44, 0.7, 0.3)
        } else if (dgm_type == 3) {
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
    
    dta <- dgm_with_treatment(sample_size, total_T, dgm_type = 2)
    summary(dta)
    # dta$A <- dta$A - dta$prob_A # action centering doesn't matter when prob_A is constant
    
    fit <- lmer(Y ~ X * A + (1 + A | userid), data = dta)
    fit
    
    fit <- lmer(Y ~ Z + X * A + (1 + A | userid), data = dta)
    fit
    
    summary(fit)$coefficients
    
    re_sd <- attr(summary(fit)$varcor$userid, "stddev")
}





##### simulation: using lmer() package #####

# if(0) {
    
max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))

set.seed(120)


nsim <- 1000

design <- expand.grid(sample_size = c(100, 500, 1000), total_T = c(10), dgm_type = 1:3)

for (idesign in 1:nrow(design)) {
    sample_size <- design$sample_size[idesign]
    total_T <- design$total_T[idesign]
    dgm_type <- design$dgm_type[idesign]
    
    # standard structure of the output for error simulations
    dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
    if (dgm_type == 1 | dgm_type == 3) {
        fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
    } else if (dgm_type == 2) {
        fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
    }
    coef <- summary(fit)$coefficients
    varcor <- summary(fit)$varcor
    for (irow in 1:nrow(coef)) {
        for (icol in 1:ncol(coef)) {
            coef[irow, icol] <- NA
        }
    }
    for (irow in 1:nrow(varcor$userid)) {
        for (icol in 1:ncol(varcor$userid)) {
            varcor$userid[irow, icol] <- NA
        }
    }
    attr(varcor, "sc") <- NA
    coef_NA_fill <- coef
    varcor_NA_fill <- varcor
    
    # start parallel jobs
    writeLines(c(""), "~/Downloads/log.txt")
    sink("~/Downloads/log.txt", append=FALSE)
    set.seed(123)
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
        if (isim %% 10 == 0) {
            cat(paste("Starting iteration",isim,"\n"))
        }
        dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
        

        
        solution <- tryCatch(
            {
                if (dgm_type == 1 | dgm_type == 3) {
                    fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
                } else if (dgm_type == 2) {
                    fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
                }
                list(coef = summary(fit)$coefficients, varcor = summary(fit)$varcor)
            },
            error = function(cond) {
                message("\nCatched error in multiroot inside weighted_centered_least_square():")
                message(cond)
                return(list(coef = coef_NA_fill, varcor = varcor_NA_fill))
            })

        output <- list(solution)
    }
    sink()
    
    saveRDS(result, file = paste0("20181025_simulation for paper_result/", idesign, ".RDS"))
}


design$alpha_0_bias <- design$alpha_0_sd <- design$alpha_0_cp <- 
    design$alpha_1_bias <- design$alpha_1_sd <- design$alpha_1_cp <- 
    design$beta_0_bias <- design$beta_0_sd <- design$beta_0_cp <- 
    design$beta_1_bias <- design$beta_1_sd <- design$beta_1_cp <- NA

for (idesign in 1:nrow(design)) {
    result <- readRDS(paste0("20181025_simulation for paper_result/", idesign, ".RDS"))
    
    dgm_type <- design$dgm_type[idesign]
    if (dgm_type == 1 | dgm_type == 3) {
        alpha_0_true <- - 1
        alpha_1_true <- - 0.3
        beta_0_true <- 0.5
        beta_1_true <- 0.1
        sigma_b0_true <- 2
        sigma_b1_true <- 0
        sigma_b2_true <- 1
        sigma_b3_true <- 0
        sigma_eps_true <- 1
    } else if (dgm_type == 2) {
        alpha_0_true <- - 1
        alpha_1_true <- - 0.3
        beta_0_true <- 0.5
        beta_1_true <- 0.1
        sigma_b0_true <- 2
        sigma_b1_true <- 0.5
        sigma_b2_true <- 1
        sigma_b3_true <- 0.5
        sigma_eps_true <- 1
    }
    
    alpha_0 <- sapply(result, function(l) l$coef["(Intercept)", "Estimate"])
    alpha_0_sd <- sapply(result, function(l) l$coef["(Intercept)", "Std. Error"])
    alpha_1 <- sapply(result, function(l) l$coef["X", "Estimate"])
    alpha_1_sd <- sapply(result, function(l) l$coef["X", "Std. Error"])
    beta_0 <- sapply(result, function(l) l$coef["A", "Estimate"])
    beta_0_sd <- sapply(result, function(l) l$coef["A", "Std. Error"])
    beta_1 <- sapply(result, function(l) l$coef["X:A", "Estimate"])
    beta_1_sd <- sapply(result, function(l) l$coef["X:A", "Std. Error"])
    
    varcor <- lapply(result, function(l) l$varcor$userid) # variance matrix for random effects; this is "G"
    
    sigma_eps <- sapply(result, function(l) attr(l$varcor, "sc"))
    
    design$alpha_0_bias[idesign] <- mean(alpha_0) - alpha_0_true
    design$alpha_0_sd[idesign] <- sd(alpha_0)
    design$alpha_0_cp[idesign] <- mean((alpha_0_true < alpha_0 + 1.96 * alpha_0_sd) & (alpha_0_true > alpha_0 - 1.96 * alpha_0_sd))
    
    design$alpha_1_bias[idesign] <- mean(alpha_1) - alpha_1_true
    design$alpha_1_sd[idesign] <- sd(alpha_1)
    design$alpha_1_cp[idesign] <- mean((alpha_1_true < alpha_1 + 1.96 * alpha_1_sd) & (alpha_1_true > alpha_1 - 1.96 * alpha_1_sd))
    
    design$beta_0_bias[idesign] <- mean(beta_0) - beta_0_true
    design$beta_0_sd[idesign] <- sd(beta_0)
    design$beta_0_cp[idesign] <- mean((beta_0_true < beta_0 + 1.96 * beta_0_sd) & (beta_0_true > beta_0 - 1.96 * beta_0_sd))
    
    design$beta_1_bias[idesign] <- mean(beta_1) - beta_1_true
    design$beta_1_sd[idesign] <- sd(beta_1)
    design$beta_1_cp[idesign] <- mean((beta_1_true < beta_1 + 1.96 * beta_1_sd) & (beta_1_true > beta_1 - 1.96 * beta_1_sd))
    
    print(idesign)
    print(apply(simplify2array(varcor), 1:2, mean))
    print(mean(sigma_eps))
}






mean(alpha_0)
mean(alpha_1)
mean(beta_0)
mean(beta_1)

apply(simplify2array(varcor), 1:2, mean)
mean(sigma_eps)

print("###############################")
print(paste0("Sample size: ", sample_size))
print(paste0("total_T: ", total_T))

mean(alpha_0) - alpha_0_true
sd(alpha_0)
mean((alpha_0_true < alpha_0 + 1.96 * alpha_0_sd) & (alpha_0_true > alpha_0 - 1.96 * alpha_0_sd))

mean(alpha_1) - alpha_1_true
sd(alpha_1)
mean((alpha_1_true < alpha_1 + 1.96 * alpha_1_sd) & (alpha_1_true > alpha_1 - 1.96 * alpha_1_sd))

mean(beta_0) - beta_0_true
sd(beta_0)
mean((beta_0_true < beta_0 + 1.96 * beta_0_sd) & (beta_0_true > beta_0 - 1.96 * beta_0_sd))

mean(beta_1) - beta_1_true
sd(beta_1)
mean((beta_1_true < beta_1 + 1.96 * beta_1_sd) & (beta_1_true > beta_1 - 1.96 * beta_1_sd))



mean(sigma_b) - sigma_b_true
sd(sigma_b)

