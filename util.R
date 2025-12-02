# basic utility functions 
sampleVar <- function(data){
  return (var(data))
}
sampleSTd <- function(data){
  return(sd(data))
}
#-------------------------------

binomExact <- function(n, p , k){
  # binomExact(times, p of win, expected)
  return(dbinom(k,n,p))
}
binomLessThanEqual <-function(n, p , k){
  # binomLessThan(times, p of win, expected)
  return(pbinom(k,n,p))
}
binomLessThan <- function(n, p, k){
  # P(X < k) = P(X <= k-1)
  return(pbinom(k - 1, size = n, prob = p))
}

# P(X > k)
binomMoreThan <- function(n, p, k){
  # P(X > k) = 1 - P(X <= k)
  return(1 - pbinom(k, size = n, prob = p))
}
binomExact(10,0.17,3)#0.1599833 5.20

geomExact <- function(p, k){
  # geomExact(p of win, trial_number)
  # P(First success is exactly on trial k)
  if (k <= 0) {
    return(0)
  }
  return(dgeom(k - 1, prob = p))
}

geomLessThanEqual <- function(p, k){
  # geomLessThan(p of win, trial_number)
  #  P(First success is on or before trial k) = P(X <= k)
  return(pgeom(k - 1, prob = p))
}
# P(X < x)
geomLessThan_Failures <- function(p, x){
  # P(X < x) = P(X <= x-1)
  # Probability of less than x failures before the first success.
  return(pgeom(x - 1, prob = p))
}

# P(X > x)
geomMoreThan_Failures <- function(p, x){
  # P(X > x) = 1 - P(X <= x)
  # Probability of strictly more than x failures before the first success.
  return(1 - pgeom(x, prob = p))
}

# poisson 
poisExact <- function(lambda, k){
  # poisExact(lambda, expected_count)
  # Calculates the probability of getting exactly k events.
  return(dpois(k, lambda = lambda))
}
poisLessThanEqual <-function(lambda, k){
  # poisLessThanEqual(lambda, expected_count)
  # Calculates the probability of getting k or fewer events.
  return(ppois(k, lambda = lambda))
}
poisLessThan <- function(lambda, k){
  # P(X<k)
  return(poisLessThanEqual(lambda,k)-poisExact(lambda,k))
}
poisMoreThan <- function(lambda,k){
  #P(x>k)
  return(1-poisLessThanEqual(lambda,k))
}
poisExact(4,5)
# fewwer than 2 admissions 
poisLessThan(4,2)
poisMoreThan(4,1)

#----------------------------------
#continuos from now 
# P(X<k) = P(X<=k)
#------------------------------
unifLessThanEqual <- function(a, b, k){
  # Calculates P(X <= k)
  return(punif(k, min = a, max = b))
}

unifLessThan <- function(a, b, k){
  # Calculates P(X < k) - Same as P(X <= k)
  return(punif(k, min = a, max = b))
}

# P(X > k)
unifMoreThan <- function(a, b, k){
  # Calculates P(X > k) = 1 - P(X <= k)
  return(1 - punif(k, min = a, max = b))
}
unifRangeInclusive <- function(a, b, x1, x2){
  # Calculates P(x1 <= X <= x2) = P(X <= x2) - P(X <= x1)
  # P(X <= x2)
  prob_up_to_x2 <- punif(x2, min = a, max = b)
  
  # P(X <= x1)
  prob_up_to_x1 <- punif(x1, min = a, max = b)
  
  return(prob_up_to_x2 - prob_up_to_x1)
}
unifRangeExclusive <- function(a, b, x1, x2){
  # Calculates P(x1 < X < x2) - Same as P(x1 <= X <= x2)
  return(unifRangeInclusive(a, b, x1, x2))
}
unifMean <- function(a, b){
  # Calculates the mean (expected value) E[X]
  return((a + b) / 2)
}

unifVariance <- function(a, b){
  # Calculates the variance Var[X]
  return(((b - a)^2) / 12)
}

#--------------------------------------------------------
#Normal distr 
#------------------------------------------------

normLessThanEqual <- function(mu, sigma, k){
  # Calculates P(X <= k)
  return(pnorm(k, mean = mu, sd = sigma))
}
normLessThan <- function(mu, sigma, k){
  # Calculates P(X < k) - Same as P(X <= k)
  return(pnorm(k, mean = mu, sd = sigma))
}

# P(X > k)
normMoreThan <- function(mu, sigma, k){
  # Calculates P(X > k) = 1 - P(X <= k)
  return(1 - pnorm(k, mean = mu, sd = sigma))
}

# P(x1 <= X <= x2) and P(x1 < X < x2) - Identical for Continuous Distributions
normRangeInclusive <- function(mu, sigma, x1, x2){
  # Calculates P(x1 <= X <= x2) = P(X <= x2) - P(X <= x1)
  
  if (x1 > x2) {
    return(0)
  }
  
  # P(X <= x2)
  prob_up_to_x2 <- pnorm(x2, mean = mu, sd = sigma)
  
  # P(X <= x1)
  prob_up_to_x1 <- pnorm(x1, mean = mu, sd = sigma)
  
  return(prob_up_to_x2 - prob_up_to_x1)
}

normRangeExclusive <- function(mu, sigma, x1, x2){
  # Calculates P(x1 < X < x2) - Same as P(x1 <= X <= x2)
  return(normRangeInclusive(mu, sigma, x1, x2))
}
normMean <- function(mu, sigma){
  # Calculates the mean (expected value) E[X]
  return(mu)
}

normVariance <- function(mu, sigma){
  # Calculates the variance Var[X]
  # Note: The input parameter is conventionally the standard deviation (sigma)
  return(sigma^2)
}

#-------------------------------------------
# in case of normalized use functions with
# mu = 0, sigma = 1 , Z = (x_mu)/sigma
#---------------------------------------------
normLessThan(6,1.5,5) # check this, here I got a different result 

# critical value and alpha 
# tail: A string indicating the tail ("upper" for P(X > k) = alpha, 
#       "lower" for P(X < k) = alpha). Defaults to "upper".



#-------------------------- 
# approximatiosn left, will we need it in the exam ? 
# -----------------------------

#----------------------------------
# Chi-Squared 
#-----------------------------------

chisqLessThanEqual <- function(df, k){
  # Calculates P(X <= k)
  return(pchisq(k, df = df))
}

chisqLessThan <- function(df, k){
  # Calculates P(X < k) - Same as P(X <= k)
  return(pchisq(k, df = df))
}

# P(X > k)
chisqMoreThan <- function(df, k){
  # Calculates P(X > k) = 1 - P(X <= k). This is often the P-value in tests.
  return(1 - pchisq(k, df = df))
}

# P(x1 <= X <= x2) and P(x1 < X < x2) - Identical for Continuous Distributions
chisqRangeInclusive <- function(df, x1, x2){
  # Calculates P(x1 <= X <= x2) = P(X <= x2) - P(X <= x1)
  
  if (x1 > x2) {
    return(0)
  }
  
  # P(X <= x2)
  prob_up_to_x2 <- pchisq(x2, df = df)
  
  # P(X <= x1)
  prob_up_to_x1 <- pchisq(x1, df = df)
  
  return(prob_up_to_x2 - prob_up_to_x1)
}

chisqRangeExclusive <- function(df, x1, x2){
  # Calculates P(x1 < X < x2) - Same as P(x1 <= X <= x2)
  return(chisqRangeInclusive(df, x1, x2))
}
chisqMean <- function(df){
  # Calculates the mean E[X]
  return(df)
}

chisqVariance <- function(df){
  # Calculates the variance Var[X]
  return(2 * df)
}

# 
# 
# 
# P(X <= k) or P(X < k)
tLessThanEqual <- function(df, k){
  # Calculates P(X <= k)
  return(pt(k, df = df))
}

# P(X > k)
tMoreThan <- function(df, k){
  # Calculates P(X > k) = 1 - P(X <= k)
  return(1 - pt(k, df = df))
}

# P(x1 <= X <= x2)
tRangeInclusive <- function(df, x1, x2){
  # Calculates P(x1 <= X <= x2) = P(X <= x2) - P(X <= x1)
  
  if (x1 > x2) {
    return(0)
  }
  return(pt(x2, df = df) - pt(x1, df = df))
}

# Mean and Variance
# Note: Mean is 0 if df > 1. Variance is df/(df-2) if df > 2.
tMean <- function(df){
  if (df > 1) {
    return(0)
  } else {
    return(NaN) # Undefined for df <= 1
  }
}

tVariance <- function(df){
  if (df > 2) {
    return(df / (df - 2))
  } else if (df > 0 && df <= 2) {
    return(Inf) # Infinite for 1 < df <= 2
  } else {
    return(NaN) # Undefined for df <= 0
  }
}

# 
# 
# 

# P(X <= k) or P(X < k)
fLessThanEqual <- function(df1, df2, k){
  # Calculates P(X <= k)
  return(pf(k, df1 = df1, df2 = df2))
}

# P(X > k)
fMoreThan <- function(df1, df2, k){
  # Calculates P(X > k) = 1 - P(X <= k). This is often the P-value in ANOVA.
  return(1 - pf(k, df1 = df1, df2 = df2))
}

# P(x1 <= X <= x2)
fRangeInclusive <- function(df1, df2, x1, x2){
  # Calculates P(x1 <= X <= x2) = P(X <= x2) - P(X <= x1)
  
  if (x1 > x2) {
    return(0)
  }
  
  prob_up_to_x2 <- pf(x2, df1 = df1, df2 = df2)
  prob_up_to_x1 <- pf(x1, df1 = df1, df2 = df2)
  
  return(prob_up_to_x2 - prob_up_to_x1)
}

# Mean and Variance
# Note: Mean is df2/(df2-2) if df2 > 2.

fMean <- function(df1, df2){
  if (df2 > 2) {
    return(df2 / (df2 - 2))
  } else {
    return(NaN) # Undefined for df2 <= 2
  }
}

fVariance <- function(df1, df2){
  if (df2 > 4) {
    numerator <- 2 * df2^2 * (df1 + df2 - 2)
    denominator <- df1 * (df2 - 2)^2 * (df2 - 4)
    return(numerator / denominator)
  } else {
    return(NaN) # Undefined for df2 <= 4
  }
}


############################### 
########################################
###### chapter6
#############################
########################


# distributions 
# 1. Mean : variance YES, n>=30
normSE <- function(sigma, n){
  # SE = σ / sqrt(n)
  return(sigma / sqrt(n))
}

normZStatistic <- function(x_bar, mu, sigma, n){
  # Z = (X̄ - μ) / SE
  se <- normSE(sigma, n)
  return((x_bar - mu) / se)
}

criticalZValue <- function(alpha){
  # Returns the positive Z-score (z_crit) for a two-tailed test, 
  # where P(|Z| > z_crit) = alpha.
  # Use the Standard Normal Distribution (mu=0, sigma=1)
  return(qnorm(1 - alpha/2, mean = 0, sd = 1))
}

# Mean: var NO n<30: usign tDist with n-1 df 

tEstSE <- function(S, n){
  # Est. SE = S / sqrt(n)
  return(S / sqrt(n))
}

# 2. Calculates the T-Statistic
tStatistic <- function(x_bar, mu, S, n){
  # t = (X̄ - μ) / Est. SE
  est_se <- tEstSE(S, n)
  return((x_bar - mu) / est_se)
}
# the exercise at 17 ?
denom <- tEstSE(10.85,16)
t1 <- 8/denom
tRangeInclusive(15,-t1,t1)

#########3
## also other distributions : I left them 
###########3

#-----------------------------
# confidence intervals 
#-----------------------------
intervalMean <- function(x_bar, n, alpha, sigma = NULL, S = NULL) {
  
  # --- Input Validation ---
  if (is.null(sigma) && is.null(S)) {
    stop("Error: Must provide either the population standard deviation (sigma) or the sample standard deviation (S).")
  }
  
  conf_level <- 1 - alpha
  df <- n - 1
  
  # --- Determine Case ---
  
  # CASE 1: Population Variance KNOWN (Always Z-Interval)
  if (!is.null(sigma)) {
    
    # 1. Known Variance: I = [ X̄ ± z_α/2 * (σ / √n) ]
    
    # Standard Error (SE)
    se <- sigma / sqrt(n)
    
    # Critical Value (Z-score)
    crit_value <- qnorm(1 - alpha / 2) 
    method <- "Z-Interval (Sigma Known)"
    
  } 
  
  # CASE 2 & 3: Population Variance UNKNOWN (Use S)
  else if (!is.null(S)) {
    
    # Estimated Standard Error (Est. SE)
    se <- S / sqrt(n)
    
    # Sub-Case 2: Large Sample (n > 30) -> Z-approximation (from your image)
    if (n > 30) {
      crit_value <- qnorm(1 - alpha / 2) 
      method <- "Z-Interval (Sigma Unknown, n > 30)"
      
    } 
    
    # Sub-Case 3: Small Sample (n <= 30) -> T-Interval
    else {
      # 2. Unknown Variance, Small Sample: I = [ X̄ ± t_α/2, n-1 * (S / √n) ]
      crit_value <- qt(1 - alpha / 2, df = df)
      method <- paste0("T-Interval (Sigma Unknown, n <= 30, df=", df, ")")
    }
  }
  
  # --- Calculate Interval ---
  
  # Margin of Error (ME)
  me <- crit_value * se
  
  lower_bound <- x_bar - me
  upper_bound <- x_bar + me
  
  # --- Return Result ---
  result <- list(
    Method = method,
    Confidence_Level = conf_level,
    Sample_Mean = x_bar,
    Margin_of_Error = me,
    Interval = c(Lower = lower_bound, Upper = upper_bound)
  )
  
  return(result)
}
# Calculates the Confidence Interval for the Population Variance (σ^2).
# Parameters: S2 (Sample Variance S^2), n (Sample Size), alpha (Significance Level)
intervalVariance <- function(S2, n, alpha){
  df <- n - 1
  
  # 1. Find Critical Chi-Squared Values
  # Lower Critical Value (Area to the left = alpha/2)
  chisq_lower_crit <- qchisq(alpha / 2, df = df)
  
  # Upper Critical Value (Area to the left = 1 - alpha/2)
  chisq_upper_crit <- qchisq(1 - alpha / 2, df = df)
  
  # 2. Calculate Numerator (n-1) * S^2
  numerator <- df * S2
  
  # 3. Calculate Interval Bounds (Note the inversion in the denominator)
  # Lower Bound uses the Upper Critical Value:
  lower_bound <- numerator / chisq_upper_crit
  
  # Upper Bound uses the Lower Critical Value:
  upper_bound <- numerator / chisq_lower_crit
  
  return(c(Lower = lower_bound, Upper = upper_bound))
}

# Calculates the Confidence Interval for the Population Proportion (p).
# Parameters: p_hat (Sample Proportion), n (Sample Size), alpha (Significance Level)
intervalProportion <- function(p_hat, n, alpha){
  
  # 1. Calculate Estimated Standard Error of Proportion (Est. SE_phat)
  # This uses p_hat as the estimate for p in the SE formula.
  est_se_phat <- sqrt(p_hat * (1 - p_hat) / n)
  
  # 2. Find Critical Z-Value (z_alpha/2)
  # qnorm uses the area to the left, which is 1 - alpha/2.
  z_crit <- qnorm(1 - alpha / 2, mean = 0, sd = 1)
  
  # 3. Calculate Margin of Error (ME)
  me <- z_crit * est_se_phat
  
  # 4. Calculate Interval Bounds
  lower_bound <- p_hat - me
  upper_bound <- p_hat + me
  
  return(c(Lower = lower_bound, Upper = upper_bound))
}



data <- c (2.2 ,2.66 ,2.74 ,3.41 ,2.46 ,2.96 ,3.34 ,
           2.16 ,2.46 ,2.71 ,2.04 ,3.74 ,3.24 ,3.92 ,2.38 ,
           2.82 ,2.2 , 2.42 ,2.82 ,2.84 ,4.22 ,3.64 ,1.77 ,
           3.44 ,1.53)
mn <- mean(data)
n <- length(data)
S <- sampleSTd(data)
intervalMean(mn,n,0.05,NULL,S)
# Attention: for mean intercval, the other that is NULL must input NULL 


####-----------------------------------
# confidence interval for 2 samples: 
# problem here 
###------------------------------------


# CI for the difference between 2 means 
intervalDiffMean <- function(mean1, mean2, n1, n2, alpha, 
                             sigma1 = NULL, sigma2 = NULL, 
                             S1 = NULL, S2 = NULL, 
                             equal_variance = FALSE) {
  
  # --- Setup ---
  diff_mean <- mean1 - mean2
  conf_level <- 1 - alpha
  
  # --- Input Validation ---
  if (is.null(sigma1) && is.null(sigma2) && is.null(S1) && is.null(S2)) {
    stop("Error: Must provide EITHER population SDs (sigma1, sigma2) OR sample SDs (S1, S2).")
  }
  
  # --- CASE 1: Population Variances KNOWN (Always Z-Interval) ---
  if (!is.null(sigma1) && !is.null(sigma2)) {
    
    # 1. Standard Error of the Difference (SE_diff)
    se_diff <- sqrt((sigma1^2 / n1) + (sigma2^2 / n2))
    
    # 2. Critical Value (Z-score)
    crit_value <- qnorm(1 - alpha / 2) 
    df <- Inf 
    method <- "Z-Interval (Variances Known)"
    
  } 
  
  # --- CASES 2-4: Population Variances UNKNOWN (Use S1, S2) ---
  else if (!is.null(S1) && !is.null(S2)) {
    
    # Determine the Standard Error and Degrees of Freedom (df) based on pooling
    
    if (equal_variance) {
      # CASE 2: Variances UNKNOWN but ASSUMED EQUAL (Pooled T-Interval)
      df <- n1 + n2 - 2
      
      # 1a. Calculate Pooled Sample Variance (S_p^2)
      Sp2 <- ((n1 - 1) * S1^2 + (n2 - 1) * S2^2) / df
      
      # 1b. Standard Error of the Difference (SE_diff) - Pooled
      se_diff <- sqrt(Sp2 * (1/n1 + 1/n2))
      
      method <- "T-Interval (Variances Unknown, Pooled)"
      
    } else {
      # CASE 3/4: Variances UNKNOWN and ASSUMED UNEQUAL (Unpooled T/Z-Interval)
      
      # 1. Standard Error of the Difference (SE_diff) - Unpooled
      se_diff <- sqrt((S1^2 / n1) + (S2^2 / n2))
      
      # 2. Degrees of Freedom (Welch-Satterthwaite Approximation)
      numerator_welch <- (S1^2/n1 + S2^2/n2)^2
      denominator_welch <- (S1^2/n1)^2/(n1 + 1) + (S2^2/n2)^2/(n2 + 1)
      df <- floor(numerator_welch / denominator_welch)-2 # Round down to nearest integer
      
      method <- "T-Interval (Variances Unknown, Unpooled/Welch)"
    }
    
    # Determine Critical Value (T or Z, using the determined df)
    if (n1 > 30 && n2 > 30) {
      crit_value <- qnorm(1 - alpha / 2) 
      method <- paste(method, "(Approx. Z)")
    } else {
      crit_value <- qt(1 - alpha / 2, df = df) 
    }
    
  } else {
    stop("Fatal Error: Logic failure in checking sigma/S parameters.")
  }
  
  # --- Calculate Interval ---
  
  # Margin of Error (ME)
  me <- crit_value * se_diff
  
  lower_bound <- diff_mean - me
  upper_bound <- diff_mean + me
  
  # --- Return Result ---
  result <- list(
    Method = method,
    Confidence_Level = conf_level,
    Degrees_of_Freedom = df,
    Difference_of_Means = diff_mean,
    Margin_of_Error = me,
    Interval = c(Lower = lower_bound, Upper = upper_bound)
  )
  
  return(result)
}

data1 <- c(39.3,39.7,39.9,40,38.9,38.9,39.7,39.3)
data2 <- c(38.1,38,38.3,37.9,38.7,38.3,37)
mean1 <- mean(data1)
mean2 <- mean(data2)
n1 <- length(data1)
n2 <- length(data2)
S1 <- sampleSTd(data1)
S2 <- sampleSTd(data2)
intervalDiffMean(mean1,mean2,n1,n2,0.1,NULL,NULL,S1,S2)
intervalDiffMean(mean2,mean1,n2,n1,0.1,NULL,NULL,S2,S1)

############################################
###Hypothesis Testing 
#######################################

