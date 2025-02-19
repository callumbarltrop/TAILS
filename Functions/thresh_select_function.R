#' Threshold selection method for univariate extremes
#'
#' 'thresh_select_function' selects a constant threshold which best captures the tail of the data via a Generalised Pareto distribution.
#'
#' @author Tom Collings et al. 
#'
#' @param data A numeric vector. The data for which to select a threshold.
#' @param thresh A numeric vector of proposed thresholds to test.
#' @param threshold_probabilities A numeric vector of the corresponding threshold probabilities. Can be null
#' @param baseline_probability A probability value above which we define the tail of the data.
#' @param ppy A positive value denoting the number of observations per year
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#' @param parallel A TRUE/FALSE value specifying whether or not to run the code in parallel. 
#' @param cores A positive integer denoting the number of cores for running in parallel. 
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.
#'
#' @examples
#' set.seed(12345)
#' data_test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
#' thresholds1 <- quantile(data_test1,seq(0,0.95,by=0.05))
#' baseline_probability <- 0.5
#' ppy <- 365
#' (example <- thresh_select_function(data = data_test, thresh = thresholds,baseline_probability = baseline_probability,ppy = ppy))

thresh_select_function <- function(data, thresh, threshold_probabilities=NULL, baseline_probability,ppy,k = 100, m = 500,parallel = FALSE,cores=NULL){

  # Check inputs are valid
  if (!is.numeric(data)) stop("data must be a vector")
  if (!is.numeric(thresh)) stop("thresh to be tested needs to be a vector")
  if (!is.numeric(threshold_probabilities) & !is.null(threshold_probabilities)) stop("threshold_probabilities need to be a vector or NULL")
  if (!is.numeric(baseline_probability) | length(baseline_probability)!=1 | baseline_probability < 0 | baseline_probability > 1) stop("baseline_probability must be a single value between 0 and 1")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples (k) must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities (m) must be a positive integer")
  if (!is.numeric(ppy) | length(ppy) != 1) stop("ppy must be a single numeric value")
  if (!is.logical(parallel)) stop("parallel can only take the values 'TRUE' or 'FALSE'")
  
  #length of data 
  n = length(data)
  
  if(is.null(threshold_probabilities)){
    ##### convert thresholds to empirical exceedance probabilities ######
    # Create the empirical CDF from the data
    ecdf_func <- ecdf(data)
    
    # Evaluate the empirical CDF at the selected quantiles
    threshold_probabilities <- ecdf_func(thresh)
  }
  
  #select probabilities of interest
  if (baseline_probability>1 - (10/n)) stop("Baseline event probability is too large - please decrease.")
  quantile_probs <- seq(baseline_probability,1 - (10/n),length=m) 
  
  if(parallel == FALSE){ #if parallel is set to false, we run in a for loop, which takes longer 
    
    meandistances <- xis <- sigmas <- num_excess <- numeric(length(thresh))
    
    for (i in 1:length(thresh)) { #for each threshold
      u <- thresh[i] #select threshold
      excess <- data[data > u] - u #exceedences over threshold
      num_excess[i] <- length(excess) #number of exceedences over threshold
      threshold_probability <- threshold_probabilities[i] #probability of the threshold compared to the empirical distribution
      
      if (threshold_probability <= 1 - (1/ppy) ) { #We only fit the model if the threshold is less than the 1 in 1 year empirical return level
        mle0 <- mean(excess)
        init.fit <- optim(GPD_LL, z = excess, par = c(mle0,0.1), control = list(fnscale = -1))
        xis[i] <- init.fit$par[[2]]
        sigmas[i] <- init.fit$par[[1]]
        distances <- numeric(k)
        
        #computing theoretical quantiles. This remains constant over bootstrapping 
        exceedance_probs <- 1 - (1-quantile_probs)/(1-threshold_probability)
        
        # mask to remove any exceedance probabilities that are less than 0 or equal to 1 (i.e., not valid)
        # note that this accounts for the cases when the baseline event is lower than the threshold - we just estimate quantiles above the threshold in t
        mask <- exceedance_probs > 0 & exceedance_probs < 1
        
        for (j in 1:k) { #for each bootstrap
          X <- sample(excess, num_excess[i], replace = TRUE) #sample with replacement from the exceedences
          mle <- mean(X)
          ifelse(xis[i] < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
          gpd.fit <- optim(GPD_LL, z = X, par = pars_init, control = list(fnscale = -1)) #fit the gpd to the sample
          
          ############ using matched quantiles (from callum's code)
          
          theoretical_quantiles <- qgpd(exceedance_probs[mask], scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
          empirical_quants <- quantile(X,probs=exceedance_probs[mask]) #to get the empirical quantiles for this iteration
          
          #find the errors between empirical and theoretical quantiles
          errors <- abs(empirical_quants - theoretical_quantiles) 
          
          #calculate mean errors - EQD
          distances[j] <- (1/length(errors)) * sum(errors) # errors over the threshold is inheriant in the quantile selection
        }
        
        meandistances[i] <- mean(distances) #find mean across all bootstraps
        
      } else{
        meandistances[i] <- NA
      }
    }
    
    chosen_index <- which.min(meandistances)
    chosen_threshold <- thresh[chosen_index]
    chosen_threshold_prob <- threshold_probabilities[chosen_index]
    xi <- xis[chosen_index]
    sigma <- sigmas[chosen_index]
    len <- num_excess[chosen_index]
    
  } else {
    
    if(is.null(cores)){
      cores = floor(0.5*detectCores())
    } else {
      if (!is.numeric(cores) | length(cores) != 1 | cores %% 1 != 0) stop("cores must be a single integer value, or NULL")
      if (cores > detectCores()-1) stop("cores is too high for your CPU setup. Please reduce")
    }
    
    cl <- makeCluster(cores)
    
    clusterExport(cl, list("data","GPD_LL","quantile_probs","k","m","qgpd","ppy"), envir = environment())
    
    TS_info <- parApply(cl, rbind(thresh,threshold_probabilities),2, function(xcol) {
      
      u <- xcol[1]
      excess <- data[data > u] - u #exceedences over threshold
      num_excess <- length(excess) #number of exceedences over threshold
      threshold_probability <- xcol[2] #probability of the threshold compared to the empirical distribution
      
      if (threshold_probability <= 1 - (1/ppy)) { #We only fit the model if the threshold is less than the 1 in 1 year empirical return level
        mle0 <- mean(excess)
        init.fit <- optim(GPD_LL, z = excess, par = c(mle0,0.1), control = list(fnscale = -1))
        xi <- init.fit$par[[2]]
        sigma <- init.fit$par[[1]]
        distances <- numeric(k)
        
        #computing theoretical quantiles. This remains constant over bootstrapping 
        exceedance_probs <- 1 - (1-quantile_probs)/(1-threshold_probability)
        
        # mask to remove any exceedance probabilities that are less than 0 or equal to 1 (i.e., not valid)
        # note that this accounts for the cases when the baseline event is lower than the threshold - we just estimate quantiles above the threshold in t
        mask <- exceedance_probs > 0 & exceedance_probs < 1
        
        for (j in 1:k) { #for each bootstrap
          X <- sample(excess, num_excess, replace = TRUE) #sample with replacement from the exceedences
          mle <- mean(X)
          ifelse(xi < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(sigma, xi) )
          gpd.fit <- optim(GPD_LL, z = X, par = pars_init, control = list(fnscale = -1)) #fit the gpd to the sample
          
          ############ using matched quantiles (from callum's code)
          
          theoretical_quantiles <- qgpd(exceedance_probs[mask], scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
          empirical_quants <- quantile(X,probs=exceedance_probs[mask]) #to get the empirical quantiles for this iteration
          
          #find the errors between empirical and theoretical quantiles
          errors <- abs(empirical_quants - theoretical_quantiles) 
          
          #calculate mean errors - EQD
          distances[j] <- (1/length(errors)) * sum(errors) # errors over the threshold is inheriant in the quantile selection
          
        }
        
        meandistance <- mean(distances) #find mean across all bootstraps
      } else{
        meandistance <- NA
        xi <- NA
        sigma <- NA
        num_excess <- NA
      }
      
      return(c(num_excess,sigma,xi,meandistance))
    })
    
    # Stop parallel processing
    stopCluster(cl)
    
    TS_info = t(TS_info)
    
    colnames(TS_info) = c("n_exc","sigmas","xis","dists")
    
    TS_info = as.data.frame(TS_info)
    
    meandistances <- TS_info$dists
    chosen_index <- which.min(meandistances)
    chosen_threshold <- thresh[chosen_index]
    chosen_threshold_prob <- threshold_probabilities[chosen_index]
    xi <- TS_info$xis[chosen_index]
    sigma <- TS_info$sigmas[chosen_index]
    len <- TS_info$n_exc[chosen_index]
    
  }

  result <- list(
    thresh = chosen_threshold, 
    threshold_prob = chosen_threshold_prob,
    par = c(sigma,xi), 
    num_excess = len, 
    dists = meandistances
   )
  
  return(result)

}
