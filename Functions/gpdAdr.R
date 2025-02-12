gpdAdGen <- function(n, theta) {
  data1 <- rgpd(n, shape = theta[2], scale = theta[1])
  mle <- mean(data1)
  ifelse(theta[2] < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(theta) )
  
  fit1 <- tryCatch(optim(GPD_LL, z = data1, par = pars_init, control = list(fnscale = -1)), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit1)) {
    teststat <- NA
  } else {
    scale1 <- fit1$par[[1]]
    shape1 <- fit1$par[[2]]
    newdata1 <- pgpd(data1, scale = scale1, shape = shape1)
    newdata1 <- sort(newdata1)
    i <- seq(1, n, 1)
    teststat <- n/2 + -2*sum(newdata1)-sum((2-(2*i-1)/n)*log(1-newdata1)) # right-sided AD test
  }
  teststat
}


#' Generalized Pareto Distribution Anderson-Darling Test
#'
#' Anderson-Darling goodness-of-fit test for the Generalized Pareto (GPD) distribution.
#' @param thresh_data Data should be in vector form, assumed to be from the GPD.
#' @param bootstrap Should bootstrap be used to obtain p-values for the test? By default, a table of critical values is used via interpolation. See details.
#' @param bootnum Number of replicates if bootstrap is used.
#' @param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @references Choulakian, V., & Stephens, M. A. (2001). Goodness-of-fit tests for the Generalized Pareto distribution. Technometrics, 43(4), 478-484.
#' @examples
#' ## Generate some data from GPD
#' x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
#' gpdAd(x)
#' @details Edited function from the eva package to get AD test stats. Edited to only include a right-sided AD test.
#' A table of critical values were generated via Monte Carlo simulation for shape
#' parameters -0.5 to 1.0 by 0.1, which provides p-values via log-linear interpolation from
#' .001 to .999. For p-values below .001, a linear equation exists by regressing -log(p-value)
#' on the critical values for the tail of the distribution (.950 to .999 upper percentiles). This
#' regression provides a method to extrapolate to arbitrarily small p-values.
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Estimated value of theta for the initial data.}
#' \item{effective_bootnum}{Effective number of bootstrap replicates if bootstrap
#' based p-value is used (only those that converged are used).}
#' @import parallel
#' @export
gpdAd <- function(thresh_data,paras, bootstrap = FALSE, bootnum = NULL, allowParallel = FALSE, numCores = 1) {
  
  if(bootstrap == TRUE & is.null(bootnum))
    stop("Must specify some number of boostrap samples")
  
  n <- length(thresh_data)
  
  scale <- paras[1]
  shape <- paras[2]
  theta <- c(scale, shape)
  
  if(bootstrap == FALSE & shape > 1)
    stop("Estimated parameters are outside the table range, please use the bootstrap version")
  
  newdata <- pgpd(thresh_data, scale = scale, shape = shape) #generate cdf
  
  newdata <- sort(newdata)
  
  i <- seq(1, n, 1)
  
  #do AD test
  stat <- n/2 + -2*sum(newdata)-sum((2-(2*i-1)/n)*log(1-newdata)) # right-sided AD test
  
  #get p value by computing AD test on bootstrap samples
  if(bootstrap == TRUE) {
    if(allowParallel == TRUE) {
      cl <- makeCluster(numCores)
      fun <- function(cl) {
        parSapply(cl, 1:bootnum, function(i,...) {gpdAdGen(n, theta)})
      }
      teststat <- fun(cl)
      stopCluster(cl)
    } else {
      teststat <- replicate(bootnum, gpdAdGen(n, theta))
    }
    teststat <- teststat[!is.na(teststat)]
    eff <- length(teststat)
    p <- (sum(teststat > stat) + 1) / (eff + 2)
    
  } 
  
  names(theta) <- c("Scale", "Shape")
  if(!bootstrap) {
    out <- list(as.numeric(stat), as.numeric(p), theta)
    names(out) <- c("statistic", "p.value", "theta")
  } else {
    out <- list(as.numeric(stat), as.numeric(p), theta, eff)
    names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  }
  out
}
