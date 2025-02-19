# Clear current workspace 
rm(list=ls())

# Load all required functions
source('Functions/helper_functions.R')
source('Functions/thresh_select_function.R')

# Checking for required packages. This function will install any required packages if they are not already available
packages = c("parallel") 
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Load in data from csv. Make sure your directory is set to the Github repo home.
data = read.csv(file="ExampleDataset/dauphin_island_detren_declust.csv")

# Extract sea level variable
data = data$sea_level

# Specify length of data (in years)
record_length = 40.42

# Specify baseline event return period. We recommend keeping this value fixed at 0.25 years
baseline_event_RP = 0.25 

# Compute the number of observations observed each year (on average)
points_per_year = length(data)/record_length 

# Compute non-exceedance probability corresponding to the baseline event
baseline_probability = 1-1/(baseline_event_RP*points_per_year)

# Specify non-exceedance quantile probabilities at which to define thresholds
threshold_probabilities = seq(0.5,0.999,by=0.001)

# Evaluate empirical quantiles from data
thresholds = unname(quantile(data,probs = threshold_probabilities))

# Apply the threshold selection algorithm to the data
thresh_select = thresh_select_function(data = data, # Specify the data vector
                                       thresh = thresholds, # Specify candidate thresholds 
                                       threshold_probabilities = threshold_probabilities, # Specify threshold non-exceedance probabilities. This can be left blank
                                       baseline_probability = baseline_probability, # Specify non-exceedance probability of baseline event
                                       ppy=points_per_year) # Specify the number of points per year

# Specify graphical plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

# Plot the distance metric over the candidate thresholds
plot(thresholds,thresh_select$dists,
     xlim=c(min(thresholds),max(thresholds[!is.na(thresh_select$dists)])),
     type="l",lwd=3,col="grey",
     xlab=expression(u),ylab=expression(tilde(d)(u)),
     main="Distance metric vs threshold",cex.lab=1.3, 
     cex.axis=1.5, cex.main=1.5)

# Highlight threshold choice 
abline(v=thresh_select$thresh,lwd=4,col="grey1")

# Plot the distance metric over the candidate threshold non-exceedance probabilities
plot(threshold_probabilities,thresh_select$dists,
     xlim=c(min(threshold_probabilities),max(threshold_probabilities[!is.na(thresh_select$dists)])),
     type="l",lwd=3,col="grey",
     xlab=expression(p),ylab=expression(tilde(d)(u)),
     main="Distance metric vs non-exceedance probability",cex.lab=1.3, 
     cex.axis=1.5, cex.main=1.5)

# Highlight non-exceedance probability choice 
abline(v=thresh_select$threshold_prob,lwd=4,col="grey1")

# Specify graphical plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

# Extract the exceedances of the chosen threshold
data_over_u = data[data > thresh_select$thresh]-thresh_select$thresh

# Compute number of exceedances
m = length(data_over_u)

# Compute observed order statistics 
observed_quants = sort(data_over_u)

# Compute corresponding theoretical (model) quantiles 
theoretical_quants = qgpd((1:m)/(m+1),shape=thresh_select$par[2],scale = thresh_select$par[1])

# Plot the observed and theoretical quantiles against eachother, check for good agreement 
plot(theoretical_quants,observed_quants,
     xlim=range(theoretical_quants,observed_quants),
     ylim=range(theoretical_quants,observed_quants),
     pch=16,
     col=1,
     ylab="Observed quantiles",
     xlab="Theoretical quantiles",
     main="GPD QQ plot",
     cex.lab=1.3, 
     cex.axis=1.2,
     cex.main=1.8, 
     cex=0.5)

abline(a=0,b=1,lwd=3,col=2)

points(theoretical_quants,observed_quants,pch=16,col=1, cex=1)

# Threshold select algorithm with different number of bootstraps and steps
thresh_select_diff_k_and_m = thresh_select_function(data = data, # Specify the data vector
                                       thresh = thresholds, # Specify candidate thresholds 
                                       threshold_probabilities = threshold_probabilities, # Specify threshold non-exceedance probabilities. This can be left blank
                                       baseline_probability = baseline_probability, # Specify non-exceedance probability of baseline event
                                       ppy=points_per_year, # Specify the number of points per year
                                       k = 200, # Specify the number of bootstrap iterations for the algorithm. Default is 100 
                                       m = 250) # Specify the number of quantile/probability levels over which to evaluate GPD fits. We DO NOT recommend changing this value 


# Specify graphical plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

# Plot the distance metric over the candidate thresholds
plot(thresholds,thresh_select$dists,xlim=c(min(thresholds),max(thresholds[!is.na(thresh_select$dists)])),type="l",lwd=3,col="grey",xlab=expression(u),ylab=expression(tilde(d)(u)),main="Distance metric vs threshold",cex.lab=1.3, cex.axis=1.5, cex.main=1.5)

# Add a line illustrating the same distances for different k and m values 
lines(thresholds,thresh_select_diff_k_and_m$dists,type="l",lwd=3,col="lightblue")

# Highlight threshold choice 
abline(v=thresh_select$thresh,lwd=4,col="grey1")

# Highlight threshold choice for different k and m values 
abline(v=thresh_select_diff_k_and_m$thresh,lwd=4,col=4)

# Plot the distance metric over the candidate threshold non-exceedance probabilities
plot(threshold_probabilities,thresh_select$dists,xlim=c(min(threshold_probabilities),max(threshold_probabilities[!is.na(thresh_select$dists)])),type="l",lwd=3,col="grey",xlab=expression(p),ylab=expression(tilde(d)(u)),main="Distance metric vs non-exceedance probability",cex.lab=1.3, cex.axis=1.5, cex.main=1.5)

# Add a line illustrating the same distances for different k and m values 
lines(threshold_probabilities,thresh_select_diff_k_and_m$dists,type="l",lwd=3,col="lightblue")

# Highlight non-exceedance probability choice 
abline(v=thresh_select$threshold_prob,lwd=4,col="grey1")

# Highlight non-exceedance probability choice for different k and m values 
abline(v=thresh_select_diff_k_and_m$threshold_prob,lwd=4,col=4)

legend("bottomleft",legend=c("k = 100, m = 500", "k = 200, m = 250"),lwd=4,col=c("grey","lightblue"),cex=1.5,bg="white")

# Threshold select algorithm in parallel
thresh_select_parallel = thresh_select_function(data = data, # Specify the data vector
                                                thresh = thresholds, # Specify candidate thresholds 
                                                threshold_probabilities = threshold_probabilities, # Specify threshold non-exceedance probabilities. This can be left blank
                                                baseline_probability = baseline_probability, # Specify non-exceedance probability of baseline event
                                                ppy=points_per_year, # Specify the number of points per year
                                                parallel = TRUE, # Specify that you wish for the algorithm to be run in parallel
                                                cores = 4) # Specify the number of cores for parallel computation
