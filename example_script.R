rm(list=ls())

#Load in helper functions and call required packages 
source('Functions/helper_functions.R')
source('Functions/TailEQD.R')


# Example data set --------------------------------------------------------

#We load in an example sea level data set from Dauphin Island, Alabama
data = read.csv(file="ExampleDataset/dauphin_island_detren_declust.csv")

#Extract sea level
data = data$sea_level

# Candidate thresholds ----------------------------------------------------

#Specify non-exceedance quantile probabilities at which to define thresholds
threshold_probabilities = seq(0.5,0.999,by=0.001)

#Evaluate empricial quantiles from data
thresholds = unname(quantile(data,probs = threshold_probabilities))


# Specifying the baseline event ------------------------------------------

#Specify length of data (in years)
record_length = 40.42

#Specify baseline event return period. We recommend keeping this value fixed at 0.25 years
baseline_event_RP = 0.25 

#Compute the number of observations observed each year (on average)
points_per_year = length(data)/record_length 

#Compute non-exceedance probability corresponding to the baseline event
baseline_probability = 1-1/(baseline_event_RP*points_per_year)


# Selecting a threshold ---------------------------------------------------

thresh_select = thresh_select_function(data = data,thresh = thresholds,threshold_probabilities = threshold_probabilities,baseline_probability = baseline_probability,ppy=points_per_year,k=100)

plot(thresholds,thresh_select$dists,xlim=c(min(thresholds),max(thresholds[!is.na(thresh_select$dists)])),type="l",lwd=3,col="grey",xlab="Threshold",ylab="Tail EQD",main="Updated threshold selection",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)

abline(v=thresh_select$thresh,lwd=4,col=2)

plot(threshold_probabilities,thresh_select$dists,xlim=c(min(threshold_probabilities),max(threshold_probabilities[!is.na(thresh_select$dists)])),type="l",lwd=3,col="grey",xlab="Threshold",ylab="Tail EQD",main="Updated threshold selection",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)

abline(v=threshold_probabilities[which.min(thresh_select$dists)],lwd=4,col=2)

data_over_u = data[data > thresh_select$thresh]-thresh_select$thresh

m = length(data_over_u)

observed_quants = sort(data_over_u)

theoretical_quants = qgpd((1:m)/(m+1),shape=thresh_select$par[2],scale = thresh_select$par[1])

plot(theoretical_quants,observed_quants,xlim=range(theoretical_quants,observed_quants),
     ylim=range(theoretical_quants,observed_quants),pch=16,col=1,ylab="Empirical",xlab="Model",
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)
abline(a=0,b=1,lwd=3,col=2)
points(theoretical_quants,observed_quants,pch=16,col=1, cex=1)

# Parallelised code  ------------------------------------------------------

threshold_probabilities = seq(0.5,0.999,by=0.001)

thresholds = unname(quantile(data,probs = threshold_probabilities))

record_length = 40.42

baseline_event_RP = 0.25 

points_per_year = length(data)/record_length #average number of observations per year  

baseline_probability = 1-1/(baseline_event_RP*points_per_year) #probability of baseline event

start <- Sys.time()
thresh_select <- TailEQD(data = data,thresh = thresholds,threshold_probabilities = threshold_probabilities,baseline_probability = baseline_probability,ppy=points_per_year,k=200) 
end <- Sys.time()
unpara_time <- end - start 

start <- Sys.time()
thresh_select_para <- TailEQDpara(data = data,thresh = thresholds,threshold_probabilities = threshold_probabilities,baseline_probability = baseline_probability,ppy=points_per_year,k=200)
end <- Sys.time()
para_time <- end - start 

print(paste0("Parallelised code is ", format((as.numeric(unpara_time)/as.numeric(para_time)),digits=3)," times quicker than unparallised"))

#To demonstrate we get the same output, subject to bootstrap variability, we plot the distance functions against the quantile levels

plot(threshold_probabilities,thresh_select$dists,xlim=c(min(threshold_probabilities),max(threshold_probabilities[!is.na(thresh_select$dists)])),type="l",lwd=3,col="grey",xlab="Threshold",ylab="Tail EQD",main="Updated threshold selection",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
lines(threshold_probabilities,thresh_select_para$dists,lwd=3,col="navy")
legend("topleft", legend = c("Unparallelised", "Parallelised"), col = c("grey", "navy"), lty = 1,lwd=4,cex=1.3,bg="white")
abline(v=threshold_probabilities[which.min(thresh_select$dists)],lwd=4,col=2)
abline(v=threshold_probabilities[which.min(thresh_select_para$dists)],lwd=4,col="darkred")

#Extracting optimal threshold choices 

data_over_u = data[data > thresh_select_para$thresh]-thresh_select_para$thresh

m = length(data_over_u)

observed_quants = sort(data_over_u)

theoretical_quants = qgpd((1:m)/(m+1),shape=thresh_select_para$par[2],scale = thresh_select_para$par[1])

plot(theoretical_quants,observed_quants,xlim=range(theoretical_quants,observed_quants),
     ylim=range(theoretical_quants,observed_quants),pch=16,col=1,ylab="Empirical",xlab="Model",
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)
abline(a=0,b=1,lwd=3,col=2)
points(theoretical_quants,observed_quants,pch=16,col=1, cex=1)
