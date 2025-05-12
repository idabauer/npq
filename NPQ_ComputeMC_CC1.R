
################################################
###	Basic settings => user-input required
################################################

###
###	Set work directory and clear workspace
###

#	setwd("<your-directory-path-here>")

rm(list = ls())







###
###	Simulation basics and data
###

logarithmize	<- FALSE



R	<- 100


set.seed(1)


for(T in c(100, 200, 400)){
#	T	<- 200
#	T	<- 400
#	T	<- 100


dgp.file		<- "NPQ_MC_DGP-Cai-Chen.R"



###	DGP specification

dat.name		<- "DGP_CC1"	# DGP of Cai and Cheng with AR1 errors, lownoise (0.1), and trend

alpha.constant	<- FALSE	# Set to TRUE for no trend

eps.sd	<- 0.1		#	Amount of noise

e.AR1.coef	<- 0.9		#	AR1-errors vs. iid-errors





###	Compute realizations of Monte Carlo DGP

source(dgp.file)





###
###	Sample splitting proportions, and folder name for saving all results
###

train.perc	<- 90		# What percentage of observations should be used the estimate the models (NPS, ISS, SAR)
min.perc.bw	<- 60		# What minimal percentage of observations should be used to estimate the NPS-bandwidths
#	Note: train.perc=80 and min.perc.bw=60 yield a 60:20:20 split of the data for NPS, as is common in Machine Learning
#	cover.perc	<- 80		# Nominal prediction interval coverage (percentage)
extending.window		<- FALSE



#tauL	<- 0.05
tauL	<- 0.1
tauU	<- 1 - tauL



###	Determine sub-interval lengths

sample.range.start	<- floor(T*min.perc.bw/100)
insample.end		<- floor(T*train.perc/100)





###
###	NPS configuration
###

###	Prespecify bandwidth vector or grid 



###	BW grid (for trend just the equidistant part)

bt.seq.equi		<- seq(from = -1, to = 1, by = 0.5)

bs.seq		<- c(0.1, 1:3/4, 0.9)	
#	bs.seq		<- c(0.01, 0.05, 0.1, 0.15, 1:3/4, 0.9)	




###	NPS configuration remainders 

trend.year.set	<- c("year")			#	trend.year.set	<- c("year", "t")
#	trend.year.set	<- c("t")			#	trend.year.set	<- c("year", "t")


#	crit.set		<- c("AICc", "LSCV")		#	crit.set		<- c("MSEh", "AICc", "LSCV")
crit.set		<- c("all")
#	regtype.set		<- c("LL")				#	regtype.set		<- c("LL", "LC")
regtype.temp		<- "LL"
#	krnl.s.set		<- c("krnl.SEN")			#	krnl.s.set		<- c("krnl.RLO", "krnl.SEN", "krnl.SEL", "krnl.SEK", "krnl.SEG")
#	krnl.s.set		<- c("krnl.RLO", "krnl.SEN")
krnl.s.temp		<- "krnl.SEN"

trend.year		<- trend.year.set[1]






###
###	Packages, functions to be sourced, names of functions to be sourced later, as well as folder for saving all results
###

#	install.packages("quantreg")
library(quantreg)

#	install.packages("forecast")
library(forecast)


source("NPQ_Functions.R")		# Functions for automated forecasting

examplecontrol.file	<- "NPQ_ControlMC.R"	# Functions used for point prediction computation and evaluation



day.of.computations	<- as.character(Sys.Date())

#folder.name		<- paste("QResults_", min.perc.bw, "-", train.perc - min.perc.bw, "-", 100-train.perc, "_", dat.name, "_", day.of.computations, sep = "")
folder.name		<- paste("MC_", dat.name, "_", min.perc.bw, "-", train.perc - min.perc.bw, "-", 100-train.perc, "_", day.of.computations, sep = "")

if(!file.exists(folder.name)){dir.create(folder.name)}











################################################
###	Computation of the Simulation
################################################

###
###	Start time measurement
###

time.computation.starts	<- Sys.time()



###
###	Core of the simulation
###

for(r in 1:R){
#	for(r in 1:2){
  #	r	<- 1		# which replication

dat.input	<- ts(data = x_t[r, ], frequency = 4)

h.set		<- c(1, frequency(dat.input))
#	h.set		<- 1:frequency(dat.input)
#	h.set		<- 1







###
###	Computation of interval predictions
###

#	res.cover	<- matrix(data = NA, nrow = length(h.set), ncol = 2, dimnames = list(h.set, c("ISV", "ISVC")))



for(h in h.set){
#	h	<- 1
#	h	<- frequency(dat.input)

sam	<- insample.end

  dat		<- window(dat.input, end = time(dat.input)[sam])

#  sim		<- FALSE

  source(examplecontrol.file)

  assign(
    x		= paste("res.T", T, ".r", r, ".h", h, sep = ""),
    value	= np.results
  )

}	# h-loop
#	}	# tau-loop




source("NPQ_ComputeISSforMC.R")



print(paste("Run ", r, " of ", R, " done.", sep = ""))
}	# r-loop


save(list = c(ls(pattern = "res.T"), ls(pattern = "iss.T")), file = paste(folder.name, "/", dat.name, "_T", T, ".RData", sep = ""))



time.computation.ends	<- Sys.time()

difftime(time.computation.ends, time.computation.starts, unit = "mins")

}


















################################################
###	Evaluation
################################################

day.of.computations	<- as.character(Sys.Date())

dat.name		<- "DGP_CC1"
#	dat.name		<- "DGP_CC2"
#	dat.name		<- "DGP_CC3"

folder.name		<- paste("MC_", dat.name, "_", min.perc.bw, "-", train.perc - min.perc.bw, "-", 100-train.perc, "_", day.of.computations, sep = "")



###
###	MC
###

h.set	<- c(1, 4)



var.temp	<- "intervalscore"; xlimits		<- c(0, 10)
#	var.temp	<- "coverage"; xlimits		<- c(0, 100)



pdf(file = paste(folder.name, "_", var.temp, ".pdf", sep = ""), width = 5, height = 5)


for(h in h.set){

for(T in c(2, 4)*100){
#	for(T in 200){
#	T	<- 100
#	T	<- 200
#	T	<- 400

#	R

load(file = paste(folder.name, "/", dat.name, "_T", T, ".RData", sep = ""))



#for(h in h.set){
#	h	<- 1

  res	<- NULL

  for(r in 1:R){
  #	r	<- 1

    res	<- rbind(
      res,
      c(
        unlist(get(paste("res.T", T, ".r", r, ".h", h, sep = ""))[["results.test"]][var.temp]), 
        get(paste("iss.T", T, ".r", r, sep = ""))[[paste("results.test.h", h, sep = "")]][var.temp]
      )
    )
  }

  colnames(res)	<- c("ISV", "ISVC", "ISS")
  assign(x = paste("res.h", h, sep = ""), value = res)
#}

#var.temp
#apply(X = res.h1, MARGIN = 2, FUN = median)
#apply(X = res.h4, MARGIN = 2, FUN = median)



col.set	<- colorRampPalette(c("red", "darkorange1", "limegreen", "steelblue", "blue"))(3)



#	h	<- 1
#	h	<- 4



par(mfrow = c(1, 1))

res	<- get(paste("res.h", h, sep = ""))
plot(ecdf(res[, 1]), xlim = xlimits, col = col.set[1], main = paste(var.temp, " for h=", h, " and T=", T, sep = ""), xlab = var.temp, ylab = paste("F(", var.temp, ")", sep = ""), verticals = TRUE, pch = NA)
lines(ecdf(res[, 2]), col = col.set[2], verticals = TRUE, pch = NA)
lines(ecdf(res[, 3]), col = col.set[3], verticals = TRUE, pch = NA)
if(var.temp == "coverage"){abline(v = (tauU - tauL)*100)}
legend(x = "bottomright", legend = colnames(res), col = col.set, lty = 1, box.col = "white", bg = "white")
box()

}
}

dev.off()

