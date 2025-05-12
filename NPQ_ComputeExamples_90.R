
################################################
###	Basic settings => user-input required
################################################

###
###	Set work directory and clear workspace
###

#	setwd("<your-directory-path-here>")

rm(list = ls())




###
###	Load data set and define naming and whether its a simulation 
###

#	name.modifier		<- ""
name.modifier		<- "90"
#	name.modifier	<- "-season-off"
#	name.modifier		<- "zhou-portnoy"; 	mod.tau	<- TRUE
#	name.modifier		<- "short-ts"

logarithmize		<- FALSE
#	logarithmize	<- TRUE



###	Alternative 1: Federal reserve industrial production monthly data set 

dat.input   <- read.csv(file = "./Data/IPG2211A2N.csv")      # Whole ts
dat.input	<- dat.input[dat.input$DATE >= "1973-01-01", ]
dat.input	<- ts(
    data = if(logarithmize){log(dat.input$IPG2211A2N)} else {dat.input$IPG2211A2N},
    start	= c(1973, 1),
    end	= c(2022, 12),
#	end	= c(1998, 12),
    frequency = 12
)
dat.name		<- paste("inpro73m", name.modifier, ifelse(test = logarithmize, yes = "-log", no = ""), sep = "")
#	dat.input	<- window(dat.input, end = c(1998,12))




###	Alternative 2: Federal reserve industrial production quarterly data set 

dat.input   <- read.csv(file = "./Data/IPG2211A2NQ.csv")      # Whole ts
dat.input	<- dat.input[dat.input$DATE >= "1973-01-01", ]
dat.input	<- ts(
    data = if(logarithmize){log(dat.input$IPG2211A2NQ)} else {dat.input$IPG2211A2NQ},
    start 	= c(1973, 1),
    end	= c(2022, 4),
    frequency = 4
)
dat.name		<- paste("inpro73q", name.modifier, ifelse(test = logarithmize, yes = "-log", no = ""), sep = "")



if(FALSE){
###	Alternative 1b: "passengers" data set
#
dat.input		<- ts(
  data		= read.csv2(file = "passengers.csv")$passengers, 
  start		= c(2000, 1), 
  frequency		= 12
)
dat.name		<- "passengers"
}









###
###	Prediction horizon, sample splitting proportions, and folder name for saving all results
###

h.set		<- 1:frequency(dat.input)
#	h.set		<- 1


train.perc	<- 70		# What percentage of observations should be used the estimate the models (NPS, ISS, SAR)
min.perc.bw	<- 40		# What minimal percentage of observations should be used to estimate the NPS-bandwidths
#	Note: train.perc=80 and min.perc.bw=60 yield a 60:20:20 split of the data for NPS, as is common in Machine Learning
#	cover.perc	<- 80		# Nominal prediction interval coverage (percentage)
extending.window		<- FALSE



tauL	<- 0.05
#	tauL	<- 0.1
tauU	<- 1 - tauL





###	BW grid (for trend just the equidistant part)

bt.seq.equi		<- seq(from = -1, to = 1, by = 0.5)
#	bt.seq.equi		<- seq(from = -1, to = 1, by = 0.1)
#	bs.seq		<- c(0.05, 1:3/4, 0.95)	
bs.seq		<- c(0.1, 1:3/4, 0.9)	
#	bs.seq		<- seq(from = 0, to = 1, by = 0.05)
#	bs.seq		<- 1:3/4

#	bt.seq.equi		<- seq(from = -1, to = 1, by = 0.1)
#	bs.seq		<- seq(from = 0.05, to = 0.95, by = 0.1)
#	bs.seq	<- 0:4/4

#	bs.seq	<- 0:4/20






################################################
###	Additional settings => no user-input required
################################################

###
###	Packages, functions to be sourced, names of functions to be sourced later, as well as folder for saving all results
###

#	install.packages("quantreg")
library(quantreg)

source("NPQ_Functions.R")		# Functions for automated forecasting

examplecontrol.file	<- "NPQ_Control.R"	# Functions used for point prediction computation and evaluation



day.of.computations	<- as.character(Sys.Date())

folder.name		<- paste("QResults_", min.perc.bw, "-", train.perc - min.perc.bw, "-", 100-train.perc, "_", dat.name, "_", day.of.computations, sep = "")



###
###	Determine time series length and sub-interval lengths
###

T	<- length(dat.input)

sample.range.start	<- floor(T*min.perc.bw/100)
insample.end		<- floor(T*train.perc/100)


if(FALSE){
if(mod.tau){
  alp			<- (tauL + (1 - tauU))/2
  delta.tauL	<- qnorm(p = 1 - alp)*sqrt(tauL*(1 - tauL))/sample.range.start
  delta.tauU	<- qnorm(p = 1 - alp)*sqrt(tauU*(1 - tauU))/sample.range.start
  tauL	<- tauL - delta.tauL
  tauU	<- tauU + delta.tauU
}
}




###
###	NPS configuration
###

###	Prespecify bandwidth vector or grid 






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
###	Start time measurement
###

time.computation.starts	<- Sys.time()












################################################
###	Computation of interval predictions
################################################

res.cover	<- matrix(data = NA, nrow = length(h.set), ncol = 2, dimnames = list(h.set, c("ISV", "ISVC")))


#	for(tau in tau.set){
#	tau	<- 0.5
for(h in h.set){



#	h	<- 1
#	h	<- frequency(dat.input)



#	sam.set	<- insample.end:(T - h)
#	for(sam in sam.set){
sam	<- insample.end

  dat		<- window(dat.input, end = time(dat.input)[sam])

  sim		<- FALSE

  source(examplecontrol.file)

  assign(
    x		= paste("res.h", h, sep = ""),
    value	= np.results
  )



#  print(paste("Computation for ", sam, " of ", max(sam.set), " done.", sep = ""))
#	}	# sam-loop
}	# h-loop
#	}	# tau-loop



library(forecast)

source("NPQ_ComputeISS.R")




save(list = c(ls(pattern = "res.h"), "iss.results"), file = paste(dat.name, ".RData", sep = ""))













################################################
###	Evaluation
################################################

###
###	Components of Tables 1 & 2 of paper
###

#	load(file = paste(dat.name, ".RData", sep = ""))
#	load(file = "inpro73m90.RData")

#	load(file = "inpro73q80.RData")
#	load(file = "inpro73m80.RData")


h	<- 1
#	h	<- 4
#	h	<- 12


res	<- get(paste("res.h", h, sep = ""))

#	names(res)

column.set	<- c("overpred.pct", "coverage", "underpred.pct", "overpred.pen", "avg.int.length", "underpred.pen", "intervalscore")



###	Training sample results

data.frame(cov.nominal = (tauU - tauL)*100, h = h, round(res$results.train[res$train.opt.rows, column.set], 2), crit = res$results.test[, "train.opt"])



###	Test sample results

data.frame(cov.nominal = (tauU - tauL)*100, h = h, round(res$results.test[, column.set], 2), crit = res$results.test[, "train.opt"])

data.frame(cov.nominal = (tauU - tauL)*100, h = h, t(round(iss.results[[paste("results.test.h", h, sep = "")]][column.set], 2)), method = "iss")






###
###	Figure 4 of paper
###

config.set	<- c("m90", "q90", "m80", "q80")


for(cf in config.set){

  S		<- 12*(substr(x = cf, start = 1, stop = 1) == "m") + 4*(substr(x = cf, start = 1, stop = 1) == "q")

  res		<- data.frame(array(data = NA, dim = c(S, 4), dimnames = list(NULL, c("h", "ISV", "ISVC", "ISS"))))
  res$h	<- 1:S


  load(file = paste("inpro73", cf, ".RData", sep = ""))

  for(h in 1:S){
    res[h, c("ISV", "ISVC")]	<- get(paste("res.h", h, sep = ""))$results.test$coverage
    res[h, "ISS"]			<- iss.results[[paste("results.test.h", h, sep = "")]]["coverage"]
  }

  assign(x = paste("res_", cf, sep = ""), value = res)
}




col.set	<- c("darkorange", "navy", "limegreen")

#	par(mfrow = c(2, 2), mgp = c(2, 1, 0), mai = c(0.6, 0.6, 0.1, 0.1))
for(cf in config.set){

#	pdf(file = paste(cf, "_Coverage.pdf", sep = ""), width = 5.5, height = 3.3)

par(mfrow = c(1, 1), mgp = c(2, 1, 0), mai = c(0.6, 0.6, 0.1, 0.1))

res	<- get(paste("res_", cf, sep = ""))
nomin	<- as.numeric(substr(x = cf, start = 2, stop = 3))


matplot(res, xaxt = "n", type = "n", ylim = c(45, 99), frame.plot = FALSE, xlab = "h", ylab = paste("Coverage (nominal: ", nomin, "%)", sep = ""))

axis(side = 1, at = res$h, labels = res$h)
abline(h = nomin, col = 1, lty = 2)
matlines(res[, -1], type = "b", pch = 19, col = col.set, lty = 1)
text(x = min(res$h), y = res[1, -1] + c(-2, +2), labels = c("IS", expression(IS[c]), "ISS"), col = col.set)

#	dev.off()
}








