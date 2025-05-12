
#################################################
###	Basic functions
#################################################

###
###	Function for determining evaluation sample
###

f.sample.range.set	<- function(
  dat,
  est.end.set	= (max(20, 4*frequency(dat))):(length(dat) - h),
  h	= 1,
  start.period = 1,
  extending	= FALSE
){
  if(extending){
    est.start.set		<- rep(x = start.period, times = length(est.end.set))
  } else {
    est.start.set		<- (start.period - 1) + 1:length(est.end.set)
  }
  assign(x = "sample.range.set", value = cbind(est.start.set, est.end.set), envir = .GlobalEnv)
}


#	f.sample.range.set(dat = dat, h = 12, extending = FALSE)
#	f.sample.range.set(dat = dat, h = 1, extending = TRUE)




###
###	Check function and ATWE for quantile regression evaluation
###

check		<- function(
  vec,
  tau.check
){
  (tau.check - as.numeric(vec < 0))*vec
}


ATWE	<- function(
  vec.obs,
  vec.est,
  tau.ATWE
){
  mean(check(vec = vec.obs - vec.est, tau.check = tau.ATWE))
}













#################################################
###	Functions for kernels, bandwidth estimation and nonparametric forecasting
#################################################

###
###	Established Weighting-/Kernel-Functions
###

###	Second order Gaussian kernel

krnl.GAU	<- function(
  x0,		# season/covariate position, where the regression function is estimated at
  x,		# observed season/covariate position
  b,		# bandwidth
  ...
){
  w		<- 1/b*dnorm(x = ((x0 - x)/b), mean = 0, sd = 1)		# weight
  return(w)
}



###	Unordered kernel of Racine and Li (2004)

krnl.RLU	<- function(
  x0,		# season/covariate position, where the regression function is estimated at
  x,		# observed season/covariate position
  b,		# bandwidth
  ...
){
  w		<- b^(x != x0)		# weight
  return(w)
}



###	Ordered kernel of Racine and Li (2004)

krnl.RLO	<- function(
  x0,		# season/covariate position, where the regression function is estimated at
  x,		# observed season/covariate position
  b,		# bandwidth
  ...
){
  d0		<- abs(x - x0)	# absolute distance
  w		<- b^d0		# weight
  return(w)
}

#	sapply(X = 1:12, FUN = krnl.LRO, x0 = 11, b = 1)




###
###	Seasonal kernel
###

###	Seasonal kernel based on RLO

krnl.SEN	<- function(	# former: SEA
  x0,		# season/covariate position, where the regression function is estimated at
  x,		# observed season/covariate position
  S,		# number of seasons
  b		# bandwidth
){
  d0		<- abs(x - x0)		# absolute distance
  d0.star	<- d0 - (d0 > (S/2))*(2*d0 - S)	# modified absolute distance
  w		<- b^d0.star		# weight
  return(w)
}











###
###	Weight all rows of ndat w.r.t. b
###


weight.x	<- function(
  x0.temp,			#	= ndat[ro, "s"],
  x.f0	= ndat$s,
  S.f0	= max(ndat$s),
  krnl.f0	= krnl.SEN,
  b.f0	= bw.grid["bs"]
){
  sapply(
    X		= x.f0, 
    FUN	= krnl.f0, 
    x0	= x0.temp, 
    S		= S.f0, 
    b		= b.f0
  )
}







###
###	Functions for future-directed regression
###

fdr.LL	<- function(
  frc.row,
  sample.range.set.f1	= sample.range.set,
  h.f1			= 1,	# h-step ahead forecast
  W.f1			= W,
  ndat.f1			= ndat,
  tau.f1			= tau,
  ...
){
  est.T.start	<- sample.range.set.f1[frc.row, "est.start.set"]
  est.T.end		<- sample.range.set.f1[frc.row, "est.end.set"]

  W.temp	<- W.f1[est.T.start:est.T.end, est.T.end + h.f1]
  if(sum(W.temp) != 0){
#    W0	<- diag(W.temp)
#    X0	<- as.matrix(data.frame(
#      int		= rep(1, times = length(W.temp)), 
#      trnd		= ndat.f1$t[est.T.start:est.T.end] - ndat.f1$t[est.T.end + h.f1]		# only local linear w.r.t. trend
#    ))

y.hat.LL	<- tryCatch(
  predict(
    object	= rq(
      y ~ t,
      tau		= tau.f1,
      data		= ndat.f1[est.T.start:est.T.end, ],
      weights	= W.temp, 
      method	= "fn"
    ), 
    newdata	= ndat.f1[est.T.end + h.f1, ]
  ),
  error = function(e){NA}
)
  } else {
    y.hat.LL	<- NA
  }
  return(y.hat.LL)
}







###
###	Objective functions for bandwidth estimation
###

np.reg	<- function(
  ro,
  bw.grid.f3		= bw.grid,
  ndat.f3			= ndat,
  tauL.f3			= tauL,
  tauU.f3			= tauU,
  sample.range.set.f3	= sample.range.set,
#  crit.f3			= "all",	# one of c("all", "QEV", "ISV")
  regtype.f3		= "LL",	# one of c("LC", "LL")
  h.f3			= 1,
  krnl.s.f3			= krnl.SEN,
  krnl.t.f3			= krnl.GAU,
  show.progress.f3	= FALSE,
  ...
){

#	ro	<- 105
  b	<- bw.grid.f3[ro, ]
  W.s	<- sapply(
    X		= ndat.f3$s, 
    FUN	= weight.x,
    x.f0	= ndat.f3$s,
    S.f0	= max(ndat.f3$s),
    krnl.f0	= krnl.s.f3,
    b.f0	= as.numeric(b["bs"])
  )
  W.t	<- sapply(
    X		= ndat.f3$t, 
    FUN	= weight.x,
    x.f0	= ndat.f3$t,
    S.f0	= max(ndat.f3$t),
    krnl.f0	= krnl.t.f3,
    b.f0	= as.numeric(b["bt"])
  )
  W	<- W.s*W.t

if(identical(W, diag(diag(W)))){	# computability
  obj.value		<- NA
} else {

  if(regtype.f3 == "LL"){

##    y.obs		<- ndat.f3[sample.range.set.f3[, "est.end.set"] + h.f3, "y"]
    y.frcL		<- sapply(X = 1:nrow(sample.range.set.f3), FUN = fdr.LL, sample.range.set.f1 = sample.range.set.f3, h.f1 = h.f3, W.f1 = W, ndat.f1 = ndat.f3, tau.f1 = tauL.f3)
    y.frcU		<- sapply(X = 1:nrow(sample.range.set.f3), FUN = fdr.LL, sample.range.set.f1 = sample.range.set.f3, h.f1 = h.f3, W.f1 = W, ndat.f1 = ndat.f3, tau.f1 = tauU.f3)






##    avg.int.length	<- mean(y.frcU - y.frcL)
##    pen.obsL	<- if(any(y.obs < y.frcL)){mean((1/tauL.f3)*(y.frcL - y.obs)[y.obs < y.frcL])} else {0}
##    pen.obsU	<- if(any(y.obs > y.frcU)){mean((1/(1 - tauU.f3))*(y.obs - y.frcU)[y.obs > y.frcU])} else {0}
##    ISV		<- avg.int.length + pen.obsL + pen.obsU

##    QEVL		<- ATWE(vec.obs = y.obs, vec.est = y.frcL, tau.ATWE = tauL.f3)
##    QEVU		<- ATWE(vec.obs = y.obs, vec.est = y.frcU, tau.ATWE = tauU.f3)

##    all.int.lengths.pos	<- all((y.frcU - y.frcL) > 0)
##    int.cover.pct	<- sum((y.obs > y.frcL) & (y.obs < y.frcU))/length(y.obs)*100

##    obj.value	<- c(avg.int.length, pen.obsL, pen.obsU, ISV, QEVL, QEVU, all.int.lengths.pos, int.cover.pct)



}	# regtype LL

}	# computability

#  if(show.progress.f3){print(paste("Nonparametric computation using bandwidths of row ", ro, " of ", nrow(bw.grid.f3), sep = ""))}
  if(show.progress.f3){print(paste("NPS computation using bandwidths of row ", ro, " of ", nrow(bw.grid.f3), " for ", regtype.f3, " h =", h.f3, "done.", sep = " "))}

##  return(obj.value)
##  return(list(y.frcL = y.frcL, y.frcU = y.frcU))
  return(c(y.frcL, y.frcU))
}



	



###
###	Nonparametric forecasts
###

np.frc	<- function(	# function as able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro,
  bw.grid.f4		= bw.grid,
  ndat.f4			= ndat,
  tau.f4			= tau,
  sample.range.set.f4	= sample.range.set,
#  crit.f4			= "MSEh",	# one of c("MSEh", "LSCV", "AICc")
  regtype.f4		= "LC",	# one of c("LC", "LL")
  h.f4			= 1,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU,
  ...
){

#	ro	<- 105
  b	<- bw.grid.f4[ro, ]
  W.s	<- sapply(
    X		= ndat.f4$s, 
    FUN	= weight.x,
    x.f0	= ndat.f4$s,
    S.f0	= max(ndat.f4$s),
    krnl.f0	= krnl.s.f4,
    b.f0	= as.numeric(b["bs"])
  )
  W.t	<- sapply(
    X		= ndat.f4$t, 
    FUN	= weight.x,
    x.f0	= ndat.f4$t,
    S.f0	= max(ndat.f4$t),
    krnl.f0	= krnl.t.f4,
    b.f0	= as.numeric(b["bt"])
  )
  W	<- W.s*W.t

  if(regtype.f4 == "LL"){
    y.frc		<- sapply(X = 1:nrow(sample.range.set.f4), FUN = fdr.LL, sample.range.set.f1 = sample.range.set.f4, h.f1 = h.f4, W.f1 = W, ndat.f1 = ndat.f4, tau.f1 = tau.f4)
#    obj.value	<- mean((ndat.f4[sample.range.set.f4[, "est.end.set"] + h.f4, "y"] - y.frc)^2)
#    obj.value	<- ATWE(vec.obs = ndat.f4[sample.range.set.f4[, "est.end.set"] + h.f4, "y"], vec.est = y.frc, tau.ATWE = tau.f4)
  }
  return(y.frc)
}








###
###	Nonparametric bandwidth estimation and forecasting
###

f.np		<- function(
  bw.grid.f5		= bw.grid,
  ndat.f5			= ndat,
  tauL.f5			= tauL,
  tauU.f5			= tauU,
  sample.range.set.f5	= sample.range.set,
  h.f5			= 1,
#  crit.f5			= "all",
  regtype.f5		= "LL",
  krnl.s.f5			= krnl.SEN,
  krnl.t.f5			= krnl.GAU,
  ...
){


###	Bandwidth estimation

  obj.temp	<- sapply(
#  bw.grid.f5$obj.value	<- sapply(
#  bw.grid.f5[, 3:6]	<- sapply(
    X				= 1:nrow(bw.grid.f5), 
#    X				= 1, 
    FUN			= np.reg, 
    bw.grid.f3		= bw.grid.f5,		# 2021-01-12
    ndat.f3			= ndat.f5,
    tauL.f3			= tauL.f5,
    tauU.f3			= tauU.f5,
    sample.range.set.f3	= sample.range.set.f5,
    h.f3			= h.f5,
#    crit.f3			= crit.f5, 
    regtype.f3		= regtype.f5, 
    krnl.s.f3		= krnl.s.f5, 
    krnl.t.f3		= krnl.t.f5, 
    show.progress.f3	= TRUE
  )
  obj.temp	<- as.data.frame(t(obj.temp))
##  colnames(obj.temp)	<- c("avg.int.length", "pen.obsL", "pen.obsU", "ISV", "QEVL", "QEVU", "all.int.lengths.pos", "int.cover.pct")







if(FALSE){		### COMMENTED OUT - START
###	Nonparametric forecast

  ISV.min.row	<- as.numeric(rownames(obj.temp[obj.temp[, "all.int.lengths.pos"] == 1, ][which.min(obj.temp[obj.temp[, "all.int.lengths.pos"] == 1, "ISV"]), ]))
  np.forecastsL.ISV	<- np.frc(
    ro			= ISV.min.row, 
    bw.grid.f4		= bw.grid.f5,
    ndat.f4			= ndat.f5,
    tau.f4			= tauL.f5,
    sample.range.set.f4	= sample.range.set.f5,
    h.f4			= h.f5,
    regtype.f4		= regtype.f5, 
    krnl.s.f4		= krnl.s.f5, 
    krnl.t.f4		= krnl.t.f5
  )

  np.forecastsU.ISV	<- np.frc(
    ro			= ISV.min.row, 
    bw.grid.f4		= bw.grid.f5,
    ndat.f4			= ndat.f5,
    tau.f4			= tauU.f5,
    sample.range.set.f4	= sample.range.set.f5,
    h.f4			= h.f5,
    regtype.f4		= regtype.f5, 
    krnl.s.f4		= krnl.s.f5, 
    krnl.t.f4		= krnl.t.f5
  )

  np.forecastsL.QEV	<- np.frc(
    ro			= which.min(obj.temp[, "QEVL"]), 
    bw.grid.f4		= bw.grid.f5,
    ndat.f4			= ndat.f5,
    tau.f4			= tauL.f5,
    sample.range.set.f4	= sample.range.set.f5,
    h.f4			= h.f5,
    regtype.f4		= regtype.f5, 
    krnl.s.f4		= krnl.s.f5, 
    krnl.t.f4		= krnl.t.f5
  )

  np.forecastsU.QEV	<- np.frc(
    ro			= which.min(obj.temp[, "QEVU"]), 
    bw.grid.f4		= bw.grid.f5,
    ndat.f4			= ndat.f5,
    tau.f4			= tauU.f5,
    sample.range.set.f4	= sample.range.set.f5,
    h.f4			= h.f5,
    regtype.f4		= regtype.f5, 
    krnl.s.f4		= krnl.s.f5, 
    krnl.t.f4		= krnl.t.f5
  )


###	to return

  to.return	<- list(
    np.bw.grid		= data.frame(bw.grid.f5, obj.temp),
    bw.ISV			= bw.grid.f5[ISV.min.row, ],
    bw.QEVL			= bw.grid.f5[which.min(obj.temp[, "QEVL"]), ],
    bw.QEVU			= bw.grid.f5[which.min(obj.temp[, "QEVU"]), ],
    np.forecastsL.ISV	= np.forecastsL.ISV,
    np.forecastsU.ISV	= np.forecastsU.ISV,
    np.forecastsL.QEV	= np.forecastsL.QEV,
    np.forecastsU.QEV	= np.forecastsU.QEV
  )
}	### COMMENTED OUT - END

  to.return		<- list(
    obs	= ndat.f5[sample.range.set.f5[, "est.end.set"] + h.f5, "y"],
    prd.L	= cbind(bw.grid.f5, obj.temp[, 1:(ncol(obj.temp)/2)]), 
    prd.U	= cbind(bw.grid.f5, obj.temp[, (ncol(obj.temp)/2 + 1):ncol(obj.temp)])
  )
  return(to.return)
}



