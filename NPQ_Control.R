
###
###	(Unit for) time measurement
###

time.unit	<- "mins"

time.A	<- Sys.time()








###
###	Create estimation windows, nonparametric training data set, and bandwidth grid
###

###	Determine estimation windows for training sample

f.sample.range.set(dat = dat, h = h, extending = extending.window, est.end.set = sample.range.start:(length(dat) - h))
#	sample.range.set




###	Create training-ndat for nonparametric regression, season and trend are separate columns

#	for(trend.year in trend.year.set){

#	trend.year	<- "year"		# trend as year (first year = 1)
#	trend.year	<- "t"		# trend as t

  ndat	<- data.frame(
    y	= as.numeric(dat), 	# changed: 2021-01-12
    t	= if(trend.year == "year"){as.integer(floor(time(dat)) - floor(min(time(dat))) + 1)} else {1:length(dat)}
  )
  ndat$t	<- ndat$t/ndat$t[insample.end]

  ndat$s	<- (start(dat)[2] + -1:(length(dat) - 2)) %% (frequency(dat) ) + 1





###	Create bandwidth grid

#	bt.default	<- 1.06*sd(ndat$t[1:insample.end])*insample.end^(-0.2)

bt.default	<- 1.06*sd(ndat$t[1:sample.range.start])*sample.range.start^(-0.2)
bt.seq	<- exp(bt.seq.equi)*bt.default

bw.grid	<- expand.grid(
  bt	= bt.seq,
  bs	= bs.seq
)







###
###	Compute and evaluate nonparametric quantile regressions on training sample
###

###	Compute nonparametric quantile regression (lower and upper quantile) for every row of bandwidth grid

#  for(crit.temp in crit.set){
#    for(regtype.temp in regtype.set){
#      for(krnl.s.temp in krnl.s.set){

#        time.A.np	<- Sys.time()
        np.results		<- f.np(
          ndat.f5			= ndat, 
          sample.range.set.f5	= sample.range.set, 
          tauL.f5			= tauL,
          tauU.f5			= tauU,
          bw.grid.f5		= bw.grid, 
#          bw.grid.f5		= bw.grid[1:2, ], 		#### XXX: to change
          h.f5			= h,
#          crit.f5			= crit.temp,
          regtype.f5		= regtype.temp,
          krnl.s.f5		= get(krnl.s.temp),
          krnl.t.f5		= krnl.GAU
        )









###	Evaluate all combinations of upper and lower bandwidth grid row

results.train	<- expand.grid(bt.seq, bs.seq, bt.seq, bs.seq)
colnames(results.train)	<- c("bt.L", "bs.L", "bt.U", "bs.U")
results.train	<- cbind(
  results.train, 
  matrix(
    data = NA, 
    nrow = nrow(results.train), 
    ncol = 12, 
    dimnames = list(
      NULL, c(
        "coverage", "overpred.pct", "underpred.pct", "total.pct", "pct.int.length.pos", "avg.int.length", 
        "overpred.pen", "underpred.pen", "intervalscore", "atwe.L", "atwe.U", "atwe.LpU"
      )
    )
  )
)

results.test		<- results.train[1:2, ]



for(ro in 1:nrow(results.train)){
  b.L		<- as.numeric(results.train[ro, 1:2])
  b.U		<- as.numeric(results.train[ro, 3:4])

  row.L	<- which(b.L[1] == np.results$prd.L[, 1] & b.L[2] == np.results$prd.L[, 2])
  row.U	<- which(b.U[1] == np.results$prd.U[, 1] & b.U[2] == np.results$prd.U[, 2])

  obs		<- if(logarithmize){exp(np.results$obs)} else {np.results$obs}
  limit.L	<- as.numeric(if(logarithmize){exp(np.results$prd.L[row.L, -(1:2)])} else {np.results$prd.L[row.L, -(1:2)]})
  limit.U	<- as.numeric(if(logarithmize){exp(np.results$prd.U[row.U, -(1:2)])} else {np.results$prd.U[row.U, -(1:2)]})

  coverage		<- sum((obs >= limit.L) & (obs <= limit.U))/length(obs)*100
  overpred.pct	<- sum(obs < limit.L)/length(obs)*100
  underpred.pct	<- sum(obs > limit.U)/length(obs)*100
  total.pct		<- overpred.pct + coverage + underpred.pct
  pct.int.length.pos	<- sum(limit.U - limit.L > 0)/length(obs)*100
  avg.int.length	<- mean(limit.U - limit.L)
  overpred.pen	<- if(any(obs < limit.L)){sum((1/tauL)*(limit.L - obs)[obs < limit.L])/length(obs)} else {0}
  underpred.pen	<- if(any(obs > limit.U)){sum((1/(1 - tauU))*(obs - limit.U)[obs > limit.U])/length(obs)} else {0}
#  overpred.pen	<- if(any(obs < limit.L)){mean((1/tauL)*(limit.L - obs)[obs < limit.L])} else {0}
#  underpred.pen	<- if(any(obs > limit.U)){mean((1/(1 - tauU))*(obs - limit.U)[obs > limit.U])} else {0}
  intervalscore	<- avg.int.length + overpred.pen + underpred.pen
  atwe.L		<- ATWE(vec.obs = obs, vec.est = limit.L, tau.ATWE = tauL)
  atwe.U		<- ATWE(vec.obs = obs, vec.est = limit.U, tau.ATWE = tauU)
  atwe.LpU		<- atwe.L + atwe.U

  results.train[ro, 5:ncol(results.train)]		<- c(
    coverage, overpred.pct, underpred.pct, total.pct, pct.int.length.pos, avg.int.length, 
    overpred.pen, underpred.pen, intervalscore, atwe.L, atwe.U, atwe.LpU
  )
#  print(paste("Row ", ro, " of ", nrow(results.train), " done.", sep = ""))
}

np.results$results.train	<- results.train





###	Determine best combination of upper/lower bandwidths w.r.t. IS and ISc

ISV.min.row		<- which.min(results.train$intervalscore)
atleast.coverage.rows	<- results.train$coverage > 100*(tauU - tauL)
if(sum(atleast.coverage.rows) > 0){
  ISVC.min.row	<- which(results.train$intervalscore == min(results.train$intervalscore[atleast.coverage.rows]))
} else {
  ISVC.min.row	<- which.max(results.train$coverage)
}

np.results$train.opt.rows	<- c(ISV = ISV.min.row, ISVC = ISVC.min.row)







###
###	Compute and evaluate nonparametric quantile regression for best upper/lower bandwidths on test sample
###

#	time.A	<- Sys.time()


f.sample.range.set(dat = dat.input, h = h, extending = extending.window, est.end.set = insample.end:(T - h), 
  start.period = insample.end - (sample.range.start - 1)	# check with h = 1, 3, etc.
)





###	Create ndat object

dat.temp	<- dat.input
ndat	<- data.frame(
  y	= as.numeric(dat.temp),
  t	= if(trend.year == "year"){as.integer(floor(time(dat.temp)) - floor(min(time(dat.temp))) + 1)} else {1:length(dat.temp)}
)
ndat$t	<- ndat$t/ndat$t[insample.end]
ndat$s	<- (start(dat.temp)[2] + -1:(length(dat.temp) - 2)) %% (frequency(dat.temp) ) + 1



###	Determine time series values of test-sample (to be forecast)

obs	<- ndat[sample.range.set[, "est.end.set"] + h, "y"]



###	Compute test sample interval limits of ISV-best upper/lower bandwidth combination

prd.L.ISV		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISV.min.row, "bt.L"], bs = results.train[ISV.min.row, "bs.L"]),
  ndat.f4			= ndat,
  tau.f4			= tauL,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",	# one of c("LC", "LL")
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)
prd.U.ISV		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISV.min.row, "bt.U"], bs = results.train[ISV.min.row, "bs.U"]),
  ndat.f4			= ndat,
  tau.f4			= tauU,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",	# one of c("LC", "LL")
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)



###	Compute test sample interval limits of ISVcorrected-best upper/lower bandwidth combination

prd.L.ISVC		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISVC.min.row, "bt.L"], bs = results.train[ISVC.min.row, "bs.L"]),
  ndat.f4			= ndat,
  tau.f4			= tauL,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",	# one of c("LC", "LL")
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)
prd.U.ISVC		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISVC.min.row, "bt.U"], bs = results.train[ISVC.min.row, "bs.U"]),
  ndat.f4			= ndat,
  tau.f4			= tauU,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",	# one of c("LC", "LL")
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)




###	Adjust results object

results.test$train.opt	<- c("ISV", "ISVC")
results.test[1, 1:4]	<- results.train[ISV.min.row, 1:4]
results.test[2, 1:4]	<- results.train[ISVC.min.row, 1:4]





###
#	time.B	<- Sys.time()
#	difftime(time.B, time.A, unit = "mins")
###






###	 Compute test sample measures (i.e., using 80%-100%-part of ts)

obs		<- if(logarithmize){exp(obs)} else {obs}

for(ro in 1:nrow(results.test)){

  prd.L	<- if(ro == 1){prd.L.ISV} else {prd.L.ISVC}
  prd.U	<- if(ro == 1){prd.U.ISV} else {prd.U.ISVC}

  limit.L	<- if(logarithmize){exp(prd.L)} else {prd.L}
  limit.U	<- if(logarithmize){exp(prd.U)} else {prd.U}


  coverage		<- sum((obs >= limit.L) & (obs <= limit.U))/length(obs)*100
  overpred.pct	<- sum(obs < limit.L)/length(obs)*100
  underpred.pct	<- sum(obs > limit.U)/length(obs)*100
  total.pct		<- overpred.pct + coverage + underpred.pct
  pct.int.length.pos	<- sum(limit.U - limit.L > 0)/length(obs)*100
  avg.int.length	<- mean(limit.U - limit.L)
  overpred.pen	<- if(any(obs < limit.L)){sum((1/tauL)*(limit.L - obs)[obs < limit.L])/length(obs)} else {0}
  underpred.pen	<- if(any(obs > limit.U)){sum((1/(1 - tauU))*(obs - limit.U)[obs > limit.U])/length(obs)} else {0}
#  overpred.pen	<- if(any(obs < limit.L)){mean((1/tauL)*(limit.L - obs)[obs < limit.L])} else {0}
#  underpred.pen	<- if(any(obs > limit.U)){mean((1/(1 - tauU))*(obs - limit.U)[obs > limit.U])} else {0}
  intervalscore	<- avg.int.length + overpred.pen + underpred.pen
  atwe.L		<- ATWE(vec.obs = obs, vec.est = limit.L, tau.ATWE = tauL)
  atwe.U		<- ATWE(vec.obs = obs, vec.est = limit.U, tau.ATWE = tauU)
  atwe.LpU		<- atwe.L + atwe.U

  results.test[ro, 5:(ncol(results.test) - 1)]		<- c(
    coverage, overpred.pct, underpred.pct, total.pct, pct.int.length.pos, avg.int.length, 
    overpred.pen, underpred.pen, intervalscore, atwe.L, atwe.U, atwe.LpU
  )

}




###
###	Save test sample objects
###


###	Mimic np.results$prd.L/U-objects

np.results$prd.L.test		<- cbind(
  rbind(results.train[ISV.min.row, c("bt.L", "bs.L")], results.train[ISVC.min.row, c("bt.L", "bs.L")]),
  rbind(prd.L.ISV, prd.L.ISVC), 
  train.opt = c("ISV", "ISVC")
)

np.results$prd.U.test		<- cbind(
  rbind(results.train[ISV.min.row, c("bt.U", "bs.U")], results.train[ISVC.min.row, c("bt.U", "bs.U")]),
  rbind(prd.U.ISV, prd.U.ISVC), 
  train.opt = c("ISV", "ISVC")
)



np.results$results.test		<- results.test







###
###	Compute real forecasts (i.e., out-of-sample/-support)
###

###	Create data set and sample start/end

dat.tmp	<- ts(NA, start = start(dat.input), end = end(dat.input) + c(1, 0), frequency = frequency(dat.input))
 ndat	<- data.frame(
    y	= c(as.numeric(dat.input), rep(NA, times = frequency(dat.input))),
    t	= if(trend.year == "year"){as.integer(floor(time(dat.tmp)) - floor(min(time(dat.tmp))) + 1)} else {1:(length(dat.input) + frequency(dat.input))}
  )
  ndat$t	<- ndat$t/ndat$t[insample.end]
  ndat$s	<- (start(dat.tmp)[2] + -1:(length(dat.tmp) - 2)) %% (frequency(dat.tmp) ) + 1

sample.range.set	<- matrix(data = c(if(extending.window){1} else {T - sample.range.start + 1}, T), nrow = 1, ncol = 2, dimnames = list(NULL, c("est.start.set", "est.end.set")))







###	Compute out-of-sample interval limits of ISV-best upper/lower training-bandwidth combination

oos.prd.L.ISV		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISV.min.row, "bt.L"], bs = results.train[ISV.min.row, "bs.L"]),
  ndat.f4			= ndat,
  tau.f4			= tauL,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)

oos.prd.L.ISVC		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISVC.min.row, "bt.L"], bs = results.train[ISVC.min.row, "bs.L"]),
  ndat.f4			= ndat,
  tau.f4			= tauL,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)




oos.prd.U.ISV		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISV.min.row, "bt.U"], bs = results.train[ISV.min.row, "bs.U"]),
  ndat.f4			= ndat,
  tau.f4			= tauU,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)

oos.prd.U.ISVC		<- np.frc(	# function is able to conduct real forecasts (i.e. not restricted to nrow(ndat) - h), then, sample.range.set.f4 has to be adjusted
  ro	= 1,
  bw.grid.f4		= data.frame(bt = results.train[ISVC.min.row, "bt.U"], bs = results.train[ISVC.min.row, "bs.U"]),
  ndat.f4			= ndat,
  tau.f4			= tauU,
  sample.range.set.f4	= sample.range.set,
  regtype.f4		= "LL",
  h.f4			= h,
  krnl.s.f4			= krnl.SEN,
  krnl.t.f4			= krnl.GAU
)





###	Save oos-results

np.results$prd.L.oos		<- cbind(
  rbind(results.train[ISV.min.row, c("bt.L", "bs.L")], results.train[ISVC.min.row, c("bt.L", "bs.L")]),
  rbind(oos.prd.L.ISV, oos.prd.L.ISVC), 
  train.opt = c("ISV", "ISVC")
)


np.results$prd.U.oos		<- cbind(
  rbind(results.train[ISV.min.row, c("bt.U", "bs.U")], results.train[ISVC.min.row, c("bt.U", "bs.U")]),
  rbind(oos.prd.U.ISV, oos.prd.U.ISVC), 
  train.opt = c("ISV", "ISVC")
)



#	Add computation time

time.B	<- Sys.time()
np.results[["computation.time"]]	<- difftime(time.B, time.A, unit = time.unit)



###	Drop obs-object

np.results$obs	<- NULL




