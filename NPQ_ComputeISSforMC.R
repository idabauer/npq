
time.iss.A		<- Sys.time()

iss.model	<- ets(
  y		= dat,
  model	= "ZZZ", 
  allow.multiplicative.trend	= TRUE
#,
#  opt.crit	= "mae"
)

time.iss.B		<- Sys.time()
#iss.model.time	<- difftime(time.iss.B, time.iss.A, unit = time.unit)


iss.results		<- list(
  iss.model			= iss.model,
  computation.time	= rep(NA, length(h.set) + 1)
)
names(iss.results$computation.time)	<- c("time.model", paste("time.h", h.set, sep = ""))

iss.results[["computation.time"]]["time.model"]	<- difftime(time.iss.B, time.iss.A, unit = time.unit)


for(h in h.set){
#	h	<- 4	

  time.iss.C		<- Sys.time()

  f.sample.range.set(dat = dat.input, h = h, extending = extending.window, est.end.set = insample.end:(T - h), 
    start.period = insample.end - (sample.range.start - 1)	# check with h = 1, 3, etc.
  )

  iss.frc.test	<- NULL

  for(ro in 1:nrow(sample.range.set)){
  #	ro	<- 1
  #	ro	<- 2

    iss.test	<- ets(
      y		= window(x = dat.input, start = time(dat.input)[sample.range.set[ro, ][1]], end = time(dat.input)[sample.range.set[ro, ][2]]),
      model		= paste(iss.model$components[1:3], collapse = ""), 
      damped	= as.logical(iss.model$components[4])
#,
#      opt.crit	= "mae"
    )

    iss.frc.test	<- rbind(iss.frc.test, as.data.frame(forecast(object = iss.test, h = h, level = tauU - tauL))[h, ])

  }

  iss.results[[paste("prd.test.h", h, sep = "")]]		<- iss.frc.test

  time.iss.D		<- Sys.time()

  iss.results[["computation.time"]][paste("time.h", h, sep = "")]	<- difftime(time.iss.D, time.iss.C, unit = time.unit)





  obs	<- window(
    x		= dat.input, 
    start	= time(dat.input)[sample.range.set[1, ][2] + h],
    end	= time(dat.input)[sample.range.set[ro, ][2] + h]
  )
  obs		<- if(logarithmize){exp(obs)} else {obs}
  limit.L	<- if(logarithmize){exp(iss.frc.test[, 2])} else {iss.frc.test[, 2]}
  limit.U	<- if(logarithmize){exp(iss.frc.test[, 3])} else {iss.frc.test[, 3]}


  coverage		<- sum((obs >= limit.L) & (obs <= limit.U))/length(obs)*100
  overpred.pct	<- sum(obs < limit.L)/length(obs)*100
  underpred.pct	<- sum(obs > limit.U)/length(obs)*100
  total.pct		<- overpred.pct + coverage + underpred.pct
  pct.int.length.pos	<- sum(limit.U - limit.L > 0)/length(obs)*100
  avg.int.length	<- mean(limit.U - limit.L)
  overpred.pen	<- if(any(obs < limit.L)){sum((1/tauL)*(limit.L - obs)[obs < limit.L])/length(obs)} else {0}
  underpred.pen	<- if(any(obs > limit.U)){sum((1/(1 - tauU))*(obs - limit.U)[obs > limit.U])/length(obs)} else {0}
  intervalscore	<- avg.int.length + overpred.pen + underpred.pen
  atwe.L		<- ATWE(vec.obs = obs, vec.est = limit.L, tau.ATWE = tauL)
  atwe.U		<- ATWE(vec.obs = obs, vec.est = limit.U, tau.ATWE = tauU)
  atwe.LpU		<- atwe.L + atwe.U

  iss.results[[paste("results.test.h", h, sep = "")]]		<- c(
    coverage		= coverage, 
    overpred.pct		= overpred.pct, 
    underpred.pct		= underpred.pct, 
    total.pct		= total.pct, 
    pct.int.length.pos	= pct.int.length.pos, 
    avg.int.length	= avg.int.length, 
    overpred.pen		= overpred.pen, 
    underpred.pen		= underpred.pen, 
    intervalscore		= intervalscore, 
    atwe.L			= atwe.L, 
    atwe.U			= atwe.U, 
    atwe.LpU		= atwe.LpU
  )

}


  assign(
    x		= paste("iss.T", T, ".r", r, sep = ""),
    value	= iss.results
  )

