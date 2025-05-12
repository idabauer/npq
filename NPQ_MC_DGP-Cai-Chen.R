
##################################################
###	Create data set of Cai & Chen (2006)
##################################################

###
###	Basic configuration
###

#	R	<-	1000

#	T	<- 200	# Number of time periods (w.r.t. 4 quarters per year)

S	<- 4

n	<- T/S	# Number of years -> Cai & Chen (2006) use n \in {50, 100, 300}



###
###	Error term (Note: this is the only stochastic part in the DGP)
###

#	eps.sd	<- 0.1
eps_t		<- array(data = rnorm(n = R*T, mean = 0, sd = eps.sd), dim = c(R, T))


###	AR1-error 

#	e.AR1.coef	<- 0.9	# value 0 yields iid error

e_t	<- array(data = NA, dim = c(R, T))
e_t[, 1]	<- eps_t[, 1]
for(t in 2:T){
  e_t[, t]	<- e.AR1.coef*e_t[, t-1] + eps_t[, t]
}





###
###	Functions for creating regression function parts alpha and beta
###

f.alpha	<- function(x){exp(-0.7 + 3.5*x)}


f.beta_1	<- function(x){-3.1*x^2 + 17.1*x^4 - 28.1*x^5 + 15.7*x^6}
f.beta_2	<- function(x){-0.5*x^2 + 15.7*x^6 - 15.2*x^7}
f.beta_4	<- function(x){-0.2 + 4.8*x^2 - 7.7*x^3}
f.beta_3	<- function(x){-f.beta_1(x) - f.beta_2(x) - f.beta_4(x)}





###
###	Create object "dat" that contains the time series components and systematics
###

###	t, i, and j

dat	<- data.frame(
  t	= 1:T,			# time period t
  i	= rep(1:n, each = 4),	# year i
  j	= rep(1:4, times = n)	# quarter j
)



###	alpha value for year i

#	alpha.constant	<- TRUE
if(alpha.constant){
  dat$alpha	<- exp(1)
} else {
  dat$alpha	<- f.alpha(dat$i/n)
}


###	beta value for year i and quarter j

for(ro in 1:nrow(dat)){
  dat$beta[ro]	<- do.call(
    what	= get(paste("f.beta_", dat[ro, "j"], sep = "")),
    args	= list(dat[ro, "i"]/n)
  )
}


alpha.beta	<- matrix(data = dat$alpha + dat$beta, nrow = R, ncol = T, byrow = TRUE)

x_t	<- alpha.beta + e_t



