library(fixest)
library(ggplot2)

source("SetParameters.R")

set.seed(123)

dt = data.frame(group=rep(1:ngroups,each=ntimes),time=rep(1:ntimes,ngroups))

dt$c = dt$time*slopes[dt$group]
dt$x = dt$c + as.vector(replicate(ngroups,as.vector(arima.sim(model=list(order=c(1,0,0),ar=theta_AR),n=ntimes,sd=sd_AR))))

reg = feols(c~x|group+time,data=dt)
bias = unname(reg$coefficients) # bias of the estimator of beta conditionally on x
df_t = degrees_freedom(reg,type="t",vcov="iid")
df_t_tw = degrees_freedom(reg,type="t",vcov="twoway")

onesimul_xfixed = function(){
# one simulation with x fixed

	dt$y = dt$c + rnorm(ngroups*ntimes,0,sd_error)
	reg = feols(y~x|group+time,data=dt)
	est = unname(reg$coefficients)
	se_iid = sqrt(vcov(reg,vcov="iid")[1,1])
	se_hc = sqrt(vcov(reg,vcov="hetero")[1,1])
	se_tw = sqrt(vcov(reg,vcov="twoway")[1,1])
	return(c(est=est,se_iid=se_iid,se_hc=se_hc,se_tw=se_tw))
}

onesimul_xvar = function(){
# one simulation with x random

	dt$x = dt$c + as.vector(replicate(ngroups,as.vector(arima.sim(model=list(order=c(1,0,0),ar=theta_AR),n=ntimes,sd=sd_AR))))
	dt$y = dt$c + rnorm(ngroups*ntimes,0,sd_error)
	reg = feols(y~x|group+time,data=dt)
	est = unname(reg$coefficients)
	se_iid = sqrt(vcov(reg,vcov="iid")[1,1])
	se_hc = sqrt(vcov(reg,vcov="hetero")[1,1])
	se_tw = sqrt(vcov(reg,vcov="twoway")[1,1])
	bias = unname(feols(c~x|group+time,data=dt)$coefficients)
	return(c(est=est,se_iid=se_iid,se_hc=se_hc,se_tw=se_tw,bias=bias))
}

nsimul = 1000 # number of simulations
simul_xfixed = replicate(nsimul,onesimul_xfixed())
simul_xvar = replicate(nsimul,onesimul_xvar())
bias_var_est = mean(simul_xvar["bias",]) # approximate unconditional bias of the estimator of beta

# compute p-values for simulations with x fixed, testing H0:beta=0, with various covariance matrices
pval_xfixed_iid = 2*pt(abs(simul_xfixed["est",])/simul_xfixed["se_iid",],df=df_t,lower.tail=FALSE)
pval_xfixed_hc = 2*pt(abs(simul_xfixed["est",])/simul_xfixed["se_hc",],df=df_t,lower.tail=FALSE)
pval_xfixed_tw = 2*pt(abs(simul_xfixed["est",])/simul_xfixed["se_tw",],df=df_t_tw,lower.tail=FALSE)

# compute p-values for simulations with x fixed, testing H0:beta=bias, with various covariance matrices
pval_xfixed_bias_iid = 2*pt(abs(simul_xfixed["est",]-bias)/simul_xfixed["se_iid",],df=df_t,lower.tail=FALSE)
pval_xfixed_bias_hc = 2*pt(abs(simul_xfixed["est",]-bias)/simul_xfixed["se_hc",],df=df_t,lower.tail=FALSE)
pval_xfixed_bias_tw = 2*pt(abs(simul_xfixed["est",]-bias)/simul_xfixed["se_tw",],df=df_t_tw,lower.tail=FALSE)

# compute p-values for simulations with x random, testing H0:beta=0, with various covariance matrices
pval_xvar_iid = 2*pt(abs(simul_xvar["est",])/simul_xvar["se_iid",],df=df_t,lower.tail=FALSE)
pval_xvar_hc = 2*pt(abs(simul_xvar["est",])/simul_xvar["se_hc",],df=df_t,lower.tail=FALSE)
pval_xvar_tw = 2*pt(abs(simul_xvar["est",])/simul_xvar["se_tw",],df=df_t_tw,lower.tail=FALSE)

# compute p-values for simulations with x random, testing H0:beta=bias, with various covariance matrices
pval_xvar_bias_iid = 2*pt(abs(simul_xvar["est",]-bias_var_est)/simul_xvar["se_iid",],df=df_t,lower.tail=FALSE)
pval_xvar_bias_hc = 2*pt(abs(simul_xvar["est",]-bias_var_est)/simul_xvar["se_hc",],df=df_t,lower.tail=FALSE)
pval_xvar_bias_tw = 2*pt(abs(simul_xvar["est",]-bias_var_est)/simul_xvar["se_tw",],df=df_t_tw,lower.tail=FALSE)

# cumulative distribution functions for p-values
F_xfixed_bias_iid = ecdf(pval_xfixed_bias_iid)
F_xfixed_bias_hc = ecdf(pval_xfixed_bias_hc)
F_xfixed_bias_tw = ecdf(pval_xfixed_bias_tw)
F_xvar_bias_iid = ecdf(pval_xvar_bias_iid)
F_xvar_bias_hc = ecdf(pval_xvar_bias_hc)
F_xvar_bias_tw = ecdf(pval_xvar_bias_tw)

# dataframe for plotting cumulative distribution functions of p-values for tests with H0:beta=bias
xaxis = seq(from=0,to=1,by=0.001)
plot_xfixed = data.frame(xaxis,F_iid=sapply(xaxis,F_xfixed_bias_iid),F_hc=sapply(xaxis,F_xfixed_bias_hc),F_tw=sapply(xaxis,F_xfixed_bias_tw))
plot_xvar = data.frame(xaxis,F_iid=sapply(xaxis,F_xvar_bias_iid),F_hc=sapply(xaxis,F_xvar_bias_hc),F_tw=sapply(xaxis,F_xvar_bias_tw))

# when testing H0:beta=0 all p-values are tiny, so we just compute the maximum values and show in output file
sink("output.txt")
cat("Maximum p-value in ",nsimul, "simulations\n\n")
cat("Testing H0:beta=0, x fixed\n")
cat("vcov=iid: ",max(pval_xfixed_iid),"\n")
cat("vcov=hetero: ",max(pval_xfixed_hc),"\n")
cat("vcov=twoway: ",max(pval_xfixed_tw),"\n")
cat("\n")
cat("Testing H0:beta=0, x random\n")
cat("vcov=iid: ",max(pval_xvar_iid),"\n")
cat("vcov=hetero: ",max(pval_xvar_hc),"\n")
cat("vcov=twoway: ",max(pval_xvar_tw),"\n")
sink()

ggplot(data=plot_xfixed)+geom_line(aes(x=xaxis,y=F_iid,col="iid"))+
	geom_line(aes(x=xaxis,y=F_hc,col="hc"))+
	geom_line(aes(x=xaxis,y=F_tw,col="tw"))+
	geom_abline(intercept=0,slope=1,col="red")+
	scale_x_continuous(limits=c(0,1),name="significance level")+scale_y_continuous(limits=c(0,1),name="type I error")+
	scale_color_manual(values=c("iid"="black","hc"="green","tw"="orange"),
		labels=c("iid"="iid","hc"="hetero","tw"="twoway"),name="covariance matrix")+
	labs(title=paste(nsimul,"simulations with x fixed, testing H0:beta=bias"))
ggsave(filename="pval_xfixed.jpg",device="jpeg",width=14.23,height=6.77)

ggplot(data=plot_xvar)+geom_line(aes(x=xaxis,y=F_iid,col="iid"))+
	geom_line(aes(x=xaxis,y=F_hc,col="hc"))+
	geom_line(aes(x=xaxis,y=F_tw,col="tw"))+
	geom_abline(intercept=0,slope=1,col="red")+
	scale_x_continuous(limits=c(0,1),name="significance level")+scale_y_continuous(limits=c(0,1),name="type I error")+
	scale_color_manual(values=c("iid"="black","hc"="green","tw"="orange"),
		labels=c("iid"="iid","hc"="hetero","tw"="twoway"),name="covariance matrix")+
	labs(title=paste(nsimul,"simulations with x random, testing H0:beta=bias"))
ggsave(filename="pval_xvar.jpg",device="jpeg",width=14.23,height=6.77)
