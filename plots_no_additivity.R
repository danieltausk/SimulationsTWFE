library(fixest)
library(ggplot2)

source("SetParameters.R")

set.seed(123)

dt = data.frame(group=rep(1:ngroups,each=ntimes),time=rep(1:ntimes,ngroups))

dt$c = dt$time*slopes[dt$group]
dt$x = dt$c + as.vector(replicate(ngroups,as.vector(arima.sim(model=list(order=c(1,0,0),ar=theta_AR),n=ntimes,sd=sd_AR))))
dt$y = dt$c + rnorm(ngroups*ntimes,0,sd_error)

reg = feols(y~x|group+time,data=dt)
reg_exp = feols(c~x|group+time,data=dt)
dt$res = reg$residuals
dt$res_exp = reg_exp$residuals # expected values of residuals

pdf("xy_no_additivity.pdf",width=14.23,height=6.77)
for(i in 1:ngroups)
{
	print(ggplot(data=dt[dt$group == i,])+geom_line(aes(x=time,y=y,col="y"))+geom_line(aes(x=time,y=x,col="x"))+
			scale_y_continuous(limits=c(0,ntimes),name=NULL)+
			labs(title=paste("Group",i)))
}
dev.off()

pdf("residuals_no_additivity.pdf",width=14.23,height=6.77)
for(i in 1:ngroups)
{
	print(ggplot(data=dt[dt$group == i,])+geom_point(aes(x=time,y=res,col="residuals"))+
			geom_line(aes(x=time,y=res_exp,linetype="expected values"),col="red")+
			scale_linetype_manual(values="solid",name=NULL)+
			scale_color_manual(values="black",name=NULL)+
			labs(title=paste("Group",i),y=NULL))
}
dev.off()
