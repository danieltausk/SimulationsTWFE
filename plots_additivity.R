library(fixest)
library(ggplot2)

source("SetParameters.R")

set.seed(123)

dt = data.frame(group=rep(1:ngroups,each=ntimes),time=rep(1:ntimes,ngroups))

dt$c = dt$time+intercepts[dt$group]
dt$x = dt$c + as.vector(replicate(ngroups,as.vector(arima.sim(model=list(order=c(1,0,0),ar=theta_AR),n=ntimes,sd=sd_AR))))
dt$y = dt$c + rnorm(ngroups*ntimes,0,sd_error)

reg = feols(y~x|group+time,data=dt)
dt$res = reg$residuals

pdf("xy_additivity.pdf",width=14.23,height=6.77)
for(i in 1:ngroups)
{
	print(ggplot(data=dt[dt$group == i,])+geom_line(aes(x=time,y=y,col="y"))+geom_line(aes(x=time,y=x,col="x"))+
			scale_y_continuous(limits=c(-20,120),name=NULL)+
			labs(title=paste("Group",i)))
}
dev.off()

pdf("residuals_additivity.pdf",width=14.23,height=6.77)
for(i in 1:ngroups)
{
	print(ggplot(data=dt[dt$group == i,])+geom_point(aes(x=time,y=res,col="residuals"))+
			geom_hline(yintercept=0,col="red")+
			scale_color_manual(values="black",name=NULL)+
			labs(title=paste("Group",i),y=""))
}
dev.off()
