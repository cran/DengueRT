
###packages
require(nlme)
require(ggplot2)
require(drc)

            

####Richards model funtions

Richardsm<-function(inc, time, start=NULL){
if(is.null(start)){
tt<-data.frame(inc,time)
gi<-stats::getInitial(cumsum(inc) ~ SSlogis(time, Asym, xmid, scal), data = tt)
c<-as.data.frame(gi)
pars <- unname(coef(nls(cumsum(inc) ~ (1+exp((xmid - time)/scal))^(-exp(-lpow)),
                            data = tt, c(xmid = c$gi[2], scal = c$gi[3],
                                         lpow = 0.001), algorithm = "plinear")))
start=list(alpha=c$gi[1],k=pars[3],gamma=1/c$gi[3],eta=c$gi[2])
start=start}
gnls.richards<-gnls(cumsum(inc) ~ alpha*((1+k*exp(-gamma*k*(time-eta)))^(-1/k)),start=start)
p<-summary(gnls.richards);
return(p)}


####3P logistic model 
logistic3pm <-function(inc, time, start=NULL){
if(is.null(start)){
tt<-data.frame(inc,time)
gi<-stats::getInitial(cumsum(inc) ~ SSlogis(time, Asym, xmid, scal), data = tt)
c<-as.data.frame(gi)
start=list(alpha=c$gi[1],gamma=1/c$gi[3],eta=c$gi[2])
start=start}
gnls.3plogistic<-gnls(cumsum(inc)~alpha/(1+exp(-gamma*(time-eta))),start=start)
p<-summary (gnls.3plogistic);
return(p)}


####Sigmoidal Emax Model
SigmEmaxm <-function(inc, time, start=NULL){
if(is.null(start)){
tt<-data.frame(inc,time)
gi<-drm(cumsum(inc) ~time, data=tt,fct=LL.4())
start=list(alpha=gi$fit$par[3],eta=gi$fit$par[4],beta=gi$fit$par[2],n=-gi$fit$par[1])
start=start}
gnls.Emax<-gnls(cumsum(inc)~beta + (time^n)*(alpha-beta)/(time^n + eta^n),start=start)
p<-summary (gnls.Emax);
return(p)}



####Gompertz model 4 parametrs
gompertzm <-function(inc, time, start=NULL){
if(is.null(start)){
tt<-data.frame(inc,time)
gi<-stats::getInitial(cumsum(inc) ~ SSgompertz(time, Asym, b2, b3), data = tt)
c<-as.data.frame(gi)
start=list(alpha=c$gi[1],eta=-log(c$gi[2])/log(c$gi[3]),beta=0, gamma=-log(c$gi[3]))
start=start}
gnls.gompertz<-gnls(cumsum(inc)~beta + (alpha-beta)*exp(-exp(-gamma*(time-eta))),start=start)
p<-summary (gnls.gompertz);
return(p)}



####Modelo de Weibull
weibullm <-function(inc, time, start=NULL){
if(is.null(start)){
tt<-data.frame(inc,time)
gi<-stats::getInitial(cumsum(inc) ~ SSweibull(time, Asym, Drop, lrc, pwr), data = tt)
c<-as.data.frame(gi)
start=list(alpha=c$gi[1],eta=exp(-c$gi[3]/c$gi[4]),beta=c$gi[1]-c$gi[2], k=c$gi[4])
start=start}
gnls.weibull<-gnls(cumsum(inc)~beta + (alpha-beta)*exp(-(time/eta)^-k),start=start)
p<-summary(gnls.weibull);
return(p)}


####5p logistic model
logistic5pm <-function(inc, time, start=NULL){
if(is.null(start)){
tt<-data.frame(inc,time)
gi<-drm(cumsum(inc) ~time, data=tt,fct=LL.5())
gi2<-drm(cumsum(inc) ~time, data=tt,fct=LL.4())
start=list(alpha=gi$fit$par[2],beta=gi$fit$par[3],g=gi$fit$par[5],eta=gi2$fit$par[4], k=-gi$fit$par[1])
start=start}
gnls.5plogistic<-gnls(cumsum(inc)~beta+ (alpha-beta)/(1+(2^(1/g)-1)*(time/eta)^-k)^g,start=start)
p<-summary (gnls.5plogistic);
return(p)}


#####model averaging weight function
weight<-function(inc , time, start=NULL){
aic<-AIC(Richardsm(inc , time, start[[1]]),logistic3pm(inc , time, start[[2]]),SigmEmaxm (inc , time, start[[3]]),
     gompertzm(inc , time, start[[4]]),weibullm(inc , time, start[[5]]),logistic5pm (inc , time, start[[6]]))
minAIC<-min(aic[,2])
deltaAIC<-aic[,2]-minAIC
w<-signif(exp(-deltaAIC/2)/sum(exp(-deltaAIC/2)),6)
names(w)<-c("Richards","3P Logistic","Sigmoidal Emax","Gompertz","Weibull","5P Logistic");
return(w)}



####Model averaging estimator Final size
estimFSizeMA<-function(inc , time, start=NULL){
alpha.MA<-weight(inc , time, start)[1]*coef(Richardsm(inc , time, start[[1]]))[1] + weight(inc , time, start)[2]*coef(logistic3pm(inc , time, start[[2]]))[1] + weight(inc , time, start)[3]*coef(SigmEmaxm (inc , time, start[[3]]))[1] +weight(inc , time, start)[4]*coef(gompertzm(inc , time, start[[4]]))[1] +weight(inc , time, start)[5]*coef(weibullm(inc , time, start[[5]]))[1] + weight(inc , time, start)[6]*coef(logistic5pm (inc , time, start[[6]]))[1]
stderror.alpha.MA<-weight(inc , time, start)[1]*sqrt(Richardsm(inc , time, start[[1]])$tTable[1, 2]^2 + (coef(Richardsm(inc , time, start[[1]]))[1]-alpha.MA)^2) + 
weight(inc , time, start)[2]*sqrt(logistic3pm(inc , time, start[[2]])$tTable[1, 2]^2 + (coef(logistic3pm(inc , time, start[[2]]))[1]-alpha.MA)^2) + 
weight(inc , time, start)[3]*sqrt(SigmEmaxm (inc , time, start[[3]])$tTable[1, 2]^2 + (coef(SigmEmaxm (inc , time, start[[3]]))[1]-alpha.MA)^2) + 
weight(inc , time, start)[4]*sqrt(gompertzm(inc , time, start[[4]])$tTable[1, 2]^2 + (coef(gompertzm(inc , time, start[[4]]))[1]-alpha.MA)^2) + 
weight(inc , time, start)[5]*sqrt(weibullm(inc , time, start[[5]])$tTable[1, 2]^2 + (coef(weibullm(inc , time, start[[5]]))[1]-alpha.MA)^2) + 
weight(inc , time, start)[6]*sqrt(logistic5pm (inc , time, start[[6]])$tTable[1, 2]^2 + (coef(logistic5pm (inc , time, start[[6]]))[1]-alpha.MA)^2)
low<-alpha.MA-qt(0.975,5)*stderror.alpha.MA
up<-alpha.MA+qt(0.975,5)*stderror.alpha.MA
names(alpha.MA)<-c("est.")
names(low)<-c("lower")
names(up)<-c("upper")
p<-c(low,alpha.MA,up);
return(p)}


####Model averaging estimator Turning point 

estimTpointMA<-function(inc , time, start=NULL){
eta.MA<-weight(inc , time, start)[1]*coef(Richardsm(inc , time, start[[1]]))[4] + weight(inc , time, start)[2]*coef(logistic3pm(inc , time, start[[2]]))[3] + weight(inc , time, start)[3]*coef(SigmEmaxm (inc , time, start[[3]]))[2] +weight(inc , time, start)[4]*coef(gompertzm(inc , time, start[[4]]))[2] +weight(inc , time, start)[5]*coef(weibullm(inc , time, start[[5]]))[2] + weight(inc , time, start)[6]*coef(logistic5pm (inc , time, start[[6]]))[4]
stderror.eta.MA<-weight(inc , time, start)[1]*sqrt(Richardsm(inc , time, start[[1]])$tTable[4, 2]^2 + (coef(Richardsm(inc , time, start[[1]]))[4]-eta.MA)^2) + 
weight(inc , time, start)[2]*sqrt(logistic3pm(inc , time, start[[2]])$tTable[3, 2]^2 + (coef(logistic3pm(inc , time, start[[2]]))[3]-eta.MA)^2) + 
weight(inc , time, start)[3]*sqrt(SigmEmaxm (inc , time, start[[3]])$tTable[2, 2]^2 + (coef(SigmEmaxm (inc , time, start[[3]]))[2]-eta.MA)^2) + 
weight(inc , time, start)[4]*sqrt(gompertzm(inc , time, start[[4]])$tTable[2, 2]^2 + (coef(gompertzm(inc , time, start[[4]]))[2]-eta.MA)^2) + 
weight(inc , time, start)[5]*sqrt(weibullm(inc , time, start[[5]])$tTable[2, 2]^2 + (coef(weibullm(inc , time, start[[5]]))[2]-eta.MA)^2) + 
weight(inc , time, start)[6]*sqrt(logistic5pm (inc , time, start[[6]])$tTable[4, 2]^2 + (coef(logistic5pm (inc , time, start[[6]]))[4]-eta.MA)^2)
low<-eta.MA-qt(0.975,5)*stderror.eta.MA
up<-eta.MA+qt(0.975,5)*stderror.eta.MA
names(eta.MA)<-c("est.")
names(low)<-c("lower")
names(up)<-c("upper")
p<-c(low,eta.MA,up);
return(p)
}


#### final size Confidence interval function
CIFinalsize<-function(inc , time, start=NULL, model){
switch(model,
       Richards = {intervals(Richardsm(inc , time, start))$coef[1,]},
       logistic3P = {intervals(logistic3pm(inc , time, start))$coef[1,]},
       SigmEmax = {intervals(SigmEmaxm (inc , time, start))$coef[1,]},
       Gompertz = {intervals(gompertzm(inc , time, start))$coef[1,]},
       Weibull = {intervals(weibullm(inc , time, start))$coef[1,]},
       logistic5P = {intervals(logistic5pm (inc , time, start))$coef[1,]},
       all = {p<-rbind(intervals(Richardsm(inc , time, start[[1]]))$coef[1,],
                  intervals(logistic3pm(inc , time, start[[2]]))$coef[1,],
                  intervals(SigmEmaxm (inc , time, start[[3]]))$coef[1,],
                  intervals(gompertzm(inc , time, start[[4]]))$coef[1,],
                  intervals(weibullm(inc , time, start[[5]]))$coef[1,],
                  intervals(logistic5pm (inc , time, start[[6]]))$coef[1,],
                  estimFSizeMA(inc , time, start))
                  rownames(p)<-c("Richards","3P logistic","SigmEmax","Gompertz",
                                  "Weibull","5P logistic","Model averaged")
                  p},
       {cat("This model name is not correct\n")})
}


#### turning point Confidence interval  function
CIturningpoint<-function(inc , time, start=NULL, model){
switch(model,
       Richards = {intervals(Richardsm(inc , time, start))$coef[4,]},
       logistic3P = {intervals(logistic3pm(inc , time, start))$coef[3,]},
       SigmEmax = {intervals(SigmEmaxm (inc , time, start))$coef[2,]},
       Gompertz = {intervals(gompertzm(inc , time, start))$coef[2,]},
       Weibull = {intervals(weibullm(inc , time, start))$coef[2,]},
       logistic5P = {intervals(logistic5pm (inc , time, start))$coef[4,]},
       all = {p<-rbind(intervals(Richardsm(inc , time, start[[1]]))$coef[4,],
                  intervals(logistic3pm(inc , time, start[[2]]))$coef[3,],
                  intervals(SigmEmaxm (inc , time, start[[3]]))$coef[2,],
                  intervals(gompertzm(inc , time, start[[4]]))$coef[2,],
                  intervals(weibullm(inc , time, start[[5]]))$coef[2,],
                  intervals(logistic5pm (inc , time, start[[6]]))$coef[4,],
                  estimTpointMA(inc , time, start))
                  rownames(p)<-c("Richards","3P logistic","SigmEmax","Gompertz",
                                  "Weibull","5P logistic","Model averaged")
                  p},
       {cat("This model name is not correct\n")})
}


allmodels<-function(inc , time, start=NULL, model){
switch(model,
       Richards = {w<-Richardsm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-Richardsm(inc , time, start)$tTable
                   r2<-intervals(Richardsm(inc , time, start))$coef[1,]
                   r3<-intervals(Richardsm(inc , time, start))$coef[4,]
                   r4<-predict(Richardsm(inc , time, start))
                   pp1<-D(expression(alpha*((1+k*exp(-gamma*k*(t-eta)))^(-1/k))),"t")
                   alpha<-coef(Richardsm(inc , time, start))[1]
                   k<-coef(Richardsm(inc , time, start))[2]
                   gamma<-coef(Richardsm(inc , time, start))[3]
                   eta<-coef(Richardsm(inc , time, start))[4]
                   t<-time
                   r5<-eval(pp1)
			output = list(Incidence=inc, Time=time,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodels",model.type="Richards")
			class(output) = "dengue"	
			return(output)
                  },
      logistic3P = {w<-logistic3pm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-logistic3pm(inc , time, start)$tTable
                   r2<-intervals(logistic3pm(inc , time, start))$coef[1,]
                   r3<-intervals(logistic3pm(inc , time, start))$coef[3,]
                   r4<-predict(logistic3pm(inc , time, start))
                   pp2<-D(expression(alpha/(1+exp(-gamma*(t-eta)))),"t")
                   alpha<-coef(logistic3pm(inc , time, start))[1]
                   gamma<-coef(logistic3pm(inc , time, start))[2]
                   eta<-coef(logistic3pm(inc , time, start))[3]
                   t<-time
                   r5<-eval(pp2)
			output = list(Incidence=inc, Time=time,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodels",model.type="logistic3P")
			class(output) = "dengue"	
			return(output)
                  },
       SigmEmax = {w<-SigmEmaxm (inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-SigmEmaxm (inc , time, start)$tTable
                   r2<-intervals(SigmEmaxm (inc , time, start))$coef[1,]
                   r3<-intervals(SigmEmaxm (inc , time, start))$coef[2,]
                   r4<-predict(SigmEmaxm (inc , time, start))
                   pp3<-D(expression(e0 + (t^n)*(alpha-e0)/(t^n + eta^n)),"t")
                   alpha<-coef(SigmEmaxm (inc , time, start))[1]
                   eta<-coef(SigmEmaxm (inc , time, start))[2]
                   e0<-coef(SigmEmaxm (inc , time, start))[3]
                   n<-coef(SigmEmaxm (inc , time, start))[4]
                   t<-time
                   r5<-eval(pp3)
			output = list(Incidence=inc, Time=time,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodels",model.type="SigmEmax")
			class(output) = "dengue"	
			return(output)
                  },
       Gompertz = {w<-gompertzm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-gompertzm(inc , time, start)$tTable
                   r2<-intervals(gompertzm(inc , time, start))$coef[1,]
                   r3<-intervals(gompertzm(inc , time, start))$coef[2,]
                   r4<-predict(gompertzm(inc , time, start))
                   pp4<-D(expression(e0 + (alpha-e0)*exp(-exp(-k*(t-eta)))),"t")
                   alpha<-coef(gompertzm(inc , time, start))[1]
                   eta<-coef(gompertzm(inc , time, start))[2]
                   e0<-coef(gompertzm(inc , time, start))[3]
                   k<-coef(gompertzm(inc , time, start))[4]
                   t<-time
                   r5<-eval(pp4)
			output = list(Incidence=inc, Time=time,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodels",model.type="Gompertz")
			class(output) = "dengue"	
			return(output)
                  },
       Weibull =  {w<-weibullm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-weibullm(inc , time, start)$tTable
                   r2<-intervals(weibullm(inc , time, start))$coef[1,]
                   r3<-intervals(weibullm(inc , time, start))$coef[2,]
                   r4<-predict(weibullm(inc , time, start))
                   pp5<-D(expression(e0 + (alpha-e0)*exp(-(t/eta)^-gamma)),"t")
                   alpha<-coef(weibullm(inc , time, start))[1]
                   eta<-coef(weibullm(inc , time, start))[2]
                   e0<-coef(weibullm(inc , time, start))[3]
                   gamma<-coef(weibullm(inc , time, start))[4]
                   t<-time
                   r5<-eval(pp5)
			output = list(Incidence=inc, Time=time,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodels",model.type="Weibull")
			class(output) = "dengue"	
			return(output)
                  },
     logistic5P = {w<-logistic5pm (inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-logistic5pm (inc , time, start)$tTable
                   r2<-intervals(logistic5pm (inc , time, start))$coef[1,]
                   r3<-intervals(logistic5pm (inc , time, start))$coef[4,]
                   r4<-predict(logistic5pm (inc , time, start))
                   pp6<-D(expression(d+ (alpha-d)/(1+(2**(1/g)-1)*(t/eta)**-b)**g),"t")
                   alpha<-coef(logistic5pm (inc , time, start))[1]
                   d<-coef(logistic5pm (inc , time, start))[2]
                   g<-coef(logistic5pm (inc , time, start))[3]
                   eta<-coef(logistic5pm (inc , time, start))[4]
                   b<-coef(logistic5pm (inc , time, start))[4]
                   t<-time
                   r5<-eval(pp6)
			output = list(Incidence=inc, Time=time,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodels",model.type="logistic5P")
			class(output) = "dengue"	
			return(output)
                  },
       all = {w<-AIC(Richardsm(inc , time, start[[1]]),logistic3pm(inc , time, start[[2]]),SigmEmaxm (inc , time, start[[3]]),
                     gompertzm(inc , time, start[[4]]),weibullm(inc , time, start[[5]]),logistic5pm (inc , time, start[[6]]))[,2]
                 names(w)<-c("Richards","3P Logistic","Sigmoidal Emax","Gompertz","Weibull","5P Logistic")
               p0<-w
               p1<-weight(inc , time, start)
               p2<-CIFinalsize(inc , time, start,model)
               p3<-CIturningpoint(inc , time, start, model)
               p41<-predict(Richardsm(inc , time, start[[1]]))
               p42<-predict(logistic3pm(inc , time, start[[2]]))
               p43<-predict(SigmEmaxm (inc , time, start[[3]]))
               p44<-predict(gompertzm(inc , time, start[[4]]))
               p45<-predict(weibullm(inc , time, start[[5]]))
               p46<-predict(logistic5pm (inc , time, start[[6]]))
               p4<-list(p41,p42,p43,p44,p45,p46)
               p5<-weight(inc , time, start)[1]*p41 + weight(inc , time, start)[2]*p42 + weight(inc , time, start)[3]*p43 +weight(inc , time, start)[4]*p44 +weight(inc , time, start)[5]*p45 + weight(inc , time, start)[6]*p46
               pp1<-D(expression(alpha*((1+k*exp(-gamma*k*(t-eta)))^(-1/k))),"t")
               alpha<-coef(Richardsm(inc , time, start[1]))[1]
               k<-coef(Richardsm(inc , time, start[1]))[2]
               gamma<-coef(Richardsm(inc , time, start[1]))[3]
               eta<-coef(Richardsm(inc , time, start[1]))[4]
               t<-time
               p61<-eval(pp1)
               pp2<-D(expression(alpha/(1+exp(-gamma*(t-eta)))),"t")
               alpha<-coef(logistic3pm(inc , time, start[2]))[1]
               gamma<-coef(logistic3pm(inc , time, start[2]))[2]
               eta<-coef(logistic3pm(inc , time, start[2]))[3]
               p62<-eval(pp2)
               pp3<-D(expression(e0 + (t^n)*(alpha-e0)/(t^n + eta^n)),"t")
               alpha<-coef(SigmEmaxm (inc , time, start[3]))[1]
               eta<-coef(SigmEmaxm (inc , time, start[3]))[2]
               e0<-coef(SigmEmaxm (inc , time, start[3]))[3]
                n<-coef(SigmEmaxm (inc , time, start[3]))[4]
               p63<-eval(pp3)
               pp4<-D(expression(e0 + (alpha-e0)*exp(-exp(-k*(t-eta)))),"t")
               alpha<-coef(gompertzm(inc , time, start[4]))[1]
               eta<-coef(gompertzm(inc , time, start[4]))[2]
               e0<-coef(gompertzm(inc , time, start[4]))[3]
               k<-coef(gompertzm(inc , time, start[4]))[4]
               p64<-eval(pp4)
               pp5<-D(expression(e0 + (alpha-e0)*exp(-(t/eta)^-gamma)),"t")
               alpha<-coef(weibullm(inc , time, start[5]))[1]
               eta<-coef(weibullm(inc , time, start[5]))[2]
               e0<-coef(weibullm(inc , time, start[5]))[3]
               gamma<-coef(weibullm(inc , time, start[5]))[4]
               p65<-eval(pp5)
               pp6<-D(expression(d+ (alpha-d)/(1+(2**(1/g)-1)*(t/eta)**-b)**g),"t")
               alpha<-coef(logistic5pm (inc , time, start[6]))[1]
               d<-coef(logistic5pm (inc , time, start[6]))[2]
               g<-coef(logistic5pm (inc , time, start[6]))[3]
               eta<-coef(logistic5pm (inc , time, start[6]))[4]
               b<-coef(logistic5pm (inc , time, start[6]))[5]
               p66<-eval(pp6)
               p6<-list(p61,p62,p63,p64,p65,p66)
               p7<-weight(inc , time, start)[1]*p61 + weight(inc , time, start)[2]*p62 + weight(inc , time, start)[3]*p63 +weight(inc , time, start)[4]*p64 +weight(inc , time, start)[5]*p65 + weight(inc , time, start)[6]*p66
               output = list(Incidence=inc, Time=time,AIC=p0,Weights=p1,FinalSize=p2,TurningPoint=p3,Predict=p4,PredictMA=p5,PredInc=p6, PredMAinc=p7,function.type="allmodels",model.type="all")
               class(output) = "dengue"
		   return(output)},
     {stop("This model name is not correct")})
}


### SUMMARY FUNCTION
summary.dengue <- function(object,...){
   if(object$function.type=="allmodels"){
           if(object$model.type!="all"){
                if(object$model.type=="Richards")
                {cat("\nRichards model\n")}
                if(object$model.type=="logistic3P")
                {cat("\n3P logistic model\n")}
                if(object$model.type=="SigmEmax")
                {cat("\nSigmoid Emax model\n")}
                if(object$model.type=="Gompertz")
                {cat("\nGompertz model\n")}
                if(object$model.type=="Weibull")
                {cat("\nWeibull model\n")}
                if(object$model.type=="logistic5P")
                {cat("\n5P logistic model\n")}
                cat("\nAIC")
                print(object$AIC)
                cat("\nParameter estimate\n")
                print(object$tTable)
                cat("\nFinal size estimate:\n")
                print(object$FinalSize)
                cat("\nTurning point estimate:\n")
                print(object$TurningPoint)}
          else {cat("\nAIC\n")
                print(object$AIC)
                cat("\nModel weights\n")
                print(object$Weights)
                cat("\nModel specific and model average estimate of the final size \n")
                print(object$FinalSize)
                cat("\nModel specific and model average estimate of the turning point \n")
                print(object$TurningPoint)
                }
          }
       if(object$function.type=="allmodelpredict"){
            if(object$model.type!="all"){
                if(object$model.type=="Richards")
                {cat("\nRichards model\n")}
                if(object$model.type=="logistic3P")
                {cat("\n3P logistic model\n")}
                if(object$model.type=="SigmEmax")
                {cat("\nSigmoid Emax model\n")}
                if(object$model.type=="Gompertz")
                {cat("\nGompertz model\n")}
                if(object$model.type=="Weibull")
                {cat("\nWeibull model\n")}
                if(object$model.type=="logistic5P")
                {cat("\n5P logistic model\n")}
                cat("\nAIC")
                print(object$AIC)
                cat("\nPrediction of the final size at the end of epidemic:\n")
                print(object$FinalSize)
                cat("\nPrediction of turning point at the end of epidemic:\n")
                print(object$TurningPoint)
                cat("\nPrediction of the incidence at the time point ",length(object$PredTime)," \n")
                print(object$PredInc[[length(object$PredTime)]])
                cat("\nPrediction of the cummulative number of cases at the time point ",length(object$PredTime)," \n")
                print(object$Predict[[length(object$PredTime)]])}
          else {cat("\nAIC\n")
                print(object$AIC)
                cat("\nModel weights\n")
                print(object$Weights)
                cat("\nModel specific and model average prediction of the final size at the end of epidemic\n")
                print(object$FinalSize)
                cat("\nModel specific and model average prediction of the turning point at the end of epidemic\n")
                print(object$TurningPoint)
                cat("\nModel averaged prediction of the incidence at the time point ",length(object$PredTime)," \n")
                print(object$PredMAinc[[length(object$PredTime)]])
                cat("\nModel averaged prediction of the cummulative number of cases at the time point ",length(object$PredTime)," \n")
                print(object$PredictMA[[length(object$PredTime)]])}
          }
       if(object$function.type=="changetimeFSTP"){
             cat("\nChanges over time of the parameter estimates for the final size\n")
              print(object$FSchangetime)
              cat("\nChanges over time of the parameter estimates for the turning point\n")
              print(object$TPchangetime)
           }
}


#### PLOT FUNCTION
plot.dengue= function(x,which=c(1,2),xlab="",...){
             L<-NULL
             U<-NULL
             t<-as.character(c("Richards","3P logistic","Sigmoidal Emax","Gompertz","Weibull","5p logistic","Model averaged"))
             t<-factor(t,levels=c("Richards","3P logistic","Sigmoidal Emax","Gompertz","Weibull","5p logistic","Model averaged"))
     if( 1 %in% which){
         if(x$function.type=="changetimeFSTP"){
           stop("Cannot produce plot with option which=1, because you should use allmodels or allmodelpredict objects",call.=FALSE)
         }
         if(x$function.type=="allmodels"){
             dev.new()
		        if(x$model.type == "Richards"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict,col=2,lwd=2,lty=4)
              legend("right",legend=c("Real data","Richards model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,2), inset=0.01,cex=0.65)
		        }
            if(x$model.type == "logistic3P"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict,col=3,lwd=2,lty=4)
              legend("right",legend=c("Real data","3P Logistic model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,3), inset=0.01,cex=0.65)
            }
            if(x$model.type == "SigmEmax"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict,col=4,lwd=2,lty=4)
              legend("right",legend=c("Real data","Sigmoid Emax model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,4), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Gompertz"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict,col=5,lwd=2,lty=4)
              legend("right",legend=c("Real data","Gompertz model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,5), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Weibull"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict,col=6,lwd=2,lty=4)
              legend("right",legend=c("Real data","Weibull model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,6), inset=0.01,cex=0.65)
            }
            if(x$model.type == "logistic5P"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict,col=7,lwd=2,lty=4)
              legend("right",legend=c("Real data","5P Logistic model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,7), inset=0.01,cex=0.65)
            }
            if(x$model.type == "all"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(cumsum(x$Incidence))+10),...)
              lines(x$Time,x$Predict[[1]],col=2,lwd=2,lty=4)
              lines(x$Time,x$Predict[[2]],col=3,lwd=2,lty=4)
              lines(x$Time,x$Predict[[3]],col=4,lwd=2,lty=4)
              lines(x$Time,x$Predict[[4]],col=5,lwd=2,lty=4)
              lines(x$Time,x$Predict[[5]],col=6,lwd=2,lty=4)
              lines(x$Time,x$Predict[[6]],col=7,lwd=2,lty=4)
              lines(x$Time,x$PredictMA,col=8,lwd=3,lty=1)
              legend("right",legend=c("Real data","Richards","3P logistic","Sigmoidal Emax","4P Gompertz","Weibull","5P logistic","Model averaging"),pch=c(19,-1,-1,-1,-1,-1,-1,-1),lwd=c(-1,2,2,2,2,2,2,2), lty=c(-1,4,4,4,4,4,4,1),col=c(1,2,3,4,5,6,7,8),inset=0.01,cex=0.60)
            }
          }
         if(x$function.type=="allmodelpredict"){
            dev.new()
            if(x$model.type == "Richards"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[1]],cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict,col=2,lwd=2,lty=4)
              legend("right",legend=c("Real data","Richards model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,2), inset=0.01,cex=0.65)
            }
            if(x$model.type == "logistic3P"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[2]],cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict,col=3,lwd=2,lty=4)
              legend("right",legend=c("Real data","3P Logistic model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,3), inset=0.01,cex=0.65)
            }
            if(x$model.type == "SigmEmax"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[3]],cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict,col=4,lwd=2,lty=4)
              legend("right",legend=c("Real data","Sigmoid Emax model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,4), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Gompertz"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[4]],cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict,col=5,lwd=2,lty=4)
              legend("right",legend=c("Real data","Gompertz model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,5), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Weibull"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[5]],cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict,col=6,lwd=2,lty=4)
              legend("right",legend=c("Real data","Weibull model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,6), inset=0.01,cex=0.65)
            }
            if(x$model.type == "logistic5P"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[6]],cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict,col=7,lwd=2,lty=4)
              legend("right",legend=c("Real data","5P Logistic model"),pch=c(19,-1),lwd=c(-1,2), lty=c(-1,4),col=c(1,7), inset=0.01,cex=0.65)
            }
            if(x$model.type == "all"){
              plot(x$Time,cumsum(x$Incidence),pch=19,xlab=xlab,ylab="Cumulative number of cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$Predict[[1]],x$Predict[[2]],x$Predict[[3]],x$Predict[[4]],x$Predict[[5]],x$Predict[[6]],x$PredictMA,cumsum(x$Incidence))+10),...)
              lines(x$PredTime,x$Predict[[1]],col=2,lwd=2,lty=4)
              lines(x$PredTime,x$Predict[[2]],col=3,lwd=2,lty=4)
              lines(x$PredTime,x$Predict[[3]],col=4,lwd=2,lty=4)
              lines(x$PredTime,x$Predict[[4]],col=5,lwd=2,lty=4)
              lines(x$PredTime,x$Predict[[5]],col=6,lwd=2,lty=4)
              lines(x$PredTime,x$Predict[[6]],col=7,lwd=2,lty=4)
              lines(x$PredTime,x$PredictMA,col=8,lwd=3,lty=1)
              legend("right",legend=c("Real data","Richards","3P logistic","Sigmoidal Emax","4P Gompertz","Weibull","5P logistic","Model averaging"),pch=c(19,-1,-1,-1,-1,-1,-1,-1),lwd=c(-1,2,2,2,2,2,2,2), lty=c(-1,4,4,4,4,4,4,1),col=c(1,2,3,4,5,6,7,8),inset=0.01,cex=0.60) 
          }      
     }
     }
       if( 2 %in% which){
          if(x$function.type=="changetimeFSTP"){
            stop("Cannot produce plot with option which=2, because you should use allmodels or allmodelpredict objects",call.=FALSE) 
          }
           if(x$function.type=="allmodels"){
             dev.new()
		        if(x$model.type == "Richards"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc,col=2,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Richards model"),lwd=c(6,2), lty=c(1,4),col=c(1,2), inset=0.01,cex=0.65)
		        }
            if(x$model.type == "logistic3P"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc,col=3,lwd=2,lty=4)
              legend("topright",legend=c("Real data","3P Logistic model"),lwd=c(6,2), lty=c(1,4),col=c(1,3), inset=0.01,cex=0.65)
            }
            if(x$model.type == "SigmEmax"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc,col=4,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Sigmoid Emax model"),lwd=c(6,2), lty=c(1,4),col=c(1,4), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Gompertz"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc,col=5,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Gompertz model"),lwd=c(6,2), lty=c(1,4),col=c(1,5), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Weibull"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc,col=6,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Weibull model"),lwd=c(6,2), lty=c(1,4),col=c(1,6), inset=0.01,cex=0.65)
            }
            if(x$model.type == "logistic5P"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc,col=7,lwd=2,lty=4)
              legend("topright",legend=c("Real data","5P Logistic model"),lwd=c(6,2), lty=c(1,4),col=c(1,7), inset=0.01,cex=0.65)
            }
            if(x$model.type == "all"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$Incidence)+3),ylim=c(0,max(x$Incidence)+5),...)
              lines(x$Time,x$PredInc[[1]],col=2,lwd=2,lty=4)
              lines(x$Time,x$PredInc[[2]],col=3,lwd=2,lty=4)
              lines(x$Time,x$PredInc[[3]],col=4,lwd=2,lty=4)
              lines(x$Time,x$PredInc[[4]],col=5,lwd=2,lty=4)
              lines(x$Time,x$PredInc[[5]],col=6,lwd=2,lty=4)
              lines(x$Time,x$PredInc[[6]],col=7,lwd=2,lty=4)
              lines(x$Time,x$PredMAinc,col=8,lwd=3,lty=1)
              legend("topright",legend=c("Real data","Richards","3P logistic","Sigmoidal Emax","4P Gompertz","Weibull","5P logistic","Model averaging"),lwd=c(6,2,2,2,2,2,2,2), lty=c(1,4,4,4,4,4,4,1),col=c(1,2,3,4,5,6,7,8),inset=0.01,cex=0.60)
            }
          }
          if(x$function.type=="allmodelpredict"){
            dev.new()
            if(x$model.type == "Richards"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[1]],x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc,col=2,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Richards model"),lwd=c(6,2), lty=c(1,4),col=c(1,2), inset=0.01,cex=0.65)
            }
            if(x$model.type == "logistic3P"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[2]],x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc,col=3,lwd=2,lty=4)
              legend("topright",legend=c("Real data","3P Logistic model"),lwd=c(6,2), lty=c(1,4),col=c(1,3), inset=0.01,cex=0.65)
            }
            if(x$model.type == "SigmEmax"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[3]],x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc,col=4,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Sigmoid Emax model"),lwd=c(6,2), lty=c(1,4),col=c(1,4), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Gompertz"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[4]],x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc,col=5,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Gompertz model"),lwd=c(6,2), lty=c(1,4),col=c(1,5), inset=0.01,cex=0.65)
            }
            if(x$model.type == "Weibull"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[5]],x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc,col=6,lwd=2,lty=4)
              legend("topright",legend=c("Real data","Weibull model"),lwd=c(6,2), lty=c(1,4),col=c(1,6), inset=0.01,cex=0.65)
            }
            if(x$model.type == "logistic5P"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[6]],x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc,col=7,lwd=2,lty=4)
              legend("topright",legend=c("Real data","5P Logistic model"),lwd=c(6,2), lty=c(1,4),col=c(1,7), inset=0.01,cex=0.65)
            }
            if(x$model.type == "all"){
              plot(x$Time,x$Incidence,type="h",lwd=6,xlab=xlab,ylab="Reported cases",xlim=c(0,length(x$PredTime)+3),ylim=c(0,max(x$PredInc[[1]],x$PredInc[[2]],x$PredInc[[3]],x$PredInc[[4]],x$PredInc[[5]],x$PredInc[[6]],x$PredMAinc,x$Incidence)+5),...)
              lines(x$PredTime,x$PredInc[[1]],col=2,lwd=2,lty=4)
              lines(x$PredTime,x$PredInc[[2]],col=3,lwd=2,lty=4)
              lines(x$PredTime,x$PredInc[[3]],col=4,lwd=2,lty=4)
              lines(x$PredTime,x$PredInc[[4]],col=5,lwd=2,lty=4)
              lines(x$PredTime,x$PredInc[[5]],col=6,lwd=2,lty=4)
              lines(x$PredTime,x$PredInc[[6]],col=7,lwd=2,lty=4)
              lines(x$PredTime,x$PredMAinc,col=8,lwd=3,lty=1)
              legend("topright",legend=c("Real data","Richards","3P logistic","Sigmoidal Emax","4P Gompertz","Weibull","5P logistic","Model averaging"),lwd=c(6,2,2,2,2,2,2,2), lty=c(1,4,4,4,4,4,4,1),col=c(1,2,3,4,5,6,7,8),inset=0.01,cex=0.60)
            }
         }
       } 
     if( 3 %in% which){
              if(x$model.type!="all") {
                stop("Cannot produce plot with option which=3, because you should use all models",call.=FALSE)
		  } 
             else {
             dev.new()
		 df <- data.frame(t = t,
                 F = x$FinalSize[,2],
                 L = x$FinalSize[,1],
                 U = x$FinalSize[,3])
              print(ggplot(df, aes(x = t, y = F)) +
               geom_point(size = 4) +  
               geom_errorbar(aes(ymax = U, ymin = L))+ 
               coord_flip(xlim = c(0,8), ylim = c(min(x$FinalSize[,1])-3,max(x$FinalSize[,3])+3)) +
               xlab("Models") + ylab("Final Size"))           
             }
         }
     if( 4 %in% which){
               if(x$model.type!="all") {
                 stop("Cannot produce plot with option which=4, because you should use all models",call.=FALSE)
		     } 
               else {
		   dev.new()
               
              df<- data.frame(t =t,
                 F =x$TurningPoint[,2],
                 L =x$TurningPoint[,1],
                 U =x$TurningPoint[,3])
                print(ggplot(df, aes(x = t, y = F)) +
                 geom_point(size = 4) +  
                 geom_errorbar(aes(ymax = U, ymin = L))+ 
                  coord_flip(xlim = c(0,8), ylim = c(min(x$TurningPoint[,1])-0.1,max(x$TurningPoint[,3])+0.1)) +
                  xlab("Models") + ylab("Turning Point"))
                }
               } 
    if( 5 %in% which){
           if(x$function.type!="changetimeFSTP") {
                 stop("Cannot produce plot with option which=5, because you should use a changetimeFSTP object",call.=FALSE)
		     } 
               else {
		   dev.new()
            plot(x$Period,x$FSchangetime[,1],type="l",col=2,lwd=2,lty=4,xlab=xlab,ylab="Final size of outbreak",xlim=c(min(x$Period)-1,max(x$Period)+1),ylim=c(min(x$FSchangetime)-10,max(x$FSchangetime)+10),...)
              lines(x$Period,x$FSchangetime[,2],col=3,lwd=2,lty=4)
              lines(x$Period,x$FSchangetime[,3],col=4,lwd=2,lty=4)
              lines(x$Period,x$FSchangetime[,4],col=5,lwd=2,lty=4)
              lines(x$Period,x$FSchangetime[,5],col=6,lwd=2,lty=4)
              lines(x$Period,x$FSchangetime[,6],col=7,lwd=2,lty=4)
              lines(x$Period,x$FSchangetime[,7],col=8,lwd=3,lty=1)
		         legend("topright",legend=c("Richards","3P logistic","Sigmoidal Emax","4P Gompertz","Weibull","5P logistic","Model averaging"),lwd=c(2,2,2,2,2,2,3), lty=c(4,4,4,4,4,4,1),col=c(2,3,4,5,6,7,8),inset=0.01,cex=0.60)
               }  
         } 
    if( 6 %in% which){
           if(x$function.type!="changetimeFSTP") {
                 stop("Cannot produce plot with option which=6, because you should use a changetimeFSTP object",call.=FALSE)
		     } 
               else {
		   dev.new()
            plot(x$Period,x$TPchangetime[,1],type="l",col=2,lwd=2,lty=4,xlab=xlab,ylab="Turning point of outbreak",xlim=c(min(x$Period)-1,max(x$Period)+1),ylim=c(min(x$TPchangetime)-3,max(x$TPchangetime)+3),...)
              lines(x$Period,x$TPchangetime[,2],col=3,lwd=2,lty=4)
              lines(x$Period,x$TPchangetime[,3],col=4,lwd=2,lty=4)
              lines(x$Period,x$TPchangetime[,4],col=5,lwd=2,lty=4)
              lines(x$Period,x$TPchangetime[,5],col=6,lwd=2,lty=4)
              lines(x$Period,x$TPchangetime[,6],col=7,lwd=2,lty=4)
              lines(x$Period,x$TPchangetime[,7],col=8,lwd=3,lty=1)
		       legend("topright",legend=c("Richards","3P logistic","Sigmoidal Emax","4P Gompertz","Weibull","5P logistic","Model averaging"),lwd=c(2,2,2,2,2,2,3), lty=c(4,4,4,4,4,4,1),col=c(2,3,4,5,6,7,8),inset=0.01,cex=0.60)
               }  
         } 
                
}

 
     
allmodelpredict<-function(inc , time, pred, start=NULL, model){
switch(model,
       Richards = {w<-Richardsm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-Richardsm(inc , time, start)$tTable
                   r2<-intervals(Richardsm(inc , time, start))$coef[1,]
                   r3<-intervals(Richardsm(inc , time, start))$coef[4,]
                   pp1<-D(expression(alpha*((1+k*exp(-gamma*k*(t-eta)))^(-1/k))),"t")
                   alpha<-coef(Richardsm(inc , time, start))[1]
                   k<-coef(Richardsm(inc , time, start))[2]
                   gamma<-coef(Richardsm(inc , time, start))[3]
                   eta<-coef(Richardsm(inc , time, start))[4]
                   t<-c(1:pred)
                   r5<-eval(pp1)
                   r4<-alpha*(1+k*exp(-gamma*k*(t-eta)))^(-1/k)
			output = list(Incidence=inc, Time=time,PredTime=t, AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodelpredict",model.type="Richards")
			class(output) = "dengue"	
			return(output)
                  },
      logistic3P = {w<-logistic3pm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-logistic3pm(inc , time, start)$tTable
                   r2<-intervals(logistic3pm(inc , time, start))$coef[1,]
                   r3<-intervals(logistic3pm(inc , time, start))$coef[3,]
                   pp2<-D(expression(alpha/(1+exp(-gamma*(t-eta)))),"t")
                   alpha<-coef(logistic3pm(inc , time, start))[1]
                   gamma<-coef(logistic3pm(inc , time, start))[2]
                   eta<-coef(logistic3pm(inc , time, start))[3]
                   t<-c(1:pred)
                   r4<-alpha/(1+exp(-gamma*(t-eta)))
                   r5<-eval(pp2)
			output = list(Incidence=inc, Time=time,PredTime=t,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodelpredict",model.type="logistic3P")
			class(output) = "dengue"	
			return(output)
                  },
       SigmEmax = {w<-SigmEmaxm (inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-SigmEmaxm (inc , time, start)$tTable
                   r2<-intervals(SigmEmaxm (inc , time, start))$coef[1,]
                   r3<-intervals(SigmEmaxm (inc , time, start))$coef[2,]
                   pp3<-D(expression(e0 + (t^n)*(alpha-e0)/(t^n + eta^n)),"t")
                   alpha<-coef(SigmEmaxm (inc , time, start))[1]
                   eta<-coef(SigmEmaxm (inc , time, start))[2]
                   e0<-coef(SigmEmaxm (inc , time, start))[3]
                   n<-coef(SigmEmaxm (inc , time, start))[4]
                   t<-c(1:pred)
                   r5<-eval(pp3)
                   r4<-e0 + (t^n)*(alpha-e0)/(t^n + eta^n)
			output = list(Incidence=inc, Time=time,PredTime=t,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodelpredict",model.type="SigmEmax")
			class(output) = "dengue"	
			return(output)
                  },
       Gompertz = {w<-gompertzm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-gompertzm(inc , time, start)$tTable
                   r2<-intervals(gompertzm(inc , time, start))$coef[1,]
                   r3<-intervals(gompertzm(inc , time, start))$coef[2,]
                   pp4<-D(expression(e0 + (alpha-e0)*exp(-exp(-k*(t-eta)))),"t")
                   alpha<-coef(gompertzm(inc , time, start))[1]
                   eta<-coef(gompertzm(inc , time, start))[2]
                   e0<-coef(gompertzm(inc , time, start))[3]
                   k<-coef(gompertzm(inc , time, start))[4]
                   t<-c(1:pred)
                   r5<-eval(pp4)
                   r4<-e0 + (alpha-e0)*exp(-exp(-k*(t-eta)))
			output = list(Incidence=inc, Time=time,PredTime=t,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodelpredict",model.type="Gompertz")
			class(output) = "dengue"	
			return(output)
                  },
       Weibull =  {w<-weibullm(inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-weibullm(inc , time, start)$tTable
                   r2<-intervals(weibullm(inc , time, start))$coef[1,]
                   r3<-intervals(weibullm(inc , time, start))$coef[2,]
                   pp5<-D(expression(e0 + (alpha-e0)*exp(-(t/eta)^-gamma)),"t")
                   alpha<-coef(weibullm(inc , time, start))[1]
                   eta<-coef(weibullm(inc , time, start))[2]
                   e0<-coef(weibullm(inc , time, start))[3]
                   gamma<-coef(weibullm(inc , time, start))[4]
                   t<-c(1:pred)
                   r5<-eval(pp5)
                   r4<-e0 + (alpha-e0)*exp(-(t/eta)^-gamma)
			output = list(Incidence=inc, Time=time,PredTime=t,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodelpredict",model.type="Weibull")
			class(output) = "dengue"	
			return(output)
                  },
     logistic5P = {w<-logistic5pm (inc , time, start)$AIC
                   names(w)<-c("")
                   r0<-w
                   r1<-logistic5pm (inc , time, start)$tTable
                   r2<-intervals(logistic5pm (inc , time, start))$coef[1,]
                   r3<-intervals(logistic5pm (inc , time, start))$coef[4,]
                   pp6<-D(expression(d+ (alpha-d)/(1+(2**(1/g)-1)*(t/eta)**-b)**g),"t")
                   alpha<-coef(logistic5pm (inc , time, start))[1]
                   d<-coef(logistic5pm (inc , time, start))[2]
                   g<-coef(logistic5pm (inc , time, start))[3]
                   eta<-coef(logistic5pm (inc , time, start))[4]
                   b<-coef(logistic5pm (inc , time, start))[4]
                   t<-c(1:pred)
                   r5<-eval(pp6)
                   r4<-d+ (alpha-d)/(1+(2**(1/g)-1)*(t/eta)**-b)**g
			output = list(Incidence=inc, Time=time,PredTime=t,AIC=r0,tTable=r1,FinalSize=r2,TurningPoint=r3,Predict=r4,PredInc=r5,function.type="allmodelpredict",model.type="logistic5P")
			class(output) = "dengue"	
			return(output)
                  },
       all = {w<-AIC(Richardsm(inc , time, start[[1]]),logistic3pm(inc , time, start[[2]]),SigmEmaxm (inc , time, start[[3]]),
                     gompertzm(inc , time, start[[4]]),weibullm(inc , time, start[[5]]),logistic5pm (inc , time, start[[6]]))[,2]
                 names(w)<-c("Richards","3P Logistic","Sigmoidal Emax","Gompertz","Weibull","5P Logistic")
               p0<-w
               p1<-weight(inc , time, start)
               p2<-CIFinalsize(inc , time, start,model)
               p3<-CIturningpoint(inc , time, start, model)
               pp1<-D(expression(alpha*((1+k*exp(-gamma*k*(t-eta)))^(-1/k))),"t")
               alpha<-coef(Richardsm(inc , time, start[1]))[1]
               k<-coef(Richardsm(inc , time, start[1]))[2]
               gamma<-coef(Richardsm(inc , time, start[1]))[3]
               eta<-coef(Richardsm(inc , time, start[1]))[4]
               t<-c(1:pred)
               p61<-eval(pp1)
               p41<-alpha*((1+k*exp(-gamma*k*(t-eta)))^(-1/k))
               pp2<-D(expression(alpha/(1+exp(-gamma*(t-eta)))),"t")
               alpha<-coef(logistic3pm(inc , time, start[2]))[1]
               gamma<-coef(logistic3pm(inc , time, start[2]))[2]
               eta<-coef(logistic3pm(inc , time, start[2]))[3]
               p62<-eval(pp2)
               p42<-alpha/(1+exp(-gamma*(t-eta)))
               pp3<-D(expression(e0 + (t^n)*(alpha-e0)/(t^n + eta^n)),"t")
               alpha<-coef(SigmEmaxm (inc , time, start[3]))[1]
               eta<-coef(SigmEmaxm (inc , time, start[3]))[2]
               e0<-coef(SigmEmaxm (inc , time, start[3]))[3]
                n<-coef(SigmEmaxm (inc , time, start[3]))[4]
               p63<-eval(pp3)
               p43<-e0 + (t^n)*(alpha-e0)/(t^n + eta^n)
               pp4<-D(expression(e0 + (alpha-e0)*exp(-exp(-k*(t-eta)))),"t")
               alpha<-coef(gompertzm(inc , time, start[4]))[1]
               eta<-coef(gompertzm(inc , time, start[4]))[2]
               e0<-coef(gompertzm(inc , time, start[4]))[3]
               k<-coef(gompertzm(inc , time, start[4]))[4]
               p64<-eval(pp4)
               p44<-e0 + (alpha-e0)*exp(-exp(-k*(t-eta)))
               pp5<-D(expression(e0 + (alpha-e0)*exp(-(t/eta)^-gamma)),"t")
               alpha<-coef(weibullm(inc , time, start[5]))[1]
               eta<-coef(weibullm(inc , time, start[5]))[2]
               e0<-coef(weibullm(inc , time, start[5]))[3]
               gamma<-coef(weibullm(inc , time, start[5]))[4]
               p65<-eval(pp5)
               p45<-e0 + (alpha-e0)*exp(-(t/eta)^-gamma)
               pp6<-D(expression(d+ (alpha-d)/(1+(2**(1/g)-1)*(t/eta)**-b)**g),"t")
               alpha<-coef(logistic5pm (inc , time, start[6]))[1]
               d<-coef(logistic5pm (inc , time, start[6]))[2]
               g<-coef(logistic5pm (inc , time, start[6]))[3]
               eta<-coef(logistic5pm (inc , time, start[6]))[4]
               b<-coef(logistic5pm (inc , time, start[6]))[5]
               p66<-eval(pp6)
               p46<-d+ (alpha-d)/(1+(2**(1/g)-1)*(t/eta)**-b)**g
               p6<-list(p61,p62,p63,p64,p65,p66)
               p7<-weight(inc , time, start)[1]*p61 + weight(inc , time, start)[2]*p62 + weight(inc , time, start)[3]*p63 +weight(inc , time, start)[4]*p64 +weight(inc , time, start)[5]*p65 + weight(inc , time, start)[6]*p66
               p4<-list(p41,p42,p43,p44,p45,p46)
               p5<-weight(inc , time, start)[1]*p41 + weight(inc , time, start)[2]*p42 + weight(inc , time, start)[3]*p43 +weight(inc , time, start)[4]*p44 +weight(inc , time, start)[5]*p45 + weight(inc , time, start)[6]*p46
               output = list(Incidence=inc, Time=time,PredTime=t,AIC=p0,Weights=p1,FinalSize=p2,TurningPoint=p3,Predict=p4,PredictMA=p5,PredInc=p6, PredMAinc=p7,function.type="allmodelpredict",model.type="all")
               class(output) = "dengue"
		   return(output)},
     {stop("This model name is not correct")})
}


changetimeFSTP<-function(inc, time,ini,start=NULL){
           y1<-NULL  
           y2<-NULL
           y3<-NULL
           y4<-NULL
           y5<-NULL
           y6<-NULL
           y7<-NULL 
           t1<-NULL  
           t2<-NULL
           t3<-NULL
           t4<-NULL
           t5<-NULL
           t6<-NULL
           t7<-NULL  
           p<-NULL
           for (i in ini:length(time))
           {p[[i]]<-allmodels(inc[1:i],time[1:i], start=start, model="all")
           y1[[i-ini+1]]<-p[[i]]$FinalSize[1,2]
           y2[[i-ini+1]]<-p[[i]]$FinalSize[2,2]
           y3[[i-ini+1]]<-p[[i]]$FinalSize[3,2]
           y4[[i-ini+1]]<-p[[i]]$FinalSize[4,2]
           y5[[i-ini+1]]<-p[[i]]$FinalSize[5,2]
           y6[[i-ini+1]]<-p[[i]]$FinalSize[6,2]
           y7[[i-ini+1]]<-p[[i]]$FinalSize[7,2]
           t1[[i-ini+1]]<-p[[i]]$TurningPoint[1,2]
           t2[[i-ini+1]]<-p[[i]]$TurningPoint[2,2]
           t3[[i-ini+1]]<-p[[i]]$TurningPoint[3,2]
           t4[[i-ini+1]]<-p[[i]]$TurningPoint[4,2]
           t5[[i-ini+1]]<-p[[i]]$TurningPoint[5,2]
           t6[[i-ini+1]]<-p[[i]]$TurningPoint[6,2]
           t7[[i-ini+1]]<-p[[i]]$TurningPoint[7,2]}
           y<-c(y1,y2,y3,y4,y5,y6,y7)
           t<-c(t1,t2,t3,t4,t5,t6,t7)
           d<-matrix(y,nrow=length(time)-ini+1, ncol=7, byrow=FALSE, dimnames=list(paste("1",1:length(time),sep="-")[ini:length(time)],c("Richards","logistic3P","SigmEmax","Gompertz","Weibull","logistic5P","Model averaged")))
           tt<-matrix(t,nrow=length(time)-ini+1, ncol=7, byrow=FALSE, dimnames=list(paste("1",1:length(time),sep="-")[ini:length(time)],c("Richards","logistic3P","SigmEmax","Gompertz","Weibull","logistic5P","Model averaged")))
           g<-c(ini:length(time))
           output = list(Incidence=inc, Time=time,Period=g,FSchangetime=d,TPchangetime=tt,function.type="changetimeFSTP",model.type="all")
           class(output) = "dengue"
	     return(output)
}

