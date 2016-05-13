\name{allmodels}
\alias{allmodels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak.
}
\description{
This function gives the parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak. Also this function, when all the built-in models are used,  gives the AIC of each model, the model averaged weights and predicted incidence and cumulative cases. This function is used retrospectively, that is, when all the data are available.
}
\usage{
## Object of the S3 class dengue
allmodels(inc,time,start=NULL,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list with the starting values of the model to be used for fitting the data. If model= "all" the imput must be a list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters. By default, the initial values are provided by self-starting functions.
\itemize{
\item 3P Logistic model uses SSlogis selfstart function from the stats package.
\item Sigmoid Emax model uses LL.4 selfstart function from the drc package.
\item Gompertz model uses SSgompertz selfstart function from the stats package.
\item Weibull model uses SSweibull selfstart function from stats package.
\item 5P Logistic model uses LL.5 and LL.4 selfstart functions from drc package.
\item Richards model uses SSlogis selfstart function and some ideas of SSRichards  selfstart function from NRAIA package now archived.
}
}
\item{model}{
The nonlinear model to be used for fitting the data. Built-in models are  "Richards", "Logistic3P", "SigmEmax", "Gompertz", "Weibull" and "Logistic5P". If model = "all" the parameter estimate will be done taking into account all built-in models via model averaging.For all built-in model, C(t) represents the cumulative number of reported cases at time t and also the turning point (eta) and the final size of epidemic (alpha) are parameters in the models. The model expressions of each built-in model are given as follow:
\itemize{
\item The Richards model is given by the expression:\cr\cr
C(t)=alpha/[1+k.exp(-k.gamma.(t-eta))]^(1/k)\cr\cr
Where \cr
gamma:  Per capita growth rate or intrinsic growth rate.\cr
k: Exponent of the deviation from standard logistic curve.
\item The three parameter logistic model (3P Logistic ), proposed by Verhults, is given by the expression:\cr\cr
C(t)=alpha/[1+exp(-gamma.(t-eta))]\cr\cr
Where \cr
gamma:  Per capita growth rate or intrinsic growth rate.\cr
\item The Sigmoid Emax model is given by the expression:\cr\cr
C(t)=beta + (t^n)(alpha-beta)/(t^n + eta^n)\cr\cr
Where \cr
n: Slope factor or Hill factor.\cr
beta: Lower asymptote.
\item The Gompertz model is given by the expression :\cr\cr
C(t)=beta+(alpha-beta).exp(-exp(-gamma(t-eta)))\cr\cr
Where \cr
gamma:  Per capita growth rate or intrinsic growth rate.\cr
beta: Lower asymptote.
\item The Weibull model is given by the expression:\cr\cr
C(t)=alpha + (beta-alpha).exp(-(t/eta)^k)\cr\cr
Where \cr
k: Shape factor.\cr
beta: Lower asymptote.
\item The five parameter logistic model (5P Logistic ) is given by the expression:\cr\cr
C(t)=alpha + (beta-alpha)/[1+(2^(1/g)-1).(t/eta)^k]^g \cr\cr
Where \cr
g: Asymmetry factor. \cr
k: Shape factor. \cr
beta: Lower asymptote.
}
}
}

\value{
 An object with the parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak. It is a list:
\item{Incidence }{All the available incidences}
\item{Time}{All the available time points}
\item{AIC}{The AIC for each built-in model and model averaged}
\item{tTable}{A table with parameter estimates and t test. It is not available when all the model are used.}
\item{Weights}{Model averaged weights. It is not availabe when is used only one built-in model}
\item{FinalSize}{95\% confidence interval and point estimate of the final size of outbreak for each built-in model and model averaged estimate}
\item{TurningPoint}{95\% confidence interval and point estimate of the turning point of outbreak for each built-in model and model averaged estimate}
\item{Predict}{Predicted cumulative cases for each built-in model}
\item{PredictMA}{Predicted cumulative cases for model averaged}
\item{PredInc}{Predicted incidence for each built-in model}
\item{PredMAInc}{Predicted incidence for model averaged}
\item{function.type}{Name of the function}
\item{model.type}{models used to estimate}
Generic functions such as plot and summary have methods to show the results of the fit
}
\references{
K. Burnham, D. R. Anderson, Model Selection and Multimodel Inference: A Practical Information-theoretic Approach,
2nd Edition, Springer-Verlag, New York, 2002.\cr
J. MacDougall, Analysis of dose responses Studies: Emax model, in: N. Ting (Ed.), Dose Finding in Drug Development, Statistics for Biology and Health, Springer New York,  pp. 127, 2006.\cr
G. Claeskens, N. L. Hjort, Model selection and model averaging, Cambridge University Press, 2008.\cr
D. Ratkowsky, Handbook of nonlinear regression models, Marcel Dekker, New York, 1990.\cr
F. Richards, A flexible growth function for empirical use, Journal of Experimental Botany 10 (29), pp 290-301, 1959.\cr
Y.H. Hsieh, Temporal trend and regional variability of 2001-2002 multiwave DENV-3 epidemic in Havana City: did Hurricane Michelle contribute to its severity?, Tropical Medicine and International Health, Vol. 18, no. 7, pp 830-838, 2013.\cr
A. Tsoularis, J. Wallace, Analysis of logistic growth models, Mathematical Biosciences, Vol. 179, no. 1, pp 21-55, 2002.\cr
J. Liao, R. Liu, Re-parameterization of five-parameter logistic function, Journal of Chemometrics 23 (5), pp 248-253, 2009.
}

\author{
Carlos Sebrango, Lizet Sanchez, Ziv Shkedy
}

\examples{
## (data example 1)
data("dengueoutbreak1")
\dontrun{
## Parameter estimate for Richards model (for Incidence data example 1)
allmodels(dengueoutbreak1$Incidence,dengueoutbreak1$Time, 
start=list(alpha=375,k=2.38,gamma=0.76,eta=16), model="Richards")}
## or
p1<-allmodels(dengueoutbreak1$Incidence,dengueoutbreak1$Time, model="Richards")
## summary function for a allmodels object 
summary(p1)
## plot function for a allmodels object
plot(p1,which=c(1,2))

\dontrun{
## Parameter estimate for each built-in model and model averaged 
##estimate for final size and turning point of outbreak 

allmodels(dengueoutbreak1$Incidence,dengueoutbreak1$Time, 
start=list(list(alpha=375,k=2.38,gamma=0.76,eta=16), 
list(alpha=375,gamma=1,eta=16),list(alpha=380,eta=13,beta=7,n=5),
list(alpha=380,eta=20,beta=0,gamma=1),list(alpha=410,eta=12,beta=11,k=3),
list(alpha=380,beta=4,g=1,eta=13,k=15)), model="all")
## or
allmodels(dengueoutbreak1$Incidence,dengueoutbreak1$Time,model="all")

## (data example 2)
data("dengueoutbreak2")
# Parameter estimate for 3P Logistic model 
allmodels(dengueoutbreak2$Incidence,dengueoutbreak2$Time, 
start=list(alpha=375,gamma=1,eta=16), model="logistic3P")
## or
allmodels(dengueoutbreak2$Incidence,dengueoutbreak2$Time,model="logistic3P")

## Parameter estimate for each built-in model and model averaged estimate 
##for final size and turning point of outbreak 
##for Incidence data example 2
allmodels(dengueoutbreak2$Incidence,dengueoutbreak2$Time, 
start=list(list(alpha=355,k=1,gamma=1,eta=14),
list(alpha=355,gamma=1,eta=14), list(alpha=355,eta=13,beta=10,n=6),
list(alpha=355,eta=11,beta=20,gamma=1),list(alpha=355,eta=12,beta=22,k=3),
list(alpha=355,beta=15,g=1,eta=13,k=10)),model="all")
## or
allmodels(dengueoutbreak2$Incidence,dengueoutbreak2$Time,model="all")}
}
