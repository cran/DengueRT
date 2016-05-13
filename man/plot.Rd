\name{plot}
\alias{plot.dengue}
\title{
Plots for a visualization of the results
}
\description{
plot is a generic function used to produce plots for a visualization of the results  of the functions: allmodels, allmodelpredict and changetimeFSTP from the S3 class dengue. The function invokes particular methods of the S3 class dengue. 
}
\usage{
\method{plot}{dengue}(x, which = c(1, 2), xlab = "", \dots)
}
\arguments{
  \item{x}{
An object for which a plot is desired. It could be any of the following objects: allmodels, allmodelpredict and changetimeFSTP
}
  \item{which}{
Different graph options, from 1 to 6. 
\itemize{
\item Option 1 gives a plot with the cumulative number of cases and fitted models for each built-in model and model averaged. This option is only available for allmodels and allmodelpredict objects.
\item Option 2 gives a plot with the incidence and fitted models for each built-in model and model averaged. This option is only available for allmodels and allmodelpredict objects.
\item Option 3 gives a plot with the final size 95\% confidence interval and point estimate for each built-in model and model averaged. This option is only available for allmodels and allmodelpredict objects when all models are used.
\item Option 4 gives a plot with the turning point 95\% confidence interval and point estimate for each built-in model and model averaged. This option is only available for allmodels and allmodelpredict objects when all models are used.
\item Option 5 gives a plot with the change over time of the final size estimates for each built-in model and model averaged. This option is only available for changetimeFSTP objects.
\item Option 6 gives a plot with the change over time of the turning point estimates for each built-in model and model averaged. This option is only available  for changetimeFSTP objects.
}

}
  \item{xlab}{
A title for the x axis. It is not necesarry for options 3 y 4.
}
  \item{\dots}{
Additional arguments affecting the plot produced. Arguments to be passed to methods, such as graphical 
parameters (lines, segment, etc) 
}
}
\details{
When the which options 1, 2, 5 and 6 are used, the plot displays all fitted models with the following colours: (Richards (red), 3P logistic (green), Sigmoid Emax (dark blue), Gompertz (clear blue), Weibull (pink) and 5P logistic (yelow))and the fitted model via model averaging (gray colour). \cr
The which options 1, 2, 3 and 4 are only used with the allmodels and allmodelpredict objects. The which options 5 and 6 are only used with the changetimeFSTP object.
}

\references{
K. Burnham, D. R. Anderson, Model Selection and Multimodel Inference: A Practical Information-theoretic Approach,
2nd Edition, Springer-Verlag, New York, 2002.\cr
G. Claeskens, N. L. Hjort, Model selection and model averaging, Cambridge University Press, 2008.
}
\author{
Carlos Sebrango, Lizet Sanchez, Ziv Shkedy
}

\examples{
\dontrun{
## (data example 1)
data("dengueoutbreak1")

## plot function for a allmodels object 
##Richards model 
ee<-allmodels(dengueoutbreak1$Incidence,dengueoutbreak1$Time, model="Richards")
plot(ee,which=c(1,2))

## (data example 2)
data("dengueoutbreak2")
##using all built-in models 
ee2<-allmodels(dengueoutbreak2$Incidence,dengueoutbreak2$Time, model="all")
plot(ee2,which=c(1,2,3,4))

##plot function for allmodelpredict object
##3P logistic model (Incidence data example 2)
## information avalaible to the time point 20
## observation to predict: 30 
ee3<-allmodelpredict(dengueoutbreak2$Incidence[1:20],dengueoutbreak2$Time[1:20],30, 
model="logistic3P")
plot(ee3,which=c(1,2))

##using all built-in models (Incidence data example 2)
## information avalaible to the time point 22
## observation to predict: 30 
ee4<-allmodelpredict(dengueoutbreak2$Incidence[1:22],dengueoutbreak2$Time[1:22],30, 
model="all")
plot(ee4,which=c(1,2,3,4))}

##plot function for changetimeFSTP object
## See example changetimeFSTP function 

}

