\name{changetimeFSTP}
\alias{changetimeFSTP}
\title{
Changes over time of the final size and turning point estimates for each built-in model and model averaged.
}
\description{
This function provides the changes over time of the final size and turning point estimates (for specific models and model average) from the time point
required until the last time point available. The output of the function are two tables, one for final size and the other for turning point. the rows represent the data used for the estimation of the model parameters.
}
\usage{
## Object of the S3 class dengue
changetimeFSTP(inc, time, ini, start = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
 \item{ini}{
Time point from which is required to see the changes over time of the final size and turning point estimates
}
  \item{start}{
A list with the starting values of the model to be used for fitting the data. If model="all" the imput must be a list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters. By default, the initial values are provided by self-starting functions.
}
}
\value{
 An object with the changes over time of the final size and turning point estimates for specific models and model average. It is a list:
\item{Incidence }{All the available incidences}
\item{Time}{All the available time points}
\item{Period}{Time period required to see the changes over time of the final size and turning point estimates}
\item{FSchangetime}{A Table with the changes over time of the final size for the time period required}
\item{TPchangetime}{A Table with the changes over time of the turning point for the time period required}
\item{function.type}{Name of the function}
\item{model.type}{models used to estimate}
Generic functions such as plot and summary have methods to show the results of the fit.
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
## (data example 2)
data("dengueoutbreak2")

##Time point from which is required to see the changes over time 
##of the final size and turning point estimates:20

t2<-changetimeFSTP(dengueoutbreak2$Incidence,dengueoutbreak2$Time,20)

##summary function for changetimeFSTP object
summary(t2)

## plot function for changetimeFSTP object 
plot(t2,which=c(5,6))
}
}

