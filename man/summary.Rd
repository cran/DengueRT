\name{summary}
\alias{summary.dengue}
\title{
Result summaries of the functions
}
\description{
summary is a generic function used to produce result summaries of the functions: allmodels, allmodelpredict and changetimeFSTP from the S3 class dengue. The function invokes particular methods of the S3 class dengue. 
}
\usage{
\method{summary}{dengue}(object, \dots)
}

\arguments{
  \item{object}{
An object for which a summary is desired. It could be any of the following objects: allmodels, allmodelpredict and changetimeFSTP.
}
\item{\dots}{
Additional arguments affecting the summary produced.
}
}

\value{
The form of the value returned by summary depends on the object of the S3 class dengue. In the case of allmodels returns the AIC, model averaged Weights and model averaged estimates for the final size and turning point of outbreak. If the object is allmodelpredict to the returned value is added a prediction of the incidence and cumulative cases in the required future observation. When the object is changetimeFSTP function the returned values are two tables with the changes over time of the final size and turning point for the time period required.
}
\references{
K. Burnham, D. R. Anderson, Model Selection and Multimodel Inference: A Practical Information-theoretic Approach,
2nd Edition, Springer-Verlag, New York, 2002.\cr
G. Claeskens, N. L. Hjort, Model selection and model averaging, Cambridge University Press, 2008.

}
\author{
Carlos Sebrango, Lizet Sanchez, Ziv Shkedy, Ewoud DE Troyer
}

\examples{
\dontrun{
## (data example 1)
data("dengueoutbreak1")

## summary function for a allmodels object 
##Richards model (Incidence data example 1)
ee<-allmodels(dengueoutbreak1$Incidence,dengueoutbreak1$Time, model="Richards")
summary(ee)

## (data example 2)
data("dengueoutbreak2")
##using all built-in models (Incidence data example 2)
ee2<-allmodels(dengueoutbreak2$Incidence,dengueoutbreak2$Time, model="all")
summary(ee2)

##summary function for allmodelpredict object
##3P logistic model (Incidence data example 2)
## information avalaible to the time point 20
## observation to predict: 30 
ee3<-allmodelpredict(dengueoutbreak2$Incidence[1:20],dengueoutbreak2$Time[1:20],30, 
model="logistic3P")
summary(ee3)

##using all built-in models (Incidence data example 1)
## information avalaible to the time point 22
## observation to predict: 30 
ee4<-allmodelpredict(dengueoutbreak1$Incidence[1:22],dengueoutbreak1$Time[1:22],30, 
model="all")
summary(ee4)}

##summary function for changetimeFSTP object
## See example changetimeFSTP function 



}


