\name{allmodelpredict}
\alias{allmodelpredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Real-time predictions of the final size and the turning point at the end of the epidemic, the incidence and the cumulative number of cases in future observations.
}
\description{
This function provides real-time predictions of the final size and the turning point at the end of the epidemic for each built-in model and model averaged, as well as the incidence and the cumulative number of cases in future observations. Also this function, when all the built-in models are used,  gives the AIC of each model, the model averaged weights and predicted incidence and cumulative cases.
}
\usage{
## Object of the S3 class dengue
allmodelpredict(inc,time,pred,start=NULL,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
\item{pred}{
Number of observation in which the incidence and the cumulative number of cases will be predicted.
}

  \item{start}{
A list with the starting values of the model to be used for fitting the data. If model="all" the imput must be a list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters. By default, the initial values are provided by self-starting functions (see argument start in \code{allmodels} function).
}
\item{model}{
The nonlinear model to be used for fitting the data. Built-in models are  "Richards"", "Logistic3P", "SigmEmax", "Gompertz", "Weibull" and "Logistic5P" (see argument model in \code{allmodels} function). If model = "all" the parameter estimate will be done taking into account all built-in models via model averaging.
}
}

\value{
 An object with the parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak. It is a list:
\item{Incidence }{All the available incidences}
\item{Time}{All the available time points}
\item{PredTime}{Period of time for which the prediction is required}
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
K. Burnham, D. R. Anderson, Model Selection and Multimodel Inference: A Practical Information - theoretic Approach,
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

## Using only the information until time point 20
## Prediction of the final size and turning point at the end of epidemic,  
## the incidence and the cumulative number of cases in the observation  
## number 36 using only the Richards model
allmodelpredict(dengueoutbreak1$Incidence[1:20],dengueoutbreak1$Time[1:20],36, 
model = "Richards")

## Using only the information until time point 22
## Now using all built-in model, Prediction of the final size and turning point  
## at the end of epidemic, the incidence and the cumulative number of cases 
## in the observation number 30 
allmodelpredict(dengueoutbreak1$Incidence[1:22],dengueoutbreak1$Time[1:22],30, 
model = "all")}

## (data example 2)
data("dengueoutbreak2")
## Using only the information until time point 18, Prediction of the final size and  
## turning point at the end of epidemic, the incidence and the cumulative number  
## of cases in the observation number 31 using only the 3P logistic model

allmodelpredict(dengueoutbreak2$Incidence[1:18],dengueoutbreak2$Time[1:18],31,
model = "logistic3P")

\dontrun{
## Using only the information until time point 20
## Now using all built-in model, Prediction of the final size and turning point 
## at the end of  epidemic, the incidence and the cumulative number of cases  
## in the observation number 30
allmodelpredict(dengueoutbreak2$Incidence[1:20],dengueoutbreak2$Time[1:20],34, 
model = "all")}

}