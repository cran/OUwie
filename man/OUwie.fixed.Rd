\name{OUwie.fixed}
\alias{OUwie.fixed}
\title{Generalized Hansen model likelihood calculator}
\description{Allows the user to calculate the likelihood given a specified set of parameter values}
\usage{
OUwie.fixed(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),simmap.tree=FALSE,
root.station=TRUE, alpha=NULL, sigma.sq=NULL, theta=NULL, clade=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing species information (see Details)}
\item{model}{models to fit to comparative data (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated.}
\item{alpha}{a numeric vector giving the values of \eqn{\alpha}{alpha} for each selective regime}
\item{sigma.sq}{a numeric vector giving the values of \eqn{\sigma^2}{sigma^2} for each selective regime} 
\item{theta}{a numeric vector giving the values of \eqn{\theta}{theta} for each selective regime}
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest.}
}
\details{
The input is a tree and a data file. The tree must be of class \dQuote{phylo} and must contain the ancestral selective regimes as internal node labels. The data file is a dataframe that must have column entries in the following order: [,1] species names and [,2] their current selective regime. The user specifies the  parameter values (i.e. \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, and \eqn{\theta}{theta}).}
\value{
\item{loglik}{the maximum log-likelihood.}
\item{AIC}{Akaike information criterion.}
\item{AICc}{Akaike information criterion corrected for sample-size.}
\item{Param.est}{a matrix containing \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{theta}{a matrix containing \eqn{\theta}{theta}. Note: if \code{model=OU1} and \code{root.station=TRUE}, then \code{theta}=\eqn{\theta_0}{theta_0}.}
\item{ahat}{the estimate of the root for models based on a Brownian motion process (i.e., \code{model=BM1} or \code{model=BMS}).}
}
\examples{
data(tworegime)

#Calculate the likelihood based on known values of
#alpha, sigma^2, and theta:
alpha=c(1.0,0.5)
sigma.sq=c(0.45,0.9)
theta=c(1.0,2.0)

OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, clade=NULL, alpha,sigma.sq,theta)
}
\author{Jeremy M. Beaulieu}
