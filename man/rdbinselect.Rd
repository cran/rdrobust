\name{rdbinselect}
\alias{rdbinselect}
\alias{print.rdbinselect}
\alias{summary.rdbinselect}
\title{Data-Driven Optimal Length Selector of Evenly-Spaced and Quantile-Based Bins Employed to Approximate the Underlying Regression Functions in RD Estimation}
\description{
 \code{\link{rdbinselect}} implements IMSE-optimal data-driven evenly-spaced and quantile-spaced RD plots.

    NOTE: this command is no longer mantained, and it is included only for back compatibility. Please use the new command \code{\link{rdplot}}.

    See Calonico, Cattaneo and Titiunik (2015a) and Calonico, Cattaneo and Titiunik (2014b) for a discussion on the results implemented by these commands.
    
    A companion \code{Stata} package is described in Calonico, Cattaneo and Titiunik (2014b).
}
\usage{
rdbinselect(y, x, subset = NULL, c = 0, p = 4, 
          numbinl = NULL, numbinr = NULL, binselect = "es", 
          lowerend = NULL, upperend = NULL, scale = 1, hide = FALSE, 
          par=NULL, title = NULL, x.label = NULL, y.label = NULL, 
          x.lim = NULL, y.lim = NULL)

}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the running variable (a.k.a. score or forcing variable).}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{c}{specifies the RD cutoff in \code{x}; default is \code{c = 0}.}
  \item{p}{specifies the order of the global-polynomial used to approximate the population conditional mean functions for control and treated units; default is \code{p = 4}.} 
  \item{numbinl}{specifies the number of bins used to the left of the cutoff, denoted \eqn{J_-}. If not specified, \eqn{J_-} is estimated using the method and options chosen below.}
    \item{numbinr}{specifies the number of bins used to the right of the cutoff, denoted \eqn{J_+}. If not specified, \eqn{J_+} is estimated using the method and options chosen below.}
\item{binselect}{specifies the partition length selection procedure to be used. This option is available only if \eqn{J_-} and \eqn{J_+} are not set manually. Options are:  

      \code{es}: evenly-spaced method using spacings estimators; default option. 
      
      \code{espr}: evenly-spaced method with polynomial regression. 
      
      \code{esdw}:  density-weighted evenly-spaced method. 
      
      \code{qs}: quantile-spaced method with spacings estimators. 
      
      \code{qspr}: quantile-spaced method using polynomial regression. 
      
      \code{qsdw}:  density-weighted quantile-spaced method.}
    \item{lowerend}{specifies the lower bound for \code{x} to the left of the cutoff. Default is the minimum value in sample.}
    \item{upperend}{specifies the upper bound for \code{x} to the right of the cutoff. Default is the maximum value in sample.}
    \item{scale}{specifies a multiplicative factor to be used with the optimal numbers of bins selected. Specifically, the number of bins used for the treatment and control groups will be \code{scale}\eqn{\times \hat{J}_+} and \code{scale}\eqn{\times \hat{J}_-}, where \eqn{\hat{J}_\cdot} denotes the estimated optimal numbers of bins originally computed for each group; default is \code{scale = 1}.}
  \item{hide}{logical. If \code{TRUE}, it omits the RD plot; default is \code{hide = FALSE}.}
  \item{par}{optional graphical parameters.}
  \item{title}{optional title for the RD plot.}
  \item{x.label}{optional label for the x-axis of the RD plot.}
  \item{y.label}{optional label for the y-axis of the RD plot.}
  \item{x.lim}{optional setting for the range of the x-axis in the RD plot.}
  \item{y.lim}{optional setting for the range of the x-axis in the RD plot.}
}
\value{
  \item{method}{method used to compute the optimal number of bins.}
  \item{results}{matrix containing the main estimation results (sample size, polynomial order, number of bins, scales and bin lengths).}
  \item{coef}{matrix containing the coefficients of the \eqn{p^{th}} order global polynomial estimated both sides of the cutoff.}
}
\references{
        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014a. "Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs." \emph{Econometrica} 82(6): 2295-2326, November 2014.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}.

        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014b. "Robust Data-Driven Inference in the Regression-Discontinuity Design." \emph{Stata Journal} 14(4): 909-946, 4th Quarter 2014.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_Stata.pdf}.

        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2015a. "Optimal Data-Driven Regression Discontinuity Plots." \emph{Journal of the American Statistical Association}, forthcoming.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2015_JASA.pdf}.

        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2015b. "rdrobust: An R Package for Robust Inference in Regression-Discontinuity Designs." \emph{R Journal}, forthcoming. \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2015_R.pdf}.

}
\author{
Sebastian Calonico, University of Miami, Coral Gables, FL.  \email{scalonico@bus.miami.edu}.

Matias D. Cattaneo, University of Michigan, Ann Arbor, MI.  \email{cattaneo@umich.edu}.

Rocio Titiunik, University of Michigan, Ann Arbor, MI.  \email{titiunik@umich.edu}.
}
\examples{
x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdbinselect(y,x)
}
\keyword{ regression discontinuity }
\keyword{ RD plots }
\keyword{ binning }
\keyword{ partitioning}
\keyword{ tuning parameter selection}
\seealso{
  \code{\link{rdbwselect}}, \code{\link{rdrobust}}, \code{\link{rdplot}}
}