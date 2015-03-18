\name{rdplot}
\alias{rdplot}
\alias{print.rdplot}
\alias{summary.rdplot}
\title{Data-Driven Regression Discontinuity Plots}
\description{
 \code{\link{rdplot}} implements several data-driven RD plots, using either evenly-spaced or quantile-spaced partitioning.  Two type of RD plots are constructed:  (i) RD plots with binned sample means tracing out the underlying regression function (IMSE-optimal selectors), and (ii) RD plots with binned sample means mimicking the underlying variability of the data. For all technical and methodological details see Calonico, Cattaneo and Titiunik (2015a).
 
    Calonico, Cattaneo and Titiunik (2015b) provides an introduction to this command. For a review on RD methods see Imbens and Lemieux (2008), Lee and Lemieux (2010), Calonico, Cattaneo and Titiunik (2014a), and references therein.
    
    For local polynomial inference methods see the commands \code{\link{rdrobust}} and \code{\link{rdbwselect}}.

    A companion \code{Stata} package is described in Calonico, Cattaneo and Titiunik (2014b).
}
\usage{
rdplot(y, x, subset = NULL, c = 0, p = 4, 
          numbinl = NULL, numbinr = NULL, binselect = "esmv", 
          lowerend = NULL, upperend = NULL, scale = 1, scalel = 1, scaler = 1, 
          hide = FALSE, par=NULL, title = NULL, x.label = NULL, y.label = NULL, 
          x.lim = NULL, y.lim = NULL, col.dots = NULL, col.lines = NULL, 
          type.dots = NULL,...)
}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the running variable (a.k.a. score or forcing variable).}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{c}{specifies the RD cutoff in \code{x}; default is \code{c = 0}.}
  \item{p}{specifies the order of the global-polynomial used to approximate the population conditional mean functions for control and treated units; default is \code{p = 4}.} 
  \item{numbinl}{specifies the number of bins used to the left of the cutoff, denoted \eqn{J_-}. If not specified, \eqn{J_-} is estimated using the method and options chosen below.}
    \item{numbinr}{specifies the number of bins used to the right of the cutoff, denoted \eqn{J_+}. If not specified, \eqn{J_+} is estimated using the method and options chosen below.}
\item{binselect}{specifies the procedure to select the number of bins. This option is available only if \eqn{J_-} and \eqn{J_+} are not set manually. Options are:  

        \code{es}: IMSE-optimal evenly-spaced method using spacings estimators.
        
        \code{espr}: IMSE-optimal evenly-spaced method using polynomial regression.
        
        \code{esmv}: mimicking variance evenly-spaced method using spacings estimators. This is the default option.
        
        \code{esmvpr}: mimicking variance evenly-spaced method using polynomial regression.
        
        \code{qs}: IMSE-optimal quantile-spaced method using spacings estimators.
        
        \code{qspr}: IMSE-optimal quantile-spaced method using polynomial regression.
        
        \code{qsmv}: mimicking variance quantile-spaced method using spacings estimators.
        
        \code{qsmvpr}: mimicking variance quantile-spaced method using polynomial regression.
     }
    \item{lowerend}{specifies the lower bound for \code{x} to the left of the cutoff. Default is the minimum value in sample.}
    \item{upperend}{specifies the upper bound for \code{x} to the right of the cutoff. Default is the maximum value in sample.}
    \item{scale}{specifies a multiplicative factor to be used with the optimal numbers of bins selected. Specifically, the number of bins used for the treatment and control groups will be \code{scale}\eqn{\times \hat{J}_+} and \code{scale}\eqn{\times \hat{J}_-}, where \eqn{\hat{J}_\cdot} denotes the estimated optimal numbers of bins originally computed for each group; default is \code{scale = 1}.}
  \item{scalel}{same as scale but only for the data to the left of the cutoff.}
  \item{scaler}{same as scale but only for the data to the right of the cutoff.}
  \item{hide}{logical. If \code{TRUE}, it omits the RD plot; default is \code{hide = FALSE}.}
  \item{par}{optional graphical parameters.}
  \item{title}{optional title for the RD plot.}
  \item{x.label}{optional label for the x-axis of the RD plot.}
  \item{y.label}{optional label for the y-axis of the RD plot.}
  \item{x.lim}{optional setting for the range of the x-axis in the RD plot.}
  \item{y.lim}{optional setting for the range of the y-axis in the RD plot.}
  \item{col.dots}{optional setting for the color of the dots in the RD plot.}
  \item{col.lines}{optional setting for the color of the lines in the RD plot.}
  \item{type.dots}{optional setting for the type of lines in the RD plot.}
  \item{...}{other graphical parameters.}
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

    Imbens, G., and T. Lemieux. 2008. "Regression Discontinuity Designs: A Guide to Practice." \emph{Journal of Econometrics} 142(2): 615-635.

    Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in Economics." \emph{Journal of Economic Literature} 48(2): 281-355.
    
        
}
\author{
Sebastian Calonico, University of Miami, Coral Gables, FL.  \email{scalonico@bus.miami.edu}.

Matias D. Cattaneo, University of Michigan, Ann Arbor, MI.  \email{cattaneo@umich.edu}.

Rocio Titiunik, University of Michigan, Ann Arbor, MI.  \email{titiunik@umich.edu}.
}
\examples{
x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdplot(y,x)
}
\keyword{ regression discontinuity }
\keyword{ RD plots }
\keyword{ binning }
\keyword{ partitioning}
\keyword{ tuning parameter selection}
\seealso{
  \code{\link{rdbwselect}}, \code{\link{rdrobust}}
}