\name{rdbinselect}
\alias{rdbinselect}
\alias{print.rdbinselect}
\alias{summary.rdbinselect}
\alias{print.summary.rdbinselect}
\title{Data-driven optimal length selector of evenly-spaced and quantile-based bins employed to approximate the underlying regression functions in RD estimation}
\description{
 \code{\link{rdbinselect}} employs the results in Calonico, Cattaneo and Titiunik (2014d) to implement several novel data-driven optimal choices of evenly-spaced and quantile-based bins used to approximate the underlying regression functions by local sample means, which is then used to construct the RD plots usually found in empirical applications. 

    Calonico, Cattaneo and Titiunik (2014b) provides an introduction to this command. For a review on RD methods see Imbens and Lemieux (2008), Lee and Lemieux (2010), Dinardo and Lee (2011),  Calonico, Cattaneo and Titiunik (2014a), and references therein.
    
    A companion \code{Stata} package is described in Calonico, Cattaneo and Titiunik (2014c).
}
\usage{
rdbinselect(y, x, data, subset = NULL, c = 0, p = 4, 
          numbinl = NULL, numbinr = NULL, binselect = "es", 
          lowerend = NULL, upperend = NULL, scale = 1, hide = FALSE, 
          title = NULL, x.label = NULL, y.label = NULL, 
          x.lim = NULL, y.lim = NULL, model = FALSE, frame = FALSE)

}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the running variable (a.k.a. score or forcing variable).}
  \item{data}{an optional data frame.}
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
  \item{title}{optional title for the RD plot.}
  \item{x.label}{optional label for the x-axis of the RD plot.}
  \item{y.label}{optional label for the y-axis of the RD plot.}
  \item{x.lim}{optional setting for the range of the x-axis in the RD plot.}
  \item{y.lim}{optional setting for the range of the x-axis in the RD plot.}
  \item{frame}{logical. If \code{TRUE}, the data frame used in model fitting will be returned.}
  \item{model}{logical. If \code{TRUE}, the model object  will be returned.}
}
\value{
  \item{method}{method used to compute the optimal number of bins.}
  \item{results}{matrix containing the main estimation results (sample size, polynomial order, number of bins, scales and bin lengths).}
  \item{coef}{matrix containing the coefficients of the \eqn{p^{th}} order global polynomial estimated both sides of the cutoff.}
}
\references{
        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014a. "Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs." \emph{Econometrica}, forthcoming.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}.
        
      Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014b. "rdrobust: An R Package for Robust Inference in Regression-Discontinuity Designs." University of Michigan, Department of Economics. \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_JSS.pdf}.

      Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014c. "Robust Data-Driven Inference in the Regression-Discontinuity Design." University of Michigan, Department of Economics. \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_Stata.pdf}.

      Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014d. "Optimal Data-Driven Regression Discontinuity Plots." University of Michigan, Department of Economics. \url{http://www-personal.umich.edu/~cattaneo/papers/RD-binselect.pdf}.

    Dinardo, J., and D. S. Lee. 2011. "Program Evaluation and Research Designs." In \emph{Handbook of Labor Economics}, ed. O. Ashenfelter and D. Card, vol. 4A, 463-536. Elsevier Science B.V.

    Imbens, G., and T. Lemieux. 2008. "Regression Discontinuity Designs: A Guide to Practice." \emph{Journal of Econometrics} 142(2): 615-635.

    Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in Economics." \emph{Journal of Economic Literature} 48(2): 281-355.
    
        
}
\author{
Sebastian Calonico, University of Michigan, Ann Arbor, MI.  \email{calonico@umich.edu}.

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
  \code{\link{rdbwselect}}, \code{\link{rdrobust}}
}