\name{rdrobust}
\alias{rdrobust}
\alias{print.rdrobust}
\alias{summary.rdrobust}
\title{Local-Polynomial RD Estimation with Robust Confidence Intervals}
\description{
    \code{\link{rdrobust}} implements a local-polynomial RD estimator with robust confidence intervals as proposed in Calonico, Cattaneo and Titiunik (2014a). It also computes alternative procedures available in the literature.

    Calonico, Cattaneo and Titiunik (2015b) provides an introduction to this command. Additional details for conventional approaches to conduct inference in the RD design can be found in Imbens and Lemieux (2008), Lee and Lemieux (2010), and references therein.
    
    For data-driven RD plots use the command \code{\link{rdplot}}, which implements the results in Calonico, Cattaneo and Titiunik (2015a).
    
    A companion \code{Stata} package is described in Calonico, Cattaneo and Titiunik (2014b).
}
\usage{rdrobust(y, x, subset = NULL, c = 0,  p = 1, q = 2, deriv = 0, 
            fuzzy = NULL, h = NULL, b = NULL, rho = NULL, scalepar = 1,
            kernel = "tri", bwselect = "CCT",  scaleregul = 1, delta = 0.5, 
            cvgrid_min = NULL, cvgrid_max = NULL, cvgrid_length = NULL, 
            cvplot = FALSE, vce = "nn", matches = 3, level = 95,
            all = FALSE)}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the running variable (a.k.a. score or forcing variable).}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{c}{specifies the RD cutoff in \code{x}; default is \code{c = 0}.}
  \item{p}{specifies the order of the local-polynomial used to construct the point-estimator; default is \code{p = 1} (local linear regression).}
  \item{q}{specifies the order of the local-polynomial used to construct the bias-correction; default is \code{q = 2} (local quadratic regression).}
    \item{deriv}{specifies the order of the derivative of the regression function to be estimated; default is \code{deriv = 0} (Sharp RD, or Fuzzy RD if \code{fuzzy} is also specified). Setting it equal to 1 results in estimation of a Kink RD design (or Fuzzy Kink RD if \code{fuzzy} is also specified).}
  \item{fuzzy}{specifies the treatment status variable used to implement Fuzzy RD estimation (or Fuzzy Kink RD if \code{deriv = 1} is also specified).
Default is Sharp RD design and hence this option is not used. 
\emph{Note}: For fuzzy RD designs, bandwidths are estimated using sharp RD bandwidth selectors for the reduced form outcome equation.}

  \item{h}{directly sets the main bandwidth. If not specified, it is computed by the companion command \code{\link{rdbwselect}}.}
  \item{b}{directly sets the pilot bandwidth. If not specified, it is computed by the companion command \code{\link{rdbwselect}}.}
  \item{rho}{sets the value of \code{rho}, so that the pilot bandwidth \code{b} equals \code{h/rho}. Default is \code{rho = 1} if \code{h} is specified but \code{b} is not.}
    \item{scalepar}{specifies scaling factor for RD parameter of interest. This option is useful when the population parameter of interest involves a known multiplicative factor (e.g., sharp kink RD).  Default is \code{scalepar = 1} (no scaling).}
    \item{kernel}{is the kernel function used to construct the local-polynomial estimator(s). Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.}
    
  \item{bwselect}{selects the bandwidth selection procedure to be used. By default it computes both \code{h} and \code{b}, unless \code{rho} is specified, in which case it only computes \code{h} and sets \code{b} = \code{h/rho}. Options are:
  
  \code{CCT} for bandwidth selector proposed by Calonico, Cattaneo and Titiunik (2014a)  (default option).
  
  \code{IK} for bandwidth selector proposed by Imbens and Kalyanaraman (2012) (only available for Sharp RD design).
  
  \code{CV} for cross-validation method proposded by Ludwig and Miller (2007) (only available for Sharp RD design).}
  
  \item{scaleregul}{specifies scaling factor for the regularization terms of \code{CCT} and \code{IK} bandwidth selectors. Setting \code{scaleregul = 0} removes the regularization term from the bandwidth selectors; default is \code{scaleregul = 1}.}

  \item{delta}{sets the quantile that defines the sample used in the cross-validation procedure. This option is used only  if \code{bwselect(CV)} is specified. See companion command \code{\link{rdbwselect}} for more details; default is \code{delta = 0.5}, that is, the median of the control and treated samples.}
  \item{cvgrid_min}{sets the minimum value of the bandwidth grid used in the cross-validation procedure. This option is used only if \code{bwselect = "CV"} is specified. See companion command \code{\link{rdbwselect}} for more details.}
  \item{cvgrid_max}{sets the maximum value of the bandwidth grid used in the cross-validation procedure. This option is used only if \code{bwselect = "CV"} is specified. See companion command \code{\link{rdbwselect}} for more details.}
  \item{cvgrid_length}{sets the bin length of the (evenly-spaced) bandwidth grid used in the cross-validation procedure. This option is used only if \code{bwselect = "CV"} is specified. See companion command \code{\link{rdbwselect}} for more details.}
  \item{cvplot}{generates a graph of the CV objective function. This option is used only if \code{bwselect = "CV"} is specified. See companion command \code{\link{rdbwselect}} for more details.}

\item{vce}{specifies the procedure used to compute the variance-covariance matrix estimator. This option is used only if \code{CCT} or \code{IK} bankdwith procedures are employed. Options are:

  \code{nn} for nearest-neighbor matches residuals using \code{matches} number of matches. This is the default option (with \code{matches = 3}, see below).

  \code{resid} for estimated plug-in residuals using \code{h} bandwidth.}
  \item{matches}{specifies the number of matches in the nearest-neighbor based variance-covariance matrix estimator. This options is used only when nearest-neighbor matches residuals are employed; default is \code{matches = 3}.}
    
  \item{level}{sets the confidence level for confidence intervals; default is \code{level = 95}.}
  \item{all}{if specified, \code{\link{rdrobust}} reports three different procedures:
  
  (i) conventional RD estimates with conventional standard errors.
  
  (ii) bias-corrected estimates with conventional standard errors.
  
  (iii) bias-corrected estimates with robust standard errors.}
}
\value{
  \item{c}{cutoff value.}
  \item{N_l}{sample size used to the left of the cutoff.}
  \item{N_r}{sample size used to the right of the cutoff.}
  \item{N}{overall sample size.}
  \item{p}{order of the polynomial used for estimation of the regression function.}
  \item{q}{order of the polynomial used for estimation of the bias of the regression function.}
  \item{h}{bandwidth used for estimation of the regression function.}
  \item{b}{bandwidth used for estimation of the bias of the regression function estimator.}
  \item{coef}{vector containing conventional and bias-corrected local-polynomial RD estimates.}
  \item{se}{vector containing conventional and robust standard errors of the local-polynomial RD estimates.}
  \item{pv}{vector containing the p-values associated with conventional, bias-corrected and robust local-polynomial RD estimates.}
  \item{ci}{matrix containing the confidence intervals associated with conventional, bias-corrected and robust local-polynomial RD estimates.}
}
\references{
        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014a. "Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs." \emph{Econometrica} 82(6): 2295-2326, November 2014.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}.

        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2014b. "Robust Data-Driven Inference in the Regression-Discontinuity Design." \emph{Stata Journal} 14(4): 909-946, 4th Quarter 2014.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_Stata.pdf}.

        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2015a. "Optimal Data-Driven Regression Discontinuity Plots." \emph{Journal of the American Statistical Association}, forthcoming.   \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2015_JASA.pdf}.

        Calonico, S., Cattaneo, M. D., and R. Titiunik. 2015b. "rdrobust: An R Package for Robust Inference in Regression-Discontinuity Designs." \emph{R Journal}, forthcoming. \url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2015_R.pdf}.
        
    Imbens, G., and T. Lemieux. 2008. "Regression Discontinuity Designs: A Guide to Practice." \emph{Journal of Econometrics} 142(2): 615-635.

    Imbens, G. W., and K. Kalyanaraman. 2012. "Optimal Bandwidth Choice for the Regression Discontinuity Estimator." \emph{Review of Economic Studies} 79(3): 933-959.
    
    Lee, D. S., and T. Lemieux. 2010. "Regression Discontinuity Designs in Economics." \emph{Journal of Economic Literature} 48(2): 281-355.

    Ludwig, J., and D. L. Miller. 2007. "Does Head Start Improve Children's Life Chances? Evidence from a Regression Discontinuity Design." \emph{Quarterly Journal of Economics} 122(1): 159-208.
       
}
\author{
Sebastian Calonico, University of Miami, Coral Gables, FL.  \email{scalonico@bus.miami.edu}.

Matias D. Cattaneo, University of Michigan, Ann Arbor, MI.  \email{cattaneo@umich.edu}.

Rocio Titiunik, University of Michigan, Ann Arbor, MI.  \email{titiunik@umich.edu}.
}
\examples{
x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdrobust(y,x)
}
\keyword{ RDD }
\keyword{ Robust Estimation }
\seealso{
  \code{\link{rdbwselect}}, \code{\link{rdplot}}
}
