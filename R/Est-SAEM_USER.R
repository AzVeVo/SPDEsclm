#' ML estimation of spatial censored linear models via the SAEM algorithm using SPDE
#'
#' It fits the left, right, or interval spatial censored linear model using the
#' Stochastic Approximation EM (SAEM) algorithm by approximating the spatial process through
#' a Gaussian Markov random field obtained as a solution of a stochastic partial differential equation (SPDE).
#' It provides estimates and standard errors of the parameters and supports missing values on the dependent variable.
#'
#' @param y vector of responses of length \eqn{n}.
#' @param x design matrix of dimensions \eqn{n\times q}, where \eqn{q} is the number
#' of fixed effects, including the intercept.
#' @param ci vector of censoring indicators of length \eqn{n}. For each observation:
#' \code{1} if censored/missing, \code{0} otherwise.
#' @param lcl,ucl vectors of length \eqn{n} representing the lower and upper bounds
#' of the interval, which contains the true value of the censored observation. Default
#' \code{=NULL}, indicating no-censored data. For each observation: \code{lcl=-Inf} and
#' \code{ucl=c} (left censoring); \code{lcl=c} and \code{ucl=Inf} (right censoring); and
#' \code{lcl} and \code{ucl} must be finite for interval censoring. Moreover, missing data
#' could be defined by setting \code{lcl=-Inf} and \code{ucl=Inf}.
#' @param coords 2D spatial coordinates of dimensions \eqn{n\times 2}.
#' @param phi0 initial value for the spatial scaling parameter.
#' @param gamma0 initial value for the proportion of the total variance related to the spatial process.
#' @param lower,upper vectors of lower and upper bounds for the optimization method.
#' If unspecified, the default is \code{c(0.01,0.01)} for lower and \code{c(30,0.99)} for upper.
#' @param iedge the largest proportion allowed triangle edge. By default \code{iedge=0.30}.
#' @param MaxIter maximum number of iterations of the SAEM algorithm. By default \code{=300}.
#' @param M number of Monte Carlo samples for stochastic approximation. By default \code{=20}.
#' @param pc percentage of initial iterations of the SAEM algorithm with no memory.
#' It is recommended that \code{50<MaxIter*pc<100}. By default \code{=0.20}.
#' @param error maximum convergence error. By default \code{=1e-4}.
#' @param show_se logical. It indicates if the standard errors
#' should be estimated by default \code{=TRUE}.
#'
#' @details The spatial Gaussian model is given by
#'
#' \eqn{Y = X\beta + \xi},
#'
#' where \eqn{Y} is the \eqn{n\times 1} response vector, \eqn{X} is the \eqn{n\times q}
#' design matrix, \eqn{\beta} is the \eqn{q\times 1} vector of regression coefficients
#' to be estimated, and \eqn{\xi} is the error term which is normally distributed with
#' zero-mean and covariance matrix \eqn{\Sigma=\sigma^2[ \gamma R(\phi) + (1-\gamma) I_n]}.
#' The parameter \eqn{\sigma^2} denotes the total variance, \eqn{\gamma} denotes the proportion of the total variance
#' that is related to the spatial process, \eqn{\phi} represents the range parameter, and \eqn{R} is the correlation matrix,
#' which is computed from the Matérn correlation function with smooth parameter \eqn{\nu = 1}.
#' We assume that \eqn{\Sigma} is non-singular and \eqn{X} has full rank \insertCite{diggle2007springer}{SPDEcensSpatial}.
#'
#' The estimation process is performed via the SAEM algorithm, initially proposed by
#' \insertCite{delyon1999convergence;textual}{SPDEcensSpatial}. The spatial censored
#' (SAEM) algorithm was previously proposed by \insertCite{lachos2017influence;textual}{SPDEcensSpatial} and
#' \insertCite{ordonez2018geostatistical;textual}{SPDEcensSpatial} and is available in packages \code{CensSpatial}
#' and \code{RcppCensSpatial}.
#'
#' @note The SAEM final estimates correspond to the estimates obtained at the last iteration
#' of the algorithm.
#'
#' To fit a regression model for non-censored data, just set \code{ci} as a vector of zeros.
#'
#' @return An object of class "sclmSPDE". Generic functions \code{print} and \code{summary} have
#' methods to show the results of the fit. The function \code{plot} can extract
#' convergence graphs for the parameter estimates.
#'
#' Specifically, the following components are returned:
#' \item{Theta}{estimated parameters in all iterations, \eqn{\theta = (\beta, \sigma^2, \phi, \gamma)}.}
#' \item{theta}{final estimation of \eqn{\theta = (\beta, \sigma^2, \phi, \gamma)}.}
#' \item{beta}{estimated \eqn{\beta}.}
#' \item{sigma2}{estimated \eqn{\sigma^2}.}
#' \item{phi}{estimated \eqn{\phi}.}
#' \item{gamma}{estimated \eqn{\gamma}.}
#' \item{EY}{stochastic approximation of the first conditional moment.}
#' \item{EYY}{stochastic approximation of the second conditional moment.}
#' \item{SE}{vector of standard errors of \eqn{\theta = (\beta, \sigma^2, \phi, \gamma)}.}
#' \item{InfMat}{observed information matrix.}
#' \item{loglik}{log-likelihood for the SAEM method.}
#' \item{AIC}{Akaike information criterion.}
#' \item{BIC}{Bayesian information criterion.}
#' \item{Iter}{number of iterations needed to converge.}
#' \item{time}{processing time.}
#' \item{call}{\code{SPDEspatial} call that produced the object.}
#' \item{tab}{table of estimates.}
#' \item{critFin}{selection criteria.}
#' \item{range}{effective range.}
#' \item{ncens}{number of censored/missing observations.}
#' \item{MaxIter}{maximum number of iterations for the SAEM algorithm.}
#'
#' @author Pablo Zúñiga, Katherine L. Valeriano, and Lourival Lima
#'
#' @examples
#' library("RcppCensSpatial")
#' library("SPDEcensSpatial")
#'
#' # Example 1: 8% of left-censored observations
#' set.seed(1000)
#' n = 50   # Test with another values for n
#' coords = round(matrix(runif(2*n,0,15),n,2), 5)
#' x = cbind(rnorm(n), rnorm(n))
#' data = rCensSp(c(4,-2), 2, 3, 0.50, x, coords, "left", 0.08, 0,
#'                "matern", 1)
#'
#' fit = SPDEsclm(y=data$y, x=x, ci=data$ci, lcl=data$lcl, ucl=data$ucl,
#'                coords, phi0=2, lower=c(0.01,0.01), upper=c(30,0.99),
#'                iedge=0.30, MaxIter=100, M=10, pc=0.18)
#' fit
#' \donttest{
#' # Example 2: censored and missing observations
#' set.seed(123)
#' n = 200
#' coords = round(matrix(runif(2*n,0,20),n,2), 5)
#' x = cbind(runif(n), rnorm(n), rexp(n))
#' data = rCensSp(c(1,4,-1), 2, 3, 0.50, x, coords, "left", 0.05, 0,
#'                "matern", 1)
#' data$y[c(10,120)] = NA
#' data$ci[c(10,120)] = 1
#' data$ucl[c(10,120)] = Inf
#'
#' # Fitting the model
#' fit2 = SPDEsclm(data$y, data$x, data$ci, data$lcl, data$ucl, coords, phi0=2,
#'                 lower=c(0.01,0.01), upper=c(30,0.99), iedge=0.30, MaxIter=300,
#'                 M=20, pc=0.2, error=1e-4, show_se=TRUE)
#' fit2$tab
#' plot(fit2) }
#' @references \insertAllCited
#'
#' @export SPDEsclm
#'
#' @importFrom ggplot2 aes geom_line ggplot labs theme_bw
#' @importFrom gridExtra grid.arrange
#' @importFrom mvtnorm dmvnorm pmvnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom INLA inla.mesh.2d inla.mesh.fem inla.spde.make.A
#' @importFrom graphics points
#' @importFrom Rdpack reprompt
#' @importFrom RcppCensSpatial rCensSp


SPDEsclm = function(y, x, ci, lcl=NULL, ucl=NULL, coords, phi0, gamma0=0.50, lower=c(0.01,0.01), upper=c(30,0.99),
                    iedge=0.30, MaxIter=300, M=20, pc=0.2, error=1e-4, show_se=TRUE){
  n = length(c(y))
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x = as.matrix(x)
  if (!is.matrix(y)) y = as.matrix(y)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")

  #No data
  if ((length(x) == 0) | (length(y) == 0) | (length(ci) == 0)) stop("All parameters must be provided")

  #Validating if exists NA's
  if (sum(ci%in%c(0,1)) < length(ci)) stop("The elements of the vector ci must be 0 or 1")
  if (sum(is.na(x)) > 0) stop("There are some NA values in x")
  if (!all(c(is.finite(x)))) stop("x must contain only finite values.")
  if (sum(is.na(ci)) > 0) stop("There are some NA values in ci")
  miss = which(is.na(y))
  if (sum(ci[miss]) != length(miss)) stop ("NA values in y must be specified through arguments ci, lcl, and ucl")
  if (!all(c(is.finite(coords)))) stop("coords must contain only finite values.")

  #Validating dims data set
  if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if (ncol(as.matrix(ci)) > 1) stop("ci must have just one column")
  if (nrow(as.matrix(x)) != n) stop("x does not have the same number of lines than y")
  if (length(ci) != n) stop("ci does not have the same length than y")
  if (nrow(coords)!=n | ncol(coords)!=2) stop("Non-conformable dimensions between coords and y.")

  if (sum(ci) > 0){
    if (is.null(lcl) | is.null(ucl)) stop("lcl and ucl must be provided for censored data")
    if (!is.numeric(lcl) | !is.numeric(ucl)) stop("lcl and ucl must be numeric vectors")
    if (length(miss)>0){
      censor = (ci==1 & !is.na(y))
      if (any(is.infinite(lcl[censor])) & any(is.infinite(ucl[censor]))) stop("lcl or ucl must be finite for censored data")
    } else {
      if (any(is.infinite(lcl[ci==1])) & any(is.infinite(ucl[ci==1]))) stop("lcl or ucl must be finite for censored data")
    }
    if (length(lcl) != n) stop("lcl does not have the same length than y")
    if (length(ucl) != n) stop("ucl does not have the same length than y")
    if (ncol(as.matrix(lcl)) > 1) stop("lcl must have just one column")
    if (ncol(as.matrix(ucl)) > 1) stop("ucl must have just one column")
    if (sum(is.na(lcl))>0 | sum(is.na(ucl))>0) stop("There are some NA values in lcl or ucl")
    if (!all(lcl[ci==1]<ucl[ci==1])) stop ("lcl must be smaller than ucl")
  }

  #Validating supports
  if (length(c(phi0))>1 | !is.numeric(phi0)) stop("Initial value for phi must be provided")
  if (phi0 <= 0) stop("phi0 must be non-negative")
  if (length(c(gamma0))>1 | !is.numeric(gamma0)) stop("Initial value for gamma0 must be provided")
  if (gamma0<=0 | gamma0>=1) stop("gamma0 must be a number in (0, 1)")
  if (length(c(lower))!=2 | length(c(upper))!=2) stop("lower and upper must be vectors of length 2")
  if (any(is.na(lower)) | !is.numeric(lower)) stop("lower is not specified or contains NA")
  if (any(is.na(upper)) | !is.numeric(upper)) stop("upper is not specified or contains NA")
  if (any(lower>=upper)) stop("lower must be smaller than upper")
  if (any(lower<=0)) stop("Values in lower must be non-negative")
  if (!all(is.finite(upper))) stop("upper must contain only finite values")
  if (length(c(iedge))>1 | !is.numeric(iedge)) stop("Value for iedge must be a number in (0, 0.50)")
  if (iedge<=0 | iedge>=0.50) stop("iedge must be a number in (0, 0.50)")
  if (length(c(MaxIter))>1 | !is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
  if (MaxIter<=0 | MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  if (length(c(M))>1 | !is.numeric(M)) stop("M must be a positive integer value")
  if (M<=1 | M%%1!=0) stop("M must be a positive integer value (greater than 1)")
  if (length(c(pc))>1 | !is.numeric(pc)) stop("pc must be a real number in [0,1]")
  if (pc>1 | pc<0) stop("pc must be a real number in [0,1]")
  if (length(c(error))>1 | !is.numeric(error)) stop("error must be specified")
  if (error <= 0) stop("error must be a positive value (suggested to be small)")
  if (!is.logical(show_se)) stop("show_se must be logical (TRUE/FALSE).")

  ci  = as.matrix(ci)
  lcl = as.matrix(lcl)
  ucl = as.matrix(ucl)
  coords = as.matrix(coords)

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  out.Sp = SPDE_Spatial(y, x, ci, lcl, ucl, coords, phi0, gamma0, lower, upper,
                        iedge, MaxIter, M, pc, error, show_se)
  # Estimates
  q = ncol(x)
  lab = numeric(q + 3)
  if (sum(abs(x[,1])) == nrow(x)){ for (i in 1:q) lab[i] = paste('beta',i-1,sep='')
  } else { for (i in 1:q) lab[i] = paste('beta',i,sep='') }
  lab[q+1] = 'sigma2';  lab[q+2] = 'phi';  lab[q+3] = 'gamma'

  if (show_se) {
    tab = round(rbind(c(out.Sp$theta), c(out.Sp$SE)), 4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(c(out.Sp$theta)), 4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  # Information criteria
  critFin = c(out.Sp$loglik, out.Sp$AIC, out.Sp$BIC)
  critFin = round(t(as.matrix(critFin)), digits=3)
  dimnames(critFin) = list(c("Value"), c("Loglik", "AIC", "BIC"))

  out.Sp$call = match.call()
  out.Sp$tab = tab
  out.Sp$critFin = critFin
  out.Sp$ncens = sum(ci)
  out.Sp$MaxIter = MaxIter

  class(out.Sp) <- "sclmSPDE"
  return(out.Sp)
}

#' @export
summary.sclmSPDE = function(object, ...){
  cat('----------------------------------------------------------------\n')
  cat('     Censored Linear Spatial Regression Model using SPDE   \n')
  cat('----------------------------------------------------------------\n')
  cat("Call:\n")
  print(object$call)
  cat('\nEstimated parameters:\n')
  print(object$tab)
  cat('\n')
  cat('\nModel selection criteria:\n')
  print(object$critFin)
  cat('\nDetails:\n')
  cat('Number of censored/missing values:', object$ncens, '\n')
  cat('Convergence reached?:', (object$Iter < object$MaxIter), '\n')
  cat('Iterations:', object$Iter,'/',object$MaxIter, '\n')
  cat('Processing time:', round(object$time,4), units(object$time), '\n')
}


#' @export
print.sclmSPDE = function(x, ...){
  cat('----------------------------------------------------------------\n')
  cat('     Censored Linear Spatial Regression Model using SPDE   \n')
  cat('----------------------------------------------------------------\n')
  cat("Call:\n")
  print(x$call)
  cat('\nEstimated parameters:\n')
  print(x$tab)
  cat('\n')
  cat('\nModel selection criteria:\n')
  print(x$critFin)
  cat('\nDetails:\n')
  cat('Number of censored/missing values:', x$ncens,'\n')
  cat('Convergence reached?:',(x$Iter < x$MaxIter),'\n')
  cat('Iterations:',x$Iter,'/',x$MaxIter,'\n')
  cat('Processing time:',round(x$time,4),units(x$time),'\n')
}


#' @export
plot.sclmSPDE = function(x, ...){
  plotconvergence(x)
}


#' Prediction in spatial models with censored/missing responses
#'
#' It performs spatial prediction in a set of new \code{S} spatial locations.
#'
#' @param object object of class \code{'sclmSPDE'} given as output of \code{\link{SPDEsclm}} function.
#' @param coord.new matrix of coordinates for which prediction is performed.
#' @param x.new matrix of covariates for which prediction is performed.
#' @param ... further arguments passed to or from other methods.
#'
#' @details This function predicts using the mean squared error (MSE) criterion, which
#' takes the conditional expectation E(Y|X) as the best linear predictor.
#'
#' @return The function returns a list with:
#' \item{coord}{matrix of coordinates.}
#' \item{predValues}{predicted values.}
#' \item{sdPred}{predicted standard deviations.}
#'
#' @author Pablo Zúñiga, Katherine L. Valeriano, and Lourival Lima
#'
#' @seealso \code{\link{SPDEsclm}}
#' @examples \donttest{library("RcppCensSpatial")
#' library("SPDEcensSpatial")
#'
#' set.seed(1000)
#' n = 120
#' coords = round(matrix(runif(2*n,0,15),n,2), 5)
#' x = cbind(rbinom(n,1,0.50), rnorm(n), rnorm(n))
#' data = rCensSp(c(1,4,-1), 2, 3, 0.50, x, coords, "left", 0.10, 20, "matern", 1)
#'
#' ## Estimation
#' data1 = data$Data
#'
#' # Estimation: SAEM algorithm
#' fit = SPDEsclm(y=data1$y, x=data1$x, ci=data1$ci, lcl=data1$lcl,
#'                ucl=data1$ucl, coords=data1$coords, phi0=2.5,
#'                lower = c(0.01, 0.01), upper = c(30, 0.99), iedge = 0.3)
#' fit$tab
#'
#' # Prediction
#' data2 = data$TestData
#' pred = predict(fit, data2$x, data2$coords)
#'
#' # Cross-validation
#' mean((data2$y - pred$predValues)^2)}
#' @export
predict.sclmSPDE = function(object, x.new, coord.new, ...){

  if (is.null(object)) stop("object must be specified")

  if (!is.null(coord.new) & !is.null(x.new)){
    coord.new = as.matrix(coord.new)
    x.new     = as.matrix(x.new)
    if (!all(c(is.finite(coord.new)))) stop("coord.new must contain only finite values")
    if (!all(c(is.finite(x.new)))) stop("x.new must contain only finite values")
    if(ncol(x.new)!=length(c(object$beta))) stop("Non-conformable dimensions between x.new and beta")
    if (nrow(coord.new)!=nrow(x.new) | ncol(coord.new)!=2) stop("Non-conformable dimensions between coord.new and x.new")
  } else { stop("coord.new and x.new must be specified") }

  mats = object$mats
  cpred  = c(rep(0, nrow(object$EY)), rep(1, nrow(x.new)))
  A.pred = inla.spde.make.A(mesh = mats$mesh, loc = as.matrix(coord.new))
  output = prediccion(object, mats, A.pred, x.new, coord.new, cpred)
  return(output)
}
