#' A Constrained functional additive model (main function)
#'
#'
#' \code{cfam} is the main function for fitting the constrained functional additive model (CFAM).
#' \code{cfam} uses pre-treatment functional regressors, X1D, and scalar regressors, Xsc, for modeling a scalar-valued treatment outcome y.
#' For each of the functional regressors, \code{cfam} reduces its dimension via a data-driven linear projection to define a scalar-valued index variable,
#' and fits a treatment-specific additive model over those index variables and the scalar regressors simultaneousely.
#' The differential treatment effect is estimated by these treatment-specific nonparametrically-defined link functions on the indices and the scalar regressors.
#' For simultaneous variable selection for the treatment effect modifiers,
#'  \code{cfam} estimates a sparse combination of the components of the additive model via a \eqn{L1} regularization.
#' The resulting \code{cfam} object can be used to estimate an optimal treatment selection rule
#' for a new patient with pre-treatment clinical information.
#'
#' The sequence of the model coefficients implied by the tuning paramters \code{lambda} is fit by a coordinate descent algorithm.
#'
#'
#'
#' @param y   treatment outcomes, n-by-1 vector
#' @param Tr  treatment indicators, n-by-1 vector; each element represents one of the K available treatment options
#' @param X1D   a list of p pretreatment functional regressors; each of the functional regressors is evaluated on a grid (in [0,1]), giving a matrix (or a data.frame)
#' @param Xsc   a matrix of q pretreatment scalar-regressors, n-by-q matrix (or a n-by-q data.frame)
#' @param nbasis.t   a length K+1 vector; each element specifies the number of B-spline basis funtions for approximating the treatment-specific link function; the last element is for the "main effect" link function; the default is \code{nbasis.t=NULL}, and will be determined depending on the sample size.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters the p-splines in the link function estimation
#' @param eps   a value specifying the convergence criterion of algorithm for estimating coef (linear prjections associated with the functional regressors)
#' @param max.iter  an integer value specifying the maximum number of iterations for estimating coef (linear prjections associated with the functional regressors)
#' @param cd.eps  a value specifying the convergence criterion of coordinate descent algorithm for estimating the link functions
#' @param cd.max.iter  an integer value specifying the maximum number of coordinate descent iterations for estimating the link funtions
#' @param type  can choose bewteen \code{"AIC"}, \code{"BIC"}, and \code{"GCV"}, for the sparsity tuninig parameter selection.
#' @param lambda.grid   a grid of the sparsity tuning parameters, \code{lambda}, used in the \eqn{L1} regularization.
#' @param main.effect  if \code{TRUE}, first regress out the main effect of the regressors on the outcome, fitted by \code{fam}; the default is \code{FALSE}.
#' @param lambda.grid2  when \code{main.effect=TRUE}, a grid of the sparsity tuning parameters to be used for the main effect estimation
#' @param warm.start if \code{TRUE}, fit a \code{cmim} with some given penalty parameters to obtain a reasonable starting value for the index coeffiicent vectors.
#' @param warm.start.rho  a smoothness penalty to be used in obtaining an initial estimate.
#' @param warm.start.lam  a sparsity penalty to be used in obtaining an initial estimate.
#' @param csim.ini if \code{TRUE}, use a \code{csim} (a constrained single-index model) estimate as a warm start initial estimate  of the index coeffiicent vectors.
#' @param linear.link  a length (p+q) vector of 1's or 0's, indicating whether to restrict the jth link functions to be linear or not; if \code{NULL}, will take the vector of 0's.
#' @param trace if \code{TRUE}, show the trace of the fitting procedure; the default is \code{FALSE}.
#' @param plots if \code{TRUE}, produce plots of the estiamated link functions and the estiamated coefficient functions
#' @param scale.X  if \code{TRUE}, scale the FPC scores to have a unit variance for the parameter estimation (and will be scaled back in the final estimate)
#' @param sparse.coef  if \code{TRUE}, estimate spasre coef vectors associated with the FPC scores, using lasso with \code{glmnet}; the default is \code{FALSE}.
#' @param pve  the proportion of variance explained in determining the number of FPCs for the functional regressors
#'
#'
#' @return a list of information of the fitted model including
#' \item{cmim.obj}{a \code{cmim} (a constrained multiple-index model) object containing information about the fitted interaction-effect model, based on the FPC (scalar) scores and the scalar regressors.}
#' \item{mim.obj}{a \code{mim} (a multiple-index model) object containing information about the fitted main-effect model, based on the FPC (scalar) scores and the scalar regressors, if \code{main.effect=TRUE}.}
#' \item{scale.param}{a length (p+q) vector representing the scale parameters associated with the (p+q) regression component functions; 0 indicates the associated component function is shrinked to 0 (due to the L1 regularization).}
#' \item{scale.param.path}{a path of \code{scale.param} over a decreasing sequence of the sparsity parameters, \code{lambda}.}
#' \item{lambda.opt}{a value indicating the estimated optimal sparsity parameter, \code{lambda}.}
#' \item{coef}{a length p list of the estimated coef (the linear projection coefficients) associated with the FPC scores.}
#' \item{link.fn.fit}{a list of information about the estimated (p+q) link functions, including the knot sequences used in the B-spline approximations and the smoother matrices.}
#' \item{beta.t.coef}{a list of the estimated B-spline coefficient vectors associated with the (p+q) link functions.}
#' \item{beta.0.coef}{a list of the estimated B-spline coefficient vectors associated with the (p+q) link functions, used in orthogonalzing the link functions against the main-effect.}
#' \item{coef.fn}{a length p list of the estimated coefficient functions associated with the unctional regressors.}
#' \item{eigen.fn}{a length p list of the matrices of eigenfunctions associated with the p functional regressors.}
#' \item{mean.fn}{a length p list of the mean functions associated with the p functional regressors.}
#' \item{rho.opt}{a length (p+q) vector of the estimated optimal penalty parameters used in the P-splines, associated with the (p+q) link functions}
#' \item{link.fn.plot}{a plot object for the (p+q) link functions.}
#' \item{coef.fn.plot}{a plot object for the p coefficient functions.}
#'

#' @author Park, Petkova, Tarpey, Ogden
#' @import glmnet ggplot2 csim
#' @importFrom splines splineDesign
#' @importFrom plyr dlply
#' @importFrom refund fpca.face
#' @importFrom mgcv gam
#' @importFrom MASS lm.ridge
#' @importFrom magic adiag
#' @seealso \code{pred.cfam},  \code{cmim},  \code{fam}
#' @export
#' @examples
#' ## generate a dataset
#' dat <- dataGnFn(n=400, p=10, q=3, contrast ="nonlinear")
#' y <- dat$y   # a vector of scalar-valued outcomes
#' Tr <- dat$Tr  # a vector of treatment indicators
#' X1D <- dat$X  # a list of 1-D functional regressors
#' Xsc <- dat$Z  # a matrix of scalar regressors
#'
#' # a grid of tuning parameters for the treatment effect modifier selection
#' lambda.grid <- seq(0.1, 0.25, length.out = 10);
#' cfam.obj <- cfam(y=y, Tr=Tr, X1D=X1D, Xsc=Xsc,lambda.grid=lambda.grid, type="AIC",
#' nbasis.t= c(5,5,7), main.effect = FALSE)
#' cfam.obj$lambda.opt
#' #cfam.obj$coef  # coefficient associated with the functoinal covariates
#' # a length (p+q) vector of the shrinkage factors:
#' cfam.obj$scale.param  # 0 indicates that the corresponding is not selected as a treatment effect modifier
#' #cfam.obj$scale.param.path
#' cfam.obj$link.fn.plot[[1]]
#' cfam.obj$coef.fn.plot[[1]]
#'
#' # can also fit FAM to model the covariates' main effects;
#' # due to the orthogonality between the main and interaction effect components,
#' # the main effect model can be fitted seprately from the interaction effect model,
#' # and regressed out before fitting the interaction effect model.
#'
#' # a grid of tuning parameters for variable selection for the main effect model;
#' # this is to avoid over-fitting
#' lambda.grid2 <- seq(0.3, 0.6, length.out = 10);
#' cfam.obj2 <- cfam(y=y, Tr=Tr, X1D=X1D, Xsc=Xsc, lambda.grid=lambda.grid,
#' nbasis.t= c(5,5,7), main.effect = TRUE, lambda.grid2 =lambda.grid2)
#' cfam.obj2$lambda.opt
#' cfam.obj2$scale.param  # can identify the selected treatment effect modifiers
#' cfam.obj2$link.fn.plot[[1]]
#' cfam.obj2$coef.fn.plot[[1]]
#' cfam.obj2$mim.obj$scale.param  # can also identify the selected main-effect covariates
#' cfam.obj2$mim.obj$link.fn.plot[[1]]


## 1/13/2019
## constrained funcitonal additive models
#library(plyr)
#library(glmnet)
#library(magic)  # for adiag()
#library(ggplot2)
#library(mgcv)
#library(refund)
#library(MASS)
# X1D is a list of functional regressors
# Xsc is a matrix of scalar regressors
cfam <- function(y, Tr, X1D = NULL, Xsc = NULL,   # data
                 nbasis.t = NULL, # number of b-spline basis for the link functions
                 rho.grid = c(0.1, 0.3, 0.5), # a grid of penalty parameters for the p-splines in the link function estimation
                 eps = 0.05, # alpha estimation stopping rule
                 max.iter = 20, # max number of alpha estimation iterations
                 cd.eps = 0.005,  # coordinate descent stopping rule
                 cd.max.iter = 100,  # max number of coordinate descent iterations for fitting the link funtions
                 type = "AIC",  # "BIC","GCV"
                 lambda.grid = seq(0.05, 0.25, length.out =10), # lambda is the sparsity tuning parameters for variable selection
                 main.effect = FALSE, # if TRUE, then regress out the main effect estimated by fam().
                 lambda.grid2 = seq(0.1, 0.5, length.out =10),  # when main.effect = TRUE, the sparsity parameters to be used in the main effect estimation
                 warm.start = TRUE,  # if TRUE, fit a cmim() with some given penalty parameters to obtain a reasonable starting value for the index coeffiicents.
                 warm.start.rho = 0.1, # some reasonable smoothness penalty
                 warm.start.lam = 0.15,  # some reasonable sparsity penalty
                 csim.ini = TRUE,
                 linear.link = NULL, # a length (p+q) vector of 1's or 0's, indicating whether to restrict the jth link functions to be linear or not; if NULL, will take the vector of 0's.
                 trace = FALSE, plots = TRUE,
                 scale.X = TRUE,  # this indicates whether we scale the FPC scores to have unit variance or not (will be scaled back)
                 sparse.coef = FALSE,
                 pve = 0.95)  # the proportion of variance explained in determining the number of FPCs.
{

  if(!is.list(X1D))
    stop("X1D must be of class `list'");
  n <- length(y);
  p <- length(X1D);
  q <- ncol(Xsc);
  if(is.null(q))  q <- 0;
  K <- length(unique(Tr));

  # X.list is a list of data matrices associted with the p+q number of regressors
  X.list <- vector("list", length=p+q);

  # represent the p functional regressors in terms of (centered) FPC scores and the associetd eigenfunctions
  eigen.fn = mean.fn <- vector("list", length= p);
  if(p!=0)
  {
    for(j in 1:p)
    {
      tmp <- fpca.face(X1D[[j]], center = TRUE, pve=pve);
      eigen.fn[[j]] <- tmp$efunctions;
      mean.fn[[j]] <- tmp$mu;
      X.list[[j]] <- as.matrix(tmp$scores);
      rm(tmp);
    }
  }

  if(q!=0)
  {
    for(k in 1:q)
    {
      X.list[[p+k]] <- as.matrix(Xsc[,k]);
    }
  }

  if(is.null(linear.link))  linear.link <- rep(0, p+q);  # the default is nonparametrically-defined link.

  mim.obj <- NULL;
  if(main.effect)
  {
    mim.obj <- mim(y, X.list, nbasis=5, lambda.grid = lambda.grid2, rho.grid = rho.grid, scale.X = scale.X, sparse.coef = sparse.coef, plots=plots, type=type);
    y <- y - pred.mim(mim.obj);  # regress out the main effect
  }

  coef.warm.start <- NULL;
  if(warm.start)
  {
    cmim.warm.start <- cmim(y=y, Tr=Tr, X.list=X.list,
                            nbasis.t = c(rep(5, K), 7), rho.grid = warm.start.rho,
                            lambda.grid = warm.start.lam, type= type,
                            csim.ini = csim.ini,
                            eps = eps, max.iter=max.iter/2,
                            cd.eps = cd.eps, cd.max.iter = cd.max.iter,
                            linear.link = linear.link, trace = trace,
                            scale.X = scale.X, sparse.coef = FALSE, plots= FALSE)
    coef.warm.start <- cmim.warm.start$coef;
    scale.param.warm.start <- cmim.warm.start$scale.param;
    rm(cmim.warm.start);
  }

  # fit the constrained multiple index models (over a grid of the sparsity tuning parameters, lambda.grid)
  cmim.obj <- cmim(y=y, Tr=Tr, X.list=X.list,
                   nbasis.t = nbasis.t, rho.grid = rho.grid,
                   lambda.grid = lambda.grid, type= type,
                   coef = coef.warm.start,
                   csim.ini = csim.ini,
                   eps = eps, max.iter=max.iter,
                   cd.eps = cd.eps, cd.max.iter = cd.max.iter,
                   linear.link = linear.link, trace = trace,
                   scale.X = scale.X, sparse.coef = sparse.coef, plots= plots)
  link.fn.fit <- cmim.obj$link.fn.fit;
  cmim.obj$scale.param

  # obtain the coefficient functions associated with the functional regressors
  coef.fn = coef <- vector("list", length=p);
  if(p!=0)
  {
    for(j in 1:p)
    {
      if(is.null(cmim.obj$sd.X[[j]]))
      {
        coef[[j]] <- cmim.obj$coef[[j]];
      }else{
        tmp <- cmim.obj$coef[[j]]/ cmim.obj$sd.X[[j]];
        coef[[j]] <- tmp/sqrt(sum(tmp^2));
        rm(tmp);
      }
      coef.fn[[j]] <- eigen.fn[[j]] %*% coef[[j]];
    }
  }

  ###################
  ## visualization ##
  ###################
  # plots
  link.fn.plot = coef.fn.plot <- NULL
  if(plots)
  {
    # 1) link function plots
    link.fn.plot <- cmim.obj$link.fn.plot;
    # 2) coefficient function plot
    if(p!=0)
    {
      for(j in 1:p)
      {
        tmp <- seq(0, 1, length.out = nrow(coef.fn[[j]]));
        tmp2 <- data.frame(x= tmp, y= coef.fn[[j]]);
        coef.fn.plot[[j]] <- ggplot(data = tmp2, aes(x = x, y= y)) + geom_line(lwd=0.7, aes(colour = 2)) +
          scale_x_continuous(name = " ") + ylab(" ") + ggtitle(bquote(alpha[.(j)])) +theme(legend.position="none") + theme_bw() + theme(legend.position="none");
        rm(tmp); rm(tmp2);
      }
    }
  }

  results <-  list(
    cmim.obj = cmim.obj,
    mim.obj = mim.obj,
    scale.param = cmim.obj$scale.param,
    scale.param.path = cmim.obj$scale.param.path,
    coef.path = cmim.obj$coef.path,
    lambda.opt = cmim.obj$lambda.opt,
    coef = coef,
    link.fn.fit = link.fn.fit,
    beta.t.coef = cmim.obj$beta.t.coef,
    beta.0.coef = cmim.obj$beta.0.coef,
    intercept.y = cmim.obj$intercept.y,
    coef.fn = coef.fn,
    eigen.fn = eigen.fn,
    mean.fn = mean.fn,
    rho.opt = cmim.obj$rho.opt,
    link.fn.plot = link.fn.plot,
    coef.fn.plot = coef.fn.plot,
    X1D= X1D, Xsc=Xsc,
    p=p, q=q, n=n, K=K,
    coef.warm.start=coef.warm.start,
    scale.param.warm.start=scale.param.warm.start)

  class(results) <- c("cfam", "list")
  return(results)
}





#' A Constrained multiple-index model, for estimating interaction effects
#'
#' \code{cmim} is a function for fitting the constrained multiple-index model (CMIM).
#' Suppose there are V pre-specified groups of scalar-valued covariates for a treatment outcome variable, where each group consists of  multiple (or single) scalar-valued covariates.
#' \code{cmim}  fits a multiple-index model for estimating treatment-by-covariates interactions, where each index is defined to be a (data-driven) linear combination
#' of the covariates within group.
#' The function optimizes the within-group linear combinations, giving V scalar-valued index variables, and simultaneousely
#' fits a treatment-specific additive model over these scalar-valued index variables (i.e., a multiple-index model).
#' The differential treatment effect is estimated by these treatment-specific nonparametrically-defined link functions on the indices.
#' For simultaneous group selection for the treatment effect modifying groups,
#'  \code{cmim} estimates a sparse combination of the components of the additive model via a \eqn{L1} regularization.
#' The resulting \code{cmim} object can be used to estimate an optimal treatment selection rule
#' for a new patient with pre-treatment clinical information.
#'
#' The sequence of the model coefficients implied by the tuning paramters \code{lambda} is fit by a coordinate descent algorithm.
#'
#'
#'
#' @param y   treatment outcomes, n-by-1 vector
#' @param Tr  treatment indicators, n-by-1 vector; each element represents one of the K available treatment options
#' @param X.list  a list of V matrices, in which each matrix is the data matrix corresponding to each of the V groups of the scalar-valued covariates.
#' @param nbasis.t   a length K+1 vector; each element specifies the number of B-spline basis funtions for approximating the treatment-specific link function; the last element is for the "main effect" link function; the default is \code{nbasis.t=NULL}, and will be determined depending on the sample size.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters the P-splines in the link function estimation
#' @param eps   a value specifying the convergence criterion of algorithm for estimating the index coefficient vectors
#' @param max.iter  an integer value specifying the maximum number of iterations for estimating the index coefficient vectors
#' @param cd.eps  a value specifying the convergence criterion of coordinate descent algorithm for estimating the link functions
#' @param cd.max.iter  an integer value specifying the maximum number of coordinate descent iterations for estimating the link funtions
#' @param type  can choose bewteen \code{"AIC"}, \code{"BIC"}, and \code{"GCV"}, for the sparsity tuninig parameter selection.
#' @param lambda.grid   a grid of the sparsity tuning parameters, \code{lambda}, used in the \eqn{L1} regularization.
#' @param coef  a list of length p; can explicitely specify an initial estimate of the index coef vectors; the default is \code{NULL}.
#' @param csim.ini if \code{TRUE} and \code{coef} is not explictely given, use a \code{csim} (a constrained single-index model) estimate as an initial estimate of the index coefficient vectors.
#' @param linear.link  a length V vector of 1's or 0's, indicating whether to restrict the jth link functions to be linear or not; if \code{NULL}, will take the vector of 0's.
#' @param trace if \code{TRUE}, show the trace of the fitting procedure; the default is \code{FALSE}.
#' @param plots if \code{TRUE}, produce plots of the estiamated link functions and the estiamated coefficient functions
#' @param scale.X  if \code{TRUE}, scale the covariates to have a unit variance for the parameter estimation (and will be scaled back for the final estimate)
#' @param sparse.coef  if \code{TRUE}, estimate spasre index coefficient vectors, using lasso with \code{glmnet}; the default is \code{FALSE}.
#'
#'
#' @return a list of information of the fitted constrained multiple-index model including
#' \item{coef}{a length V list of the estimated index coefficient vectors}
#' \item{link.fn.fit}{a list of information about the estimated V link functions, including the knot sequences used in the B-spline approximations and the smoothers}
#' \item{scale.param}{a length V vector representing the scale parameters associated with the V regression component functions; 0 indicates the associated component function is shrinked to 0 (due to the L1 regularization).}
#' \item{scale.param.path}{a path of \code{scale.param} over a decreasing sequence of the sparsity parameters, \code{lambda}.}
#' \item{lambda.opt}{a value indicating the estimated optimal sparsity parameter, \code{lambda}.}
#' \item{beta.t.coef}{a length V list of the estimated B-spline coefficient vectors associated with the link functions.}
#' \item{beta.0.coef}{a length V list of the estimated B-spline coefficient vectors associated with the link functions, used in orthogonalzing the link functions against the main-effect.}
#' \item{rho.opt}{a length V vector of the estimated optimal penalty parameters used in the P-splines, associated with the V link functions}
#' \item{link.fn.plot}{a plot object for the (p+q) link functions.}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{pred.cmim},  \code{cfam},  \code{mim}, \code{fam}
#' @export

# X.list is a list of "group" covariates, each of them is a n x p_k matrices; p_k is the number of covariates for the kth group.
# Z.list is a list of "singleton" covariates, each of them a n x 1 vector
# this function optimizes the constrained additive model w.r.t. the interactions with the treatment, and performs trt effect modifier selection
cmim <- function(y, Tr, X.list,
                 nbasis.t = NULL,
                 rho.grid = c(0.1, 0.3, 0.5),
                 lambda.grid = seq(0.1, 0.3, length.out = 10),
                 type = "AIC", #"BIC"
                 coef = NULL,
                 csim.ini = FALSE,
                 cd.eps = 0.001, cd.max.iter=100,
                 linear.link = NULL,
                 eps=0.01, max.iter=20,
                 trace = FALSE, scale.X = TRUE,
                 sparse.coef=TRUE, plots = TRUE)
{

  ###################
  ## pre-procssing ##
  ###################
  # sort the observations by the order of the treatment indicators
  y <- data.preprocess(Tr=Tr, y=y)$y
  X.list <- data.preprocess(Tr=Tr, X1D= X.list)$X1D
  Tr <- data.preprocess(Tr=Tr)$Tr

  n <- length(y);
  K <- length(unique(Tr));
  V <- length(X.list)
  ind.groups <- which(sapply(X.list, ncol) > 1)

  # center and scale X
  sd.X = mean.X <- vector("list", length = V);
  if(scale.X)
  {
    for(j in 1:V)
    {
      tmp <- scale(X.list[[j]]);
      mean.X[[j]] <- attr(tmp, "scaled:center");
      sd.X[[j]] <- attr(tmp, "scaled:scale");
      X.list[[j]] <- tmp;
      rm(tmp);
    }
  }

  # center y within each treatment, and store the treatment-specific mean at intercept.y
  dat <- data.frame(Tr=Tr, y=y);
  dat.list <- dlply(dat, .(dat$Tr), function(x) {as.matrix(x[,-1])});
  intercept.y <- vector("list", K);
  yc <- NULL;
  for(t in 1:K)
  {
    dat.list[[t]] <- scale(dat.list[[t]], center = TRUE, scale = FALSE);
    intercept.y[[t]] <- attr(dat.list[[t]], "scaled:center");
    yc <- c(yc, dat.list[[t]])
  }
  rm(dat); rm(dat.list);

  if(is.null(nbasis.t))
  {
    nt <- summary(as.factor(Tr));
    for(t in 1:K)
    {
      nbasis.t[t] <- floor(nt[t]^{1/5.5})  + 4;
    }
    nbasis.t[K+1] <- floor(n^{1/5.5}) + 6;
  }

  if(is.null(linear.link)) linear.link <- rep(0, V);
  nonpar.link.thresh <- sum(nbasis.t[1:K])*2; # in practice, we require at least this many number of unique data points, to use a nonparametric link

  # initialize the index coefficients (coef)
  if(is.null(coef))
  {
    coef <- vector("list", V);
    for(j in 1:V)
    {
      if(ncol(X.list[[j]]) > 1)
      {
        if(csim.ini)
        {
          tmp <- drop(fit.csim(yc, Tr, X = X.list[[j]], nbasis.t = c(rep(5, K), 6), rho.grid =0.1, sparse = F, it.max= 2)$coef);
          if(sparse.coef) tmp[tmp < 0.01] <- 0;
          coef[[j]] <-  tmp/sqrt(sum(tmp^2));
          rm(tmp);
        }else{
          coef[[j]] <- c(1, rep(0, ncol(X.list[[j]])-1));
        }
      }else{
        coef[[j]] <- 1
      }
    }
  }

  # initialize the indices (u.list), and the smoothers (smoother)
  u.list = smoother  <- vector("list", V);
  for(j in 1:V)
  {
    # compute the jth index
    u.list[[j]] <-  drop(X.list[[j]] %*% coef[[j]]);
    # use a linear link, if not enough number of unique design points for nonparametric regrssion
    if(length(unique(u.list[[j]])) < nonpar.link.thresh)  linear.link[j] <- 1;
    # initialize the smoother matrices for fitting the link functions
    smoother[[j]] <- smoother.fn(Tr, u.list[[j]], nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link= linear.link[j]);
  }

  # over the path of lambdas, we will store the sets of fitted link functions (and the sets of scale parameters of the regressors);
  # will start from the most sparse fit (corresponding to the largest lambda),
  # to the most non-sparse fit (corresponding to the smallest lambda).
  lambda.grid <- sort(lambda.grid, decreasing =TRUE);
  link.fn.path = scale.param.path <- vector("list", length(lambda.grid));
  AIC.path = BIC.path = GCV.path <- vector("numeric", length(lambda.grid));

  # initialize the link funcion (yhat) and the scale parameters of the regressors (scale.param)
  yhat <- NULL;
  scale.param <- rep(1, V)

  # iterate until convergence (i.e., the projection coefficients change less than eps) or reaching max.iter
  sdiff <- 10^10;
  iter <- 0;
  diff <- sdiff;

  while(diff > eps)
  {
    iter <- iter + 1;
    coef.old <- coef;
    # optimize the constrained additive model via coordinate descent (for each lambda)
    for(lambda.index in seq_along(lambda.grid))
    {
      link.fn.fit  <- fit.link.fns(yc, smoother = smoother,
                                   yhat = yhat,
                                   lambda = lambda.grid[lambda.index],
                                   cd.max.iter = cd.max.iter, cd.eps = cd.eps,
                                   linear.link= linear.link,
                                   trace = trace)
      yhat <- link.fn.fit$yhat;    # this wil be used as a warm start for the subsequent iteration (and for the subsequent value of lambda)
      link.fn.path[[lambda.index]]  <- link.fn.fit;
      scale.param.path[[lambda.index]]  <- link.fn.fit$scale.param;
      AIC.path[lambda.index] <- link.fn.fit$AIC;
      BIC.path[lambda.index] <- link.fn.fit$BIC;
      GCV.path[lambda.index] <- link.fn.fit$GCV;
      #rm(link.fn.fit);
      #print(lambda.index)
    }

    # choose the "best" lambda
    if(type=="AIC")  lambda.index.opt <- which.min(AIC.path);
    if(type=="BIC")  lambda.index.opt <- which.min(BIC.path);
    if(type=="GCV")  lambda.index.opt <- which.min(GCV.path);

    # based on the "best" set of the link functions (that corresponds to the "best" lambda), optimize the projection coefficients
    coef <- fit.coefs(link.fn.fit = link.fn.path[[lambda.index.opt]],
                      X.list = X.list, coef = coef,
                      ind.groups = ind.groups, sparse.coef = sparse.coef)

    # based on the updated projection coefficients, update the indices and the smoother matrices
    for(j in 1:V)
    {
      if(ncol(X.list[[j]]) > 1)
      {
        if(scale.param.path[[lambda.index.opt]][j]!=0)
        {
          u.list[[j]] <- drop(X.list[[j]] %*% coef[[j]]);    # update the indices
          if(length(unique(u.list[[j]])) < nonpar.link.thresh)  linear.link[j] <- 1;
          smoother[[j]] <-  smoother.fn(Tr, u.list[[j]], nbasis.t = nbasis.t, rho.grid = rho.grid,
                                        linear.link = linear.link[j])   # update the smoother matrix
        }
      }
    }

    # keep track of the iterative procedure for estimating the projection coefficients
    ediff <- 0;
    for(j in 1:V)  ediff <- ediff + max(abs(coef[[j]] - coef.old[[j]]));
    diff <- abs(sdiff - ediff);
    sdiff <- ediff;
    #cat("iter", iter, "diff", diff, "\n");
    #cat("scale.param ", round(scale.param.path[[lambda.index.opt]], 2), "\n");
    if(iter >= max.iter | length(ind.groups) < 1)  break;
  }

  # once the algorithm converges, choose the final model
  if(type=="AIC") lambda.index.opt <- which.min(AIC.path);
  if(type=="BIC") lambda.index.opt <- which.min(BIC.path);
  if(type=="GCV") lambda.index.opt <- which.min(GCV.path);

  lambda.opt <- lambda.grid[lambda.index.opt];
  #lambda.opt
  link.fn.fit <- link.fn.path[[lambda.index.opt]];
  rm(link.fn.path);
  scale.param <- scale.param.path[[lambda.index.opt]];

  rho.opt <- numeric(V);
  for(j in 1:V)
  {
    rho.opt[j] <- link.fn.fit$rho.grid[link.fn.fit$rho.index.opt[j]];
  }


  # B-spline coefficients for the link functions
  beta.0.coef <- link.fn.fit$beta.0.coef;
  beta.t.coef <- vector("list", length=V);
  for(j in 1:V)
  {
    t.ind  <- unlist(lapply(1:K, function(x) rep(x, link.fn.fit$smoother[[j]]$nbasis.t[-(K+1)][x])));
    beta.t.coef[[j]] <- split(link.fn.fit$beta.t.coef[[j]], t.ind);
    rm(t.ind);
  }

  link.fn.plot <- NULL;
  if(plots)
  {
    # link function plots
    y.means <- rep(0, n);
    if(is.factor(Tr))
    {
      for(t in 1:K)  y.means <- y.means + intercept.y[[t]]*(Tr==attributes(Tr)$levels[t]);
    }else{
      for(t in 1:K)  y.means <- y.means + intercept.y[[t]]*(Tr==unique(Tr)[t]);
    }
    for(j in 1:V)
    {
      dat.j <- data.frame(y = y.means + link.fn.fit$resid.list[[j]]
                          - link.fn.fit$smoother[[j]]$B0 %*% link.fn.fit$beta.0.coef[[j]],
                          x = link.fn.fit$smoother[[j]]$u.t[[K+1]],
                          #Treatment = factor(Tr, labels=c("Placebo","Active drug")) )
                          Treatment = factor(Tr));

      link.fn.plot.j  <- ggplot(dat.j, aes(x = x, y = y, color=Treatment, shape=Treatment, linetype= Treatment))+
        geom_point(aes(color= Treatment, shape =Treatment), size= 1,  fill="white") +
        scale_colour_brewer(palette = "Set1", direction = -1) + theme( axis.title.x=element_text(size=15,face="bold")) +
        theme(title =element_text(size=12)) + ylab("(Adjusted) Partial residuals") + theme_bw(base_size = 14);

      if(linear.link[j])
      {
        link.fn.plot.j = link.fn.plot.j +
          geom_smooth(method=lm, se=TRUE, fullrange=FALSE, alpha = 0.35)
      }else{
        tmp1 <- 0;  for(t in 1:K)  tmp1 <- tmp1 + link.fn.fit$smoother[[j]]$nbasis.t[t];
        link.fn.plot.j <- link.fn.plot.j +
          geom_smooth(method=gam, formula = y ~ s(x, bs = "ps", k= floor(tmp1/K), sp= rho.opt[j]),
                      se=TRUE, fullrange=FALSE, alpha = 0.35)
      }

      if(j %in% ind.groups)  link.fn.plot.j <- link.fn.plot.j + xlab(bquote("<"~alpha[.(j)]~","~  x[.(j)]~">") );
      link.fn.plot[[j]] <- link.fn.plot.j;
      rm(link.fn.plot.j); rm(dat.j);
    }
  }

  results <- list(coef = coef,
                  link.fn.fit = link.fn.fit,
                  beta.t.coef = beta.t.coef,
                  beta.0.coef = beta.0.coef,
                  scale.param = scale.param,
                  intercept.y = intercept.y,
                  scale.param.path = scale.param.path,
                  lambda.grid = lambda.grid,
                  lambda.opt = lambda.opt,
                  AIC.path = AIC.path, BIC.path = BIC.path, GCV.path = GCV.path,
                  rho.grid = rho.grid,
                  rho.opt = rho.opt,
                  nbasis.t = nbasis.t,
                  y = y, Tr = Tr, X.list = X.list,
                  mean.X = mean.X, sd.X = sd.X, scale.X = scale.X,
                  link.fn.plot = link.fn.plot);
  class(results) <- c("cmim", "list");

  return(results)
}







# a set of optimal smoothing parameters, rho.opt, is chosen by minimizing GCV.
# a subfunction to fit a set of the (B-spline approximated) link functions using coordinate descent, given the list of the smoothers.
fit.link.fns <- function(yc, smoother = NULL,
                         yhat = NULL, scale.param = NULL, lambda = 1,
                         cd.max.iter = 100, cd.eps = 0.001,
                         linear.link= FALSE,
                         trace= FALSE, ortho.constr=TRUE)
{

  V <- length(smoother);  # this includes all scalar and functional predictors
  n <- length(yc);
  sqrt.n <- sqrt(n);

  # initialize the link functions
  if(is.null(yhat))
  {
    for(j in 1:V)
    {
      yhat[[j]] <- numeric(n);
    }
  }

  if(is.null(scale.param))  scale.param <- rep(1, V);

  # compute the initial residuals
  residuals <- yc;
  for(j in 1:V)
  {
    if(scale.param[j]!=0)  residuals <- residuals - yhat[[j]]
  }

  # perform a coordinate descent procedure to optimize the link functions, until convergence
  sdiff <- 100
  diff <- sdiff
  iter <- 0

  rho.grid <- smoother[[1]]$rho.grid;
  gcv  <- numeric(length= length(rho.grid));
  denom =  rho.index.opt <- rep(1, V);
  proj.Vt = proj.V0 = beta.0.coef = beta.t.coef <- vector("list", length= V);

  while(diff > cd.eps & diff <= sdiff & iter < cd.max.iter)
  {
    iter <- iter + 1;
    yhat.old <- yhat;

    for(j in 1:V)
    {
      # compute the jth partial residuals
      residuals <- residuals + yhat[[j]];

      # we will apply a ridge-type regularization (equivalent to an OLS with added 0s) to smooth the link functions (P-splines)
      residuals.aug <- c(residuals, rep(0, smoother[[j]]$ncol.Bt-2));

      # pick an optimal regularization (smoothing) parameter by generalized cross-validation (gcv)
      if(length(smoother[[j]]$rho.grid) >1)
      {
        gcv <- numeric();
        for(s in seq_along(smoother[[j]]$rho.grid))
        {
          gcv[s] <- sum((residuals - qr.fitted(smoother[[j]]$Bt.qr[[s]], residuals.aug)[1:n])^2) /(1 - smoother[[j]]$edf[[s]]/n)^2;
        }
        rho.index.opt[j] = s.opt <- which.min(gcv);
      }else{
        rho.index.opt[j] <- 1;
      }

      # compute the treatment-specific regression fit (without any group-wise penalization)
      proj.Vt[[j]]  <- qr.fitted(smoother[[j]]$Bt.qr[[rho.index.opt[j]]], residuals.aug);

      # compute the main effect
      if(ortho.constr)
      {
        beta.0.coef[[j]] <- qr.coef(smoother[[j]]$B0.qr, proj.Vt[[j]][1:n]);
        beta.0.coef[[j]][is.na(beta.0.coef[[j]])] <- 0;
      }else{
        beta.0.coef[[j]] <- rep(0, ncol(smoother[[j]]$B0));
      }
      # the fitted main effect
      proj.V0[[j]] <- drop(smoother[[j]]$B0 %*% beta.0.coef[[j]]);

      # regress out the main effect, proj.V0[[j]], from the fit, proj.Vt[[j]]:
      hold  <- proj.Vt[[j]] - c(proj.V0[[j]], rep(0, smoother[[j]]$ncol.Bt-2));

      # compute the sparsity scale parameter of the jth regressor
      denom[j] <- sqrt(sum(hold[1:n]^2));
      scale.param[j] <- (1 - lambda*sqrt.n/denom[j]);
      scale.param[j] <- as.numeric(scale.param[j] * (scale.param[j] > 0));

      beta.t.tmp <- qr.coef(smoother[[j]]$Bt.qr[[rho.index.opt[j]]], hold);
      beta.t.tmp[is.na(beta.t.tmp)] <- 0;
      beta.t.coef[[j]] <- scale.param[j]*beta.t.tmp;

      # compute the jth component regression function
      tmp  <- scale.param[j] * hold[1:n];
      yhat[[j]] <- tmp - mean(tmp);
      residuals <- residuals - yhat[[j]];

      rm(tmp); rm(hold); rm(beta.t.tmp);
    }

    diff <- 0;
    for(j in 1:V) diff <- diff + max(abs(yhat[[j]] - yhat.old[[j]]));

    if(diff > sdiff)
    {
      lambda <- lambda + 0.1;
      yhat <- yhat.old;
    }

    if(trace)
    {
      cat("MSE", var(residuals), "diff", diff, "\n")
      cat("L2 norm of components: ", denom, "\n")
    }
  }

  resid.list <- vector("list", V);
  EDF <- 0;
  working.residuals <- residuals;
  for(j in 1:V)
  {
    resid.list[[j]] <-  residuals + yhat[[j]]   # compute the partial-residuals for fitting the jth component.
    working.residuals <- working.residuals -  (scale.param[j]!=0) *proj.V0[[j]];
    EDF  <- EDF + (scale.param[j]!=0) *smoother[[j]]$edf[[rho.index.opt[j]]];   # effective degrees of freedom of the model
  }

  MSE <- mean(residuals^2);
  AIC <- log(MSE) + 2*EDF/n;
  BIC <- log(MSE) + log(n)*EDF/n;
  GCV <- MSE/ (1 - EDF/n)^2

  return(list(MSE=MSE, EDF=EDF,
              resid.list = resid.list,
              residuals = residuals,
              working.residuals=working.residuals,
              iter= iter, scale.param = scale.param,
              AIC= AIC, BIC=BIC, GCV=GCV,
              lambda = lambda,
              smoother = smoother,
              yhat = yhat, proj.Vt=proj.Vt, proj.V0=proj.V0,
              beta.t.coef = beta.t.coef, beta.0.coef = beta.0.coef,
              rho.grid = rho.grid, rho.index.opt = rho.index.opt,
              n=n, V=V)
  )
}


# a subfunction to compute the first derivative of the link functions, evaluated at the current single index variable u = alpha'X.
deriv_link.fns <- function(link.fn.fit, ind.groups = NULL)
{
  smoother <- link.fn.fit$smoother
  V <- length(smoother)
  K <- smoother[[1]]$K;
  d.link.fn <- vector("list", V);
  if(is.null(ind.groups))  ind.groups <- 1:V;

  # we only need to perform optimzation for the non-singleton groups; for the singleton group, the coef is simply 1.
  for(j in ind.groups)
  {
    if(link.fn.fit$scale.param[j]!=0)
    {
      u.t <- smoother[[j]]$u.t;
      knots.t <- smoother[[j]]$knots.t;
      t.ind  <- unlist(lapply(1:K, function(x) rep(x, smoother[[j]]$nbasis.t[-(K+1)][x])));
      beta.t.coef <- split(link.fn.fit$beta.t.coef[[j]], t.ind);

      d.design.t <- vector("list", K);
      for(t in 1:K)
      {
        if(smoother[[j]]$linear.link)
        {
          d.link.fn[[j]] <- c(d.link.fn[[j]], rep(beta.t.coef[[t]][2], length(u.t[[t]])) )
        }else{
          # compute the 1st derivative of the design functions
          d.design.t[[t]] <- splineDesign(knots.t[[t]], x=u.t[[t]], derivs=rep(1, length(u.t[[t]])), outer.ok = TRUE)
          d.link.fn[[j]] <- c(d.link.fn[[j]], d.design.t[[t]] %*% beta.t.coef[[t]])
        }
      }
      rm(d.design.t);
    }
  }

  return(list(d.link.fn=d.link.fn))
}



# for updating the projection coefficients, we only need to consider the groups with size > 1.
# for a singleton group, the projection coefficient is simply 1.
#ind.groups <- which(sapply(X.list, ncol) > 1)

# coef is the list of initial values for coef.
fit.coefs <- function(link.fn.fit, X.list, coef = NULL, ind.groups = NULL, sparse.coef= TRUE)
{

  V <- length(X.list);
  scale.param <- link.fn.fit$scale.param;
  #residuals <- link.fn.fit$residuals;
  working.residuals <- link.fn.fit$working.residuals;
  if(is.null(ind.groups))
  {
    ind.groups <- which(sapply(X.list, ncol) > 1);
  }
  if(is.null(coef))
  {
    for(j in 1:V)
    {
      coef[[j]] <- c(1, rep(0, ncol(X.list[[j]])-1));
    }
  }
  ## perform a linear approximation of the objective function w.r.t. coef, and minimize w.r.t. coef.
  # compute the 1st derivatives of the link functions
  d.link.fn.fit <- deriv_link.fns(link.fn.fit, ind.groups);
  # only need to perform optimzation for the non-singleton groups (i.e., ind.groups); for the singleton group, the coef is simply 1.
  for(j in ind.groups)
  {
    if(scale.param[j]!=0)
    {
      tmp <- d.link.fn.fit$d.link.fn[[j]] * X.list[[j]];
      y.star <- working.residuals + tmp %*%coef[[j]];
      if(sparse.coef)
      {
        glmnet.fit <- cv.glmnet(x= tmp, y= y.star, alpha=1, nfolds=5);
        lambda.index.opt <- max(which.min(glmnet.fit$cvm), min(which(glmnet.fit$nzero > 0)));
        lamda.opt <- glmnet.fit$lambda[lambda.index.opt];
        coef[[j]] <- round(coef(glmnet(x= tmp, y=y.star, lambda = lamda.opt, alpha=1)), 1)[-1];
      }else{
        coef[[j]] <- lm.ridge(y.star~ tmp-1, lambda = 0.01)$coef;
      }
      coef[[j]] <- coef[[j]] / sqrt(sum(coef[[j]]^2));
      names(coef[[j]]) <- colnames(X.list[[j]]);
      rm(tmp);
    }
  }
  return(coef)
}



#' \code{cfam} prediction function
#'
#' This function makes predictions from the constrained functional additive model, given a \code{cfam} object and pretreatment covariates, at which predictions are to be made.
#' The function returns predicted outcomes for each treatment and treatment selection rules.
#'
#' @param cfam.obj  a \code{cfam} object
#' @param X1D  a list of new values for the p pretreatment functional regressors (each of them evaluated on the same grid as the training data), at which predictions are to be made.
#' @param Xsc  a matrix of new values for the q pretreatment scalar regressors, at which predictions are to be made.
#'
#' @return
#' \item{pred.new}{a matrix of predicted values; column represents each treatment.}
#' \item{trt.rule}{a vector of suggested treatments}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{cfam}, \code{cmim}, \code{pred.cmim}
#' @export
#'
pred.cfam <- function(cfam.obj, Xsc=NULL, X1D =NULL)
{
  if(!inherits(cfam.obj, "cfam"))   # checks input
    stop("Object must be of class `cfam'");

  if(is.null(Xsc)) Xsc <- cfam.obj$Xsc;
  if(is.null(X1D)) X1D <- cfam.obj$X1D;
  if(!is.list(X1D)) stop("X1D must be of class `list'");

  p <- cfam.obj$p;
  q <- cfam.obj$q;
  K <- cfam.obj$K;

  X.list <- vector("list", length=p+q);
  if(p!=0)
  {
    n <- nrow(X1D[[1]]);
    for(j in 1:p)
    {
      m.fn  <- matrix(rep(cfam.obj$mean.fn[[j]], n), nrow = n, byrow = T);
      X.list[[j]]  <- drop(as.matrix(X1D[[j]] - m.fn) %*% cfam.obj$eigen.fn[[j]]);
    }
  }
  if(q!=0)
  {
    Xsc <- as.matrix(Xsc);
    for(k in c(1:q +p))
    {
      X.list[[k]]  <- Xsc[,k-p];
    }
  }

  results <- pred.cmim(cfam.obj$cmim.obj, X.list);
  if(!is.null(cfam.obj$mim.obj))  results$predicted  <- results$predicted + pred.mim(cfam.obj$mim.obj, X.list);

  return(results)
}


#' \code{cmim} prediction function
#'
#' This function makes predictions from the constrained multiple-index model, given a \code{cmim} object and pretreatment covariates, at which predictions are to be made.
#' The function returns predicted outcomes for each treatment and treatment selection rules.
#'
#' @param cmim.obj  a \code{cmim} object
#' @param X.list  a list of V matrices, in which each matrix corresponds to new values for each of the V groups of the covariates, at which predictions are to be made.
#'
#' @return
#' \item{pred.new}{a matrix of predicted values; column represents each treatment.}
#' \item{trt.rule}{a vector of suggested treatments}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{cmim}, \code{cfam}, \code{pred.cfam}
#' @export
#'
# prediction function for a constrained additive index model with multiple links
pred.cmim <- function(cmim.obj, X.list = NULL)
{
  if(!inherits(cmim.obj, "cmim"))   # checks input
    stop("Object must be of class `cmim'")

  if(is.null(X.list))  X.list <- cmim.obj$X.list;
  V <- length(X.list);
  n <- nrow(as.matrix(X.list[[1]]));
  K <- length(cmim.obj$intercept.y);

  index.score <- vector("list", length=V);
  for(j in 1:V)
  {
    if(!is.null(cmim.obj$sd.X[[j]]))
    {
      tmp <-  scale(X.list[[j]], cmim.obj$mean.X[[j]], cmim.obj$sd.X[[j]]);
    }else{
      tmp <-  scale(X.list[[j]], cmim.obj$mean.X[[j]], scale = FALSE);
    }
    index.score[[j]] <- as.matrix(tmp) %*% cmim.obj$coef[[j]];
    rm(tmp);
  }

  # compute treatment-specific predicted value
  predicted <- matrix(0, n, K);
  for(t in 1:K)
  {
    predicted.t <- cmim.obj$intercept.y[[t]]
    for(j in 1:V)
    {
      index.score[[j]][index.score[[j]] < cmim.obj$link.fn.fit$smoother[[j]]$u.min ] <- cmim.obj$link.fn.fit$smoother[[j]]$u.min
      index.score[[j]][index.score[[j]] > cmim.obj$link.fn.fit$smoother[[j]]$u.max ] <- cmim.obj$link.fn.fit$smoother[[j]]$u.max
      if(cmim.obj$link.fn.fit$smoother[[j]]$linear.link)
      {
        uj.expanded <- cbind(1, index.score[[j]])
      }else{
        uj.expanded <- splineDesign(cmim.obj$link.fn.fit$smoother[[j]]$knots.t[[t]], x= index.score[[j]], outer.ok = TRUE);
      }
      predicted.t  <- predicted.t + uj.expanded %*% cmim.obj$beta.t.coef[[j]][[t]]
      rm(uj.expanded);
    }
    predicted[, t] <- predicted.t
  }

  # compute treatment assignment
  trt.rule <- apply(predicted, 1, which.max)  # WLLG, we assume a higher y is prefered.
  if(K==2)  colnames(predicted) <- c("Tr1", "Tr2")

  results <- list(trt.rule = trt.rule,  pred.new =  predicted)
  return(results)
}





######################################################################
# a subfunction to construct (B-spline) smoother matrices, given the current index u = alpha'X.
# return: the QR decomposed design matrices (and the knot sequences used in constructing the B-spline design matrices);
######################################################################
smoother.fn <- function(Tr, u, nbasis.t=c(6,6,8), rho.grid = c(0, 0.25, 0.5), linear.link = FALSE)
{

  K <- length(unique(Tr));
  u.t = design.t <- vector("list", K+1);
  # create a list, dat.list, grouped by the treatment indicator
  dat <- data.frame(Tr=Tr, u =u);
  dat_list <- dlply(dat, .(dat$Tr), function(dat) as.matrix(dat[,-1]));
  for(t in 1:K)
  {
    u.t[[t]] <-  dat_list[[t]][,1];     # data points from the tth treatment group
  }
  u.t[[K+1]] <- u;
  u.min <- max(sapply(u.t, min));
  u.max <- min(sapply(u.t, max));
  # construct treatment group-specific design matrices
  design.t = knots.t  <- vector("list", K+1)
  for(t in 1:(K+1))
  {
    if(linear.link)   # if linear.link==TRUE, construct the linear model design matrix
    {
      nbasis.t[t] <- 2
      design.t[[t]] <- cbind(1, u.t[[t]])
    }else{
      #knots.t[[t]] <- c(rep(u.min, 3),  quantile(u.t[[t]], probs = seq(0, 1, length = nbasis.t[t] -2)),  rep(u.max,3))
      knots.t[[t]] <- seq(u.min, u.max, length.out= nbasis.t[t]+4);
      design.t[[t]] <-  splineDesign(knots.t[[t]], x= u.t[[t]], outer.ok = TRUE)
    }
  }
  # construct the block-diagonal matrix consist of the treatment-specific design matrices, to approximate E(Y| u, T)
  design.t.block <- NULL;
  for(t in 1:K)
  {
    design.t.block <- c( design.t.block, list(design.t[[t]]))
  }
  Bt <- Reduce(adiag,  design.t.block)  # Bt is the block-diagonal design matrix
  # QR decomposition of the design matrix Bt, given each value of the smoothness tuning parameter, rho
  Bt.qr <- vector("list", length(rho.grid))
  ncol.Bt <- ncol(Bt)
  D <- diff(diag(ncol.Bt), differences = 2)
  for(r in seq_along(rho.grid)) # a ridge-type smoothing (equivalent to a regular least squares with added observations)
  {
    #Bt.qr[[r]] <- qr(rbind(Bt, diag(sqrt(rho.grid[r]), ncol.Bt)))
    Bt.qr[[r]] <- qr(rbind(Bt, sqrt(rho.grid[r])*D))
  }
  # compute effective degrees of freedom of smoothers, so that later we use GCV to select an optimal smoothing parameter
  edf <- vector("list", length=length(rho.grid))
  svd.Bt <- svd(Bt)
  for(r in seq_along(rho.grid))
  {
    edf[[r]] <-  sum(svd.Bt$d[svd.Bt$d>0]^2/(svd.Bt$d[svd.Bt$d>0]^2 +rho.grid[r] )) #/K
  }

  # QR decomposition of the design matrix B0
  B0 <- design.t[[K+1]]
  B0.qr  <- qr(B0)
  results <- list(Bt.qr= Bt.qr, B0.qr= B0.qr, Bt=Bt, B0= B0,
                  u.t = u.t, u.min = u.min, u.max = u.max,
                  edf= edf, rho.grid = rho.grid, K=K,
                  knots.t = knots.t, nbasis.t = nbasis.t, ncol.Bt = ncol.Bt,
                  linear.link=linear.link)
  return(results)
}




########################################################
# A pre-processing  function that can be  used to create a data frame that consists of the observations
# in the order of the treatment indicators;
# i.e., (y, Tr, X) with Tr==1 come in the first rows of the data frame,
# (y, Tr, X) with Tr==2 second, ..., and (y, Tr, X) with Tr==K come last.
########################################################
data.preprocess <- function(y=NULL, Tr, Xsc=NULL, X1D=NULL, optTr=NULL)
{

  K <- length(unique(Tr))
  n <- length(y)
  p <- length(X1D)

  if(!is.null(Xsc))
  {
    q <- dim(Xsc)[2]
    dat.temp <- data.frame(Tr=Tr, X = Xsc)   # organize your data in a dataframe
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
    Xsc.new <- NULL
    for(t in 1:K) Xsc.new <- rbind(Xsc.new, dat.temp.list[[t]])
    Xsc <- Xsc.new
  } #else{
  #  q <- 0
  #}

  if(!is.null(X1D))
  {
    for(j in 1:p){
      X1D.j <- X1D[[j]]
      dat.temp <- data.frame(Tr=Tr, X = X1D.j)   # organize your data in a dataframe
      dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
      X1D.new <- NULL
      for(t in 1:K) X1D.new <- rbind(X1D.new, dat.temp.list[[t]])
      X1D[[j]] <- X1D.new
    }
  } #else{
  #  p <- 0
  #}

  if(!is.null(y))
  {
    dat.temp <- data.frame(Tr=Tr, y = y)   # organize your data in a dataframe
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x)  as.matrix(x[,-1]) )
    y.new <- NULL
    for(t in 1:K) y.new <- rbind(y.new, dat.temp.list[[t]])
    y <- y.new
  }

  if(!is.null(optTr))  # this is only avaialbe for a simulated dataset where we know the true model
  {
    optTr <- optTr
    dat.temp <- data.frame(Tr=Tr, optTr)
    dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
    optTr.new <- NULL
    for(t in 1:K) optTr.new <- rbind(optTr.new, dat.temp.list[[t]])
    optTr <- optTr.new  # this is in (-1, 1) values
    optTr <- 0.5*optTr +1.5  # this is in (1, 2) values
  }

  dat.temp <- data.frame(Tr=Tr, Tr )   # organize your data in a dataframe
  dat.temp.list <- dlply(dat.temp, .(dat.temp$Tr), function(x) {as.matrix(x[,-1])})
  Tr.new <- NULL
  for(t in 1:K) Tr.new <- rbind(Tr.new, dat.temp.list[[t]])
  if(is.factor(Tr))
  {
    Tr <- as.factor(Tr.new)
  }else{
    Tr <- as.vector(Tr.new)
  }
  #for(i in 1:n) if(Tr[i,]==-1) Tr[i,] <- 0

  return(list(y=y, Tr=Tr, Xsc=Xsc, X1D=X1D, optTr=optTr))
}



#####################################################################
# A performance mesure function;
# performance_measure assesses the performance of any given treatment selection rule,
# in terms of the "Value" of the treatment decision rule
# and the proportion of correct decisions (PCD).
#####################################################################
performance_measure <- function(pred.test, # predicted value for the testing data
                                y, Tr, # these are the testing data
                                value.opt= NULL, optTr =NULL)
{

  Tr.coded <- as.numeric(as.factor(Tr))

  if(is.null(value.opt))  value.opt <- 1

  #n <- length(y)
  ## Prediction error estimation
  #predicted <- numeric(n)
  #for(i in 1:n)   predicted[i] <- pred.test[i, Tr[i]]
  #MSE <-  mean((y - predicted)^2)
  #MAD <- mean(abs(y - predicted))
  #MSE.sd <- sd((y - predicted)^2) /sqrt(n)

  # Value estimation
  regime <- apply(pred.test, 1, which.max)
  right.trt.index <- regime == Tr.coded

  value <- sum(y[right.trt.index])/sum(right.trt.index)
  value.s <- value/value.opt

  #percent correct decision
  if(is.null(optTr)) optTr <-  numeric(length(y))
  pcd <- mean(regime == optTr)

  results <- list(#MSE=MSE, MAD=MAD,
    value =value, value.s= value.s, pcd = pcd)
  return(results)
}



#' A dataset simulation function
#'
#' \code{dataGenerationFn} generates an example dataset under a model that contains a main effect component, a treatment-by-covariates interaction effect component, and a random noise component.
#'
#' @param n  sample size.
#' @param p  number of functional-valued pretreatment covariates.
#' @param q  number of scalar-valued pretreatment covariates.
#' @param sigma.eps  standard deviation of the random noise term for the outcome model.
#' @param sigma.u  standard deviation of the measurement error for the functional covariates.
#' @param correlationZ  correlation among the scalar covariates.
#' @param sigmaZ   standard deviation of the scalar covariates.
#' @param contrast  indicates the shape of the treatment-specific link function that defines the interaction effect component.
#' \describe{
#' \item{"linear"}{a linear conntrast}
#' \item{"nonlinear"}{a nonlinear contrast; cosine functions}
#' }
#' @param delta  controls the intensity of the main effect.
#' \describe{
#' \item{\code{delta=1}}{moderate main effect}
#' \item{\code{delta=2}}{big main effect}
#' }
#' @param V  the number of support (evaluation) points for functional covariates
#' @param beta.fn  a list of index coefficient functions associated with the main effect; if \code{NULL}, the coefficients are randomly generated (and scaled to have a unit norm).
#' @param sim.set specifies the shape of the index coefficient functions associated with the interaction effect.
#'
#' @return
#' \item{y}{a n x 1 vector of treatment outcomes.}
#' \item{Tr}{a n x 1 vector of treatment indicators.}
#' \item{X}{a length p list of functional covariates; each element of the p list is a n x V matrix.}
#' \item{Z}{a n x q matrix of scalar covariates.}
#' \item{SNR}{the "signal" (interaction effects) to "nuisance" (main effects + noise) variance ratio (SNR) of the dataset.}
#' \item{optTr}{a n x 1 vector of treatments indicating the optimal treatment selections.}
#' \item{value.opt}{the "Value" of the optimal treatment selection rule, \code{optTr}.}
#' @export
#'
## a function for generating simulation data
dataGnFn <- function(n = 200, p = 5, q= 5,
                     sigma.eps= 0.4,  sigma.u= 0.4,
                     correlationZ= 0.1, sigmaZ= pi/2,
                     contrast ="linear",
                     V =100, # the number of support (evaluation) points for functional covariates
                     beta.fn =NULL, sim.set="A", delta=1)
{

  # define the interaction effect "link" functions
  if(contrast=="linear")
  {
    g.fn <- function(u)   0.5 * u
  }
  if(contrast=="nonlinear")
  {
    g.fn <- function(u)  cos(u)
  }

  # define the main effect "link" functions
  mu.fn <- function(u)  sin(u);

  # set up a set of functional basis
  s <- seq(0, 1, length.out =V)  # a grid of support points
  b1 <- sqrt(2)*sin(2*pi*s);
  b2 <- sqrt(2)*cos(2*pi*s);
  b3 <- sqrt(2)*sin(4*pi*s);
  b4 <- sqrt(2)*cos(4*pi*s);
  B <- cbind(b1, b2, b3, b4);   # a V by Nx matrix

  ## specify the alpha.fn:
  alpha.fn <- list(); Nx <- ncol(B);
  theta  <- matrix(rep(0, Nx*p), Nx, p);
  if(sim.set =="A") # simulation set A
  {
    theta[,1] <- c(1, 1, 1, 1)/sqrt(4)
    theta[,2] <- c(1, -1, 1, -1)/sqrt(4)
  }
  if(sim.set =="B") # simulation set B
  {
    theta[,1] <- c(1, 1, 0, 0)/sqrt(2)
    theta[,2] <-c(-1, -1, 0, 0)/sqrt(2)
  }
  for(j in 1:p)
  {
    alpha.fn[[j]] <- drop(B %*% theta[ ,j]);
    #plot(alpha.fn[[j]])
  }


  if(is.null(beta.fn))
  {
    # generate the main effect index functions (betas) (randomly generated coefficients)
    beta.fn <- list();
    eta  <- matrix(rnorm(Nx*p), Nx, p);
    eta <- apply(eta, 2, function(s)  s/sqrt(sum(s^2) ));
    for(j in 1:p)
    {
      beta.fn[[j]] <- drop(B %*% eta[ ,j]);
    }
  }

  # generate scalar covariates
  Psix <- sigmaZ*(diag(1 - correlationZ, nrow = q, ncol = q) + matrix(correlationZ, nrow = q, ncol = q));   # Z covariance matrix.
  ePsix <- eigen(Psix);
  Z <- sapply(1:q, function(x) rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors);

  # generate functional covarites
  X = xi <- vector("list", length= p);
  for(j in 1:p)
  {
    xi[[j]] <- matrix(rnorm(Nx*n, 0, sigmaZ), n, Nx);
    X[[j]] <- xi[[j]] %*% t(B) + matrix(rnorm(n*V, 0, sigma.u), n, V);
  }

  # generate the treatment indicators:
  Tr <- drop(rbinom(n, 1, 0.5) + 1)

  # generate the interaction effects:
  contrasts <- g.fn(Z[,1]) - g.fn(Z[,2]) + g.fn(drop(X[[1]]%*% alpha.fn[[1]])/V) - g.fn(drop(X[[2]]%*% alpha.fn[[2]])/V);
  interaction.effects <- contrasts*(-1)^(Tr);

  # generate the main effects:
  main.effects <- mu.fn(Z[,1]) + mu.fn(Z[,2]) + mu.fn(Z[,3]) +
    mu.fn(drop(X[[1]] %*% beta.fn[[1]])/V) + mu.fn(drop(X[[2]] %*% beta.fn[[2]])/V) + mu.fn(drop(X[[3]] %*% beta.fn[[3]])/V);

  # generate the noise:
  noise <- rnorm(n, 0, sigma.eps);

  # genreate the outcomes:
  Ey <- 1 + delta * main.effects  + interaction.effects;
  y <- Ey + noise;

  # some information about the datset
  var.interaction <- var(interaction.effects)
  var.main <- delta^2*var(main.effects)
  var.noise <- var(noise)
  SNR <- var.interaction/ (var.main + var.noise);
  SNR

  Y1 <- -contrasts;
  Y2 <-  contrasts;
  optTr <- as.numeric(Y2 > Y1) + 1;
  value.opt <- mean(Ey[Tr == optTr ])
  value.opt

  return(list(y=y, Tr=Tr, X = X, Z= Z, optTr= optTr, value.opt = value.opt, SNR=SNR))
}



######################################################################
## END OF THE FILE
######################################################################


