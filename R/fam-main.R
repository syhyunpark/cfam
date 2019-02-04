#' A functional additive model, for estimating main effects
#'
#'
#' \code{fam} is the main function for fitting the functional additive model (FAM) of Fan et al. (2015).
#' \code{fam} uses functional regressors, X1D, and scalar regressors, Xsc, for modeling a scalar-valued outcome y.
#' For each of the functional regressors, \code{fam} reduces its dimension via a data-driven linear projection to define a scalar-valued index variable,
#' and fits an additive model over those index variables and the scalar regressors simultaneousely.
#' For simultaneous variable selection, \code{fam} estimates a sparse combination of the components of the additive model via a \eqn{L1} regularization.
#'
#' The sequence of the model coefficients implied by the tuning paramters \code{lambda} is fit by a coordinate descent algorithm.
#'
#'
#' @param y   treatment outcomes, n-by-1 vector
#' @param X1D   a list of p pretreatment functional regressors; each of the functional regressors is evaluated on a grid (in [0,1]), giving a matrix (or a data.frame)
#' @param Xsc   a matrix of q pretreatment scalar-regressors, n-by-q matrix (or a n-by-q data.frame)
#' @param nbasis   an interger value specifying the number of B-spline basis funtions for approximating the link functions; the default is \code{nbasis.t=NULL}, and will be determined depending on the sample size.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters the p-splines in the link function estimation
#' @param eps   a value specifying the convergence criterion of algorithm for estimating coef (linear prjections associated with the functional regressors)
#' @param max.iter  an integer value specifying the maximum number of iterations for estimating coef (linear prjections associated with the functional regressors)
#' @param cd.eps  a value specifying the convergence criterion of coordinate descent algorithm for estimating the link functions
#' @param cd.max.iter  an integer value specifying the maximum number of coordinate descent iterations for estimating the link funtions
#' @param type  can choose bewteen \code{"AIC"}, \code{"BIC"}, and \code{"GCV"}, for the sparsity tuninig parameter selection.
#' @param lambda.grid   a grid of the sparsity tuning parameters, \code{lambda}, used in the \eqn{L1} regularization.
#' @param warm.start if \code{TRUE}, fit a \code{mim} with some given penalty parameters to obtain a reasonable starting value for the index coeffiicent vectors.
#' @param warm.start.rho  a smoothness penalty to be used in obtaining an initial estimate.
#' @param warm.start.lam  a sparsity penalty to be used in obtaining an initial estimate.
#' @param sim.ini if \code{TRUE}, use a single-index model estimate as an initial estimate of the index coefficient vector.
#' @param linear.link  a length (p+q) vector of 1's or 0's, indicating whether to restrict the jth link functions to be linear or not; if \code{NULL}, will take the vector of 0's.
#' @param trace if \code{TRUE}, show the trace of the fitting procedure; the default is \code{FALSE}.
#' @param plots if \code{TRUE}, produce plots of the estiamated link functions and the estiamated coefficient functions
#' @param scale.X  if \code{TRUE}, scale the FPC scores to have a unit variance and perform estimation (and will be scaled back in the final estimate)
#' @param sparse.coef  if \code{TRUE}, estimate spasre coef vectors associated with the FPC scores, using lasso with \code{glmnet}; the default is \code{FALSE}.
#' @param pve  the proportion of variance explained in determining the number of FPCs for the functional regressors
#'
#'
#'
#' @return a list of information of the fitted model including
#' \item{mim.obj}{a \code{mim} (a multiple-index model) object containing information about the fitted main-effect model, based on the FPC (scalar) scores and the scalar regressors.}
#' \item{scale.param}{a length (p+q) vector representing the scale parameters associated with the (p+q) regression component functions; 0 indicates the associated component function is shrinked to 0 (due to the L1 regularization).}
#' \item{scale.param.path}{a path of \code{scale.param} over a decreasing sequence of the sparsity parameters, \code{lambda}.}
#' \item{lambda.opt}{a value indicating the estimated optimal sparsity parameter, \code{lambda}.}
#' \item{coef}{a length p list of the estimated coef (the linear projection coefficients) associated with the FPC scores.}
#' \item{link.fn.fit}{a list of information about the estimated (p+q) link functions, including the knot sequences used in the B-spline approximations.}
#' \item{beta.0.coef}{a list of the estimated B-spline coefficient vectors associated with the (p+q) link functions.}
#' \item{smoother}{an object containing information about the estimated link functions, including the knot sequences used in the B-spline approximation \code{knots.t.}}
#' \item{coef.fn}{a length p list of the estimated coefficient functions associated with the unctional regressors.}
#' \item{eigen.fn}{a length p list of the matrices of eigenfunctions associated with the p functional regressors.}
#' \item{mean.fn}{a length p list of the mean functions associated with the p functional regressors.}
#' \item{rho.opt}{a length (p+q) vector of the estimated optimal penalty parameters used in the P-splines, associated with the (p+q) link functions}
#' \item{link.fn.plot}{a plot object for the (p+q) link functions.}
#' \item{coef.fn.plot}{a plot object for the p coefficient functions.}
#'

#' @author Park, Petkova, Tarpey, Ogden
#' @import refund plyr glmnet ggplot2
#' @seealso \code{pred.fam},  \code{mim},  \code{cfam}
#' @export
#' @examples
#' ## generate a dataset
#' dat <- dataGnFn(n=300, p=10, q=3, contrast ="nonlinear")
#' y <- dat$y   # a vector of scalar-valued outcomes
#' X1D <- dat$X  # a list of 1-D functional regressors
#' Xsc <- dat$Z  # a matrix of scalar regressors
#'
#' lambda.grid <- seq(0.1, 0.25, length.out = 10);  # a grid of tuning parameters for variable selection
#' fam.obj <- fam(y=y, X1D=X1D, Xsc=Xsc, lambda.grid=lambda.grid, type="AIC", sim.ini = TRUE,  sparse.coef=FALSE)
#' fam.obj$lambda.opt
#' #fam.obj$coef  # coefficient associated with the functoinal covariates
#' fam.obj$scale.param   # the shrinkage factors:  0 indicates that the corresponding variale is unselected
#' #fam.obj$scale.param.path
#' fam.obj$link.fn.plot[[1]]
#' fam.obj$coef.fn.plot[[1]]


## 1/13/2019
## funcitonal additive models
## the following code is for fitting functional additive index model to estimate the main effect.
#library(plyr)
#require(fda)
#library(glmnet)
#library(ggplot2)
#library(mgcv)
#library(refund)
# Xsc is a matrix of scalar covariates
# X1D is a list of functional covariates
fam <- function(y, X1D=NULL, Xsc=NULL,  # data
                nbasis = NULL, # number of b-spline basis for the link functions
                rho.grid = c(0.1, 0.2, 0.3), # grid of the ridge-type smoothing parameter for the link functions
                eps= 0.05, # alpha estimation stopping rule
                max.iter = 20, # max number of alpha estimation iterations a
                cd.eps = 0.005,  # coordinate descent stopping rule
                cd.max.iter = 100,  # max number of  coordinate descent iterations
                lambda.grid = seq(0.05, 0.25, length.out =10), # lambda is the tuning parameters for the sparsity penalty
                type = "AIC",  # "BIC", "GCV"
                warm.start = TRUE,  # if TRUE, fit a mim() with some given penalty parameters to obtain a reasonable starting value for the index coeffiicents.
                warm.start.rho = 0.1, # some reasonable smoothness penalty
                warm.start.lam = 0.15,  # some reasonable sparsity penalty
                sim.ini = TRUE,   # if sim.ini = TRUE, initial solution for the index coefficients is determined by fitting a single-index model
                linear.link = NULL,  # whether we restrict the jth component's link functions to be linear; a vector of 1's or 0's; if NULL, take a vector of 0's.
                trace = FALSE, plots = TRUE,
                scale.X = TRUE,  # whether we scale the FPC scores X.scores in performing estimation (will be scaled back)
                sparse.coef = TRUE,
                pve=0.99)
{

  if(!is.list(X1D))
    stop("X1D must be of class `list'");
  n <- length(y);
  p <- length(X1D);
  q <- ncol(Xsc);
  if(is.null(q))  q <- 0;
  # X.list is a list of data matrices associted with the p+q number of covariates
  X.list <- vector("list", length=p+q);

  # for the p funcitonal covariates, express them in terms of centered FPC scores and the associetd eigenfunctions
  eigen.fn = mean.fn <- vector("list", length= p);
  if(p!=0)
  {
    for(j in 1:p)
    {
      fpca.obj.j <- fpca.face(X1D[[j]], center = TRUE, pve=pve)
      eigen.fn[[j]] <- fpca.obj.j$efunctions
      mean.fn[[j]] <- fpca.obj.j$mu
      X.list[[j]] <- as.matrix(fpca.obj.j$scores);
    }
  }

  if(q!=0)
  {
    for(k in 1:q)
    {
      X.list[[p+k]] <- as.matrix(Xsc[,k]);
    }
  }

  if(is.null(linear.link))  linear.link <- rep(0, p+q);  # so, the default is nonparametrically-defined link.

  coef.warm.start <- NULL;
  if(warm.start)
  {
    mim.warm.start <- mim(y=y, X.list=X.list,
                          nbasis=5,
                          rho.grid = warm.start.rho,
                          lambda.grid = warm.start.lam,
                          type= type,
                          sim.ini = sim.ini,
                          eps = eps, max.iter=max.iter,
                          cd.eps = cd.eps, cd.max.iter = cd.max.iter,
                          linear.link = linear.link, trace = trace,
                          scale.X = scale.X, sparse.coef = FALSE, plots= FALSE);
    coef.warm.start <- mim.warm.start$coef;
    scale.param.warm.start <- mim.warm.start$scale.param;
    rm(mim.warm.start);
  }


  # fit the multiple-index model (over a grid of the sparsity tuning parameters, lambda.grid)
  mim.obj <- mim(y=y, X.list=X.list,
                 nbasis = nbasis, rho.grid = rho.grid,
                 lambda.grid = lambda.grid, type= type,
                 coef = coef.warm.start,
                 sim.ini = sim.ini,
                 eps = eps, max.iter=max.iter,
                 cd.eps = cd.eps, cd.max.iter = cd.max.iter,
                 linear.link = linear.link, trace = trace,
                 scale.X = scale.X, sparse.coef = sparse.coef, plots= plots);

  link.fn.fit <- mim.obj$link.fn.fit;

  # obtain the coefficient functions associated with the functional covariates
  coef.fn = coef <- vector("list", length=p);
  if(p!=0)
  {
    for(j in 1:p)
    {
      if(is.null(mim.obj$sd.X[[j]]))
      {
        coef[[j]] <- mim.obj$coef[[j]];
      }else{
        tmp <- mim.obj$coef[[j]]/ mim.obj$sd.X[[j]];
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
    link.fn.plot <- mim.obj$link.fn.plot
    # 2) coefficient function plot
    if(p!=0)
    {
      for(j in 1:p)
      {
        xj <- seq(0, 1, length.out = nrow(coef.fn[[j]]));
        dat2.j <- data.frame(x= xj, y= coef.fn[[j]]);
        coef.fn.plot.j <- ggplot(data = dat2.j, aes(x = x, y= y)) + geom_line(lwd=0.7, aes(colour = 2)) +
          scale_x_continuous(name = " " ) + ylab(" ") + ggtitle(bquote(alpha[.(j)])) +theme(legend.position="none") + theme_bw() + theme(legend.position="none");
        coef.fn.plot[[j]] <- coef.fn.plot.j;
        rm(coef.fn.plot.j)
      }
    }
  }

  results <- list(
    mim.obj= mim.obj,
    scale.param = mim.obj$scale.param,
    scale.param.path = mim.obj$scale.param.path,
    coef.path = mim.obj$coef.path,
    lambda.opt = mim.obj$lambda.opt,
    coef = coef,
    link.fn.fit = link.fn.fit,
    beta.0.coef = mim.obj$beta.0.coef,
    intercept.y = mim.obj$intercept.y,
    coef.fn = coef.fn,
    eigen.fn = eigen.fn,
    mean.fn = mean.fn,
    rho.grid = mim.obj$rho.grid,
    rho.opt = mim.obj$rho.opt,
    link.fn.plot = link.fn.plot,
    coef.fn.plot = coef.fn.plot,
    X1D= X1D, Xsc=Xsc,
    p=p, q=q, n=n,
    coef.warm.start=coef.warm.start,
    scale.param.warm.start=scale.param.warm.start)

  class(results) <- c("fam", "list")
  return(results)
}





#' A multiple-index model, for estimating main effects
#'
#' \code{mim} is a function for fitting the multiple-index model (MIM).
#' Suppose there are V pre-specified groups of scalar-valued regressors, in which each group consists of single or multiple scalar-valued regressor(s).
#' \code{mim} fits a multiple-index model, where each of the indices is defined as a (data-driven) unknown linear combination of the regressors within each group.
#' \code{mim} optimizes these within-group linear combinations, resulting in V scalar-valued index variables,
#' and simultaneousely
#' fits an additive model over these scalar-valued indices (i.e., a multiple-index model).
#' For simultaneous group selection,
#'  \code{mim} estimates a sparse combination of the components of this additive model via a \eqn{L1} regularization.
#'
#' The sequence of the model coefficients implied by the tuning paramters \code{lambda} is fit by a coordinate descent algorithm.
#'
#'
#'
#' @param y   treatment outcomes, n-by-1 vector
#' @param X.list  a list of V matrices, in which each matrix is the data matrix corresponding to each of the V groups of the scalar-valued regressors.
#' @param nbasis   an integer value spefifying the number of B-spline basis funtions for approximating the link function; the default is \code{nbasis=NULL}, and will be determined depending on the sample size.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters the P-splines in the link function estimation
#' @param eps   a value specifying the convergence criterion of algorithm for estimating the index coefficient vectors
#' @param max.iter  an integer value specifying the maximum number of iterations for estimating the index coefficient vectors
#' @param cd.eps  a value specifying the convergence criterion of coordinate descent algorithm for estimating the link functions
#' @param cd.max.iter  an integer value specifying the maximum number of coordinate descent iterations for estimating the link funtions
#' @param type  can choose bewteen \code{"AIC"}, \code{"BIC"}, and \code{"GCV"}, for the sparsity tuninig parameter selection.
#' @param lambda.grid   a grid of the sparsity tuning parameters, \code{lambda}, used in the \eqn{L1} regularization.
#' @param coef  a list of length p; can explicitely specify an initial estimate of the index coef vectors; the default is \code{NULL}.
#' @param sim.ini if \code{TRUE}, use a sigle-index model estimate as an initial estimate for the index coefficient vectors.
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
# fit a multiple index model
mim  <- function(y, X.list,
                 nbasis = NULL,
                 rho.grid = c(0, 0.25, 0.5),
                 lambda.grid = seq(0.1, 1, length.out = 30),
                 type = "AIC", #"BIC"
                 coef = NULL,
                 sim.ini = TRUE,
                 cd.eps = 0.001, cd.max.iter=100,
                 linear.link = NULL,
                 eps=0.01, max.iter=20,
                 trace = FALSE, scale.X = TRUE, sparse.coef=TRUE, plots = TRUE)
{

  n <- length(y);
  V <- length(X.list);
  ind.groups <- which(sapply(X.list, ncol) > 1)

  # center y
  intercept.y <- mean(y);
  yc <- y - intercept.y;
  # center and scale X
  sd.X = mean.X <- vector("list", length = V);
  if(scale.X)
  {
    for(j in 1:V)
    {
      tmp <- scale(X.list[[j]]);
      mean.X[[j]] <- attr(tmp, "scaled:center");
      sd.X[[j]] <- attr(tmp, "scaled:scale");
      X.list[[j]] <- tmp; rm(tmp);
    }
  }

  if(is.null(nbasis))
  {
    nbasis <- floor(n^{1/5.5}) + 4
  }

  if(is.null(linear.link)) linear.link <- rep(0, V);
  nonpar.link.thresh <- nbasis*2;

  # initialize the index coefficients (coef)
  if(is.null(coef))
  {
    coef <- vector("list", V);
    for(j in 1:V)
    {
      if(ncol(X.list[[j]]) > 1)
      {
        if(sim.ini)
        {
          tmp <- ppr(x=X.list[[j]], y = yc, nterms=1, sm.method= "gcvspline");
          coef[[j]] <- tmp$alpha; rm(tmp);
          #lasso.temp <- cv.glmnet(x= X.list[[j]], y= yc, nfolds = 5)
          #lambda.index.opt <- max(which.min(lasso.temp$cvm), min(which(lasso.temp$nzero > 0)))
          #lamda.opt <- lasso.temp$lambda[lambda.index.opt]
          #lasso.fit <- glmnet(x= X.list[[j]], y= yc, lambda = lamda.opt, alpha = 1)
          #alpha.coef <- as.numeric(lasso.fit$beta)
          #coef[[j]] <- alpha.coef/sqrt(sum(alpha.coef^2))
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
    smoother[[j]] <- smoother.fn2(u.list[[j]], nbasis = nbasis, rho.grid = rho.grid, linear.link= linear.link[j]);
  }


  # over the path of lambdas, will store the sets of estimated link functions (and the associated scale parameters)
  lambda.grid <- sort(lambda.grid, decreasing =TRUE); # start from the most sparse fit (the largest lambda) to the most non-sparse fit (the smallest lambda)
  link.fn.path = scale.param.path <- vector("list", length(lambda.grid));
  AIC.path = BIC.path = GCV.path <- vector("numeric", length(lambda.grid));

  # initialize the link funcion (yhat) and the scale parameters (scale.param)
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

    # start from the most non-sparse fit (the smallest lambda) to the most sparse fit (the largest lambda)
    # optimize the constrained additive model via coordinate descent (for each lambda)
    for(lambda.index in seq_along(lambda.grid))
    {
      link.fn.fit  <- fit.link.fns2(yc, smoother = smoother,
                                    yhat = yhat,
                                    lambda = lambda.grid[lambda.index],
                                    cd.max.iter = cd.max.iter, cd.eps = cd.eps,
                                    linear.link= linear.link,
                                    trace = trace)

      yhat <- link.fn.fit$yhat    # this wil be used as a warm start for the subsequent iteration
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
    coef <- fit.coefs2(link.fn.fit = link.fn.path[[lambda.index.opt]],
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
          smoother[[j]] <-  smoother.fn2(u.list[[j]], nbasis = nbasis, rho.grid = rho.grid, linear.link = linear.link[j])   # update the smoother matrix
        }
      }
    }

    # keep track of the iterative procedure for estimating the projection coefficients
    ediff <- 0;
    for(j in 1:V)  ediff <- ediff + max(abs(coef[[j]] - coef.old[[j]]));
    diff <- abs(sdiff - ediff);
    sdiff <- ediff;
    #cat("iter", iter, "diff", diff, "\n");
    #cat("scale.param ", round(scale.param.path[[lambda.index.opt]], 2), "\n")
    if(iter >= max.iter | length(ind.groups) < 1)  break;
  }

  # once the algorithm is converged, choose the final model
  if(type=="AIC")  lambda.index.opt <- which.min(AIC.path);
  if(type=="BIC")  lambda.index.opt <- which.min(BIC.path);
  if(type=="GCV")  lambda.index.opt <- which.min(GCV.path);

  lambda.opt <- lambda.grid[lambda.index.opt];
  link.fn.fit <- link.fn.path[[lambda.index.opt]];
  rm(link.fn.path);
  scale.param <- scale.param.path[[lambda.index.opt]];

  rho.opt <- numeric(V);
  for(j in 1:V)
  {
    rho.opt[j] <- link.fn.fit$rho.grid[link.fn.fit$rho.index.opt[j]];
  }

  # b-spline coefficients for the link functions
  beta.0.coef <- link.fn.fit$beta.0.coef;

  link.fn.plot <- NULL;
  if(plots)
  {
    # link function plots
    for(j in 1:V)
    {
      dat.j <- data.frame(y = link.fn.fit$resid.list[[j]], x = link.fn.fit$smoother[[j]]$u);

      link.fn.plot.j  <- ggplot(dat.j, aes(x = x, y = y)) + geom_point(size= 1,  fill="white") +
        scale_colour_brewer(palette = "Set1", direction = -1) + theme( axis.title.x=element_text(size=15,face="bold")) +
        theme(title =element_text(size=12)) + ylab("(Adjusted) Partial residuals") + theme_bw(base_size = 14);

      if(linear.link[j])
      {
        link.fn.plot.j = link.fn.plot.j + geom_smooth(method=lm, se=TRUE, fullrange=FALSE, alpha = 0.35)
      }else{
        link.fn.plot.j <- link.fn.plot.j + geom_smooth(method=gam, se=TRUE, fullrange=FALSE, alpha = 0.35,
                                                       formula = y ~ s(x, bs = "ps", k= link.fn.fit$smoother[[j]]$nbasis,
                                                                       sp= rho.opt[j]));
      }

      if(j %in% ind.groups)  link.fn.plot.j <- link.fn.plot.j + xlab(bquote("<"~alpha[.(j)]~","~  x[.(j)]~">") );
      link.fn.plot[[j]] <- link.fn.plot.j;
      rm(link.fn.plot.j);
      rm(dat.j);
    }
  }

  results <- list(coef = coef,
                  link.fn.fit = link.fn.fit,
                  beta.0.coef= beta.0.coef,
                  scale.param = scale.param,
                  intercept.y = intercept.y,
                  #link.fn.path = link.fn.path,
                  scale.param.path = scale.param.path,
                  lambda.grid = lambda.grid,
                  lambda.opt = lambda.opt,
                  AIC.path = AIC.path, BIC.path= BIC.path, GCV.path= GCV.path,
                  rho.opt = rho.opt,
                  nbasis = nbasis,
                  y=y, X.list=X.list, intercept.y = intercept.y,
                  mean.X=mean.X, sd.X=sd.X, scale.X=scale.X,
                  link.fn.plot=link.fn.plot)
  class(results) <- c("mim", "list")

  return(results)
}



######################################################################
# a subfunction to construct (B-spline) smoother matrices, given the current index u = alpha'X.
# return: the QR decomposed design matrices (and the knot sequences used in constructing the B-spline design matrices);
######################################################################
smoother.fn2 <- function(u, nbasis = 6, rho.grid = c(0, 0.25, 0.5), linear.link = FALSE)
{

  u.min <- min(u);
  u.max <- max(u);
  # construct a design matrix B0
  if(linear.link)   # if linear.link==TRUE, construct the linear model design matrix
  {
    nbasis <- 2
    B0  <- cbind(1, u)
  }else{
    #knots <- c(rep(u.min, 3),  quantile(u, probs = seq(0, 1, length = nbasis -2)),  rep(u.max,3))
    knots <- seq(u.min, u.max, length.out = nbasis+4)
    B0 <-  splineDesign(knots, x= u, outer.ok = TRUE)
  }

  # QR decomposition of the design matrix B0, given each value of the smoothness tuning parameter, rho
  B0.qr <- vector("list", length(rho.grid))
  ncol.B0 <- ncol(B0)
  D <- diff(diag(ncol.B0), differences = 2)
  for(r in seq_along(rho.grid)) # a ridge-type smoothing (equivalent to a regular least squares with added observations)
  {
    #B0.qr[[r]] <- qr( rbind(B0, diag(sqrt(rho.grid[r]), ncol.B0) ) )
    B0.qr[[r]] <- qr( rbind(B0, sqrt(rho.grid[r])*D) )
  }


  # compute effective degrees of freedom of smoothers, so that later we use GCV to select an optimal smoothing parameter
  edf <- vector("list", length=length(rho.grid))
  svd.B0 <- svd(B0)
  for(r in seq_along(rho.grid))
  {
    edf[[r]] <-  sum(svd.B0$d[svd.B0$d>0]^2/(svd.B0$d[svd.B0$d>0]^2 +rho.grid[[r]] )) #/K
  }

  results <- list(B0.qr= B0.qr, B0= B0, ncol.B0 =ncol.B0,
                  u = u, u.min = u.min, u.max = u.max,
                  edf= edf, rho.grid = rho.grid,
                  knots = knots, nbasis = nbasis,
                  linear.link=linear.link)
  return(results)
}




# a set of optimal smoothing parameters, rho.opt, is chosen by minimizing GCV.
# a subfunction to fit a set of the (B-spline approximated) link functions using coordinate descent, given the list of the smoothers.
fit.link.fns2 <- function(yc, smoother = NULL,
                          yhat = NULL,
                          scale.param = NULL, lambda = 1,
                          cd.max.iter = 100, cd.eps = 0.001,
                          linear.link= FALSE,
                          trace= FALSE, ortho.constr=TRUE)
{

  V <- length(smoother);  # this includes all scalar and functional predictors
  n <- length(yc);
  sqrt.n <- sqrt(n);
  # initialize the link functions
  if(is.null(yhat)) for(j in 1:V) yhat[[j]] <- numeric(n);
  if(is.null(scale.param))  scale.param <- rep(1,V);

  # compute the initial residuals
  residuals <- yc;
  for(j in 1:V) if(scale.param[j]!=0)  residuals <- residuals - yhat[[j]]

  # perform a coordinate descent procedure to optimize the link functions, until convergence
  sdiff <- 100
  diff <- sdiff
  iter <- 0

  rho.grid <- smoother[[1]]$rho.grid;
  gcv  <- numeric(length= length(rho.grid));
  denom =  rho.index.opt <- rep(1, V);
  proj.V0 = beta.0.coef <- vector("list", length= V);

  while(diff > cd.eps & diff <= sdiff & iter < cd.max.iter)
  {
    iter <- iter + 1
    yhat.old <- yhat

    for(j in 1:V)
    {
      # compute the jth partial residuals
      residuals <- residuals + yhat[[j]];

      # apply a ridge-type regularization (equivalent to an OLS with added 0s) for estimating the interaction effect
      residuals.aug <- c(residuals, rep(0, smoother[[j]]$ncol.B0-2));

      # pick an optimal regularization (smoothing) parameter by generalized cross-validation (gcv)
      if(length(smoother[[j]]$rho.grid) >1)
      {
        gcv <- numeric()
        for(s in seq_along(smoother[[j]]$rho.grid))
        {
          gcv[s] <- sum((residuals - qr.fitted(smoother[[j]]$B0.qr[[s]], residuals.aug)[1:n])^2) /(1 - smoother[[j]]$edf[[s]]/n)^2
        }
        rho.index.opt[j] = s.opt <- which.min(gcv)
      }else{
        rho.index.opt[j] <- 1
      }

      # compute the main effect
      beta.0.tmp <- qr.coef(smoother[[j]]$B0.qr[[rho.index.opt[j]]], residuals.aug);
      beta.0.tmp[is.na(beta.0.tmp)] <- 0;
      proj.V0[[j]] <- drop(smoother[[j]]$B0 %*% beta.0.tmp);

      # compute the sparsity scale parameter
      denom[j] <- sqrt(sum(proj.V0[[j]]^2));
      scale.param[j] <- (1 - lambda*sqrt.n/denom[j]);
      scale.param[j] <- as.numeric(scale.param[j] * (scale.param[j] > 0));

      beta.0.coef[[j]] <- scale.param[j] *beta.0.tmp; rm(beta.0.tmp);

      # compute the jth regression function
      yhat.tmp <- scale.param[j] * proj.V0[[j]];
      yhat[[j]] <- yhat.tmp - mean(yhat.tmp); rm(yhat.tmp);
      residuals <- residuals - yhat[[j]];
    }

    diff <- 0
    for(j in 1:V)   diff <- diff + max(abs(yhat[[j]] - yhat.old[[j]]));

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
  for(j in 1:V)
  {
    resid.list[[j]] <-  residuals + yhat[[j]]   # compute the partial-residuals for fitting the jth component.
    EDF  <- EDF + (scale.param[j]!=0) *smoother[[j]]$edf[[rho.index.opt[j]]];   # effective degrees of freedom of the model
  }

  MSE <- mean(residuals^2);
  AIC <- log(MSE) + 2*EDF/n;
  BIC <- log(MSE) + log(n)*EDF/n;
  GCV <- MSE/ (1 - EDF/n)^2

  return(list(MSE=MSE, EDF=EDF,
              resid.list = resid.list,
              residuals = residuals,
              iter= iter, scale.param = scale.param,
              AIC= AIC, BIC=BIC, GCV=GCV,
              lambda = lambda,
              smoother = smoother,
              yhat = yhat,
              proj.V0=proj.V0,
              beta.0.coef = beta.0.coef,
              rho.grid = rho.grid, rho.index.opt = rho.index.opt,
              n=n, V=V)
  )
}



# for updating the projection coefficients, we only need to consider the groups with size > 1.
# for a singleton group, the projection coefficient is simply 1.
#ind.groups <- which(sapply(X.list, ncol) > 1)

# coef is the list of initial values for coef.
fit.coefs2 <- function(link.fn.fit, X.list, coef = NULL, ind.groups = NULL, sparse.coef= TRUE)
{

  V <- length(X.list);
  scale.param <- link.fn.fit$scale.param;
  residuals <- link.fn.fit$residuals;
  if(is.null(ind.groups))  ind.groups <- which(sapply(X.list, ncol) > 1);
  if(is.null(coef))
  {
    for(j in 1:V)
    {
      coef[[j]] <- c(1, rep(0, ncol(X.list[[j]])-1));
    }
  }

  # perform a linear approximation of the objective function w.r.t alpha.j's, and minimize w.r.t. alpha.j's
  # compute the first derivatives of the link functions
  d.link.fn.fit <- deriv_link.fns2(link.fn.fit, ind.groups);
  # we only need to perform optimzation for the non-singleton groups; for the singleton group, the coef is simply 1.
  for(j in ind.groups)
  {
    if(scale.param[j])
    {
      tmp <- as.vector(d.link.fn.fit$d.link.fn[[j]]) * X.list[[j]];
      y.star <- residuals + tmp %*%coef[[j]];
      if(sparse.coef)
      {
        glmnet.fit <- cv.glmnet(x= tmp, y= y.star, alpha=1, nfolds=5);
        lambda.index.opt <- max(which.min(glmnet.fit$cvm), min(which(glmnet.fit$nzero > 0)));
        lamda.opt <- glmnet.fit$lambda[lambda.index.opt];
        coef[[j]] <- round(coef(glmnet(x= tmp, y=y.star, lambda = lamda.opt, alpha=1)), 1)[-1];
      }else{
        coef[[j]] <- lm.ridge(y.star ~ tmp-1, lambda = 0.01)$coef;
      }
      coef[[j]] <- coef[[j]] / sqrt(sum(coef[[j]]^2));
      names(coef[[j]]) <- colnames(X.list[[j]]);
      rm(tmp);
    }
  }

  return(coef)
}


deriv_link.fns2 <- function(link.fn.fit, ind.groups = NULL)
{
  smoother <- link.fn.fit$smoother;
  beta.0.coef <- link.fn.fit$beta.0.coef;
  V <- length(smoother);
  d.link.fn <- vector("list", V);
  if(is.null(ind.groups))  ind.groups <- 1:V;
  # we only need to perform optimzation for the non-singleton groups; for the singleton group, the coef is simply 1.
  for(j in ind.groups)
  {
    if(link.fn.fit$scale.param[j])
    {
      u <- smoother[[j]]$u
      knots <- smoother[[j]]$knots
      if(smoother[[j]]$linear.link)
      {
        d.link.fn[[j]] <-  rep(beta.0.coef[[j]][2], length(u))
      }else{
        d.link.fn[[j]] <- splineDesign(knots, x=u, derivs=rep(1,length(u)), outer.ok = TRUE) %*% beta.0.coef[[j]];
      }
    }
  }
  return(list(d.link.fn=d.link.fn))
}



#' \code{fam} prediction function
#'
#' This function makes predictions from the functional additive model, given a \code{fam} object and covariates, at which predictions are to be made.
#' The function returns predicted outcomes.
#'
#' @param fam.obj  a \code{fam} object
#' @param X1D  a list of new values for the p functional regressors (each of them evaluated on the same grid as the training data), at which predictions are to be made.
#' @param Xsc  a matrix of new values for the q scalar regressors, at which predictions are to be made.
#'
#' @return
#' \item{predicted}{a vector of predicted values.}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{fam}, \code{mim}, \code{pred.mim}
#' @export
#'
pred.fam <- function(fam.obj, Xsc=NULL, X1D =NULL)
{
  if(!inherits(fam.obj, "fam"))   # checks input
    stop("Object must be of class `fam'");

  if(is.null(Xsc)) Xsc <- fam.obj$Xsc;
  if(is.null(X1D)) X1D <- fam.obj$X1D;
  if(!is.list(X1D)) stop("X1D must be of class `list'");

  p <- fam.obj$p;
  q <- fam.obj$q;

  X.list <- vector("list", length=p+q);
  if(p!=0)
  {
    n <- nrow(X1D[[1]]);
    for(j in 1:p)
    {
      m.fn  <- matrix(rep(fam.obj$mean.fn[[j]], n), nrow = n, byrow = T);
      X.list[[j]]  <- drop(as.matrix(X1D[[j]] - m.fn) %*% fam.obj$eigen.fn[[j]]);
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

  results <- pred.mim(mim.obj = fam.obj$mim.obj, X.list);
  return(results)
}



#' \code{mim} prediction function
#'
#' This function makes predictions from the multiple-index model, given a \code{mim} object and pretreatment covariates.
#' The function returns predicted outcomes.
#'
#' @param mim.obj  a \code{mim} object
#' @param X.list  a list of V matrices, in which each matrix corresponds to new values for each of the V groups of the covariates, at which predictions are to be made.
#'
#' @return
#' \item{predicted}{a vector of predicted values}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{mim}, \code{fam}, \code{pred.fam}
#' @export
#'
# prediction function for a constrained additive index model with multiple links
pred.mim <- function(mim.obj, X.list = NULL) #Xsc.new=NULL, X1D.new=NULL)
{
  if(!inherits(mim.obj, "mim"))   # checks input
    stop("Object must be of class `mim'")

  if(is.null(X.list))  X.list <- mim.obj$X.list;

  V <- length(X.list);
  n <- nrow(as.matrix(X.list[[1]]));
  K <- length(mim.obj$intercept.y);

  index.score <- vector("list", length=V);
  for(j in 1:V)
  {
    if(!is.null(mim.obj$sd.X[[j]]))
    {
      tmp <-  scale(X.list[[j]], mim.obj$mean.X[[j]], mim.obj$sd.X[[j]]);
    }else{
      tmp <-  scale(X.list[[j]], mim.obj$mean.X[[j]], scale = FALSE);
    }
    index.score[[j]] <- as.matrix(tmp) %*% mim.obj$coef[[j]];
    rm(tmp);
  }

  # compute the predicted values
  predicted <- numeric(n);
  predicted <- predicted + mim.obj$intercept.y;
  for(j in 1:V)
  {
    if(mim.obj$scale.param[j])
    {
      index.score[[j]][index.score[[j]] < mim.obj$link.fn.fit$smoother[[j]]$u.min ] <- mim.obj$link.fn.fit$smoother[[j]]$u.min
      index.score[[j]][index.score[[j]] > mim.obj$link.fn.fit$smoother[[j]]$u.max ] <- mim.obj$link.fn.fit$smoother[[j]]$u.max
      if(mim.obj$link.fn.fit$smoother[[j]]$linear.link)
      {
        uj.expanded <- cbind(1, index.score[[j]])
      }else{
        uj.expanded <- splineDesign(mim.obj$link.fn.fit$smoother[[j]]$knots, x= index.score[[j]], outer.ok = TRUE);
      }
      predicted  <- predicted + uj.expanded %*% mim.obj$beta.0.coef[[j]]
      rm(uj.expanded);
    }
  }

  return(drop(predicted))
}
