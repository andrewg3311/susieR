#' @title SUm of Single Effect (SuSiE) Logistic Regression of Y on X
#' @details Performs Bayesian logistic regression of Y on X.
#' That is, this function
#' fits the regression model Y ~ Bern(sum_l Xb_l), where
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0, var = prior_variance).
#' @param X an n by p matrix of covariates
#' @param Y an n vector (binary)
#' @param L maximum number of non-zero effects
#' @param prior_variance the prior variance (vector of length L, or scalar. In latter case gets repeated L times). The prior variance on each non-zero element of b is set to be prior_variance.
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting.
#' Note that `prior_variance` specifies the prior on the coefficients of X *after* standardization (if performed).
#' If you do not standardize you may need
#' to think more carefully about specifying
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE). The latter is generally not recommended.
#' @param Z an n by q vector of covariates to include in the model (e.g. top 10 PCs, gender, etc). The estimated effects
#' for these covarites are not penalized, and there is no prior on the effects, the algorithm simply maximizes the ELBO w.r.t. the effects.
#' Do not include a constant column for the intercept, this is handled separately.
#' @param estimate_prior_variance indicates whether to estimate prior (initial estimate given by prior_variance. If prior_variance is
#' a vector, then a different prior variance is estimated for each of the L effects)
#' @param s_init a previous susie fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations (NOT CURRENTLY IMPLEMENTD IN LOGISTIC VERSION)
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Xr}{an n vector of fitted values, equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = rbinom(1000, 1, 1 / (1 + exp(-X %*% beta)))
#' res = susie_logistic(X,y,L=10)
#' coef(res)
#' plot(y,predict(res))
#'
#' @importFrom stats var
#' @importFrom utils modifyList
#'
#' @export
#'
susie_logistic = function(X, Y, L = min(10, ncol(X)), prior_variance = 1, prior_weights = NULL, 
                 standardize = TRUE, intercept = TRUE, Z = NULL, estimate_prior_variance = FALSE,
                 s_init = NULL, coverage=0.95, min_abs_corr=0.5,
                 max_iter = 100, tol = 1e-3, verbose = FALSE, track_fit = FALSE) {
  
  
  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix")) {
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  }
  p = ncol(X)
  n = nrow(X)
  if (is.null(prior_weights)) {
    prior_weights = rep(1 / p, p)
  }
  
  X = set_X_attributes(X, center = FALSE, scale = standardize) # in logistic case, no need to center (intercept handled differently)
  attr(X, "X2") = X^2 # store the matrix w/ all entries squared (compute once and store)
  
  
  # Check Z
  if (is.null(Z)) {
    if (intercept == TRUE) {
      Z = matrix(1, nrow = n, ncol = 1)
    } else {
      Z = matrix(0, nrow = n, ncol = 1) # if no intercept and no control covariates, set to 0
    }
  } else {
    col_variances = apply(Z, MARGIN = 2, var)
    if (any(col_variances == 0)) { # is constant column in Z matrix
      stop("Matrix 'Z' cannot have a constant column")
    }
    Z = cbind(matrix(1, nrow = n, ncol = 1), Z) # add intercept column
  }
  
  if (is.null(prior_weights)) {
    prior_weights = rep(1 / p, p)
  }
  
  # WILL CLEAN THIS UP LATER, NEED TO BE MORE CAREFUL IN GENERAL
  if (!is.null(s_init)) { # if starting from a previous fit, initialize there
    
    if (intercept) {
      if (ncol(Z) == 1) {
        delta = s_init$intercept
      } else {
        delta = c(s_init$intercept, s_init$delta)
      }
    } else {
      if (identical(Z, matrix(0, nrow = n, ncol = 1))) {
        delta = 0
      } else {
        delta = s_init$delta
      }
    }
    
    Alpha = t(s_init$alpha)
    Mu = t(s_init$mu)
    Sigma2 = t(s_init$mu2 - s_init$mu^2)
    
  } else { # initialize regularly
    # initialize: could think of something better
    #delta = glm(Y ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution
    if (all(Z == 0)) { # if no covariates and no intercept
      delta = 0
    } else {
      delta = glm(as.numeric(Y) ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution w/ just Z (no X)
    }
    Alpha = matrix(prior_weights, nrow = p, ncol = L)
    #Alpha = t(MCMCpack::rdirichlet(L, prior_weights)) # alternate initialization method
    Mu = matrix(0, nrow = p, ncol = L)
    Sigma2 = matrix(prior_variance, nrow = p, ncol = L, byrow = T)
    
  }
  
  xi = update_xi(X, Sigma2, Mu, Alpha, Z, delta)

  post_info = list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi, V = prior_variance)
  
  beta_post = post_info$Alpha * post_info$Mu
  
  elbo = numeric(max_iter + 1)
  elbo[1] = -Inf
  elbo[2] = calc_ELBO(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, post_info$V, prior_weights, Z, post_info$delta, post_info$xi)
  
  iter = 1
  while(abs(elbo[iter+1] - elbo[iter]) > tol) { # repeat until ELBO increase is negligible
    post_info = update_all(X, Y, post_info$xi, prior_weights, post_info$V, post_info$Sigma2, post_info$Mu, post_info$Alpha, Z, post_info$delta, estimate_prior_variance)
    
    iter = iter + 1
    elbo[iter + 1] = calc_ELBO(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, post_info$V, prior_weights, Z, post_info$delta, post_info$xi)
    
    if(verbose){
      print(paste("objective:", elbo[iter + 1]))
    }
    
    if (iter > max_iter) {
      stop("Maximum number of iterations reached")
    }
  }
  
  beta_post = colSums(t(post_info$Alpha) * t(post_info$Mu))
  lin_pred = (X %*% beta_post) + (Z %*% delta)
  
  # change output format to match the GLM version's output
  int = NULL
  delta = post_info$delta
  if (intercept == TRUE) {
    int = post_info$delta[1]
    delta = post_info$delta[-1]
  }
  if (length(delta) == 0) {
    delta = 0
  }
  
  # make V a vector, if it's a scalar (needed to make CSs)
  if (length(post_info$V) == 1) {
    post_info$V = rep(post_info$V, L)
  }
  
  # change output to match linear version
  s = list(alpha = t(post_info$Alpha), mu = t(post_info$Mu), mu2 = t(post_info$Sigma2 + post_info$Mu^2),
           Xr = rep(NA, n), KL = rep(NA, L), lbf = rep(NA, L), sigma2 = NA, V = post_info$V, pi = prior_weights, null_index = 0)
  class(s) = 'susie'
  

  elbo = elbo[2:(iter+1)] #remove trailing NAs and leading -Inf
  s$elbo <- elbo
  s$niter <- iter

  if(intercept){
    s$intercept = int
  } else {
    s$intercept = 0
  }

  s$delta = delta
  
  s$fitted = 1 / (1 + exp(-lin_pred))

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage = coverage, X = X, min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s, prune_by_cs = FALSE)
  }
  return(s)
}
