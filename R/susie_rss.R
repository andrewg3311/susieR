#' @title Bayesian sum of single-effect (susie) linear regression using z scores
#' @details Performs sum of single-effect (susie) linear regression with z scores.
#' The summary data required are the p by p correlation matrix R, the p vector z. The summary stats should come from the same individuals.
#' This function fits the regression model z = sum_l Rb_l + e, where e is N(0,residual_variance * R) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=prior_variance).
#' @param z a p vector of z scores.
#' @param R a p by p symmetric and positive semidefinite correlation matrix.
#' @param maf_thresh threshold for MAF
#' @param maf Minor Allele Frequency
#' @param L maximum number of non-zero effects
#' @param lambda fudge factor
#' @param prior_variance the prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance, a scaler between 0 and 1
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param restrict whether to restrict the resiudal variance between 0 and 1
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior
#' @param estimate_prior_method The method used for estimating prior variance, 'optim' or 'EM'
#' @param max_iter maximum number of iterations to perform
#' @param s_init a previous susie fit with which to initialize
#' @param intercept_value a value to assign to the intercept (since the intercept cannot be estimated from centered summary data). This
#' value will be used by coef.susie() to assign an intercept value, for consistency with the non-summary-statistic version of this function \code{susie}.
#' Set to NULL if you want coef.susie() not to include an intercept term (and so only return a p vector).
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' @param tol convergence tolerance based on alpha
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param check_input whether to skip the checks for R and z, the checks are:
#'
#' 1. Check whether R is positive semidefinite
#'
#' 2. Check whether z in space spanned by the non-zero eigenvectors of R
#'
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Rr}{an p vector of t(X) times fitted values, the fitted values equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#'
#' @export
susie_rss = function(z, R, maf_thresh=0, maf=NULL,
                     L=10, lambda = 0,
                     prior_variance=50,residual_variance=NULL,
                     r_tol = 1e-08,
                     prior_weights = NULL, null_weight = NULL,
                     restrict = TRUE,
                     estimate_residual_variance = TRUE,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = c("optim","EM"),
                     max_iter=100,s_init = NULL, intercept_value=0,
                     coverage=0.95, min_abs_corr=0.5,
                     tol=1e-3, verbose=FALSE, track_fit = FALSE, check_input = FALSE){

  warning('This method is under active development, so it should not be considered stable.')

  if(L > 1){
    warning('The maximum number of non-zero effects is greater than 1, this feature is experimental.')
  }
  estimate_prior_method <- match.arg(estimate_prior_method)

  # MAF filter
  if(!is.null(maf)){
    if(length(maf) != length(z)){
      stop(paste0('The length of maf does not agree with expected ', length(z)))
    }
    id = which(maf > maf_thresh)
    R = R[id, id]
    z = z[id]
  }

  # replace NA in z with 0
  if (any(is.na(z))){
    warning('z scores contain NA, it is replaced with 0.')
    z[is.na(z)] = 0
  }
  if(any(is.infinite(z))){
    stop('z contains infinite value.')
  }
  # Check NA in R
  if(any(is.na(R))){
    stop('R matrix contains NA.')
  }
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix.")

  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R)*(1-null_weight), ncol(R)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    R = cbind(rbind(R, 0),0)
    z = c(z, 0)
  }

  p = ncol(R)

  # Check input R.
  if(nrow(R) != length(z)) {
    stop(paste0('The dimension of R (', nrow(R), ' by ', ncol(R), ') does not agree with expected (', length(z), ' by ', length(z), ')'))
  }
  if(!is_symmetric_matrix(R)){
    stop('R is not a symmetric matrix.')
  }
  ## eigen decomposition
  semi_pd = check_semi_pd(R, r_tol)
  R = semi_pd$matrix

  ## check whether z in space spanned by the non-zero eigenvectors of R
  if(check_input){
    # Check whether R is positive semidefinite
    if(semi_pd$status == FALSE){
      stop('R is not a positive semidefinite matrix.')
    }
    # Check whether z in space spanned by the non-zero eigenvectors of R
    proj = check_projection(R, z)
    if(proj$status == FALSE)
      warning('z does not lie in the space of non-zero eigenvectors of R')
  }
  attr(R, 'eigen')$values = semi_pd$eigenvalues
  attr(R, 'd') <- diag(R)
  attr(R, 'scaled:scale') <- rep(1, length = p)

  # initialize susie fit
  s = init_setup_rss(p,L,prior_variance,residual_variance,prior_weights,null_weight,1)

  if (!missing(s_init)) {
    s = modifyList(s, s_init)
    s = init_finalize_rss(s, R=R)
  } else {
    s = init_finalize_rss(s)
  }

  # intialize elbo to NA
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  attr(R, 'lambda') = lambda
  Sigma = update_Sigma(R, s$sigma2, z) # sigma2*R + lambda I

  # alpha_new = s$alpha
  for(i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    # alpha_old = alpha_new
    s = update_each_effect_rss(R, z, s, Sigma, estimate_prior_variance,estimate_prior_method)
    # alpha_new = s$alpha

    if(verbose){
      print(paste0("before estimate sigma2 objective:",get_objective_rss(R, z, s)))
    }
    if(estimate_residual_variance){

      if(lambda == 0){
        tmp = s
        tmp$sigma2 = 1
        est_sigma2 = (1/sum(attr(R, 'eigen')$values!=0))* get_ER2_rss(R,z,tmp)
        if(est_sigma2 < 0){
          stop('Estimating residual variance failed: the estimated value is negative')
        }
        if(restrict){
          if(est_sigma2 > 1){
            est_sigma2 = 1
          }
        }
      }else{
        if(restrict){
          sigma2 = seq(0.1,1,by=0.1)
        }else{
          sigma2 = seq(0.1,6,by=0.1)
        }
        tmp = s
        obj = numeric(length(sigma2))
        for(j in 1:length(sigma2)){
          tmp$sigma2 = sigma2[j]
          obj[j] = Eloglik_rss(R, z, tmp)
        }
        est_sigma2 = sigma2[which.max(obj)]
      }

      s$sigma2 = est_sigma2

      if(verbose){
        print(paste0("after estimate sigma2 objective:", get_objective_rss(R, z, s)))
      }

      Sigma = update_Sigma(R, s$sigma2, z)
    }
    elbo[i+1] = get_objective_rss(R, z, s)

    if((elbo[i+1]-elbo[i])<tol) {
      s$converged = TRUE
      break;
    }
  }
  elbo = elbo[1:(i+1)] #remove trailing NAs
  s$elbo <- elbo
  s$niter <- i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in", max_iter, "iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$fitted = s$Rz

  s$X_column_scale_factors = attr(R,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(R)
    s$sets = susie_get_cs(s, coverage=coverage, Xcorr=R, min_abs_corr=min_abs_corr)
    s$pip = susie_get_pip(s, prune_by_cs = FALSE)
  }

  return(s)
}

update_Sigma = function(R, sigma2, z){
  Sigma = sigma2*R + attr(R, 'lambda') *diag(length(z))
  eigenS = attr(R, 'eigen')
  eigenS$values = sigma2*eigenS$values + attr(R, 'lambda')

  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  attr(Sigma, 'eigenS') = eigenS

  # Sigma^(-1) R_j
  attr(Sigma, 'SinvRj') = eigenS$vectors %*% (Dinv*attr(R, 'eigen')$values * t(eigenS$vectors))

  if(attr(R, 'lambda')==0){
    attr(Sigma, 'RjSinvRj') = attr(R, 'd')/sigma2
  }else{
    tmp = t(eigenS$vectors)
    attr(Sigma, 'RjSinvRj') = colSums(tmp * (Dinv*(attr(R, 'eigen')$values^2) * tmp))
  }

  return(Sigma)
}

