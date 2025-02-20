### Functions required to run mircoSLAM, a R package for performing metagenome-wide association studies (MWAS)
###
### The methodology involves fitting generalized linear mixed effects models (GLMMs) for trait associations.
###
### The GLMM fitting code (REML AI) is adapted from the SAIGE package (GNU GPL):
###   https://saigegit.github.io/SAIGE-doc/
### The methodology and equations are detailed in the supplementary
### material of the following paper:
### Zhou, W., Nielsen, J. B., Fritsche, L. G., Dey, R., Gabrielsen, M. E., Wolford, B. N., LeFaive, J., et al.
### "Efficiently Controlling for Case-Control Imbalance and Sample Relatedness in Large-Scale Genetic Association Studies."
### Nature Genetics 50 (2018): 1335–1341. https://doi.org/10.1038/s41588-018-0184-y.
###
### The saddle point approximation (SPA) code are from the SPAtest package (GNU GPL):
###   https://github.com/leeshawn/SPAtest


#' calculate_grm
#'
#' This function alculates genetic relatedness matrix (GRM) from sample-by-gene matrix
#'
#' @param sample_by_gene_matrix must have first column as sample names and other columns as gene names
#' @return grm genetic relatedness matrix for tau and beta test
#' @export
calculate_grm <- function(sample_by_gene_matrix) {
  check_gene_matrix(sample_by_gene_matrix)
  freq_mat_dist_man = parDist(as.matrix(sample_by_gene_matrix[, -1]), method = "manhattan") / (ncol(sample_by_gene_matrix[, -1]) - 1)
  freq_mat_dist_man = as.matrix(freq_mat_dist_man)
  freq_mat_grm_man = 1 - freq_mat_dist_man
  grm = as.matrix(freq_mat_grm_man)
  colnames(grm) = as.matrix(sample_by_gene_matrix[, 1])
  rownames(grm) = colnames(grm)
  return(grm)
}

gen_sigma <- function(W, var_vec, grm) {
  # Generate covariance matrix sigma for GLMM
  # sigma = phi * diag(1 / W) + tau * GRM
  # Adapted from getDiagOfSigma and get_sp_Sigma from SAIGE package  
  grm = as.matrix(grm)
  dtkin = W^-1 * (var_vec[1]) 
  new_grm = grm * var_vec[2]
  diag_new_grm = diag(new_grm) + dtkin # update diag
  diag(new_grm) = diag_new_grm
  new_grm[new_grm < 1e-4] = 1e-4 #check nothing is 0
  return(as.matrix(new_grm))
}

get_coef_inner <- function(Y, X, W, var_vec, grm) {
  # This function calculate fixed and random effect coefficients Equation 3 and Equation 4
  # Y is working vector Y=alpha X + b
  # X: sample-by-covariate matrix
  # W coefficient of variation
  # var_vec is the variance component [tau, phi]
  # Adapted from cpp::getCoefficients in SAIGE package
  simga = gen_sigma(W, var_vec, grm) # V
  Y = as.vector(Y)
  sigmai = Rfast::spdinv(simga)

  if(all(is.na(sigmai))){
    print("While fitting Sigma, a transformation of the GRM became singular or was not symmetrical, this can happen when many columns of the GRM are highly correlated to the y, or when there is not enough varaiblity in the GRM.")
    stop()
  }

  sigmai_Y = sigmai %*% Y # V^-1 Y
  sigmai_X = sigmai %*% X # V^-1 X
  cov_var = Matrix::solve(forceSymmetric(t(X) %*% sigmai_X), sparse = TRUE, tol = 1e-10) # (Xt V^-1 X)^-1
  sigmai_Xt = t(sigmai_X) # t(V^-1 X)
  sigmaiXtY = sigmai_Xt %*% Y # XtV^-1Y
  # Equation 3
  alpha = cov_var %*% sigmaiXtY # (Xt V X)^-1 XtVY
  epsilon = var_vec[1] *(t(sigmai_Y) - t(sigmai_X %*% alpha)) / as.vector(W) # phi to act on W
  # Updated working response
  eta = as.vector(Y - epsilon) # Y-var_vec \sigma (Y-X\alpha)
  b = eta - X %*% alpha
  coef_list = list("sigmai_Y" = sigmai_Y, "sigmai_X" = sigmai_X, "cov_var" = cov_var, "alpha" = alpha, "eta" = eta, "b" = b, "epsilon" = epsilon)
  return(coef_list)
}

get_alpha <- function(y, X, var_vec, grm, family, alpha0, eta0, offset, maxiter, tol.coef = tol) {
  # Adapted from Get_Coef from SAIGE package
  mu = family$linkinv(eta0)
  mu_eta = family$mu.eta(eta0)

  Y = eta0 - offset + (y - mu) / mu_eta
  sqrt_W = mu_eta / sqrt(family$variance(mu)) # weight matrix
  W = sqrt_W^2

  for (i in 1:maxiter) {
    alpha.obj = get_coef_inner(Y, X, W, var_vec, grm)

    alpha = as.matrix(alpha.obj$alpha)
    eta = as.matrix(alpha.obj$eta + offset)

    mu = family$linkinv(eta)
    mu_eta = family$mu.eta(eta)

    Y = eta - offset + (y - mu) / mu_eta
    sqrt_W = mu_eta / sqrt(family$variance(mu))
    W = sqrt_W^2

    if (max(abs(alpha - alpha0) / (abs(alpha) + abs(alpha0) + tol.coef)) < tol.coef) {
      break
    }
    alpha0 = alpha
  }

  alpha_result = list("Y" = Y, "alpha" = alpha, "eta" = eta, "W" = W, "cov_var" = alpha.obj$cov_var, 
                      "sqrt_W" = sqrt_W, "sigmai_Y" = alpha.obj$sigmai_Y, "sigmai_X" = alpha.obj$sigmai_X,
                      "mu" = mu, "eta_2" = alpha.obj$eta_2, "b" = alpha.obj$b)
  return(alpha_result)
}

get_AI_score <- function(Y, X, grm, W, var_vec, sigmai_Y, sigmai_X, cov_var) {
  # Get AI score used to estimate var_vec from the supplementary materials of the SAIGE paper.
  # Adapted from getAIScore from SAIGE package
  sigma = gen_sigma(W, var_vec, grm)

  sigmai =  Rfast::spdinv(sigma)
  sigmai_Xt = t(sigmai_X) # transpose sigma inverse times X
  
  P = sigmai - sigmai_X %*% cov_var %*% sigmai_Xt # solve for P
  
  PY1 = P %*% Y # \hat{Y}-\hat(X) (Xt V X)^-1 PY
  
  APY = grm %*% PY1 # grm (\hat{Y}-\hat(X) (Xt V X)^-1)
  YPAPY = t(PY1) %*% APY # dot product
  YPAPY = YPAPY[1] 
  
  PA = P %*% grm
  trace_P_grm = sum(diag(PA))
  
  score1 = YPAPY - trace_P_grm
  PAPY = P %*% APY
  
  AI = (t(PAPY) %*% APY) # AI=t(Y)%*%P%*%grm%*%P%*%grm%*%P%*%Y
  return(list("YPAPY" = YPAPY, "PY" = PY1, "trace_P_grm" = trace_P_grm, "score1" = score1, "AI" = AI[1]))
}

get_AI_score_quant <- function(Y, X, grm, W, var_vec, sigmai_Y, sigmai_X, cov_var) {
  ## Compute AI 2*2 matrix for quantative traits
  ## https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0184-y/MediaObjects/41588_2018_184_MOESM1_ESM.pdf
  # Adapted from getAIScore_q from SAIGE package
  
  # Compute sigma
  sigma = gen_sigma(W, var_vec, grm)
  # Compute projection matrix P
  sigmai = Rfast::spdinv(sigma)
  sigmai_Xt = t(sigmai_X) # transpose X
  P = sigmai - sigmai_X %*% cov_var %*% sigmai_Xt
  PY1 = P %*% Y # \hat{Y}-\hat(X) (Xt V X)^-1 PY

  APY = grm %*% PY1
  YPAPY = t(PY1) %*% APY # dot product
  YPAPY = YPAPY[1] 
  
  PA = P %*% grm
  trace_P_grm = (sum(sigmai * grm) - sum(sigmai_X * crossprod(grm, t(cov_var %*% sigmai_Xt))))

  score1 = YPAPY - trace_P_grm # score 1

  ## Quantitative
  # compute wPY, this is equivalent to scale PY by the inverse of the phenotypic residual variance
  wPY = PY1 / (W)
  YPwPY = t(PY1) %*% wPY
  YPwPY = YPwPY[1]

  # compute SCORE statsitcs
  diag_P = diag(P) / (W)
  trace_PW = sum(diag_P)
  score0 = YPwPY - trace_PW # score 0 

  score_vector = as.matrix(c(score0[1], score1[1]))
  
  PwPY = P %*% wPY
  PAPY = P %*% APY
  
  AI_11 = (t(PAPY) %*% APY)
  AI_00 = (t(PwPY) %*% wPY)
  AI_01 = (t(PAPY) %*% wPY)
  AI_mat = matrix(c(AI_00[1], AI_01[1], AI_01[1], AI_11[1]), 2, 2)
  
  return(list("YPAPY" = YPAPY, "PY" = PY1, "YPwPY" = YPwPY, "trace_P_grm" = trace_P_grm, "trace_PW" = trace_PW, "AI" = AI_mat, "score_vector" = score_vector))
}

fit_vars <- function(Y_vec, X_mat, grm, w_vec, var_vec, sigmai_Y, sigmai_X, cov_var, tol, quant = FALSE, verbose, write_log) {
  # Fitting process for tau and phi
  # Adapted from fitglmmaiRPCG_q + fitglmmaiRPCG from SAIGE package
  if (!quant) {
    re.AI = get_AI_score(Y_vec, X_mat, grm, w_vec, var_vec, sigmai_Y, sigmai_X, cov_var)
    Dtau = re.AI$score1 / re.AI$AI
    var_vec0 = var_vec
    var_vec[2] = var_vec0[2] + Dtau # update tau
    step = 1.0
    while (var_vec[2] < 0) {
      step = step * 0.5
      var_vec[2] = var_vec0[2] + step * Dtau
    }
  } else {
    re.AI = get_AI_score_quant(Y_vec, X_mat, grm, w_vec, var_vec, sigmai_Y, sigmai_X, cov_var)
    Dtau = Matrix::solve(re.AI$AI, re.AI$score_vector) #=-
    var_vec0 = var_vec
    var_vec = var_vec0 + Dtau
    step = 1.0
    while (any(var_vec < 0)) { # if var_vec is less than zero step through
      step = step * 0.5
      var_vec = var_vec0 + step * Dtau
    }
  }

  if (var_vec[1] < tol) {
    print("Warning! The first variance component parameter estimate is set at the tolarance")
    var_vec[1] = tol * 10
  }
  if (var_vec[2] < tol) {
    var_vec[2] = 0
  }
  return(list("var_vec" = var_vec))
}


simulate_tau_inner <- function(glm.fit0, grm, species_id = "s_id", tau0, phi0) {
  family_to_fit = glm.fit0$family
  data_new = glm.fit0$data
  formulate_to_fit = glm.fit0$formula

  data_new_shuffled = data_new[sample(nrow(data_new)), ]
  data_new_shuffled$sample_name = rownames(grm)
  
  refit0 = glm(formulate_to_fit, data = data_new_shuffled, family = family_to_fit)

  fit.glmm <- tryCatch(
    fit_tau_test(
      refit0, grm,
      tau0 = tau0, phi0 = phi0,
      verbose = FALSE,
      species_id = species_id,
      log_file = NA
    ),
    error = function(e) {
      message("Error in fitting GLMM: ", e$message)
      return(NULL)
    }
  )

  # Initialize default results
  tau = 0
  t = 0
  phi = 0
  if (length(fit.glmm$t) > 0) {
    t = sum(fit.glmm$b^2, na.rm = TRUE) / length(fit.glmm$sample_names)
    tau = fit.glmm$var_vec[2]
    phi = fit.glmm$var_vec[1]
  } else {
    message("GLMM fitting failed or produced no results.")
  }

  return(data.frame("tau" = tau, "t = t", "phi" = phi))
}


for_beta <- function(mu, mu2, y, X) {
  ## score test for var_vec test used to fit beta test
  ## inputs: fitted mu, original y, covariates X
  # Adapted from ScoreTest_NULL_Model_binary from SAIGE package
  V = as.vector(mu2)
  res = as.vector(y - mu)
  
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = Matrix::solve(XVX)
  
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  
  S_a = colSums(X * res)
  for_beta = list("XV" = XV, "XVX" = XVX, "XXVX_inv" = XXVX_inv, "XVX_inv" = XVX_inv, "S_a" = S_a, "XVX_inv_XV" = XVX_inv_XV, "V" = V)
  return(for_beta)
}


fit_beta_one_gene <- function(glmm_fit, glm_fit0, grm, one_gene_df, SPA = FALSE) {
  # Function: fit_beta_one_gene
  # Args:
  # - glmm_fit: Output from the fit_tau_test function
  # - glm_fit0: Output from the glm model for data in the same order
  # - grm: Genetic relatedness matrix
  # - one_gene_df: The presence or absence for ONE gene in a long format
  #            Must have columns for sample_name, gene_id, and gene_value!
  # - SPA: Whether to run saddle point approximation for p-values (mostly for binomial case)
  # Returns:
  # - List of values for each gene examined
  # Adapted from SAIGE package

  check_grm(grm, glm_fit0, verbose = TRUE)
  check_beta(glmm_fit, glm_fit0, grm, gene_df, SPA)
  
  forbeta.obj = glmm_fit$forbeta.obj
  family = glm_fit0$family
  eta = glmm_fit$linear_predictors
  mu = glmm_fit$fitted_values
  mu_eta = family$mu.eta(eta)
  sqrt_W = mu_eta / sqrt(glm_fit0$family$variance(mu))
  W1 = sqrt_W^2 

  G0 = as.vector(one_gene$gene_value)
  G_tilde = G0 - glmm_fit$forbeta.obj$XXVX_inv %*% (glmm_fit$forbeta.obj$XV %*% G0) # G1 is X adjusted
  Y = eta + (glmm_fit$y - mu) / mu_eta
  t_score = t(G_tilde) %*% (glmm_fit$y - mu)
  m1 = sum(mu * G_tilde)
  var1 = sum(W1 * G_tilde^2)
  t_adj_2 = (t_score^2) / var1
  beta = t_score / var1
  pval = (pchisq(t_adj_2, lower.tail = FALSE, df = 1, log.p = FALSE))
  z = (qnorm(pval / 2, log.p = F, lower.tail = F))
  se_beta = abs(beta) / sqrt(abs(z))
  qtilde = t_score + m1

  # Generate the result data frame
  ret_basic = data.frame(
        "gene_id" = colnames(one_gene_df)[2], 
        "tau" = as.numeric(glmm_fit$var_vec[2]), 
        "cor_to_y" = cor(G0, glmm_fit$y), 
        "cor_to_b" = cor(as.numeric(glmm_fit$b), G0), 
        "z" = z, 
        "var1" = var1)
  
  if (SPA) {
    if (var1 < 0) {
      ret = cbind(ret_basic, data.frame("beta" = NA, "se_beta" = NA, "t_adj" = NA, "SPA_pvalue" = NA, 
                                      "spa_score" = NA, "SPA_zvalue" = NA, "pvalue_noadj" = NA, "converged" = NA))
    } else {
      out1 = saddle_prob(q = qtilde, mu = mu, g = G_tilde, var1, cutoff = 2, log.p = FALSE)
      ret = cbind(ret_basic, data.frame("beta" = beta, "se_beta" = se_beta, "t_adj" = t_adj_2, "SPA_pvalue" = out1$p.value, 
                  "spa_score" = out1$score, "SPA_zvalue" = out1$z_value, "pvalue_noadj" = out1$pvalue_noadj, "converged" = out1$converged))
    }
  } else {
    # no SPA
    ret = cbind(ret_basic, data.frame("beta" = beta, "se_beta" = se_beta, "t_adj" = t_adj_2, "SPA_pvalue" = NA, 
                  "spa_score" = NA, "SPA_zvalue" = NA, "pvalue_noadj" = pval,  "converged" = NA)) #<-- DOUBLE CHECK converged
  }
  return(ret)
}


#' fit_tau_test
#'
#' Fit the base model for population structure finding random effects
#'
#' @param glm.fit0 glm model. Model output with no sample relatedness accounted for
#' @param grm Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param species_id for tracking species
#' @param tau0 inital tau estimate (variance on the population strucutre of genetic relatedness)
#' @param phi0 inital phi estimate (variance on means of outputs will be 1 for logit or binonimal dispersion)
#' @param maxiter maximum iterations to fit the glmm model
#' @param verbose whether outputting messages in the process of model fitting
#' @param log_file log file to write to
#' @return model output for the tau test on population structure 
#' @export
fit_tau_test <- function(glm.fit0, grm, species_id, tau0 = 1, phi0= 1, maxiter = 100, verbose = TRUE, tol = .0001, log_file = NA) {
  # Fits the null generalized linear mixed model for a binary trait
  # Adapted from glmmkin.ai_PCG_Rcpp_Binary from SAIGE package
  # Args:
  #  glm.fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  grm: Genetic Relatedness Matrix  in the same sample order as glm.fit0!
  #  species_id: Species ID of the species for record
  #  tau0: initial values for the variance component parameter tau
  #  phi0: initial values for the variance component parameter phi (only for quant)
  #  maxiter: maximum iterations to fit the glmm model
  #  verbose: whether outputting messages in the process of model fitting
  #  tol: tolance for varaince estimates
  #  log_file: whether to write to a log file and what that would be
  # Returns:
  #  model object pop.struct.glmm

  log_file = ifelse(is.na(log_file), file.path(tempdir(), "microSLAM.log"), log_file)

  t_begin = proc.time()
  
  log_open(log_file, logdir = F, show_notes = F)
  log_print("====fit_tau_test::start::====", console = verbose)
  log_print(paste("====Fixed-effect coefficients:"), console = verbose)
  log_print(glm.fit0, console = verbose)
  
  check_inputs_tau(glm.fit0, grm, species_id, tau0, phi0, maxiter, tol, verbose, write_log, log_file)

  y = glm.fit0$y
  n = length(y)
  X = model.matrix(glm.fit0)
  Xorig = model.matrix(glm.fit0)

  offset = glm.fit0$offset
  if (is.null(offset)) {
    offset = rep(0, n)
  }

  family = glm.fit0$family
  eta = glm.fit0$linear.predictors
  mu = glm.fit0$fitted.values
  mu_eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu) / mu_eta
  alpha0 = glm.fit0$coef
  eta0 = eta

  sample_ids = colnames(grm)
  
  ### check for quantitiative or not from glm obj
  if (family$family %in% c("poisson", "binomial")) {
    phi0 = 1
    quant = FALSE
  } else {
    quant = TRUE
  }
 
  var_vec0 = c(phi0, tau0) ## var_vec is a vector of phi and tau for the variance estimates for tau test
  if (var_vec0[1] <= 0) {
    stop("\nERROR! The first variance component parameter estimate is 0\n")
  }

  var_vec=var_vec0
  log_print(paste("====get_alpha::initial var_vec:: ", paste(round(var_vec, 4), collapse = ", "), sep = ""), console = verbose)

  ## Fixed-effect coefficients => alpha.obj$alpha
  alpha.obj = get_alpha(y, X, var_vec0, grm, family, alpha0, eta0, offset, maxiter = maxiter, tol.coef = tol)
 
  ### Initial Update of variance components phi and tau 
  if (quant) {
    re = get_AI_score_quant(alpha.obj$Y, X, grm, alpha.obj$W, var_vec0, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
    var_vec[2] = max(0, as.numeric(var_vec0[2] + var_vec0[2]^2 * (re$YPAPY - re$trace_P_grm) / n))
    var_vec[1] = max(tol*10, as.numeric(var_vec0[1] + var_vec0[1]^2 * (re$YPwPY - re$trace_PW) / n))
  } else {
    re = get_AI_score(alpha.obj$Y, X, grm, alpha.obj$W, var_vec0, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
    var_vec[2] = max(0, as.numeric(var_vec0[2] + var_vec0[2]^2 * ((re$YPAPY - re$trace_P_grm)) / n)) # first fit vars
  }
  log_print(paste("====get_AI_score::var_vec::", paste(round(var_vec, 4), collapse = ", ")), console = verbose)

  t_begin_fit_tau = proc.time()

  tol_limt = FALSE
  list_of_log_messages <- list()
  for (i in seq_len(maxiter)) {
    alpha0 = alpha.obj$alpha
    var_vec0 = var_vec
    eta0 = eta
    rss_0 = sum((y - mu)^2)
    
    ## get alpha
    alpha.obj = get_alpha(y, X, var_vec, grm, family, alpha0, eta0, offset, maxiter = maxiter, tol.coef = tol)
    ## fit tau and phi
    fit.obj = fit_vars(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var, tol = tol, verbose = verbose, write_log = write_log, quant = quant)
    
    ## update all params
    var_vec = as.numeric(fit.obj$var_vec)
    tau = as.numeric(fit.obj$var_vec[2])
    cov_var = alpha.obj$cov_var
    alpha = alpha.obj$alpha
    eta = alpha.obj$eta
    Y = alpha.obj$Y
    mu = alpha.obj$mu
    res = y - mu

    list_of_log_messages[[i]] = paste("  ===fit_vars::iteration ", i, "::var_vec ", paste(round(var_vec, 4), collapse = ", "), sep = "")
    
    #if(sum(res^2)/n < tol) break
    if (var_vec[1] <= 0) {
      stop("\nERROR! The first variance component parameter estimate is 0\n")
    }
    if (var_vec[2] == 0) break

    var_condition = max(abs(var_vec - var_vec0) / (abs(var_vec) + abs(var_vec0) + tol)) < tol 
    if (var_vec[1] <= tol*10 & var_vec0[1] <= tol*10) {
      tol_limt = TRUE
      print("Warning! The first variance component parameter estimate is set at the tolarance breaking iterations")
      break
    }

    if (var_condition) break
    if (max(var_vec) > tol^(-2)) {
      tol_limt = TRUE
      i = maxiter
      print("Warning! max(var_vec) > tol^(-2)")
      break
    }
  }
  t_end_fit_tau = proc.time()

  for (msg in list_of_log_messages) {
    log_print(msg, console = verbose)
  }
  
  if (max(var_vec) > tol^(-2) | i == maxiter) {
    log_print("Model not converged", console = verbose)
  }
  
  alpha.obj = get_alpha(y, X, var_vec, grm, family, alpha, eta, offset, maxiter = maxiter, tol.coef = tol)
  
  cov_var = alpha.obj$cov_var
  alpha = alpha.obj$alpha
  names(alpha) = names(glm.fit0$coefficients)
  eta = alpha.obj$eta
  Y = alpha.obj$Y
  mu = alpha.obj$mu

  fit.final = fit_vars(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var, tol = tol, verbose = verbose, write_log = write_log, quant = quant)
  var_vec = as.numeric(fit.final$var_vec)
  names(var_vec) = c("phi", "tau")

  if (FALSE) {
    # This was the same with the initialization step. But, we should update the tau and phi just like the for loop.
    if (quant) {
      if(tol_limt){
        fit.final = get_AI_score_quant(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
        var_vec[2] = max(0, var_vec0[2] + var_vec0[2]^2 * (fit.final$YPAPY - fit.final$trace_P_grm) / n)
        var_vec[1] = max(tol*10, var_vec0[1] + var_vec0[1]^2 * (fit.final$YPwPY - fit.final$trace_PW) / n)
        i = maxiter
        
      }else{
        fit.final = get_AI_score_quant(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
        var_vec[2] = max(0, var_vec0[2] + var_vec0[2]^2 * (fit.final$YPAPY - fit.final$trace_P_grm) / n)
        var_vec[1] = max(tol*10, var_vec0[1] + var_vec0[1]^2 * (fit.final$YPwPY - fit.final$trace_PW) / n)
      }
      
    } else {
      fit.final = get_AI_score(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
      var_vec[2] = max(0, as.numeric(var_vec0[2] + var_vec0[2]^2 * ((fit.final$YPAPY - fit.final$trace_P_grm)) / n)) # tau + Dtau 
    }

  }
  
  if (quant) {
    mu2 = rep(1 / var_vec[1], length(y))
  } else {
    mu2 = mu * (1 - mu)
  }
  forbeta.obj = for_beta(mu, mu2, y, Xorig)
  
  log_print(paste("====fit_tau_phi::final var_vec::", paste(round(var_vec, 4), collapse = ", ")), console = verbose)
  log_print(paste("====fit_tau_phi elapsed time:"), console = verbose)
  log_print(t_end_fit_tau - t_begin_fit_tau, console = verbose)
  
  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu
  rss = sum(res^2)
  
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
  ss = sum((y - mean(y))^2)
  var_fixed = var(Xorig %*% alpha)
  var_random = var(as.vector(alpha.obj$b))
  var_error = var(res)

  if (!quant) {
    model_metrics = list(
      "S" = sd(res),
      "AUC" = suppressMessages(suppressWarnings(auc(glm.fit0$y, as.vector(mu)))),
      "R-sq" = (1 - rss / ss),
      "R-sq(marginal)" = var_fixed / (var_fixed + var_random + var_error),
      "r-sq(conditional)" = (var_fixed + var_random) / (var_fixed + var_random + var_error)
    )
  } else {
    model_metrics = list(
      "S" = sd(res),
      "R-sq" = (1 - rss / ss),
      "R-sq(marginal)" = var_fixed / (var_fixed + var_random + var_error),
      "r-sq(conditional)" = (var_fixed + var_random) / (var_fixed + var_random + var_error)
    )
  }

  glmm_result = list(
    "tau" = var_vec[2],
    "var_vec" = var_vec,
    "coefficients" = alpha, 
    "b" = alpha.obj$b, 
    "t" = sum(alpha.obj$b^2)/length(sample_ids),
    "linear_predictors" = eta, 
    "linear_model" = Xorig %*% alpha + alpha.obj$b,
    "fitted_values" = mu, 
    "var_mu" = mu2, 
    "Y" = Y, 
    "residuals" = res,
    "cov_var" = cov_var, 
    "converged" = converged,
    "sample_names" = sample_ids,
    "forbeta.obj" = forbeta.obj,
    "y" = y, 
    "X" = Xorig,
    "trait_type" = glm.fit0$family,
    "formula" = paste0(glm.fit0$formula, " + b"),
    "iter_finised" = i,
    "model_metrics" = model_metrics, 
    "species_id" = species_id #,#grm = grm #<- we dont' need to carry GRM around
  )
  class(glmm_result) = "pop.struct.glmm"
  t_end = proc.time()

  log_print(paste("====total elapsed time:"), console = verbose)
  log_print(t_end - t_begin, console = verbose)
  log_print("====fit_tau_test::end::====", console = verbose)
  for (msg in summary.pop.struct.glmm(glmm_result)) {
    log_print(msg, console = verbose)
  }
  log_close()

  return(glmm_result)
}


#' run_tau_test
#'
#' take output from population structure test and test how significant the tau is for that set of data
#'
#' @param glm.fit0 glm model. Model output with no sample relatedness accounted for
#' @param grm Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param n_tau number of tau to simulate
#' @param species_id species id for bacterial species
#' @param tau0 starting tau
#' @param phi0 starting phi
#' @return df of values of T for tau for different permutations of the covarites matrix
#' @export
run_tau_test <- function(glm.fit0, grm, n_tau, species_id = "s_id", tau0, phi0, seed=1) {
  set.seed(seed)
  list_of_tau = lapply(seq(1, n_tau), function(x) simulate_tau_inner(glm.fit0, grm, species_id = species_id, tau0, phi0))
  df_of_tau = do.call(rbind, list_of_tau)
  return(df_of_tau)
}


#' fit_beta
#'
#' Fit a beta for each genes including the population structure model with the random effects
#'
#' @param pop.struct.glmm output of fit_tau_test; GLMM of species with grm accounted for
#' @param glm.fit0 glm model. Model output with no sample relatedness accounted for
#' @param grm Genetic Relatedness Matrix (from calculate_grm or user) NxN matrix of sample relatedness
#' @param sample_by_gene_df long data frame with, gene_id, sample_name, and gene_value
#' @param SPA whether to run Saddle point approximation for pvalues (will slow down output)
#' @return dataframe of beta estimates for all genes tested
#' @export
fit_beta <- function(pop.struct.glmm, glm.fit0, grm, sample_by_gene_df, SPA = FALSE) {
  # Args:
  # pop.struct.glmm: output from fit_tau_test
  # glm.fit0: output from glm model for data in same order
  # grm: genetic relatedness matrix
  # sample_by_gene_df: the presance or absence for each gene in a long format 
  # sample_by_gene_df must have column for sample; column for gene_id; column for gene_value
  # SPA: whether to run saddle point approximation for pvalues (mostly for binomial case)
  # Returns:
  # list of values for each gene examined
  # Adapted from SAIGE package
  

  list_of_samples_from_glmm = glmm_fit$sample_names
  stopifnot(all.equal(rownames(grm), list_of_samples_from_glmm))
  
  sample_by_gene_df <- sample_by_gene_df %>% 
    mutate(sample_name = factor(sample_name, levels = list_of_samples_from_glmm)) %>%
    arrange(sample_name)
  stopifnot(all.equal(rownames(sample_by_gene_df), as.vector(rownames(grm))))
  
  t_begin = proc.time()
  iter = 1
  list_of_genes <- setdiff(colnames(sample_by_gene_df), "sample_name")
  list_of_rets <- list()
  
  for (my_gene in list_of_genes) {
    if (iter %% 1000 == 0) {
      cat(paste("number of genes done ", iter, "\n"))
      cat("time past:")
      t_now = proc.time()
      cat(t_now - t_begin)
      cat("\n")
    }
    
    one_gene_df <- sample_by_gene_df %>% select(sample_name, all_of(my_gene))
    ret <- fit_beta_one_gene(pop.struct.glmm, glm.fit0, grm, one_gene_df, SPA=SPA)
    list_of_rets[[my_gene]] <- ret
    iter = iter + 1
  }
  cat("total time past:")
  t_end = proc.time()
  cat(t_end - t_begin)
  cat("\n")

  genes_test_df <- bind_rows(list_of_rets)

  return(genes_test_df)
}


#' summary.pop.struct.glmm 
#'
#' Summarize the output of fit tau test
#'
#' @param x a pop.struct.glmm objecct the output of fit_tau_test; GLMM of species with grm accounted for
#' @export summary.pop.struct.glmm
#' @export 
summary.pop.struct.glmm <- function(x, ...) {
  cat("Species ID: ", x$species_id, "\n")
  cat("Formula: ", x$formula, "\n")
  cat("family: ", paste(x$trait_type[1:2]), "\n")
  cat("Fixed-effect covariates estimates: \n", names(x$coefficients), "\n", round(x$coefficients,3), "\n")
  cat("Converged: ", x$converged, "\n")
  cat("Number of iterations:",x$iter_finised,"\n")
  cat("Tau: ", round(x$var_vec[2],3), "\n")
  cat("Phi: ", round(x$var_vec[1],3), "if logit or binomail should be 1", "\n")
  cat("T value of tau:", round(x$t,3),"\n")
  cat("Number of Samples:", length(x$sample_names),"\n")
}

#' ... all the usual documentation for summary() ...
#' @export
summary <- function(x, ...) {
  UseMethod("summary")
}


filter_pop_obj <- function(pop.struct.glmm, sample_indexs) {
  pop.struct.glmm$gene_value = sample_indexs$gene_value
  pop.struct.glmm$residuals = pop.struct.glmm$residuals[sample_indexs$index]
  pop.struct.glmm$b = pop.struct.glmm$b[sample_indexs$index]
  pop.struct.glmm$linear.predictors = pop.struct.glmm$linear.predictors[sample_indexs$index]
  pop.struct.glmm$fitted.values = pop.struct.glmm$fitted.values[sample_indexs$index]
  pop.struct.glmm$forbeta.obj$XXVX_inv = pop.struct.glmm$forbeta.obj$XXVX_inv[sample_indexs$index, ]
  pop.struct.glmm$forbeta.obj$XV = pop.struct.glmm$forbeta.obj$XV[, sample_indexs$index]
  pop.struct.glmm$sigma = pop.struct.glmm$sigma[sample_indexs$index, sample_indexs$index]
  pop.struct.glmm$X = pop.struct.glmm$X[sample_indexs$index, ]
  pop.struct.glmm$y = pop.struct.glmm$y[sample_indexs$index]
  pop.struct.glmm$sigmai_X = pop.struct.glmm$sigmai_Xt[sample_indexs$index, ]
  pop.struct.glmm$sigmai_Y = pop.struct.glmm$sigmai_Y[sample_indexs$index]
  return(pop.struct.glmm)
}


########
saddle_prob <- function(q, mu, g, var1, cutoff = 2, log.p = FALSE) {
  #### taken from ‘SPAtest’ with a few changes for use case
  m1 = sum(mu * g)
  var1 = sum(mu * (1 - mu) * g^2)
  p1 = NULL
  p2 = NULL

  score = q - m1

  qinv = -sign(q - m1) * abs(q - m1) + m1
  t_adj_sq = ((q - m1)^2) / var1
  pval_noadj = pchisq(t_adj_sq, lower.tail = FALSE, df = 1, log.p = FALSE)
  converged = TRUE

  if (is.na(abs(q - m1)) || is.na(var1) || var1 < 0) {
    converged = FALSE
    pval = pval_noadj
  }
  if (isTRUE(abs(q - m1) / sqrt(var1) < cutoff)) {
    pval = pval_noadj
  } else {
    out_uni1 = get_root_K1(0, mu = mu, g = g, q = q)
    out_uni2 = get_root_K1(0, mu = mu, g = g, q = qinv)
    if (out_uni1$converged == TRUE && out_uni2$converged == TRUE && converged == TRUE) {
      p1 = tryCatch(get_saddle_prob(out_uni1$root, mu, g, q, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval_noadj - log(2))
        } else {
          return(pval_noadj / 2)
        }
      })
      p2 = tryCatch(get_saddle_prob(out_uni2$root, mu, g, qinv, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval_noadj - log(2))
        } else {
          return(pval_noadj / 2)
        }
      })
      if (log.p) {
        pval = add_logp(p1, p2)
      } else {
        pval = abs(p1) + abs(p2)
      }
      converged = TRUE
    } else {
      cat("Error_Converge")
      pval = pval_noadj
      converged = FALSE
    }
  }
  z_value = qnorm(pval / 2, log.p = F, lower.tail = F) * sign(q - m1)
  if (pval != 0 && pval_noadj / pval > 10^3) {
    return(saddle_prob(q, mu, g, var1, cutoff = cutoff * 2, log.p = log.p))
  } else if (pval == 0) {
    return(list(p.value = pval, pvalue_noadj = pval_noadj, z_value = z_value, converged = FALSE, score = score))
  } else {
    return(list(p.value = pval, pvalue_noadj = pval_noadj, z_value = z_value, converged = converged, score = score))
  }
}

get_root_K1 <- function(init, mu, g, q, m1, tol = .0001, maxiter = 1000) {
  #### taken from ‘SPAtest package’
  g_pos = sum(g[which(g > 0)])
  g_neg = sum(g[which(g < 0)])
  if (q >= g_pos || q <= g_neg) {
    return(list(root = Inf, n_iter = 0, converged = TRUE))
  } else {
    t = init
    K1_eval = K1_adj(t, mu, g, q)
    prevJump = Inf
    rep = 1
    repeat    {
      K2_eval = K2(t, mu, g)
      tnew = t - K1_eval / K2_eval
      if (is.na(tnew)) {
        conv = FALSE
        break
      }
      if (abs(tnew - t) < tol) {
        conv = TRUE
        break
      }
      if (rep == maxiter) {
        conv = FALSE
        break
      }

      newK1 = K1_adj(tnew, mu, g, q)
      if (sign(K1_eval) != sign(newK1)) {
        if (abs(tnew - t) > prevJump - tol) {
          tnew = t + sign(newK1 - K1_eval) * prevJump / 2
          newK1 = K1_adj(tnew, mu, g, q)
          prevJump = prevJump / 2
        } else {
          prevJump = abs(tnew - t)
        }
      }

      rep = rep + 1
      t = tnew
      K1_eval = newK1
    }
    return(list(root = t, n_iter = rep, converged = conv))
  }
}

Korg <- function(t, mu, g) {
  #### From ‘SPAtest’ package
  n.t = length(t)
  out = rep(0, n.t)

  for (i in 1:n.t) {
    t1 = t[i]
    temp = log(1 - mu + mu * exp(g * t1))
    out[i] = sum(temp)
  }
  return(out)
}

get_saddle_prob <- function(zeta, mu, g, q, log.p = FALSE) {
  #### From ‘SPAtest' package
  k1 = Korg(zeta, mu, g)
  k2 = K2(zeta, mu, g)

  if (is.finite(k1) && is.finite(k2)) {
    temp1 = zeta * q - k1


    w = sign(zeta) * (2 * temp1)^{
      1 / 2
    }
    v = zeta * (k2)^{
      1 / 2
    }

    Z.test = w + 1 / w * log(v / w)


    if (Z.test > 0) {
      pval = pnorm(Z.test, lower.tail = FALSE, log.p = log.p)
    } else {
      pval = -pnorm(Z.test, lower.tail = TRUE, log.p = log.p)
    }
  } else {
    if (log.p) {
      pval = -Inf
    } else {
      pval = 0
    }
  }

  return(pval)
}

K1_adj <- function(t, mu, g, q) {
  #### From ‘SPAtest' package
  n.t = length(t)
  out = rep(0, n.t)

  for (i in 1:n.t) {
    t1 = t[i]
    temp1 = (1 - mu) * exp(-g * t1) + mu
    temp2 = mu * g
    out[i] = sum(temp2 / temp1) - q
  }
  return(out)
}

K2 <- function(t, mu, g) {
  n.t = length(t)
  out = rep(0, n.t)

  for (i in 1:n.t) {
    t1 = t[i]
    temp1 = ((1 - mu) * exp(-g * t1) + mu)^2
    temp2 = (1 - mu) * mu * g^2 * exp(-g * t1)
    out[i] = sum(temp2 / temp1, na.rm = TRUE)
  }
  return(out)
}



########
check_gene_matrix <- function(gene_matrix) {
  ### function to check gene matrix is as expected for creation of GRM"
  # gene_matrix NxM+1  matrix with column sample_name, N is number of samples, M is number of genes
 if(!colnames(gene_matrix)[1] == "sample_name"){
   warning("first column in gene matrix is not sample_name this might cause unexpected behavior")
 }
  if(typeof(gene_matrix[,1]) != "character"){
    warning("first column in gene matrix is expected to be sample_name, if this is not the sample name and is a value of a gene, this will not be included in the computation.")
  }
  if(any(duplicated(gene_matrix[,1]))){
    stop("gene matrix should not have any duplicated sample names, sample names should be unique.")
  }
  if(any(is.na(gene_matrix[,1]))){
    stop("gene matrix should not have any NA sample names, sample names should be unique.")
  }
  if(!all(apply(gene_matrix[,-1],2,function(x) all(!is.na(as.numeric(x)))))){
    stop("all columns besides the first column should be numeric")
  }
  if((ncol(gene_matrix)) < 3 | ncol(gene_matrix) < nrow(gene_matrix)){
    warning("gene matrix does not have as many columns as rows, gene matrix should be a gene X sample matrix, if there are not many genes, the GRM will less precise.")
  }
  if(!is.data.frame(gene_matrix)){
    warning("gene_matrix is not a data frame, could cause unexpected behavior")
  }
}

check_inputs_tau <- function(glm.fit0,grm,species_id,tau0,phi0,maxiter,tol,verbose,write_log,log_file) {
  ## check inputs are as expected for tau test
  # glm.fit0 baseline glm from glm function
  # grm output from create_grm or NxN matrix
  # species_id id to denote species
  # tau0 numeric > 0 
  # phi0 numeric > 0 
  # maxiter maximum iterations
  # tol tolarnce for fitting usually <1
  # verbose whether to write updates as fitting occurs
  # write_log whether to write log file
  # log_file location to write log file
  if(!is.logical(verbose)){
    stop("verbose should be a logical")
  }
  if(write_log){
    if(!file_test("-x",log_file)){
      warning("log file is not writing to a file, check path, running without log file")
      write_log = FALSE
    }
  }
  if(!is.numeric(tol)){
    stop("tolarance should be numeric")
  }
  if(tol>.01){
    warning("tolarance is usually smaller than .01")
  }
  if(!is.numeric(maxiter)){
    stop("maxiter should be numeric")
  }
  if(maxiter<3){
    warning("maxiter is usually larger than 3")
  }
  if(maxiter>100){
    warning("large maxiter, might take a long time to converge")
  }
  if(!is.numeric(tau0)){
    stop("tau0 must be numeric")
  }
  if(!is.numeric(phi0)){
    stop("phi0 must be numeric")
  }
  if(tau0 <= 0 | phi0 <= 0){
    stop("tau0 and phi0 must be a postive number")
  }
  check_grm(grm,glm.fit0,verbose)
  
}

check_beta <- function(pop.struct.glmm, glm.fit0, grm, gene_df, SPA = FALSE) {
  ## check inputs for beta test for expected input
  # pop.struct.glmm output from tau test
  # glm.fit0 baseline glm from glm function
  # grm output from create_grm or NxN matrix
  # gene_df dataframe of gene_id, sample_name, and gene_value
  if(!is.logical(SPA)){
    stop("SPA should be a logical")
  }
  if(!(class(pop.struct.glmm)=="pop.struct.glmm")){
    stop("pop.struct.glmm should be an output of fit_tau_test")
  }
  if(all(pop.struct.glmm$y != glm.fit0$y)){
    warning("y for pop.struct.glmm does not match y for glm.fit0, if this is not intentional double check data.")
  }
  if(!is.data.frame(gene_df)){
    stop("gene_df is expected to be a data frame")
  }
  if(!any("gene_id" == colnames(gene_df))){
    stop("Expected column named gene_id, please use a gene_df with a column named gene_id")
  }
  
  if(!any("sample_name" == colnames(gene_df))){
    stop("Expected column named sample_name, please use a gene_df with a column named sample_name with sample names")
  }
  if(!any("gene_value" == colnames(gene_df))){
    stop("Expected column named gene_value, please use a gene_df with a column named gene_value with gene value")
  }
  sample_genes = unique(gene_df$gene_id)
  one_gene =  gene_df[which(gene_df$gene_id == sample_genes[1]),]
  if(!all(one_gene$sample_name == pop.struct.glmm$sample_names)){
    warning("sample names for gene_df do not match sample names for  pop.struct.glmm, if this is not intentional double check data for order of samples and number of samples.")
  }
  if(ncol(gene_df)>3){
    warning("expected stacked data frame of 3 columns, gene_id, sample_name, gene_value, data frame has more than 3 columns, check that data frame is stacked.")
  }
}

check_grm <- function(grm, glm.fit0, verbose) {
  ### function to check grm and glm.fit0 that they are the expected inputs
  # grm, NxN matrix
  # glm.fit0, output from glm, baseline glm
  if(!is.matrix(grm)){
    stop("grm is expected to be a matrix")
  }
  if(!is.numeric(grm) |  any(is.na(grm))){
    stop("grm is expected to have all numeric values")
  }
  if(nrow(grm) != ncol(grm)){
    stop("grm is expected to be a NxN matrix, with N as the number of samples")
  }
  if(!all(colnames(grm) == rownames(grm))){
    warning("column names do not match row names, grm should be a sample by sample matrix")
  }
  if(class(glm.fit0)[1] != "glm"){
    stop("Expected a glm object for glm.fit0, Example: glm_fit0=glm(`y ~ age  + 1`, data = exp_metadata, family = `binomial`)")
  }
  if(nrow(glm.fit0$data) != nrow(grm)){
    stop("number of samples from baseline glm does not match number of samples from grm, these must match")
  }
  if ("sample_name" %in% colnames(glm.fit0$data)) {
    if (verbose) {
      cat("\nchecking sample names of grm and glm fit match\n")
    }
    if (!all(glm.fit0$data$sample_name == rownames(grm))) {
      stop("\nERROR! the sample names for glm and grm do not match")
    } else {
      if (verbose) {
        cat("check complete")
      }
    }
  } else {
    warning("not running check on sample names ensure grm sample order and glm fit sample order match! Will only check if data from glm has column named sample_name")
  }
}
