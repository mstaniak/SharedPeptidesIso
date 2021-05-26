#' Fit a model based on model design
#' @param model_design object of class IsoAPQModelDesign
#' @return list
#' @importFrom stats optim
#' @keywords internal
fit_iso_full_model = function(model_design, ...) {
    grad_m = getIsoAnalyticalGradient(model_design)
    loglik_m = getIsoLogLikelihood(model_design)
    num_random_effects = getIsoNumRandom(model_design)
    num_fixed_effects = getIsoNumFixed(model_design)
    num_proteins = getIsoNumProteins(model_design)
    starting_values = c(0.1, rep(0.5, num_random_effects),
                        rep(12, num_proteins),
                        rep(1, num_fixed_effects))
    # TODO: getIsoStartingValues
    optimized_loglik = optim(starting_values, loglik_m, grad_m, method = "BFGS",
                             ...)
    optimized_loglik
}


#' Extract model information from optim() output and model design
#' @param fitted_model output of optim()
#' @param model_design object of class IsoAPQModelDesign
#' @return list
#' @keywords internal
get_full_model_information = function(fitted_model, model_design) {
    response = getIsoResponse(model_design)
    num_fixed = getIsoNumFixed(model_design)
    num_random = getIsoNumRandom(model_design)
    num_proteins = getIsoNumProteins(model_design)

    if (fitted_model[["convergence"]] == 0) {
        pars = fitted_model[["par"]]
        loglik = -fitted_model[["value"]]
        coef = c(exp(pars[1]), pars[-1])
        sigma = coef[1]
        random_effects_var = (sigma ^ 2) * coef[2:(num_random + 1)]
        names(random_effects_var) = names(getIsoEffectsCounts(model_design))
        protein_abundance = coef[(num_random + 2):(num_random + num_proteins + 1)]
        fixed_effects = coef[(num_random + num_proteins + 2):length(coef)]

        fitted = as.numeric(log(getIsoProteinDesign(model_design) %*% exp(protein_abundance)) +
                                getIsoFixedDesign(model_design) %*% fixed_effects)
        residuals = as.numeric(response - fitted)
        message("Computing hessian, please wait...")
        hessian = get_hessian(model_design, pars)
        conv_message = "model converged"
    } else {
        coef = c(NA, num_fixed + num_random + num_proteins)
        loglik = NA
        sigma = NA
        random_effects_var = rep(NA, num_random)
        protein_abundance = rep(NA, num_proteins)
        fixed_effects = rep(NA, num_fixed)

        fitted = rep(NA, getIsoNumObservations(model_design))
        residuals = rep(NA, getIsoNumObservations(model_design))
        hessian = matrix()
        conv_message = "model did not converge"
    }
    model_information = list(coef = coef, sigma = sigma,
                             random_effects_var = random_effects_var,
                             protein_abundance = protein_abundance,
                             fixed_effects = fixed_effects,
                             fitted = fitted, residuals = residuals,
                             hessian = hessian,
                             loglik = loglik, message = conv_message)
}


#' Get hessian based on model design and parameters
#' @param model_design object of class IsoAPQModelDesign
#' @param pars numeric
#' @return matrix
#' @keywords internal
get_hessian = function(model_design, pars) {
    grad_m = getIsoAnalyticalGradient(model_design)
    hessian_num = numDeriv::jacobian(grad_m, pars)
    row.names(hessian_num) = c("sigma", names(getIsoEffectsCounts(model_design)),
                               colnames(getIsoProteinDesign(model_design)),
                               colnames(getIsoFixedDesign(model_design)))
    colnames(hessian_num) = row.names(hessian_num)
    hessian_num
}
