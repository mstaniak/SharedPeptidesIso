fit_iso_full_model = function(model_design) {
    grad_m = getIsoAnalyticalGradient(model_design)
    loglik_m = getIsoLogLikelihood(model_design)
    num_random_effects = getIsoNumRandom(model_design)
    num_fixed_effects = getIsoNumFixed(model_design)
    num_proteins = getIsoNumProteins(model_design)
    starting_values = c(0.1, rep(0.5, num_random_effects),
                        rep(12, num_proteins),
                        rep(1, num_fixed_effects))
    # TODO: getIsoStartingValues
    optimized_loglik = optim(starting_values, loglik_m, grad_m, method = "BFGS")
    optimized_loglik
}


setClass("IsoAPQModel",
         slots = c(model_design = "IsoAPQModelDesign",
                   fitted_model_raw = "ANY",
                   coef = "numeric",
                   sigma = "numeric",
                   random_effects_var = "numeric",
                   protein_abundance = "numeric",
                   fixed_effects = "numeric",
                   fitted = "numeric",
                   residuals = "numeric",
                   hessian = "matrix",
                   loglik = "numeric",
                   message = "character"))

IsoAPQModel = function(model_formula, model_data) {
    message("Preparing model design...")
    model_design = IsoAPQModelDesign(model_formula, model_data = model_data)
    response = getIsoResponse(model_design)
    message("Fitting nonlinear mixed effects model, please wait...")
    fitted_model = fitIsoModel(model_design)
    if (inherits(fitted_model, "list")) {
        model_information = get_full_model_information(fitted_model,
                                                       model_design)
    } else {
        stop("Not implemented yet")
    }

    new("IsoAPQModel",
        model_design = model_design,
        fitted_model_raw = fitted_model,
        coef = model_information[["coef"]],
        sigma = model_information[["sigma"]],
        random_effects_var = model_information[["random_effects_var"]],
        protein_abundance = model_information[["protein_abundance"]],
        fixed_effects = model_information[["fixed_effects"]],
        fitted = model_information[["fitted"]],
        residuals = model_information[["residuals"]],
        hessian = model_information[["hessian"]],
        loglik = model_information[["loglik"]],
        message = model_information[["message"]])
}


setGeneric("getIsoModelDesign", function(x) standardGeneric("getIsoModelDesign"))
setMethod("getIsoModelDesign", "IsoAPQModel", function(x) x@model_design)
setGeneric("getIsoResponse", function(x) standardGeneric("getIsoResponse"))
setMethod("getIsoResponse", "IsoAPQModel", function(x) getIsoResponse(x@model_design))
setGeneric("getIsoModelFormula", function(x) standardGeneric("getIsoModelFormula"))
setMethod("getIsoModelFormula", "IsoAPQModel", function(x) getIsoFormula(x@model_design))
setGeneric("getIsoRawFit", function(x) standardGeneric("getIsoRawFit"))
setMethod("getIsoRawFit", "IsoAPQModel", function(x) x@fitted_model_raw)
setMethod("coef", "IsoAPQModel", function(object, ...) object@coef)
setMethod("fitted", "IsoAPQModel", function(object, ...) object@fitted)
setMethod("residuals", "IsoAPQModel", function(object, ...) object@residuals)
setMethod("sigma", "IsoAPQModel", function(object, ...) object@sigma)

setGeneric("getIsoRandomVariances", function(x) standardGeneric("getIsoRandomVariances"))
setMethod("getIsoRandomVariances", "IsoAPQModel", function(x) x@random_effects_var)
setGeneric("getIsoProteinAbundances", function(x) standardGeneric("getIsoProteinAbundances"))
setMethod("getIsoProteinAbundances", "IsoAPQModel", function(x) x@protein_abundance)
setGeneric("getIsoFixedEffects", function(x) standardGeneric("getIsoFixedEffects"))
setMethod("getIsoFixedEffects", "IsoAPQModel", function(x) x@fixed_effects)
setGeneric("getIsoHessian", function(x) standardGeneric("getIsoHessian"))
setMethod("getIsoHessian", "IsoAPQModel", function(x) x@hessian)
setGeneric("getIsoLoglik", function(x) standardGeneric("getIsoLoglik"))
setMethod("getIsoLoglik", "IsoAPQModel", function(x) x@loglik)
setGeneric("getIsoFitStatus", function(x) standardGeneric("getIsoFitStatus"))
setMethod("getIsoFitStatus", "IsoAPQModel", function(x) x@message)


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


get_hessian = function(model_design, pars) {
    grad_m = getIsoAnalyticalGradient(model_design)
    hessian_num = numDeriv::jacobian(grad_m, pars)
    row.names(hessian_num) = c("sigma", names(getIsoEffectsCounts(model_design)),
                               colnames(getIsoProteinDesign(model_design)),
                               colnames(getIsoFixedDesign(model_design)))
    colnames(hessian_num) = row.names(hessian_num)
    hessian_num
}
