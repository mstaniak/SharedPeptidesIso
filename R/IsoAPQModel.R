#' S4 class to represent a nonlinear mixed effects model based on isotopic distribution
#'
#' @slot model_design object of class IsoAPQModelDesign that describes
#' experimental design reflected in the model
#' @slot fitted_model_raw usually list (output of stats::optim), potentially
#' could be a different object with raw outputs of function used to fit a model
#' @slot coef numeric vector of model coefficient in the following order:
#' standard deviation of the error, parameters that describe random effects,
#' protein abundances, fixed effects
#' @slot sigma standard deviation of the random error
#' @slot random_effects_var numeric vector of estimated variances of random effects
#' @slot protein_abundance numeric vector of estimated protein abundances
#' @slot fixed_effects numeric vector of estimated fixed effects
#' @slot fitted numeric vector of fitted values
#' @slot residuals numeric vector of residuals
#' @slot hessian hessian of the negative log-likelihood
#' @slot loglik log-likelihood calculated at the solution
#' @slot message message indicating if the model converged
#'
#' @include IsoAPQModelDesign.R
#' @rdname isoapqmodel
#' @export
#'
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


#' Fit a nonlinear mixed effects model based on isotopic distribution
#'
#' @param model_formula formula that specifies the model. Special notation `..`
#' denotes all protein variables
#' @param model_data input data - data.table with columns: intensity,
#' iso_prob or isoDistr, sequence, rep, precursor_scan, charge_in_scan and
#' protein indicator columns.
#' @param additional parameters to the model fitting function (currently stats::optim)
#'
#' @return object of class IsoAPQModel
#' @rdname isoapqmodel
#' @export
#'
IsoAPQModel = function(model_formula, model_data, ...) {
    message("Preparing model design...")
    model_design = IsoAPQModelDesign(model_formula, model_data = model_data)
    response = getIsoResponse(model_design)
    message("Fitting nonlinear mixed effects model, please wait...")
    fitted_model = fitIsoModel(model_design, ...)
    if (inherits(fitted_model, "list")) {
        if (hasIsoFixedEffects(model_design) & hasIsoRandomEffects(model_design)) {
            model_information = get_full_model_information(fitted_model,
                                                           model_design)
        } else if (!hasIsoFixedEffects(model_design) & hasIsoRandomEffects(model_design)) {
            model_information = get_random_model_information(fitted_model,
                                                             model_design)
        } else if (hasIsoFixedEffects(model_design) & !hasIsoRandomEffects(model_design)) {
            model_information = get_fixed_model_information(fitted_model,
                                                            model_design)
        } else {
            stop("Not implemented yet")
            # get_iso_protein_loglikelihood(x)
        }
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


#' Extract model design from IsoAPQModel object
#'
#' @param x object of class IsoAPQModel
#'
#' @rdname apqmodel-methods
#' @return object of class IsoAPQModelDesign
#' @export
#'
setGeneric("getIsoModelDesign", function(x) standardGeneric("getIsoModelDesign"))


#' Extract raw fit from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return list or other type returned by the model fitting function such as optim
#' @export
#'
setGeneric("getIsoRawFit", function(x) standardGeneric("getIsoRawFit"))


#' Extract variances of random effects from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return numeric, vector of estimated variances of random effects
#' @export
#'
setGeneric("getIsoRandomVariances", function(x) standardGeneric("getIsoRandomVariances"))

#' Extract protein abundances from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return numeric, vector of values of estimated protein abundances
#' @export
#'
setGeneric("getIsoProteinAbundances", function(x) standardGeneric("getIsoProteinAbundances"))

#' Extract fixed effects from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return numeric, vector of values of estimated fixed effects
#' @export
#'
setGeneric("getIsoFixedEffects", function(x) standardGeneric("getIsoFixedEffects"))

#' Extract hessian of the negative log-likelihood from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return matrix - hessian of the negative log-likelihood evaluated at the solution
#' @export
#'
setGeneric("getIsoHessian", function(x) standardGeneric("getIsoHessian"))

#' Extract fixed effects from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return numeric- value of log-likelihood at the solution
#' @export
#'
setGeneric("getIsoLoglik", function(x) standardGeneric("getIsoLoglik"))

#' Extract fit status (message) from IsoAPQModel object
#'
#' @inheritParams getIsoModelDesign
#' @rdname apqgmodel-methods
#' @return character - message indicating if model converged
#' @export
#'
setGeneric("getIsoFitStatus", function(x) standardGeneric("getIsoFitStatus"))


setMethod("getIsoModelDesign", "IsoAPQModel", function(x) x@model_design)
setMethod("getIsoResponse", "IsoAPQModel", function(x) getIsoResponse(x@model_design))
setMethod("getIsoFormula", "IsoAPQModel", function(x) getIsoFormula(x@model_design))
setMethod("getIsoRawFit", "IsoAPQModel", function(x) x@fitted_model_raw)
setMethod("getIsoRandomVariances", "IsoAPQModel", function(x) x@random_effects_var)
setMethod("getIsoProteinAbundances", "IsoAPQModel", function(x) x@protein_abundance)
setMethod("getIsoFixedEffects", "IsoAPQModel", function(x) x@fixed_effects)
setMethod("getIsoHessian", "IsoAPQModel", function(x) x@hessian)
setMethod("getIsoLoglik", "IsoAPQModel", function(x) x@loglik)
setMethod("getIsoFitStatus", "IsoAPQModel", function(x) x@message)

#' @export
setMethod("coef", "IsoAPQModel", function(object, ...) object@coef)
#' @export
setMethod("fitted", "IsoAPQModel", function(object, ...) object@fitted)
#' @export
setMethod("residuals", "IsoAPQModel", function(object, ...) object@residuals)
#' @export
setMethod("sigma", "IsoAPQModel", function(object, ...) object@sigma)
