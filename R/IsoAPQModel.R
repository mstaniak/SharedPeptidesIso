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
