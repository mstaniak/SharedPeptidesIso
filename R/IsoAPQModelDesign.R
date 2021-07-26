setClassUnion("matrixOrNull", members = c("matrix", "NULL"))

#' Create a description of model design for nonlinear mixed effects model
#'
#' @slot protein protein design matrix - 0-1 matrix with a column for each protein in the dataset
#' @slot fixed design matrix for fixed effects
#' @slot random design matrix for random effects
#' @slot response response - numeric vector
#' @slot counts counts of levels of random effects
#' @slot formula model formula
#' @slot independent_unit character vector that specifiec independent blocks in the data
#' @slot num_proteins number of proteins
#' @slot num_random_effects number of random effects
#' @slot num_fixed_effects number of fixed effects
#' @slot num_observations number of observations
#'
#' @rdname isoapqmodeldesign
#' @export
#'
setClass("IsoAPQModelDesign",
         slots = c(protein = "matrixOrNull",
                   fixed = "matrixOrNull",
                   random = "matrixOrNull",
                   response = "numeric",
                   counts = "integer",
                   formula = "formula",
                   independent_unit = "character",
                   num_proteins = "integer",
                   num_random_effects = "integer",
                   num_fixed_effects = "integer",
                   num_observations = "integer")
)


#' Extract protein model design from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return matrix
#' @export
#'
setGeneric("getIsoProteinDesign", function(x) standardGeneric("getIsoProteinDesign"))

#' Extract fixed effects model design from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return matrix
#' @export
#'
setGeneric("getIsoFixedDesign", function(x) standardGeneric("getIsoFixedDesign"))

#' Extract random model design from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return matrix
#' @export
#'
setGeneric("getIsoRandomDesign", function(x) standardGeneric("getIsoRandomDesign"))

#' Extract response vector from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return numeric vector
#' @export
#'
setGeneric("getIsoResponse", function(x) standardGeneric("getIsoResponse"))

#' Extract counts for levels of random effects from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return numeric vector
#' @export
#'
setGeneric("getIsoEffectsCounts", function(x) standardGeneric("getIsoEffectsCounts"))

#' Extract vector that describes independent blocks from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return character vecotr
#' @export
#'
setGeneric("getIsoIndependentUnit", function(x) standardGeneric("getIsoIndependentUnit"))

#' Check if model described by IsoAPQModelDesign has fixed effects
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return logical
#' @export
#'
setGeneric("hasIsoFixedEffects", function(x) standardGeneric("hasIsoFixedEffects"))

#' Check if model described by IsoAPQModelDesign has random effects
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return logical
#' @export
#'
setGeneric("hasIsoRandomEffects", function(x) standardGeneric("hasIsoRandomEffects"))

#' Extract model formula from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return formula
#' @export
#'
setGeneric("getIsoFormula", function(x) standardGeneric("getIsoFormula"))

#' Extract number of random effects from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return integer
#' @export
#'
setGeneric("getIsoNumRandom", function(x) standardGeneric("getIsoNumRandom"))

#' Extract number of fixed effects from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return integer
#' @export
#'
setGeneric("getIsoNumFixed", function(x) standardGeneric("getIsoNumFixed"))

#' Extract number of proteins from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return integer
#' @export
#'
setGeneric("getIsoNumProteins", function(x) standardGeneric("getIsoNumProteins"))

#' Extract number of observations from IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return integer
#' @export
#'
setGeneric("getIsoNumObservations", function(x) standardGeneric("getIsoNumObservations"))

#' Check if model described by IsoAPQModelDesign is linear
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return logical
#' @export
#'
setGeneric("isIsoLinear", function(x) standardGeneric("isIsoLinear"))

#' Get analytical gradient for a IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return function of a single parameter
#' @export
#'
setGeneric("getIsoAnalyticalGradient",
           function(x) standardGeneric("getIsoAnalyticalGradient"))

#' Get negative log-likelihood for a IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return function of a single parameter
#' @export
#'
setGeneric("getIsoLogLikelihood",
           function(x) standardGeneric("getIsoLogLikelihood"))

#' Fit a nonlinear mixed effects model based on a IsoAPQModelDesign object
#'
#' @param x object of class IsoAPQModelDesign
#'
#' @rdname apqmodeldesign-methods
#' @return object of class IsoAPQModel
#' @export
#'
setGeneric("fitIsoModel", function(x, ...) standardGeneric("fitIsoModel"))


#' Create design for nonlinear mixed effects model based on isotopic distribution
#'
#' @param inheritParams IsoAPQModel
#' @param model_design_list optional list of elements of the model design
#'
#' @rdname isoapqmodeldesign
#' @return object of class IsoAPQModelDesign
#' @export
#'
IsoAPQModelDesign = function(model_formula, model_design_list = NULL,
                             model_data = NULL) {
    if (is.null(model_design_list)) {
        model_design_list = parse_formula(model_formula, model_data)
    }
    new("IsoAPQModelDesign",
        protein = model_design_list[["protein"]],
        fixed = model_design_list[["fixed"]],
        random = model_design_list[["random"]],
        response = model_design_list[["response"]],
        counts = model_design_list[["counts"]],
        formula = model_formula,
        independent_unit = model_design_list[["independent_unit"]],
        num_proteins = ncol(model_design_list[["protein"]]),
        num_random_effects = length(model_design_list[["counts"]]),
        num_fixed_effects = ncol(model_design_list[["fixed"]]),
        num_observations = nrow(model_design_list[["protein"]]))
}


setMethod("getIsoProteinDesign", "IsoAPQModelDesign", function(x) x@protein)
setMethod("getIsoFixedDesign", "IsoAPQModelDesign", function(x) x@fixed)
setMethod("getIsoRandomDesign", "IsoAPQModelDesign", function(x) x@random)
setMethod("getIsoResponse", "IsoAPQModelDesign", function(x) x@response)
setMethod("getIsoEffectsCounts", "IsoAPQModelDesign", function(x) x@counts)
setMethod("getIsoIndependentUnit", "IsoAPQModelDesign", function(x) x@independent_unit)
setMethod("hasIsoFixedEffects", "IsoAPQModelDesign", function(x) !is.null(x@fixed) & !(ncol(x@fixed) == 0))
setMethod("hasIsoRandomEffects", "IsoAPQModelDesign", function(x) {
    if (is.null(x@random)) {
        FALSE
    } else {
        ncol(x@random) != 0
    }
})
setMethod("getIsoFormula", "IsoAPQModelDesign", function(x) x@formula)
setMethod("getIsoNumRandom", "IsoAPQModelDesign", function(x) x@num_random_effects)
setMethod("getIsoNumFixed", "IsoAPQModelDesign", function(x) x@num_fixed_effects)
setMethod("getIsoNumProteins", "IsoAPQModelDesign", function(x) x@num_proteins)
setMethod("getIsoNumObservations", "IsoAPQModelDesign", function(x) x@num_observations)
setMethod("isIsoLinear", "IsoAPQModelDesign",
          function(x) !grepl("log(..)", as.character(getIsoFormula(x))[3],
                             fixed = TRUE))
# setGeneric("getIsoDesignSummary", function(x) standardGeneric("getIsoProteinDesign"))
# TODO: ^

setMethod("getIsoAnalyticalGradient", "IsoAPQModelDesign",
          function(x) {
              if (hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_full_gradient(x)
              } else if (!hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_random_gradient(x)
              } else if (hasIsoFixedEffects(x) & !hasIsoRandomEffects(x)) {
                  get_iso_fixed_gradient(x)
              } else {
                  get_iso_protein_gradient(x)
              }
          })

setMethod("getIsoLogLikelihood", "IsoAPQModelDesign",
          function(x) {
              if (hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_full_loglikelihood(x)
              } else if (!hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_random_loglikelihood(x)
              } else if (hasIsoFixedEffects(x) & !hasIsoRandomEffects(x)) {
                  get_iso_fixed_loglikelihood(x)
              } else {
                  get_iso_protein_loglikelihood(x)
              }
          })


setMethod("fitIsoModel", "IsoAPQModelDesign", function(x, ...) {
    if (hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
        fit_iso_full_model(x, ...)
    } else if (!hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
        fit_iso_random_model(x, ...)
    } else if (hasIsoFixedEffects(x) & !hasIsoRandomEffects(x)) {
        fit_iso_fixed_model(x, ...)
    } else {
        stop("Not implemented yet")
    }
})
