setClassUnion("matrixOrNull", members = c("matrix", "NULL"))
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

setGeneric("getIsoProteinDesign", function(x) standardGeneric("getIsoProteinDesign"))
setMethod("getIsoProteinDesign", "IsoAPQModelDesign", function(x) x@protein)
setGeneric("getIsoFixedDesign", function(x) standardGeneric("getIsoFixedDesign"))
setMethod("getIsoFixedDesign", "IsoAPQModelDesign", function(x) x@fixed)
setGeneric("getIsoRandomDesign", function(x) standardGeneric("getIsoRandomDesign"))
setMethod("getIsoRandomDesign", "IsoAPQModelDesign", function(x) x@random)
setGeneric("getIsoResponse", function(x) standardGeneric("getIsoResponse"))
setMethod("getIsoResponse", "IsoAPQModelDesign", function(x) x@response)
setGeneric("getIsoEffectsCounts", function(x) standardGeneric("getIsoEffectsCounts"))
setMethod("getIsoEffectsCounts", "IsoAPQModelDesign", function(x) x@counts)
setGeneric("getIsoIndependentUnit", function(x) standardGeneric("getIsoIndependentUnit"))
setMethod("getIsoIndependentUnit", "IsoAPQModelDesign", function(x) x@independent_unit)
setGeneric("hasIsoFixedEffects", function(x) standardGeneric("hasIsoFixedEffects"))
setMethod("hasIsoFixedEffects", "IsoAPQModelDesign", function(x) !is.null(x@fixed))
setGeneric("hasIsoRandomEffects", function(x) standardGeneric("hasIsoRandomEffects"))
setMethod("hasIsoRandomEffects", "IsoAPQModelDesign", function(x) !is.null(x@random))
setGeneric("getIsoFormula", function(x) standardGeneric("getIsoFormula"))
setMethod("getIsoFormula", "IsoAPQModelDesign", function(x) x@formula)
setGeneric("getIsoNumRandom", function(x) standardGeneric("getIsoNumRandom"))
setMethod("getIsoNumRandom", "IsoAPQModelDesign", function(x) x@num_random_effects)
setGeneric("getIsoNumFixed", function(x) standardGeneric("getIsoNumFixed"))
setMethod("getIsoNumFixed", "IsoAPQModelDesign", function(x) x@num_fixed_effects)
setGeneric("getIsoNumProteins", function(x) standardGeneric("getIsoNumProteins"))
setMethod("getIsoNumProteins", "IsoAPQModelDesign", function(x) x@num_proteins)
setGeneric("getIsoNumObservations", function(x) standardGeneric("getIsoNumObservations"))
setMethod("getIsoNumObservations", "IsoAPQModelDesign", function(x) x@num_observations)
setGeneric("isIsoLinear", function(x) standardGeneric("isIsoLinear"))
setMethod("isIsoLinear", "IsoAPQModelDesign",
          function(x) !grepl("log(..)", as.character(getIsoFormula(x))[3],
                             fixed = TRUE))
# setGeneric("getIsoDesignSummary", function(x) standardGeneric("getIsoProteinDesign"))
# TODO: ^

setGeneric("getIsoAnalyticalGradient",
           function(x) standardGeneric("getIsoAnalyticalGradient"))
setMethod("getIsoAnalyticalGradient", "IsoAPQModelDesign",
          function(x) {
              if (hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_full_gradient(x)
              } else if (!hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  stop("Not implemented yet")
                  # get_iso_random_gradient(x, response)
              } else if (hasIsoFixedEffects(x) & !hasIsoRandomEffects(x)) {
                  stop("Not implemented yet")
                  # get_iso_fixed_gradient(x, response)
              } else {
                  stop("Not implemented yet")
                  # get_iso_protein_gradient(x, response)
              }
          })

setGeneric("getIsoLogLikelihood",
           function(x) standardGeneric("getIsoLogLikelihood"))
setMethod("getIsoLogLikelihood", "IsoAPQModelDesign",
          function(x) {
              if (hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_full_loglikelihood(x)
              } else if (!hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  stop("Not implemented yet")
                  # get_iso_random_loglikelihood(x, response)
              } else if (hasIsoFixedEffects(x) & !hasIsoRandomEffects(x)) {
                  stop("Not implemented yet")
                  # get_iso_fixed_loglikelihood(x, response)
              } else {
                  stop("Not implemented yet")
                  # get_iso_protein_loglikelihood(x, response)
              }
          })
