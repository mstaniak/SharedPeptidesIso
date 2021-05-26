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

get_iso_full_gradient = function(x) {
    if (isIsoLinear(x)) {
        # get_iso_full_linear_loglikelihood(x, response)
        stop("Not implemented yet")
    } else {
        get_iso_full_nonlinear_gradient(x)
    }
}

get_iso_full_nonlinear_gradient = function(model_design) {
    response = getIsoResponse(model_design)
    random_effects_counts = getIsoEffectsCounts(model_design)
    num_random_effects = getIsoNumRandom(model_design)
    proteins_design = getIsoProteinDesign(model_design)
    num_proteins = getIsoNumProteins(model_design)
    fixed_design = getIsoFixedDesign(model_design)
    num_fixed_effects = getIsoNumFixed(model_design)
    random_design = getIsoRandomDesign(model_design)
    independent_unit = getIsoIndependentUnit(model_design)

    gradient = get_gradient_full(random_design, proteins_design, fixed_design,
                                 response, num_random_effects, num_fixed_effects,
                                 num_proteins, random_effects_counts,
                                 independent_unit)
    gradient
}

# get_iso_random_gradient = function(x) {
#     if (isIsoLinear(x)) {
#         get_iso_random_linear_model(x)
#     } else {
#         get_iso_random_nonlinear_model(x)
#     }
# }
#
# get_iso_fixed_loglikelihood = function(x) {
#     if (isIsoLinear(x)) {
#         get_iso_fixed_linear_loglikelihood(x)
#     } else {
#         get_iso_fixed_nonlinear_loglikelihood(x)
#     }
# }
#
# get_iso_protein_loglikelihood = function(x) {
#     if (isIsoLinear(x)) {
#         get_iso_protein_linear_loglikelihood(x)
#     } else {
#         get_iso_protein_nonlinear_loglikelihood(x)
#     }
# }
