setGeneric("getIsoLogLikelihood",
           function(x, response) standardGeneric("getIsoLogLikelihood"),
           signature = "x")
setMethod("getIsoLogLikelihood", "IsoAPQModelDesign",
          function(x, response) {
              if (hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_full_loglikelihood(x, response)
              } else if (!hasIsoFixedEffects(x) & hasIsoRandomEffects(x)) {
                  get_iso_random_loglikelihood(x, response)
              } else if (hasIsoFixedEffects(x) & !hasIsoRandomEffects(x)) {
                  get_iso_fixed_loglikelihood(x, response)
              } else {
                  get_iso_protein_loglikelihood(x, response)
              }
          })

get_iso_full_loglikehood = function(x, response) {
    if (isIsoLinear(x)) {
        get_iso_full_linear_loglikelihood(x, response)
    } else {
        get_iso_full_nonlinear_loglikelihood(x, response)
    }
}

get_iso_full_nonlinear_loglikelihood = function(model_design, response) {
    random_effects_counts = getIsoEffectsCounts(model_design)
    num_random_effects = length(random_effects_counts)
    proteins_design = getIsoProteinDesign(model_design)
    num_proteins = ncol(proteins_design)
    fixed_design = getIsoFixedDesign(model_design)
    num_fixed_effects = ncol(fixed_design)
    random_design = getIsoRandomDesign(model_design)
    independent_unit = getIsoIndependentUnit(model_design)

    neg_loglik = get_neglog_full(random_design, proteins_design, fixed_design,
                                 response, num_random_effects, num_fixed_effects,
                                 num_proteins, random_effects_counts,
                                 independent_unit)
    neg_loglik
}

get_iso_random_loglikelihood = function(x) {
    if (isIsoLinear(x)) {
        get_iso_random_linear_model(x)
    } else {
        get_iso_random_nonlinear_model(x)
    }
}

get_iso_fixed_loglikelihood = function(x) {
    if (isIsoLinear(x)) {
        get_iso_fixed_linear_loglikelihood(x)
    } else {
        get_iso_fixed_nonlinear_loglikelihood(x)
    }
}

get_iso_protein_loglikelihood = function(x) {
    if (isIsoLinear(x)) {
        get_iso_protein_linear_loglikelihood(x)
    } else {
        get_iso_protein_nonlinear_loglikelihood(x)
    }
}
