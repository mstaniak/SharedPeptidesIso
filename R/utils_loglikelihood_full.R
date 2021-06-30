#' Negative log-likelihood for a full model based on IsoAPQModelDesign object
#' @param model_design object of class IsoAPQModelDesign
#' @param function of a single parameter
#' @keywords internal
get_iso_full_loglikelihood = function(model_design) {
    response = getIsoResponse(model_design)
    random_effects_counts = getIsoEffectsCounts(model_design)
    num_random_effects = getIsoNumRandom(model_design)
    proteins_design = getIsoProteinDesign(model_design)
    num_proteins = getIsoNumProteins(model_design)
    fixed_design = getIsoFixedDesign(model_design)
    num_fixed_effects = getIsoNumFixed(model_design)
    random_design = getIsoRandomDesign(model_design)
    independent_unit = getIsoIndependentUnit(model_design)

    neg_loglik = get_neglog_full(random_design, proteins_design, fixed_design,
                                 response, num_random_effects, num_fixed_effects,
                                 num_proteins, random_effects_counts,
                                 independent_unit)
    neg_loglik
}

#' Negative log-likelihood for a full model
#' @inheritParams get_gradient_full
#' @return function of a single parameter that returns negative log-likelihood
#' @keywords internal
get_neglog_full = function(random_design, proteins_design, fixed_design,
                           response, num_random_effects, num_fixed_effects,
                           num_proteins, random_effects_counts,
                           independent_unit) {
    function(params) {
        sigma = params[1]
        random_parameters = params[2:(num_random_effects + 1)]
        proteins = params[(num_random_effects + 2):(num_random_effects + 2 + num_proteins - 1)]
        fixed = params[(num_random_effects + 2 + num_proteins):length(params)]

        loglik_value = 0
        independent_values = unique(independent_unit)
        for (scan in independent_values) {
            id = independent_unit == scan
            loglik = get_neglog_full_unit(random_design[id, , drop = FALSE],
                                          proteins_design[id, , drop = FALSE],
                                          fixed_design[id, , drop = FALSE],
                                          response[id],
                                          num_random_effects,
                                          num_fixed_effects,
                                          num_proteins,
                                          sum(id),
                                          random_effects_counts)
            loglik_value = loglik_value + loglik(params)
        }
        loglik_value
    }
}


#' Negative log-likelihood within a single independent unit
#' @inheritParams get_gradient_full
#' @param num_rows number of observations
#' @return function a single parameter
#' @keywords internal
get_neglog_full_unit = function(random_design, proteins_design, fixed_design,
                                response, num_random_effects, num_fixed_effects,
                                num_proteins, num_rows, random_effects_counts) {
    function(params) {
        sigma = exp(params[1])
        random_parameters = exp(params[2:(num_random_effects + 1)])
        proteins = params[(num_random_effects + 2):(num_random_effects + 2 + num_proteins - 1)]
        fixed = params[(num_random_effects + 2 + num_proteins):length(params)]
        random_diagonal = rep(random_parameters, times = random_effects_counts)

        V = diag(1, num_rows) + random_design %*% diag(random_diagonal) %*% t(random_design)
        V_inv = solve(V)
        r = (response - log(proteins_design %*% exp(proteins)) - fixed_design %*% fixed)

        as.numeric(log(det((sigma ^ 2) * V)) + (1 / (sigma ^ 2)) * t(r) %*% V_inv %*% r)
    }
}
