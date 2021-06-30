#' Negative log-likelihood for a full model based on IsoAPQModelDesign object
#' @param model_design object of class IsoAPQModelDesign
#' @param function of a single parameter
#' @keywords internal
get_iso_fixed_loglikelihood = function(model_design) {
    response = getIsoResponse(model_design)
    proteins_design = getIsoProteinDesign(model_design)
    num_proteins = getIsoNumProteins(model_design)
    fixed_design = getIsoFixedDesign(model_design)
    num_fixed_effects = getIsoNumFixed(model_design)

    neg_loglik = get_neglog_fixed(proteins_design, fixed_design,
                                  response, num_fixed_effects,
                                  num_proteins)
    neg_loglik
}

#' Negative log-likelihood for a full model
#' @inheritParams get_gradient_full
#' @return function of a single parameter that returns negative log-likelihood
#' @keywords internal
get_neglog_fixed = function(proteins_design, fixed_design,
                            response, num_fixed_effects,
                            num_proteins) {
    num_rows = nrow(fixed_design)
    function(params) {
        sigma = exp(params[1])
        proteins = params[2:(num_proteins + 1)]
        fixed = params[(2 + num_proteins):length(params)]

        V = diag(1, num_rows)
        V_inv = solve(V)
        r = (response - log(proteins_design %*% exp(proteins)) - fixed_design %*% fixed)

        as.numeric(log(det((sigma ^ 2) * V)) + (1 / (sigma ^ 2)) * t(r) %*% V_inv %*% r)
    }
}
