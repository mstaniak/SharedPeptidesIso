#' Analytical gradient for the full model
#' @param model_design object of class IsoAPQModelDesign
#' @keywords internal
get_iso_fixed_gradient = function(model_design) {
    response = getIsoResponse(model_design)
    proteins_design = getIsoProteinDesign(model_design)
    num_proteins = getIsoNumProteins(model_design)
    fixed_design = getIsoFixedDesign(model_design)
    num_fixed_effects = getIsoNumFixed(model_design)

    gradient = get_gradient_fixed(proteins_design, fixed_design,
                                  response, num_fixed_effects,
                                  num_proteins)
    gradient
}


#' Nonlinear mixed effects model based on isotopic distribution
#' @param proteins_design design matrix for protein parameters
#' @param fixed_design design matrix for fixed effects
#' @param response response vector
#' @param num_fixed_effects number of fixed effects
#' @param num_proteins number of proteins
#' @param random_effects_counts counts of groups for each random effect
#' @param independent_unit vector of values of the independent unit
#' @return function of a single parameter that returns vector of gradient values
#' @keywords internal
get_gradient_fixed = function(proteins_design, fixed_design,
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

        # PROTEIN PARAMETERS GRADIENT
        grad_prot = (-2 / (sigma ^ 2) )* (
            t(proteins_design) %*% ((V_inv %*% r) /
                                        (proteins_design %*% exp(proteins)))) * exp(proteins)
        # SIGMA GRADIENT
        grad_sigma_logdet = 2 * sigma * sum(diag(V %*% solve((sigma ^ 2)*V)))
        grad_sigma_main = (-2 / (sigma ^ 3))* as.numeric(t(r) %*% V_inv %*% r)
        grad_sigma = grad_sigma_logdet + grad_sigma_main
        # FIXED GRADIENT
        fixed_gradient = -2 * (sigma ^ (-2)) * t(fixed_design) %*% V_inv %*% r
        c(sigma * grad_sigma, grad_prot, fixed_gradient)

    }
}
