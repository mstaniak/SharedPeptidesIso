

#' Nonlinear mixed effects model based on isotopic distribution
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @return function of a single parameter that returns vector of gradient values
#' @keywords internal
get_gradient_full = function(random_design, proteins_design, fixed_design,
                             response, num_random_effects, num_fixed_effects,
                             num_proteins, random_effects_counts,
                             independent_unit) {
    function(params) {
        sigma = params[1]
        random_parameters = params[2:(num_random_effects + 1)]
        proteins = params[(num_random_effects + 2):(num_random_effects + 2 + num_proteins - 1)]
        fixed = params[(num_random_effects + 2 + num_proteins):length(params)]

        i = 1
        independent_values = unique(independent_unit)
        gradient_values = vector("list", length(independent_values))
        for (scan in independent_values) {
            id = independent_unit == scan
            gradient = get_gradient_full_unit(random_design[id, , drop = FALSE],
                                              proteins_design[id, , drop = FALSE],
                                              fixed_design[id, , drop = FALSE],
                                              response[id],
                                              num_random_effects,
                                              num_fixed_effects,
                                              num_proteins,
                                              sum(id),
                                              random_effects_counts)
            gradient_values[[i]] = gradient(params)
            i = i + 1
        }
        sapply(1:length(params), function(i) {
            sum(sapply(gradient_values, function(x) x[i]))
        })
    }
}


get_gradient_full_unit = function(random_design, proteins_design, fixed_design,
                                  response, num_random_effects, num_fixed_effects,
                                  num_proteins, num_rows, random_effects_counts) {
    function(params) {
        sigma = exp(params[1])
        random_parameters = params[2:(num_random_effects + 1)]
        proteins = params[(num_random_effects + 2):(num_random_effects + 2 + num_proteins - 1)]
        fixed = params[(num_random_effects + 2 + num_proteins):length(params)]
        random_diagonal = rep(random_parameters, times = random_effects_counts)

        V = diag(1, num_rows) + random_design %*% diag(random_diagonal) %*% t(random_design)
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
        # THETA GRADIENT
        grad_theta = numeric(num_random_effects)
        for (i in seq_len(num_random_effects)) {
            random_effect_values_tmp = rep(0, num_random_effects)
            random_effect_values_tmp[i] = 1
            random_diagonal_tmp = rep(random_effect_values_tmp,
                                      times = random_effects_counts)

            grad_theta_logdet_1 = sum(diag(V_inv %*% (
                random_design %*% diag(random_diagonal_tmp) %*% t(random_design))))
            grad_theta_main_1 = -(1 / (sigma ^ 2)) *
                t(r) %*% (V_inv %*% (
                    random_design %*% diag(random_diagonal_tmp) %*% t(random_design)) %*% V_inv) %*% r
            grad_theta_1 = as.numeric(grad_theta_logdet_1 + grad_theta_main_1)
            grad_theta[i] = grad_theta_1
        }
        # FIXED GRADIENT
        fixed_gradient = -2 * (sigma ^ (-2)) * t(fixed_design) %*% V_inv %*% r
        c(sigma * grad_sigma, grad_theta, grad_prot, fixed_gradient)
    }
}
