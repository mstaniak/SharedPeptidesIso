# FULL MODEL
get_neglog_full_unit = function(random_design, proteins_design, fixed_design,
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

        as.numeric(log(det((sigma ^ 2) * V)) + (1 / (sigma ^ 2)) * t(r) %*% V_inv %*% r)
    }
}

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
