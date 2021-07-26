#' Get design matrices and response from formula and input data
#' @param model_formula formula
#' @param model_data data.table
#' @return list
#' @keywords internal
parse_formula = function(model_formula, model_data) {
    formula_chr = as.character(model_formula)
    has_random = any(grepl("|", model_formula, fixed = TRUE))
    has_fixed = !all(grepl("|", setdiff(attr(terms(model_formula),
                                             "term.labels"), ".."), fixed = TRUE))
    candidate_formula = get_candidate_formula(formula_chr)
    design_matrices = get_design_matrices(candidate_formula, has_random, has_fixed, model_data)
    design_matrices
}


#' Get updated formula for use in lmer
#' @param formula_chr as.character(formula)
#' @keywords internal
get_candidate_formula = function(formula_chr) {
    if (grepl("log", formula_chr[3])) {
        formula_chr[3] = gsub("log(..)", "1", formula_chr[3], fixed = TRUE)
    } else {
        formula_chr[3] = gsub("..", "1", formula_chr[3], fixed = TRUE)
    }
    formula_chr = paste(formula_chr[c(2, 1, 3)], collapse = " ")
    as.formula(formula_chr)
}


#' Get design matrices
#' @param candidate_formula formula
#' @param has_random if TRUE, there are random effects
#' @param has_fixed, if TRUE, there are fixed effects
#' @param model_data data.table
#' @keywords internal
get_design_matrices = function(candidate_formula, has_random,
                               has_fixed, model_data) {
    if (has_random & has_fixed) {
        design_matrices = get_full_design(candidate_formula, model_data)
    } else if (!has_random & !has_fixed) {
        design_matrix = get_empty_design(model_data)
    } else if (!has_random & has_fixed) {
        design_matrices = get_fixed_design(candidate_formula, model_data)
    } else {
        design_matrices = get_random_design(candidate_formula, model_data)
    }
    design_matrices = add_protein_design(design_matrices, model_data)
    design_matrices = add_effect_counts(design_matrices, candidate_formula,
                                        model_data)
    design_matrices = add_independent_unit(design_matrices, model_data)
    design_matrices = add_response(design_matrices, candidate_formula,
                                   model_data)
    design_matrices
}


#' Design matrices for full model
#' @inheritParams get_design_matrices
#' @keywords internal
get_full_design = function(candidate_formula, model_data) {
    lmer_model = lme4::lmer(candidate_formula, data = model_data)
    fixed_design = lme4::getME(lmer_model, "X")
    random_design = lme4::getME(lmer_model, "Z")
    list(fixed = fixed_design,
         random = as.matrix(random_design))
}


#' Design matrices for model with just proteins
#' @inheritParams get_design_matrices
#' @keywords internal
get_empty_design = function(model_data) {
    list(fixed = NULL,
         random = NULL)
}

#' Get design matrices for model with no random effects
#' @inheritParams get_design_matrices
#' @keywords internal
get_fixed_design = function(candidate_formula, model_data) {
    fixed_design = model.matrix(candidate_formula, data = model_data)
    list(fixed = fixed_design,
         random = NULL)
}


#' Get design matrices for model with no fixed effects
#' @inheritParams get_design_matrices
#' @keywords internal
get_random_design = function(candidate_formula, model_data) {
    lmer_model = lme4::lmer(candidate_formula, data = model_data)
    random_design = lme4::getME(lmer_model, "Z")
    list(fixed = NULL,
         random = as.matrix(random_design))
}


#' Add design matrix for proteins
#' @param design_matrices list of design matrices
#' @inheritParams get_design_matrices
#' @keywords internal
add_protein_design = function(design_matrices, model_data) {
    protein_design = get_protein_matrix(model_data)
    design_matrices[["protein"]] = protein_design
    design_matrices
}


#' Create design matrix for proteins
#' @inheritParams get_design_matrices
#' @keywords internal
get_protein_matrix = function(model_data) {
    protein_cols = setdiff(colnames(model_data),
                           c("intensity", "y", "isoDistr", "iso_prob",
                             "precursor_scan", "charge_in_scan", "charge",
                             "rep", "sequence", "iso_id", "sequence_in_rep"))
    protein_formula = paste("intensity ~ 0 +",
                            paste(protein_cols, collapse = " + "))
    model.matrix(as.formula(protein_formula), data = model_data)
}


#' Add counts of levels of random effects to list of design matrices
#' @inheritParams add_protein_design
#' @inheritParams get_design_matrices
#' @import data.table
#' @keywords internal
add_effect_counts = function(design_matrices, candidate_formula, model_data) {
    if (!is.null(design_matrices[["random"]])) {
        column_names = colnames(design_matrices$random)
        all_effects = strsplit(as.character(candidate_formula)[3],
                               " + ", fixed = TRUE)[[1]]
        random_effects = all_effects[grepl("|", all_effects, fixed = TRUE)]
        random_effects = gsub("[\\(\\)1 \\|]", "", random_effects)
        random_effects


        random_effects_order = sapply(random_effects, function(x) {
            min(which(column_names %in% unique(model_data[[x]])))
        })
        random_effects = random_effects[order(random_effects_order)]
        random_effects_counts = model_data[, lapply(.SD, data.table::uniqueN),
                                           .SDcols = random_effects]
        random_effects_counts = unlist(random_effects_counts, use.names = TRUE)
    } else {
        random_effects_counts = integer()
    }
    design_matrices[["counts"]] = random_effects_counts
    design_matrices
}


#' Add vector with independent unit to list of design matrices
#' @inheritParams add_protein_design
#' @keywords internal
add_independent_unit = function(design_matrices, model_data) {
    # TODO: make it general?
    if (is.null(design_matrices[["random"]])) {
        unit = character(0)
    } else {
        random_effects_names = names(design_matrices[["counts"]])
        counts = sapply(model_data[, random_effects_names, with = FALSE],
                        data.table::uniqueN)
        independent_unit = random_effects_names[which.min(counts)]
        unit = as.character(model_data[[independent_unit]])
    }
    design_matrices[["independent_unit"]] = unit
    design_matrices
}


#' Add response to list of design matrices
#' @inheritParams add_protein_design
#' @param model_formula formula
#' @keywords internal
add_response = function(design_matrices, model_formula, model_data) {
    response_chr = as.character(model_formula)[2]
    response = eval(str2expression(response_chr), envir = model_data)
    design_matrices[["response"]] = response
    design_matrices
}
