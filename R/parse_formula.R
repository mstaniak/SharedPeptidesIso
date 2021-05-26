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

#' @keywords internal
get_candidate_formula = function(formula_chr) {
    if (grepl("log", formula_chr[3])) {
        formula_chr[3] = gsub("log(..)", "0", formula_chr[3], fixed = TRUE)
    } else {
        formula_chr[3] = gsub("..", "0", formula_chr[3], fixed = TRUE)
    }
    formula_chr = paste(formula_chr[c(2, 1, 3)], collapse = " ")
    as.formula(formula_chr)
}

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

#' @keywords internal
get_full_design = function(candidate_formula, model_data) {
    lmer_model = lme4::lmer(candidate_formula, data = model_data)
    fixed_design = lme4::getME(lmer_model, "X")
    random_design = lme4::getME(lmer_model, "Z")
    list(fixed = fixed_design,
         random = as.matrix(random_design))
}

#' @keywords internal
get_empty_design = function(model_data) {
    list(fixed = NULL,
         random = NULL)
}

#' @keywords internal
get_fixed_design = function(candidate_formula, model_data) {
    fixed_design = model.matrix(candidate_formula, data = model_dat)
    list(fixed = fixed_design,
         random = NULL)
}

#' @keywords internal
get_random_design = function(candidate_formula, model_data) {
    lmer_model = lme4::lmer(candidate_formula, data = model_data)
    random_design = lme4::getME(lmer_model, "Z")
    list(fixed = NULL,
         random = as.matrix(random_design))
}

#' @keywords internal
add_protein_design = function(design_matrices, model_data) {
    protein_design = get_protein_matrix(model_data)
    design_matrices[["protein"]] = protein_design
    design_matrices
}

#' @keywords internal
get_protein_matrix = function(model_data) {
    protein_cols = setdiff(colnames(model_data),
                           c("intensity", "y", "isoDistr", "iso_prob",
                             "precursor_scan", "charge_in_scan", "charge",
                             "rep", "sequence", "iso_id"))
    protein_formula = paste("intensity ~ 0 +", paste(protein_cols, collapse = " + "))
    model.matrix(as.formula(protein_formula), data = model_data)
}

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

#' @keywords internal
add_independent_unit = function(design_matrices, model_data) {
    # TODO: make it general?
    if (is.null(design_matrices[["random"]])) {
        unit = character(0)
    } else {
        unit = as.character(model_data[["precursor_scan"]])
    }
    design_matrices[["independent_unit"]] = unit
    design_matrices
}

#' @keywords internal
add_response = function(design_matrices, model_formula, model_data) {
    response_chr = as.character(model_formula)[2]
    response = eval(str2expression(response_chr), envir = model_data)
    design_matrices[["response"]] = response
    design_matrices
}
