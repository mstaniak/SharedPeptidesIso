input = data.table::as.data.table(expand.grid(
    iso_prob = 1,
    isoDistr = 1,
    sequence = letters[1:3],
    rep = paste0("r", 1:2),
    precursor_scan = 1:2,
    charge_in_scan = 1:2,
    iso_id = 1:3
))[order(rep, sequence, precursor_scan, charge_in_scan, iso_id)]
input[, precursor_scan := paste(sequence, rep, precursor_scan, sep = "_")]
input[, charge_in_scan := paste(precursor_scan, charge_in_scan, sep = "_")]
input
input$intensity = 15 + runif(nrow(input))
input[, y := intensity]
input[, protein_1 := as.integer(sequence %in% c("a", "b"))]
input[, protein_2 := as.integer(sequence %in% c("c", "b"))]
input

model_formula = intensity ~ log(..) + sequence + rep + (1|precursor_scan) + (1|charge_in_scan)

design_matrices_test = parse_formula(model_formula, model_data = input)
x = IsoAPQModelDesign(model_formula, design_matrices_test)

model_test = get_iso_full_nonlinear_loglikelihood(x)
params = c(0.1, rep(0.5, getIsoNumRandom(x)), rep(13, getIsoNumProteins(x)),
           rep(1, getIsoNumFixed(x)))
model_test(params)
