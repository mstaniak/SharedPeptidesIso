library(SharedPeptidesIso)
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

x = IsoAPQModelDesign(model_formula, model_data = input)
loglik = getIsoLogLikelihood(x)
agrad = getIsoAnalyticalGradient(x)

point = c(1, rep(0.5, getIsoNumRandom(x)), rep(12, getIsoNumProteins(x)),
          rep(2, getIsoNumFixed(x)))

analytical_1 = agrad(point)
numerical_1 = numDeriv::grad(loglik, point)

tinytest::expect_equal(analytical_1, numerical_1)

point_2 = c(0.1, rep(-0.5, getIsoNumRandom(x)), rep(12, getIsoNumProteins(x)),
          rep(1.5, getIsoNumFixed(x)))

analytical_2 = agrad(point_2)
numerical_2 = numDeriv::grad(loglik, point_2)

tinytest::expect_equal(analytical_2, numerical_2)

random_points = lapply(1:20, function(i) stats::runif(getIsoNumRandom(x) + getIsoNumProteins(x) + getIsoNumFixed(x) + 1))
for (point in random_points) {
    print(tinytest::expect_equal(agrad(point), numDeriv::grad(loglik, point), tolerance = 1e-6))
}

