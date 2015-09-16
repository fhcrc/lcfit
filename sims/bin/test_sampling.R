script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "rejection_sampling.R"))
source(file.path(script_dir, "log_sampling.R"))

set.seed(0)

lambda <- 0.1
N <- 1000

# small model
m.s <- list(c = 5, m = 8, r = 1, b = 0.5)

results.s <- test_model(m.s, lambda, N)
#inv_vs_rejR.s <- compare_results(results.s$inv, results.s$rejR)
inv_vs_rejC.s <- compare_results(results.s$inv, results.s$rejC)
#rejR_vs_rejC.s <- compare_results(results.s$rejR, results.s$rejC)

print(inv_vs_rejC.s$plot)

# medium model
m.m <- list(c = 1100, m = 800, r = 2, b = 0.5)

results.m <- test_model(m.m, lambda, N)
#inv_vs_rejR.m <- compare_results(results.m$inv, results.m$rejR)
inv_vs_rejC.m <- compare_results(results.m$inv, results.m$rejC)
#rejR_vs_rejC.m <- compare_results(results.m$rejR, results.m$rejC)

print(inv_vs_rejC.m$plot)

# long model
m.l <- list(c = 2000, m = 500, r = 2, b = 0.5)

results.l <- test_model(m.l, lambda, N)
#inv_vs_rejR.l <- compare_results(results.l$inv, results.l$rejR)
inv_vs_rejC.l <- compare_results(results.l$inv, results.l$rejC)
#rejR_vs_rejC.l <- compare_results(results.l$rejR, results.l$rejC)

print(inv_vs_rejC.l$plot)
