# install.packages("slendr")
# setup_env()

library(slendr)
init_env(quiet = TRUE)

anc <- population("ancestor", N = 10000, time = 1000000, remove = 649000)
afr <- population("AFR", parent = anc, N = 10000, time = 650000)
nea <- population("NEA", parent = anc, N = 2000, time = 650000)
eur <- population("EUR", parent = afr, N = 5000, time = 80000)

gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 55000 - 30)

model <- compile_model(
  populations = list(anc, afr, nea, eur), gene_flow = gf,
  generation_time = 30
)

samples <- schedule_sampling(model, times = 0, list(eur, 100))

# plot_model(model, proportions = TRUE,
#            order = c("AFR", "EUR", "ancestor", "NEA"))

ts <- msprime(model, sequence_length = 50e6, recombination_rate = 1e-8, samples = samples, random_seed = 42)

tracts <- ts_tracts(ts, census = 55000, quiet = TRUE) %>%
  dplyr::select(-node_id, -pop, -source_pop, -source_pop_id) %>%
  dplyr::rename(individual = name)

n <- 50
tracts$bin <- cut(tracts$length, breaks = n, labels = FALSE)

write.table(tracts, file = "tracts.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)