# install.packages("slendr")
# setup_env()

library(slendr)
init_env(quiet = TRUE)

anc <- population("ancestors", N = 5000, time = 6500000, remove = 649000)
chimp <- population("Chimpanzees", N = 5000, time = 6000000, parent = anc)
afr <- population("Africans", parent = anc, N = 10000, time = 650000)
nea <- population("Neanderthals", parent = anc, N = 2000, time = 650000)
eur <- population("Europeans", parent = afr, N = 3000, time = 60000)

gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 50000)

true_model <- compile_model(
  populations = list(chimp, anc, afr, nea, eur), gene_flow = gf,
  generation_time = 30
)

samples <- schedule_sampling(
  true_model, times = 0,
  list(eur, 4), list(afr, 5), list(nea, 2), list(chimp, 1)
)

# plot_model(true_model, proportions = TRUE,
#            order = c("Africans", "Europeans", "ancestors", "Neanderthals", "Chimpanzees"))

ts <- msprime(true_model, sequence_length = 200e6, recombination_rate = 1e-8, samples = samples, random_seed = 1234567) %>%
  ts_mutate(mutation_rate = 1e-8, random_seed = 42)

gt <- ts_genotypes(ts)

# turn genotypes into haploid by only taking one chromosome
gt <- dplyr::select(gt, -dplyr::contains("chr2"))

# remove the _chr1 suffix
colnames(gt) <- gsub("_chr1", "", colnames(gt))

# rename individuals to make them easier to work with
gt <- dplyr::rename(
  gt, African = Africans_1, Neanderthal = Neanderthals_1,
  another_Neanderthal = Neanderthals_2,
  Chimp = Chimpanzees_1,
  A = Africans_2, B = Africans_3, C = Africans_4, D = Africans_5,
  E = Europeans_1, F = Europeans_2, G = Europeans_3, H = Europeans_4
) %>% dplyr::select(pos, African, Neanderthal, Chimp, dplyr::everything())

# genotypes without "another_Neanderthal"
# saveRDS(gt[, -ncol(gt)], "genotypes_ex1.rds") # doesn't load in the ancient KU R
write.table(gt[, -ncol(gt)], file = "genotypes_ex1.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# genotypes with "another_Neanderthal"
# saveRDS(gt, "genotypes_ex2.rds") # doesn't load in the ancient KU R
write.table(gt, file = "genotypes_ex2.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
