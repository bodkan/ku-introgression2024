# install.packages("slendr")
# setup_env()

library(slendr)
init_env(quiet = TRUE)

anc <- population("ancestors", N = 10000, time = 6500000, remove = 649000)
chimp <- population("Chimpanzees", N = 20000, time = 6000000, parent = anc)
afr <- population("Africans", parent = anc, N = 10000, time = 650000)
nea <- population("Neanderthals", parent = anc, N = 2000, time = 650000)
eur <- population("Europeans", parent = afr, N = 5000, time = 60000)

gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 45000)

true_model <- compile_model(
  populations = list(chimp, anc, afr, nea, eur), gene_flow = gf,
  generation_time = 30
)

samples <- schedule_sampling(
  true_model, times = 0,
  list(eur, 4), list(afr, 5), list(nea, 1), list(chimp, 1)
)

# plot_model(true_model, proportions = TRUE,
#            order = c("Africans", "Europeans", "ancestors", "Neanderthals", "Chimpanzees"))

ts <- msprime(true_model, sequence_length = 100e6, recombination_rate = 1e-8, samples = samples, random_seed = 42) %>%
  ts_mutate(mutation_rate = 1e-8, random_seed = 42)

gt <- ts_genotypes(ts)

# turn genotypes into haploid by only taking one chromosome
gt <- dplyr::select(gt, -dplyr::contains("chr2"))

# remove the _chr1 suffix
colnames(gt) <- gsub("_chr1", "", colnames(gt))

# rename individuals to make them easier to work with
gt <- dplyr::rename(
  gt, African = Africans_1, Neanderthal = Neanderthals_1, Chimp = Chimpanzees_1,
  A = Africans_2, B = Africans_3, C = Africans_4, D = Africans_5,
  E = Europeans_1, F = Europeans_2, G = Europeans_3, H = Europeans_4
) %>% dplyr::select(pos, African, Neanderthal, Chimp, dplyr::everything())

# subsample genotypes to a more manageable size
set.seed(12345)
gt <- dplyr::sample_n(gt, 300000) %>% dplyr::arrange(pos)

invariant <- dplyr::select(gt, -pos) %>% rowMeans() %>% {. == 0 | . == 1}
fixed <- dplyr::select(gt, -pos, -Chimp) %>%
  rowMeans() %>%
  {. == 0 | . == 1} %>%
  {(. - gt$Chimp) != 0}
gt <- gt[!invariant & !fixed, ]

set.seed(19)
gt <- dplyr::sample_n(gt, 100000) %>% dplyr::arrange(pos)

saveRDS(gt, "genotypes.rds")
