# Archaic introgression lecture exercises (February 2024)

## Exercise 1 -- D statistics / introgression test

You sequenced genomes of two some Africans and Eurasians (four genomes from each population). Unfortunately, there's been a mix up in the lab and you don't know which one is which! You only know that they are labeled A, B, C, ..., H. What a disaster!

Fortunately, you also have three other genomes which you know are from an African individual, a Neanderthal, and a Chimpanzee outgroup. This means you are able to compute a D statistic which will test for evidence of Neanderthal introgression in a given sample.

Can you save the day and determine which of the A, B, C, ..., H genomes are African and which are Eurasian based on the following D statistic test?

$$
D(\textrm{African}, X; \textrm{Neanderthal}, \textrm{Chimp}).
$$

(Recall that only Eurasian genomes are expected to have appreciable amounts of Neanderthal ancestry but Africans don't.)

First load the genotype table into R:

```         
gt <- readRDS(url("https://github.com/bodkan/ku-introgression2024/raw/main/genotypes.rds"))
```

The `gt` data set is a plain R data frame where each individual's column contains the genotype of that individual's chromosome (`0` - ancestral allele, `1` - derived allele).

Familiarize yourself with the data by running this R command which shows information from only the first few genotypes:

```         
head(gt)
```

You can extract *all the genotypes* of a given individual by using the `$` or `[[` subsetting operators of R data frames like this:

```         
gt$African

gt[["African"]]
```

A useful trick for comparing two chromosomes in their entirety is to rely in the fact that R can perform *vectorized operations* (operations performed on multiple elements of a vector at once). For instance, if this gives us the genotypes of an African and Neanderthal chromosome:

```         
gt[["African"]]

gt[["Neanderthal"]]
```

Then we can use this to sum up how many positions of those two chromosomes carry the same allele:

```         
gt[["African"]] == gt[["Neanderthal"]] # this gives us TRUE/FALSE values 

# sum() treats TRUE as 1 and FALSE as 0, so we can sum everything up!
sum(gt[["African"]] == gt[["Neanderthal"]])
```

On the other hand, this would count how many alleles are _different_ between our African and Chimpanzee chromosome:

```
sum(gt[["African"]] != gt[["Neanderthal"]]) # note the != instead of ==
```

Inside the `sum()` function we can compose multiple logical conditions to create more complex comparison operations using the `&` operator (AND operation in mathematical logic).

Armed with this knowledge, we can compute the BABA and ABBA counts using this bit of R code:

```
X = "A"  # or "B", or "C", ..., or "H"

abba <- sum(
  (gt[["African"]] == gt[["Chimp"]]) &         # filters for   A**A
  (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for   A*B*
  (gt[[X]]         == gt[["Neanderthal"]])     # filters for   *BB*
)                                              # together then ABBA

baba <- sum(
  (gt[["African"]] != gt[["Chimp"]]) &         # filters for   B**A
  (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for   B*B*
  (gt[[X]]         == gt[["Chimp"]])           # filters for   *A*A
)                                              # together then BABA

baba - abba
```

**You know that if `X` is a African, you expect to see roughly the same count of `BABA` and `ABBA` site patterns, so the difference should "be about zero". Compute this for all of your mixed up chromosomes A, B, C, ..., H and note down the `baba - abba` values you got for each -- which ones are most likely African and which ones are Eurasian?**

**[If you are more familiar with R, compute the counts automatically in a loop of some kind and make a figure.]**

**What does it mean for a test statistic to "be about zero"? What are we missing to truly use this as a statistical significance test?**

_Solution:_

```
X <- c("A", "B", "C", "D", "E", "F", "G", "H")

D_values <- sapply(X, function(x) {
  abba <- sum(
    (gt[["African"]] == gt[["Chimp"]]) &         # filters for   A**A
    (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for   A*B*
    (gt[[x]]         == gt[["Neanderthal"]])     # filters for   *BB*
  )                                              # together then ABBA
  
  baba <- sum(
    (gt[["African"]] != gt[["Chimp"]]) &         # filters for   B**A
    (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for   B*B*
    (gt[[x]]         == gt[["Chimp"]])           # filters for   *A*A
  )                                              # together then BABA
  
  (baba - abba) / (baba + abba)
})

plot(D_values, xaxt = "n", xlab = "test sample", ylab = " D(African, X; Neanderthal, Chimp)")
abline(h = 0, lty = 2, col = "red")
axis(side = 1, at = seq_along(X), labels = X)

# In this simple example we're missing confidence intervals -- those would allow
# us to do a proper statistical test to determine for which samples we really cannot
# reject a null hypothesis of no gene flow from Neanderthals. Real-world software
# such as ADMIXTOOLS (https://github.com/DReichLab/AdmixTools) computes confidence
# intervals using a so-called bootstrap procedure across windows along a genome
# (https://en.wikipedia.org/wiki/Bootstrapping_(statistics)).
```

