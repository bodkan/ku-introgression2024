# Archaic introgression lecture exercises (February 2024)

## Exercise 1 -- D statistics / introgression test

### Introduction

*You sequenced the genomes of four Africans and four Eurasians and got genotypes from a single chromosome from each of them (so you have genotypes of four African and four Eurasian chromosomes). Unfortunately, there's been a mix up in the lab and you don't know which one is which! You only know that they are labeled A, B, C, ..., H. What a disaster!*

*Fortunately, you also have genotypes from three other individuals whose identity you know for certain: an African, a Neanderthal, and a Chimpanzee. This means you are able to compute a D statistic which will test for evidence of Neanderthal introgression in a given individual.*

*Can you save the day and determine which of the A, B, C, ..., H samples are African and which are Eurasian based on the following D statistic test?*

$$
D(\textrm{African}, X; \textrm{Neanderthal}, \textrm{Chimp}).
$$

*(Recall that only Eurasians are expected to have appreciable amounts of Neanderthal ancestry but Africans don't.)*

### Moving over to R

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
# -- this gives us the number of positions at which an African carries the same
#    allele as the Neanderthal
sum(gt[["African"]] == gt[["Neanderthal"]])
```

On the other hand, this would count how many alleles are *different* between our African and Chimpanzee chromosome:

```         
sum(gt[["African"]] != gt[["Neanderthal"]]) # note the != instead of ==
```

Inside the `sum()` function we can compose multiple logical conditions to create more complex comparison operations using the `&` operator (AND operation in mathematical logic).

Armed with this knowledge, we can compute the BABA and ABBA counts using this bit of R code:

```         
X = "A"  # or "B", or "C", ..., or "H"

abba <- sum(
  (gt[["African"]] == gt[["Chimp"]]) &         # filters for A**A sites
  (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for A*B* sites
  (gt[[X]]         == gt[["Neanderthal"]])     # filters for *BB* sites
)                                              # together then ABBA

baba <- sum(
  (gt[["African"]] != gt[["Chimp"]]) &         # filters for B**A sites
  (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for B*B* sites
  (gt[[X]]         == gt[["Chimp"]])           # filters for *A*A sites
)                                              # together then BABA

baba - abba
```

**You know that if `X` is a African, you expect to see roughly the same count of `BABA` and `ABBA` site patterns, so the difference should "be about zero". Compute this for all of your mixed up chromosomes A, B, C, ..., H and note down the `baba - abba` values you got for each -- which ones are most likely African and which ones are Eurasian?**

**[If you are more familiar with R, compute the counts automatically in a loop of some kind and make a figure.]**

**What does it mean for a test statistic to "be about zero"? What are we missing to truly use this as a statistical significance test?**

*Solution:*

```         
X <- c("A", "B", "C", "D", "E", "F", "G", "H")

D_values <- sapply(X, function(x) {
  abba <- sum(
    (gt[["African"]] == gt[["Chimp"]]) &         # filters for A**A sites
    (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for A*B* sites
    (gt[[x]]         == gt[["Neanderthal"]])     # filters for *BB* sites
  )                                              # together then ABBA
  
  baba <- sum(
    (gt[["African"]] != gt[["Chimp"]]) &         # filters for B**A sites
    (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for B*B* sites
    (gt[[x]]         == gt[["Chimp"]])           # filters for *A*A sites
  )                                              # together then BABA
  
  (baba - abba) / (baba + abba)
})

plot(D_values, xaxt = "n", xlab = "test sample", ylab = " D(African, X; Neanderthal, Chimp)")
abline(h = 0, lty = 2, col = "red")
axis(side = 1, at = seq_along(X), labels = X)
```

We can see that the samples A-D are consistent with a D statistic value of about 0, meaning that the BABA and ABBA counts were about the same. This is what we would expect for African samples who are not expected to be closer to a Neanderthal genome than another African.

On the other hand, samples E-H show a much more negative value of the D statistic, which is consistent with an access of ABBA site patterns -- which arise with an increased sharing of derived alleles between the sample X and a Neanderthal genome.

**Important:** In this simple example we're missing confidence intervals -- those would allow
us to do a proper statistical test to determine for which samples we really cannot
reject a null hypothesis of no gene flow from Neanderthals. This way we can
avoid the vague and statistically unsatisfying talk about a value being "almost zero",
and some other value being "much more negative".

Real-world software such as ADMIXTOOLS (https://github.com/DReichLab/AdmixTools) computes confidence
intervals using a so-called bootstrap procedure across windows along a genome
(https://en.wikipedia.org/wiki/Bootstrapping_(statistics)).
------------------------------------------------------------------------

If you want to take a closer look at how the genotype data was prepared (it was simulated!), you can see the complete code [here](generate_genotypes.R).

## Exercise 2 -- dating Neanderthal admixture

### Introduction

You managed to sequence a bunch more Eurasian genomes (100 of diploid individuals) and are interested in double-checking the time of Neanderthal introgression into the ancestors of all Eurasians that's in the literature. To be able to do this, you ran an inference software which gives you the exact coordinates of Neanderthal DNA tracts present in every Eurasian genome that you sequenced. This of course means that you also know the lengths of each of those tracts.

You the distribution of the Neanderthal tract lengths in your population to estimate the time of Neanderthal introgression!

### Moving over to R

First load the table with coordinates of all Neanderthal tracts into R:

```         
tracts <- readRDS(url("https://github.com/bodkan/ku-introgression2024/raw/main/tracts.rds"))
```

The `gt` data set is a plain R data frame where each individual's column contains the genotype of that individual's chromosome (`0` - ancestral allele, `1` - derived allele).

Familiarize yourself with the data by running this R command which shows information from only the first few genotypes:

```         
head(tracts)
```

For how many individuals do we have information about Neanderthal DNA tracts that they carry?

```         
length(unique(tracts$individual))
```

It looks like the inference software (or a helpful bioinformatician) binned each tract according to its length. What does the distribution of Neanderthal tract lengths looks like in your data? Knowing that recombination has acted on the introgressed Neanderthal DNA over time, each generation, suggests that the distribution should look exponential -- do you see this in the data? To answer this, plot the proportion of tracts in each bin.

```         
# get the bin numbers
bins <- sort(unique(tracts$bin))

# count the tracts in each bin and compute the proportion of tracts in each bin
counts <- as.integer(table(tracts$bin))
props <- counts / sum(counts)

plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion")
```

The distribution does, indeed, look quite exponential. Let's try to use this to date the Neanderthal introgression using information encoded in the distribution of tract lengths!

As we know, over time since admixture, recombination breaks up longer haplotypes into shorter one, regularly almost like a clock. And it turns out that the distribution of tract lengths after time $t$ follows exponential decay, leading to the distribution of tract lengths $x$ to have the following form:

$$
f(\textrm{x}) \sim \exp(-\lambda x) = \exp(-r t x)
$$

Where the $\lambda$ parameter determines the _rate_ of exponential decay and, under some simplifying assumption can be computed as the product of the recombination rate (traditionally in humans with value of about $1e-8$) and $t$ which is the time since admixture -- **the latter is our unknown we're trying to compute here**.

It also turns out that the [expected value](https://en.wikipedia.org/wiki/Exponential_distribution#Mean,_variance,_moments,_and_median) of this exponential distribution (which can be computed simply as $1 / \lambda$) gives us the expected tract length after time $t$. Starting from our data, we can compute this expected length simply by computing the average introgressed tract length in our data like this:

```         
L <- mean(tracts$length)
L # length in units of base pairs
```

Taking all the math together and doing a little algebra, we can write this:

$$
\textrm{average tract length} L = \frac{1}{\lambda} = \frac{1}{rt}
$$

But because we know $L$ (average tract length) and $r$ (recombination rate), this means we can get an estimate of the time since the admixture like this by simply rearranging the equation to get $t$:

$$
t = \frac{1}{rL}
$$

Note that this time estimate will be in units of generations, so we'll have to multiply this quantity by the generation time (roughly 30 years for humans) to get time in years before the present.

We can use this equation to compute the estimate of the admixture time now:

```         
r <- 1e-8 # crossovers per bp per generation

t <- 1 / (r * L)
t

# convert the time of introgression to years before present assuming
# generation time of 30 years
t * 30
```

You should get a value will be quite close to ~55 thousand years ago, an estimate which is often found in the literature as the time when Neanderthals and anatomically modern humans interbred!

As a last sanity check, if we use this time to compute the rate of exponential decay $\lambda$, we should get a nice fit of the theoretical exponential decay curve over the empirical counts of tract lengths in each bin. As a reminder, this is the decay:

```
plot(bins, props)
```

Let's try if we can plot the theoretical exponential decay from the estimated time of admixture.

First, because the exponential plot above shows the decay of the size of introgressed tracts in bins, the $\lambda$ parameter determines the decay with each sucessive bin. However, our recombination rate $r$ which features in the equation to compute $\lambda$ above describes recombination in units of base pairs, not bins. First compute the average increase in tract length as we move from bin to bin:

```
# compute the average length in each bin
average_bins <- aggregate(length ~ bin, data = tracts, FUN = mean)

# compute the difference between bins to get "average bin step size"
bin_step <- mean(diff(average_bins$length))

bin_step
```

```         
r <- 1e-8 # crossovers per bp per generation
t <- 1800 # time of admixture (in generations) we computed above

lambda <- r * t * bin_step
y <- dexp(bins, rate = lambda)

plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion")
lines(bins, y, col = "red")
```
