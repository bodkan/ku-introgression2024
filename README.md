Archaic introgression lecture exercises (February 2024)
================

## Exercise 1 – $f_4$ statistic / introgression test

### Introduction

*You sequenced the genomes of four Africans and four Eurasians and got
genotypes from a single chromosome from each of them (so you have
genotypes of four African and four Eurasian chromosomes). Unfortunately,
there’s been a mix up in the lab and you don’t know which one is which!
You only know that they are labeled A, B, C, …, H. What a disaster!*

*Fortunately, you also have genotypes from three other individuals whose
identity you know for certain: an African, a Neanderthal, and a
chimpanzee. This means that you are able to compute an $f_4$ statistic
which will test for evidence of Neanderthal introgression in a given
sample $X$.*

*Can you save the day and determine which of the A, B, C, …, H samples
are African and which are Eurasian based on the following $f_4$
statistic test?*

$$
f_4(\textrm{African}, X; \textrm{Neanderthal}, \textrm{Chimp}).
$$

*Recall that only Eurasians are expected to have appreciable amounts of
Neanderthal ancestry but Africans don’t.*

### Moving over to R

First **load the genotype table into R**:

``` r
gt <- readRDS(url("https://github.com/bodkan/ku-introgression2024/raw/main/genotypes_ex1.rds"))
```

The `gt` data set is a plain R data frame where each column contains the
genotype of that individual (`0` - ancestral allele, `1` - derived
allele).

**Familiarize yourself with the data** by running this R command which
shows information from only the first few genotypes:

``` r
head(gt)
```

**How many positions of the genome do we have genotyped?**

``` r
nrow(gt)
```

You can extract *all the genotypes* of a given individual by using the
`$` or `[[` subsetting operators of R data frames like this:

``` r
gt$African

gt[["African"]]
```

A useful trick for comparing two chromosomes in their entirety is to
rely in the fact that R can perform *vectorized operations* (operations
performed on multiple elements of a vector at once). For instance, if
this gives us the genotypes of an African and Neanderthal chromosome:

``` r
gt[["African"]]

gt[["Neanderthal"]]
```

then we can sum up positions at which those two samples carry the same
allele:

``` r
gt[["African"]] == gt[["Neanderthal"]] # this gives us TRUE/FALSE values 

# sum() treats TRUE as 1 and FALSE as 0, so we can sum everything up!
# -- this gives us the number of positions at which an African carries the same
#    allele as the Neanderthal
sum(gt[["African"]] == gt[["Neanderthal"]])
```

On the other hand, this would count how many alleles are *different*
between our African and chimpanzee chromosome:

``` r
sum(gt[["African"]] != gt[["Chimp"]]) # note the != instead of ==
```

Inside the `sum()` function we can compose multiple logical conditions
to create more complex comparison operations using the `&` operator (AND
operation in mathematical logic).

Armed with this knowledge, we can compute the BABA and ABBA counts using
this bit of R code:

``` r
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
```

From these counts we can get an idea about whether one or the other are
more frequently appearing in the data:

``` r
baba - abba
```

Finally, we can compute an $f_4$ statistic like this, which simply
normalizes the raw difference by how many observations of either BABA or
ABBA site patterns we observed:

``` r
(baba - abba) / nrow(gt)
```

**You know that if `X` is a African, you expect to see roughly the same
count of `BABA` and `ABBA` site patterns, so the difference should “be
about zero”. Use the code above to compute the $f_4$ statistic for all
of your mixed up samples A, B, C, …, H and note down the values you got
for each – which samples are most likely African and which ones are
Eurasian?**

**\[If you are more familiar with R, compute the counts automatically in
a loop of some kind and make a figure.\]**

**What does it mean for this test statistic to “be about zero”? What are
we missing to truly use this as a statistical significance test?**

------------------------------------------------------------------------

If you’re not comfortable with R, feel free to run this in full and
answer the questions based on the results you get:

``` r
X <- c("A", "B", "C", "D", "E", "F", "G", "H")

f4_values <- sapply(X, function(x) {
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
  
  (baba - abba) / nrow(gt)
})

plot(f4_values, xaxt = "n", xlab = "test sample", ylab = "f4(African, X; Neanderthal, Chimp)")
abline(h = 0, lty = 2, col = "red")
axis(side = 1, at = seq_along(X), labels = X)
```

We can see that the samples A-D are consistent with an $f_4$ statistic
“value of about 0”“, meaning that the BABA and ABBA counts were”about
the same”. This is what we would expect for African samples who are not
expected to be closer to a Neanderthal genome than another African.

On the other hand, samples E-H show a much “more negative value of the D
statistic”“, which is consistent with an access of ABBA sites compared
to BABA sites – which arise with an increased sharing of derived alleles
between the sample X and a Neanderthal genome, just as we would expect
when X is of Eurasian ancestry.

**Important:** In this simple example we’re missing confidence intervals
– those would allow us to do a proper statistical test to determine for
which samples we really cannot reject a null hypothesis of no gene flow
from Neanderthals. This would allow us to avoid the vague and
statistically unsatisfying talk about some value being “almost zero”,
and some other value being “much more negative” than that. The
confidence interval for a given $f_4$ statistic would either intersect
the 0 null hypothesis or not. For an example, see Figure 3 in [this
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6485383/).

Real-world software such as
[ADMIXTOOLS](https://github.com/DReichLab/AdmixTools) computes
confidence intervals using a so-called
[bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
procedure across windows along a genome.

------------------------------------------------------------------------

If you want to take a closer look at how the genotype data was prepared
(it was simulated!), you can see the complete code
[here](generate_genotypes.R).

## Exercise 2 – estimating the proportion of Neanderthal ancestry

### Introduction

*Having saved the day by identifying which of the A-H samples are of
African or Eurasian origin by testing which of them appear to carry
evidence of Neanderthal introgression, you now want to estimate how much
of their genome derives from the Neanderthals. In order to do this, you
need to compute a ratio of $f_4$ values as described in the lecture.*

*Of course, in order to compute this $f_4$-ratio estimate, you will need
“another Neanderthal” genome! Luckily, we now have genomes of several
Neanderthals so this is not an issue and a local friendly
bioinformatician merged you `gt` genotype table with the genotypes of
“another Neanderthal”.*

*Estimate the proportion of Neanderthal ancestry in each of your A-H
samples!*

### Moving over to R

**You will be using the same genotype table as in the previous exercise,
with one additional column called `another_Neanderthal`. You can load it
again like this:**

``` r
gt <- readRDS(url("https://github.com/bodkan/ku-introgression2024/raw/main/genotypes_ex2.rds"))
```

As always, **verify that the format of the data set you have matches
what you expect:**

``` r
head(gt)

nrow(gt)
```

From the lecture you know that we can get an estimate for the proportion
of Neanderthal ancestry in a sample $X$ by dividing the rate of allele
sharing between $X$ and a Neanderthal genome (one $f_4$ statistic) by
the rate of allele sharing expected between two Neanderthals (another
$f_4$ statistic).

**Take the $f_4$ values you computed in Exercise 1 for all samples A-H,
and estimate the proportion of Neanderthal ancestry in these samples by
dividing those values by newly computed
$f_4(\textrm{African, another Neanderthal; Neanderthal, Chimp})$.**

``` r
# we can compute the f4 values for everyone (A-H samples as well as
# "another_neanderthal") using the same code as above
X <- c("A", "B", "C", "D", "E", "F", "G", "H",
       "another_Neanderthal") # <---- we added this to our loop

f4_values <- sapply(X, function(x) {
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
  
  (baba - abba) / nrow(gt)
})


# to arrive at the estimate of Neanderthal ancestry, we divide f4 values for
# samples A-H by f4 value comparing the two Neanderthals
proportions <- f4_values / f4_values["another_Neanderthal"]

proportions
```

To make the results clearer to see, let’s visualize them:

``` r
plot(proportions[-length(proportions)] * 100,
     xlab = "test individual", ylab = "proportion of Neanderthal ancestry [%]",
     ylim = c(0, 10))

abline(h = 3, lty = 2, col = "red")
```

**How much Neanderthal ancestry did you estimate in Africans vs
Eurasians? Do those numbers fit what you’ve learned from the lecture?**

**We did we not plot the proportion of Neanderthal ancestry in the very
last item of the `proportions` variable? What does that last element of
the vector `proportions` contain and why?**

------------------------------------------------------------------------

If you want to take a closer look at how the genotype data was prepared
(it was simulated!), you can see the complete code
[here](generate_genotypes.R).

## Exercise 3 – dating Neanderthal admixture

### Introduction

*You sequenced 100 diploid genomes from a Eurasian population and are
interested in estimating the time of Neanderthal introgression into the
ancestors of this population. The literature suggests that the
introgression happened around 55 thousand years ago – does this also
apply to the population that you sequenced?*

*To be able to do this, you ran an inference software which gives you
the exact coordinates of Neanderthal DNA tracts present in every
Eurasian genome that you sequenced. This of course means that you also
know the lengths of each of those tracts.*

*Use the distribution of the Neanderthal tract lengths in your
population to estimate the time of Neanderthal introgression!*

### Moving over to R

**First load the table with coordinates of all Neanderthal tracts** into
R:

``` r
tracts <- readRDS(url("https://github.com/bodkan/ku-introgression2024/raw/main/tracts.rds"))
```

The `gt` data set is a plain R data frame where each individual’s column
contains the genotype of that individual’s chromosome (`0` - ancestral
allele, `1` - derived allele).

**Familiarize yourself with the data** by running this R command which
shows information from only the first few genotypes:

``` r
head(tracts)
```

**For how many individuals do we have information about Neanderthal DNA
tracts that they carry?**

``` r
length(unique(tracts$individual))
```

It looks like the inference software (or a helpful bioinformatician)
binned each tract according to its length. **What does the distribution
of Neanderthal tract lengths looks like in your data?** Knowing that
recombination has acted on the introgressed Neanderthal DNA over time,
each generation, suggests that the distribution should look exponential
– do you see this in the data? To answer this, plot the proportion of
tracts in each bin.

``` r
# get the bin numbers
bins <- sort(unique(tracts$bin))

# count the tracts in each bin and compute the proportion of tracts in each bin
counts <- as.integer(table(tracts$bin))
props <- counts / sum(counts)

plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion of tracts")
```

The distribution does, indeed, look quite exponential. **Let’s try to
use this to date Neanderthal introgression using the information encoded
in the distribution of tract lengths!**

As we know, over time since admixture, recombination breaks up longer
haplotypes into shorter ones, regularly almost like a clock. And it
turns out that the distribution of tract lengths after time $t$ follows
exponential decay, leading to the distribution of tract lengths $x$ to
have the following form:

$$
\textrm{tract length of }~x \sim \exp(-\lambda x) = \exp(-r t x),
$$

where the $\lambda$ parameter determines the *rate* of exponential decay
which can be computed, under some simplifying assumptions, as the
product of the recombination rate (traditionally in humans with value of
about $10^{-8}$) and $t$ which is the time since admixture:

$$
\lambda = rt
$$

The $t$ in this equation latter is our unknown we’re trying to compute
in this exercise.

It also turns out that the expected value of this exponential
distribution (which [can be
computed](https://en.wikipedia.org/wiki/Exponential_distribution#Mean,_variance,_moments,_and_median)
simply as $1 / \lambda = 1 / rt$) gives us the theoretical expression
for the expected tract length after time $t$.

Of course, you can also compute this expectation from the data – **get
an estimate of the expected length of an introgressed fragment after
(unknown) time $t$ by computing the average introgressed tract length
observed in data:**

``` r
L <- mean(tracts$length)
L # length in units of base pairs
```

**What is the average length of a chunk of a Neanderthal DNA in a
Eurasian genome?**

Taking the simple math above together and doing a little algebra, we can
express the average expected length of an introgressed fragment after
time $t$ using this formula:

$$
\textrm{average tract length}~L = \frac{1}{\lambda} = \frac{1}{rt}
$$

**But because we know $L$ (average tract length as you just computed)
and $r$ (recombination rate of $10^{-8}$), this means we can estimate
the time since the admixture by a simple rearrangement of the above
equation to separate the time of introgression $t$:**

$$
t = \frac{1}{rL}
$$

Note that this time estimate will be in units of generations, so we’ll
have to multiply this quantity by the generation time (roughly 30 years
for humans) to get time in years before the present.

**Use this simple equation to compute the estimate of the admixture
time:**

``` r
r <- 1e-8 # crossovers per bp per generation
L <- mean(tracts$length) # average tract length after time t

t <- 1 / (r * L)
t

# convert the time of introgression to years before present assuming
# generation time of 30 years
t * 30
```

You should get a value will be quite close to ~55 thousand years ago, an
estimate which is often found in the literature as the time when
Neanderthals and anatomically modern humans interbred!

------------------------------------------------------------------------

**Bonus:**

As a last sanity check, if we use this time to compute the rate of
exponential decay $\lambda$, we should get a nice fit of the theoretical
exponential decay curve over the empirical counts of tract lengths in
each bin. As a reminder, this is the decay:

``` r
plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion of tracts")
```

Let’s try if we can plot the theoretical exponential decay from the
estimated time of admixture.

First, because the exponential tract decay you plotted above shows the
decay of the length of introgressed tracts in bins, the $\lambda$
parameter determines the decay expected for tracts falling into to
successive bins (so decay of longer and longer tracts). However, our
recombination rate $r$ which features in the equation to compute
$\lambda$ above describes recombination in units of base pairs, not
bins. In order to be able to plot the theoretical exponential decay
curve, we therefore first need to compute the average increase in tract
lengths from bin to bin:

``` r
# compute the average length in each bin
average_bins <- aggregate(length ~ bin, data = tracts, FUN = mean)

# compute the difference between bins to get "average bin step size"
bin_step <- mean(diff(average_bins$length))

bin_step
```

With this, we can overlay the theoretical exponential decay expectation
across the empirical distribution of tract lengths in different bins:

``` r
r <- 1e-8 # crossovers per bp per generation
t <- 1800 # time of admixture (in generations) we computed above

# plot the empirical proportions of tracts in each bin
plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion of tracts")

# overlay theoretical exponential density curve assuming recombination rate r
# and time of introgression t generations ago
lambda <- r * bin_step * t
y <- dexp(bins, rate = lambda)
lines(bins, y, col = "red", lty = 2)
```

------------------------------------------------------------------------

If you want to take a closer look at how the tracts data was prepared
(it was simulated!), you can see the complete code
[here](generate_tracts.R).
