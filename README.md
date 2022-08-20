# gwas.winners.curse: a tool for probabilistic correction of Winner's Curse in two-stage GWAS

<!-- badges: start -->
[![R-CMD-check](https://github.com/cpalmer718/gwas-winners-curse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cpalmer718/gwas-winners-curse/actions/workflows/R-CMD-check.yaml)

[![codecov](https://codecov.io/gh/cpalmer718/gwas-winners-curse/branch/default/graph/badge.svg?token=IIPN9LBZNQ)](https://codecov.io/gh/cpalmer718/gwas-winners-curse)
<!-- badges: end -->



## Overview

This is an `R` package designed to adjust GWAS discovery results to correct
for the [Winner's Curse](https://en.wikipedia.org/wiki/Winner%27s_curse) or
regression to the mean phenomenon. The resultant corrected values can more
appropriately be used to evaluate success or failure of replication studies,
or to compute power to replicate signals _a priori_.

This implementation will serve as a total replacement for an old C++ program
that attempted the same goal. More documentation will be added as this package
comes into existence.

## Documentation

### Overview

This package replaces a C++ program accompanying an [old publication on the Winner's Curse](https://doi.org/10.1371/journal.pgen.1006916).
If you're interested in the original implementation, please see [this tag](https://github.com/cpalmer718/gwas-winners-curse/tree/1.0.0-beta.2);
but trust me, newer is better in this case.

I encourage anyone who's interested in the Winner's Curse (or regression to the mean) to read the aforementioned citation. I can
also suggest some really great papers on the topic.

- [Zhong and Prentice](https://doi.org/10.1093%2Fbiostatistics%2Fkxn001) and [Xiao and Boehnke](https://doi.org/10.1002/gepi.20398)
  were the inspiration and statistical engine behind the original paper. I wholeheartedly recommend both papers.
- There has been plenty of interesting work on the Winner's Curse in GWAS in the intervening years, more than I could possibly enumerate
  here. I will mention in particular [Zou _et al._](https://www.biorxiv.org/content/10.1101/856898v1.full.pdf) who take the work
  of my manuscript and improve it by addressing something we only mentioned in passing: that certainly the Winner's Curse is only one
  reason you might observe discordances between discovery and replication. They go past what we did and actually attempt to model
  such effects. It is a pleasure to see someone make use of any of the information I worked on back then.
- As was the case in the original publication, I will mention that I surveyed several hundred GWAS papers for that original manuscript.
  I really appreciate the effort that everyone involved in all of those papers put into their work, and how it enabled my work to happen.
  The full list of citations is [here](https://doi.org/10.1371/journal.pgen.1006916.s012) for anyone who is curious; and the parsed 
  variant summary information from those papers is available in the supplement of the original manuscript, and is also now stored in 
  `tests/testthat/testthat_files/correct_winners_curse` for posterity.

### Usage

I recommend installing the package with `devtools::install_github()`:

```r
library(devtools)
devtools::install_github("cpalmer718/gwas-winners-curse", ref = "default")
library(gwas.winners.curse)
```

To actually run Winner's Curse correction, prepare an input file containing the following columns:

- variant regression coefficient (discovery)
- variant regression standard error (discovery)
- variant sample size (discovery; please make it as accurate as possible)
- variant allele frequency (discovery; as with sample size, accuracy is important)
- (optional) phenotype distribution mean; only has to be a column if it differs per-row
- (optional) p-value threshold for discovery; required if phenotype distribution mean is also specified per-row
- variant regression coefficient (replication)
- variant regression standard error (replication)
- variant sample size (replication)
- variant allele frequency (replication)

For the moment, only the discovery and mean/threshold entries are used. Replication values will be used in some
other functions in the near future.

Once this is prepared, run Winner's Curse correction as follows:

```r
input.file <- "my_input.tsv"
output.file <- "my_output.tsv"
trait.mean <- NA ## can be specified as a parameter if missing from input file
p.threshold <- NA ## can be specified as a parameter if missing from input file
gwas.winners.curse::correct.winners.curse(input.file,
                                          output.file,
                                          trait.mean,
                                          p.threshold,
                                          header = TRUE,
                                          sep = "\t")
```

The results file will contain:

- the input data, possibly with modified header
- trait mean and p-value threshold columns appended, if absent from input
- variant regression coefficient (adjusted with MLE method)
- variant regression lower 95% confidence interval (adjusted with MLE method)
- variant regression upper 95% confidence interval (adjusted with MLE method)
- variant regression coefficient (adjusted with MSE method)
- variant regression lower 95% confidence interval (adjusted with MSE method)
- variant regression upper 95% confidence interval (adjusted with MSE method)

A few things to note for anyone comparing this information to runs with the old program:

- the MLE standard error is no longer reported. the program uses the input standard error.
  this is more appropriate for linear regression than what it was doing before.
- the MSE coefficient will differ in most cases from the original program's. MSE was a random
  addition that didn't ever seem to work very well originally, and so it was actually disabled
  and just copied the MLE version. that behavior is reset to the correct behavior in the R package;
  my apologies for any confusion.

## Version History
See [changelog](CHANGELOG.md) for more information.
 * 19 Aug 2022: remove old implementation, begin R package

## How to contribute to development

### Step 1: Set up a development environment (OSX and Linux only)

- If needed, install miniconda by following the steps [here](https://docs.conda.io/en/latest/miniconda.html).
- If needed, install [mamba](https://github.com/mamba-org/mamba): `conda install mamba`
- Clone a copy of this repo: 

```
git clone https://github.com/cpalmer718/gwas-winners-curse.git
```

- Navigate into the repo directory: `cd gwas-winners-curse`
- Create a conda environment with, minimally, the dependencies defined in `r-dev.yaml`.  Make sure to activate your dev environment whenever you are writing/committing code!

```
# create the env
mamba env create -f r-dev.yaml

# activate the env
conda activate r-dev
```

- Install [commitizen](https://github.com/commitizen/cz-cli) as follows

```
npm install -g commitizen cz-conventional-changelog
commitizen init cz-conventional-changelog --save-dev --save-exact
```

- Set up pre-commit hook scripts.  This will apply linting and check for some common code formatting errors every time you commit.  See https://pre-commit.com/ for more details.  

```
pre-commit install
```

- Install pre-commit in R (either in an R terminal or in Rstudio):

```{r}
install.packages("precommit")
```

### Step 2: Select an issue to work on, or submit one that you'd like to address

See the current [issues](https://github.com/cpalmer718/gwas-winners-curse/issues) for this project.

### Step 3: Contribute code

- All development work should branch off of the `dev` branch.  Make sure you're on the right branch: `git checkout dev`
- Make sure your repository is up-to-date with the remote: `git pull origin dev`
- Create a new branch named for the feature you're going to work on: `git checkout <feature_branch>`
- Write code and commit often!
    - Stage changes with `git add .`
    - Commit code with `git cz`; make sure to cite the issue number you're working on
    - Push your changes to the remote repository with `git push origin <feature_branch>`
- When you're all done, submit a pull request [here](https://github.com/cpalmer718/gwas-winners-curse/pulls).  Other developers will review your code, make comments, and merge in your changes when ready!
