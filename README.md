# gwas.winners.curse: a tool for probabilistic correction of Winner's Curse in two-stage GWAS

<!-- badges: start -->
[![R-CMD-check](https://github.com/lightning-auriga/gwas-winners-curse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lightning-auriga/gwas-winners-curse/actions/workflows/R-CMD-check.yaml)
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

TBD

## Version History
See [changelog](CHANGELOG.md) for more information.
 * 19 Aug 2022: remove old implementation, begin R package

## How to contribute to development

### Step 1: Set up a development environment (OSX and Linux only)

- If needed, install miniconda by following the steps [here](https://docs.conda.io/en/latest/miniconda.html).
- If needed, install [mamba](https://github.com/mamba-org/mamba): `conda install mamba`
- Clone a copy of this repo: 

```
git clone https://github.com/lightning-auriga/gwas-winners-curse.git
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

See the current [issues](https://github.com/lightning-auriga/gwas-winners-curse/issues) for this project.

### Step 3: Contribute code

- All development work should branch off of the `dev` branch.  Make sure you're on the right branch: `git checkout dev`
- Make sure your repository is up-to-date with the remote: `git pull origin dev`
- Create a new branch named for the feature you're going to work on: `git checkout <feature_branch>`
- Write code and commit often!
    - Stage changes with `git add .`
    - Commit code with `git cz`; make sure to cite the issue number you're working on
    - Push your changes to the remote repository with `git push origin <feature_branch>`
- When you're all done, submit a pull request [here](https://github.com/lightning-auriga/gwas-winners-curse/pulls).  Other developers will review your code, make comments, and merge in your changes when ready!
