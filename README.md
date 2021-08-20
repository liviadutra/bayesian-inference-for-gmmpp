# Algorithms for examples in "Exact and computationally efficient Bayesian inference for generalized Markov modulated Poisson processes"

The algorithms are coded in Ox

Please, download Ox Console (including OxEdit) from the link https://www.doornik.com/download.html.
Ox is free for academic use.

For all the examples, there is a main file, named main_code.ox, which is the one to be run and that contains all the required specifications to be set. For each of the examples provide, there are a few specifications that can be modified by the user by following the comments in the code.

The datasets used in the manuscript are saved in files named data.mat inside the respective folder of each example.
The first row of the data file contains its number of rows and columns (Ox format).

We are attaching five of the examples presented in the paper. The codes for any of the other examples in the manuscript are available upon request to the authors.


## Scenario A2 - Section 4.1

Folder: [Section_4.1](https://github.com/liviadutra/bayesian-inference-for-gmmpp/tree/main/Section_4.1)


## Example 3 - Section 4.3

Folder: [Section_4.3](https://github.com/liviadutra/bayesian-inference-for-gmmpp/tree/main/Section_4.3)


## Coal mining disasters - Section 5.1

Folder: [Section_5.1](https://github.com/liviadutra/bayesian-inference-for-gmmpp/tree/main/Section_5.1)


## BRLxUSD exchange rate - Section 5.2

Folder: [Section_5.2](https://github.com/liviadutra/bayesian-inference-for-gmmpp/tree/main/Section_5.2)


## COVID-19 Romania - Appendix E

Folder: [Appendix_E](https://github.com/liviadutra/bayesian-inference-for-gmmpp/tree/main/Appendix_E)


## The output files

The output files for each example are:

- "lambda_mean.mat": posterior mean of the IF for a grid of 1000 equaly spaced points in \[0,S\].
- "lambda_ci.mat": 95% pointwise credibility interval for the IF for the same grid with 1000 points. Row 1: quantile(0.025), row 2: quantile(0.975).
- "psi.mat": MCMC chain of the \psi parameters - one parameter per column.
- "logposterior.mat": MCMC chain of the log-posterior density (logarithm of the non-normalized posterior density)
- "integral.mat": MCMC chain of the integrated IF in \[0,S\].
- "theta.mat": MCMC chain of the \theta parameters - one parameter per column. No output if the Q-matrix is fixed.

If prediction is performed:

- "Npredic.mat": MCMC chain of the number of observetions in the predicted interval.
- "integralpredic.mat": MCMC chain of the integrated IF in the predicted interval.
- "timepredict.mat": MCMC chain of the time instant when the IF hits a predefined value in the epidemic model.
