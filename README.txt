This repository contains code used for the ppaer ``Shared Subspace Models for Multi-Group Covariance Estimation''
by Franks and Hoff. (https://arxiv.org/abs/1607.03045)

The EM algorithm used to estimate the maximum marginal likelihood of a shared subsace
(described in Section 3.1 of the paper) can be run using the optimV function in helper.R

Code to run Bayesian inference for the projected covariance matrices (Section 3.2)
is provided in fit-subspace.R

The code which was used for the simulation studies (Section 4) can be found in
simulationFigures, bayes-simulation.R, bayes-coverage.R, Sdim-test.R

The code used for the figures on asymptotic bias (Section 5)
is in asymptotics.R

The code used for the two applied examples from Section 5 can be found
in fit-leukemia.R and fit-dmelan.R


