BELA
====

Bayesian Estimation of Latent Abundance for microbiome data

A Variational Inference suite for NB Ordination Analysis as described in [Ren et al (2016)](http://arxiv.org/abs/1601.05156).

For MCMC version of the same model, see Boyu Ren's [Repo](https://github.com/boyuren158/DirichletFactor)

Documentation is located [here](https://www.dropbox.com/s/92z4gncbmvrxbgs/LRVI.pdf?dl=0)

Observations
---
* In pollution example: 
    * Factor is identifiable when sample Y by row, with mean(Y) = log(Theta)
    * It is not identifiable when sample by column, or mean(Y) = Theta