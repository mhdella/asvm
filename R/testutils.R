
# anomaly SVMs
# Wikum Dinalankara
# 2015-01-07

# testing functions
# ==========================

# create simulation data
make_simulation_data = function(n, p, m){

    # make normals
    n_means = rnorm(m, mean=3.5, sd=1)
    n_sds = abs(rnorm(m, mean=2, sd=1))
    
    Z = t(sapply(1:m, function(j) rnorm(p, mean=n_means[j], sd=n_sds[j]) ))
    
    low_means = jitter(n_means + 1)
    high_means = jitter(n_means + 1.5)

    low_sds = jitter(n_sds + 2)
    high_sds = jitter(n_sds + 3)

    L = t(sapply(1:m, function(j) rnorm(n/2, mean=low_means[j], sd=low_sds[j]) ))

    H = t(sapply(1:m, function(j) rnorm(n/2, mean=high_means[j], sd=high_sds[j]) ))

    y = rep(c(-1, 1), c(dim(L)[2], dim(H)[2]))
    
    list(Z=Z, X=cbind(L, H), y=y)

}

# make training and testing subsets: take p portion of each class




# divide data for cross validation

# tuning nu-paramater for linear kernel

# sigma and nu tuning for rbf kernel


