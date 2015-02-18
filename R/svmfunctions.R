
# anomaly SVMs
# Wikum Dinalankara
# 2015-01-07

# svm optimization functions
# ==========================

INFI = Inf

# convert svmpath kernel functions for use with row-feature matrices

poly_kernel = function(x, y=x, param.kernel=1){

    poly.kernel(t(x), t(y), param.kernel)
}

radial_kernel = function(x, y=x, param.kernel=1){

    radial.kernel(t(x), t(y), param.kernel)
}

# classification svms
#======================================================

# 

#' C-svm for classification, qp solver
#' 
#' @param y vector of class labels as -1 or 1
#' @param X matirx of input data, columns as samples
#' @param kernel kernel function, polynomial kernel by default
#' @param param.kernel kernel parameter, 1 by default
#' @param c cost parameter, 1 by default
#' @param epsilon threshold for selecting support vectors, 1e-6 by default
#' @export

fit_c_svm_qp = function(y, X, kernel=poly_kernel, param.kernel=1, c=1, epsilon=1e-6){

    if(length(y) != dim(X)[2]) stop("incorrect length of y vector")

    n = dim(X)[2]

    K = kernel(X, X, param.kernel)

    Dmat = diag(y) %*% K %*% diag(y)
    dvec = rep(1, n)
    
    Amat = t(rbind(y, diag(rep(1, n)), diag(rep(-1, n))))
    bvec = c(0, rep(0, n), rep(-c, n))
    meq = 1

    sol = solve.QP(Dmat, dvec, Amat, bvec, meq, factorized=FALSE)

    alphas=sol$solution
    bvals = c(0)
    tryCatch({
        bvals = sapply(which(alphas > epsilon), function(j) 1 - sum(alphas * y * K[, j]))
    }, error = function(e){ print(e) })

    list(alphas=alphas, optim=sol, K=K, bvals=bvals, b=mean(bvals))
}

# nu-svm for classification, qp solver
fit_nu_svm_qp = function(y, X, kernel=poly_kernel, param.kernel=1, c=1, epsilon=1e-6){

    if(length(y) != dim(X)[2]) stop("incorrect length of y vector")

    n = dim(X)[2]

    K = kernel(X, X, param.kernel)

    K = as.matrix(nearPD(K)$mat)

    Dmat = diag(y) %*% K %*% diag(y)

    Dmat = as.matrix(nearPD(Dmat)$mat)

    dvec = rep(0, n)
    
    Amat = t(rbind(y, rep(1, n), diag(rep(1, n)), diag(rep(-1, n))))
    bvec = c(0, nu, rep(0, n), rep(-1/n, n))
    meq = 1

    sol = solve.QP(Dmat, dvec, Amat, bvec, meq, factorized=FALSE)

    alphas=sol$solution

    bvals = c(0)
    tryCatch({
        bvals = sapply(intersect(which(alphas > epsilon), which(alphas < ((1/n)-epsilon))) , function(j) 1 - sum(alphas * y * K[, j]))
    }, error = function(e){ print(e) })

    list(alphas=alphas, optim=sol, K=K, bvals=bvals, b=mean(bvals))
}

# nu-svm for classification, ipop solver
fit_c_svm_ipop = function(y, X, kernel=poly_kernel, param.kernel=1, c=1, epsilon=1e-6){

    if(length(y) != dim(X)[2]) stop("incorrect length of y vector")

    n = dim(X)[2]

    K = kernel(X, X, param.kernel)

    K = as.matrix(nearPD(K)$mat)

    c = 1/n

    cc = rep(0, n)
    H = diag(y) %*% K %*% diag(y)

    H = as.matrix(nearPD(H)$mat)

    A = rbind(matrix(y, 1, n), rep(1, n))
    b = c(0, nu)
    r = c(0, INFI)
    l = rep(0, n)
    u = rep(c, n)

    sol = ipop(cc, H, A, b, l, u, r, margin = epsilon)

    alphas=primal(sol)

    bvals = c(0)
    tryCatch({
        bvals = sapply(intersect(which(alphas > epsilon), which(alphas < ((1/n)-epsilon))) , function(j) 1 - sum(alphas * y * K[, j]))
    }, error = function(e){ print(e) })

    list(alphas=alphas, optim=sol, K=K, bvals=bvals, b=mean(bvals))
}

# exapnded anomaly nu-svm for classification, qp solver
fit_exapnded_nu_svm_qp = function(y, X, Z=NULL, kernel=poly_kernel, param.kernel=1, nu=0.1, epsilon=1e-6){

    if(length(y) != dim(X)[2]) stop("incorrect length of y vector")

    n = dim(X)[2]
    p = dim(Z)[2]
    m = dim(X)[1]

    # K is a p*n x p*n matrix

    X1 = matrix(0, m, p*n)
    Y = rep(0, p*n)
    for(i in 1:n){
        for(j in 1:p){
            X1[, (p*(i-1))+j] = X[, i] - Z[, j]
            Y[(p*(i-1))+j] = y[i]
        }
    }

    K = kernel(X1, X1, param.kernel)

    Dmat = diag(Y) %*% K %*% diag(Y)
    # d should be zero for nu-svm
    dvec = rep(0, n*p)
    Amat = t(rbind(Y, rep(1, n*p), diag(rep(1, n*p)), diag(rep(-1, n*p))))
    bvec = c(0, nu, rep(0, n*p), rep(-1/(n*p), n*p))
    meq = 1

    sol = solve.QP(Dmat, dvec, Amat, bvec, meq, factorized=FALSE)    

    alphas=sol$solution

    bvals = c(0)
    tryCatch({
        bvals = sapply(intersect(which(alphas > epsilon), which(alphas < ((1/n)-epsilon))) , function(j) 1 - sum(alphas * y * K[, j]))
    }, error = function(e){ print(e) })

    list(alphas=alphas, optim=sol, K=K, bvals=bvals, b=mean(bvals), X1=X1, Y=Y)
}

# one class svms
#======================================================

# one-class svm for regression, qp solver
fit_one_svm_qp = function(y, X, kernel=poly_kernel, param.kernel=1, nu=0.1, epsilon=1e-6){

    n = dim(Z)[2]

    K = kernel(Z, Z, param.kernel)

    K = as.matrix(nearPD(K)$mat)

    c = 1/(nu * n)

    Dmat = K

    Dmat = as.matrix(nearPD(Dmat)$mat)

    dvec = rep(0, n)
    Amat = t(rbind(rep(1, n), diag(rep(1, n)), diag(rep(-1, n))))
    bvec = c(1, rep(0, n), rep(-c, n))
    meq = 1

    sol = solve.QP(Dmat, dvec, Amat, bvec, meq, factorized=FALSE)

    alphas=sol$solution
    bvals = c(0)
    tryCatch({
        bvals = sapply(intersect(which(alphas > epsilon), which(alphas < (c-epsilon))) , function(j) 1 - sum(alphas * K[, j]))
    }, error = function(e){ print(e) })

    list(alphas=alphas, optim=sol, K=K, bvals=bvals, b=mean(bvals))
}

# one-class svm for regression, ipop solver
fit_one_svm_ipop = function(y, X, kernel=poly_kernel, param.kernel=1, nu=0.1, epsilon=1e-6){

    n = dim(Z)[2]

    K = kernel(Z, Z, param.kernel)

    K = as.matrix(nearPD(K)$mat)

    c = 1/(nu * n)

    cc = rep(0, n)
    H = K

    H = as.matrix(nearPD(H)$mat)

    A = matrix(1, 1, n)
    b = c(1)
    r = c(0)
    l = rep(0, n)
    u = rep(c, n)

    sol = ipop(cc, H, A, b, l, u, r, margin = epsilon)

    alphas=primal(sol)
    bvals = c(0)
    tryCatch({
        bvals = sapply(intersect(which(alphas > epsilon), which(alphas < (c-epsilon))) , function(j) 1 - sum(alphas * K[, j]))
    }, error = function(e){ print(e) })

    list(alphas=alphas, optim=sol, K=K, bvals=bvals, b=mean(bvals))
}

# regression svms
#======================================================

make_regression_K = function(X, kernel_function, kernel_param){

    # columns are samples

    K = kernel_function(X, X, kernel_param)
    K = as.matrix(nearPD(K)$mat)

    Kstar = rbind(cbind(K, -K), cbind(-K, K))

    Kstar = as.matrix(nearPD(Kstar)$mat)

    list(Kstar=Kstar, K=K)
}

# nu-svm for regression, qp solver
fit_nu_regression_svm_qp = function(y, X, kernel=poly_kernel, param.kernel=1, nu=0.1, u = 0.001, epsilon=1e-6){

    if(length(y) != dim(X)[2]) stop("incorrect length of y vector")

    u = 0.1
    tryCatch({
        du = diff(sort(y))
        u = abs( min(du[du != 0]) )/2
    }, error = function(e){
        print(e)
    })

    u = 0.001
 
    n = dim(X)[2]

    va = 1 / (nu * n)

    # K = kernel(X, X, param.kernel)

    Kmats = make_regression_K(X, kernel, param.kernel)
    Dmat = Kmats$Kstar

    K=Kmats$K

    dvec = c(y-u, -y-u)

    Amat = t(rbind( c(rep(1, n), rep(-1, n) ),
              diag(rep(1, n+n)),
              diag(rep(-1, n+n))
            ))
    bvec = c(0, rep(0, n+n), -c( rep(va, n+n) )) 

    meq = 1

    sol = solve.QP(Dmat, dvec, Amat, bvec, meq, factorized=FALSE)

    # calculating alphas and b

    full_alphas=sol$solution

    alpha_i = full_alphas[1:n]
    alpha_i_star = full_alphas[(n+1):(n+n)]

    alphas = alpha_i - alpha_i_star

    bvals = c(0)
    tryCatch({
        sel = intersect(which(full_alphas > epsilon), which(full_alphas < (va-epsilon)))
        bfx = sapply(sel, function(j) sum(alphas * Kmats$Kstar[, j]))
        bvals = (c(y, y) + rep(c(u, -u), each=n))[sel] - bfx 
    }, error = function(e){ 
        #print(e) 
    })

    list(alphas=alphas, optim=sol, Kstar=Kmats$Kstar, K=Kmats$K, bvals=bvals, b=mean(bvals),
        alpha_i=alpha_i, alpha_i_star=alpha_i_star, y=y)
}



