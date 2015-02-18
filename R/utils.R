
# anomaly SVMs
# Wikum Dinalankara
# 2015-01-07

# utility functions
# ==================

# get hat matrix and projection of anomalies
# Z : normal matrix
# X : anomaly matrix
# Note: columns are samples
get_projections = function(Z, X=NULL){

    Qz = qr(Z)
    H = qr.Q(Qz) %*% t(qr.Q(Qz))

    P = NULL
    if(! is.null(X))
        P = H %*% X

    list(H=H, X=P)

}

# check whether a matrix is positive semi-definite
isPD = function(X){ sum(eigen(X)$values < 0) == 0 }

# convert a binary vector of {0, 1} to {-1, 1}
y_conv = function(y){  
    if( (length(unique(y)) == 2) && (min(y) == 0) ){
        y[y == 0] = -1
    }
    y
}

