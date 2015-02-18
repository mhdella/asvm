#' asvm.
#'
#' @name asvm
#' @docType package
#' @include svmfunctions.R

classificationSvm = setRefClass("classificationSvm",

	slots = list(
			
			svm_function = "function",
			kernel_function = "function",
			params = "list",
			optimization = "character",
			type = "character",
			penalization_type = "character",
			penalty = "numeric",
			y = "vector",
			X = "matrix"

		),

	prototype = list(

			svm_function = fit_nu_svm_qp,
			kernel_function = poly_kernel,
			params = list(param=1),
			optimization = "qp", # options: qp or ipop
			type = "regular", # options: regular, anomaly, one-class, anomaly-one-class, expanded
			penalization_type = "nu", # options: nu or c
			penalty = .1

		)

	# methods:
	# 


)