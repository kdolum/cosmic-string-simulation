;;;Interface to Gnu Scientific Library in Lisp

(in-package "CL-USER")

(declaim (inline besselk gsl-sf-bessel-knu))

;;Direct access much faster than gsl:cylindrical-bessel-K
(define-alien-routine "gsl_sf_bessel_Knu" double
    (nu double)
    (z double))

(defun besselk (nu z)
  (gsl-sf-bessel-knu nu z))

(define-alien-routine "gsl_sf_bessel_Jnu" double
    (nu double)
    (z double))

(defun besselj (nu z)
  (declare (double-float z))
  (setq nu (double-float nu))
  (handler-case (gsl-sf-bessel-jnu nu z)
    ;;Give zero if underflow detected by gsl.  Changing floating point modes can't fix this, because
    ;;the error is detected explicitly by the code.
    (gsll:underflow () 0.0)))

;;Again direct access faster
(declaim (inline ci gsl-sf-ci))
(define-alien-routine "gsl_sf_Ci" double
    (x double))

(defun ci (x)
  (gsl-sf-ci x))

;;There is no user interface to constant folding in sbcl.  Rather than using compiler internals,
;;I'm doing it by hand here

;;We allow the constant to have any type, so we can use 1/3, for example, but if it is not
;;a constant we don't emit code to convert it.
;;We ignore the lexical environment here.
(defmacro gamma (z)
  (if (constantp z) (gsl:gamma (double-float (eval z)))
    `(gsl:gamma ,z)))

;;Interface to Gaussian adaptive integration
(defun integrate (function a b
			   &key
			   (absolute-error 1e-6)
			   (relative-error 1e-6)
			   (limit 1000))
  ;;Handle singularities because we often have discontinuities at the edges of the caching routine ranges.
  (gsl:integration-qags function a b absolute-error relative-error limit))

;;This is needed to produce the right kind of "grid" for integration limits
(setq grid:*default-grid-type* 'grid:foreign-array)

;;Tell what's going on in integration
(defun debug-integrate (function a b
			   &key
			   (absolute-error 1e-6)
			   (relative-error 1e-6)
			   (limit 1000))
  (let ((table (make-array 1000 :adjustable t :fill-pointer 0)))
   (unwind-protect
       (gsl:integration-qags
	#'(lambda (x)
	    (let ((result (funcall function x)))
	      (vector-push-extend (list x result) table)
	      result))
	a b absolute-error relative-error limit)
     (plot-data-list (sort (coerce table 'list) #'< :key #'car) "integrand" :styles :points :reuse t)
     (break "plot")
     )))
 
;;Eigenvalues and eigenvectors of real	a symmetric matrix
(defun eigensystem-real-symmetric (m)
  (multiple-value-bind (eigenvalues eigenvectors)
      (gsl:eigenvalues-eigenvectors (grid:copy m :grid-type 'grid:foreign-array))
    (values (grid:copy eigenvalues :grid-type 'array)
	    (grid:copy eigenvectors :grid-type 'array))))

(defun eigenvalues-nonsymmetric (m)
  (grid:copy (gsl:eigenvalues-nonsymm (grid:copy m :grid-type 'grid:foreign-array))
	     :grid-type 'array))
    
;;Find eigensystem of nonsymmetric matrix.  Eigenvectors are the columns of the returned array.
(defun eigensystem-nonsymmetric (m)
  (multiple-value-bind (eigenvalues eigenvectors)
      (gsl:eigenvalues-eigenvectors-nonsymm (grid:copy m :grid-type 'grid:foreign-array))
    (values (grid:copy eigenvalues :grid-type 'array)
	    (grid:copy eigenvectors :grid-type 'array))))