(in-package "CL-USER")
;;Modified with declarations
(declaim (inline signp dfloat))

;;Common-lisp compatibility for Numerical Recipes

;;Fortran-offset array-reference.  Fortran arrays start at 0, lisp
;;arrays at 1.
(defmacro fref (array &rest indicies)
  `(aref ,array
	 ,@(mapcar #'(lambda (index) `(1- ,index)) indicies)))

;;Return X with the sign of Y.  Various routines like this exist
;;in the book, but they differ about what to do when Y is 0.
;;This one gives X rather than -X.
(defun signp (x y)
  (if (minusp y) (- x) x))


;;Convert a number to a double-float
(defun dfloat (x)
  (coerce x 'double-float))

;;Check that a number returned from a user function is a double-float
(defun dfloat-check (x)
  (check-type x double-float)
  x)

