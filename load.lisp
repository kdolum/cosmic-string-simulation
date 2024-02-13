(in-package "CL-USER")

(eval-when (:compile-toplevel)
  (error "Don't compile this file.  Load the source."))

(defparameter simulation-files
  '("definitions"
    "run-program"
    "gnuplot"
    "nr-header"				;Numerical recipes header
    "lu" "fft" "svdcmp" "jacobi" "interpolation" "integration"	;Numerical recipes code
    ;; "bessel"				;NR Bessel function code.  GSL is better
    ;; "rk" ; Not used except in deleted smoothing code
    "dop853"
    "gsl"				;Interface to GSLL
    "gc"
    "heap"
    "thread"
    "resources"
    "vector"
    "matrix"
    "handle"
    "geometry"
    "tag"
    "diamond"
    "junction"
    "successors"
    "cells"
    "communication"
    "calendar"
    "output"
    "initial-conditions"
    "intersections"
    "input"
    "solve"
    "evolve"
    "test"
    "peer"
    "submit"
    "b-splines"
    "nufft"
    "interpolating-function"
    "plot"
    "recreate"
    "gravitational"
    "gw-near-cusp"
    "paraview"
    "manager"
    "worker"
    "backreaction"
    "parallel-backreaction"
    "cusps"
    "debug"
    "blackholes"
    ))
 
;;Everything including files that are not compiled: this file, system, local-modifications
(defparameter all-simulation-files
  (list* "load" "system" "local-modifications" simulation-files))  

(load (merge-pathnames "system.lisp" *load-pathname*))

;;Some operating settings
;;(proclaim '(optimize debug))		;Allow debugging at cost of speed

;;Some internals depend on sbcl version
(when (char> (char (lisp-implementation-version) 0) #\1)
  (pushnew :sbcl2 *features*)
  (when (find-symbol "CLASSOID-WRAPPER" "SB-KERNEL")
    (pushnew :sbcl-wrappers *features*)))

(setq *compile-print* nil)		;Don't tell about each fn compiled
(setq *load-verbose* t)			;Do tell about files loaded
(setq *print-right-margin* (max (or *print-right-margin* 0) 100)) ;print 4 floats on a line
(setq *read-default-float-format* 'double-float)
;;This gives us double floats by default even in the debugger.  It's a kluge because this
;;is not a print variable, but in fact this (documented but sbcl-specific" variable allows you
;;to bind anything.
(push '(*read-default-float-format* . double-float) *debug-print-variable-alist*)
;;(setq *print-array* nil)
(setq *debug-beginner-help-p* nil) ;shorter debugger messages

;;Make debugging easier
(import '(sb-debug:arg sb-debug:var sb-kernel:make-lisp-obj))

(proclaim '(sb-ext:disable-package-locks double-float single-float intersection))

;;Set up for GNU scientific library in Lisp.  Must be here to create packages.
;;sbcl must be at least version 1.2.5, and you must say module load libffi.

;;Load from each user's home directory.  This is probably a bad idea, but it is how quicklisp is set up
(load "~/quicklisp/setup")
(ql:quickload "gsll")

;;All traps are disabled by GSLL.  Why?  Put them back.
(sb-int:set-floating-point-modes :traps '(:overflow :divide-by-zero :invalid))

(defvar *recompile* nil)

;;Compile and load the system as needed.
(load-system *load-pathname* simulation-files *recompile*)

(let ((local-modifications (merge-pathnames "local-modifications.lisp" *load-pathname*)))
  (when (probe-file local-modifications)
    (with-open-file (stream local-modifications)
      (loop for first = t then nil
	    for form = (read stream nil)
	    while form
	    when first do (format t "~&******************** Local modifications ********************~%")
	    do (write form :level 3 :length 10) (terpri)
	    do (eval form)
	    finally (unless first (format t "~&******************** ******************* ********************~%"))
	    ))))
