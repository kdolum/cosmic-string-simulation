;;Centered cardinal B-splines
(in-package "CL-USER")

;;See Fourier.tex.  The indices to the array are m, l if m is even or l-1/2 if odd, and k
(defvar b-spline-coefficients)
(declaim (type (simple-array double-float (* * *)) b-spline-coefficients))
(defvar b-spline-max-order nil)		;Order to which coefficients have been computed

(defun compute-b-spline-coefficients (max-order)
  (unless (and b-spline-max-order (>= b-spline-max-order max-order)) ;Unless already done
    (let ((coefficients (make-array (list (1+ max-order) (ceiling max-order 2) max-order))))
      (setf (aref coefficients 1 0 0) 1) ;Do computation in rationals
      (loop for m from 2 to max-order	 ;Compute order m
	    do (loop for l below (ceiling m 2)
		     do (loop for k below m
			      do (setf (aref coefficients m l k)
				       (compute-b-spline-coefficient coefficients m l k)))))
      (setq b-spline-max-order nil)	;Invalidate old array
      (setq b-spline-coefficients (make-array (array-dimensions coefficients) :element-type 'double-float))
      (dotimes (i (array-dimension coefficients 0)) ;Convert to double-floats for rapid computation later
	(dotimes (j (array-dimension coefficients 1))
	  (dotimes (k (array-dimension coefficients 2))
	    (setf (aref b-spline-coefficients i j k) (double-float (aref coefficients i j k) )))))
      (setq b-spline-max-order max-order)))) ;Say coefficients installed
		   
;;Compute coefficients at order m from coefficients at lower orders
(defun compute-b-spline-coefficient (coefficients m l k)
  ;;If m is even, then l represents l+1/2.  The same l for m-1 represents the left
  ;;previous coefficient.  But if m is odd, l is just l, and l-1 represents the previous
  ;;left coefficient.
  (let ((left-l (if (evenp m) l (1- l)))
	(right-l (cond ((= l (1- (ceiling m 2))) nil) ;No right predecessor
		       ((evenp m) (1+ l))
		       (t l))))
    (* (/ 1 (1- m))			;Common prefactor
       (loop for j from (max (1- k) 0) below (1- m)
	     for prev = (+ (if right-l (aref coefficients (1- m) right-l j) 0) ;Right if any
			   (if (minusp left-l)				       ;Need l=-1/2 coeff?
			       (* (expt -1 k) (aref coefficients (1- m) 0 j)) ;(-)^j times l=1/2
			     (* (expt -1 (- j k)) (aref coefficients (1- m) left-l j)))) ;Left
	     sum (* prev
		    (expt 2 (- k j 1))
		    (+ (if (< j k) 0	;No first term when j = k-1
			 (* m (binomial j k)))
		       (if (plusp k)	;No second term if k = 0
			   (binomial j (1- k))
			 0)))))))

(defun factorial (n)
  (assert (typep n '(integer 0 *)))
  (if (zerop n) 1
    (* n (factorial (1- n)))))

(defun binomial (n m)
  (/ (factorial n) (factorial m) (factorial (- n m))))

(declaim (inline b-spline))

;;Compute b-spline of order m.  Coefficients must be set up.
(defun b-spline (m x)
  (let* ((x (abs x))		       ;Functions are even
	 (l (if (evenp m) (fixnum-floor x)		       ;Find which piece to compute.  Even m: break at integers
	      (fixnum-floor (+ x 0.5)))))		       ;break at half-integers.  0...1/2 -> l=0
    (if (>= l (ceiling m 2))				       ;Too big?
	0.0
      (loop with result double-float = 0.0
	    for k from (1- m) downto 0
	    do (setq result (+ (* result x) (aref b-spline-coefficients m l k)))
	    finally (return result)))))
