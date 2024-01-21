;; Code to compute the gravitational radiation from the near cusp region.
(in-package "CL-USER")

;;;Numerical calculations for analytic approximation

(declaim (inline besselK13 besselK23 h-body))

;;Modified Bessel Function K_{1/3,x}
(defun besselK13 (x)
  (if (> x 100.0) 0.0			;Avoid underflow
    (besselk (/ 1.0 3.0) x)))

;;Modified Bessel Function K_{2/3,x}
(defun besselK23 (x)
  (if (> x 100.0) 0.0			;Avoid underflow
    (besselk (/ 2.0 3.0) x)))

;;The function H(a)
(defun h-body (z sa plusp)
  (declare (double-float z sa))
  (let ((zsmall (/ z sa))
	(zbig (* z sa)))
    (declare (double-float zbig zsmall))
    (* (expt z 2)
       (+ (* (+ (expt (besselK13 zbig) 2) (expt (besselK23 zbig) 2))
	     (+ (expt (besselK13 zsmall) 2) (expt (besselK23 zsmall) 2)))
	  (* 4 (if plusp 1 -1)
	     (besselK13 zbig)
	     (besselK13 zsmall)
	     (besselK23 zbig)
	     (besselK23 zsmall))))))

;;Explicit integration
(defun H2-direct (a y plusp  &key (low-limit 0.0))
  (if (zerop y) 0.0
    (let ((sa (if (< a 1) (sqrt (/ 1 a))
		(sqrt a))))
      (declare (double-float sa))
      (gsl:integration-qagp
       #'(lambda (z)
	   (locally			;Avoid warning about consing return value
	       (declare (optimize speed))
	     (h-body z sa plusp)))
       (grid:grid low-limit y)))))

(defun H-direct (a plusp &key (Bessel-function-argument-limit 20)) ;K(1/3,20) and K(2/3,20) are about 5e-10
  (H2-direct a (/ Bessel-function-argument-limit (sqrt a)) plusp))

;;Caching H(a)
(defparameter h-cache-a-max 1000.0)	;Maximum cached value
(defparameter h-cache-a-interval 0.1)	;Step size
(defparameter h-cache-a-count		;Count of steps.  Fudge to make sure max is covered.
  (1+ (ceiling (/ (- (* h-cache-a-max (+ 1 1e-15)) 1.0) h-cache-a-interval))))

(defun make-h-cache (plusp)
  (make-interpolation-1
   #'(lambda (a)
       (h-direct a plusp))
   1.0
   h-cache-a-interval
   h-cache-a-count))

(defun write-h-caches ()
  (write-interpolation-data (make-h-cache t) "Hplus.dat")
  (write-interpolation-data (make-h-cache nil) "Hminus.dat"))

(defvar *hplus-data* (read-interpolation-data (merge-pathnames "Hplus.dat" *load-pathname*)))
(defvar *hminus-data* (read-interpolation-data (merge-pathnames "Hminus.dat" *load-pathname*)))

(defun h-cached (a plusp)
  (interpolate (if plusp *hplus-data* *hminus-data*) a))

;; We defined the functions H+ and H- from Eq. 50 of the notes 
;; We split the definition in 2 parts the interpolated one
;; and the asymptotic form.
(defun hpm (a plusp)
  (when (< a 1) (setq a (/ 1 a)))
  (if (< a h-cache-a-max)			;Do by linear interpolation
      (h-cached a plusp)
    (asymptotic-h a plusp)))			;Large value: use asymptotic form

(defun asymptotic-h (a plusp)
  (+ (asymptotic-hcommon a) (* (if plusp 1 -1) (asymptotic-hdifference a))))

;;Parts of H(a) when a>>1.  This means that we're concerned only with z/sqrt{a} << 1, because otherwise we have
;;z sqrt{a} >> 1, giving exponential suppression.  By expanding the small-argument Bessel functions in power series
;;we get a definite integral we can do.
(defun asymptotic-hcommon (a)
  (declare (double-float a))
  (+ (/ (* 2 (expt 3.0 -1/2) (expt a -1/6) (expt pi 7/2))
	(gamma 1/6) (gamma 1/3))
     (/ (* 1/2 (expt 3.0 -1/2) (expt a -5/6) (expt pi 5/2) (gamma 7/6))
	(gamma 2/3))
     (* -1/8 (expt 3.0 1/2) (expt a -3/2) (expt pi 3))))

(defun asymptotic-hdifference (a)
  (declare (double-float a))
  (+ (* 1/3 (expt 3.0 -1/2) (expt a -1/2) (expt pi 3))
     (/ (* -1/2 (expt 3.0 1/2) (expt a -7/6) (expt pi 7/2))
	(gamma 1/6) (gamma 1/3))))

;;H(a,y), the integral only from 0 to y.

;;Small-argument approximations in both Bessel functions, valid when y sqrt{a} << 1.
;;Returns list of coefficients starting at y^(1/3) with steps of y^(2/3)
;;First omitted term y^(11/3).
;;For h23, power of y is higher by 1/3.  Thus if p=2/3, 4/3, 6/3,... is that new power,
;;the coefficient resulting from integration should be multiplied by (p-1/3)/p
(defun h2small-common-coefficients (a)
  (declare (double-float a))
  (list
   (* 3/2 (expt 2.0 -1/3) (expt (gamma 2/3) 4)) ;y^(1/3)
   (* 1/3 (+ 1 (expt a 2/3)) (expt a -1/3) (expt pi 2)) ;y
   (* 3/20 (expt 2 -2/3)				;y^(5/3) 
      (+ (expt (gamma 1/3) 4)
	 (* -6 (expt 3 1/2) (expt a -2/3) (+ 1 (expt a 4/3)) pi (expt (gamma 2/3) 2))))
   (* 1/56 (expt 2.0 -1/3) (expt a -1)	;y^(7/3)
      (+ (* -18 (expt 3.0 1/2) (+ 1 (expt a 2/3)) (expt a 2/3) pi (expt (gamma 1/3) 2))
	 (* 5 (+ 1 (expt a 2)) (expt (gamma -1/3) 2) (expt (gamma 2/3) 2))))
   (* 1/48 (expt a -4/3) (expt Pi 2)	;y^3
      (+ 7 (* 20 (expt a 2/3)) (* 108 (expt a 4/3)) (* 20 (expt a 2)) (* 7 (expt a 8/3))))))


;;First omitted term y^(11/3).  There must be the same number of terms here as above.
(defun h2small-difference-coefficients (a)
  (declare (double-float a))
  (list 0.0				;y^(1/3)
	(* 4/3 (expt Pi 2))		;y
	(* -3/5 (expt 2.0 1/3) (expt 3.0 1/2) (+ 1 (expt a 2/3)) (expt a -1/3) Pi (expt (Gamma 2/3) 2));y^(5/3)
	(*  (+ (* 1/42 (expt 2.0 -1/3) (expt (gamma -1/3) 4))					       ;y^(7/3)
	       (* 18/7 (expt a -2/3) (+ 1 (expt a 4/3)) (expt pi 5/2) (expt (gamma -1/6) -1))))
	(* 1/2 (expt a -1) (+ 2 (expt a 2/3) (expt a 4/3) (* 2 (expt a 2))) (expt Pi 2)))) ;y^3

;;h+(a,y) if plusp, h-(a,y) if not
(defun h2small (a y plusp)
  (loop for common-coefficient in (h2small-common-coefficients a)
	for difference-coefficient in (h2small-difference-coefficients a)
	for power-cubed from 1 by 2
	sum (* (expt y (/ power-cubed 3.0))
	       (if plusp (+ common-coefficient difference-coefficient)
		 (- common-coefficient difference-coefficient)))))

;;(1/theta) \int_0^theta dtheta' theta' h(a,y).
;;y ~ theta^3, so terms are multiplied by 1/2, 1/4, 1/6, ...
(defun h2small-integrated (a y plusp)
  (loop for common-coefficient in (h2small-common-coefficients a)
	for difference-coefficient in (h2small-difference-coefficients a)
	for power-cubed from 1 by 2
	sum (* (/ (expt y (/ power-cubed 3.0)) (1+ power-cubed))
	       (if plusp (+ common-coefficient difference-coefficient)
		 (- common-coefficient difference-coefficient)))))

;;Parts of H(a,y) when a>>1.  As for H(a) in this case, we can expand the small-argument Bessel functions in power
;;series, but now we get an indefinite integral we cannot do.  But we can cache the results
;;Common parts.
;;All these multiply the given power of a, which is is -(4*step+1)/6.
;;First omitted term a^-13/6.  Since we only use this for a>1000, this is plenty
(defun h2biga-common-integrand (step x)
  (* (ecase step
       (0 (* (expt 2.0 -2/3) (expt x 2/3) (expt (Gamma 2/3) 2)))     ;a^-1/6
       (1 (* 1/2 (expt 2.0 -1/3) (expt x 4/3) (expt (Gamma 1/3) 2)))) ;a^-5/6
     (+ (expt (BesselK13 x) 2) (expt (BesselK23 x) 2))))

(defun h2biga-common-integral (step limit)
  (if (zerop limit) 0.0
    (integrate #'(lambda (x) (h2biga-common-integrand step x))
	       0.0 limit)))

;;Difference parts.  Power of a is -(4*step+3)/6.  First omitted term a^-11/6
(defun h2biga-difference-integrand (step x)
  (* (ecase step
       (0 (* 2 (Gamma 1/3) (Gamma 2/3) x))			     ;a^-3/6
       (1 (* (expt 2.0 1/3) (expt x 5/3) (Gamma -1/3) (Gamma 2/3)))) ;a^-7/6
     (BesselK13 x) (BesselK23 x)))

(defun h2biga-difference-integral (step limit)
  (if (zerop limit) 0.0
    (integrate #'(lambda (x) (h2biga-difference-integrand step x))
	       0.0 limit)))

(defparameter h2biga-steps 2)

;;Range of y parameters that will be cached.  At a=1, y=4 is within one part in 10^7 of the y=infinity value
;;and for larger a it is even closer.
(defparameter h2-cache-y-max 4)
(defparameter h2-cache-y-steps 40)
(defun h2-cache-q-interval (max steps)		;Interval in q for caching
  (+ (/ (expt max (/ 1.0 3.0)) steps)
     1e-16))		;Fudge factor makes sure that y=h2-cache-y-steps is covered

;;We cache h(a,q^3) for q evenly spaced
(defun make-h2-cache (plusp)
  (make-interpolation
   #'(lambda (a q)
       (h2-direct a (expt q 3) plusp))
   '(1.0 0.0)				;a = 1... and q=0...
   (list h-cache-a-interval (h2-cache-q-interval h2-cache-y-max h2-cache-y-steps))
   (list h-cache-a-count
	 (1+ h2-cache-y-steps))))

(defun write-h2-caches ()
  (write-interpolation-data (make-h2-cache t) "H2plus.dat")
  (write-interpolation-data (make-h2-cache nil) "H2minus.dat"))

;;These cache the function f(a,q) = h(a,y=q^3)
(defvar *h2plus-data* (read-interpolation-data (merge-pathnames "H2plus.dat" *load-pathname*)))
(defvar *h2minus-data* (read-interpolation-data (merge-pathnames "H2minus.dat" *load-pathname*)))

(defun h2-cached (a y plusp)
  (interpolate (if plusp *h2plus-data* *h2minus-data*) a (expt y 1/3)))

;;We cache for y sqrt{a}=q^3, q evenly spaced
(defun make-h2biga-cache (step differencep)
  (make-interpolation-1
   #'(lambda (limit)			;y sqrt{a}
       (if differencep (h2biga-difference-integral step (expt limit 3))
	 (h2biga-common-integral step (expt limit 3))))
   0.0
   ;;Max y doubled to give max sqrt{a} y.  In h2 cache, we have (sqrt{a}+1/sqrt{a}) y,
   ;;which is always at least 2 y, so y can be less.
   (h2-cache-q-interval (* h2-cache-y-max 2) h2-cache-y-steps)
   (1+ h2-cache-y-steps)))

(defun write-h2biga-caches ()
  (loop for differencep in '(t nil)
	do (dotimes (step h2biga-steps)
	     (write-interpolation-data (make-h2biga-cache step differencep)
				       (format nil "H2biga~A~D.dat" (if differencep "difference" "common") step)))))

(defvar *h2bigacommon-data*
  (loop for step below h2biga-steps
	collect (read-interpolation-data (merge-pathnames (format nil "H2bigacommon~D.dat" step) *load-pathname*))))
(defvar *h2bigadifference-data*
  (loop for step below h2biga-steps
	collect (read-interpolation-data (merge-pathnames (format nil "H2bigadifference~D.dat" step) *load-pathname*))))

(defun h2biga (a y plusp)
  (loop for step below h2biga-steps
	for common in *h2bigacommon-data*
	for difference in *h2bigadifference-data*
	sum (* (expt a (- (/ (+ 1 (* step 4)) 6.0)))
	       (interpolate common (expt (* y (sqrt a)) (/ 1.0 3.0))))
	sum (* (if plusp 1 -1) (expt a (- (/ (+ 3 (* step 4)) 6.0)))
	       (interpolate difference (expt (* y (sqrt a)) (/ 1.0 3.0))))))

(defun h2pm (a y plusp)
  (when (< a 1) (setq a (/ 1 a)))
  (let ((sa (sqrt a)))
    (cond ((>= (* (+ sa (/ 1 sa)) y)	;The overall argument to exp is twice this.
	       (* 2 h2-cache-y-max))	;Out of cache for a=1.
	   (hpm a plusp))		;The fractional difference from y=infinity is less than 10^-7
	  ((< (* sa y) 0.02)		;Small zbig limit
	   (h2small a y plusp))
	  ((<= a h-cache-a-max)		;In cached range of y.  In range of a also?
	   (h2-cached a y plusp))	;Use cache
	  (t
	   (h2biga a y plusp))		;Large a approximation.
	  )))

;;This is like h except the integrand has an extra 1/3 power of z
(defun h3-body (z sa plusp)
  (* (expt z (/ 1.0 3.0)) (h-body z sa plusp)))

(defun h23-direct (a y plusp)
  (if (zerop y) 0.0
    (let ((sa (if (< a 1) (sqrt (/ 1 a))
		(sqrt a))))
      (declare (double-float sa))
      (gsl:integration-qagp
       #'(lambda (z)
	   (locally			;Avoid warning about consing return value
	       (declare (optimize speed))
	     (h3-body z sa plusp)))
       (grid:grid 0.0 y)))))

(defun H3-direct (a plusp &key (Bessel-function-argument-limit 20))
  (H23-direct a (/ Bessel-function-argument-limit (sqrt a)) plusp))

(defun make-h3-cache (plusp)
  (make-interpolation-1
   #'(lambda (a)
       (h3-direct a plusp))
   1.0
   h-cache-a-interval
   h-cache-a-count))

(defun write-h3-caches ()
  (write-interpolation-data (make-h3-cache t) "H3plus.dat")
  (write-interpolation-data (make-h3-cache nil) "H3minus.dat"))

(defvar *h3plus-data* (read-interpolation-data (merge-pathnames "H3plus.dat" *load-pathname*)))
(defvar *h3minus-data* (read-interpolation-data (merge-pathnames "H3minus.dat" *load-pathname*)))

(defun h3-cached (a plusp)
  (interpolate (if plusp *h3plus-data* *h3minus-data*) a))

(defun h3pm (a plusp)
  (when (< a 1) (setq a (/ 1 a)))
  (if (< a h-cache-a-max)			;Do by linear interpolation
      (h3-cached a plusp)
    (asymptotic-h3 a plusp)))			;Large value: use asymptotic form

(defun asymptotic-h3 (a plusp)
  (+ (asymptotic-h3common a) (* (if plusp 1 -1) (asymptotic-h3difference a))))

(defun asymptotic-h3common (a)
  (declare (double-float a))
   (+ (* 2/5 (expt a -1/3) (expt Pi 3/2) (Gamma 11/6)) 
      (* -6/5 (expt a -5/3) (expt Pi 3/2) (Gamma 11/6)) 
      (* 9/7 (expt 2 -1/3) (expt a -1) (expt (* 3 Pi) 1/2) (expt (Gamma 4/3) 2) (Gamma 13/6))
      (* 2/3 (expt a -7/3) (expt Pi 3/2) (Gamma 11/6)))) 

(defun asymptotic-h3difference (a)
  (declare (double-float a))
  (+ (* -1/3 (expt 2 1/3) (expt 3 -1/2) (expt a -2/3) Pi (Gamma -1/3) (Gamma 2/3)) 
     (* 2/9 (expt 2 1/3) (expt 3 -1/2) (expt a -4/3) Pi (Gamma -1/3) (Gamma 2/3))))

;;LIke h2small, but coefficients differ.  See h2small-common-coefficients.
(defun h23small (a y plusp)
  (loop for common-coefficient in (h2small-common-coefficients a)
	for difference-coefficient in (h2small-difference-coefficients a)
	for power-cubed from 2 by 2
	for power = (/ power-cubed 3.0)
	sum (* (expt y power)
	       (/ (- power 1/3) power)	;Correct h2 to h23
	       (if plusp (+ common-coefficient difference-coefficient)
		 (- common-coefficient difference-coefficient)))))

(defparameter h23-cache-y-steps 200)	;Need more steps because power only 1.5

(defun h23-cache-q-interval (max steps)		;Interval in q=y^(3/2) for caching
  (+ (/ (expt max (/ 2.0 3.0)) steps)
     1e-16))		;Fudge factor makes sure that y=h2-cache-y-steps is covered

;;We cache h(a,q^1.5) for q evenly spaced
(defun make-h23-cache (plusp)
  (make-interpolation
   #'(lambda (a q)
       (h23-direct a (expt q 1.5) plusp))
   '(1.0 0.0)				;a = 1... and q=0...
   (list h-cache-a-interval (h23-cache-q-interval h2-cache-y-max h23-cache-y-steps))
   (list h-cache-a-count
	 (1+ h23-cache-y-steps))))

(defun write-h23-caches ()
  (write-interpolation-data (make-h23-cache t) "H23plus.dat")
  (write-interpolation-data (make-h23-cache nil) "H23minus.dat"))

(defvar *h23plus-data* (read-interpolation-data (merge-pathnames "H23plus.dat" *load-pathname*)))
(defvar *h23minus-data* (read-interpolation-data (merge-pathnames "H23minus.dat" *load-pathname*)))

(defun h23-cached (a y plusp)
  (interpolate (if plusp *h23plus-data* *h23minus-data*) a (expt y 2/3)))

;;Like h2biga-common-integrand, but extra 1/3 power of x and -1/6 power of a.
(defun h23biga-common-integrand (step x)
  (* (expt x 1/3) (h2biga-common-integrand step x)))

(defun h23biga-common-integral (step limit)
  (if (zerop limit) 0.0
    (integrate #'(lambda (x) (h23biga-common-integrand step x))
	       0.0 limit)))

(defun h23biga-difference-integrand (step x)
  (* (expt x 1/3) (h2biga-difference-integrand step x)))

(defun h23biga-difference-integral (step limit)
  (if (zerop limit) 0.0
    (integrate #'(lambda (x) (h23biga-difference-integrand step x))
	       0.0 limit)))

;;We cache for y sqrt{a}=q^1.5, q evenly spaced
(defun make-h23biga-cache (step differencep)
  (make-interpolation-1
   #'(lambda (limit)			;y sqrt{a}
       (if differencep (h23biga-difference-integral step (expt limit 1.5))
	 (h23biga-common-integral step (expt limit 1.5))))
   0.0
   (h23-cache-q-interval (* h2-cache-y-max 2) h23-cache-y-steps)
   (1+ h23-cache-y-steps)))

(defun write-h23biga-caches ()
  (loop for differencep in '(t nil)
	do (dotimes (step h2biga-steps)
	     (write-interpolation-data (make-h23biga-cache step differencep)
				       (format nil "H23biga~A~D.dat" (if differencep "difference" "common") step)))))

(defvar *h23bigacommon-data*
  (loop for step below h2biga-steps
	collect (read-interpolation-data (merge-pathnames (format nil "H23bigacommon~D.dat" step) *load-pathname*))))
(defvar *h23bigadifference-data*
  (loop for step below h2biga-steps
	collect (read-interpolation-data (merge-pathnames (format nil "H23bigadifference~D.dat" step) *load-pathname*))))

(defun h23biga (a y plusp)
  (loop for step below h2biga-steps
	for common in *h23bigacommon-data*
	for difference in *h23bigadifference-data*
	sum (* (expt a (- (/ (+ 1 (* step 2)) 3.0)))
	       (interpolate common (expt (* y (sqrt a)) (/ 2.0 3.0))))
	sum (* (if plusp 1 -1) (expt a (- (/ (+ 2 (* step 2)) 3.0)))
	       (interpolate difference (expt (* y (sqrt a)) (/ 2.0 3.0))))))

(defparameter h-cache-a-smoothing-range 0) ;feature disabled

(defun h23pm (a y plusp)
  (when (< a 1.0) (setq a (/ 1 a)))
  (let ((sa (sqrt a)))
    (cond ((>= (* (+ sa (/ 1 sa)) y)	;The overall argument to exp is twice this.
	       (* 2 h2-cache-y-max))	;Out of cache for a=1.
	   (h3pm a plusp))		;The fractional difference from y=infinity is less than 10^-7
	  ((< (* sa y) 0.05)		;Small zbig limit.  More range in h23 than h2.
	   (h23small a y plusp))
	  ((<= a (- h-cache-a-max h-cache-a-smoothing-range))		;In cached range of y.  In range of a also?
	   (h23-cached a y plusp))	;Use cache
	  ;;Kluge: smooth out transition at edge of cached region to avoid problems in integration routines
	  ((< a h-cache-a-max)
	   (let ((cached (h23-cached a y plusp))
		 (approximation (h23biga a y plusp))
		 (parameter (/ (- a (- h-cache-a-max h-cache-a-smoothing-range))
			       h-cache-a-smoothing-range))) ;interpolation parameter 0...1
	     (+ (* cached (- 1 parameter)) (* approximation parameter))))
	  (t
	   (h23biga a y plusp))		;Large a approximation.
	  )))


;;Get information for cusp radiation calculations.
;;Returns values alpha-, alpha+, phi-, phi+, theta, sign((sin alpha+)(sin alpha-))
(defun cusp-observer-info (cusp-info observer-direction)
  (mirror-image-let*
      ((cusp-direction (cusp-info-direction cusp-info))
       (a-dd (cusp-info-a-pp cusp-info))	  ;Dimensionless second derivatives L a'', b''
       (a-dd-normalized (3vector-normalize a-dd)) ;Unit vectors in the directions of a'', b''
       (a-len (3vector-length a-dd))		  ;Magnitude of the vector a'', b''
       ;;The y-vector of the orthonormal basis whose z direction is the cusp
       (cusp-y-direction (3vector-normalize (3vector-cross-product cusp-direction observer-direction))) ;
       ;;The x-vector of the orthonormal basis whose z direction is the cusp
       (cusp-x-direction (3vector-cross-product cusp-y-direction cusp-direction))
       ;;Angle between a'' and x and b'' and x.  Between -pi and pi.  Thus sign(sin phi-a) = sign(phi-a)
       (phi-a (atan (3vector-dot a-dd-normalized cusp-y-direction) (3vector-dot a-dd-normalized cusp-x-direction)))
       (theta (spherical-angle observer-direction cusp-direction))) ;Angle between the cusp and observation directions
    (values a-len b-len phi-a phi-b theta (* (signum phi-a) (signum phi-b)))))

;;Angular power density in a single mode in units of G\mu^2.  The mode number does not need to be an integer.
;;Look at the notes on gw from cusps for the definition of the parameters.
(defun mode-angular-power-distribution-from-cusp (cusp-info observer-direction mode-number)
  (multiple-value-bind (aminus aplus phiminus phiplus theta sign)
      (cusp-observer-info cusp-info observer-direction)
    (let ((ximinus (/ (* 2 pi mode-number (expt theta 3.0) (expt (abs (sin phiminus)) 3.0)) (* 3 aminus)))
	  (xiplus (/ (* 2 pi mode-number (expt theta 3.0) (expt (abs (sin phiplus)) 3.0)) (* 3 aplus)))
	  (coefficient (/ (* 128 pi (expt mode-number 2) (expt theta 8) (expt (* (sin phiminus) (sin phiplus)) 4)) 
			  (* 9 (expt (* aminus aplus) 2)))))
    (* coefficient
     (+ (* (+ (expt (besselK13 xiplus) 2.0) (expt (besselK23 xiplus) 2.0))
	   (+ (expt (besselK13 ximinus) 2.0) (expt (besselK23 ximinus) 2.0)))
	(* 4.0 sign
	   (besselK13 xiplus)
	   (besselK13 ximinus)
	   (besselK23 xiplus)
	   (besselK23 ximinus)))))))

;;Spectral angular power distribution dP/(dOmega domega) for omega=4 pi n/L
(defun spectral-angular-power-distribution-from-cusp (cusp-info observer-direction loop-length mode-number)
  (/ (* loop-length (mode-angular-power-distribution-from-cusp cusp-info observer-direction mode-number))
     4 pi))


;; If we want to compare to the Burden loop this is an example.
;;(* 2.0 (spectral-angular-power-distribution-from-cusp (nth 0 cusps-from-loop) (3vector-normalize (make-3vector (* (cos 0.0) (sin 0.001)) (sin 0.0) (cos 0.001))) (* 2.0 pi) 10000.0))

;;(burden-radiation-analytic 1 1 10000 (/ pi 2.0) 0.001 0.0)
;; This calculates the angular power distribution for gravitational waves (dp/dOmega)*theta
;; from a cusp in G\mu^2 units. It is basically the integral over frequencies of the
;; spectral angular power distribution.
;; This depends only on the azimuthal angle phi and not on theta because we have (dp/dOmega)*theta
;; Look at the notes on gw from cusps for the definition of the parameters.
;;If n given, give only power above that harmonic
(defun angular-power-distribution-from-cusp (cusp-info observer-direction &optional n)
  (multiple-value-bind (aminus aplus phiminus phiplus theta)
      (cusp-observer-info cusp-info observer-direction)
    (angular-power-distribution-from-cusp-1 aminus aplus phiminus phiplus theta n)))

;;Cusp radiation power theta dP/dOmega (which doesn't depend on theta).
;;aminus and aplus are the magnitude of the second derivatives times L.
;;See Eq. 49 in the notes
;;If n given, give only power above that harmonic, in which case we need theta as well
(defun angular-power-distribution-from-cusp-1 (aminus aplus phiminus phiplus &optional theta n)
  (let* ((a (* (/ aminus aplus) (expt (abs (/ (sin phiplus) (sin phiminus))) 3)))
	 (coefficient (/ 48 (expt pi 2) (sqrt (abs (* aminus aplus (sin phiplus) (sin phiminus))))))
	 (plusp (plusp (* (sin phiplus) (sin phiminus))))
	 (integral (hpm a plusp))) ;Integral 0...infinity
;;    (format t "Coefficient ~S, integral ~S" coefficient integral)
    (when (and n (plusp n))		;Must remove first modes
      (let ((z (* 2/3 pi (sqrt (abs (/ (expt (* (sin phiplus) (sin phiminus)) 3) aminus aplus)))
		  (expt theta 3)
		  (+ n 0.5))))		;For power starting with harmonic n+1, start integral at n+0.5
	(let ((firstn (h2small a z plusp)))	;Remove power below threshold
;;	  (format t "decrement ~S~%" firstn)
	  (decf integral firstn))))
    (* coefficient integral)))

;;Returns the cusp approximation to (1/theta) \int_0^theta dtheta' theta' dP/dOmega from 0 to theta, considering
;;the first n harmonics only.  Including all harmonics, this integral would be theta dP/dOmega, which is what
;;is returned by angular-power-distribution-from-cusp-1.
;;This is not a good approximation.  The purpose is to subtract it off.
(defun integrated-angular-power-distribution-from-cusp-n (aminus aplus phiminus phiplus theta n)
  (let* ((a (* (/ aminus aplus) (expt (abs (/ (sin phiplus) (sin phiminus))) 3)))
	 (coefficient (* (/ 48 (expt pi 2) (sqrt (abs (* aminus aplus (sin phiplus) (sin phiminus)))))))
	 (plusp (plusp (* (sin phiplus) (sin phiminus))))
	 (z (* 2/3 pi (sqrt (abs (/ (expt (* (sin phiplus) (sin phiminus)) 3) aminus aplus)))
	       (expt theta 3)
	       (+ n 0.5))))		;For power starting with harmonic n+1, start integral at n+0.5
    (* coefficient (h2small-integrated a z plusp))))

;;First n only, for testing
(defun angular-power-distribution-from-cusp-n (cusp-info observer-direction n)
  (multiple-value-bind (aminus aplus phiminus phiplus theta plusp)
      (cusp-observer-info cusp-info observer-direction)
    (let ((a (* (/ aminus aplus) (expt (abs (/ (sin phiplus) (sin phiminus))) 3)))
	  (coefficient (* (/ 48 (expt pi 2) (sqrt (abs (* aminus aplus (sin phiplus) (sin phiminus))))))))
      (let ((z (* 2/3 pi (sqrt (abs (/ (expt (* (sin phiplus) (sin phiminus)) 3) aminus aplus)))
		  (expt theta 3)
		  (+ n 0.5))))		;For power up through harmonic n, end integral at n+0.5
	(* coefficient (h2small a z plusp))))))

;; We can also integrate the angular power distribution in the azimuthal angle to get
;; the total power from the cusp.
;; In reality this is \int{d\phi [\theta] \int{ dw {dP/dw dOmega}}.
;; We take some observer direction, but this is just to get the energy in a direction close to the cusp.
;; This is missing both a factor of 1/L and an overall factor \theta.
(defun total-power-from-cusp (cusp-info observer-direction)
  (multiple-value-bind
   (aminus aplus phiminus phiplus theta sign)
   (cusp-observer-info cusp-info observer-direction)
   (declare (ignore theta))
   (declare (ignore sign))
   (let ((coeff (* (/ 48.0 (expt pi 2)) (/ 1.0 (sqrt (* aminus aplus))))))
     ;;(a (* (/ aminus aplus) (expt (abs (/ (sin phiplus) (sin phiminus))) 3))))
     (flet ((f (x)
	       (* coeff (/ 1.0d0 (sqrt (abs (* (sin (- phiplus x)) (sin (- phiminus x))))))
		  (hpm (* (/ aminus aplus) (expt (abs (/ (sin (- phiplus x)) (sin (- phiminus x)))) 3))
		       (plusp (signum (* (sin (- phiplus x)) (sin (- phiminus x)))))))))
	   (qromb #'f 1.0d-12 (* 2.0d0 pi) :eps 1e-9)) ;; Here we perform the integration over the azimuthal angle.
     )))

;; behaves as TOTAL-POWER-FROM-CUSP, but assumes that DISCRETE-CUSP-INFO is a list of the form
;; ((x-a x-b) cusp-direction a'' b'')
;; which is then passed to GET-DISCRETE-CUSP-PARAMETERS in order to find the apm, phipm, etc.
;; L is the length of the loop the cusp resides upon.
;; :THETA is the angle (in radians) out to which we measure the cusp power. It should always be small, s.t.
;;     sin(THETA)~THETA, and includes modes out to n~1/THETA^3
(defun total-power-from-cusp-discrete (discrete-cusp-info L &key (theta 0.1))
  (let ((observer-direction ;we need a direction close to the cusp direction, but not exactly
	 (3vector-normalize (3vector+ (nth 1 discrete-cusp-info) 
				      (3vector-normalize (nth 2 discrete-cusp-info) 0.1)))))
    (multiple-value-bind
     (aminus aplus phiminus phiplus fake-theta sign)
     (get-discrete-cusp-parameters discrete-cusp-info observer-direction)
     (declare (ignore fake-theta))
     (declare (ignore sign))
     (let ((coeff (* theta (/ 48.0 (expt pi 2) L) (/ 1.0 (sqrt (* aminus aplus))))))
       ;;(a (* (/ aminus aplus) (expt (abs (/ (sin phiplus) (sin phiminus))) 3))))
       (flet ((f (x)
		 (* coeff (/ 1.0d0 (sqrt (abs (* (sin (- phiplus x)) (sin (- phiminus x))))))
		    (hpm (* (/ aminus aplus) (expt (abs (/ (sin (- phiplus x)) (sin (- phiminus x)))) 3))
			 (plusp (signum (* (sin (- phiplus x)) (sin (- phiminus x)))))))))
	     (qromb #'f 1.0d-12 (* 2.0d0 pi) :eps 1e-9)) ;; Here we perform the integration over the azimuthal angle.
       ))))

;; operates similarly to CUSP-OBSERVER-INFO with respect to what it returns, but instead takes as its argument a list of
;; the form
;; ((x-a x-b) cusp-direction a'' b'')
(defun get-discrete-cusp-parameters (discrete-cusp-info observer-direction)
  (let*
   ((cusp-direction (nth 1 discrete-cusp-info))
    (a-dd (nth 2 discrete-cusp-info))
    (a-dd-normalized (3vector-normalize a-dd))
    (b-dd (nth 3 discrete-cusp-info))
    (b-dd-normalized (3vector-normalize b-dd))
    (cusp-y-direction (3vector-normalize (3vector-cross-product cusp-direction observer-direction)))
    (cusp-x-direction (3vector-normalize (3vector-cross-product cusp-y-direction cusp-direction)))
    (phi-a (atan (3vector-dot a-dd-normalized cusp-y-direction) (3vector-dot a-dd-normalized cusp-x-direction)))
    (phi-b (atan (3vector-dot b-dd-normalized cusp-y-direction) (3vector-dot b-dd-normalized cusp-x-direction)))
    (theta (spherical-angle observer-direction cusp-direction)))
   (values (3vector-length a-dd) (3vector-length b-dd) phi-a phi-b theta (* (signum phi-a) (signum phi-b)))))

		       

;;********************************************************************************************************************************
;;********************************************************************************************************************************
;;; The following code is meant to test some of the functions above
;;********************************************************************************************************************************
;;********************************************************************************************************************************

#| 

(load "load")
;;Initialize for infinite-volume
(initialize :total-size nil )



(setq *loop-preservation-threshold* 0)
(setq *initial-time* 0.0)

;;Do not allow for intersections
(setq *intersection-probability* 0.0)

;;Must start late enough to not be in other regions
;; This code sets up the Burden loops family parametrizes by m, n and psi
(defun setup-burden-loop (&optional (n 1) (m 2) (psi (/ pi 2.0)) (points 512))
  (make-test-ab
   #'(lambda (sigma)
       (incf sigma (/ pi points m))	;Offset half step to avoid cusp at kink
       (make-3vector (/ (sin (* m (- 0.0 sigma))) m)
		     0.0
		     (/ (cos (* m (- 0.0 sigma))) m)))
   #'(lambda (sigma) 
       (incf sigma (/ pi points n))	;Offset half step to avoid cusp at kink
       (make-3vector (* (/ (sin (* n sigma)) n) (cos psi))
		     (* (/ (sin (* n sigma)) n) (sin psi))
		     (/ (cos (* n sigma)) n)
		     ))
   (* 2 pi)
   points
))

(setup-burden-loop 1 1 (/ pi 2.0))



;; Load the loop using the smoothing code with the lowest possible smoothing.
;; 1024 number of kinks
(let*((initial-smoothing-fraction  (nth 0 *smoothing-fraction-sequence*))
	(diamond-selected-loop (longest-loop))
	(rest-frame-loop-diamonds (rest-frame-diamonds diamond-selected-loop
						       (coerce (current-time) 'double-float) 
						       t)))
    (load-loop-in-simulation  rest-frame-loop-diamonds initial-smoothing-fraction 1024 0.0 200.0)
)

;;Evolve and check the properties of the loop
(evolve-loop-until-all-in-time) 
(print-loop-velocity)
(loop-sizes-kinks)
(evolve-loop-and-record-energy)


(setq cusps-from-loop (let*((diamond (longest-loop)))
			(multiple-value-bind (a-hats a-sigmas) (get-a-data-sigmas diamond)
			  (multiple-value-bind (b-hats b-sigmas) (get-b-data-sigmas diamond)
			    (let*((a-amplitudes (amplitudes-from-hats (lorentzian-smooth a-hats 0.005 :sigmas a-sigmas :n 512)))	
				  (b-amplitudes (amplitudes-from-hats (lorentzian-smooth b-hats 0.005 :sigmas b-sigmas :n 512))))
			      (find-smooth-cusps a-amplitudes b-amplitudes))))))
  

(* 2.0 (spectral-angular-power-distribution-from-cusp (nth 0 cusps-from-loop) (3vector-normalize (make-3vector (sin 0.001) 0.0 (cos 0.001))) (* 2.0 pi) 10000.0))

(angular-power-distribution-from-cusp (nth 0 cusps-from-loop) (3vector-normalize (make-3vector (* (sin 0.001) (cos 0.0)) (* (sin 0.001) (sin 0.0)) (cos 0.001))) (* 2.0 pi))

 
(total-power-from-cusp  (nth 0 cusps-from-loop) (3vector-normalize (make-3vector (* (sin 0.001) (cos 0.0)) (* (sin 0.001) (sin 0.0)) (cos 0.001))) (* 2.0 pi))

|#

