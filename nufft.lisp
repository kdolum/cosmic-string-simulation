;;Fourier transform with non-uniform samples.  See Fourier.tex.

(in-package "CL-USER")

;;X is an array of positions x_i in [0,1]
;;F is an array of values f(x_i)
;;Computes F(j) = sum_i exp{-2 pi i j x_i} f(x_i)
;;Returns an array from REALFT, so element 0 is freq 0, element 1 is n/2, and
;;elements i and i+1 are the real and imaginary parts of frequency i
;;If the number of points in X and F is a power of two, this code does not do any oversampling, 
;;i.e., the oversampling factor alpha = 1.  As a result, the last small fraction of the returned
;;transform is wrong.  My understanding is that oversampling in this type of transform just means
;;computing more frequencies than you need.  Callers of this function generally use the first half only,
;;which corresponds to oversampling factor 2.  Probably this could be improved.
(defun nufft (x f m &optional n)
  (compute-b-spline-coefficients (* 2 m)) ;Set up
  (let* ((nn (length x))		;Number of samples N
	 (n (or n (expt 2 (1+ (ceiling (log nn 2)))))) ;Size of FFT: At least 2N, rounded up to power of 2
	 (g (make-array n :element-type 'double-float :initial-element 0.0))
	 (m2 (* 2 m)))			;Cache this
    (locally (declare (optimize speed)
		      (fixnum n nn m)
		      (type (signed-byte 32) m2) ;helps optimize expt
		      (type (simple-array double-float (*)) x f g))
      ;;Compute convolution g(l/n) for l=0...n-1
      (loop with nf = (double-float n)
	    for xi double-float across x ;I don't know why type declarations necessary
	    for fi double-float across f
	    for nxi double-float = (* nf xi)
	    do (locally (declare (optimize (safety 0)))	;Avoid check-fixnum calls
		 (loop for l fixnum from (- (fixnum-floor (+ nxi m))) repeat (* 2 m) ;2m points are nonzero
		       do (incf (aref g (mod l n)) (* fi (b-spline m2 (+ nxi l)))))))
      (realft g (/ n 2))		;Real Fourier transform of g
      (loop for i below n
	    for k double-float = (* pi (/ (floor i 2) (double-float n)))
	    ;;Divide by n and by Phi(k) = (1/n) sinc^{2m}(pi j/n) for j = 1..n to get f(j)
	    do (setf (aref g i)
		     (/ (aref g i)
			(expt (cond ((zerop i) 1.0) ;First elt is freq 0: sinc = 1
				    ((= i 1) (/ 2 pi)) ;This is frequency n/2, sinc = 2/pi
				    ;;Otherwise i and i+1 have freq j = i/2
				    (t (the (double-float 0.0) (/ (sin k) k))))
			      m2))))
      g)))				;Return g, now really F

    

;;compare NUFFT result with direct integration
(defun nufft-check (x f m &optional n)
  (let ((result (nufft x f m n)))
    (setq n (length result))		;See how many actually used in case defaulted
    ;;Return values have frequencies up through n/2, but some towards the end are not accurate
    (loop for j from 0 to (/ n 4)
	  ;;It's unclear if more than N (number of samples) should be expected to be right, but it seems to be OK
	  for fft = (cond ((zerop j) (aref result 0))   ;First elt is freq 0, real
;;			  ((= j (/ n 2)) (aref result 1)) ;Last is real and in a funny place 
			  (t (complex (aref result (* j 2)) (aref result (1+ (* j 2))))))
	  for ft = (loop for xi across x
			 for fi across f
			 sum (* fi (exp (* -2 pi (complex 0 1) j xi))))
	  for error = (abs (- ft fft))
;;	  do (format t "NUFFT: ~F; FT: ~F, error ~F~%" fft ft error)
	  maximize error into max-error
	  sum (abs ft) into total
	  finally (format t "~&Worst error ~G, total magnitude of true answers ~F~%" max-error total))
    result))


(defun check-nufft (m &key samples x f n print)
  (cond ((or x f)
	 (unless (and x f) (error "Giving only one of X and F does not make sense"))
	 (when samples (error "Giving sample number and supplying data doesn't make sense")))
	(t
	 (unless samples (error "Must give a number of samples to generate"))
	 (setq x (make-array samples :element-type 'double-float)
	       f (make-array samples :element-type 'double-float))
	 (dotimes (i samples)
	   (setf (aref x i) (random 1.0)
		 (aref f i) (- (random 1.0) 0.5)))))
  (when print
    (format t "X = ~S~%F = ~S~%" x f))
  (nufft-check x f m n))
