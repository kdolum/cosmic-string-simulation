;;;Gravitational radiation
;;Convention for function and variable names in this file:
;; wave = I or J (possibly rescaled)
;; ij = |I_perp|^2 and maybe Im(I_x I_y^*) and the same for J
;; spectrum = power spectrum (i.e., the angular power density for each harmonic) in a certain direction
;; power = the power in a given direction
;; total-spectrum, total-power = the above quantities integrated over solid angle
;; binned means the spectrum summed in a series of logarithmic bins in harmonic number
(in-package "CL-USER")

;;Check that sigmas are not erroneously dsigmas
(defun check-sigmas (sigmas)
  (unless (loop for index from 1 below (length sigmas)
		always (< (aref sigmas (1- index)) (aref sigmas index)))
    (error "sigmas are not in order, so probably they are dsigmas by mistake")))

;;Check that dsigmas are not erroneously sigmas
(defun check-dsigmas (dsigmas)
  (unless (loop for dsigma across dsigmas
		always (<= dsigma 1.0))
    (error "Some dsigmas are > 1, so probably they are sigmas by mistake")))

(defun total-sigma (dsigmas)
  (loop for dsigma across dsigmas
	sum dsigma))

;;Average of a and b.  Of course they should be the same...
(defun average-total-sigma (a-dsigmas b-dsigmas)
  (/ (+ (total-sigma a-dsigmas) (total-sigma b-dsigmas)) 2.0))

;;Compute gravitational radiation given arrays of a' (or b') and sigmas.
;;Return two arrays indexed by harmonic number k, with the k = 0 slot zero:
;; |I_perp^k|^2 and Im(I_x I_y^*), where x, y, and Omega form a righthand coordinate system and
;; I^(k)(Omega) = (1/L) int_0^L dv a'(v) exp(2 pi i k (v-a_z(v))/L)
;;See Fourier.tex 
;;If you want them computed in the rest frame, go to the rest frame first
;;Omega is a unit vector pointing toward the observer.  n is the number of modes that you want
(defun gravitational-wave (hats sigmas omega n)
  (check-sigmas sigmas)
  (multiple-value-bind (wavex wavey n) (gravitational-wave-1 hats sigmas omega n)
    (multiple-value-bind (iperp ixy) (process-gravitational-wave wavex wavey n)
      (values iperp ixy))))

;;Return array of loci = (sigma_i - Omega . a_i)/l
(defun gravitational-wave-locations (hats sigmas omega)
  (loop with n = (length sigmas)
	with l = (aref sigmas (1- n))
	with locations = (make-array n :element-type 'double-float)
	for index below n
	for pos = (make-zero-3vector) then (prog1 (3vector+ pos (3vector-scale hat (- next-sigma sigma)))
					     (deallocate 3vectors pos))
	for hat = (aref hats index)
	for sigma double-float = 0.0 then next-sigma
	for next-sigma double-float = (aref sigmas index)
	do (setf (aref locations index) (/ (- sigma (3vector-dot omega pos)) l))
	finally (return locations)))

;;Omega is the observation direction, direction1 and direction2 are the two perpendicular directions
;;Returns arrays of locations (phases) x_i and values f1_i and f2_i for non-uniform FFT
;; I_c =  1/(2 pi n i) a'_c / (1-Omega.a') [exp{-2 pi j (sigma_{i+1} - Omega.a_{i+1}) - exp{-2 pi j (sigma_i - Omega.a_i)}] for c=1,2
;;     =  1/(2 pi n i) sum_i fc_i exp{-2 pi i n x_i}
(defun gravitational-wave-data (hats sigmas omega direction1 direction2)
  (declare (simple-vector hats) (type (simple-array double-float (*)) sigmas))
  (let* ((nn (length sigmas))
         (loci (gravitational-wave-locations hats sigmas omega)) ;array of (sigma_i - Omega.a_i)
         (fvali1 (make-array nn :element-type 'double-float :initial-element 0.0))
         (fvali2 (make-array nn :element-type 'double-float :initial-element 0.0)))
    (declare (type (simple-array double-float (*)) loci fvali1 fvali2) (fixnum nn))
    ;;Now f_i
    (loop for index below nn
	  for next-index = (mod (1+ index) nn)
	  for hat = (aref hats index)
	  for denominator double-float = (/ (3vector-squared-length (3vector- hat omega)) 2) ;1-Omega.a'
	  for coefficient1 double-float = (/ (3vector-dot hat direction1) denominator)	     ;a'_1/(1-Omega.a')
	  for coefficient2 double-float = (/ (3vector-dot hat direction2) denominator)
	  do
	  (decf (aref fvali1 index) coefficient1)		 ;negatively in this slot
	  (incf (aref fvali1 next-index) coefficient1)		 ;positively in next slot
	  (decf (aref fvali2 index) coefficient2)
	  (incf (aref fvali2 next-index) coefficient2)
	  )
    (values loci fvali1 fvali2)))

;;Gather the data and do the Fourier transform and return (2 pi i n) I_x and similarly I_y, in the format returned by REALFT,
;;and the valid number of frequencies n
(defun gravitational-wave-1 (hats sigmas omega n)
  ;;Number of modes for NUFFT to compute.  It must be a power of 2.
  ;;Some toward the end will be wrong, so we compute twice as many as needed
  ;;Probably it could be a lot smaller.
  (let* ((nufft-modes (* n 2)))	;Minimum to compute
    (setq nufft-modes (expt 2 (ceiling (log nufft-modes 2))) ;Up to power of 2
	  n (/ nufft-modes 2))				     ;Number we will actually return
      (multiple-value-bind (unitx unity) (two-perpendicular-vectors omega)
	(multiple-value-bind (loc fvalx fvaly)
	    (gravitational-wave-data hats sigmas omega unitx unity)
	  (values (nufft loc fvalx *nufft-m* (* nufft-modes 2))
		  (nufft loc fvaly *nufft-m* (* nufft-modes 2))
		  n)))))

;;Process data from gravitational-wave-1, return |I_perp^k|^2 and Im(I_x I_y^*).  The argument n is the number of valid
;;frequencies.  The arrays may be larger.
(defun process-gravitational-wave (wavex wavey n)
  (declare (type (simple-array double-float (*)) wavex wavey))
  (let ((perp2spec (make-array n :element-type 'double-float))
	(xyspec (make-array n :element-type 'double-float)))
    (setf (aref perp2spec 0) 0.0
	  (aref xyspec 0) 0.0)
    (loop for j from 1 below n
	  for jj = (* 2 j)
	  do (setf (aref perp2spec j)	;find Fourier spectra for |I_perp|^2 (or J) 
		   (/ (+ (expt (aref wavex jj) 2) (expt (aref wavex (1+ jj)) 2)
			 (expt (aref wavey jj) 2) (expt (aref wavey (1+ jj)) 2))
		      4 pi pi j j)
		   ;;find Fourier spectra for Im(I_xI_y^*) (or J)
		   (aref xyspec j)
		   (/ (- (* (aref wavex jj) (aref wavey (1+ jj)))
			 (* (aref wavex (1+ jj)) (aref wavey jj))) 4 pi pi j j)))
    (values perp2spec xyspec)))

;;Compute one |I_perp|^2 or |J_perp|^2 by brute force using complex 3-vectors.
(defun slow-gravitational-wave (hats sigmas omega k)
  (check-sigmas sigmas)
  (let ((result (make-array 3 :initial-element 0.0)))
    (loop with count = (length sigmas)
	  with l = (aref sigmas (1- count)) ;Length of loop
	  with expt = (* (complex 0 -1) (/ (* 2 pi k) l))
	  for index below count
	  for sigma = 0.0 then next-sigma    ;Start of segment
	  for pos = zero-3vector then next-pos   ;Start of segment
	  for next-sigma = (aref sigmas index) ;end of segment
	  for hat = (aref hats index)      ;direction of segment
	  for next-pos = (3vector+ pos (3vector-scale hat (- next-sigma sigma))) ;end of segment
	  for scale = (/ (- (exp (* expt (- next-sigma (3vector-dot omega next-pos))))
			    (exp (* expt (- sigma (3vector-dot omega pos)))))
			 (- 1 (3vector-dot omega hat)))
	  do (dotimes (coordinate 3)
	       (incf (aref result coordinate) (* (aref hat coordinate) scale)))
	  finally
	  (unless (< (3vector-length next-pos) (* l 1d-10))
	    (warn "This loop does not appear to be in the rest frame.  Delta = ~F" (3vector-length next-pos)))
	  (return		;|I_perp|^2 = |I|^2 - |I.Omega|^2
	   (/ (- (loop for coordinate below 3
		       sum (expt (abs (aref result coordinate)) 2))
		 (expt (abs (loop for coordinate below 3
				  sum (* (aref result coordinate) (aref omega coordinate))))
		       2))
	      4 pi pi k k)))))

;;Number of bins for N harmonics
(defun bin-count (n bin-size) 
 (round (real-log n bin-size)))	;Full round so that we get remainder also

(declaim (inline bin-start-harmonic bin-end-harmonic))

;;First harmonic in bin.  If the bin edge falls between integers, take the next larger ones.  But if
;;it is a just a tiny amount above an integer, use that integer.  This avoids problems with, e.g,
;;(expt (expt 2.0 1/8) 8) = 2 + 4e-16.
(defun bin-start-harmonic (bin-size bin)
  (declare (optimize speed)
	   (type (double-float 0.0) bin-size)
	   (type (signed-byte 32) bin))	;Makes expt faster
  (let ((result (fixnum-ceiling (* (expt bin-size bin) (- 1 1e-14)))))
    (without-compiler-notes result)))	;Avoid warning about consing float when not inline

;;First harmonic beyond bin end
(defun bin-end-harmonic (bin-size bin)
  (bin-start-harmonic bin-size (1+ bin)))

;;Number of harmonics in bin
(defun bin-harmonic-count (bin-size bin)
  (- (bin-end-harmonic bin-size bin) (bin-start-harmonic bin-size bin)))

(declaim (inline harmonic-number-approximate sum-cos-n-small-x sum-cos-n-large-xn sum-cos-n))

(without-compiler-notes-return
;;Approximate H_n. = sum_1^n (1/n)  The error is less than 10^-8 for n>10.
(defun harmonic-number-approximate (n)
  (declare (optimize speed)
	   (type (unsigned-byte 62) n))
  (let* ((nf (double-float n))
	 (result (+ (log nf) euler-constant (/ 1 2 nf) (/ -1 12 (expt nf 2)) (/ 1 120 (expt nf 4)))))
    result)))

(without-compiler-notes-return
;;Calculate sum_{k=n_1}^infty cos(2 pi n x)/n in case where x is small using Euler-Maclaurin.  See Lerch.tex.
(defun sum-cos-n-small-x (x n)
  (declare (optimize speed)
	   (double-float x)
	   (fixnum n))
  (let* ((nn (double-float n))
	 (phase (* 2 pi x nn)))
    (+ (- (ci phase))				     ;cos integral
       (/ (+ (* (cos phase) (/ (+ 1 (/ 1 6 nn)) 2))  ;cos/2n + cos/12n^2 + sin/6n
	     (without-compiler-notes		     ;Avoid complaint involving reciprocal 6
	      (/ (* (sin phase) pi x) 6)))
	  nn)))))

(without-compiler-notes-return
;;Calculate sum_{k=n_1}^infty cos(2 pi n x)/n in case where x n is large using asymptotic
;;series for Lerch Phi in large n.  See Lerch.tex.
(defun sum-cos-n-large-xn (x n)
  (declare (optimize speed)
	   (double-float x)
	   (fixnum n))
  (let* ((nn (double-float n))
	 (phase1 (* 2 pi x))
	 (cos1 (cos phase1))
	 (cosn (cos (* phase1 nn)))
	 (cosn-1 (cos (* phase1 (1- nn))))
	 (cosn+1 (cos (* phase1 (1+ nn))))
	 (d (* 2 (- 1.0 cos1) nn)))
    (/ (- (- cosn cosn-1)
	  (/ (- (+ cosn+1 cosn-1) (* 2 cosn)) d)
	  )
       d))))

;;Return sum_{k=n_1}^infty cos(2 pi n x)/n
;;except that if x=0 instead of the divergent sum_{k=n_1}^infty 1/n, return -sum_{k=1}^{1-n} 1/n
;;Thus (- (sum-cos-n-0 n1) (sum-cos-n-0 n2)) = sum_{k=n_1}^{n_2-1} 1/n regardless
;;See Lerch.tex
(defun sum-cos-n (x n)
  (declare (type fixnum-size-float x) (fixnum n))
  (setq x (mod x 1))			;Mod out 2pi in phase
  (when (> x 0.5) (setq x (- 1.0 x)))	;-x and x are the same.
  (cond ((zerop x)
	 (- (harmonic-number-approximate (1- n))))
	((< (* n (expt x 3)) 0.18)
	 (sum-cos-n-small-x x n))
	(t (sum-cos-n-large-xn x n))))

(defun test-sum-cos-n (x n1 n2)
  (let ((result (- (sum-cos-n x n1) (sum-cos-n x n2)))
	(direct (loop for n from n1 below n2
		      sum (/ (cos (* 2 pi n x)) n))))
    (format t "sum-cos-n ~S ~D ~D: our code: ~S direct sum: ~S, diff ~S~%" x n1 n2 result direct (- result direct))
    result))
		 
;;Increasing this above 1 causes more pairs of a'/b' values to be included when computing high frequencies.
(defparameter sum-cos-n-threshold-factor 1)

;;Suppose n >> 1, n even, and x = 0.5.  Then (sum-cos-n x n infinity) = 1/n - 1/(n+1) + 1/(n+2) - ... 
;;= 1/(n(n+1)) + 1/((n+2)(n+3) + ... = int_n^infty dk 1/(2k^2) = 1/(2n) approximately.
;;For x<<1, the individual cosines are close together, so we can approximate by int_n^infty dk cos(2 pi n x)/k
;;= -Ci(2 pi n x), whose magnitude is at most 1/(2 pi n x).  The values of sum-cos-n oscillate, but the envelope
;;gradually declines.  If there are N elements, and there will be N^2 contributions to the power, but their
;;signs are random, so we expect about N/(2 pi n x) in total.  Thus we will not consider contributions from
;;x > N/(2 pi n), because their total contribution will not be significant.
(defun sum-cos-n-threshold (element-count n)
  (min 0.5 (/ (* element-count sum-cos-n-threshold-factor) (* 2 pi n))))

;;Compute spectrum of a loop in a given direction
(defun binned-gravitational-spectrum (a-hats a-sigmas b-hats b-sigmas omega direct-n &key (total-n direct-n) (bin-size 2.0))
  (mirror-images (check-sigmas a-sigmas))
  (multiple-value-bind (bins error) (bin-count total-n bin-size) ;Total number of bins
    (unless (< (abs error) 1e-12) (error "Bin size ~S did not go evenly into total frequencies ~D" bin-size total-n))
    (multiple-value-bind (iperp ixy) (gravitational-wave a-hats a-sigmas omega direct-n)
      (setq direct-n (length iperp))	;Actual number directly computed
      (multiple-value-bind (direct-bins error) (bin-count direct-n bin-size) ;Number of bins directly computed
	(unless (< (abs error) 1e-12) (error "Bin size ~S did not go evenly into directly computed frequencies ~D"
				     bin-size direct-n))
	(multiple-value-bind (jperp jxy) (gravitational-wave b-hats b-sigmas omega direct-n)
	  (let* ((spectrum (spectrum-from-ij iperp ixy jperp jxy)) ;spectrum through direct-n
		 (power (bin-gravitational-spectrum spectrum total-n bin-size))) ;Bin low-f data
	    (when (> total-n direct-n)						 ;Extension?
	      (extend-gravitational-spectrum power a-hats a-sigmas b-hats b-sigmas omega direct-bins bins bin-size))
	    power))))))

(defvar *bpd-split-factor* 4)

(defun extend-gravitational-spectrum (power a-hats a-sigmas b-hats b-sigmas omega direct-bins bins bin-size
					    &optional (split-factor *bpd-split-factor*))
  ;;Split up bins for better accuracy
  (let* ((first-bin-split (* direct-bins split-factor))
	 (bins-split (* bins split-factor))
	 (bin-size-split (expt bin-size (/ 1 (double-float split-factor))))
	 (iperp-extension (binned-perp-direct a-hats a-sigmas omega first-bin-split bins-split bin-size-split))
	 (jperp-extension (binned-perp-direct b-hats b-sigmas omega first-bin-split bins-split bin-size-split)))
    ;;Extend power spectrum using I and J computed individually from a' and b'
    ;;iperp-extension has the sum of n |I_perp|^2.  To get the average we should divide by
    ;;of the number of harmonics in the bin.  Then we multiply by jperp-extension to get the
    ;;average of n^2 |I_perp|^2 |J_perp|^2.  To get the sum we multiply by the number of harmonics,
    ;;but only once, so we must divide once.
    (loop for bin from direct-bins below bins
	  do (setf (aref power bin)
		   (loop for split below split-factor
			 for sub = (+ (* bin split-factor) split)
			 ;;Add up results from sub-bins, processed separately
			 sum (/ (* 8 pi (* (aref iperp-extension sub) (aref jperp-extension sub)))
				(bin-harmonic-count bin-size-split sub)))))))

(defun test-binned-gravitational-spectrum (a-hats a-sigmas b-hats b-sigmas omega direct-n
						       &key  (total-n direct-n) (bin-size 2.0))
  (let ((start-time (get-internal-run-time))
	(power (binned-gravitational-spectrum a-hats a-sigmas b-hats b-sigmas omega direct-n
						   :total-n total-n :bin-size bin-size))
	(extended-time (get-internal-run-time))
	(direct (binned-gravitational-spectrum a-hats a-sigmas b-hats b-sigmas omega total-n :bin-size bin-size))
	(end-time (get-internal-run-time)))
    (format t "~D a', ~D b', ~S seconds extended, ~S seconds direct" (length a-sigmas) (length b-sigmas)
	    (/ (double-float (- extended-time start-time)) internal-time-units-per-second)
	    (/ (double-float (- end-time extended-time)) internal-time-units-per-second))
    (gnuplot 3 (length power)
	     #'(lambda (plot point)
		 (if (eq point :title) (nth plot '("directly computed" "extended" "abs error"))
		   (ecase plot
		     (0 (values (* (log bin-size 2.0) point) (aref direct point)))
		     (1 (values (* (log bin-size 2.0) point) (aref power point)))
		     (2 (and (> (expt bin-size point) direct-n)
			     (values (* (log bin-size 2.0) point)
				     (abs (- (aref power point) (aref direct point)))))))))
	     :logscale :y :styles :lines
	     :prelude (format nil "set format '%h'~%")
	     )))

;;Return an array of (phase coefficient1 coefficient2)
;;I^k = sum (coefficient1 direction1 + coefficient2 direction2) / (2 pi i k) exp(-2 pi i k phase)
;;returned phases are in order.
(defun gravitational-wave-elements (hats sigmas omega direction1 direction2)
  (loop with n = (length sigmas)
	with result = (make-array (* n 2))
	with locations = (gravitational-wave-locations hats sigmas omega)
	for index below n
	for next-index = (1+ index)
	for hat = (aref hats index)
	for denominator double-float = (/ (3vector-squared-length (3vector- hat omega)) 2)
	for coefficient1 double-float = (/ (3vector-dot hat direction1) denominator)
	for coefficient2 double-float = (/ (3vector-dot hat direction2) denominator)
	do (setf (aref result (1+ (* index 2)))	;Slot 1 has first loc when it is this phase.  N-1 has last loc
		 (list (aref locations index)	;phase at start of segment
		       (- coefficient1) (- coefficient2))) ;negatively for this phase
	if (< next-index n)				   ;Normal case
	do (setf (aref result (+ 2 (* index 2)))	   ;Even slots have next phase data
		  (list (aref locations next-index) ;phase at end of segment
			coefficient1 coefficient2)) ;positively for this phase
	else do (setf (aref result 0)		    ;Slot 0 has next-phase data for last location
		      (list (aref locations 0)
			    coefficient1 coefficient2))
	finally (return result)))

;;Account spectral contribution of a single pair of I_perp elements
;;Returns T if we did anything, NIL if first threshold was not satisfied.
(defun binned-perp-direct-1 (loc fvalx fvaly index1 index2 thresholds result first-bin bins bin-size wrap)
  (declare (optimize speed)
	   (type (simple-array double-float (*)) loc fvalx fvaly thresholds result)
	   (float bin-size)
	   (fixnum index1 index2 first-bin bins))
  (let ((phase (- (aref loc index2) (aref loc index1)))) ;relative phase
    (when wrap					 ;We have wrapped.  x2 < x1
      (setq phase (+ phase 1)))			 ;mod into 0...1
    (loop for bin from first-bin below bins
	  for threshold double-float = (aref thresholds bin) ;Max value of phase to use in this bin
	  until (and (>= phase threshold)		;Stop if phase is larger than threshold
		     (or wrap (> phase threshold)))	;Probably silly, but prevents duplication if phase=0.5 exactly
	  for this-sum-cos double-float = (sum-cos-n phase (bin-start-harmonic bin-size first-bin)) then next-sum-cos ;Sum from start of this bin
	  for next-sum-cos double-float = (sum-cos-n phase (bin-end-harmonic bin-size bin)) ;Sum from end of this bin
	  do (incf (aref result bin) (* 2 (+ (* (aref fvalx index1) (aref fvalx index2)) (* (aref fvaly index1) (aref fvaly index2)))
					(- this-sum-cos next-sum-cos)))
	  ;;If we exited on the first test of threshold, then bin = first-bin.  Otherwise it has been incremented (even if there's only one bin)
	  finally (return (> bin first-bin)))))

;;Special case for element interacting with itself
(defun binned-perp-direct-1-self (fvalx fvaly index result first-bin bins bin-size)
  (declare (optimize speed)
	   (type (simple-array double-float (*)) fvalx fvaly result)
	   (float bin-size)
	   (fixnum index first-bin bins))
  (loop for bin from first-bin below bins
	for this-sum-cos double-float = (sum-cos-n 0.0 (bin-start-harmonic bin-size first-bin)) then next-sum-cos ;Sum from start of this bin
	for next-sum-cos double-float = (sum-cos-n 0.0 (bin-end-harmonic bin-size bin)) ;Sum from end of this bin
	do (incf (aref result bin) (* (+ (expt (aref fvalx index) 2) (expt (aref fvaly index) 2))
				      (- this-sum-cos next-sum-cos)))
	))

;;Return the sum of n |I_perp|^2 (or n |J_perp|^2) in each bin
(defun binned-perp-direct (hats sigmas omega first-bin bins bin-size)
  (multiple-value-bind (direction1 direction2) (two-perpendicular-vectors omega)
    (multiple-value-bind (loc fvalx fvaly)
	(gravitational-wave-data hats sigmas omega direction1 direction2)
      (let ((result (make-array bins :element-type 'double-float))
	    (thresholds (make-array bins :element-type 'double-float))
	    (count (length loc)))
	(loop for bin from first-bin below bins
	      do (setf (aref thresholds bin) (sum-cos-n-threshold count (bin-start-harmonic bin-size bin))))
	(loop for index1 below (length loc)
	      do (binned-perp-direct-1-self fvalx fvaly index1 result first-bin bins bin-size) ;Self term
	      do (loop for index2 fixnum from (1+ index1) below (length loc) ;Start with next element
		       ;;Account this pair of elements.  If phase the difference too large for even first threshold, stop scanning.
		       while (binned-perp-direct-1 loc fvalx fvaly index1 index2 thresholds result first-bin bins bin-size nil))
	      ;;Now do elements that have wrapped, i.e., they are still after our location (mod 1) by less than threshold.
	      do (loop for index2 fixnum from 0
		       while (binned-perp-direct-1 loc fvalx fvaly index1 index2 thresholds result first-bin bins bin-size t)))
	(loop for bin from first-bin below bins
	      when (minusp (aref result bin)) do (warn "Negative result in ~S" 'binned-perp-direct)
	      do (setf (aref result bin) (/ (aref result bin) 4 pi pi)))
	result))))



;;Compute spectrum from I_perp, Im(I_x I_y^*), and the same for J.
(defun spectrum-from-ij (iperp ixy jperp jxy)
  ;;To avoid aliasing, it may be necessary to compute more frequencies in, say, A, than will ever be used because
  ;;they are strongly suppressed in B.
  (let* ((frequencies (min (length iperp) (length jperp)))
	 (spectrum (make-array frequencies :element-type 'double-float)))
    (setf (aref spectrum 0) 0.0)
    (loop for harmonic from 1 below frequencies
	  do (setf (aref spectrum harmonic)
		   (* 8 pi (expt harmonic 2)
		      (+ (* (aref iperp harmonic) (aref jperp harmonic))
			 (* 4 (aref ixy harmonic) (aref jxy harmonic))
			 ))))
    spectrum))

;;Separate directly computed power into bins.  Bin-size must go evenly into n.
(defun bin-gravitational-spectrum (power n bin-size)
  (declare (type (simple-array double-float (*)) power)
	   (double-float bin-size)
	   (fixnum n))
  (let* ((count (length power))
	 (bins (round (bin-count n bin-size)))
	 (data (make-array bins :element-type 'double-float :initial-element 0.0)))
    (locally (declare (optimize speed))
      (loop for i from 1 below count
	    do (incf (aref data (fixnum-floor (/ (real-log (double-float i)) (real-log bin-size))))
		     (aref power i))))
    data))

;;Power in gravitational waves from a loop in a given direction, in units of G mu^2
(defun gravitational-power (a-hats a-sigmas b-hats b-sigmas omega direct-n &key (total-n direct-n) (bin-size 2.0))
  (let ((spectrum (binned-gravitational-spectrum a-hats a-sigmas b-hats b-sigmas omega direct-n :total-n total-n :bin-size bin-size)))
    (loop for k below (length spectrum)
          sum (aref spectrum k))))

#|
;;Wrong calculation!
(defun slow-gravitational-wave-power (diamond omega k)
  (multiple-value-bind (a-hats a-sigmas) (get-a-data-sigmas diamond)
    (multiple-value-bind (b-hats b-sigmas) (get-b-data-dsigmas diamond)
      (let ((i (slow-gravitational-wave-a a-hats a-sigmas omega k)) ;I_perp^2
	    (j (slow-gravitational-wave-b b-hats b-sigmas omega k))) ;J_perp^2
	(* (* 8 pi) k k i j)))))
|#

;;Total gravitational power emitted by a loop
(defun total-gravitational-power (a-hats a-sigmas b-hats b-sigmas direct-n &key (total-n direct-n) (bin-size 2.0) (split-levels 0))
  (spherical-integral 
   #'(lambda (omega) (gravitational-power a-hats a-sigmas b-hats b-sigmas omega direct-n :total-n total-n :bin-size bin-size))
   split-levels))

;;Returns total power and momentum in gravitational waves as a 4vector
(defun gravitational-momentum-4vector (a-hats a-sigmas b-hats b-sigmas direct-n
					     &key (total-n direct-n) (bin-size 2.0) (split-levels 0))
  (let ((result (make-zero-4vector)))
    (flet ((integrate-power (omega area)
	     (let ((power (* area ;Power in this direction
			     (gravitational-power a-hats a-sigmas b-hats b-sigmas omega
						  direct-n :total-n total-n :bin-size bin-size)))) 
	       (4vector-incf result (4vector-scale (3to4vector omega 1.0) power)))))
      (spherical-integral-1 #'integrate-power split-levels)
      result)))

;;Backward compatibility: energy and magnitude of momentum as separate values
(defun gravitational-energy-momentum (&rest args)
  (let ((momentum (apply #'gravitational-momentum-4vector args)))
    (values (4vector-t momentum) (3vector-length momentum))))

;;Binned 4-momentum spectrum integrated over directions 
(defun total-gravitational-4vector-spectrum (a-hats a-sigmas b-hats b-sigmas direct-n
						 &key (total-n direct-n) (bin-size 2.0) (split-levels 0))
  (let* ((bins (bin-count total-n bin-size))
	 (momentum (coerce (loop repeat bins collect (make-zero-4vector)) 'vector)))
    (flet ((integrate-momentum (omega area)
	     (loop with spectrum = (binned-gravitational-spectrum a-hats a-sigmas b-hats b-sigmas omega direct-n
								  :total-n total-n :bin-size bin-size)
		   for k below bins
		   do (4vector-incf (aref momentum k) 
				    (4vector-scale (3to4vector omega 1.0)
						   (* (aref spectrum k) area))))))
      (spherical-integral-1 #'integrate-momentum split-levels)
      momentum)))
				  
;;Just energy
(defun total-gravitational-spectrum (&rest args)
  (map 'vector #'(lambda (x) (4vector-t x)) (apply #'total-gravitational-4vector-spectrum args)))


;;Spherical integration

;;Do a spherical integral by averaging faces in in a triangulation.
(defun spherical-integral (function split-levels)
  (let ((result 0.0))
    (spherical-integral-1 #'(lambda (center area) ;Add up area times function value at center
			      (incf result (* area (funcall function center))))
			  split-levels)
    result))

;;Call function with center and area for each triangle
(defun spherical-integral-1 (function split-levels &optional (triangles (triangulate-sphere split-levels)))
  (format t "Total of ~D triangles: " (length triangles)) (force-output)
  (loop for triangle across triangles
	for count from 1
	for center = (triangle-center triangle)
	for area = (apply #'spherical-triangle-area triangle)
	do (funcall function center area)
	when (zerop (mod count 10)) do (format t "~D " count) (force-output))
  (terpri))

;;Adaptive spherical integral, following Boal and Sayas, "Adaptive Numerical Integration on Spherical Triangles"

(defstruct (spherical-integration-element
	    (:conc-name SI-))
  a b c					;Outer triangle
  ab bc ca				;Midpoints of segments
  i0					;Integral over outer triangle
  ia ib ic icenter			;Integrals over 4 subtriangles
  level					;Split level of smaller triangles
  cost					;Cost of doing the 4 integrations
  )

;;Estimated (absolute) error due to this triangle is 4/3 (see Boal and Sayas) times the difference between using the 4
;;subtriangles and using only the full triangle.
;;We could make an approximation by just comparing the outer subtriangles to the center, but when we split we already have the
;;outer triangle integral, so we don't need to recompute it.
(defun si-error (si)
  (abs (* 4/3 (- (si-value si) (si-i0 si)))))

(defun si-error-larger-p (si1 si2)
  (> (si-error si1) (si-error si2)))

;;See if si1 is the better candidate to split than si2
(defun si-better-split-p (si1 si2)
  (> (/ (si-error si1) (si-cost si1))
     (/ (si-error si2) (si-cost si2))))

;;Contribution to integral from this element
(defun si-value (si)
  (+ (si-ia si) (si-ib si) (si-ic si) (si-icenter si)))
 
(defvar my-heap)
(defvar my-si)

(defun adaptive-spherical-integral (function initial-split-levels tolerance max-recursion)
  (let* ((heap (setup-initial-si function initial-split-levels))
	 (integral 0.0)
	 (error 0.0)
	 (cost 0.0))
    (setq my-heap heap)
    (heap-map heap #'(lambda (si)
		       (incf integral (si-value si))
		       (incf error (si-error si))
		       (incf cost (si-cost si))))
    (loop with last-count = 0
	  with si
	  when (>= (heap-count heap) (* 2 last-count))
	  do (format t "~&Error ~S with ~D samples and cost ~S~%" error (* 4 (setq last-count (heap-count heap))) cost)
	  when (<=  error (* tolerance (abs integral))) return t
	  do (setq si (heap-remove heap))	;Triangle with largest error
	  when (= (si-level si) max-recursion)
	    do (setq my-si si)
	    and do (error "Max recursion ~D reached without achieving requested tolerance ~S.  Value ~S, error estimate ~S"
		      max-recursion tolerance integral error)
	  do (decf integral (si-value si)) ;Remove from result
	  do (decf error (si-error si))
	  do (decf cost (si-cost si))
	  do (dolist (sub (si-split si heap function))	;Split into 4 SI's
		 (incf integral (si-value sub))	;Add them in to  result
		 (incf error (si-error sub))
		 (incf cost (si-cost sub))))
    (format t "Integration succeeded with ~D samples, error ~S, and cost ~S" (* 4 (heap-count heap)) error cost)
    integral))
	  

(defun setup-initial-si (function initial-split-levels)
  (setq initial-split-levels (max 0 (1- initial-split-levels)))
  (format t "~&Starting with ~D samples" (* 20 (expt 4 (1+ initial-split-levels))))
  (loop with heap = (create-heap #'si-better-split-p)
	for (a b c) across (triangulate-sphere initial-split-levels)
	do (create-si heap function a b c)
	finally (format t "initial scan done ~%")
	(return heap)
	))

;;Integrate one triangle by taking the midpoint.
(defun spherical-triangle-integrate (function a b c)
  (multiple-value-bind (value cost) (funcall function (spherical-triangle-center a b c))
    (unless cost (setq cost 1.0))
    (values (* value (spherical-triangle-area a b c)) cost)))

;;Split SI object and put subobjects into the heap.  Return list of them
(defun si-split (si heap function)
  (list (create-si heap function (si-a si) (si-ab si) (si-ca si) (1+ (si-level si)) (si-ia si))
	(create-si heap function (si-b si) (si-bc si) (si-ab si) (1+ (si-level si)) (si-ib si))
	(create-si heap function (si-c si) (si-ca si) (si-bc si) (1+ (si-level si)) (si-ic si))
	(create-si heap function (si-ab si) (si-bc si) (si-ca si) (1+ (si-level si)) (si-icenter si))))

;;Create an SI object, filling in the mid points and the sub-integrals.  Put it in the heap.
(defun create-si (heap function a b c &optional (level 1) i0)
  (let ((total-cost 0.0))
    (flet ((integrate-with-cost (a b c)
	     (multiple-value-bind (value cost) (spherical-triangle-integrate function a b c)
	       (incf total-cost cost)
	       value)))
      (let* ((ab (spherical-midpoint a b))
	     (bc (spherical-midpoint b c))
	     (ca (spherical-midpoint c a))
	     (si (make-spherical-integration-element
		  :a a :b b :c c
		  :ab ab :bc bc :ca ca
		  :ia (integrate-with-cost a ab ca)
		  :ib (integrate-with-cost b bc ab)
		  :ic (integrate-with-cost c ca bc)
		  :icenter (integrate-with-cost ab bc ca)
		  :level level
		  :i0 (or i0 (spherical-triangle-integrate function a b c))) ;i0 not included in cost
		 ))
	(setf (si-cost si) total-cost)
	(heap-insert heap si)))))

;;Triangulate sphere and split as in triangulate-sphere-separate-cusps
;;Call FUNCTION with triangle and a list of cusp-info
;;You have to add up the results yourself
(defun spherical-integral-separate-cusps (function split-levels cusp-data threshold1 threshold2)
  (let ((data (triangulate-sphere-separate-cusps split-levels cusp-data threshold1 threshold2)))
    (loop for (triangle . list) across data
	  count list into cusp-triangles
	  sum (length list) into total
	  finally (format t "Cusps: ~D, Total triangles: ~D, triangles near cusps: ~D~@[, cusps/triangle ~S~]~%"
			  (length cusp-data) (length data) cusp-triangles
			  (and (plusp cusp-triangles) (/ (double-float total) cusp-triangles))))
    (loop for (triangle . list) across data
	  for count from 1
	  do (funcall function triangle list)
	  when (zerop (mod count 10)) do (format t "~D " count) (force-output))))

;;Plots a set of triangles
;;If separate given, push each triangle in the direction of its center by this amount to separate them.
(defun plot-triangles (faces  &rest gnuplot-keys &key (separate 0.0))
  (apply #'gnuplot 1 (* (length faces) 6)
	 #'(lambda (plot point)
	     (declare (ignore plot))
	     (unless (eq point :title)
	       (multiple-value-bind (triangle corner) (floor point 6)
		 (and (< corner 4)	;double NIL to break between surfaces
		      (let ((face (aref faces triangle)))
			(values-list (3vector-list
				      (3vector+ (nth (if (= corner 3) 0 corner) face)
						(3vector-scale (triangle-center face) separate)))))))))
						   

	 :3d t
	 :styles :lines
	 gnuplot-keys))

;;Plot triangles showing which involve cusps
(defun test-triangulate-sphere-separate-cusps (levels cusp-data threshold1 threshold2
						      &rest gnuplot-keys &key no-far &allow-other-keys)
  (let ((triangles (triangulate-sphere-separate-cusps levels cusp-data threshold1 threshold2)))
    (apply #'gnuplot (1+ (* 2 (length cusp-data))) ;0 = far, 1 = first circle, 2 = first triangles, ...
	     (max 200 (* 6 (length triangles)))
	     #'(lambda (plot point)
		 (let ((direction (and (plusp plot) (cusp-info-direction (nth (floor (1- plot) 2) cusp-data)))))
		   (if (eq point :title)
		       (if (zerop plot) (unless no-far "far")
			 (format nil "~$,~$,~$" (3vector-x direction) (3vector-y direction) (3vector-z direction)))
		     (if (evenp plot)	;actual triangle set
			 (unless (and (zerop plot) no-far)
			   (multiple-value-bind (position corner) (floor point 6)
			     (and (< position (length triangles))
				  (< corner 4)				  ;double NIL to break between surfaces
				  (let ((data (aref triangles position))) ;(triangles . directions)
				    (and (if (zerop plot) (null (cdr data)) ;Triangles with no cusp
					   ;;includes this direction?
					   (member direction (cdr data) :key #'cusp-info-direction))
					 (values-list (3vector-list (nth (if (= corner 3) 0 corner)
									 (car data)))))))))
		       (and (< point 200)
			    (let ((angle (* point (/ pi 100))))	;0..2pi
			      (multiple-value-bind (x y) (two-perpendicular-vectors direction)
				(values-list (3vector-list (3vector+ (3vector-scale direction (cos threshold1))
								     (3vector-scale x (* (sin threshold1) (cos angle)))
								     (3vector-scale y (* (sin threshold1) (sin angle)))			   ))))))))))
	     :3d t
	     :styles :lines
	     gnuplot-keys)))


;;Plotting of gravitational radiation

;;Gravitational spectra plotting
(defun plot-total-gravitational-spectrum (a-hats a-sigmas b-hats b-sigmas direct-n &rest keys &key (bin-size 2.0) &allow-other-keys)
  (apply #'plot-binned-gravitational-spectrum
	 (apply #'total-gravitational-spectrum a-hats a-sigmas b-hats b-sigmas direct-n keys)
	 bin-size
	 keys))

;;Plot binned spectrum
(defun plot-binned-gravitational-spectrum (power bin-size &rest keys)
  (apply #'gnuplot 1 (length power)
	 #'(lambda (plot point)
	     (declare (ignore plot))
	     (and (numberp point)
		  (values (expt bin-size point) (aref power point))))
	 :logscale '(:x :y)
	 :prelude (format nil "set format '%h'~%")
	 keys))

;;Plot projected angular distribution of something
;;There is one point at the north pole, one at the south pole, and STEPS values of phi for each
;;of the STEPS intermediate values of theta
;;It would be better the colors corresponded to the the power at the center of the region that is colored,
;;instead of the average of the corners.
;;We could also do it by inverse mapping
(defun plot-projected-angular-distribution (function steps
						     &rest gnuplot-keys
						     &key cbmax ;Value at which to saturate the color scale.
						     (colorbox t) ;if NIL, do not show it
						     &allow-other-keys)
  (let ((npoints (+ 3 (* steps (1+ steps)))) ;2 polar points, s^2 regular points, 1+s breaks
	(dtheta (/ pi (1+ steps)))
	(dphi (/ (* 2 pi) steps)))
    (apply #'gnuplot 1 npoints
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (when (numberp point)
		 (let ((x nil) y z)
		   (cond ((zerop point) (setq x 0.0 y 0.0 z 1.0))
			 ((= point 1))	;First break
			 ((= point (1- npoints)) (setq x 0.0 y 0.0 z -1.0))
			 (t (multiple-value-bind (theta-step phi-step) (floor (- point 2) (1+ steps))
			      (unless (= phi-step steps) ;break
				(let* ((theta (* (1+ theta-step) dtheta))
				       (phi (* (- phi-step (/ (1- steps) 2.0)) dphi)))
				  (setq x (* (sin theta) (cos phi))
					y (* (sin theta) (sin phi))
					z (cos theta)))))))
		   (and x
			(multiple-value-call
			 #'values (Mollweide x y z)
			 (values (funcall function (make-3vector x y z)))
			 )))))
	   :prelude (format nil "~@[set cbrange [0:~S]~]
set view 0,0,1.8,1.8~%set pm3d ftriangles~%unset border~%unset tics~%set cbtics
unset lmargin~%unset rmargin~%unset tmargin~%unset bmargin
~:[unset colorbox~;set colorbox horizontal user origin 0.01,0.95 size 0.2, 0.05~]
"
			    cbmax colorbox)
	   :styles :pm3d	       
	   :3d t :key nil
	   gnuplot-keys)))

;;Overlay tangent vectors on distribution
;;Many features are not included
(defun plot-projected-angular-distribution-and-tangent-vectors (function steps data
						     &rest gnuplot-keys
						     &key styles
						     cbmax ;Value at which to saturate the color scale.
						     (colorbox t) ;if NIL, do not show it
						     &allow-other-keys)
  (setq data (coerce data 'vector))
  (unless styles
    (setq gnuplot-keys (list* :styles
			      (append (list :pm3d)
				      (loop for (vectors title style) across data collect style))
			      gnuplot-keys)))
  (let ((npoints (+ 3 (* steps (1+ steps)))) ;2 polar points, s^2 regular points, 1+s breaks
	(dtheta (/ pi (1+ steps)))
	(dphi (/ (* 2 pi) steps)))
    (apply #'gnuplot (+ (length data) 1)
	   (max npoints (loop for (vectors) across data maximize (1+ (length vectors))))
	   #'(lambda (plot point)
	       (when (numberp point)
		 (if (zerop plot)	;distribution
		     (let ((x nil) y z)
		       (cond ((> point npoints)) ;Not part of distribution
			     ((zerop point) (setq x 0.0 y 0.0 z 1.0))
			     ((= point 1)) ;First break
			     ((= point (1- npoints)) (setq x 0.0 y 0.0 z -1.0))
			     (t (multiple-value-bind (theta-step phi-step) (floor (- point 2) (1+ steps))
				  (unless (= phi-step steps) ;break
				    (let* ((theta (* (1+ theta-step) dtheta))
					   (phi (* (- phi-step (/ (1- steps) 2.0)) dphi)))
				      (setq x (* (sin theta) (cos phi))
					    y (* (sin theta) (sin phi))
					    z (cos theta)))))))
		       (and x
			    (multiple-value-call
			     #'values (Mollweide x y z)
			     (values (funcall function (make-3vector x y z)))
			     )))
		   ;;tangent vectors
		   (let ((vector (car (aref data (1- plot)))))
		     (when vector
		       (when (= point (length vector))
			 (setq point 0))
		       (let ((p (and (< point (length vector)) (aref vector point))))
			 (and p
			      (multiple-value-bind (sx sy) (apply #'Mollweide (3vector-list p))
				(values-list (list sx sy 0.0) ;Z coordinate
					     )))))))))
	   :prelude (format nil "~@[set cbrange [0:~S]~]
set view 0,0,1.8,1.8~%set pm3d ftriangles~%unset border~%unset tics~%set cbtics
unset lmargin~%unset rmargin~%unset tmargin~%unset bmargin
~:[unset colorbox~;set colorbox horizontal user origin 0.01,0.95 size 0.2, 0.05~]
"
			    cbmax colorbox)
	   :3d t :key nil
	   gnuplot-keys)))

;;Plot angular distribution of gravitational wave emission
(defun plot-gravitational-power (a-hats a-sigmas b-hats b-sigmas direct-n steps &rest keys &key (total-n direct-n) (bin-size 2.0)
					     &allow-other-keys)
  (apply #'plot-projected-angular-distribution
	 #'(lambda (direction)
	     (gravitational-power a-hats a-sigmas b-hats b-sigmas direction direct-n
				       :total-n total-n :bin-size bin-size))
	 steps
	 keys))

(defun plot-gravitational-power-one-bin-and-tangent-vectors (a-hats a-dsigmas b-hats b-dsigmas direct-n steps bin
						&rest keys &key (total-n direct-n) (bin-size 2.0)
						(style :points)	;Style for a' and b'
						interpolate
						&allow-other-keys)
  (mirror-image-let ((smooth-a-hats (if interpolate (interpolate-hats a-hats :split-count interpolate) a-hats)))
    (apply #'plot-projected-angular-distribution-and-tangent-vectors
	   #'(lambda (direction)
	       (aref (binned-gravitational-spectrum a-hats (dsigma-to-sigma a-dsigmas)
						    b-hats (dsigma-to-sigma b-dsigmas)
						    direction direct-n
						    :total-n total-n :bin-size bin-size)
		     bin))
	   steps
	   (append (list (list a-hats "A" style)
			 (list b-hats "B" style))
		   (and interpolate
			(list (list smooth-a-hats nil "lines lt 1")
			      (list smooth-b-hats nil "lines lt 2"))))
	   keys)))


(defun plot-ab-and-power (loop n-harmonics n-points &key reuse power-keys ab-keys)
  (setq loop (rest-frame-loop loop))
  (multiplot (:rows 2 :reuse reuse :prelude (format nil "set terminal x11 size 700,600~%"))
    (multiplot-format "set origin 0.0,0.45~%")
    (apply #'plot-gravitational-power loop n-harmonics n-points power-keys)
    (multiplot-format "set size 1.0,0.45~%")
    (apply #'plot-ab loop :rotate nil ab-keys)))

;;Different version
(defun plot-ab-and-power-1 (loop n-harmonics n-points
				 &rest gnuplot-keys &key (a-style "p lt 2") (b-style "p lt 4")
				 &allow-other-keys)
  (setq loop (rest-frame-loop loop))
  (mirror-image-let ((a-data (get-a-hats loop)))
    (destructuring-bind (power-plots power-points power-func &key prelude &allow-other-keys)
	(call-with-hats #'plot-gravitational-power loop n-harmonics n-points :return-args t)
      (assert (= power-plots 1))
      (apply  #'gnuplot 3		;power, a, b
	      (max (length a-data) (length b-data) power-points)
	      #'(lambda (plot point)
		  (cond ((eq point :title)
			 (nth plot '(nil "A" "B")))
			((zerop plot)	;Power
			 (and (< point power-points)
			      (funcall power-func plot point)))
			(t (let ((data (if (= plot 1) a-data b-data)))
			     (and (< point (length data))
				  (let ((p (aref data point)))
				    (multiple-value-bind (sx sy) (apply #'Mollweide (3vector-list p))
				      (values sx sy 0.0))))))))
	      :prelude prelude
	      :styles (list :pm3d a-style b-style)
	      :3d t
	      gnuplot-keys))))

;;Test function for sphere plot
(defun test-sphere (theta-steps &rest keys &key (phi-steps (* 2 theta-steps)) cbmax &allow-other-keys)
  (apply #'gnuplot 1 (* (1+ theta-steps) (+ 2 phi-steps))
	 #'(lambda (plot point)
	     (declare (ignore plot))
	     (when (numberp point)
	       ;;theta-step = 0...step
	       ;;phi-step = 0...phi-steps + 1, to leave room for blank record
	       (multiple-value-bind (theta-step phi-step) (floor point (+ 2 phi-steps))
		 (if (> phi-step phi-steps) ;end of circle of fixed theta
		     nil
		   (let* ((theta (* theta-step (/ pi theta-steps)))
			  (phi (* phi-step (/ (* 2 pi) phi-steps)))
			  (direction (spherical-coordinates theta phi)))
		     (values-list
		      (append
		       (3vector-list direction)
		       ;;This previous code doesn't work.  Maybe there's something wrong with the choice of corners?
		       ;; (list (+ (if (zerop phi-step) 0.0 1.0) (if (zerop theta-step) 0.0 1.0)))
		       (list (+ (/ phi-step (float phi-steps)) (/ theta-step (float theta-steps))))
		       )))))))
	     :3d t
	 :styles :pm3d
	 :prelude (format nil "~@[set cbrange [0:~S]~]~%set pm3d depthorder corners2color c2~%unset border
unset tics~%set cbtics~%unset lmargin~%unset rmargin~%unset tmargin~%unset bmargin~%" cbmax)
	 keys))


;;Batch processing
(defun compute-radiated-energy (loop-file direct-n &key total-n split-levels (output-file (merge-pathnames "energy-momentum.lisp" loop-file)))
  (initialize)
  (read-one-loop loop-file)
  (with-group-write-access
   (ensure-directories-exist output-file)
   (with-open-file (stream output-file :direction :output :if-exists :supersede)
     (multiple-value-bind (energy momentum)
	 (call-with-hats (rest-frame-loop (longest-loop))
			 #'gravitational-energy-momentum
			 direct-n :total-n total-n :split-levels split-levels)
       (format stream ";;energy, momentum~%(~S ~S)~%" energy momentum)))))

;;Loop files are directory/N/loop-N.dat
(defun run-radiated-energy (directory direct-n &key total-n split-levels)
  (load "load")
  (loop for number in (with-open-file (stream (format nil "~A/list-of-loops.dat" directory)) (read stream))
	for file = (format nil "~A/~D/loop-~D.dat" directory number number)
	while (probe-file file)
	for input-directory = (pathname-directory file)
	for output-file = (merge-pathnames (format nil "energy-momentum-~D.lisp" number)
					   (make-pathname :directory (list :absolute "cluster" "tufts" "strings" "energy-momentum"
									   (nth (- (length input-directory) 2) input-directory))))
	do (ensure-directories-exist output-file)
	unless (probe-file output-file)
	do (print number)
	  (do-submit (namestring (make-pathname :directory (pathname-directory output-file)))
		      (format nil "EM ~D" number)
		      (namestring (make-pathname :name (format nil "compute-EM-~D-" number) :type nil :defaults output-file))
		      batch-flags
		      `(compute-radiated-energy ,file ,direct-n :output-file ,(namestring output-file)
						:total-n ,total-n :split-levels ,split-levels)
		      :load-file (namestring (merge-pathnames "strings/parallel/load.lisp" (user-homedir-pathname))))))

(defun plot-radiated-energy (directory)
  (loop with slopes = '(0.1 0.2 0.3)
	for file in (directory (format nil "~A/energy-momentum-*.lisp" directory))
	for this = (with-open-file (stream file) (read stream nil))
	when this
	collect this into data
	finally (let ((max (loop for (energy) in data maximize energy)))
		  (format t "~D loops ~%" (length data))
		  (gnuplot (1+ (length slopes)) (length data)
			 #'(lambda (plot point)
			     (cond ((eq point :title)
				    (if (zerop plot) "loops"
				      (format nil "rocket fraction ~S" (nth (1- plot) slopes))))
				   ((zerop plot)
				    (values-list (nth point data)))
				   (t (case point
					(0 (values 0.0 0.0))
					(1 (values max (* max (nth (1- plot) slopes))))))))
			 :styles (cons :points (loop for entry in slopes collect :lines))
			 ))))

;;Different version
(defun plot-radiated-energy-1 (directory)
  (loop for file in (directory (format nil "~A/energy-momentum-*.lisp" directory))
	for this = (with-open-file (stream file) (read stream nil))
	when this
	collect this into data
	finally (format t "~D loops ~%" (length data))
		(gnuplot 1 (length data)
			 #'(lambda (plot point)
			     (declare (ignore plot))
			     (if (eq point :title) "loops"
			       (destructuring-bind (energy momentum) (nth point data)
				 (values energy (/ momentum energy))))))))

(defun histogram-rocket-effect (directory &rest gnuplot-keys &key (bins 100) &allow-other-keys)
  (let ((data (make-array bins :initial-element 0.0))
	(list nil))
    (loop for file in (directory (format nil "~A/energy-momentum-*.lisp" directory))
	  for (energy momentum) = (with-open-file (stream file) (read stream nil))
	  when energy
	  do  (let ((rocket (/ momentum energy)))
		(incf (aref data (floor (* rocket bins))))
		(push rocket list)))
    (setq list (sort list #'<))
    (format t "~D loops, average rocket fraction ~S, median ~S"
	    (length list)
	    (/ (reduce #'+ list) (length list))
	    (nth (floor (length list) 2) list))
    (apply #'gnuplot 1 (length data)
	     #'(lambda (plot point)
		 (declare (ignore plot))
		 (unless (eq point :title)
		   (values (/ (+ point 0.5) bins) (aref data point))))
	     :styles :boxes
	     :prelude (format nil "set style fill solid 0.7~%")
	     gnuplot-keys)))

(defun write-initial-rocket-fractions (directory)
  (with-open-file (stream (format nil "~A/initial-rocket-fractions.dat" directory)
			  :direction :output :if-exists :supersede)
   (loop for file in (directory (format nil "~A/energy-momentum-*.lisp" directory))
	  for (energy momentum) = (with-open-file (stream file) (read stream nil))
	  when energy
	  do (format stream "~S~%" (/ momentum energy)))))
	  
;;Smooth strings stored as their Fourier transforms

;;In this section, a scalar function on 0..1 is stored as its Fourier transform in the format (but not the
;;normalization)  returned by REALFT.  This consists of an array of 2*Nf real numbers representing frequencies 0...Nf-1
;;and frequencies -Nf+1...-1, which are just conjugates of the positive frequencies.
;;The first element is the zero mode, which vanishes if the function averages to 0.
;;The second element, which would be freq Nf (called N/2 by Numerical Recipes) must be 0.
;;Frequency j=1...Nf-1 has its real part in slot 2j and its imaginary part in slot 2j+1.
;;The function f(x) = sum_{j=-Nf}^Nf F_j exp{-2 pi i j x} = sum_{j=1}^Nf 2 Re (F_j exp{-2 pi i j x})
;;The inverse discrete Fourier transform defined by Numerical Recipes yields 1/(2 Nf) f(x_j) where x = j/(2 Nf).
;;Calling inverse REALFT on this data returns (1/2) f(x).

;;When we have a vector function, we store it as an array of 3 scalar amplitude arrays in the above format.

;;Compute f(x) from amplitudes as above.
(defun slow-transform-amplitudes (amplitudes x)
  (+ (aref amplitudes 0)
     (* 2 (loop for j from 1 below (floor (length amplitudes) 2)
		sum (+ (* (cos (* 2 pi j x)) (aref amplitudes (* 2 j))) 
		       (* (sin (* 2 pi j x)) (aref amplitudes (1+ (* 2 j)))))
		))))

;;Much faster: compute all sines and cosines by successive rotations.
;;I learned this idea from a web page of Inigo Quilez
(declaim (inline transform-amplitudes))
(defun transform-amplitudes (amplitudes x)
  (declare (type (simple-array double-float) amplitudes)
	   (double-float x)
	   (optimize speed)
	   (muffle-conditions compiler-note)) ;Avoid warning about consing return value
  (let ((cos 1.0)			;Go through sines and cosines by rotation
	(sin 0.0)			;Start with sin and cos of 0
	(sin1 (sin (* 2 pi x)))		;Rotate by this angle between successive steps
	(cos1 (cos (* 2 pi x)))
	(result (/ (aref amplitudes 0) 2))) ;zero mode has half amplitude relative to others
    (declare (double-float cos sin sin1 cos1 result))
    (loop for j from 1 below (floor (length amplitudes) 2)
	  do (psetq cos (- (* cos cos1) (* sin sin1)) ;Do rotation
		    sin (+ (* cos sin1) (* sin cos1)))
	  do (incf result (+ (* cos (aref amplitudes (* 2 j))) (* sin (aref amplitudes (1+ (* 2 j)))))))
    (* result 2)))
	  
;;Compute vector-valued function from amplitudes
(defun transform-vector-amplitudes (amplitudes x)
  (let ((result (make-3vector)))
    (dotimes (component 3)
      (setf (3vector-component result component) (transform-amplitudes (aref amplitudes component) x)))
    result))
  
;;Vector amplitudes are an array of amplitude vectors, but sometimes we want to get out a vector
;;for a given index
(defun vector-amplitudes-ref (amplitudes index)
  (make-3vector (aref (aref amplitudes 0) index) (aref (aref amplitudes 1) index)
		(aref (aref amplitudes 2) index)))

(defun vector-amplitudes-length (amplitudes)
  (length (aref amplitudes 0)))

;;Compute f(x_j), x=j/N, j=0...N-1 by FFT.
;;N will be rounded up to the number of real elements in amplitudes, and then to a power of 2.
(defun discrete-transform-amplitudes (amplitudes &optional (n (length amplitudes)))
  (setq n (expt 2 (ceiling (log (max n (length amplitudes)) 2))))
  (let* ((nf (floor n 2))
	 (result (make-array n :element-type 'double-float :initial-element 0.0)))
    (dotimes (i (length amplitudes))	;Install data for FFT.  If N is larger, extra frequencies have 0 amplitude
      ;;REALFT would return a result too small by factor 2, so we multiply by that in advance
      (setf (aref result i) (* (aref amplitudes i) 2)))
    (realft result nf :isign -1)
    result))

;;Accept evenly spaced hats, considering them to be samples from a smooth function a'
;;(i.e., sinc problems are not dealt with).  Fouier transform and
;;package amplitudes into a vector amplitude structure, as above
;;The loop should be in the rest frame.
(defun amplitudes-from-hats (hats)
  (let* ((n (length hats))		;Count of points, so we'll have freqs -n/2...0...n/2 with first and last equal
	 (amplitudes (make-array 3)))
    (unless (= (expt 2 (1- (integer-length n))) n)
      (error "The number of hats must be a power of two"))
    (dotimes (coordinate 3)
      (let ((data (make-array n :element-type 'double-float))) ;Array to return
	(dotimes (index n)
	  ;;The Fourier transform will involve n data points, so the inverse Fourier
	  ;;transform using Numerical Recipes' formula should have 1/n in front.  What we define above
	  ;;doesn't have this factor, so we should apply it in advance now.
	  (setf (aref data index) (/ (3vector-component (aref hats index) coordinate) n)))
	(realft data (/ n 2))		;Transform n real numbers
	(when (> (aref data 1) 1e-8)	;Amplitudes don't have frequency n/2, so it should be small
	  (warn "Ignoring freq ~D with amplitude ~S in ~S "
		(/ n 2) (aref data 1) 'amplitudes-from-hats))
	(setf (aref data 1) 0.0)
	(setf (aref amplitudes coordinate) data)))
    amplitudes))

;;Return an array of vector positions
(defun discrete-transform-vector-amplitudes (amplitudes &key n	   ;Number of points to return
							normalize) ;If set, normalize to unit vectors
  (setq n (max (or n 0) (vector-amplitudes-length amplitudes)))	   ;Must be at least number of real amplitudes
  (setq n (expt 2 (ceiling (log n 2))))	       ;Round up to power of 2
  (let ((result (make-array n)))
    (dotimes (i n)			;Create 3vectors to return
      (setf (aref result i) (make-3vector)))
    (dotimes (component 3)
      (let ((values (discrete-transform-amplitudes (aref amplitudes component) n))) ;Point values for this component
	(dotimes (i n)								    ;Copy data to 3vectors
	  (setf (3vector-component (aref result i) component) (aref values i)))))
    (when normalize
      (dotimes (i n)
	(setf (aref result i) (3vector-normalize (aref result i)))))
    result))

;;Convert amplitudes for f into amplitudes for f', by multiplying F_j by -2 pi i j
(defun differentiate-amplitudes (amplitudes)
  (let* ((n (length amplitudes))	;Count of real data in amplitudes.  # of frequencies = n/2
	 (new (make-array n :element-type 'double-float))) ;Array for results
    (setf (aref new 0) 0.0)
    (loop for j below (floor n 2)
	  for real-index = (* 2 j)
	  for imaginary-index = (1+ real-index)
	  do (setf (aref new real-index) (* 2 pi j (aref amplitudes imaginary-index)) ;Multiply by -2 pi i j
		   (aref new imaginary-index) (* -2 pi j (aref amplitudes real-index))))
    new))

;;Apply to vector of amplitudes
(defun differentiate-vector-amplitudes (amplitudes)
  (let ((new (make-array 3)))
    (dotimes (component 3)
      (setf (aref new component) (differentiate-amplitudes (aref amplitudes component))))
    new))

;;Rotate the object ascribed by the amplitudes with the given rotation matrix
(defun rotate-vector-amplitudes (amplitudes rotation)
  (let* ((frequencies (vector-amplitudes-length amplitudes))
	 (new-amplitudes (coerce (loop repeat 3 collect (make-array frequencies :element-type 'double-float))
				 'vector)))
    (loop for mode below frequencies
	  for vector = (make-array 3 :element-type 'double-float
				   :initial-contents (loop for array across amplitudes collect (aref array mode)))
	  for rotated = (dotmv rotation vector)
	  do (loop for array across new-amplitudes
		   for component across rotated
		   do (setf (aref array mode) component)))
    new-amplitudes))

;;Find a cusp, where a'=b' on a smooth string, starting from given positions xa and xb in 0..1.
;;A smooth string is given by Fourier components, and these unfortunately do not give a function which obeys |a'|=1
;;exactly.  So the best we can hope for is to find a'-hat = b'-hat.  So we normalize the output from the smooth
;;function.  Then we need the derivative of that, 
;;d(a-hat)/d(xa) = a''/|a'| - (a''.a')a'/|a'|^3. = (a'' - (a''.a-hat)a-hat)/|a'|
;;We assume we're in the rest frame, so there is no zero mode.
;;Returns cusp-info, q.v.
(defun find-smooth-cusp (a-amplitudes b-amplitudes x-a x-b &key (tolerance 1e-10) (max-steps 10))
  (mirror-image-let ((a-2-amplitudes (differentiate-vector-amplitudes a-amplitudes))) ;amplitudes for a'', b''
    (let ((m (make-array '(2 2) :element-type 'double-float))			    ;We will solve m.dx = d
	  (d (make-array 2 :element-type 'double-float))
	  (indx (make-array 2 :element-type '(unsigned-byte 16)))    ;Permutation array for LUDCMP
	  (vv (make-array 3 :element-type 'double-float))	     ;Working array for LUDCMP
	  (distance nil) difference cross
	  (steps 0))
      (loop
       (mirror-image-let*		;Get normalized a' and b' and their derivatives.  See header comment.
	   ((a-unnormalized (transform-vector-amplitudes a-amplitudes x-a))
	    (a-length (3vector-length a-unnormalized))
	    (a (3vector-scale a-unnormalized (/ 1 a-length)))	       ;a-hat
	    (a-2-raw (transform-vector-amplitudes a-2-amplitudes x-a)) ;Actual a''
	    (a-2 (3vector-scale (3vector- a-2-raw       ;d(a-hat)/d(xa).  See above.
					  (3vector-scale a (3vector-dot a-2-raw a)))
				(/ 1 a-length))))
	 (setq difference (3vector- b a))
	 (let ((new-distance (3vector-length difference)))
	   (when (< new-distance tolerance) ;Done?
	     (return (make-cusp-info
		      :x-a (mod x-a 1.0)  ;In case intersection is just at end of string
		      :x-b (mod x-b 1.0)
		      :direction (3vector-normalize (3vector-scale (3vector+ a b) 0.5)) ;Average "equal" values
		      :a-pp a-2
		      :b-pp b-2)))
	   (when (and distance		;Not first time
		      (> new-distance (* 0.9 distance)))
;;	     (format t "New distance = ~S~%" new-distance)
	     (error "Newton's method not making adequate progress")
	     )
	   (setq distance new-distance))
;;	 (format t "Xa = ~S, Xb = ~S~%a' = ~S, b' = ~S~%Distance ~S~%a'' = ~S, b'' = ~S~%" x-a x-b a b distance a-2 b-2)
	 (when (> (incf steps) max-steps) ;Should we try another step?
	   (error "Newton's method did not converge after ~D steps" max-steps))
	 (setq cross (3vector-normalize (3vector-cross-product a b))) ;unit vector perpendicular to a', b'
	 (setf (aref m 0 0) (- (/ (3vector-dot a-2 b)		      ;Effect of dxa in direction of difference
				  distance))
	       (aref m 0 1) (- (/ (3vector-dot b-2 a) ;Effect of dxb in direction of difference
				  distance))
	       (aref m 1 0) (- (3vector-dot a-2 cross)) ;Effect of dxa in crosswise direction difference
	       (aref m 1 1) (3vector-dot b-2 cross)     ;Effect of dxa in crosswise direction difference
	       (aref d 0) (- distance)			;Desired change in direction of difference
	       (aref d 1) 0.0)				;Desired change in crosswise direction
;;	 (let ((m0 (copy-matrix m)))
	   (unless (ludcmp m indx vv)	;Solve m.dx = d
	     (error "LUDCMP failed"))
	   (lubksb m indx d)		;Now d becomes desired dx
#|	   (format t "Predicted change to (dist, cross) is ~S~%" (dotmv m0 d))
	   (let ((diff1 (3vector-normalize difference))
		 (new-a (3vector+ a (3vector-scale a-2 (aref d 0))))
		 (new-b (3vector+ b (3vector-scale b-2 (aref d 1)))))
	     (mirror-images (format t "Perpendicular in ~A is ~S~%" :a (3vector-dot a a-2)))
	     (format t "Old distance in direction of difference ~S~%" (3vector-dot diff1 (3vector- b a)))
	     (format t "Predicted A= ~S~%" new-a)
	     (format t "Predicted B= ~S~%" new-b)
	     (format t "Predicted distance in direction of old difference ~S~%"
		     (3vector-dot diff1 (3vector- new-b new-a))))
	   |#
	 (incf x-a (aref d 0))		;Make changes
	 (incf x-b (aref d 1)))))))

;;Find cusps on smoothed loop.  We find N a' and N b' positions and look for crossings by brute force
;;Returns a list of (c a'' b'') -- see find-smooth-cusp.
(defun find-smooth-cusps (a-amplitudes b-amplitudes
				       &key (n (* 10 ;Must be large enough to separate crossings
						  (max (vector-amplitudes-length a-amplitudes)
						       (vector-amplitudes-length b-amplitudes))))
				       (tolerance 1e-10))
  (format t "Finding cusps...") (force-output)
  (mirror-image-let ((a (discrete-transform-vector-amplitudes a-amplitudes :n n :normalize t))
		     (result nil))
    (setq n (length a))			;What we actually got.  Might be rounded up.
    (dotimes (i n)
      (loop with a1 = (aref a i)
	    with a2 = (aref a (mod (1+ i) n))
	    for j below n
	    when (spherical-intersection a1 a2 (aref b j) (aref b (mod (1+ j) n)))
	    do (pushnew (find-smooth-cusp a-amplitudes b-amplitudes
					  (/ (+ i 0.5) n) (/ (+ j 0.5) n) ;initial guess: midpoint
					  :tolerance tolerance)
			result
			:test #'(lambda (data1 data2) ;Avoid duplicates.
				  (and (fudge=mod (cusp-info-x-a data1) (cusp-info-x-a data2)
						  1.0 (* 2 tolerance))
				       (fudge=mod (cusp-info-x-b data1) (cusp-info-x-b data2) 
						  1.0 (* 2 tolerance)))))))
    result))

;;Naively return the maximum a''.  N is the number of points to use.  The default is twice the number of real
;;amplitudes, i.e., 4 times the maximum frequency.  This means that we have at least 4 points to resolve the
;;most rapid fundamental oscillation.
(defun maximum-second-derivative (amplitudes &optional n)
  (setq n (if n (max n (vector-amplitudes-length amplitudes)) ;Always as many points as real amplitudes.
	    (* 2 (vector-amplitudes-length amplitudes))))		      ;But default is twice that
  (setq n (expt 2 (ceiling (log n 2))))	       ;Round up to power of 2
  (let ((app-amplitudes (differentiate-vector-amplitudes amplitudes))		;Amplitudes for a''
	(app (make-array n :element-type 'double-float :initial-element 0.0)))	;Construct a'' with n points
    (loop for amplitudes-1 across app-amplitudes
	  for points = (discrete-transform-amplitudes amplitudes-1 n) ;Get a''_i at N evenly spaced points
	  do (dotimes (j n)
	       (incf (aref app j) (expt (aref points j) 2)))) ;sum into |a''|^2
    (sqrt (reduce #'max app))))		;Maximum |a''|^2 at evenly spaced points.  Should adjust by sum(sinc)

;;Naively return the cosine of the closest approach of a' to Omega
;;N is the number of points to use.  See comment on maximum-second-derivative.
(defun closest-approach (amplitudes omega &optional n)
  (setq n (if n (max n (vector-amplitudes-length amplitudes)) ;Always as many points as real amplitudes.
	    (* 2 (vector-amplitudes-length amplitudes))))		      ;But default is twice that
  (setq n (expt 2 (ceiling (log n 2))))	       ;Round up to power of 2
  (let ((data (make-array n :element-type 'double-float)))			;Data for Fourier transform
    (dotimes (j (vector-amplitudes-length amplitudes))				;Make amplitudes for a' . Omega
      (setf (aref data j) (3vector-dot (vector-amplitudes-ref amplitudes j) omega)))
    (let ((az (discrete-transform-amplitudes data n))) ;Get a'_i . Omega at N evenly spaced points
      (reduce #'max az))))			     ;Biggest cosine

;;Returns scalar amplitudes for x', y', z' and z, where z = omega is the observation direction,
;;and x and y are two arbitrary perpendicular directions
;;The position of the loop is arbitrary.
(defun gravitational-wave-basis (amplitudes omega)
  (let* ((nfreq (floor (vector-amplitudes-length amplitudes) 2)) ;Frequencies 1... Nf-1 are used
	 (basis (append (multiple-value-list (two-perpendicular-vectors omega)) (list omega))) ;x, y, z unit vectors
	 ;;Amplitudes for a'_x, a'_y, and a'_z
	 (basis-amplitudes (loop repeat 3 collect (make-array (* 2 nfreq) :element-type 'double-float)))
	 (xprime-amplitudes (first basis-amplitudes)) ;Split into separate variables
	 (yprime-amplitudes (second basis-amplitudes))
	 (zprime-amplitudes (third basis-amplitudes))
	 (z-amplitudes (make-array (* 2 nfreq) :element-type 'double-float))) ;This is for a_z rather than a'_z
    (loop for vector in basis						      ;Get components in our basis separately
	  for array in basis-amplitudes
	  do (setf (aref array 0) 0.0	;Rest frame.  No zero mode.
		   (aref array 1) 0.0)
	  do (loop for index from 2 below (* 2 nfreq)
		   do (setf (aref array index)
			    (3vector-dot vector (vector-amplitudes-ref amplitudes index)))))
    (loop for j from 1 below nfreq		 ;Find amplitudes for a_z from those for a'_z by dividing by -2 pi i j,
	  for jj = (* 2 j)			 ;i.e. multiplying by i and dividing by 2 pi j
	  do (setf (aref z-amplitudes (1+ jj))	 ;Imaginary part
		   (/ (aref zprime-amplitudes jj) 2 pi j))	    ;is rescaled real
	  do (setf (aref z-amplitudes jj)			    ;Real part
		    (/ (aref zprime-amplitudes (1+ jj)) -2 pi j)))   ;rescale imaginary and change sign
    (values xprime-amplitudes yprime-amplitudes zprime-amplitudes z-amplitudes)))

;;The number of points that we can easily generate is a power of 2.  But NUFFT is supposed to oversample
;;by some factor alpha when generating the points to actually be transformed by FFT.  This also has
;;to be a power of 2, so the minimum alpha we could use is 2.  But it appears that if we don't oversample,
;;everything is fine except for the last small fraction of frequencies, so we just discard those.
(defparameter smooth-gravitational-wave-discard-fraction 0.1)

;;Accept a'(v) [not a] as array of Fourier amplitudes with the structure above.
;;Return two arrays indexed by harmonic number k, with the k = 0 slot zero:
;; |I_perp^k|^2 and Im(I_x I_y^*), where x, y, and Omega form a righthand coordinate system and
;; I^(k)(Omega) = (1/L) int_0^L dv a'(v) exp(2 pi i k (v-a_z(v))/L)
;; where z means the component in the observation direction omega.
;;FREQUENCIES is the number of frequencies that you are asking for.  We actually need to compute more than this
;;so it is not desirable for this to be a power of 2.
;;As with all Fourier transforms, you should have negligible power in frequencies that are not returned, or else
;;you'll get aliasing.
;;See Fourier.tex
(defun smooth-gravitational-wave (amplitudes omega frequencies)
  ;;Number of good frequencies will be (1-smooth-gravitational-wave-discard-fraction)N/2, where N is the
  ;;number of points to pass to NUFFT.  So N = 2*frequencies/(1-smooth-gravitational-wave-discard-fraction), rounded
  ;;up to a power of 2, because that's what we can generate in discrete-transform-amplitudes.
  (let ((n (expt 2 (ceiling (log (/ (* frequencies 2) (- 1 smooth-gravitational-wave-discard-fraction)) 2)))))
    ;;Usually we round up in the calculation above, and therefore we get more good frequencies than the user
    ;;asked for, so we might as well return them
    (setq frequencies (floor (* (floor n 2) (- 1 smooth-gravitational-wave-discard-fraction)))) ;Number to return
    (multiple-value-bind (xprime-amplitudes yprime-amplitudes zprime-amplitudes z-amplitudes)
	(gravitational-wave-basis amplitudes omega)
      (declare (ignore zprime-amplitudes))
      (let ((az (discrete-transform-amplitudes z-amplitudes n))	      ;a_z at n equally-spaced sigma
	    (gx (discrete-transform-amplitudes xprime-amplitudes n))  ;a'_x
	    (gy (discrete-transform-amplitudes yprime-amplitudes n))) ;a'_y
	(functional-gravitational-wave-1 az gx gy frequencies)))))

;;See smooth-gravitational-wave.  We accept a_z, a'_x, and a'_y samples at evenly spaced sigma
;;a constant offset in a_z is removed
(defun functional-gravitational-wave-1 (az gx gy frequencies)
  (let* ((n (length az))
	 (f (make-array n :element-type 'double-float))		      ;f = sigma - a_z(sigma)
	 (perp2 (make-array frequencies :element-type 'double-float)) ;|I_perp|^2
	 (ixy (make-array frequencies :element-type 'double-float)))  ;Im (I_x I_y*)
    (declare (optimize speed)
	     (type (unsigned-byte 31) n frequencies) ;So n^2 is a fixnum
	     (type (simple-array double-float (*)) f az gx gy))
    (dotimes (i n)
      (setf (aref f i) (- (/ (double-float i) n) (aref az i) (aref az 0)))) ;Fill in f.  Offset so f(0) = 0.
    ;;We want to compute I^n = int_0^1 g(sigma) exp(2 pi i n f(sigma), which we approximate
    ;;by (1/N) Sum g(sigma_j) exp(2 pi i n f(sigma_j) with sigma_j = j/N
    (let ((ix (nufft f gx *nufft-m*))	   ;This is complex conjugate of sum above
	  (iy (nufft f gy *nufft-m*)))	   ;The length of these arrays is n/2, but we don't want the last fraction
      (declare (type (simple-array double-float (*)) ix iy))
      ;;Combine I_x, I_y into results.
      (setf (aref perp2 0) 0.0		;Unused slots
	    (aref ixy 0) 0.0)
      (loop with nsquared = (the fixnum (expt n 2))
	    for j from 1 below frequencies
	    for jj = (* 2 j)
	    do (setf (aref perp2 j)
		     (/ (+ (expt (aref ix jj) 2) (expt (aref ix (1+ jj)) 2)	 ;|I_x|^2
			   (expt (aref iy jj) 2) (expt (aref iy (1+ jj)) 2))	 ;|I_y|^2
			nsquared)
		     (aref ixy j)
		     ;;Our NUFFT routine has opposite phase convention, so we must invert here
		     (/ (- (* (aref ix jj) (aref iy (1+ jj))) ;-Im (I_x I_y^*) = Re I_x Im I_y - Im I_x Re I_y
			   (* (aref ix (1+ jj)) (aref iy jj)))
			nsquared)))
      (values perp2 ixy))))

;;Return information like smooth-gravitational-wave.  FUNCTION returns values a(x=0...1), da/dx
(defun functional-gravitational-wave (function omega frequencies)
  (let ((n (expt 2 (ceiling (log (/ (* frequencies 2) (- 1 smooth-gravitational-wave-discard-fraction)) 2)))))
    (setq frequencies (floor (* (floor n 2) (- 1 smooth-gravitational-wave-discard-fraction))))
    (multiple-value-bind (x-direction y-direction) (two-perpendicular-vectors omega)
      (let ((az (make-array n :element-type 'double-float))
	    (gx (make-array n :element-type 'double-float))
	    (gy (make-array n :element-type 'double-float)))
	(loop for index below n
	      for x = (/ index (double-float n))
	      do (multiple-value-bind (a g) (funcall function x)
		   (setf (aref az index) (3vector-dot omega a)
			 (aref gx index) (3vector-dot x-direction g)
			 (aref gy index) (3vector-dot y-direction g))))
	(functional-gravitational-wave-1 az gx gy frequencies)))))

;;Return |I_perp^n|^2 and Im(I_x I_y^*), where x, y, and Omega form a righthand coordinate system and
;;I^(n)(Omega) = (1/L) int_0^L dv a'(v) exp(2 pi i k (v-a_z(v))/L)
;;where z means the component in the observation direction omega.
;;Computation by direct integration for a single harmonic
(defun slow-smooth-gravitational-wave (amplitudes omega n)
  (multiple-value-bind (xprime-amplitudes yprime-amplitudes zprime-amplitudes z-amplitudes)
      (gravitational-wave-basis amplitudes omega)
    (declare (ignore zprime-amplitudes))
    ;;We need to do 4 integrals, for x and y and real and imaginary
    (flet ((integral (prime-amplitudes imaginary-p)
	      (integrate-subdivide #'(lambda (x)	;x = 0...1
				       (* (transform-amplitudes prime-amplitudes x)
					  (funcall (if imaginary-p #'sin #'cos)
						   (* 2 pi n (- x (transform-amplitudes z-amplitudes x))))))
				   0.0 1.0 (ceiling n 50)
				   :limit (expt 10 5))))
      (let ((xr (integral xprime-amplitudes nil))
	    (xi (integral xprime-amplitudes t))
	    (yr (integral yprime-amplitudes nil))
	    (yi (integral yprime-amplitudes t)))
	(values (+ (expt xr 2) (expt xi 2) (expt yr 2) (expt yi 2)) ;I_perp^n|^2
		(- (* xi yr) (* xr yi)))))))			    ;Im(I_x I_y^*)

;;Sometimes the integration system erroneously decides that the integral is accurate when it has not used enough
;;points to capture all the oscillations.  Probably there is an accidental resonance between the oscillations and the
;;choice of the integration points.  By subdividing into intervals with not too many oscillations, we work around
;;this problem.
(defun integrate-subdivide (function a b subdivisions &rest keys)
  (loop with division-size = (/ (- b a) subdivisions)
	for n below subdivisions
	sum (apply #'integrate function (+ a (* division-size n)) (+ a (* division-size (1+ n)))
		   keys)))

;;Accept a'(v) as an array of Fourier amplitudes as in smooth-gravitational-wave.  Returns
;;an array of power in the different frequencies.  Slot 0 is set to 0.
;;FREQUENCIES as the minimum number to use, but for various reasons you may get more than that returned.
(defun smooth-gravitational-spectrum-2 (a-amplitudes b-amplitudes omega a-frequencies b-frequencies)
  ;;Use at least as many frequencies as the number that are used in the specification of the string
  ;;Otherwise you'll certainly get aliasing
  (mirror-images (setq a-frequencies (max a-frequencies (vector-amplitudes-length a-amplitudes))))
  (multiple-value-bind (iperp ixy) (smooth-gravitational-wave a-amplitudes omega a-frequencies) ;|I_perp|^2, Im(I_xI_y^*)
    (multiple-value-bind (jperp jxy) (smooth-gravitational-wave b-amplitudes omega b-frequencies) ;Same for J
      (spectrum-from-ij iperp ixy jperp jxy))))

(defun functional-gravitational-spectrum-1 (a-function b-function omega frequencies)
  (multiple-value-bind (iperp ixy) (functional-gravitational-wave a-function omega frequencies)
    (multiple-value-bind (jperp jxy) (functional-gravitational-wave b-function omega frequencies)
      (spectrum-from-ij iperp ixy jperp jxy))))

;;Power in one harmonic
(defun slow-smooth-gravitational-power (a-amplitudes b-amplitudes omega n)
  (multiple-value-bind (iperp ixy) (slow-smooth-gravitational-wave a-amplitudes omega n) ;|I_perp|^2, Im(I_xI_y^*)
    (multiple-value-bind (jperp jxy) (slow-smooth-gravitational-wave b-amplitudes omega n) ;Same for J
      (* 8 pi (expt n 2)
	 (+ (* iperp jperp)
	    (* 4 ixy jxy))))))

;;Maximum ratio between number of modes of a and modes of b to compute
(defparameter smooth-gravitational-max-harmonic-contrast 1) ;Disable this feature by default.  Always use name number.

;;In some cases (especially analytic things like the burden loop), there's no radiation in certain directions,
;;and this calculation would otherwise yield n=0, which we probably don't want.
(defparameter smooth-gravitational-min-max-harmonic 4)

;;If the cusp direction is this close to a corner, we don't split the triangle
(defparameter fudge-cusp-triangle-corner 1e-8)

;;If a given harmonic yields xi- and xi+ both less than this, do it
;; using the bulk code.
(defparameter cusp-bulk-max-xi 0.02)
(defparameter cusp-bulk-max-n 100)	;Limit on number of modes ever to do this way

;;Find maximum harmonic to do using the bulk code
(defun cusp-bulk-max-n (cusp-list triangle)
  (let ((result (min cusp-bulk-max-n
		     (loop for cusp-info in cusp-list
			   minimize (cusp-bulk-max-n-1 cusp-info triangle))))) ;Must be below all limits
    (when (< result 4)
      (warn "Only ~D harmonics computed with bulk code.  This might be inaccurate" result))
    result))

(defun cusp-bulk-direction-max-n (cusp-list direction)
  (let ((result (min cusp-bulk-max-n
		     (loop for cusp-info in cusp-list
			   minimize (cusp-bulk-max-n-2 cusp-info
						       (spherical-angle direction (cusp-info-direction cusp-info)))))))
    (when (< result 4)
      (warn "Only ~D harmonics computed with bulk code.  This might be inaccurate" result))
    result))

;;Max harmonics that we could do with bulk code for single cusp
(defun cusp-bulk-max-n-1 (cusp-info triangle)
  (destructuring-bind (a b c) triangle
    (let* ((cusp (cusp-info-direction cusp-info))
	   (maxtheta (max (spherical-angle cusp a) (spherical-angle cusp b) (spherical-angle cusp c))))
      (cusp-bulk-max-n-2 cusp-info maxtheta))))

(defun cusp-bulk-max-n-2 (cusp-info maxtheta)
  (mirror-image-let* ((a-pp (cusp-info-a-pp cusp-info))
		      (a-pp-length (3vector-length a-pp)))
    (floor (/ (* 3 cusp-bulk-max-xi (min a-pp-length b-pp-length))
	      pi (expt maxtheta 3)))))

;;Integrate the cusp radiation over the given triangle
;;We compute low harmonics using the bulk code, as determined by cusp-bulk-max-xi
;;and use the cusp code only for higher harmonics
(defun cusp-bulk-triangle-integral (a-amplitudes b-amplitudes cusp-list triangle
						 &key (nmax (cusp-bulk-max-n cusp-list triangle)))
  (+ (loop for cusp-info in cusp-list	;Add up contributions of all relevant cusps
	   sum (cusp-triangle-integral cusp-info triangle nmax))
     (if (plusp nmax)			;Can be zero if very large triangles, e.g. split-levels = 0
	 (* (apply #'spherical-triangle-area triangle) ;Compute low harmonics once for each triangle
	    (smooth-gravitational-power-first-n a-amplitudes b-amplitudes (triangle-center triangle) nmax))
       0)))

;;Cusp calculation without the first BULK-N modes
(defun cusp-triangle-integral (cusp-data triangle &optional bulk-n)
  (loop for (a b weight) in (cusp-triangle-decomposition cusp-data triangle)
	sum (* weight (cusp-triangle-integral-1 cusp-data a b bulk-n))))

;;Decompose a triangle into a set of triangles with weights (+/-1) with 1 corner at the cusp direction
(defun cusp-triangle-decomposition (cusp-data triangle)
  (let ((cusp (cusp-info-direction cusp-data)))		;Direction of cusp
    (destructuring-bind (a b c) triangle
      (cond ((4vector= cusp a fudge-cusp-triangle-corner)		;Near one point
	     (list (list b c 1)))
	    ((4vector= cusp b fudge-cusp-triangle-corner)
	     (list (list c a 1)))
	    ((4vector= cusp c fudge-cusp-triangle-corner)
	     (list (list a b 1)))
	    ((point-inside-spherical-triangle cusp triangle) ;Inside.  Do 3 subtriangles
	     (list (list a b 1)
		   (list b c 1)
		   (list c a 1)
		   ))
	    (t			;Outside.  Must find which of the vertices is in the middle.
	     (let ((ab (spherical-triangle-angle cusp a b))
		   (bc (spherical-triangle-angle cusp b c))
		   (ca (spherical-triangle-angle cusp c a)))
	       (if (> ab bc)
		   (when (> ab ca)	;AB biggest
		     (rotatef a b c))	;Rotate so C goes into B
		 (when (> bc ca)	;BC biggest?
		   (rotatef c b a)))	;Rotate so A goes into B
	       ;;Otherwise CA is biggest and B is already in the middle
	       ;;We integrate the power radially from the cusp out to segments ab and bc
	       ;;and also the power out to ac.  Then we subtract in one order or the other to get the desired power
	       ;;But to get the triangle points in the right order, we need separate code for the two cases
	       (if (> (spherical-triangle-angle a b cusp) ;angle at a between radius and line to b
		      (spherical-triangle-angle a c cusp)) ;angle at a between radius and line to c
		   ;;Angle to b is larger: b is beyond ac and a has smallest phi angle
		   (list (list a b 1) 
			 (list b c 1)
			 (list a c -1))
		 ;;No: angle to b is smaller: b is inside ac and c has smallest phi angle
		 (list (list c a 1)
		       (list c b -1)
		       (list b a -1)))
	       ))))))
    
;;Integrate radiation over a triangle with one point at the cusp direction.  The other two points are FIRST
;;and SECOND, in counterclockwise order (seen from outside the sphere).
;;If BULK-N only include modes above that number
(defun cusp-triangle-integral-1 (cusp-data first second &optional bulk-n)
  (let ((cusp (cusp-info-direction cusp-data)))
;;    (format t "~S, ~S, ~S, area ~S~%" cusp first second (spherical-triangle-area cusp first second))
    ;;If triangle has no area, just give zero.  This avoids a problem where QROMB is trying too hard to
    ;;get accurate results out of an integrand that fluctuates around 0 by numerical error
    (when (spherical-triangle-degenerate-p cusp first second)
      (return-from cusp-triangle-integral-1 0.0))
    (mirror-image-let* ((a-pp (cusp-info-a-pp cusp-data))
			(a-pp-length (3vector-length a-pp)) ;Separate a'' and b'' into length
			(a-pp-1 (3vector-normalize a-pp))   ;and unit vector
			;;first-phi-a is phi- when observation point at first, first-phi-b is phi+
			(first-phi-a (- (azimuthal-angle first cusp a-pp-1)))
			(travel (spherical-angle first second))) ;Distance to be integrated
      (* travel			  ;We want to integrate dx, but the actual integral below is dlambda with lambda = 0...1
	 (qtrap
	  #'(lambda (x)			;Fraction of the way from first toward second along great circle connecting them
	      (let* ((position (spherical-interpolation first second x))
		     (phi (spherical-triangle-angle cusp first position)) ;Azimuthal angle of given position
		     (theta (spherical-angle cusp position))
		     (chi (if (> x 0.5)	;Get angle at position
			     (spherical-triangle-angle position first cusp)
			   ;;If x is very small, the above will try to get the angle when one side is tiny, so instead
			   ;;compute via supplementary angle
			   (- pi (spherical-triangle-angle position second cusp)))))
;;		(break "x = ~S phi- = ~S, phi+ = ~S" x (- first-phi-a phi) (- first-phi-b phi))
		(* (sin chi)		 ;See notes
		   (/ theta (sin theta)) ;This has no effect at first order in theta, but is more accurate for testing
		   (- (angular-power-distribution-from-cusp-1 a-pp-length b-pp-length ;theta*dP/dOmega, all n
							      (- first-phi-a phi) (- first-phi-b phi))
		      (if (and bulk-n (plusp bulk-n))
			  (integrated-angular-power-distribution-from-cusp-n ;Subtract contribution of first n
			   a-pp-length b-pp-length (- first-phi-a phi) (- first-phi-b phi) theta bulk-n)
			0.0)))))
	  0.0 1.0
	  :eps 1e-4			;Can't do better considering accuracy of cache
	  )))))

;;Add spectrum, computed with the bulk code, to data.
(defun do-triangle-spectrum (a-amplitudes b-amplitudes triangle data bin-size) 
  (loop with center = (triangle-center triangle)
	with spectrum = (smooth-gravitational-spectrum-1 a-amplitudes b-amplitudes center
						   ;;Must compute even harmonics we won't use to avoid aliasing
						   :max-harmonics (expt bin-size (length data)))
	with area = (apply #'spherical-triangle-area triangle)
	for i from 1 below (min (length spectrum) (expt bin-size (length data)))
	do (incf (aref data (floor (log i bin-size))) (* area (aref spectrum i)))))

;;Add cusp radiation to spectrum being assembled
(defun do-cusp-triangle-spectrum (a-amplitudes b-amplitudes triangle cusp-list data bin-size)
  (let* ((bulk-n (cusp-bulk-max-n cusp-list triangle))			  ;Max N to do by bulk code
	 (bulk-bins (if (zerop bulk-n) 0 (floor (log bulk-n bin-size))))) ;Number of bins to do by bulk
    (when (plusp bulk-bins)
      (loop with center = (triangle-center triangle)
	    with area = (apply #'spherical-triangle-area triangle)
	    for bin below bulk-bins
	    do (incf (aref data bin)
		     ;;Include N contained in bin-size^bin-0.5 ... bin-size^(bin+1)-0.5
		     (* area (loop for harmonic from (round (expt bin-size bin)) below (round (expt bin-size (1+ bin)))
				   sum (slow-smooth-gravitational-power a-amplitudes b-amplitudes center harmonic))))))
    (when (< bulk-bins (length data))
      (dolist (cusp-info cusp-list)
	(loop for (a b weight) in (cusp-triangle-decomposition cusp-info triangle)
	      do (do-cusp-triangle-spectrum-1 cusp-info a b weight data bin-size bulk-bins))))))

;;Add spectrum due to cusp to data.  If start-bin given, only fill bins starting with that one.
(defun do-cusp-triangle-spectrum-1 (cusp-info first second weight data bin-size &optional (start-bin 0))
  (format t "C") (force-output)
  (let ((cusp (cusp-info-direction cusp-info)))
    (when (spherical-triangle-degenerate-p cusp first second)
      (return-from do-cusp-triangle-spectrum-1 0.0))
    (mirror-image-let* ((a-pp (cusp-info-a-pp cusp-info))
			(a-pp-length (3vector-length a-pp)) ;Separate a'' and b'' into length
			(a-pp-1 (3vector-normalize a-pp))   ;and unit vector
			;;first-phi-a is phi- when observation point at first, first-phi-b is phi+
			(first-phi-a (- (azimuthal-angle first cusp a-pp-1)))
			(travel (spherical-angle first second))) ;Distance to be integrated
      (loop for bin from start-bin below (length data)
	    for omega1 = (if (zerop bin) 0.0 (* 4 pi (- (expt bin-size bin) 0.5))) ;cut between harmonics
	    for omega2 = (* 4 pi (- (expt bin-size (1+ bin)) 0.5))
	    do
;;	    (format t "Bin ~D~%" bin)
	    (incf (aref data bin)
		  (* travel weight
		     (integrate
		      #'(lambda (x)	;Fraction of the way from first toward second along great circle connecting them
			  (cusp-triangle-spectrum-2 x cusp first second a-pp-length b-pp-length
						    first-phi-a first-phi-b omega1 omega2))
		      0.0 1.0 :relative-error 1e-4)))
	    ))))

(defun cusp-triangle-spectrum-2 (x cusp first second a-pp-length b-pp-length
				   first-phi-a first-phi-b omega1 omega2)
  (mirror-image-let*
      ((position (spherical-interpolation first second x))
       (phi (spherical-triangle-angle cusp first position)) ;Azimuthal angle of given position
       (theta (spherical-angle cusp position))
       (chi (if (> x 0.5)		;Get angle at position
		(spherical-triangle-angle position first cusp)
	      ;;If x is very small, the above will try to get the angle when one side is tiny, so 
	      ;;instead compute via supplementary angle
	      (- pi (spherical-triangle-angle position second cusp))))
       (phi-a (- first-phi-a phi))	;phi+/-
       (f-a (/ (expt (abs (sin phi-a)) 3) 6 a-pp-length)))
    ;;			    (break "x = ~S phi- = ~S, phi+ = ~S, theta = ~S" x phi-a phi-b theta)
    (let* ((sff (sqrt (* f-a f-b)))
	   (a (/ f-b f-a)))
      (* (/ 48 (expt pi 2) (sqrt (abs (* a-pp-length b-pp-length (sin phi-a) (sin phi-b)))))
	 ;;See notes.  We are integrating dx but we need dphi = sin xi/sin theta dx
	 (/ (sin chi) (sin theta))
	 ;;			(cusp-wedge-spectrum-1 a (* sff omega1) (* sff omega2) ;q1, q2
	 ;; theta (plusp (* (sin phi-a) (sin phi-b))))
	 (cusp-wedge-spectrum a (* sff omega1) (* sff omega2) ;q1, q2
			      0.0 theta (plusp (* (sin phi-a) (sin phi-b)))))
      )))


;;Power in wedge with angle dphi and theta = 0...theta2 in frequencies between omega1 and omega2 is this times dphi,
;;times coefficient above, with q_i = sqrt{f+ f-} omega_i
(defun cusp-wedge-spectrum (a q1 q2 theta1 theta2 plusp)
  (declare (double-float a q1 q2 theta1 theta2))
  (when (< a 1) (setq a (/ 1 a)))	;Fold to >1.
  ;;We're trying to integrate (theta_max - theta_min) h(a,z) dz
  (let ((z0 (* q1 (expt theta1 3)))	;Minimum possible z.  From here to z1, theta_min = theta_1
	(z1 (* q2 (expt theta1 3)))	;Above here, theta_min = (z/q2)^(1/3).
	(z2 (* q1 (expt theta2 3)))	;Below here, theta_max = (z/q1)^(1/3)
	(z3 (* q2 (expt theta2 3))))	;Maximum possible z.  From z2 to here, theta_max = theta_2
    (-
     ;;First term of integrand, integrated from z0 to z3
     (+ (if (zerop q1) 0.0
	  (/ (- (h23pm a z2 plusp) (h23pm a z0 plusp)) ;z0..z2 with (z/q1)^(1/3) in integrand
	     (expt q1 1/3))) 
	(* (- (h2pm a z3 plusp) (h2pm a z2 plusp)) ;z2..z3 with theta_2
	   theta2))
     ;;Now subtract second term in integrand from z0 to z3
     (+ (* (- (h2pm a z1 plusp) (h2pm a z0 plusp)) ;z0...z1 with theta_1
	   theta1)
	(/ (- (h23pm a z3 plusp) (h23pm a z1 plusp)) ;z1...z3 with (z/q2)^(1/3).
	   (expt q2 1/3))))))
	      
;;Special case for theta1 = 0.
(defun cusp-wedge-spectrum-1 (a q1 q2 theta2 plusp)
  (declare (double-float a q1 q2 theta2))
  (when (< a 1) (setq a (/ 1 a)))	;Fold to >1.
  ;;We're trying to integrate (theta_max - theta_min) h(a,z) dz
  (let ((z2 (* q1 (expt theta2 3)))	;Below here, theta_max = (z/q1)^(1/3)
	(z3 (* q2 (expt theta2 3))))	;Maximum possible z.  From z2 to here, theta_max = theta_2
    (-
     ;;First term of integrand, integrated from 0 to z3
     (+ (if (zerop q1) 0.0
	  (/ (- (h23pm a z2 plusp) (h23pm a 0.0 plusp)) ;0..z2 with (z/q1)^(1/3) in integrand
	     (expt q1 1/3))) 
	(* (- (h2pm a z3 plusp) (h2pm a z2 plusp)) ;z2..z3 with theta_2
	   theta2))
     ;;Now subtract second term in integrand from 0 to z3
     (/ (- (h23pm a z3 plusp) (h23pm a 0.0 plusp)) ;0...z3 with (z/q2)^(1/3).
	(expt q2 1/3)))))

;;Power into a single mode times n^(4/3)
;;This works for n = double-float-positive-infinity to give asymptotic power
(defun smooth-gravitational-mode-scaled (a-amplitudes b-amplitudes n
					 &key (split-levels 0) (threshold1 0.1) (threshold2 (* 2 threshold1))
					 cusp-data)
  (let ((cusp-result 0)
	(non-cusp-result 0)
	(cusp-data (or cusp-data (and (plusp threshold1)
				      (find-smooth-cusps a-amplitudes b-amplitudes))))) ;Get cusps if any
    (spherical-integral-separate-cusps
     #'(lambda (triangle list)
	 (let ((center (triangle-center triangle)))
	   (if (and list					  ;Near a cusp?
		    (> n (cusp-bulk-max-n list triangle)))	  ;Small modes do by regular code anyway
	       (incf cusp-result				   ;Do cusp calculation, accumulate in cusp result
		     (cusp-triangle-mode-scaled triangle list n))	   ;This lacks n^(-4/3) factor
	     ;;Not doing cusp.  How many harmonics do we expect to contribute?
	     (let ((max-harmonics (smooth-gravitational-max-harmonic a-amplitudes b-amplitudes center)))
	       (when (< n max-harmonics) ;Ignore if too large
		 (incf non-cusp-result	;Do non-cusp calculation, accumulate in non-cusp result
		       (* (apply #'spherical-triangle-area triangle) ;Area * power at center
			  (slow-smooth-gravitational-power a-amplitudes b-amplitudes center n)
			  (expt n 4/3)))))))) ;Remove n^(-4/3) factor
     split-levels cusp-data threshold1 threshold2)
    (values (+ cusp-result non-cusp-result) non-cusp-result cusp-result)))

;;Give maximum harmonic where significant power could be present.  See radiation-bounds.tex
(defun smooth-gravitational-max-harmonic (a-amplitudes b-amplitudes omega &key (tolerance 1e-6))
  (mirror-image-let* ((a-rmax (smooth-gravitational-max-r a-amplitudes omega))
		      (a-factor (/ pi a-rmax))) ;a must decline by e^{-factor*n}
    (let ((factor (+ a-factor b-factor))			;total power declines according to this
	  (c (* 8 pi))						;Could increase for I_x + I_y and I^2 + II*
	  (n 100))
      ;;Power is given by C n^2 e^{-f n}.  Since f<<1, we're interested in n>>1, where the total power at n
      ;;and above is roughly C (n^2/f) e^{-f n}.  We want this to be less than T, so the limit is
      ;;given by n = (1/f) ln (n^2 C/(f T)).  We approximate this by iteration.
      (setq n (/ (log (/ (* n n c) factor tolerance)) factor))
      (setq n (/ (log (/ (* n n c) factor tolerance)) factor))
      (setq n (max (ceiling n) smooth-gravitational-min-max-harmonic))
      (if (= smooth-gravitational-max-harmonic-contrast 1) ;Save trouble if if we wouldn't do anything anyway
	  (values n n)
	(let ((n-a n)
	      (n-b n))
	  ;;We now compute separate requirements for a and b, demanding that
	  ;;the decline in the separate power should obey sqrt{C} n e^{-f_i n_i} = sqrt{T}, so 
	  ;;n_i = (1/f_i) ln (n_i sqrt{C/T}).  The purpose here is to prevent aliasing by making sure
	  ;;that the power in, say, a, has declined sufficiently before it can wrap back to be erroneously combined with
	  ;;large modes of b.
	  (mirror-images
	   (setq n-a (/ (log (* n (sqrt (/ c tolerance)))) a-factor))
	   (setq n-a (/ (log (* n (sqrt (/ c tolerance)))) a-factor))
	   ;;No more than max contrast times larger.  By that point we hope that the n^{-2/3} decline in I_n will
	   ;;be enough to prevent aliasing problems, even if there is no exponential fall off.
	   (setq n-a (ceiling (min (* n smooth-gravitational-max-harmonic-contrast) n-a)))
	   )
	  (when (and (> n-a n) (> n-b n))
	    (error "Both n-a and n-b were increased in ~A" 'smooth-gravitational-max-harmonic))
	  ;;Individual requirement is normally weaker (smaller n).  We only care in the case where it is stronger.
	  (values (max n n-a) (max n n-b)))))))


;;Spectrum in a given direction, an array indexed by harmonic number giving power
(defun smooth-gravitational-spectrum-1 (a-amplitudes b-amplitudes omega &key (max-harmonics 100000))
  ;;Adjust max-harmonics to the largest number which would give the same power of 2 in
  ;;smooth-gravitational-wave, by finding n/2 by the computation there and then setting
  ;;f = n(1-smooth-gravitational-wave-discard-fraction).
  (setq max-harmonics
	(floor (* (expt 2 (ceiling (log (/ max-harmonics (- 1 smooth-gravitational-wave-discard-fraction)) 2)))
		  (- 1 smooth-gravitational-wave-discard-fraction))))
  (multiple-value-bind (a-harmonics b-harmonics)
      (smooth-gravitational-max-harmonic a-amplitudes b-amplitudes omega)
    (let ((warn nil)
	  (original (list a-harmonics b-harmonics)))
      (mirror-images
       (when (> a-harmonics max-harmonics)
	 (setq a-harmonics max-harmonics
	       warn t)))
      (when warn
       (warn "Calculation suggests ~S harmonics needed, but we're using only ~D"
	     (if (= smooth-gravitational-max-harmonic-contrast 1) ;Must be same: simplify output
		    (first original) original)
	     max-harmonics)))
    (let* ((spectrum (smooth-gravitational-spectrum-2 a-amplitudes b-amplitudes omega a-harmonics b-harmonics))
	   (harmonics (length spectrum))) ;Use what we got
      (loop with warned = nil		 ;Check for problems
	    for harmonic below harmonics
	    for power = (aref spectrum harmonic)
	    when (and (> (* harmonic 10) (* harmonics 9)) ;Last 10%
		      (> power 1d-7)		     ;Power should be small
		      (not warned))		     ;Haven't complained already
	    do (if (> power 1d-3)
		   (cerror "Continue without further checking for this Omega"
			   "Power in harmonic ~D is ~S, which is much larger than predicted" harmonic power)
		 (warn "Power in harmonic ~D is ~S, which is larger than predicted" harmonic power))
	    and do (setq warned t))		;Don't warn again for same omega
      spectrum)))

;;Power in the first N harmonics
(defun smooth-gravitational-power-first-n (a-amplitudes b-amplitudes omega n)
  (loop for harmonic from 1 to n
	sum (slow-smooth-gravitational-power a-amplitudes b-amplitudes omega harmonic)))

;;Total power and rate of momentum radiation in all directions as a 4-vector
;;Directions within threshold1 distance of a cusp direction will be handled by the cusp code.
;;If a triangle is within threshold2 of some cusp, that cusp will be included if we are
;;already handling some other cusp for this triangle.
(defun smooth-gravitational-radiation (a-amplitudes b-amplitudes &key (max-harmonics 10000) (split-levels 0)
					  (threshold1 0.1) (threshold2 (* 2 threshold1)))
  (let ((cusp-result zero-3vector)
	(non-cusp-result zero-3vector)
	(cusp-data (and (plusp threshold1)
			(find-smooth-cusps a-amplitudes b-amplitudes)))) ;Get cusps if any
    (spherical-integral-separate-cusps
     #'(lambda (triangle list)
	 (let* ((center (triangle-center triangle))
		(direction (3to4vector center 1.0))) ;4vector direction of emission for momentum-energy emitted
	   (if list			;Near a cusp?
	       (setq cusp-result	;Do cusp calculation, accumulate in cusp result
		     (4vector+ cusp-result
			       (4vector-scale direction
					      (cusp-bulk-triangle-integral a-amplitudes b-amplitudes list triangle))))
	     (setq non-cusp-result		;Do non-cusp calculation, accumulate in non-cusp result
		   (4vector+ non-cusp-result
			     (4vector-scale direction
					    (* (apply #'spherical-triangle-area triangle) ;Area * power at center
					       (reduce #'+ (smooth-gravitational-spectrum-1 a-amplitudes b-amplitudes center
										   :max-harmonics max-harmonics)))))))))
     split-levels cusp-data threshold1 threshold2)
    (values (4vector+ cusp-result non-cusp-result)
	    non-cusp-result
	    cusp-result)))

(defun functional-gravitational-power (a-function b-function frequencies &key (split-levels 0))
  (spherical-integral
   #'(lambda (omega)
       (reduce #'+
	       (functional-gravitational-spectrum-1 a-function b-function omega frequencies)))
   split-levels))

(defun smooth-gravitational-spectrum (a-amplitudes b-amplitudes
						  &key (harmonics 10000) (split-levels 0)
						  (threshold1 0.1) (threshold2 (* 2 threshold1))
						  (bin-size 2.0)
						  &allow-other-keys)
  (setq harmonics (expt 2 (ceiling (log harmonics 2)))) ;Round up
  (let* ((bins (ceiling (log harmonics bin-size)))
	 (data (make-array bins :initial-element 0.0))
	 (cusp-data (and (plusp threshold1)
			 (find-smooth-cusps a-amplitudes b-amplitudes)))) ;Get cusps if any
    (spherical-integral-separate-cusps
     #'(lambda (triangle list)		;Add spectrum to data
	 (if list
	     (do-cusp-triangle-spectrum a-amplitudes b-amplitudes triangle list data bin-size)
	   (do-triangle-spectrum a-amplitudes b-amplitudes triangle ;No cusp.  Do regular calculation
				 data bin-size)))
     split-levels cusp-data threshold1 threshold2)
    data))

(defun plot-smooth-gravitational-spectrum (a-amplitudes b-amplitudes &rest keys &key (bin-size 2.0) &allow-other-keys)
  (let ((data (apply #'smooth-gravitational-spectrum a-amplitudes b-amplitudes keys)))
    (apply #'gnuplot 1 (length data)
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (and (numberp point)
		    (values (expt bin-size (+ point 0.5)) (aref data point))))
	   :logscale '(:x :y)
	   :styles :linespoints
	   keys)
    (format t "Total power ~S~%" (reduce #'+ data))))

(defun smooth-gravitational-power (&rest args)
  (multiple-value-bind (all non-cusp cusp) (apply #'smooth-gravitational-radiation args)
    (values (3vector-t all) (3vector-t non-cusp) (3vector-t cusp))))

(defun test-smooth-gravitational-power (diamond &rest keys)
   (multiple-value-bind (a-hats a-sigmas b-hats b-sigmas) (get-ab-data-sigmas diamond)
     (apply #'smooth-gravitational-power
	    (amplitudes-from-hats (Lorentzian-smooth a-hats 0.1 :sigmas a-sigmas :n 64))
	    (amplitudes-from-hats (Lorentzian-smooth b-hats 0.1 :sigmas b-sigmas :n 64))
	    keys)))

(defun plot-smooth-gravitational-spectrum-1 (diamond smoothing-fraction omega frequencies
					       &rest keys &key harmonics (bin-size 2) &allow-other-keys)
  (multiple-value-bind (a-hats a-sigmas b-hats b-sigmas) (get-ab-data-sigmas diamond)
    (mirror-image-let ((a-amplitudes (amplitudes-from-hats ;Smooth data, resample
				      (Lorentzian-smooth a-hats smoothing-fraction :sigmas a-sigmas :n frequencies))))
      (unless harmonics
	(setq harmonics (expt 2 (ceiling (log (smooth-gravitational-max-harmonic a-amplitudes b-amplitudes omega) 2)))))
      (let* ((bins (ceiling (log harmonics bin-size)))
	     (data (make-array bins :initial-element 0.0))
	     (spectrum (smooth-gravitational-spectrum-2 a-amplitudes b-amplitudes omega harmonics harmonics)))
	(loop for i from 1 below harmonics
	      do (incf (aref data (floor (log i bin-size))) (aref spectrum i)))
	(apply #'gnuplot 1 bins
	       #'(lambda (plot point)
		   (declare (ignore plot))
		   (and (numberp point)
			(values (expt bin-size point) (aref data point))))
	       :logscale '(:x :y)
	       keys)))))

;;Plot I_n and limits.
;;Use rest frame loop.
(mirror-images
(defun plot-smooth-gravitational-wave-a (diamond smoothing-fraction omega frequencies harmonics &rest keys)
  (multiple-value-bind (hats sigmas) (get-a-data-sigmas diamond)
    (let* ((n (expt 2 (ceiling (1+ (log frequencies 2))))) ;N = # of points = 2 * desired frequencies.  Round up.
	   (amplitudes (amplitudes-from-hats (Lorentzian-smooth hats smoothing-fraction :sigmas sigmas :n n))))
      (apply #' plot-smooth-gravitational-wave-1 amplitudes omega harmonics keys)))))

(defun plot-smooth-gravitational-wave-1 (amplitudes omega harmonics &rest keys
					      &key (bin-size 2) &allow-other-keys)
  (let* ((bins (ceiling (log harmonics bin-size)))
	 (data (make-array bins :initial-element 0.0))
	 (amax (maximum-second-derivative amplitudes))
	 (biggest-cos (closest-approach amplitudes omega))
	 (hmin (- 1 biggest-cos))					 ;Minimum rate of change of phase function
	 (factor (/ (* pi (expt hmin 2)) amax))				 ;See below
	 (iperp (smooth-gravitational-wave amplitudes omega (* 2 harmonics)))) ;Don't use ixy
    (format t "Amax = ~S, Hmin = ~S, Factor = ~S~%" amax hmin factor)
    (loop for i from 1 below harmonics
	  do (incf (aref data (floor (log i bin-size))) (aref iperp i)))
    (apply #'gnuplot 2 bins
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (nth plot '("total I^2 in bin" "limit"))
		 (values (expt bin-size point) ;x coordinate
			 (if (zerop plot)
			     (aref data point)
			   ;;|I_perp| < exp(- 2 pi n hmin^2 / L a''max).  We square this and integrate over
			   ;;a bin from (bin-size)^n to (bin-size)^{n+1}
			   (/ (- (exp (* -1 factor (expt bin-size point)))
				 (exp (* -1 factor (expt bin-size (1+ point)))))
			      factor)))))
	   :logscale '(:x :y)
	   keys)))
	    
(defun smooth-gravitational-max-r (amplitudes omega &optional n)
  (setq n (max (if n (max n (vector-amplitudes-length amplitudes)) ;Always as many points as real amplitudes.
		 (* 2 (vector-amplitudes-length amplitudes)))	   ;But default is twice that
	       ;;At least some reasonable minimum to trace out the string.
	       ;;Actually doing some minimalization would really be better.
	       256))
  (setq n (expt 2 (ceiling (log n 2))))	       ;Round up to power of 2
  (let ((az-amplitudes (make-array n :element-type 'double-float)))
    (dotimes (j (vector-amplitudes-length amplitudes))		    ;Make amplitudes for a'_z = a' . Omega
      (setf (aref az-amplitudes j) (3vector-dot (vector-amplitudes-ref amplitudes j) omega)))
    (let* ((appz-amplitudes (differentiate-amplitudes az-amplitudes))      ;Amplitudes for a''_z
	   (az (discrete-transform-amplitudes az-amplitudes n))				 ;N samples of a'_z
	   (appz (discrete-transform-amplitudes appz-amplitudes n)))		 ;N samples of a''_z
      (loop for index below n
	    maximize (/ (aref appz index) (expt (- 1.0 (aref az index)) 2))) ;maximize r = a''_z/h^2
      )))

(defun plot-smooth-gravitational-power-projected (a-amplitudes b-amplitudes steps &rest keys
							 &key (max-harmonics 100000)
							 (threshold1 0.1)
							 (threshold2 (* 2 threshold1))
							 (cusp-data (find-smooth-cusps a-amplitudes b-amplitudes))
							 (cusp-bulk-max-n cusp-bulk-max-n)
							 &allow-other-keys)
  (apply #'plot-projected-angular-distribution
	 #'(lambda (direction)
	     (smooth-gravitational-power-1 direction a-amplitudes b-amplitudes :max-harmonics max-harmonics
				     :cusp-data cusp-data :threshold1 0.1 :threshold2 threshold2))
	 steps
	 keys))

;;This gives a sphere in 3d.
(defun plot-smooth-gravitational-power (a-amplitudes b-amplitudes theta-steps
					       &rest keys
					       &key (max-harmonics 100000)
					       (phi-steps (* 2 theta-steps))
					       (threshold1 0.1)
					       (threshold2 (* 2 threshold1))
					       (cusp-data (find-smooth-cusps a-amplitudes b-amplitudes))
					       cbmax ;Value at which to saturate the color scale 
					       (cusp-bulk-max-n cusp-bulk-max-n)
					       interpolate
					       &allow-other-keys)
  (apply #'gnuplot 1 (* (1+ theta-steps) (+ 2 phi-steps))
	 #'(lambda (plot point)
	     (declare (ignore plot))
	     (when (numberp point)
	       ;;theta-step = 0...step
	       ;;phi-step = 0...phi-steps + 1, to leave room for blank record
	       (multiple-value-bind (theta-step phi-step) (floor point (+ 2 phi-steps))
		 (if (> phi-step phi-steps) ;end of circle of fixed theta
		     (progn (format t "~D " theta-step) (force-output) nil)
		   (let* ((theta (* theta-step (/ pi theta-steps)))
			  (phi (* phi-step (/ (* 2 pi) phi-steps)))
			  (plot-direction (spherical-coordinates theta phi));This is the direction given to gnuplot
			  ;;Now back off half a step in both theta and phi directions.  This gives the center of the
			  ;;preceding quadrilateral, but because of "corners2color c2" this is the region that
			  ;;will be given the specified color
			  (direction (spherical-coordinates (- theta (/ pi theta-steps 2)) (- phi (/ pi phi-steps)))))
		     (values-list
		      (append
		       (3vector-list plot-direction)
		       (list (smooth-gravitational-power-1 direction a-amplitudes b-amplitudes :max-harmonics max-harmonics
						     :cusp-data cusp-data :threshold1 0.1 :threshold2 threshold2)))))))))
	 :3d t
	 :styles :pm3d
	 :prelude (format nil "~@[set cbrange [0:~S]~]~%set pm3d depthorder corners2color c2 ~@[interpolate ~D~:*,~D~]
unset border~%unset tics~%set cbtics~%unset lmargin~%unset rmargin~%unset tmargin~%unset bmargin~%set view equal xyz~%"
			  cbmax interpolate)
	 keys))

(defun smooth-gravitational-power-1 (direction a-amplitudes b-amplitudes
					 &key (max-harmonics 100000) (threshold1 0.1) (threshold2 (* 2 threshold1))
					 cusp-data
					 &allow-other-keys)
  (if (loop for info in cusp-data
	    thereis (< (spherical-angle direction (cusp-info-direction info))
		       threshold1))		     ;Close to any cusp?  Use cusp code
      (let* ((cusp-list (loop for info in cusp-data ;Close to this cusp with weaker threshold?
			      when (< (spherical-angle direction (cusp-info-direction info)) threshold2)
			      collect info))
	     (bulk-max (cusp-bulk-direction-max-n cusp-list direction)))
	(+ (if (plusp bulk-max)
	       (smooth-gravitational-power-first-n a-amplitudes b-amplitudes direction bulk-max)
	     0)
	   (loop for info in cusp-list
		 sum (/ (angular-power-distribution-from-cusp info direction bulk-max) ;This is power*theta
			(spherical-angle direction (cusp-info-direction info))))))
    ;;Not near a cusp
    (reduce #'+ (smooth-gravitational-spectrum-1 a-amplitudes b-amplitudes direction :max-harmonics max-harmonics))))
			
#| This version colors each triangle based on its power density computed by adding up the cubically spaced modes

(defun plot-smooth-gravitational-power (a-amplitudes b-amplitudes
					       &rest keys
					       &key (max-harmonics 100000)
					       (threshold1 0.1)
					       (threshold2 (* 2 threshold1))
					       (split-levels 0)
					       cusp-data
					       cbmax ;Value at which to saturate the color scale 
					       &allow-other-keys)
  (unless cusp-data (setq cusp-data (find-smooth-cusps a-amplitudes b-amplitudes)))
  (let* ((triangulation (triangulate-sphere-separate-cusps split-levels cusp-data threshold1 threshold2))
	 (modes (cubically-spaced-modes 1000 t))
	 (data (make-array (length modes) :element-type 'double-float))
	 (last-triangle nil)
	 last-power)
    (format t "~&~D triangles: " (length triangulation))
    (apply #'gnuplot 1 (* 7 (length triangulation))
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (when (numberp point)
		 (multiple-value-bind (index step) (floor point 7)
		   (when (and (zerop step) (zerop (mod index 10))) (format t "~D " index) (force-output))
		   (fill data 0.0)
		   (destructuring-bind (triangle . cusp-list) (aref triangulation index)
		     (unless (eq triangle last-triangle)
		       (setq last-power ;Cache power over the different corners of the triangle
			     (/ (smooth-gravitational-triangle-power-1 a-amplitudes b-amplitudes triangle cusp-list
								 :modes modes :array data :max-harmonics max-harmonics)
				(apply #'spherical-triangle-area triangle)))
		       (setq last-triangle triangle))
		     (flet ((4d (v) (values-list (append (3vector-list v) (list last-power)))))
		       (case step
			 (0 (4d (first triangle)))
			 (1 (4d (second triangle)))
			 (2 nil)			 ;Break between scan lines
			 ((3 4) (4d (third triangle)))	 ;twice because gnuplot wants rectangle
			 ((5 6) nil)))))))		 ;double break between triangles
	   :3d t
	   :styles :pm3d
	   :prelude (format nil "~@[set cbrange [0:~S]~]~%set pm3d depthorder~%unset border~%unset tics~%set cbtics
unset lmargin~%unset rmargin~%unset tmargin~%unset bmargin~%" cbmax)
	   keys)))
|#

;;Actual mode power including n^(-4/3) factor
(defun smooth-gravitational-mode (a-amplitudes b-amplitudes n &rest keys)
  (* (expt n -4/3) (apply #'smooth-gravitational-mode-scaled a-amplitudes b-amplitudes n keys)))

;;Return an array of mode power corresponding to the mode numbers given in the sequence MODES, each
;;multiplied by n^(4/3)
;;This works for n = double-float-positive-infinity to give asymptotic power
(defun smooth-gravitational-modes-scaled (a-amplitudes b-amplitudes modes
					 &key (split-levels 0) (threshold1 0.1) (threshold2 (* 2 threshold1))
					 cusp-data (max-harmonics 100000))
  (let ((cusp-data (or cusp-data (and (plusp threshold1)
				      (find-smooth-cusps a-amplitudes b-amplitudes)))) ;Get cusps if any
	(result (make-array (length modes) :element-type 'double-float :initial-element 0.0))) 
    (spherical-integral-separate-cusps
     #'(lambda (triangle list)
	 (smooth-gravitational-modes-scaled-1 a-amplitudes b-amplitudes triangle list modes result max-harmonics))
     split-levels cusp-data threshold1 threshold2)
    result))

;;Restore n^(-4/3) factors
(defun smooth-gravitational-modes (a-amplitudes b-amplitudes modes &rest keys)
  (let ((result (apply #'smooth-gravitational-modes-scaled a-amplitudes b-amplitudes modes keys)))
    (loop for harmonic being the elements of modes
	  for index from 0
	  do (setf (aref result index) (* (aref result index) (expt harmonic -4/3))))
    result))
		
;;Add the given modes to the array DATA
(defun smooth-gravitational-modes-scaled-1 (a-amplitudes b-amplitudes triangle list modes data max-harmonics)
  (let ((center (triangle-center triangle))
	(area (apply #'spherical-triangle-area triangle)))
    (if list				;Near a cusp?
	(loop with bulk-max = (cusp-bulk-max-n list triangle)
	      for harmonic being the elements of modes
	      for index from 0
	      do
	      (incf (aref data index)
		    (if (<= harmonic bulk-max) ;Do small modes by Fourier integration at center	
			(* area (slow-smooth-gravitational-power a-amplitudes b-amplitudes center harmonic)
			   (expt harmonic 4/3))					    ;Rescale to remove n^(-4/3) factor
		      (cusp-triangle-mode-scaled triangle list harmonic)))) ;total power in region, convert to density
      (let ((spectrum (smooth-gravitational-spectrum-1 a-amplitudes b-amplitudes center :max-harmonics max-harmonics)))
	(loop for harmonic being the elements of modes
	      for index from 0
	      while (< harmonic (length spectrum)) ;If this harmonic was calculated
	      do (incf (aref data index) (* area (aref spectrum harmonic)
					    (expt harmonic 4/3)) ;Rescale to remove n^(-4/3) factor
		       ))))))

;;Power in a single triangle.  MODES is the set of modes to use, and ARRAY is an array to collect data in.
(defun smooth-gravitational-triangle-power-1 (a-amplitudes b-amplitudes triangle cusp-list
						     &key modes array (max-harmonics 100000))
  (unless modes (setq modes (cubically-spaced-modes 1000 t)))
  (unless array (make-array (length modes) :element-type 'double-float))
  (fill array 0.0)
  (smooth-gravitational-modes-scaled-1 a-amplitudes b-amplitudes triangle cusp-list modes array max-harmonics)
  (numerical-sum-43 modes array))

;;Power in a triangle in single mode, times n^(4/3)
;;This works for n = double-float-positive-infinity to give asymptotic power
(defun cusp-triangle-mode-scaled (triangle cusp-list n)
  (loop for cusp-info in cusp-list sum
    (loop for (first second weight) in (cusp-triangle-decomposition cusp-info triangle)
	  sum (let ((cusp (cusp-info-direction cusp-info)))
		(if (spherical-triangle-degenerate-p cusp first second)
		    0.0
		  (mirror-image-let* ((a-pp (cusp-info-a-pp cusp-info))
				      (a-pp-length (3vector-length a-pp)) ;Separate a'' and b'' into length
				      (a-pp-1 (3vector-normalize a-pp))   ;and unit vector
				      ;;first-phi-a is phi- when observation point at first, first-phi-b is phi+
				      (first-phi-a (- (azimuthal-angle first cusp a-pp-1)))
				      (travel (spherical-angle first second))) ;Distance to be integrated
		    (* travel weight
		       (qtrap
			#'(lambda (x)	;Fraction of the way from first toward second along great circle connecting them
			    (cusp-triangle-mode-scaled-2 x cusp first second a-pp-length b-pp-length
						  first-phi-a first-phi-b n))
			0.0 1.0 :eps 1e-4)))
		  )))))

;;Power in a single mode integrated over a wedge is this times dphi times n^(-4/3)
;;This works for n = double-float-positive-infinity
(defun cusp-triangle-mode-scaled-2 (x cusp first second a-pp-length b-pp-length first-phi-a first-phi-b n)
  (mirror-image-let*
      ((position (spherical-interpolation first second x))
       (phi (spherical-triangle-angle cusp first position)) ;Azimuthal angle of given position
       (theta (spherical-angle cusp position))
       (chi (if (> x 0.5)		;Get angle at position
		(spherical-triangle-angle position first cusp)
	      ;;If x is very small, the above will try to get the angle when one side is tiny, so 
	      ;;instead compute via supplementary angle
	      (- pi (spherical-triangle-angle position second cusp))))
       (phi-a (- first-phi-a phi))	;phi+/-
       (f-a (/ (expt (abs (sin phi-a)) 3) 6 a-pp-length)))
;;    (format t "x = ~S, phi- = ~S, phi+ = ~S, theta = ~S, a-pp = ~S, b-pp = ~S" x phi-a phi-b theta a-pp-length b-pp-length)
    (let* ((sff (sqrt (* f-a f-b)))
	   (a (/ f-b f-a)))
      (* (/ (* 16 (expt 1.5 1/3))
	    (expt pi 7/3)		;Don't divide by factor (expt n 4/3)
	    (expt (* a-pp-length b-pp-length) 1/3)
	    (abs (* (sin phi-a) (sin phi-b))))
	 ;;See notes.  We are integrating dx but we need dphi = sin xi/sin theta dx
	 (/ (sin chi) (sin theta))
	 (if (= n double-float-positive-infinity)
	     (h3pm a (plusp (* (sin phi-a) (sin phi-b))))
	   (h23pm a (* sff (expt theta 3) 4 pi n) (plusp (* (sin phi-a) (sin phi-b)))))))))

;;Return an array of mode numbers consisting of (n/j)^3 for j=n...1, with duplicates eliminated
;;There are about (n/3)^(3/4) consecutive modes, followed by about 3^(1/4) n^(3/4) more modes with
;;progressively larger spacing, for a total of about 4 (n/3)^(3/4).
;;If requested, include infinity
(defun cubically-spaced-modes (n infinity-p)
  (coerce (loop with last = 0
		for j from n downto 0
		for this = (if (zerop j)
			       (and infinity-p double-float-positive-infinity)
			     (round (expt (/ (double-float n) (double-float j)) 3)))
		when (and this (/= this last))
		collect this and do (setq last this))
	  'vector))
		
(defun z43 (n)
  (if (= n double-float-positive-infinity) 0.0
    (gsl:hurwitz-zeta (/ 4.0 3.0) (double-float n))))
(defun z53 (n)
  (if (= n double-float-positive-infinity) 0.0
    (gsl:hurwitz-zeta (/ 5.0 3.0) (double-float n))))

(defconstant n13 (/ -1.0 3.0))

;;Give the coefficient of F_this when summing using samples previous, this, next
;;appropriate to series that decline at least as fast as n^{-4/3}.  See spectrum-today.tex
;;previous or next should be NIL at start and end of sequence respectively
;;Elements can be double-float-positive-infinity
(defun numerical-sum-coefficient-43 (previous this next)
  (- (if next
	 (/ (- (z53 this) (z53 next) (* (- (z43 this) (z43 next)) (expt next n13)))
	    (- (expt this n13) (expt next n13)))
       0)
     (if previous
	 (/ (- (z53 previous) (z53 this) (* (- (z43 previous) (z43 this)) (expt previous n13)))
	    (- (expt previous n13) (expt this n13)))
       0)))

;;Approximate the sum of n^(-4/3) Sn for n = 1...infinity from the given Sn values
(defun numerical-sum-43 (modes Sn)
  (loop with length = (length modes)
	for index below length
	sum (* (numerical-sum-coefficient-43 (and (plusp index) (elt modes (1- index))) (elt modes index)
					     (and (< index (1- length)) (elt modes (1+ index))))
	       (aref Sn index))))


;;Testing

;;Return vector amplitudes for a' and b' for a Burden loop.
(defun make-smooth-burden-loop (m n psi)
  (let ((aft (make-array 3))
	(bft (make-array 3))
	(nf (expt 2 (ceiling (log (1+ (max m n)) 2))))) ;Number of frequencies to use in arrays
    (dotimes (component 3)
      (setf (aref aft component) (make-array (* 2 nf) :element-type 'double-float :initial-element 0.0))
      (setf (aref bft component) (make-array (* 2 nf) :element-type 'double-float :initial-element 0.0)))
    ;;Now fill in non-zero components
    ;;The sign of the imaginary term of a is opposite Burden's paper, because we use a(v) = a(t-sigma), whereas he uses
    ;;a(xi) = a(sigma-t).
    ;;Thus psi=0 is the breather, psi=pi the rotator
    (setf (aref (aref aft 2) (* 2 m)) 0.5		       ;A_m_z = 0.5
	  (aref (aref aft 0) (1+ (* 2 m))) -0.5		       ;A_m_x = -0.5i
	  (aref (aref bft 2) (* 2 n)) 0.5		       ;B rotated from A by angle psi
	  (aref (aref bft 0) (1+ (* 2 n))) (* 0.5 (cos psi))   ;B_n_x
	  (aref (aref bft 1) (1+ (* 2 n))) (* 0.5 (sin psi)))   ;B_n_y
    (values aft bft)))

;;Fourier amplitudes for the Kibble-Turok loop. ALPHA sets the perturbation of a' away from circular, while PHI
;;sets the angle between a' and b'.
(defun make-smooth-kt-loop (alpha phi)
  (let ((aft (make-array 3))
	(bft (make-array 3)))
    (dotimes (component 3)
      (setf (aref aft component) (make-array 8 :element-type 'double-float :initial-element 0.0))
      (setf (aref bft component) (make-array 4 :element-type 'double-float :initial-element 0.0)))
    ;;Now fill in non-zero components
    (setf (aref (aref aft 0) 2) (* 0.5 (- 1 alpha)) ;A_1_x = 0.5 (1-alpha)
	  (aref (aref aft 1) 3) (* 0.5 (- 1 alpha)) ;A_1_y = 0.5i (1-alpha)
	  (aref (aref aft 2) 3) (sqrt (* alpha (- 1 alpha))) ;A_1_z = sqrt{alpha(1-alpha)}
	  (aref (aref aft 0) 6) (* 0.5 alpha)	       ;A_3_x = 0.5 alpha
	  (aref (aref aft 1) 7) (* 0.5 alpha)	       ;A_3_y = 0.5i alpha
	  ;;B rotated from x-y plane by phi around x-axis
	  (aref (aref bft 0) 2) 0.5		   ;B_1_x = 0.5
	  (aref (aref bft 1) 3) (* 0.5 (cos phi))  ;B_1_y = 0.5i cos phi
	  (aref (aref bft 2) 3) (* 0.5 (sin phi))) ;B_1_z = 0.5i sin phi
    (values aft bft)))

;;#|

;;Check radiation from a Burden loop against Burden's paper
(defun check-burden-radiation (m n multiple psi theta phi &optional (max-n 10000))
  (values (burden-radiation-analytic m n multiple psi theta phi)
	  (burden-radiation-numerical m n multiple psi theta phi max-n)))

;;Check find-smooth-cusp with a Burden loop
(defun check-find-smooth-cusp (m n psi xa xb)
  (multiple-value-bind (aft bft) (make-smooth-burden-loop m n psi)
    (find-smooth-cusp aft bft xa xb)))

;;dP/d Omega in units of G mu^2, according to Burden's paper, using above code
;;max-n is the number of points (and thus max-n/2 is the number of frequencies) to use in
;;smooth-gravitational-wave.  Even if the harmonic number m*n*multiple is small, this needs to be large to keep aliased power
;;from affecting your harmonic
(defun burden-radiation-numerical (m n multiple psi theta phi max-n)
  (let ((harmonic (* m n multiple))
	(omega (make-3vector (* (sin theta) (cos phi)) (* (sin theta) (sin phi)) (cos theta)))) ;observation direction
    (multiple-value-bind (aft bft) (make-smooth-burden-loop m n psi)
      (aref (smooth-gravitational-spectrum-2 aft bft omega (ceiling max-n 2) (ceiling max-n 2))
	    harmonic))))

;;dP_j/dOmega in units of G mu^2, according to Burden's paper, with j=m*n*multiple.
;;a is in the x-z plane, b in the plane making angle psi from x toward y
(defun burden-radiation-analytic (m n multiple psi theta phi)
  (let* ((pp (- psi phi))
	 (a (sqrt (+ (expt (cos theta) 2) (* (expt (sin theta) 2) (expt (cos phi) 2)))))
	 (b (sqrt (+ (expt (cos theta) 2) (* (expt (sin theta) 2) (expt (cos pp) 2)))))
	 (p (* n multiple))
	 (q (* m multiple))
	 (s1 (/ (* (sin theta) (sin phi)) (sqrt (- 1 (* (expt (sin theta) 2) (expt (sin phi) 2))))))
	 (s2 (/ (* (sin theta) (sin pp)) (sqrt (- 1 (* (expt (sin theta) 2) (expt (sin pp) 2))))))
	 (jp (besselj p (* p a)))
	 (jq (besselj q (* q b)))
	 (jpprime (/ (- (besselj (1- p) (* p a)) (besselj (1+ p) (* p a))) 2))
	 (jqprime (/ (- (besselj (1- q) (* q b)) (besselj (1+ q) (* q b))) 2)))
;;    (format t "a = ~S, b = ~S, s1 = ~S, s2 = ~S, jp = ~S, jq = ~S, jpprime = ~S, jqprime = ~S~%"
;;	    a b s1 s2 jp jq jpprime jqprime)
    (* 8 pi (expt (* m n multiple) 2)
       (+ (expt (+ (* s1 s2 jp jq) (* jpprime jqprime)) 2)
	  (expt (+ (* s1 jp jqprime) (* s2 jpprime jq)) 2)))))

;;Total power in a given harmonic
;;Alternatively, if you supply theta-max, you can get the power just around one cusp.
;;For high harmonics, integration may miss the cusp because it is so narrow, unless you use a smaller
;;theta-max.
(defun burden-power-analytic (m n multiple psi &key (phi-offset 0.01) (theta-max pi))
  (integrate #'(lambda (theta)
		 (* (sin theta)
		    (integrate #'(lambda (phi) (burden-radiation-analytic m n multiple psi theta phi))
			       phi-offset (+ (* 2 pi) phi-offset)))) ;Avoid problems at phi=0 etc.
	     0.0 theta-max))

;;Integral of power over a wedge near one of the two cusps divided by wedge angle.
;;phi-a is the angle from the observer direction to a''
(defun burden-power-analytic-phi (m n multiple psi phi-a theta-max)
  (integrate #'(lambda (theta)		;integrate only over theta
		 (* (sin theta)
		    (burden-radiation-analytic m n multiple psi theta (- phi-a))))
	     0.0 theta-max))

;;Total power in all harmonics.  Multiple-count is the n for cubically-space-modes
;;Sum over harmonics first, integrate afterward
;;The integrand then diverges in the cusp directions, but in in an integrable way
(defun burden-total-power-analytic-sum-first (m n psi multiple-count)
  (let* ((multiples (cubically-spaced-modes multiple-count nil)) ;I haven't computed the infinite limit
	 (Sn (make-array (length multiples) :element-type 'double-float)))
    (integrate #'(lambda (theta)	;integrate over solid angle
		   (format t ".") (force-output)
		   (* (sin theta)
		      (integrate #'(lambda (phi)
				     (loop for index below (length multiples)
					   for multiple = (aref multiples index)
					   do (setf (aref Sn index)
						    (* (burden-radiation-analytic m n multiple psi theta phi)
						       (expt multiple 4/3))))
				     (numerical-sum-43 multiples Sn))
				 0.01 (+ (* 2 pi) 0.01)))) ;Avoid problem at phi=pi etc.
	       0.0 pi)))

;;Total power in all harmonics.  Multiple-count is the n for cubically-space-modes
;;Integrate each harmonic separately, then sum
(defun burden-total-power-analytic (m n psi multiple-count)
  (loop with multiples = (cubically-spaced-modes multiple-count nil) ;I haven't computed the infinite limit
	with Sn = (make-array (length multiples) :element-type 'double-float)
	for index below (length multiples)
	for multiple = (aref multiples index)
	do (setf (aref Sn index) (* (burden-power-analytic m n multiple psi) (expt multiple 4/3)))
	do (format t "~D " multiple) (force-output)
	finally (return (numerical-sum-43 multiples Sn))))


(defun plot-burden-radiation (m n multiple psi theta-max phi
				&rest gnuplot-keys &key (points 100) &allow-other-keys)
  (apply #'gnuplot 2 points
	   #'(lambda (plot point)
	       (if (eq point :title) (nth plot '("analytic" "cusp approximation"))
		 (let ((theta (* (1+ point) (/ theta-max points))))
		   (values
		    theta
		    (ecase plot
		      (0 (burden-radiation-analytic m n multiple psi theta phi))
		      (1 (let ((direction (make-3vector (* (sin theta) (cos phi))
							(* (sin theta) (sin phi)) (cos theta))))
			   (mode-angular-power-distribution-from-cusp
			    (list nil nil (make-3vector 0.0 0.0 1.0) (make-3vector (- (* 2 pi)) 0.0 0.0)
				  (make-3vector 0.0 (* 2 pi) 0.0))
			    direction
			    (* multiple m n)))))))))
	   gnuplot-keys))
;; |#

;;Make a 1,1 Burden-like loop with gaps
(defun make-gap-loop-hats (gap n)
  (let ((k (* (- pi gap) 2))		;Total length of circular segments, 2 pi - 2 gaps
	(a-hats (make-array n))
	(b-hats (make-array n)))
    (loop for index from 0 below n
	  for sigma = (/ (+ index 0.5) n) ;Even spacing around gaps
	  for angle = (+ (* sigma k) (* gap (if (< sigma 0.5) 0.5 1.5)))
	  do (setf (aref a-hats index) (make-3vector (sin angle) 0.0 (cos angle))) ;xz-plane, gaps at z=1,-1
	  do (setf (aref b-hats index) (make-3vector 0.0 (sin angle) (cos angle)))) ;yz-plane, gaps at z=1,-1
    (values a-hats b-hats)))
	  
;;Two-speed loop: Travel at speed 2 pi s in gaps of size g.
(defun make-two-speed-loop-hats (gap-angle speedup n)
  (let* ((gap-speed (* 2 pi speedup))
	 (gap-time (/ gap-angle gap-speed))	;Time to cross one gap
	 (non-gap-angle (- pi gap-angle))	;Angle of one non-gap
	 (non-gap-time (- 0.5 gap-time))	;Time to cross one non-gap
	 (non-gap-speed (/ non-gap-angle non-gap-time))
	 (a-hats (make-array n))
	 (b-hats (make-array n))
	 (n2 (floor n 2)))
    (loop for index from 0 below n
	  for split = (mod index n2)		 ;0..n/2 twice
	  for sigma = (/ (double-float split) n) ;0..1/2
	  for angle = (+ (cond ((< sigma (/ gap-time 2)) (* gap-speed sigma))
			       ((> sigma (- 0.5 (/ gap-time 2))) (- pi (* gap-speed (- 0.5 sigma))))
			       (t (+ (/ gap-angle 2) (* (- sigma (/ gap-time 2)) non-gap-speed))))
			 (if (>= index n2) pi 0)) ;Second half same as first
	  do (setf (aref a-hats index) (make-3vector (sin angle) 0.0 (cos angle))) ;xz-plane, gaps at z=1,-1
	  do (setf (aref b-hats index) (make-3vector 0.0 (sin angle) (cos angle)))) ;yz-plane, gaps at z=1,-1
    (values a-hats b-hats)))

;;a, a', b, b' for the ACO loop.  N is the number of times b' goes between -z^ and z^
(defun make-functional-ACO (&optional (n 1))
  (values
   #'(lambda (x)
       (let ((cos (cos (* 2 pi x)))
	     (sin (sin (* 2 pi x))))
	 (values (make-3vector (/ cos 2 pi) (/ sin 2 pi) 0.0)
		 (make-3vector (- sin) cos 0.0))))
   #'(lambda (x)
       (setq x (mod (* x n) 1.0))
       (if (< x 0.5) (values (make-3vector 0.0 0.0 (/ x n)) ;0...1/2N
			     (make-3vector 0.0 0.0 1.0))    ;z^
	 (values (make-3vector 0.0 0.0 (/ (- 1.0 x) n))	    ;1/2N...0
		 (make-3vector 0.0 0.0 -1.0))))))	    ;-z^
