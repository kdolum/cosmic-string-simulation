(in-package "CL-USER")

;;;Explicit initial conditions of various types

;;Just make one diamond for testing
(defun make-test-diamond ()
  (let ((diamond (make-diamond :start (make-4vector 1.0 0.0 0.0 0.0)
			       :left (make-4vector 0.7 0.4 0.0 0.5)
			       :right (make-4vector 1.3 0.4 0.0 0.5))))
    (setf (diamond-end diamond) (compute-diamond-end diamond))
    diamond))

;;Make string from derivatives a'(t-sigma) and b'(sigma +t)
;;This should be rewritten to use the initial-conditions code like make-test-ab
;; Function rewritten to use the initial-conditions code like make-test-ab by ALE
(defun make-test-ab-prime (a-pfunction b-pfunction max-sigma points &key (offset zero-3vector) (string-number 0))
  (let ((starts (explicit-string-from-apbp a-pfunction b-pfunction max-sigma points :offset offset)))
    (if *test-string-explicit-output* (output-explicit-test-string starts)
      (process-explicit-string starts string-number))))

;;Make a string from functions for a and b.  MAX-SIGMA is the length of the string.
;;Functions will be called with POINTS values from 0 below MAX-SIGMA.  
;;The argument of a is tau-sigma, meaning that increasing argument goes around the string the opposite way.
(defun explicit-string-from-apbp (a-pfunction b-pfunction max-sigma points &key (offset zero-3vector))
  (mirror-image-let ((delta-sigma (/ max-sigma points))
		          (a-data (make-array points)))
      (mirror-images
       (loop for index below points
	      for sigma = (* delta-sigma index)
	       for a = (3to4vector (3vector-scale (funcall a-pfunction (+ sigma (/ delta-sigma 2))) (* 0.5 delta-sigma)) (* 0.5 delta-sigma))
	        do (setf (aref a-data index) a)))
      (explicit-string-from-discrete-apbp a-data b-data :offset offset)))

;;Accepts discrete functions ap(t-sigma) and bp(t+sigma).  These functions are lists of 3-vector values, and the index
;;in the list has no special meaning.  The implicit t+/sigma argument changes by the spatial distance at each step,
;;so the tangent vector automatically has magnitude 1.  The total argument change in the two functions must be the same.
;;We return a list of starting points for a diamond world sheet overlapping *initial-time*
(defun explicit-string-from-discrete-apbp (a b &key (offset zero-3vector))
  (loop with results = (make-array (length a))
	with start = (3to4vector offset *initial-time*)
	with start-new = start
	for index below (length a)
	do (if (< (4vector-t start) *initial-time*)  
	            (setf start-new (4vector+ start (aref b index)))
	          (setf start-new (4vector- start (aref a index))))
	do (setf (aref results index) start-new)
	do (setf start start-new)
	finally
	(return results)))

;;Create a string with a given shape, specified by position and
;;velocity functions
(defun make-test-string (x-function v-function max-sigma points &key open)
  ;;First make list of diamonds that touch at right and left
  (finish-test-string (make-test-string-1 x-function v-function max-sigma points open) open))

;;DIAMONDS is a list of diamonds linked corner to corner for the initial conditions
;;If OPEN, the string has :DELETED at the ends
;;But this does not work properly unless it is inert
(defun finish-test-string (diamonds open)
  (when *test-string-explicit-output*
    (when open (error "There are no open explicit strings"))
    (return-from finish-test-string
      (output-explicit-test-string diamonds)))
  (dolist (diamond diamonds)
    (handle-new-diamond diamond))
  ;;Now make future diamonds so they are connected
  (loop for last = (car (last diamonds)) then this
	for this in diamonds		;Every pair of adjacent diamonds
	for deleted = open then nil	;Open string starts with :deleted
	for new = (if deleted :deleted
		    (make-diamond :left (diamond-end last)
				  :right (diamond-end this)
				  :start (diamond-right last)
				  :sw last 
				  :se this
				  :tag (create-loop-tag (globalize-position (diamond-right last)))
				  :a-kink-created (diamond-a-kink-created last)
				  :b-kink-created (diamond-b-kink-created this)
				  ))
	unless (< (max (4vector-t (diamond-p this))
		       (4vector-t (diamond-q this)))
		  (/ (- diamond-span fudge-coordinates) 2))
	do (error "Diamond edge longer than half span: ~S" this)
	do (setf (diamond-ne last) new
		 (diamond-nw this) new)
	unless (eq new :deleted)
	do (setf (diamond-end new) (compute-diamond-end new))
	   (when open (setf (diamond-inertp new) t))
	   (handle-new-diamond new))
  (check-data-structures))
	  
;;This makes diamonds that are touching at the points, all of which are
;;at the present time.  We call the functions exactly POINTS times.
(defun make-test-string-1 (x-function v-function max-sigma points open)
  (loop with delta-sigma = (/ max-sigma points)
	with origin = (3to4vector (funcall x-function 0.0) *initial-time*)	;starting and ending point of loop, as 4 vector at time 0
	for index from 1 to (if open (1- points) points)			;Index of next point
	for left = origin then right	;go around segments of loop
	for right = (if (= index points) ;last point is just first
			origin
		      (3to4vector (funcall x-function (* delta-sigma index))  *initial-time*)) ;Compute next point, 4vector
	for center = (3vector-scale (3vector+ left right) 0.5)
	for delta-x = (3vector- right left) ;left-to-right distance
	for raw-v = (funcall v-function (* delta-sigma (- index 0.5))) ;midway between left and right
	;;Component perpendicular to delta-x
	for v = (3vector- raw-v (3vector-scale
				 delta-x
				 (/ (3vector-dot delta-x raw-v)
				    (3vector-dot delta-x delta-x))))
	for speed = (3vector-length v)
	for gamma = (/ 1 (sqrt (- 1 (* speed speed)))) ;boost factor
	;;Compute the real delta-sigma.  It may be different because of
	;;using a finite number of points
	for energy = (* gamma (3vector-length delta-x))
	for delta-xt = (3vector-scale v energy) ;bottom-to-top dist
	for delta-xt2 = (3vector-scale delta-xt 0.5) ;center to top dist
	for p0 = (/ energy 2)		;time component of vectors
	for start3 = (3vector- center delta-xt2)	;beginning of diamond
	for end3 = (3vector+ center delta-xt2)
	for start = (3to4vector start3 (- *initial-time* p0))
	unless (fudge= delta-sigma energy (/ delta-sigma points))
	do (warn "delta-sigma = ~S but energy = ~S.  Maybe your x is not properly parameterized"
		 delta-sigma energy)
	;;Convert everything to 4vectors and create diamond
	collect (let ((new (make-diamond :start start
					 :left left
					 :right right
					 :end (3to4vector end3 (+ *initial-time* p0))
					 :a-kink-created (globalize-position start)
					 :b-kink-created (globalize-position start)
					 :tag (create-loop-tag (globalize-position start))
					 )))
		  (when open (setf (diamond-inertp new) t))
		  new)))



(defun setup-circle (&optional (points 20))
  (make-test-string
   #'(lambda (sigma) (make-3vector (cos sigma) (sin sigma) 0.0))
   #'(lambda (sigma) (declare (ignore sigma)) zero-3vector)
   (* 2 pi)
   points))

(defun setup-links (&optional (points 16))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 1.5 (cos sigma)) (+ 2.0 (sin sigma)) 2.0)) ;position
   #'(lambda (sigma) (make-3vector (* 0.6 (sin (* 3 (+ 1 sigma)))) ;velocity
				   (* 0.6 (sin (* 4 (+ 2 sigma))))
				   (* 0.6 (sin (* 5 (+ 3 sigma))))))
   (* 2 pi)
   points)
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 2.5 (cos sigma))  2.0 (+ 2.0 (sin sigma)))) ;position
   #'(lambda (sigma) (make-3vector (* 0.6 (sin (* 3 (+ 1 sigma)))) ;velocity
				   (* 0.6 (sin (* 4 (+ 2 sigma))))
				   (* 0.6 (sin (* 5 (+ 3 sigma))))))
   (* 2 pi)
   points))

(defun setup-random-links (&key (points 80) (speed 0.1) velocities)
  (unless velocities
    (setq velocities (loop repeat points
			   collect (list (random speed) (random speed) (random speed)))))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 3.5 (cos sigma)) (+ 5.0 (sin sigma)) 5.0))
   #'(lambda (sigma) (declare (ignore sigma))
       (list-3vector (pop velocities)))
   (* 2 pi)
   points)
  (unless velocities
    (setq velocities (loop repeat points
			   collect (list (random speed) (random speed) (random speed)))))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 2.5 (cos sigma)) (+ 5.0 (/ (sin sigma) (sqrt 2))) (+ 5.0 (/ (sin sigma) (sqrt 2)))))
   #'(lambda (sigma) (declare (ignore sigma))
       (list-3vector (pop velocities)))
   (* 2 pi)
   points))
    

;;Circle with random perturbations
(defun setup-random-circle (&optional (points 8) velocities)
  (unless velocities
    (setq velocities (loop repeat points
			   collect (list (random 0.1) (random 0.1) (random 0.1)))))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 5.0 (cos sigma)) (+ 5.0 (sin sigma)) 5.0))
   #'(lambda (sigma) (declare (ignore sigma))
       (list-3vector (pop velocities)))
   (* 2 pi)
   points))

(defun setup-moving-random-circle (&optional (points 80) velocities)
  (unless velocities
    (setq velocities (loop repeat points
			   collect (list (random 0.1) (+ 0.5 (random 0.1)) (random 0.1)))))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 5.0 (cos sigma)) (+ 4.4 (sin sigma)) 5.0))
   #'(lambda (sigma) (declare (ignore sigma))
       (list-3vector (pop velocities)))
   (* 2 pi)
   points))

(defun setup-aco (&optional (points 52))
  (make-test-string
   #'(lambda (sigma) (make-3vector (/ (cos sigma) 2)
				   (/ (sin sigma) 2)
				   (/ (if (> sigma pi) (- (* 2 pi) sigma) sigma) 2)))
   #'(lambda (sigma) (make-3vector (/ (sin sigma) 2)
				   (/ (- (cos sigma)) 2)
				   (if (> sigma pi) -0.5 0.5)))
   (* 2 pi)
   points))

(defun setup-burden (&rest ignore)
  (declare (ignore ignore))
  (error "use SETUP-BURDEN-LOOP instead of SETUP-BURDEN"))

;;Sets up the Burden loop parametrized by m, n and psi
;;The created diamonds overlaps time *INITIAL-TIME* and the center is at spatial position OFFSET
;;Our functions return A and B, not A' and B'.  The tangent vectors are given by
;;A' = dA/dt = -dA/dsigma = (cos(m sigma/ampl), 0, sin(m sigma/ampl))
;;B' = dB/dt = dB/dsigma = (cos(n sigma/ampl), 0, -sin(m sigma/ampl)) rotated by psi
;;Cusps occur when A'=B'=(0,0,+/-1), so sigma_a = (j+1/2) pi ampl/m, sigma_b = (k+1/2) pi ampl/n
;;with j and k integers of opposite parity.  The A' that appears on diamond edges is the difference between two
;;successive values of A, and likewise for B.  We don't want this to be the cusp value.  Since the
;;functions are called for sigma = 2 pi ampl l/pts) with l an integer, we would not want a cusp to occur at
;;sigma = 2 pi ampl (l+1/2)/pts, which would be the difference of two successive l values.  So we need not to have
;;2 (l+1/2)/pts = (j+1/2)/m, i.e., pts/(2 m) = (l+1/2)/(j+1/2).  Thus pts/(2 m) and likewise pts/(2 n)
;;should not be the ratio of odd integers.
(defun setup-burden-loop (&optional (n 1) (m 2) (psi (/ pi 2.0)) (points 4096)
				    (amplitude 1.0) (offset zero-3vector))
  (flet ((check-ratio (m name)
	   (let ((ratio (/ points (* 2 m))))
	     (when (and (oddp (numerator ratio)) (oddp (denominator ratio)))
	       (warn "POINTS/2~A = ~S, so there will be diamonds moving at the speed of light" name ratio)))))
    (check-ratio m :m) (check-ratio n :n))
  (make-test-ab
   #'(lambda (sigma)
       (make-3vector (* amplitude (/ (sin (/ (* m (- sigma)) amplitude)) m))
		     0.0
		     (* amplitude (/ (cos (/ (* m sigma) amplitude)) m))))
   #'(lambda (sigma) 
       (make-3vector (* amplitude (/ (sin (/ (* n sigma) amplitude)) n) (cos psi))
		     (* amplitude (/ (sin (/ (* n sigma) amplitude)) n) (sin psi))
		     (* amplitude (/ (cos (/ (* n sigma) amplitude)) n))
		     ))
   (* amplitude 2 pi)
   points
   :offset offset))

;;The Garfinkel-Vachaspati loop, with 4 kinks
;;ALPHA is the angle between a' and b'
(defun setup-gv (alpha &key (offset zero-3vector))
  (let ((length 0.5))			;Diamonds must not be too big or we get errors
    (flet ((make (sigma alpha)
	     (cond ((or (zerop sigma) (fudge= sigma (* 2 length) 1e-15)) (make-zero-3vector))
		   ((fudge= sigma length 1e-15) (make-3vector (* length (cos alpha)) (* length (sin alpha)) 0.0))
		   (t (error "Unexpected sigma in setup-gv")))))
      (make-test-ab
       #'(lambda (sigma) (make sigma 0.0)) ;a always along x axis
       #'(lambda (sigma) (make sigma alpha))
       (* length 2)
       2				;Points
       :offset offset))))

;;Check gravity power against Garfinkel-Vachaspati paper
(defun gv-gravity-power (alpha)
  (let ((p (+ 1 (cos alpha)))
	(n (- 1 (cos alpha))))
    (* (/ 32 (expt (sin alpha) 2))
       (+ (* p (log (/ 2 p))) (* n (log (/ 2 n)))))))


(defun setup-loop-production (&optional (points 100))
  (make-test-string
   #'(lambda (sigma)
       (let* ((theta (- sigma (* 1.25 pi)))
	      (r (/ 1 (- (expt (* 1.30 pi) 2) (expt theta 2)))))
	 (make-3vector (* r (cos theta)) (* r (sin theta)) (* theta 0.02))))
   #'(lambda (sigma)
       (let* ((theta (- sigma (* 1.25 pi)))
;;	      (r (/ 1 (- (expt (* 1.30 pi) 2) (expt theta 2))))
	      )
	 (make-3vector (* 0.4 (cos theta)) (* 0.4 (sin theta)) (* theta -0.1))))
   (* 2.5 pi)
   points
   :open t))

;;Movie with wiggles passing to left and right
(defun setup-colliding-wiggles (&optional (points 300) (length 12.0) (amplitude 0.2))
  (flet ((xy (sigma omega)
	     (if (or (< sigma 0.5) (> sigma 1.5)) (values 0.0 0.0)
	       (let* ((angle (* pi (- sigma 0.5)))   ;0..pi in 0.5 to 1.5
		      (a (/ (* amplitude (sin angle)) omega))) ;Envelope
		 (values (* (if (= omega 1.0) (sqrt 2.0) 1.0) ;Expand in x because envelope small when this big
			    a (cos (* omega angle)))
			 (* a (sin (* omega angle))))))))
    (make-test-string
     #'(lambda (sigma)
	 (decf sigma (/ length 2))
	 (multiple-value-bind (x y) (xy (abs sigma) (if (plusp sigma) 7.0 1.0))
	   (make-3vector x y sigma)))
     #'(lambda (sigma)
	 (decf sigma (/ length 2))
	 (make-3vector 0.0 0.0 (- (signum sigma))) ;move at speed of light towards origin
	 )
     length
     points
     :open t)))
	 
;;The first 0 of J_0
(defparameter Bessel-zero 2.4048255577)

;;A loop without a cusp.  No Good.
(defun setup-no-cusp (&optional (points 32) (spread 0.4))
  (let ((main (sqrt (- 1.0 (expt spread 2))))) ;Amplitude of main part of derivative
    (flet ((ab (sigma)
	       (values (cos (* Bessel-zero (sin sigma))) ;Integrates to 0 regardless of main
		       (* Bessel-zero (sin sigma)))))
      (make-test-ab-prime
       #'(lambda (sigma) (multiple-value-bind (x y) (ab sigma)
			   (make-3vector (* main x)
					 (+ (* main y) (* spread (sin sigma)))
					 (* spread (cos sigma)))))
       #'(lambda (sigma) (multiple-value-bind (x y) (ab sigma) (make-3vector (- x) 0.0 y)))
       (* 2 pi)
       points)))) 

;;2d moving circle for periodic boundary conditions
(defun setup-moving-circle (&optional (points 8) (v 0.5))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 2.5 (cos sigma))
				   (+ 2.5 (sin sigma))
				   2.5))
   #'(lambda (sigma) (make-3vector ;;v 0.0
				   (/ v (sqrt 2)) (/ v (sqrt 2))
				   (* 0.01 (cos sigma))))	;small rotation
   (* 2 pi)
   points))

;;Slightly different version.
(defun setup-rising-circle (&optional (points 8) (v 0.5))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 2.5 (cos sigma))
				   (+ 4.0 (sin sigma))
				   2.6))
   #'(lambda (sigma) (make-3vector 0.0 (/ v (sqrt 2))
				   (* 0.01 (cos sigma))))	;small rotation
   (* 2 pi)
   points))

(defun setup-falling-circle (&optional (points 8) (v 0.5))
  (make-test-string
   #'(lambda (sigma) (make-3vector (+ 2.5 (cos sigma))
				   2.0	;X-Z plane
				   (+ 2.5 (sin sigma))))
   #'(lambda (sigma) (make-3vector 0.0 (- (/ v (sqrt 2)))
				   (* 0.01 (cos sigma))))	;small rotation
   (* 2 pi)
   points))


(defun setup-loop (points)
  (make-test-string
   #'(lambda (sigma) (make-3vector (cos sigma) (sin sigma) 0.0))
   #'(lambda (sigma) 
       (make-3vector 0.0 0.0 (* 0.5 (cos sigma))))
   (* 2 pi)
   points))



;;Get 32 bits of actually random number from /dev/urandom.  Sbcl (make-random-state t) just uses the current time,
;;which is the same in all processors
;;There is some hair involving not buffering the file, so we can only read the few bytes we need
;;We never return 0, which is used as a code for seed not supplied.
(defun get-really-random-number ()
  (loop for random = (with-open-stream (f (sb-sys:make-fd-stream (sb-unix:unix-open "/dev/urandom" sb-unix:o_rdonly 0)
								 :element-type '(unsigned-byte 32)
								 :buffering :none :file "/dev/urandom"))
		       (read-byte f))
	repeat 10
	unless (zerop random)
	do (format t "~&Random seed #x~X~%" random) (return random)
	finally (error "Failed to get a nonzero random seed")))  

;;Make a random-state from thirty-two random bits
(defun make-random-state-from-32 (random-32)
  (sb-kernel::%make-random-state (sb-kernel::init-random-state random-32)))

;;Save the random seed for reproducing runs
(defun save-random-seed (&optional (file (random-seed-file  *output-directory*)))
  (with-open-file (stream file :direction :output)
    (print *random-state* stream)))

;;Load previously stored random seed
(defun load-random-seed (&optional (file (random-seed-file (input-directory))))
  (with-open-file (stream file)
    (let ((seed (read stream)))		;get one form
      (assert (random-state-p seed))
      (format t "~&Reading random seed from ~S~%" file)
      (setq *random-state* seed)
      t)))

(defun maybe-delete-file (file)
  (when (probe-file file)
    (delete-file file)))

;;Remove previous output files.  If random-seed set, delete that as well.
(defun delete-old-output (random-seed)
  (format t "~&Deleting old output files...") (force-output)
  (mapc #'maybe-delete-file (directory (all-dump-files *job-number*))) ;Delete old dumps for this job
  (mapc #'maybe-delete-file (successor-files))
  (maybe-delete-file (job-loop-spectrum-file))
  (maybe-delete-file (job-bh-loop-spectrum-file))
  (when random-seed
    (maybe-delete-file (random-seed-file *output-directory*)))
  (format t "done~%"))

;;Gives the smallest conformal time for light to cross the simulation volume
;;Since wiggles travel roughly at half light velocity, this gives the time before 
;;correlations can contaminate evolution of the network.
(defun light-crossing-time (size)
  size)

(defparameter *x-max* 0.5)		;maximum size loop we wait for before stopping run

;;If we are logging loops we need to wait until loops produced with size x-max * total-size will
;;have been deleted.  This means the loop must oscillate twice and also the ratio of the loop size
;;to the horizon size (conformal time) must be less than *loop-preservation-threshold* (lpt).
(defun default-simulation-end (total-size initial-time lpt)
  (let ((tmax (+ initial-time (light-crossing-time total-size)))) ;Last conformal time simulated.
    ;;Last loop formed with conformal size xmax tmax
    ;;It must oscillate twice (with fixed physical period given currently by l_phys/2) to know that it is
    ;;not self-intersecting, and then one more time to reach the deletion point.
    (let ((extra (* 0.5 3 *x-max*)))
      (cond ((null lpt))		;No preservation threshold
	    ((zerop lpt)		;Never deleting loops
	     (error "Can't default duration when loops are not deleted"))
	    (t
	     (setq extra (max extra (- (/ *x-max* lpt) 1))))) ;Relative extra physical time for lpt t = tmax xmax
      (adjust-conformal-time-end tmax extra))))

;;Duration from the initial-time to the default-simulation-end time.
(defun default-duration-loops (total-size initial-time ldf *era*)
  (- (default-simulation-end total-size initial-time ldf) initial-time))

(define-timer :simulate)

;;General driver for running simulations
(defun simulate (body-function &rest args &key ((:time *report-time*)) &allow-other-keys)
  (maybe-time :overall
    (account-time :simulate
      (apply #'simulate-1 body-function args))))


(defmacro define-simulate-function (name &body body)
  `(defun ,name (body-function		;Function to call to create initial string
		 &key ,@simulate-keywords)
     ,@body))

;;The arguments to SIMULATE
(define-simulate-argument size)		;Total size of simulation volume, or NIL for infinite
(define-simulate-argument split-factor 1 fixnum) ;Number of job volumes in each dimension
(define-simulate-argument ijkl)		;Starting point in IJKL coordinates
(define-simulate-argument overwrite)  ;If set, clear output directory first.  Otherwise it must be empty or non-existent
(define-simulate-argument reproduce)	;T = reproduce previous run.  Directory: get seeds from that directory
(define-simulate-argument log)		;Write loop spectrum files
(define-simulate-argument timestamp 0.1 double-float) ;Timestamp at given interval
(define-simulate-argument start 1.0 double-float) ;Global time of start of whole simulation run.  Can be used in init.
(define-simulate-argument end nil double-float)	;Global ending time of whole simulation run.
(define-simulate-argument time)		;This argument is handled by SIMULATE instead of SIMULATE-1
(define-simulate-argument room nil)	;Report space utilization at the end
(define-simulate-argument usage-report nil) ;Report memory usage periodically
(define-simulate-argument energy-report nil) ;Report energy in initial and final string.

;;Dump control.
(define-simulate-argument dump-times nil) ;List of times at which to do full dumps (including lengths)
(define-simulate-argument length-times nil) ;List of times at which to dump lengths

;;Instead of or in addition to specifying a list, you can specify things this way:
(define-simulate-argument dump-start nil double-float) ;Global time of first dump.
(define-simulate-argument dump-interval nil double-float) ;Interval between dumps
(define-simulate-argument dump-end nil)			  ;No dumps after this global time (plus fudge factor)
(define-simulate-argument dump-diamonds t)

(define-simulate-argument length-start nil double-float) ;Same variables for writing length
(define-simulate-argument length-interval nil double-float)
(define-simulate-argument length-end nil)

(define-simulate-argument bh-number nil)
(define-simulate-argument bh-size nil)

;;bh intersection control
(define-simulate-argument bh-times nil)
(define-simulate-argument bh-start nil double-float)
(define-simulate-argument bh-interval nil double-float)
(define-simulate-argument bh-end nil)
(define-simulate-argument pointbh nil)
(define-simulate-argument all-bhs nil)
(define-simulate-argument bh-probability nil double-float)


;;Return ordered list of (TIME ONLY-LENGTH-P) for the global set of dumps.
;;OVERALL-START is the global time that the simulation began.  No dumps will be scheduled before then.
;;OUR-END is the global ending time of this simulation.  We will not return any times after that.
(defun dump-schedule (overall-start our-end dump-times dump-start dump-interval dump-end
				 length-times length-start length-interval length-end)
  (let* ((dumps (dump-schedule-1 overall-start our-end dump-times dump-start dump-interval dump-end))
	 (lengths (dump-schedule-1 overall-start our-end length-times length-start length-interval length-end))
	 (all (sort (copy-list		;SORT is destructive
		     (union dumps lengths :test #'fudge=global))
		    #'<)))
    (loop for time in all
	  collect (list time (not (member time dumps :test #'fudge=global))))))

;;Return times that are either in TIMES or in the set from START (or 0.0) by INTERVAL to END (if any)
;;and between OVERALL-START and OUR-END.  We are generous in which times are allowed, so, for example,
;; :start 4.4 :dump 0.1 will dump the initial conditions.
;;We return even times that come before the start of this job, so step number is just the position in the list.
;;Return value is not sorted.
(defun dump-schedule-1 (overall-start our-end times start interval end)
  (union				;Supplied times that are not beyond the end of the run
   (remove-if-not #'(lambda (time) (< time (+ our-end fudge-global-coordinates))) times) ;Too early checked elsewhere
   (and interval			;Generate times between start and end
	(loop for time from (or start 0.0) below (+ (if end (min end our-end) our-end) fudge-global-coordinates)
	      by interval
	      unless (< time (- overall-start fudge-global-coordinates)) ;ignore implicit too-early times
	      collect time))
   :test #'fudge=global))

(defun check-initial-time (start)
  (when (or (and (eq *era* :matter) (< start 4.0))
	    (and (eq *era* :radiation) (< start 2.0)))
    (error "Initial time is too small for ~A era" *era*)))

(defvar *simulate-dry-run* nil)		;If set, simulate does nothing but check arguments
(defvar *simulate-just-write-info* nil)	;If set, simulate does nothing but write run-info file

;;Check for problems in simulate arguments.  This happens even on a dry run.
(defun simulate-check-arguments (start end size log)
  (check-initial-time start)
  (unless (or size end)
    (error "You didn't give SIZE, so we are in infinite volume, but then you must give END"))
  (loop for last = nil then time
	for (time probability) in *intersection-probabilities*
	unless (and (numberp time) (numberp probability) (<= 0.0 probability 1.)
		    (or (null last) (> time last)))
	do (error "*INTERSECTION-PROBABILITIES* should be an ordered list of (time probability))"))
  (when *delete-unlucky-loops*
    (cond ((and (= *intersection-probability* 1.0) (null *intersection-probabilities*))
	   (warn "Setting :DELETE-UNLUCKY-LOOPS with p=1 does not make sense"))
	  (log
	   (warn "You're deleting unlucky loops in self-intersecting trajectories and also logging loop sizes, so the results will not be accurate"))))
  (when (and (not (zerop *loop-preservation-dump-start*)) ;Not the default
	     (not *loop-preservation-dump-x*))
    (error "Giving :LOOP-PRESERVATION-DUMP-START without :LOOP-PRESERVATION-DUMP-X does not make sense")))

(defun simulate-check-dump-arguments (overall-start times start interval end dump-p)
  (loop for time in times
	when (< time (- overall-start fudge-global-coordinates))
	do (error "Dump explicitly requested for time ~F before overall start time ~F" time overall-start))
  (unless (or times interval)		;Not dumping.  Shouldn't specify dump parameters
    (when (and dump-p			;Dump, not just length
	       *loop-preservation-dump-x*)
      (error "Giving :LOOP-PRESERVATION-DUMP-X without doing any dumps does not make sense"))
    (when start
      (error "Giving start time but no interval does not make sense"))
    (when end
      (error "Giving end time but no interval does not make sense"))))

(define-simulate-function simulate-1
  (declare (ignore time))		;Handled already in SIMULATE
  (simulate-coerce)			;Coerce argument types
  (simulate-check-arguments start end size log)
  (simulate-check-dump-arguments start dump-times dump-start dump-interval dump-end t)
  (simulate-check-dump-arguments start length-times length-start length-interval length-end nil)
  (when *simulate-dry-run* (return-from simulate-1 "Dry run succeeded"))
  (when pointbh
    (setf *pointbh* pointbh)
    (setf *bh-probability* bh-probability))
  (when all-bhs
    (setf *all-bhs* all-bhs))
  (when (and bh-number (null *pointbh*))
    (setf bh-times (bh-schedule bh-start bh-interval bh-end)))
  (setf *dump-diamonds* dump-diamonds)
  (if *dump-diamonds*
      (format t "Dumping Diamonds~%")
    (format t "NOT Dumping Diamonds~%"))
  (when *simulate-just-write-info*
    ;;Use supplied end or last time covered by overall run.  This is before the actual ending time of the last layer.
    ;;We only come here in the manager, so this can be figured out.
    (let* ((overall-end (or end (+ start (jobs-duration size split-factor *manager-jobs*))))
	   (dump-schedule (dump-schedule start overall-end ;List of all dumps and whether they involve only lengths
					 dump-times dump-start dump-interval dump-end
					 length-times length-start length-interval length-end))
	   (length-times (mapcar #'car dump-schedule)) ;Times when length will be dumped, i.e., any dump
	   (full-dump-times (mapcar #'car (remove-if #'second dump-schedule)))) ;Times of actual dumps
      (return-from simulate-1
	(write-run-info-file :start start
			     :end overall-end
			     :dump-times full-dump-times :length-times length-times :bh-times bh-times
			     :split-factor split-factor :total-size size
			     :bh-size bh-size :bh-number bh-number))))
  (when *job-number*
    (format t "~&********** JOB ~D **********~%" *job-number*))
  (cond (size				;Size given: normal case.  It is the periodicity distance.
	 (initialize :total-size size :split-factor split-factor :ijkl ijkl :job-number *job-number*
		     :ijkl-origin (3to4vector zero-3vector *time-offset*)
		     :start start
		     :bh-size bh-size :bh-number bh-number :bh-start bh-start)) 
	(t				;Infinite volume
	 (initialize :total-size nil)))
;  (if *bh-number* ;if number of BH is set create the blackholes.dat file
;      (if (eq *job-number* 0)
;	  (progn
;	    (sort-bhs)
;	    (format t "BHs has been created"))
;	(format t "BHs created by other job")))
  (let* ((end-local (and end (local-time end))) ;Convert to local
	 ;;Local ending time of job.  Fudge factor here is to make this go after the successor output event.
	 (job-end (+ (job-end-t) fudge-coordinates)))
    (when (> *initial-time* job-end)
      (error "Initial time ~S is after end of job volume ~S" *initial-time* job-end))
    (unless (and end-local (<= end-local job-end)) ;default or restrict end to end of cube
      (setq end-local job-end))
    (let* ((our-end (global-time end-local)) ;Global ending time of this job
	   (dump-schedule (dump-schedule start our-end ;Dumps to be done by this job
				      dump-times dump-start dump-interval dump-end
				      length-times length-start length-interval length-end))
	   (*last-dump-time* (car (last dump-times))) ;Last explicitly given dump time for preservation code, or NIL
	   (*overall-end-time* (if *last-dump-time* ;take the largest result possible
				   (if *total-size* (max *last-dump-time* (+ start *total-size*)) *last-dump-time*)
				 job-end))) ;fall back on the guaranteed number if nothing's set
      (format t "Overall end time for minimum diamond width calculation is ~S~%" *overall-end-time*)
      (when overwrite									  ;Overwriting old files?
	(delete-old-output (not reproduce)))						  ;Delete them.
      (ensure-directories-exist (job-filename *output-directory* *job-number* "")) ;Make output directory if needed
      (cond (*random-seed*							   ;Do things differently under manager
	     (setq *random-state* (make-random-state-from-32 *random-seed*)))
	    (reproduce								    ;Doing run again
	     (let ((*input-directory* (if (eq reproduce t) *input-directory*	    ;Reproduce from here
					reproduce)))				    ;or given place
	       (load-random-seed)))						    ;get random state from before
	    (t (setq *random-state* (make-random-state-from-32 (get-really-random-number)))
	       (save-random-seed)))	;Save it out for later
      (funcall body-function)           ;Now initialize the run   
      (when (and (null *pointbh*) *bh-number*) ;if number of BH is set create the blackholes.dat file
	(format t "File: ~S~%" (probe-file "blackholes.dat"))
	(format t "BHs created by other job")
	(my-bhs)
	(clean-mybhs))
      ;;Regardless of how long we will run now, put in the event to communicate strings to successors
      (main-calendar-add (job-end-t) (make-successor-output (job-end-t)))
      ;;Put in all events for dumping
      (loop with created = nil		;Create directory first time it is needed
	    with step = 0		;Count actual dumps
	    for (time only-lengths-p) in dump-schedule
	    do (when (>= (local-time time) (- *initial-time* fudge-coordinates)) ;Is this dump in our time period?
		 (request-dump (unless only-lengths-p step) time)
		 (unless (or only-lengths-p created)
		   (ensure-directories-exist (dump-file *job-number* 0)) ;If dumping create dump directory once
		   (setq created t)))
	    do (unless only-lengths-p (incf step))) ;Count step, even if dump is too early for us
      (when timestamp
	(main-calendar-add *initial-time* (make-timestamp *initial-time* timestamp)))
      (when usage-report
	(let ((first-usage-report (+ *initial-time* usage-report)))
	  (main-calendar-add first-usage-report (make-usage-report first-usage-report usage-report))))
      (let ((initial-energy (and energy-report (total-length (get-paths (+ *initial-time* .01))))))
	(setq *start-real-time* (get-internal-real-time) ;Save start times for usage reporting
	      *start-run-time* (get-internal-run-time)
	      *last-real-time* 0
	      *last-run-time* 0)
	(when log			;Writing loop file
	  (setq *loops-found* (make-array 1000 :element-type 'double-float :fill-pointer 0 :adjustable t))
	  (when *bh-start*
	    (setq *bh-loops-found* (make-array 1000 :element-type 'double-float :fill-pointer 0 :adjustable t))))
	(when (and *pointbh* *bh-start*)
	  (when (and (>= (local-time *bh-start*) (current-time)) (<= (local-time *bh-start*) end-local))
	    (evolve-until (local-time *bh-start*))
	    (format t "Creating point like BHs at ~S~%" *bh-start*)
	    (create-point-bhs)))
	(when (and (null *pointbh*) bh-times) ;if bh-time is set check that is during the evolution and analyse the intersections
	    (loop for tbh in bh-times
		  do (when (and (>= (local-time tbh) (current-time)) (<= (local-time tbh) end-local))
		       (evolve-until (local-time tbh)) 
		       (format t "Checking BH intersections at ~S~%" tbh)
		       (create-bhs))))
	(evolve-until end-local)	;Do it
	(when (and *pointbh* *bh-start*)
          (when (and (>= (local-time *bh-start*) (local-time (global-job-start *total-size* *time-offset* *job-number*))) (<= (local-time *bh-start*) end-local))
	    (writepointbhs)))
	(when log
	  (if *loops-output*		    ;Stream already open
	      (write-recorded-loops)	    ;Use it
	    (with-open-file (*loops-output* (job-loop-spectrum-file) :direction :output
					    :if-does-not-exist :create :element-type '(unsigned-byte 64))
	      (write-recorded-loops))) ;Open new file and use it
	  (when *bh-start*
	    (if *bh-loops-output*
		(write-recorded-bh-loops)
	      (with-open-file (*bh-loops-output* (job-bh-loop-spectrum-file) :direction :output
						 :if-does-not-exist :create :element-type '(unsigned-byte 64))
	      (write-recorded-bh-loops)))))
	(when *length-output*
	  (force-output *length-output*)) ;Avoid loss of buffered lengths data
	(when *report-time*
	  (format t "~&~:D diamonds advanced, ~D non-null sides could not be fixed immediately~%"
		  *advance-diamond-count* *compute-diamond-adjust-count*)
	  (format t "~D intersections performed ~@[(~D rejoinings)~], ~D rejoining~:P suppressed~@[ (~$%)~] to avoid monsters~:[~;, ~D skipped because of probability~]~%"
		  *intersections-performed* (and *count-rejoinings* *rejoinings-performed*) *suppressed-rejoinings*
		  (and (plusp *intersections-performed*)
		       (round (/ (* *suppressed-rejoinings* 100) *intersections-performed*)))	
		  (plusp *intersections-unlucky*)
		  *intersections-unlucky*))
	(format t "~&Finished with global time ~$~%" (global-time (current-time)))
	(when energy-report
	  (format t "~% Total energy in initial conditions is ~D." initial-energy)
	  (format t "~% Total energy in non-inert remains is ~D." (total-length))
	  (format t "~% Energy difference is ~D." (- initial-energy (total-length)))
	  )))
	(when room
	  (gc :full t)
	  (room))
	(format t "~&********** JOB ~D COMPLETED **********~%~%" *job-number*)))

(defun simulate-continue (job-number &rest simulate-keys)
  (apply #'simulate
	 #'(lambda ()
	     (read-initial-strings))
	   :job-number job-number
	   simulate-keys))

;;Run simulation and continue through job number N-1, all in this processor
(defun simulate-n (n function &rest simulate-keys)
  (apply #'simulate function :job-number 0 simulate-keys)
  (loop for job from 1 below n
	do (format t "~&****************************** JOB ~D ******************************~%" job)
	do (apply #'simulate-continue job simulate-keys)))
  
(defun simple-test (&rest simulate-keys &key (overwrite t) (size 20.0)
			  (start-time 5.0)
			  &allow-other-keys)
  (apply #'simulate #'setup-moving-random-circle
					;This gives the defaults or harmlessly duplicates the keywords
	 :overwrite overwrite :size size :start start-time
	 simulate-keys))

(defun test-aco (&rest simulate-keys &key (overwrite t) (size 20.0)
			  (start-time 1.0)
			  &allow-other-keys)
  (apply #'simulate
	 #'setup-aco
	 :overwrite overwrite :size size :start start-time
	 simulate-keys))

(defun test-collision (&rest simulate-keys &key (overwrite t) (size 5.0)
			  &allow-other-keys)
  (apply #'simulate
	 #'(lambda ()
	     (setup-rising-circle)
	     (setup-falling-circle))
	 :overwrite overwrite :size size
	 simulate-keys))

;; simulates links with varying resolution initial conditions (no ramdom velocities)
(defun test-links (&rest simulate-keys &key (overwrite t) (points 16) (size 4.0)
			    &allow-other-keys)
  (apply #'simulate
	 #'(lambda ()
	     (setup-links points))
	 :overwrite overwrite :size size
	 simulate-keys))

;;Vachaspati-Vilenkin initial conditions
(defun do-vv (&rest keys)
  (apply #'do-initial-1 #'setup-vv keys))

;;Explicit initial conditions.  Default start time is 0.0.
(defun do-explicit (&rest keys &key (start 0.0) &allow-other-keys)
  (apply #'do-initial-1 #'setup-explicit
	 :start start keys))

;;Driver for simulations with initial condition surface.
;;START here is what you want the clock to read at the initial condition time.
(defun do-initial-1 (setup-function &rest simulate-keys &key size (split-factor 1) (start 1.0) &allow-other-keys)
  (let ((start-offset (/ (* (/ 3.5 4.0) size (sqrt (/ 2.0 3.0))) split-factor))) ;Desired time from origin to initial
    (apply #'simulate #'(lambda () (maybe-setup-initial setup-function))
	   :time-offset (- start start-offset)
	   :start start
	   simulate-keys)))

;;Minimum size of simulation
(defun minimum-vv-simulation-size (split-factor)
  (if (> split-factor 1)
      (* split-factor 17.0)		;Size of each must be at least 12 sqrt(2).  See notes.
    20.0)				;This is about the smallest you can have without ending up with lattice overlap
  )

(defun test-vv (&rest simulate-keys &key (overwrite t)
		      size
		      (split-factor 1)
		      &allow-other-keys)
  (unless size
    (setq size (minimum-vv-simulation-size split-factor)))
  (apply #'do-vv :overwrite overwrite :size size simulate-keys))

;;Setup initial conditions if we are in one of the first layers
(defun maybe-setup-initial (setup-function)
  (cond ((initialization-surface-p)
	 (maybe-time :setup (funcall setup-function)))
	(t (format t "~&Initialization surface in the past.  Reading strings...")
	   (read-initial-strings))))


;;;Checking system

;;Check various things in our data structures
(defun check-data-structures ()
  (map-main-calendar #'check-main-calendar-entry)
  (map-final-diamonds #'check-diamond)
  (map-read-diamonds #'check-diamond)
;;Slow
;;  (check-angles)
  t)

;;Check that cells system is working.
(defun check-intersection-list (diamond)
  (let* ((seen (make-hash-table :test #'eq)) ;Creates a hash table to keep track of diamonds from end points
	 (possibilities nil))
    ;;Get possibilities old way.
    (map-main-calendar
     #'(lambda (time object)		;Find any diamond on calendar
	 (declare (ignore time))
	 (when (typep object 'diamond)	                                       ;Look at all diamonds
	   (unless (diamond-inertp object) ;Doesn't go in list if inert
	     (unless (gethash object seen)
	       (setf (gethash object seen) t) ;Set the hash table to know that we have checked this diamond already
	       (push object possibilities))))))
    ;;Get possibilities that have bounding box overlap
    (flet ((bounding-box-check (d)
	     (and (not (eq d diamond))
		  (diamond-box-overlap d diamond))))
      (setq possibilities (remove-if-not #'bounding-box-check possibilities))
      (let ((new-possibilities nil))
	(do-diamond-cells-diamonds (d diamond)
	  (when (bounding-box-check d) (push d new-possibilities)))
	;;Make sure there are no possibilities that the old system found and the new system lost.
	;;Vice versa is OK because they might belong to neighbors and not be in the calendar
	(unless (null (set-difference possibilities new-possibilities))
	  (error "Cell system failed for ~D" diamond))))))

(defun check-angles ()
  (map-main-calendar
   #'(lambda (time diamond1)
       (declare (ignore time))
       (when (diamondp diamond1)
	 (let ((p (3vector- (diamond-left diamond1) (diamond-start diamond1))))
	   (map-main-calendar
	    #'(lambda (time diamond2)
		(declare (ignore time))
		(when (diamondp diamond2)
		  (let* ((q (3vector- (diamond-right diamond2) (diamond-start diamond2)))
			 (cos (/ (3vector-dot p q) (3vector-length p) (3vector-length q))))
		    (when (> cos (- 1.0 1e-14))
		      (error "SW edge of ~S, ~S is too close to SE edge of ~S, ~S, cos ~F"
			     diamond1 p diamond2 q cos)))))))))))
      
  

(defun check-main-calendar-entry (time object)
  (assert (>= time (current-time)))
  (etypecase object
    (diamond (assert (= (4vector-t (diamond-end object)) time))
	     (check-diamond object t))
    (intersection (assert (= (4vector-t (event-location object)) time))
		  (assert (member object (diamond-pending-intersections (intersection-diamond-1 object))))
		  (when (intersection-diamond-2 object)
		    (assert (member object (diamond-pending-intersections (intersection-diamond-2 object))))))
    (dump-request
     (assert (fudge= (global-time time) (dump-request-time object) fudge-global-coordinates))) ;Check that it is in at the right time
;;    (age-objects)
    (timestamp)
    (stop-evolving)
    (successor-output)
    ))


(defun check-object-neighbors (object neighbors)
  (assert neighbors)			;shouldn't be in table if none
  (typecase object
    (diamond (check-diamond object nil))))

(defun check-position-standardized (position)
  (assert (< (4vector-Euclidian-distance position (standardize-position position)) fudge-coordinates)))

(defun check-diamond (diamond &optional from-queue-p)
  (assert (diamondp diamond))
;;  (check-position-standardized (diamond-left diamond))
;;  (check-position-standardized (diamond-start diamond))
;;  (check-position-standardized (diamond-right diamond))
  (assert (< (3vector-length (diamond-p diamond)) (/ diamond-span 2)))
  (assert (< (3vector-length (diamond-q diamond)) (/ diamond-span 2)))
  ;;Sides of diamonds should be null, but some deviation should be tolerated
  ;;when diamonds are cut.
;;Not doing this currently.  Definition of 4vector-dot was changed
;;  (assert (< (4vector-dot (diamond-p diamond) (diamond-p diamond)) 1d-5))
;;  (assert (< (4vector-dot (diamond-q diamond) (diamond-q diamond)) 1d-5))
;  (assert (> (3vector-length (diamond-p diamond)) (/ minimum-diamond-width 10)))	;allow something for redshifts
;  (assert (> (3vector-length (diamond-q diamond)) (/ minimum-diamond-width 10)))
  (mirror-images (assert (diamond-a-kink-created diamond)))
  (unless (diamond-inertp diamond)
    (assert (diamond-minima diamond))	;Make sure that bounding box has been installed
    (assert (diamond-maxima diamond)))
  (mirror-image-let ((ne (diamond-ne diamond))
		     (se (diamond-se diamond)))
    (mirror-images
     (assert (not (and ne se))) ;At most one right (left) link 
     (let ((east (or ne se)))
       (when from-queue-p ;Could be an old diamond whose advancement message has been delayed
	 (assert (or east (right-rejoining-junction diamond))))
       (when east		      ;If exists, make sure right type
	 (assert (typep east '(or diamond handle keyword))) ;east may be :deleted, a keyword
	 (loop for index below (resource-pointer diamond-resource)
	       when (eq (aref (resource-objects diamond-resource) index) east)
	       do (error "~A link goes to a diamond in the resource for reuse" :east))))
     (when (diamondp ne)
       (assert (eq (diamond-end diamond)   ;End of this diamond
		   (diamond-left ne)))	   ;the leftmost point of this one
       ;;This might not hold if the NE diamond has been cut.
       ;;(assert (eq (diamond-right diamond) (diamond-start ne)))
       (if (not (eq (diamond-sw ne) diamond))
	   (format t "~S~%" diamond))
       (assert (eq (diamond-sw ne) diamond)) ;reverse link
     )))
  (assert (not (4vector= (diamond-left diamond) (diamond-right diamond)
			 1e-14)))
  (when from-queue-p			;If in our queue, more tests
    (assert (or (diamond-finalp diamond)
		(point-mine (diamond-end diamond)))) ;Shouldn't be here unless we own it
    ;;If the end is in the past, we should have advanced it already
    (assert (<= (current-time) (4vector-t (diamond-end diamond)))
	    ()
	    "~S ends in the past" diamond)
    )
  (check-diamond-intersections diamond)
  (check-points-eq diamond)
;;Very slow
;;  (check-intersection-list diamond)
  )

;;Not used
(defun check-diamond-neighbors (diamond)
  (mirror-images
   (assert (or (diamond-w diamond)
	       (left-rejoining-junction diamond)))
   )			;mirror-images
  )


(defun check-points-eq (diamond)
  (mirror-images
   (let ((nw (diamond-nw diamond))
	 (sw (diamond-sw diamond)))
     (if (and (diamondp nw)
	      (4vector= (diamond-right nw) (diamond-end diamond) fudge-coordinates))
	 (assert (eq (diamond-right nw) (diamond-end diamond))))
     (if (and (diamondp nw)
	      (4vector= (diamond-start nw) (diamond-left diamond) fudge-coordinates))
	 (assert (eq (diamond-start nw) (diamond-left diamond))))
     (if (and (diamondp sw)
	      (4vector= (diamond-right sw) (diamond-start diamond) fudge-coordinates))
	 (assert (eq (diamond-right sw) (diamond-start diamond)))) ;failing
     (if (and (diamondp sw)
	      (4vector= (diamond-end sw) (diamond-left diamond) fudge-coordinates))
	 (assert (eq (diamond-end sw) (diamond-left diamond)))))))


(defun check-diamond-intersections (diamond)
  (let*((pending-intersection-list (diamond-pending-intersections diamond)))
    (loop for cons on pending-intersection-list
	  for intersection = (car cons)
	  do (check-pending-in-both intersection)
	  do (loop for other in (cdr cons) ;check for duplicates
		   when (4vector= (intersection-spacetime intersection)
				  (intersection-spacetime other)
				  fudge-coordinates)
		   do (error "Two intersections in the same place?")))))

(defun check-dump-junction (junction)
  (when (and (rejoining-junction-p junction)
	     (rejoining-junction-dump-time junction))
    (let* ((diamond (or (junction-left-diamond junction) (junction-right-diamond junction)))
	   (junction-time (global-time (diamond-position-time diamond :a (rejoining-junction-a junction)
							      :b (rejoining-junction-b junction)))))
      (unless (fudge= (rejoining-junction-dump-time junction) junction-time
		      fudge-global-coordinates)
	(error "Junction for dump at time ~F has time ~F" (rejoining-junction-dump-time junction) junction-time))
      ;;	(format t "~&Rescaled junction OK:~%")
      ;;	(describe junction)
      )))


(defun check-pending-in-both (intersection)
  (let* ((d1 (intersection-diamond-1 intersection))
	 (d2 (intersection-diamond-2 intersection))
	 (list1 (diamond-pending-intersections d1))
	 (list2 (diamond-pending-intersections d2))
	 (number (length (intersection list1 list2))))
    (cond ((= number 1)
	   nil)
	  ((= number 0)
	   (error "pending intersection lists not compatible: ~A" intersection))
	  (t
	   (error "two diamonds intersect more than once: ~A" intersection)))))


;;Plotting and animation.

;;Plots the past edge of the diamonds that we have, without regard to time
(defun plot-string-past (diamonds &rest gnuplot-keys &key (styles :lines) &allow-other-keys)
  (let ((data (coerce diamonds 'vector)))
    (apply #'gnuplot 1 (length data)
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (unless (eq point :title)
		 (values-list (3vector-list (diamond-start (aref data point))))))
	   :styles styles
	   gnuplot-keys)))

(defun describe-events ()
  (let ((table (make-hash-table :test #'eq)))
    (map-main-calendar
       #'(lambda (time thing)
	   (declare (ignore time))
	   (incf (gethash (type-of thing) table 0))))
    (maphash
     #'(lambda (type count)
	 (format t "~&~D event~:P of type ~S~%" count type))
     table)))


(defun do-animate (start end
			  &rest gnuplot-keys &key (delta 0.1) (sleep 0.5) filespec &allow-other-keys)
  (unless filespec (plot-time-slice :time start)) ;this prevents :reuse t from using a bad prior plot
  (loop with minimum = (job-coordinate-minimum)
	with maximum = (job-coordinate-maximum)
	with args = (append `(:xrange ,(cons (- minimum diamond-span)
					     (+ maximum diamond-span))
			      :yrange ,(cons (- minimum diamond-span)
					     (+ maximum diamond-span))
			      :zrange ,(cons (- minimum diamond-span)
					     (+ maximum diamond-span)))
				gnuplot-keys)
	for time from start by delta to end
	for step from 0
	do (evolve-until time)
	do (format t "~&~4$" time)
	do (apply #'plot-time-slice :time time :reuse (not filespec)
		  :png (and filespec (format nil filespec step))
		  :title (format nil "Time ~$" time)
		  args)
	do (unless filespec (sleep sleep))))


;;Plot CPU percentages from :usage-report
(defun plot-cpu-usage (directory &rest gnuplot-keys)
  (let ((data
	 (loop for rank from 0
	       for file = (format nil "~A/rank-~D/output" directory rank)
	       while (probe-file file)
	       collect (with-open-file (stream file)
			 (loop for line = (read-line stream nil)
			       with percentages = nil
			       while line
			       when (and (> (length line) 10)
					 (string-equal line "Time " :end1 5))
			       do (let* ((end (position #\: line))
					 (time (read-from-string line nil nil :start 5 :end end))
					 (position (position #\( line)))
				    (unless position (error "Can't find CPU percentage in ~S" line))
				    (let ((cpu (parse-integer line :start (1+ position) :junk-allowed t)))
				      (push (list time cpu) percentages)))
			       finally (return (nreverse percentages)))))))
    (format t "~&~D processors" (length data))
    (when data
      (format t ", ~D steps" (length (car data)))
      (if (car data)
	  (apply #'gnuplot
		 (length data)
		 (length (car data))	;Should all be the same
		 #'(lambda (plot point)
		     (unless (eq point :title)
		       (values-list (nth point (nth plot data)))))
		 :key nil
		 gnuplot-keys)
	(format t "~&No data found in files")
	))))

;;Plot speed data from usage report
(defun plot-speed (directory  &rest gnuplot-keys &key yrange (styles :lines) real diamond average &allow-other-keys)
  (let ((data
	 (loop for file in (directory (format nil "~A/rank-*/output" directory))
	       collect (with-open-file (stream file)
			 (loop for line = (read-line stream nil)
			       with data = nil
			       while line
			       when (and (> (length line) 10)
					 (string-equal line "Time " :end1 5))
			       do (let* ((end (position #\: line))
					 (time (read-from-string line nil nil :start 5 :end end))
					 (tag (cond (diamond "advanced, at ") (real "real speed ") (t "speed ")))
					 (position (search tag line)))
				    (unless position (error "Can't find speed in ~S" line))
				    (let ((speed (parse-integer line :start (+ position (length tag)) :junk-allowed t)))
				      (push (list time speed) data)))
			       finally (return (nreverse data)))))))
    (format t "~&~D processors" (length data))
    (when data
      (format t ", ~D steps" (length (car data)))
      (cond ((null (car data))
	     (format t "~&No data found in files"))
	    (average
	     (apply #'gnuplot 1 (length (car data))	;Should all be the same
		    #'(lambda (plot point)
			(declare (ignore plot))
			(unless (eq point :title)
			  (values (first (nth point (car data)))	;time
				  (/ (loop for processor in data sum (second (nth point processor)))
				     (length data)))))
		    :styles styles
		    gnuplot-keys))
	    (t
	     (let ((yrange
		    (or yrange
			(cons 0
			      (loop for processor in data
				    maximize (loop for (time speed) in processor
						   when (> time 2.0) ;Before this there are occasional times when there's no work, so the speed is very high
						   maximize speed))))))
	       (apply #'gnuplot
		      (length data)
		      (length (car data)) ;Should all be the same
		      #'(lambda (plot point)
			  (unless (eq point :title)
			    (values-list (nth point (nth plot data)))))
		      :key nil
		      :yrange yrange
		      :styles styles
		      gnuplot-keys)))
	))))

(defmacro monitor-allocation (&body body)
  `(monitor-allocation-1 #'(lambda () ,@body)))

(defun monitor-allocation-1 (function &key (threshold 10000) (gc t))
  (when gc (gc :full t))
  (let ((start (sb-vm::type-breakdown :dynamic)))
    (funcall function)
    (when gc (gc :full t))
    (loop for (size count type) in (sb-vm::type-breakdown :dynamic)
	  while (> size 1000)
	  for (old-size old-count) = (find type start :key #'third)
	  for bytes = (- size old-size)
	  for objects = (- count old-count)
	  sum bytes into total-bytes
	  when (> bytes threshold)
	  collect (list bytes objects type) into result
	  finally
	  (setq result (sort result #'> :key #'first))
	  (loop for (bytes objects type) in result
		do (format t "~&~:D more bytes and ~D more objects of type ~A" bytes objects type))
	  (format t "~&Total dynamic space increase ~:D bytes" total-bytes))))
	     
;;Finds number of cons cells in the values of a hash table
(defun hash-table-list-usage (table)
  (let ((conses 0))
    (maphash #'(lambda (key value)
		 (declare (ignore key))
		 (when (consp value)
		   (incf conses (length value))))
	     table)
    conses))		 

;;describe memory usage of hash tables; return byte count
(defun hash-table-memory-usage (&rest tables)
  (loop for table in tables
	for size = (hash-table-size table)
	for conses = (hash-table-list-usage table)
	sum size into total-size
	sum (hash-table-count table) into total-count
	sum (+ (* size 4 8) ;There are 2 arrays of length SIZE, and one of length 2 SIZE, each with 64-bit elements


	       (* conses 16)) into total-bytes
	finally 
	(format t "~D objects, space for ~D, ~[~:;~:*~D conses in values, ~]~:D bytes" total-count total-size conses total-bytes)
	(return total-bytes)))

;;Scan dynamic space for structures
(defun dynamic-space-structure-memory-usage ()
  (let ((data (loop for type in '(diamond handle (simple-array double-float (4)))
		    collect (list type 0 0))))
    (sb-vm::map-allocated-objects
     #'(lambda (object object-type size)
	 (declare (ignore object-type))
	 (loop for cons in data
	       when (typep object (car cons))
	       do (incf (second cons))
	          (incf (third cons) size)))
     :dynamic)
    (prog1 (loop for (type count size) in data
		 do (format t "~&~D ~Ss use ~:D bytes" count type size)
		 sum size)
      (force-output))))


;;Find our objects through our data structures
(defun structure-memory-usage ()
  (let ((diamonds 0)
	(handles (hash-table-count *object-handle*))
	(4vectors 0)
	(total-bytes 0))
    (labels ((do-point (point)
	       (declare (ignore point))
	       (incf 4vectors))
	     (do-diamond (diamond)
	       (incf diamonds)
	       (do-point (diamond-start diamond))		;Every diamond comes with a starting point
	       (unless (diamondp (diamond-ne diamond)) ;No NE diamond?
		 (do-point (diamond-right diamond)))		;Then right of this is not anyone's start
	       (unless (diamondp (diamond-w diamond)) ;No NW or SW diamond
		 (do-point (diamond-left diamond)))		;Our left is not anyone's start or right
	       (unless (or (diamondp (diamond-ne diamond)) ;No diamonds to future
			   (diamondp (diamond-nw diamond)))
		 (do-point (diamond-end diamond)))		;Count our end
	       ))
      (map-main-calendar
       #'(lambda (time object)		;Find diamonds in calendar
	   (declare (ignore time))
	   (when (diamondp object)
	     (do-diamond object))))
      (map-final-diamonds #'do-diamond)
      (map-read-diamonds #'do-diamond)
      )
    (macrolet ((do-type (type &optional size)
		 (let ((counter (intern (format nil "~AS" type) (symbol-package type)))
		       (size (or size `(defstruct-type-bytes ',type))))
		   `(progn (format t "~&~D ~Ss use ~:D bytes" ,counter ',type
				   (* ,size ,counter))
			   (incf total-bytes (* ,size ,counter))))))
      (do-type diamond)
      (do-type handle)
      (do-type 4vector (* 6 8))		;Array of 4 double-floats
      )
    total-bytes))
    

;;Number of bytes used for objects of type TYPE
(defun defstruct-type-bytes (type)
  ;;Length in layout includes slots and layout object itself.
  ;;Most add one for header, then roundup to multiple of 2 words
  (* 16 (ceiling (1+ 
		  (#-sbcl-wrappers sb-kernel:layout-length #+sbcl-wrappers sb-kernel:wrapper-length
			   (#-sbcl-wrappers sb-kernel:classoid-layout #+sbcl-wrappers sb-kernel:classoid-wrapper
				    (sb-kernel:find-classoid type))))
		 2)))




;;Describe calendar memory usage and return total
(defun calendar-memory-usage (&rest calendars)
  (loop for calendar in calendars
	sum (calendar-days calendar) into days
	sum (calendar-event-count calendar) into events
	finally (let* ((conses (* events 2)) ;get each event has (time . event) and a cons to put this in a list
		       (bytes (+ (* events 32) ;Each event requires a float w/ 2 word header and 64 bits
				 (* conses 16) ;Each cons takes 2 words
				 (* days 8)))) ;Each day is a slot in an array
		  (format t "~D events stored in ~D days, ~:D bytes"
			  events days bytes)
		  (return bytes))))

(defun memory-usage ()
  (let (structure calendar cells)
    (setq structure (structure-memory-usage))
    (format t "~&  Structures total: ~:D" structure)

    (format t "~&*CALENDAR*: ")
    (setq calendar (calendar-memory-usage *calendar*))

    (let* ((conses (reduce #'+ *cells-flat* :key #'length))
	   (elements (length *cells-flat*)))
      (setq cells (+ (* (+ 2 elements) 8) (* conses 16)))
      (format t "~&  Cells array: ~D elements, ~D conses in lists, ~:D bytes" elements conses cells)
    (format t "~&Total bytes ~:D~%" (+ structure calendar cells)))))
	
(defun find-slow-output (directory threshold-speed)
  (loop for file in (directory (format nil "~A/*/*/output" directory))
	do (with-open-file (stream file)
	     (loop for line = (read-line stream nil)
		   while line
		   for position = (search "%, speed " line)
		   when position
		   do (when (< (parse-integer line :start (+ position 9) :junk-allowed t) threshold-speed)
			(format t "~A: ~A: ~A~%" file
				(subseq line 0 (position #\: line))
				(subseq line (+ position 3) (position #\, line :start (+ position 9)))))))))

(defun largest-edges (&optional (n 5) (func #'<))
  (let ((largest nil))
    (flet ((do-edge (edge)
             (if (< (length largest) n)
		 (push edge largest)
	       (setq largest (cdr (sort (cons edge largest) func))))))
    (map-all-diamonds
     #'(lambda (diamond)
	 (do-edge (4vector-t (diamond-p diamond)))
	 (do-edge (4vector-t (diamond-q diamond)))))
    largest)))


;;;Read output files and process run-time data

;;If line starts with given string, return remainder
(defun line-starts-p (line tag)
  (when (string-equal line tag
		      :end1 (min (length line) (length tag)))
    (subseq line (length tag))))

;;Read an lsf output file and return total elapsed and CPU time in seconds, and number of runs
(defun parse-batch-output (file &optional multiple-OK)
  (with-open-file (stream file)
    (loop with starts = nil
	  with ends = nil
	  with cpus = nil
	  with this
	  for line = (read-line stream nil)
  	  while line			;Exit at EOF
	  do (cond ((setq this (line-starts-p line "Started at "))
		    (push (parse-lsf-time this) starts))
		   ((setq this (line-starts-p line "Results reported at "))
		    (push (parse-lsf-time this) ends))
		   ((setq this (line-starts-p line "    CPU time   :"))
		    (push (parse-integer this :junk-allowed t) cpus))
		   ((line-starts-p line "TERM_")
		    (error "Abnormal termination: ~A" line)))
	  finally
	  (unless starts (error "No data in ~A" file))
	  (unless (= (length starts) (length ends) (length cpus))
	    (error "unequal numbers of report lines in ~A" file))
	  (unless (or multiple-OK (= (length starts) 1))
	    (error "unexpected multiple report lines in ~A" file))
	  (return (values (reduce #'+ (mapcar #'- ends starts)) (reduce #'+ cpus) (length starts)))
	  )))

(defun parse-manager-output (file &optional multiple-OK)
  (with-open-file (stream file)
    (loop with times
	  with this
	  for line = (read-line stream nil)
  	  while line			;Exit at EOF
	  do (cond ((setq this (line-starts-p line "Run complete in "))
		    (push (parse-manager-time this) times))
		   ((line-starts-p line "Run FAILED")
		    (error "Failure: ~A" line)))
	  finally
	  (unless times (error "No data in ~A" file))
	  (unless (or multiple-OK (= (length times) 1))
	    (error "unexpected multiple report lines in ~A" file))
	  (return (values (reduce #'+ times) (length times)))
	  )))

(defun parse-manager-time (time)
  (multiple-value-bind (hours pos)
      (parse-integer time :junk-allowed t)
    (unless (char= (char time pos) #\:) (error "invalid manager time format"))
    (multiple-value-bind (minutes pos)
      (parse-integer time :start (1+ pos) :junk-allowed t)
      (unless (char= (char time pos) #\:) (error "invalid manager time format"))
      (let ((seconds (parse-integer time :start (1+ pos) :junk-allowed t)))
	(+ (* 60 (+ (* 60 hours) minutes)) seconds)))))

;;Parser for times in format Tue Oct 19 21:41:17 2010
(defparameter *month-table*
  #("Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec"))

(defun parse-lsf-time (time)
  (encode-universal-time
   (parse-integer (subseq time 17 19))	;Seconds
   (parse-integer (subseq time 14 16))	;Minutes
   (parse-integer (subseq time 11 13))	;Hours
   (parse-integer (subseq time 8 10))	;Date of month
   (1+ (position (subseq time 4 7) *month-table* :test #'string-equal))	;Month: 1..12
   (parse-integer (subseq time 20 24))	;Year
   ))
		  
(defun old-run-statistics (directories)
  (run-statistics #'parse-old-run directories))

(defun run-statistics (parser directories)
  (loop for directory in directories
	with elapsed and real and cpu
	do (multiple-value-setq (elapsed real cpu) (funcall parser (merge-pathnames directory batch-root-directory)))
	do (format t "~&~A: elapsed ~D, real ~D, CPU ~D" directory elapsed real CPU)
	count t into count
	sum elapsed into total-elapsed
	sum real into total-real
	sum cpu into total-cpu
	finally
	(format t "~&Averages of ~D run~:P:~%" count)
	(return (print-run-statistics (/ total-elapsed (* count 3600.0))
				      (/ total-real (* count 3600.0)) (/ total-cpu (* count 3600.0))))))

(defun parse-old-run (directory)
  (multiple-value-bind (elapsed cpu)
      (parse-batch-output (format nil "~A/bsub.out" directory))
    (values elapsed (* elapsed 64) cpu)))

;;Print statistics about some runs, given gain returns.  Returns given data and derived quantities
;;average number of processors, CPU utilization fraction, degree of parallelism.
(defun print-run-statistics (elapsed real cpu)
  (let ((processors (/ real elapsed))
	(utilization (/ cpu real))
	(degree (/ CPU elapsed)))
  (format t "~&Total CPU time: ~$~%Elapsed time: ~3$~%Total real time: ~$~%" CPU elapsed real)
  (format t "~&Average number of processors: ~$~%CPU utilization percentage: ~$~%Degree of parallelism: ~$~%"
	  processors (* 100 utilization) degree)
  (values elapsed real CPU processors utilization degree)))

(defun new-run-statistics (directories)
  (run-statistics #'parse-new-run directories))

(defun parse-new-run (directory)
  (loop for file in (directory (format nil "~A/worker-*/output" directory))
	with real and cpu
	do (multiple-value-setq (real cpu) (parse-batch-output file t))
	count t into count
	sum real into total-real
	sum cpu into total-cpu
	finally
	(format t "~&Read ~D worker-output files~%" count)
	(return (values (parse-manager-output (format nil "~A/manager-output" directory))
			total-real total-cpu))))

(defun compare-runs (old-directories new-directories)
  (multiple-value-bind (old-elapsed old-real old-CPU old-processors old-utilization old-degree)
      (old-run-statistics old-directories)
    (declare (ignore old-processors old-degree))
  (multiple-value-bind (new-elapsed new-real new-CPU new-processors new-utilization new-degree)
      (new-run-statistics new-directories)
    (declare (ignore new-processors new-degree))
    (format t "~%Advantages in~%")
    (format t "CPU time: ~$~%Elapsed time: ~$~%Total real-time: ~$~%CPU utilization: ~$~%"
	    (/ old-CPU new-CPU) (/ old-elapsed new-elapsed) (/ old-real new-real)
	    (/ new-utilization old-utilization)))))
  
	
;;Here's a test of Gaussian convolution

;;Smooth a square loop
(defun test-smooth-data-derivative ()
  (multiple-value-call
   #'maximum-second-derivative
   (smooth-data (coerce (list (make-3vector 1.0 0.0 0.0) (make-3vector 0.0 1.0 0.0)
			      (make-3vector -1.0 0.0 0.0) (make-3vector 0.0 -1.0 0.0)) 'vector)
		(make-array 4 :element-type 'double-float :initial-contents '(1.0 2.0 3.0 4.0))
		0.2 0.01)))
 
;;If you don't normalize the result in smooth-data, but then multiply by the factor of
;;2/L = 1/2 that is suppressed there because we will normalize anyway, you should get
;;maximum second derivative (/ 1 (sqrt pi) 0.2).

;;Because of the small smoothing interval, the smooth string is a diamond.  Normalization 
;;increases the length of the a' vectors by 2, so in a case you should get
;;maximum second derivative (/ (sqrt (/ 2 pi)) 0.2)




;; Some testing function added by ALE

;; ap and bp functions for checking the make-test-ab-prime function

(defun aprime (v)
  (make-3vector (sin v) 
		(* -1.0d0 (cos v))
		(sin v)))

(defun bprime (u)  
  (make-3vector (cos u)
		(sin u)
		(sin u)))

; Some functions to setup different initial conditions

; string using the a' and b' functions defined above
(defun setup-abp () 
  (make-test-ab-prime 'aprime 'bprime (* 2 pi) 100)) 

; string created a' and b' functions and randomly created bhs
(defun setup-abp-bhs () 
  (make-test-ab-prime 'aprime 'bprime (* 2 pi) 200)
  (create-bhs))

; string network created with vv initial conditions and randomly created bhs  
(defun setup-vvbhs () 
  ;(initialize :total-size 20.0 :start 15.0)
  (setf *bh-number* 2)
  ;(sort-bhs)
  (setup-vv)
  ;(plot-time-slice)
  ;(evolve-until 15.5)
  ;(format t "Evolved~%")
  (create-bhs)
  ;(plot-time-slice)
  ;(evolve-until 16.3)
)


; simulation test for the setups created above

;simulation test for a string created with a' and b'
(defun test-abp (&rest simulate-keys &key (overwrite t) (size 20.0) (start-time 10.0) &allow-other-keys)
  (apply #'simulate #'setup-abp :overwrite overwrite :size size :start start-time simulate-keys))

;simulation test for a string created with a' and b' and random bhs
(defun test-abp-bhs (&rest simulate-keys &key (overwrite t) (size 20.0) (start-time 10.0) &allow-other-keys)
  (apply #'simulate #'setup-abp-bhs :overwrite overwrite :size size :start start-time simulate-keys))


; some functions to save animations 

(defun animate-vvbhs () 
  (initialize :total-size 20.0 :start 15.0) 
  (setf *bh-number* 1)
  ;(sort-bhs)
  (setup-vv)
  (create-bhs)
  (do-animate 15.0 16.3 :filespec "../Animation/vvbhs-~D.png"))


;function to simulate a string network with blackholes attached to them. It needs the bh-number bh-time and bh-size if not the bh are not
;created properly
(defun do-vvbhs (&rest keys)
  (apply #'do-initial-1 #'setup-vv keys))
