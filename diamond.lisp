(in-package "CL-USER")

;;Basic datatype for diamonds
;;Slots which are initialized here must be cleared in either discard-diamond or allocate-diamond, so that
;;diamonds that our reuse have initial values.
(defstruct (diamond
	    (:predicate diamondp)     ;avoid conflict with edge name p
	    (:print-object print-diamond)
	    (:constructor construct-diamond))
  left					;left corner
  right					;right corner
  start					;starting position
  end					;ending position
  sw					;lower left neighbor
  nw					;upper left neighbor
  se					;lower right neighbor
  ne					;upper right neighbor
  tag					;Tag of event that created this loop of string
  a-kink-created			;Location where kink to se of this diamond was created
  b-kink-created			;Location where kink to sw of this diamond was created
  (minima nil)				;4vectors of coordinate minima and maxima
  maxima
  (pending-intersections nil)
  (right-rejoining-junction nil)	;Junction at which string begins in this diamond
  (left-rejoining-junction nil)		;Junction which string ends
  (created-right-dump-junctions nil)	;Junctions where we switched from dumping to writing to a successor
  (created-left-dump-junctions nil)	;Junctions where we switched from writing to a successor to dumping
  (received-right-dump-junctions nil)	;Junctions where our predecessors switched from us to dumping
  (received-left-dump-junctions nil)	;Junctions where our predecessors switched from dumping to us
  (packed 0 :type (unsigned-byte 64))	;Several packed fields
  (backreaction-data nil)		;Data used by backreaction and related functions
  (bh nil))                             ;BH attached to the diamond

;;Layout of the packed data
(defconstant diamond-inert-bit 0)	;diamond is inert
(defconstant diamond-predecessor-bit 2) ;This diamond received from predecessor
(defconstant diamond-final-bit 3)	;This diamond will be sent to successor
(defconstant diamond-sw-line 4) ;This diamond is directly SW of the observation diamond of create-initial-ring
(defconstant diamond-se-line 5) ;This diamond is directly SE of the observation diamond of create-initial-ring
(defconstant diamond-processed-bit 7)	;Used to prevent repetition in map-diamond-cells-diamonds
(define-constant diamond-countup-field (byte tag-count-bits 8)) ;number of oscillations undergone by a loop.

;;Packed data.  These slots work in MAKE-DIAMOND, because that is not the actual constructor.
;;They default to 0 or NIL.
(defmacro diamond-inertp (diamond)
  `(logbitp diamond-inert-bit (diamond-packed ,diamond)))

(defmacro diamond-processedp (diamond)
  `(logbitp diamond-processed-bit (diamond-packed ,diamond)))

(defmacro diamond-predecessorp (diamond)
  `(logbitp diamond-predecessor-bit (diamond-packed ,diamond)))

(defmacro diamond-finalp (diamond)
  `(logbitp diamond-final-bit (diamond-packed ,diamond)))

(defmacro diamond-sw-linep (diamond)
  `(logbitp diamond-sw-line (diamond-packed ,diamond)))

(defmacro diamond-se-linep (diamond)
  `(logbitp diamond-se-line (diamond-packed ,diamond)))

(defmacro diamond-countup (diamond)
  `(ldb diamond-countup-field (diamond-packed ,diamond)))

(mirror-images
 (defun diamond-right-bh-p (diamond) ;function to check if a neighbor is a BH or not
   (eq (diamond-e diamond) :BH)))

(defun print-diamond (diamond stream)
  (print-unreadable-object (diamond stream :identity t :type t)
    (cond ((diamond-start diamond)
	   (print-4vector (diamond-start diamond) stream))
	  (t (format stream "unitialized")))))

(eval-when (:compile-toplevel :load-toplevel :execute) ;Make structure available at compile time for clear-defstruct

(defstruct backreaction-data
  (center nil)
  (a nil)
  (a-half nil)
  (b nil)
  (b-half nil)
  (start-x 0.0 :type double-float)	;midline-intersection x value for diamond start
  ;;NIL = not computed, :NONE = does not exist, :FUTURE = the x position on the observer line is
  ;;in the future of the diamond start (the normal case), :PAST = the x position is in the past of the diamond start.
  (start-status nil)
  (end-x 0.0 :type double-float)	;as above for other tips
  (end-status nil)
  (left-x 0.0 :type double-float)
  (left-status nil)
  (right-x 0.0 :type double-float)
  (right-status nil)
  )

)

;;Get backreaction-data structure, creating it if needed
(defun get-backreaction-data (diamond)
  (or (diamond-backreaction-data diamond)
      (setf (diamond-backreaction-data diamond) (make-backreaction-data))))

;;Resource of saved diamonds.  Diamonds saved here should have contents reset to the default
;;Large maximum here permits workers to reuse diamonds between runs
(defvar diamond-resource (make-resource :name "diamond" :constructor #'construct-diamond :max 100000))

;;Allocate diamond and fill in fields
(defmacro make-diamond (&rest args) 
  (let ((diamond (gensym)))
    `(let ((,diamond (allocate diamond-resource))) ;Get a diamond with default contents
       ,@(loop for (key value) on args by #'cddr
	       collect `(setf (,(intern (format nil "DIAMOND-~A" key) "CL-USER") , diamond) , value))
       ,diamond)))

;;Reinitialize diamond slots and return to resource
(defun return-diamond (diamond)
  (setf (diamond-minima diamond) nil
	(diamond-maxima diamond) nil
	(diamond-pending-intersections diamond) nil
	(diamond-packed diamond) 0
	(diamond-tag diamond) nil
	(diamond-start diamond) nil	;Can't return 4vectors because it might be shared
	(diamond-end diamond) nil
	(diamond-bh diamond) nil)
  (mirror-images
   (setf (diamond-right diamond) nil
	 (diamond-ne diamond) nil
	 (diamond-se diamond) nil
	 (diamond-left-rejoining-junction diamond) nil
	 (diamond-created-left-dump-junctions diamond) nil
	 (diamond-received-left-dump-junctions diamond) nil
	 (diamond-a-kink-created diamond) nil))
  (return-backreaction-data diamond)
  (deallocate diamond-resource diamond)	;Give diamond back to resource
  )

;;Return objects stored in the back reaction data and clear the slots.  We don't
;;return the structure itself under the theory that if a diamond needed one perhaps it will need it again later.
(defun return-backreaction-data (diamond)
  (declare (optimize speed))
  (let ((data (diamond-backreaction-data diamond)))
    (when data
      (when (backreaction-data-center data) ;Deallocate cached 4vectors
	(deallocate 4vectors (backreaction-data-center data)))
      (mirror-images
       (when (backreaction-data-a data)
	 (deallocate 4vectors (backreaction-data-a data))
       (when (backreaction-data-a-half data)
	 (deallocate 4vectors (backreaction-data-a-half data)))))
      (clear-defstruct backreaction-data data)))) ;Clear slots in object
      


;;Sides of diamonds.  These are macros so that they can be processed by the nvector-operate system
(defmacro diamond-p (diamond)		;lower left edge
  (once-only (diamond)
    `(4vector- (diamond-left ,diamond) (diamond-start ,diamond))))

(defmacro diamond-q (diamond)		;lower left edge
  (once-only (diamond)
    `(4vector- (diamond-right ,diamond) (diamond-start ,diamond))))

(defmacro diamond-new-p (diamond)		;upper right edge
  (once-only (diamond)
    `(4vector- (diamond-end ,diamond) (diamond-right ,diamond))))

(defmacro diamond-new-q (diamond)		;upper left edge
  (once-only (diamond)
    `(4vector- (diamond-end ,diamond) (diamond-left ,diamond))))

(defmacro diamond-a (diamond)		;compatible with mirror-images
  `(diamond-p ,diamond))		;must be macro, not inline, for vector optimizers

(defmacro diamond-b (diamond)		;compatible with mirror-images
  `(diamond-q ,diamond))

(declaim (inline diamond-p-t diamond-q-t diamond-a-t diamond-b-t))
;;Just give time without making a 4vector
(defun diamond-p-t (diamond)
  (- (4vector-t (diamond-left diamond)) (4vector-t (diamond-start diamond))))
(defun diamond-q-t (diamond)
  (- (4vector-t (diamond-right diamond)) (4vector-t (diamond-start diamond))))
(defun diamond-a-t (diamond)
  (diamond-p-t diamond))
(defun diamond-b-t (diamond)
  (diamond-q-t diamond))

(mirror-images
;;4vector side of diamond normalized to given time
;;Wrap if needed, but run fast if not.
(defun diamond-a-wrap-normalized-1 (diamond time)
  (declare (double-float time)
	   (optimize speed))
  ;;Side should either be small or almost the periodicity distance
  ;;We only wrap if it is more than half that
  (if (and *total-size*
	   (locally (declare (double-float *total-size*)) ;since it is not NIL here
	     (> (3vector-length (diamond-p diamond)) (/ *total-size* 2))))
      (4vector-scale (diamond-a-wrap diamond) (/ time (diamond-a-t diamond)))
    (4vector-scale (diamond-a diamond) (/ time (diamond-a-t diamond))))) ;No wrap

(defun diamond-a-wrap-normalized (diamond)
  (let ((data (get-backreaction-data diamond)))
    (or (backreaction-data-a data)
	(setf (backreaction-data-a data)
	      (diamond-a-wrap-normalized-1 diamond 1.0)))))

(defun diamond-a-wrap-normalized-half (diamond)
  (let ((data (get-backreaction-data diamond)))
    (or (backreaction-data-a-half data)
	(setf (backreaction-data-a-half data)
	      (diamond-a-wrap-normalized-1 diamond 0.5)))))

)					;mirror images

;;Versions for contexts were we don't want mirror imaging
;;This works because mirror-images doesn't expand macros in its body
;;Probably this should be done with NOMIRROR now.
(defmacro diamond-p-wrap-normalized (&rest args)
  `(diamond-a-wrap-normalized ,@args))
(defmacro diamond-q-wrap-normalized (&rest args)
  `(diamond-b-wrap-normalized ,@args))
(defmacro diamond-p-wrap-normalized-half (&rest args)
  `(diamond-a-wrap-normalized-half ,@args))
(defmacro diamond-q-wrap-normalized-half  (&rest args)
  `(diamond-b-wrap-normalized-half ,@args))

;;link pointing given direction regardless of future or past
(defun diamond-w (diamond)
  (or (diamond-sw diamond) (diamond-nw diamond)))
(defun diamond-e (diamond)
  (or (diamond-se diamond) (diamond-ne diamond)))

;;Tell if bounding boxes overlap in three dimensions.
;;Two intervals overlap if the maximum of each is more than the
;;minimum of the other.  Otherwise they are disjoint.
;;Boxes overlap if they overlap in all coordinates.
;;In case of uncertainty, say yes.
(defun boxes-overlap (minima1 maxima1 minima2 maxima2)
  (loop for index below 3
	always (and (> (aref maxima1 index)
		       (- (aref minima2 index) fudge-coordinates))
		    (> (aref maxima2 index)
		       (- (aref minima1 index) fudge-coordinates)))))

;;Tell if two diamonds have overlapping bounding boxes.  Optimized version
(defun diamond-box-overlap (diamond-1 diamond-2)
  (declare (optimize speed (safety 0))
	   (type diamond diamond-1 diamond-2))
  (let ((minima1 (diamond-minima diamond-1))
	(maxima1 (diamond-maxima diamond-1))
	(minima2 (diamond-minima diamond-2))
	(maxima2 (diamond-maxima diamond-2))
	(fudge fudge-coordinates))
    (declare (type 3vector minima1 maxima1 minima2 maxima2))
    (macrolet ((body ()
		 `(and ,@(loop for index below 3
			       collect `(> (aref maxima1 ,index) (- (aref minima2 ,index) fudge))
			       collect `(> (aref maxima2 ,index) (- (aref minima1 ,index) fudge))))))
      (body))))


(defmacro add-final-diamond (diamond)
  `(setf (gethash ,diamond *final-diamonds*) t))
(defmacro delete-final-diamond (diamond)
  `(remhash ,diamond *final-diamonds*))
;;Call function for each final diamond
(defun map-final-diamonds (function)
  (maphash #'(lambda (diamond value)
	       (declare (ignore value))
	       (funcall function diamond))
	   *final-diamonds*))

(defun return-final-diamonds ()
  (map-final-diamonds #'return-diamond)
  (clrhash *final-diamonds*))

(defmacro add-read-diamond (diamond)
  `(setf (gethash ,diamond *read-diamonds*) t))
(defmacro delete-read-diamond (diamond)
  `(remhash ,diamond *read-diamonds*))	
;;Call function for each diamond read from dumps
(defun map-read-diamonds (function)
  (maphash #'(lambda (diamond value)
	       (declare (ignore value))
	       (funcall function diamond))
	   *read-diamonds*))


;;;Cosmological models

(define-simulate-variable *radiation-era-start* 18.0 double-float) ;Transition time for the :smooth-radiation model
(declaim (type double-float *radiation-era-start*))

(define-simulate-variable *smoothing-scale* 0.5 double-float)	;desired Hubble distance during smoothing 
					;time to start smoothing, so x_smoothing-scale = smoothing-scale/smoothing-start
(define-simulate-variable *smoothing-start* 180.0 double-float)
(define-simulate-variable *smoothing-end* 181.0)	;time to stop smoothing, so e-folds = (end - start)/scale
(declaim (type double-float *smoothing-scale*))
(declaim (type double-float *smoothing-start*))
(declaim (type double-float *smoothing-end*))

(define-simulate-variable *expansion-power* 1.0 double-float)
(declaim (type double-float *expansion-power*))

(declaim (inline conformal-Hubble-constant))


;;Hubble constant (da/dt)/a, where t is conformal time as we use here.
;;This is the fundamental quantity needed by the simulation
(defun conformal-Hubble-constant (time)
  (ecase *era*
    (:flat 0.0)
    (:radiation (/ 1.0 time))
    (:matter (/ 2.0 time))
    (:power (/ *expansion-power* time))
    (:string (/ 1.0 *smoothing-scale*))
    (:flat-smooth (cond ((< time *smoothing-start*) 0.0)
			((< time *smoothing-end*) (/ 1.0 *smoothing-scale*))
			(t 0.0)))
    (:radiation-smooth (cond ((< time *smoothing-start*)
			      (/ 1.0 time))
			     ((< time *smoothing-end*)
			      (/ 1.0 *smoothing-scale*))
			     (t
			      (/ 1.0 time))))
    (:smooth-radiation
     (if (< time *radiation-era-start*) (/ 1.0 *smoothing-scale*) ;Constant before transition time
       (/ 1.0 time)))			;Then radiation era
    (:matter-radiation (if (< time *radiation-era-start*)
			   (/ 2 (+ *radiation-era-start* time)) ;Matter era before radiation, with big bang at -t_r
			 (/ 1.0 time)))	;radiation era
    ))
     

;;Given conformal time, return scale factor.  It is undetermined by a constant
;;that should be the same for two calls to this function.
;;Thus we return a = exp(integral^t H dt), with any fixed lower limit.
(defun relative-scale-factor (time)
  (ecase *era*
    (:flat 1.0)				;Scale factor is constant
    (:radiation time)			;a ~ t
    (:matter (expt time 2))		;a ~ t^2
    (:power (expt time *expansion-power*)) ;a ~ t^*expansion-power*
    (:string (exp (/ time *smoothing-scale*))) ;a ~ exp(H t)
    (:flat-smooth (cond ((< time *smoothing-start*) 1.0)
			((< time *smoothing-end*)
			 (exp (/ (- time *smoothing-start*) *smoothing-scale*))) ;Expansion only during smoothing
			(t (exp (/ (- *smoothing-end* *smoothing-start*) *smoothing-scale*)))))
    (:radiation-smooth (cond ((< time *smoothing-start*)
			      time)
			     ((< time *smoothing-end*)
			      (* *smoothing-start* (exp (/ (- time *smoothing-start*)
							 *smoothing-scale*))))
			     (t
			      (* time *smoothing-start* (exp (/ (- *smoothing-end* *smoothing-start*)
							      *smoothing-scale*))))))
    (:smooth-radiation
     (if (< time *radiation-era-start*)	;At times earlier than transition, scale factor exponentially smaller
	 (exp (- (/ time *smoothing-scale*) 
		 (/ *radiation-era-start* *smoothing-scale*)))
       (/ time *radiation-era-start*)))	;At transition, a = 1, then grows linearly with time
    (:matter-radiation
     (if (< time *radiation-era-start*)	;At times earlier than transition, matter era w/ big bang at -t_r
	 (expt (/ (+ time *radiation-era-start*) 2 *radiation-era-start*) 2)
       (/ time *radiation-era-start*)))	;At transition, a = 1, then grows linearly with time
    ))

(defun scale-factor-ratio (time-2 time-1)
  (/ (relative-scale-factor time-2) (relative-scale-factor time-1)))

;;Amount by which horizon size (= conformal time, as used in the simulation, times scale factor)
;;is larger than physical time
(defun horizon-size-factor ()
  (ecase *era*
    (:matter 3.0)
    (:radiation 2.0)
    (:power (+ *expansion-power* 1.0))
    (:flat 1.0)))

;;This is used by default-simulation-end.
;;LAST is the last conformal time when we would like to detect the formation of loops
;;DELTA is the relative additional amount of conformal time that we would need to detect those loops if there were no
;;change in scale factor, i.e., a*last*delta is the amount of additional physical time that we want.
;;We return the conformal time at which to stop the simulation.
;;This doesn't work in general, only at the end
(defun adjust-conformal-time-end (last delta)
  (* last
     (expt (1+ delta)
	   (ecase *era*
	     ((:flat :flat-smooth) 1.0)		;just last+delta
	     ((:radiation :smooth-radiation :radiation-smooth :matter-radiation) 0.5)
	     (:matter (/ 1 3.0))
	     (:power (/ 1 (+ *expansion-power* 1.0))) 
	     (:string 0.0)))))

;; Compute comoving loop length at time of formation ti using li = l_phys/ai.with smeared l_phys = 2 T_phys. 
;; Compute proper period T_phys between conformal times t1 and t2 by approximating the scale factor as a degree 2 polynomial
;; This is exact in purely matter, radiation, matter+radiation, or flat eras, and is wrong by order (t2 - t1)^3 H^3 otherwise.
;; If smoothing makes a''(t) appear negative, treat a(t) as linear to prevent bad taylor expansion.
(defun loop-length-i (ti t1 t2)
  (let* ((r1 (scale-factor-ratio t1 ti))
	 (H1 (conformal-Hubble-constant t1))
	 (H2 (conformal-Hubble-constant t2))
	 (dt (- t2 t1))
	 (r1dot (* r1 H1))
	 (r1dotdot (max 0.0 
			(* r1 (+ (expt H1 2)
				 (/ (- H2 H1) dt))))))	;drop second order term if it's negative due to smoothing
    (+ (* 2.0 r1 dt)
       (* r1dot (expt dt 2.0))
       (* (/ r1dotdot 3.0) (expt dt 3.0)))))

    

;;Compute end of diamond from other 3 corners without expansion
(defun compute-diamond-end-flat (diamond)
  (4vector- (4vector+ (diamond-left diamond) (diamond-right diamond))
	    (diamond-start diamond)))

;;Expected numerical inaccuracy in computation of desired cosine.
;;If it is within 1 of this amount on either side, we just take it as 1.
(defparameter fudge-expansion-cosine 1e-13)

(declaim (inline Heron2))
;;Return twice the area of the triangle with the given sides, following
;;"Miscalculating Area and Angles of a Needle-like Triangle" by W. Kahan
(defun Heron2 (a b c diamond)
  (when (> c b) (rotatef b c))		;Reorder so that a >= b >= c
  (when (> b a) (rotatef a b))
  (when (> c b) (rotatef b c))
  (let ((product (* (+ a (+ b c)) (- c (- a b)) (+ c (- a b)) (+ a (- b c)))))
    (cond ((minusp product)		;Not a reasonable triangle
	   (locally (declare (optimize (speed 1))) ;Avoid warnings
		    (warn "Triangle inequality violated in HERON2 ~S~%" product)
		    (format t "~&A = ~F, B = ~F, C = ~F~%" a b c)
		    (describe diamond))
	   0.0)				;give 0
	  (t (/ (real-sqrt product)
		2.0)))))

;;Count of diamonds where we could not get null sides at the end
(defvar *compute-diamond-adjust-count* 0)
(declaim (fixnum *compute-diamond-adjust-count*))

;;Compute end of diamond, taking into account possible expansion of the universe.
;;A version of this function with debugging code is in old.lisp
(defun compute-diamond-end (diamond)
  (declare (optimize speed))
  (if (eq *era* :flat)
      (compute-diamond-end-flat diamond)
    (let* ((start (diamond-start diamond))
	   (left (diamond-left diamond))
	   (right (diamond-right diamond))
	   (c (conformal-Hubble-constant
	       (global-time (/ (+ (4vector-t left) (4vector-t right)) 2.0)))) ;expansion rate at center in global time
	   (p (4vector- left start))
	   (q (4vector- right start))
	   (plen (3vector-length p))	;Spatial length of vectors
	   (qlen (3vector-length q))
	   (pt (4vector-t p))	;Temporal length.  Might be different if vector is not null due to intersection
	   (qt (4vector-t q))
	   ;;Amount by which final time is less.  Compute using spatial lengths, giving smaller answer if different.
	   (shrinkage (* c (+ (* plen qlen) (3vector-dot p q))))
	   ;;Attempt to use temporal lengths
	   (newplen (- pt shrinkage))	;Length of p'
	   (newqlen (- qt shrinkage))	;Length of q'
	   (d (3vector- q p))				;Spatial distance d = q - p
	   (dlen (3vector-length d))
	   ;;We'd rather use the temporal length answer, so that the future edges will be null, but only
	   ;;if this does not decrease the theta angle.  When the edges are null, the velocity of the diamond is
	   ;;cos^2(theta/2).  If the edges are timelike, velocity will be less, but still we do not want to make
	   ;;these vectors more parallel.
	   (costheta (/ (+ (* plen plen) (* qlen qlen) (- (* dlen dlen)))
			2 plen qlen)) ;cos(theta) computed from past edges
	   (newcostheta (/ (+ (* newplen newplen) (* newqlen newqlen) (- (* dlen dlen)))
			   2 newplen newqlen))) ;cos(theta) computed from proposed future edges
;;      (format t "~%plen = ~S, qlen = ~S, C = ~S~%newplen = ~S, newqlen = ~S, sum = ~S, d = ~S
;;pt = ~S, qt = ~S, sum = ~S~%p.q = ~S, shrinkage = ~S~%"
;;	      plen qlen c newplen newplen (+ newplen newqlen)  (3vector-length d) pt qt (+ pt qt)
;;	      (/ (3vector-dot p q) plen qlen)  shrinkage)
      ;;If theta would increase we think again.  Numerical accuracy does not permit solving equations below
      ;;if the effect is too small, so we allow small increases in theta.
      (when (and (> newcostheta (+ costheta 1e-14)) ;Simple cases quickly
		 ;;The real trouble comes when we can't tell the difference between
		 ;;p q costheta and p q newcostheta by comparison with p^2 + q^2.
		 (> (/ (* (- newcostheta costheta) plen qlen) (+ (* plen plen) (* qlen qlen))) 1e-14))
	(locally (declare (optimize (safety 0))) ;Don't worry about fixnum overflow
	  (incf *compute-diamond-adjust-count*))
	(flet ((f (p q p2 q2)
		  (- (+ (* p p2) (* q q2)) (* (+ (* p q2) (* q p2)) costheta))))
	  ;;Avoid warnings about consing floats in error messages below.  This has to be early because the compiler
	  ;;generates some common code for consing lambda that is outside the scope of the declaration 
	  ;;if we put it inside the let
	  (without-compiler-notes
	   (let* ((p0 (- plen shrinkage)) ;Future lengths computed from past spatial lengths
		  (q0 (- qlen shrinkage))
		  (dp (- pt plen))
		  (dq (- qt qlen))
		  (a (f dp dq dp dq))	;Solve quadratic a lambda^2 + 2 b lambda + c = 0.  a>0 by construction
		  (b (f dp dq p0 q0))
		  (c (- (f p0 q0 p0 q0) (expt dlen 2))) ;negative always because p0,q0 decreases c0
		  (lambda				;We want larger root of quadratic
		    (cond ((plusp c)			;Can happen by numerical error
			   0.0)				;Just do old calculation
			  ((minusp b) (/ (- (real-sqrt (- (expt b 2) (* a c))) b) a)) ;usual case (-b+sqrt{b^2-ac})/a
			  (t (- (/ c (+ (real-sqrt (- (expt b 2) (* a c))) b))))))) ;more accuracy with -c/(b+sqrt{..})
	     (cond ((> lambda (+ 1 1e-14))
		    (error "Lambda = ~S > 1 in ~S" lambda 'compute-diamond-end))
		   ((> lambda 1)	;Small errors from numerics are OK
		    (setq lambda 1.0))
		   ((< lambda -1e-14)
		    (error "Lambda = ~S < 0 in ~S" lambda 'compute-diamond-end))
		   ((< lambda 0)
		    (setq lambda 0.0)))
	     (setq newplen (+ p0 (* dp lambda)))
	     (setq newqlen (+ q0 (* dq lambda)))))))
      (let* ((qparlen (/ (+ (* dlen dlen) (- (* newqlen newqlen) (* newplen newplen)))
			 2 dlen))		     ;Use law of cosines to compute length of component of q' along d
	     (qpar (3vector-scale d (/ qparlen dlen))) ;Component in direction parallel to d
	     ;;Numerical error can lead to lengths that don't reach. This can happen pretty easily for static
	     ;;or slow-moving diamonds, where the lengths only exactly reach to start with.
	     ;;Rather than getting an error in Heron2, if newplen+newqlen <= dlen, we just make the diamond static
	     (qperplen (if (>= dlen (+ newqlen newplen))
			   0.0
			 (/ (Heron2 dlen newqlen newplen diamond) dlen)))) ;length of perpendicular from Heron's formula
	(declare (type 4vector start left right))
	(let* ((end
		(if (plusp qperplen)	;Nonzero parallel component
		    (let ((f (3vector- q (3vector-scale d (/ (3vector-dot d q) (* dlen dlen)))))) ;Part of q perp to d
		      (prog1 (3vector+ left qpar (3vector-normalize f qperplen))		  ;Add both components
			(deallocate 3vectors f)))
		  (3vector+ left qpar)))) ;No perpendicular component (static string)
	  ;;If calculation above using pt and qt succeeded, we can make both vectoor null to within numerical error.
	  ;;If not because that would make p or q longer, one vector will have to be timelike.
	  ;;The use of MAX here makes one vector null and the other possibly timelike.
	  (prog1 (3to4vector end (max (+ (4vector-t left) newqlen) (+ (4vector-t right) newplen)))
	    (deallocate 3vectors p q d )))))))

(declaim (inline diamond-position))	;Speed up handle-possible-intersection-curved

;;Position inside diamond is given by interpolation, unless at the corner.  The following
;;are mathematically equivalent:
;;start + a*p + b*((1-a)q + aq')
;;start + a*((1-b)p + bp') + b*q
;;start + a*p + b*q + ab(p'-p)
;;start + a*p + b*q + ab(q'-q)
;;start*(1-a)*(1-b) + end*a*b + l*a*(1-b) + r*b*(1-a)
;;We use the last form, which is symmetrical and makes it obvious that
;;the extreme values of a and b give the corners of the diamond.
(defun diamond-position (diamond &key a b)
  (declare (muffle-conditions compiler-note) ;Don't worry about unused code when inline with const args
	   (double-float a b)
	   (optimize speed))
  (cond ((and (= a 0.0) (= b 0.0))
	 (diamond-start diamond))
	((and (= a 0.0) (= b 1.0))
	 (diamond-right diamond))
	((and (= a 1.0) (= b 1.0))
	 (diamond-end diamond))
	((and (= a 1.0) (= b 0.0))
	 (diamond-left diamond))
	(t
	 (let ((a1 (- 1.0 a))
	       (b1 (- 1.0 b)))
	   (4vector+ (4vector-scale (diamond-start diamond) (* a1 b1))
		     (4vector-scale (diamond-end diamond) (* a b))
		     (4vector-scale (diamond-left diamond) (* a b1))
		     (4vector-scale (diamond-right diamond) (* b a1)))))))



;;Just give time coordinate
(declaim (inline diamond-position-time))
(defun diamond-position-time (diamond &key a b)
  (declare (double-float a b))
  (let ((a1 (- 1.0 a))
	(b1 (- 1.0 b)))
    (+ (* (4vector-t (diamond-start diamond)) a1 b1)
       (* (4vector-t (diamond-end diamond)) a b)
       (* (4vector-t (diamond-left diamond)) a b1)
       (* (4vector-t (diamond-right diamond)) b a1))))


;;Install minimum and maximum coordinate values in 4-space in slots
(defun install-diamond-bounding-box (diamond)
  (declare (optimize speed))
  (let ((min (make-4vector))
	(max (make-4vector))
	(start (diamond-start diamond))
	(left (diamond-left diamond))
	(right (diamond-right diamond))
	(end (diamond-end diamond)))
    (declare (type 4vector start left right end))
    (loop for index below 4
	  do (setf (aref min index)
		   (min (aref start index) (aref left index)
			(aref right index) (aref end index)))
	  do (setf (aref max index)
		   (max (aref start index) (aref left index)
			(aref right index) (aref end index))))
    (setf (diamond-minima diamond) min
	  (diamond-maxima diamond) max)
    nil))

(defun diamond-bounding-box (diamond)
  (values (diamond-minima diamond) (diamond-maxima diamond)))



(defvar *discarded-objects*)		;If set, save everything here

;;Discard an object.  We never transmit discards anymore.
(defun discard-object (object)
  (etypecase object
    (diamond (discard-diamond object))
    (intersection (discard-intersection object))
    (4vector (discard-point object))
    )
  (setf (object-handle object) nil)	;Forget handle if any
  )

;;Clear data structures for point
(defun discard-point (point)
  (deallocate 4vectors point)		;Give back to resource
  ;;To be careful, we should always do this, instead of first checking for the
  ;;interior condition, because the check could yield variable results.
;;rewrite?
;;  (remhash point *point-owner*)
  )

;;Remove diamond from local data structures.
;;Probably this should only be called by DISCARD-OBJECT
(defun discard-diamond (diamond)
  (cond (*read-and-store*		;Not using calendars?
	 (delete-read-diamond diamond))	;Just remove from table
	(t
	 (maybe-delete-from-calendar diamond)
	 (maybe-delete-diamond-cells diamond)
	 (when (diamond-finalp diamond)	;If final diamond, delete from that set
	   (delete-final-diamond diamond))))
  ;;Remove from our data structures
  (mirror-image-let ((nw (diamond-nw diamond))
		     (sw (diamond-sw diamond)))
    ;;Discard reverse links, unless those diamonds were external
    (mirror-images
     (when (and nw (diamondp nw)) (setf (diamond-se nw) nil))
     (when (and sw (diamondp sw)) (setf (diamond-ne sw) nil))
     )
    ;;We would like to the deallocate the start point to the resource here.  Unfortunately,
    ;;it could have been generated by an intersection, and so be shared with someone else.
    )
  (return-diamond diamond) ;Clear and return to resource
  )
