(in-package "CL-USER")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; VARIABLE/PARAMETER DEFINITIONS
;;;; HELPER/SHORTHAND FUNCTIONS

(defvar *null-calendar*)
(defvar *observer-diamond*)
(defvar *observer-vector*)		;Unit null vector along the observer line
(defvar *observer-L* 0.0)		;Length of observer line
(declaim (double-float *observer-L*))
(defvar *other-vector*)			;Unit null vector in direction other than observer line
(defvar *most-sw*) ;to keep track of where we will end walk-ring
(defvar *most-se*) ;to keep track of where we will start walk-ring
(defparameter log-cutoff 1.0d-8)
(declaim (double-float log-cutoff))

(defmacro 4dot (v1 v2) ;To save typing.  Inline function would break vector optimizers
  `(4vector-dot-spacelike ,v1 ,v2))

;; Finds the 4dot of two vectors which are known to be null, of the same length, and which point in the
;; same temporal direction. This is done by finding the length of the difference between their 3vectors.
;; This avoids issues of numerics when the two vectors are perpendicular.
(defmacro 4dot-null-same-t (v1 v2)
  `(- (/ (3vector-squared-length (3vector- (4to3vector ,v1) (4to3vector ,v2))) 2)))

;; Finds the 4dot of two vectors which are known to be null and of different lengths by normalizing the second
;; to have the length of the first, applying 4DOT-NULL-SAME-T, and then scaling the result by the inverse of whatever
;; rescaling was necessary for the equal-length normalization.
(defmacro 4dot-null (v1 v2)
  (once-only (v1 v2)
    `(let* ((scale-factor (/ (4vector-t ,v1) (4vector-t ,v2)))
	    (dot-product (4dot-null-same-t ,v1 (4vector-scale ,v2 scale-factor))))
       (/ dot-product scale-factor))))

;; When in our pseudo basis, we want to store the components in a four-vector structure where the order of the elements
;; goes: u, v, c, d. To recover these elements later, we use 4vector-x, 4vector-y, etc., but want to be able to call
;; functions like 4vector-u, 4vector-v, etc. instead. The -a and -b versions are for mirroring compatibility.
(declaim (inline 4vector-u 4vector-v 4vector-c 4vector-d 4vector-a 4vector-b))
(defun 4vector-u (x) (4vector-x x))
(defun 4vector-v (x) (4vector-y x))
(defun 4vector-b (x) (4vector-x x)) ;same as 4vector-u, but can be mirrored
(defun 4vector-a (x) (4vector-y x)) ;same as 4vector-v, but can be mirrored
(defun 4vector-c (x) (4vector-z x))
(defun 4vector-d (x) (4vector-t x))

;; Find the 4dot of two vectors which are in the pseudo-orthonormal basis (hereafter called the pseudo basis), where
;; the diagonal components of the metric are zeta.
(defun 4dot-pseudo (v1 v2 zeta)
  (+ (* zeta (4vector-u v1) (4vector-v v2)) (* zeta (4vector-v v1) (4vector-u v2))
     (* (4vector-c v1) (4vector-c v2)) (* (4vector-d v1) (4vector-d v2))))

(defun diamond-center (diamond)
  (let ((data (get-backreaction-data diamond)))
    (or (backreaction-data-center data)
	(setf (backreaction-data-center data) (diamond-position diamond :a 0.5 :b 0.5)))))
      
(mirror-images
 (defvar *observer-diamond-center-a* 0.5)
 (coerce *observer-diamond-center-a* 'double-float))

;;It might not really be the center if we are offsetting the observer line
(defun observer-diamond-center ()
  (let ((data (get-backreaction-data *observer-diamond*)))
    (or (backreaction-data-center data)
	(setf (backreaction-data-center data) 
	      (diamond-position *observer-diamond* :a *observer-diamond-center-a* :b *observer-diamond-center-b*)))))


;;Dynamically bind commonly used resources, so that they work in threads.  This is faster than locking
;;or lock-free protocols to use a single resource.
;;We also bind some global variables used in backreaction
(defmacro with-local-resources (&body body)
  `(let ((3vectors (make-resource :name "3vector local" :constructor #'construct-3vector :max 100))
	 (4vectors (make-resource :name "4vector local" :constructor #'construct-4vector :max 100))
	 (diamond-resource (make-resource :name "diamond local" :constructor #'construct-diamond :max 10000))
	 *null-calendar* *observer-diamond* *observer-vector* 
	 (*observer-L* 0.0)
	 *other-vector* *most-sw* *most-se*)
     ,@body))

(declaim (inline log-small-ep atan-zero-ep))
;; For a term of the form (/ (log (/ (+ y (* x1 epsilon)) (+ y (* x2 epsilon)))) epsilon), we want to return
;; (/ (- x1 x2) y) if epsilon is below the (* LOG-CUTOFF y) threshold, and evaluate the expression normally otherwise.
(defun log-small-ep (y x1 x2 epsilon)
  (declare (double-float y x1 x2 epsilon log-cutoff))
  (let* ((numerator (abs (+ y (* x1 epsilon))))
	 (denominator (abs (+ y (* x2 epsilon))))
	 (tiny-n-p (< (/ numerator (+ (abs y) (abs (* x1 epsilon)))) fudge-coordinates))
	 (tiny-d-p (< (/ denominator (+ (abs y) (abs (* x2 epsilon)))) fudge-coordinates)))
    (cond
     ;; if epsilon is small relative to y, return the approximation and exit here
     ((< (abs epsilon) (abs (* y log-cutoff)))
      (/ (- x1 x2) y))
     ;; when both the numerator and denominator are tiny relative to the size of the terms in their sums, 
     ;;the answer should be zero.
     ((and tiny-n-p tiny-d-p) 0.0)
     ;; if only one is tiny...we should probably do something?...but maybe not.
     ((or tiny-n-p tiny-d-p)
      (without-compiler-notes		;Avoid warnings about consing floats
       (warn "LOG-SMALL-EP numerator (~e) and denominator (~e) out of scale" numerator denominator))
      (/ (real-log (/ numerator denominator)) epsilon))
     ;; otherwise, calculate as normal
     (t (/ (real-log (/ numerator denominator)) epsilon)))))

;; For a term of the form (/ (atan (* x epsilon)) epsilon), we want to return x if epsilon is zero, and evaluate the
;; expression normally otherwise.
(defun atan-zero-ep (x epsilon)
  (declare (double-float x epsilon))
  (let ((result (if (zerop epsilon) x (/ (atan (* x epsilon)) epsilon))))
    result))

;; Convert a list of sigmas (the sigma at which a segment ends) to dsigmas (the range of sigma of a segment)
(defun sigma-to-dsigma (sigma)
  (let ((dsigma (make-array (length sigma) :element-type 'double-float)))
    (setf (aref dsigma 0) (aref sigma 0))
    (loop for i from 1 below (length dsigma)
	  do (setf (aref dsigma i) (- (aref sigma i) (aref sigma (1- i)))))
    dsigma))

;; Convert a list of dsigmas (the range of sigma of a segment) to sigmas (the sigma at which a segment ends)
(defun dsigma-to-sigma (dsigma)
  (let ((sigma (make-array (length dsigma) :element-type 'double-float)))
    (setf (aref sigma 0) (aref dsigma 0))
    (loop for i from 1 below (length sigma)
	  do (setf (aref sigma i) (+ (aref dsigma i) (aref sigma (1- i)))))
    sigma))

;; Given a list of hats and a list of sigmas, convert into a list of four-vectors.
(defun hats-to-4vectors (hats sigmas)
  (map 'vector #'(lambda (hat dsigma) (4vector-scale (3to4vector hat 1.0) dsigma)) hats (sigma-to-dsigma sigmas)))

;; Given a list of four-vectors, convert into a list of hats and a list of sigmas.
(defun 4vectors-to-hats (four-vecs)
  (let ((hats (map 'vector #'(lambda (four-vec) (4to3vector (4vector-normalize-time four-vec))) four-vecs))
	(dsigmas (map '(array double-float (*)) #'(lambda (four-vec) (4vector-t four-vec)) four-vecs)))
    (values hats (dsigma-to-sigma dsigmas))))

;; Converts a given cartesian-basis vector to the pseudo basis. Note that u-basis and v-basis should be B'/2 and A'/2,
;; respectively. The c-basis and d-basis are the vectors returned by FIND-PSEUDO-BASIS.     
(defun cartesian-to-pseudo (cartesian-vector u-basis v-basis c-basis d-basis)
  (declare (optimize speed))
  (let* ((u-v-metric (4dot-null-same-t u-basis v-basis))
	 (u-part (/ (4dot v-basis cartesian-vector) u-v-metric))
	 (v-part (/ (4dot u-basis cartesian-vector) u-v-metric))
	 (c-part (4dot c-basis cartesian-vector))
	 (d-part (4dot d-basis cartesian-vector)))
    (make-4vector u-part v-part c-part d-part)))

;; For cartesian vectors which we know to be null, we can use a specialized version of CARTESIAN-TO-PSEUDO to more
;; accurately calculate the u and v parts of the pseudo-basis vector.
(defun null-cartesian-to-pseudo (cartesian-vector u-basis v-basis c-basis d-basis)
  (declare (optimize speed))
  (let* ((u-v-metric (4dot-null-same-t u-basis v-basis))
	 (u-part (/ (4dot-null v-basis cartesian-vector) u-v-metric))
	 (v-part (/ (4dot-null u-basis cartesian-vector) u-v-metric))
	 (c-part (4dot c-basis cartesian-vector))
	 (d-part (4dot d-basis cartesian-vector)))
    (make-4vector u-part v-part c-part d-part)))

;; Converts a given pseudo-basis vector to cartesian basis.
(defun pseudo-to-cartesian (pseudo-vector u-basis v-basis c-basis d-basis)
  (declare (optimize speed))
  (4vector+ (4vector-scale u-basis (4vector-u pseudo-vector)) (4vector-scale v-basis (4vector-v pseudo-vector))
	    (4vector-scale c-basis (4vector-c pseudo-vector)) (4vector-scale d-basis (4vector-d pseudo-vector))))

;; Finds the max or min value from a vector of number. Defaults to maximum.
(defun vector-mm (vector &key (max t))
  (let ((result (aref vector 0))
	(compare (if max #'< #'>)))
    (loop for i across vector do (if (funcall compare result i) (setf result i)))
    result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; CALCULATING THE BACKREACTION EFFECT IN GENERAL

;;Given two null vectors u and v, with equal time component, return two unit spacelike vectors c and d
;;giving a pseudo-normal basis, in such a way that the d component of the observer-vector vanishes.
;;This is called with u=a'/2, v=b'/2, but that isn't required.
(defun find-pseudo-basis (u-basis v-basis observer-vector)
  (declare (optimize speed))
  ;;First get two arbitrary spacelike vectors perpendicular to each other, u-basis, and v-basis.
  ;;We're concerned here with the possibility that the diamond is static or nearly static, so that
  ;;a'+b' = 0, or nearly so.  We're not very concerned with the diamond in which a'=b', because that one is
  ;;moving at the speed of light and has zero widths, and so will cause many problems everywhere.
  (let* ((p3 (3vector-cross-product u-basis v-basis)) ;Spatial vector perpendicular to spatial u and v
	 (p3len (3vector-length p3))
	 q3 p-vector q-vector)
    (cond ((zerop p3len)		;cross product vanished, because a'=-b'
	   (deallocate 3vectors p3)
	   (multiple-value-setq (p3 q3) (two-perpendicular-vectors u-basis)) ;Any two spatial vectors perpendicular to u
	   (setq p-vector (3to4vector p3 0.0)
		 q-vector (3to4vector q3 0.0)))
	  (t ;;Cross product nonzero.  If a' is nearly -b', this may be a tiny vector pointing in a nearly random
	     ;;direction, but it doesn't matter.
	   (setq p3 (prog1 (3vector-scale p3 (/ p3len)) ;A unit vector perpendicular to u and v
		      (deallocate 3vectors p3))
		 p-vector (3to4vector p3 0.0)	 ;No time component
		 ;;Get a spatial vector perpendicular to p and to b-a.  Thus q3.b = q3.a.
		 q3 (3vector-cross-product p3 (3vector- u-basis v-basis))
		 q-vector (3to4vector q3 (/ (3vector-dot u-basis q3) (4vector-t u-basis)))) ;Using v would be the same
	   (setq q-vector (prog1 (4vector-scale q-vector (/ (4vector-spacelike-length q-vector))) ;Normalize
			    (deallocate 4vectors q-vector)))))
    (deallocate 3vectors p3 q3)
    (let* ((o-p (4dot observer-vector p-vector)) ;Components of observer vector in the p and q directions
	   (o-q (4dot observer-vector q-vector))
	   (pq-len (real-sqrt (+ (expt o-p 2) (expt o-q 2))))) ;Length of observer vector projected into the p-q plane
      (if (zerop pq-len)				  ;Observer vector perpendicular to plane: no p or q component
	  (values q-vector p-vector)			  ;Use p and q unchanged
	(multiple-value-prog1
	    (values
	     ;;Unit spacelike vector C in direction of projection
	     (4vector+ (4vector-scale p-vector (/ o-p pq-len)) (4vector-scale q-vector (/ o-q pq-len)))
	     ;;Unit spacelike vector D in perpendicular direction
	     (4vector- (4vector-scale q-vector (/ o-p pq-len)) (4vector-scale p-vector (/ o-q pq-len))))
	  (deallocate 4vectors p-vector q-vector))))))

;;Stored versions of the x parameter from midline-intersection.
;;This defines FUNCTION with arguments (observer-diamond midline-vector L source-diamond past-ok) to
;;get and cache the x parameter associated with the given tip.  If the x parameter specifies the place in the past
;;of the tip, and past-ok is not given, we return NIL, but we still cache x.
(defmacro define-diamond-any-x (function slot)
  (let ((diamond-slot (intern (format nil "DIAMOND-~A" slot) (symbol-package slot)))
	(x-slot (intern (format nil "BACKREACTION-DATA-~A-X" slot) (symbol-package slot)))
	(status-slot (intern (format nil "BACKREACTION-DATA-~A-STATUS" slot) (symbol-package slot))))
    `(progn
       (defun ,function (source-diamond &optional past-ok)
	 (let* ((data (get-backreaction-data source-diamond))
		(x (,x-slot data))	;Look for cached result
		(status (,status-slot data)))
	   (unless status		   ;Don't know yet?
	     (multiple-value-setq (status x) ;Actually figure it out
	       (midline-intersection (,diamond-slot source-diamond)))
	     (setf (,x-slot data) x	;Store what we found
		   (,status-slot data) status))
	   (ecase status
	     (:none nil)
	     (:future x)
	     (:past (and past-ok x)))))
       (declaim (ftype (function (diamond &optional t) (or double-float null)) ,function))
       )))

(define-diamond-any-x diamond-left-x left)
(define-diamond-any-x diamond-right-x right)
(define-diamond-any-x diamond-start-x start)
(define-diamond-any-x diamond-end-x end)

(defvar *segment-direction*)		;Is this a-segment-change or b-...?

;; makes the initial intersection ring of the past lightcone of a point on the observer line (normally the beginning) 
;; with the string worldsheet.  Start-a and start-b are the indices into the hats and sigmas.  Caller must have
;; already bound *observer-vector*, *observer-L*, and *other-vector* to the proper values and *observer-diamond* to
;; anything.  We set *observer-diamond* to the first diamond we create.
;; Start-x (usually -L) is the starting point on the observer line.
(defun create-initial-ring (&key a-hats a-dsigmas b-hats b-dsigmas (start-a 0) (start-b 0) start-x open-loop)
  (create-worldsheet a-hats a-dsigmas b-hats b-dsigmas
		     #'(lambda (diamond direction)
			 (cond (*observer-diamond*
				(create-initial-ring-test diamond direction start-x))
			       ;;First time.  Set observer diamond, which should have been dynamically bound by caller
			       (t (setq *observer-diamond* diamond)
				  nil)))	;Always go SE after first diamond
		     :start-a start-a :start-b start-b :open-loop open-loop))
;;Tell whether to go north or south next.  T means north
(defun create-initial-ring-test (diamond direction start-x)
  (declare (optimize speed)
	   (double-float start-x))
  ;;We go NE when the right corner is in the past light cone of the source point, i.e., when the future light cone
  ;;of the right corner intersects the observer line at some x < start-x.
  (let ((x (choose-mirror-image	direction ;Position where future light cone of right intersects observer line, if any
	     (diamond-right-x diamond))))
    (and x (> start-x x))))    ;North if x from side of diamond would be before start of observer line

;; Create a calendar for tracking when the worldsheet and crossings need to be modified.
;;
;; The year is as long as the null parameter we are evolving over, and there are as many days as there are segments
;; associated with the null parameter we are keeping fixed. 
(defun make-null-calendar (L N)
  (setf *null-calendar* (make-calendar (* 2 L) :days N :name "null calendar"
				       :backups-allowed T :current-time (- L))))

;; Finds the x value at which the observer-line, parameterized by C+Vx/2 crosses the light cone
;; of the given source point.  Here C is the observer diamond center and V the midline-vector, which
;; should always be unit normalized.  We calculate the intersection using a straightforward method with the
;; observation-diamond center as a reference point.  But then we will find the result using the first result as a
;; reference point to avoid numerical errors.  If the first calculation yields |x| larger than 2L, we don't
;; bother to refine it.
;; Returns two values: the first is :NONE if there is no intersection, :FUTURE if the given point
;; of the observer line is in the future of the source point (the normal case), or :PAST if it is in the past.
;; The of a second is the x itself if any.
(defun midline-intersection (source-point)
  ;; C joins the source point to the midpoint of the observation diamond.  We solve (C+Vx/2)^2=0 for x.
  ;; Because V is null, this becomes (V . C)x + C^2 = 0
  (declare (optimize speed))
  (let* ((C (4vector- (observer-diamond-center) source-point))
	 (d (4dot C *observer-vector*)))
    ;;If D = V . C = 0, we're either on the observer line or spacelike separated from all of it.  Give NIL
    ;;If D>0, then the observer line crosses into the past light cone of the source.  NIL unless PAST-OK.
    (if (zerop d) :none						;Can't compute it
      (let ((x (- (/ (4vector-squared-length-spacelike C) d)))) ;tentative result
	(when (< (abs x) (* 2 *observer-L*))			;If not too big, refine it
	  (let* ((p1 (4vector- (4vector+ (observer-diamond-center)
					 (4vector-scale *observer-vector* (/ x 2)))
			       source-point))
		 (d2 (4dot p1 *observer-vector*)))	;Analytically the same as d, but perhaps it is more accurate now
	    (decf x (/ (4vector-squared-length-spacelike p1) d2)))) ;Refine result
	(without-compiler-notes				  ;Don't warn about returning float
	 (values (if (minusp d) :future :past) x))))))


;; Computes the correction to the null vector due to a past- or future-type crossing of a source diamond
;;   *observer-diamond* and source-diamond are the diamond in which the observation and source points, respectively, lie
;;   lo- and hi-x are the range of the null parameter we are evolving along for which this crossing applies
;;   future-sign is -1.0 for future, 1.0 for past.
;; NOTE: this implementation does not compute terms which are known to sum to zero when computing the total
;;   segment change, i.e., those terms which come from the total derivative term in the acceleration. We have
;;   dropped terms like F_g V^g in all instances.
(defun past-future-crossing (source-diamond lo-x hi-x future-sign)
  (declare (optimize speed)
	   (double-float lo-x hi-x future-sign))
  (when
      (or (fudge= hi-x lo-x fudge-coordinates)
	  ;;Return zero now if source and observer diamonds are parallel
	  (and (< (abs (4dot-null-same-t (diamond-a-wrap-normalized *observer-diamond*)
					 (diamond-a-wrap-normalized source-diamond)))
		  fudge-coordinates)
	       (< (abs (4dot-null-same-t (diamond-b-wrap-normalized *observer-diamond*)
					 (diamond-b-wrap-normalized source-diamond)))
		  fudge-coordinates)))
    (return-from past-future-crossing zero-4vector))
  (let* ((u-basis (diamond-b-wrap-normalized-half source-diamond))
	 (v-basis (diamond-a-wrap-normalized-half source-diamond))
	 (c-basis)
	 (d-basis)
	 (zeta (4dot-null-same-t u-basis v-basis))) ;The off-diagonal term in the metric
    (multiple-value-setq (c-basis d-basis) (find-pseudo-basis u-basis v-basis *observer-vector*))
    (let* ((observer-pseudo (null-cartesian-to-pseudo *observer-vector* u-basis v-basis c-basis d-basis))
	   (other-pseudo (null-cartesian-to-pseudo *other-vector* u-basis v-basis c-basis d-basis))
	   ;; Which source diamond tip we use depends on if this is a past- or future-type crossing
	   (source-tip (if (plusp future-sign) (diamond-start source-diamond) (diamond-end source-diamond)))
	   (tip-to-center (4vector- (observer-diamond-center) source-tip))
	   ;; The reference point on the (extended) observer line is where the source diamond tip's
	   ;; f.l.c. intersects that line.  In good conditions we will calculate from this reference point.
	   (x0 (if (plusp future-sign) (diamond-start-x source-diamond t) (diamond-end-x source-diamond t)))
	   (dx-hi 0.0) ;Half distances from x0 or center to hi-x and lo-x.  Set below depending on choosen strategy
	   (dx-lo 0.0) ;Must set this to some float, and we don't know yet which to use.
	   (neighbor (or (eq source-diamond (diamond-e *observer-diamond*))
			 (eq source-diamond (diamond-w *observer-diamond*))))
	   (hi-x-not-L (and (= hi-x (diamond-end-x source-diamond t)) (not neighbor)))
	   (lo-x-not-L (and (= lo-x (diamond-start-x source-diamond t)) (not neighbor)))
	   tip-to-reference)
      (declare (double-float dx-hi dx-lo))
      (if (and x0 (< (abs x0) (* 2 *observer-L*))) ;If x0 was returned and not too big
	  ;;Use x0 as reference point, allowing more accurate calculation.  Vector there is null.
	  (setf dx-hi (/ (- hi-x x0) 2)
		dx-lo (/ (- lo-x x0) 2)
		tip-to-reference (null-cartesian-to-pseudo
				  (4vector+ tip-to-center (4vector-scale *observer-vector* (/ x0 2)))
				  u-basis v-basis c-basis d-basis))
	;;Otherwise use center as reference point
	(setf tip-to-reference (cartesian-to-pseudo tip-to-center u-basis v-basis c-basis d-basis)
	      dx-hi (/ hi-x 2)
	      dx-lo (/ lo-x 2)))
      (let* ((tip-to-hi-c (+ (4vector-c tip-to-reference) (* dx-hi (4vector-c observer-pseudo))))
	     (tip-to-lo-c (+ (4vector-c tip-to-reference) (* dx-lo (4vector-c observer-pseudo))))
	     ;;OBSERVER-PSEUDO is zero in the d direction by definition, so for any point on the observer line:
	     (tip-to-observer-d (4dot d-basis (4vector- (observer-diamond-center)
							(diamond-center source-diamond))))
	     ;; Construct the F and f coefficients for the integration of h_ab,g with respect to x,
	     ;;where g is a pseudo direction.  Recall that F=f/V, where V is the relevant observer vector component
	     ;; When the observer line segment begins and ends due to diamond tip crossings, rewrite F-u and F-v for
	     ;;stability: the former when the crossing tips are separated by A, and the latter when by B.
	     (F-u (if (or (and lo-x-not-L (= hi-x (diamond-left-x source-diamond t)))
			  (and hi-x-not-L (= lo-x (diamond-right-x source-diamond t))))
		      (* (/ -2 (4vector-u observer-pseudo))
			 (real-log (1+ (/ (- (/ (* future-sign (4dot tip-to-center *observer-vector*))
						(* 2 zeta (4vector-u observer-pseudo) (diamond-a-t source-diamond)))
					     1)))))
		    (* 2 (- future-sign)
		       (log-small-ep (4vector-u tip-to-reference) dx-hi dx-lo (4vector-u observer-pseudo)))))
	     (F-v (if (or (and lo-x-not-L (= hi-x (diamond-right-x source-diamond t)))
			  (and hi-x-not-L (= lo-x (diamond-left-x source-diamond t))))
		      (* (/ -2 (4vector-v observer-pseudo)) 
			 (real-log (1+ (/ (- (/ (* future-sign (4dot tip-to-center *observer-vector*))
						(* 2 zeta (4vector-v observer-pseudo) (diamond-b-t source-diamond)))
					     1)))))
		    (* 2 (- future-sign)
		       (log-small-ep (4vector-v tip-to-reference) dx-hi dx-lo (4vector-v observer-pseudo)))))
	     (f-d (* 4 future-sign
		     (atan  (* (- tip-to-hi-c tip-to-lo-c) tip-to-observer-d)
			    (+ (expt tip-to-observer-d 2)
			       (* tip-to-hi-c tip-to-lo-c)))))
	     ;; Construct the velocity correction in each pseudo direction, then transform back to cartesian.
	     (velocity-correction-u (- (* (/ F-v zeta) (4vector-c observer-pseudo) (4vector-c other-pseudo))))
	     (velocity-correction-v (- (* (/ F-u zeta) (4vector-c observer-pseudo) (4vector-c other-pseudo))))
	     (velocity-correction-c (+ (* F-u (4vector-u other-pseudo) (4vector-c observer-pseudo))
				       (* F-v (4vector-v other-pseudo) (4vector-c observer-pseudo))
				       (* f-d (4vector-d other-pseudo))))
	     (velocity-correction-d (- (* f-d (4vector-c other-pseudo))))
	     (velocity-correction-pseudo (make-4vector velocity-correction-u velocity-correction-v
						       velocity-correction-c velocity-correction-d))
	     (velocity-correction (pseudo-to-cartesian velocity-correction-pseudo u-basis v-basis c-basis d-basis)))
	(deallocate 4vectors c-basis d-basis observer-pseudo other-pseudo tip-to-center tip-to-reference
		    velocity-correction-pseudo) ;these are cached: observer-center source-center u-basis v-basis
	velocity-correction))))

(defun past-crossing (source-diamond lo-x hi-x)
  (past-future-crossing source-diamond lo-x hi-x 1.0))
(defun future-crossing (source-diamond lo-x hi-x)
  (past-future-crossing source-diamond lo-x hi-x -1.0))

;; Computes the correction to the null vector due to a b-type crossing of a source diamond (that is, a crossing which
;; connects edges of constant b). Mirror-image produces a function which computes the correction to the null vector due
;; to an a-type crossing.
;; NOTE: this implementation does not compute terms which are known to sum to zero when computing the total
;;   segment change, i.e., those terms which come from the total derivative term in the acceleration. We have
;;   dropped terms like F_g V^g in all instances.
(mirror-images
(defun b-crossing-general (source-diamond lo-x hi-x)
  (declare (optimize speed)
	   (double-float lo-x hi-x))
  (let* ((L-b (diamond-b-t source-diamond))
	 (a-basis (diamond-a-wrap-normalized-half source-diamond)) ;e_(u), regardless of mirroring
	 (b-basis (diamond-b-wrap-normalized-half source-diamond)) ;e_(v).
	 (c-basis)
	 (d-basis)
	 (zeta (4dot-null-same-t a-basis b-basis)))
    (multiple-value-setq (c-basis d-basis) (find-pseudo-basis (nomirror b-basis a-basis) *observer-vector*))
    (let* ((observer-part-b (/ (4dot-null *observer-vector* a-basis) zeta))
	   (center-link-b (/ (4dot (4vector- (observer-diamond-center) (diamond-center source-diamond)) a-basis)
			     zeta))
	   (end-to-center-b (- center-link-b L-b))
	   (start-to-center-b (+ center-link-b L-b))
	   (F-b (* 2 (+ (log-small-ep end-to-center-b (/ hi-x 2) (/ lo-x 2) observer-part-b)
				  (log-small-ep start-to-center-b (/ lo-x 2) (/ hi-x 2) observer-part-b))))
	   (observer-pseudo (null-cartesian-to-pseudo *observer-vector* (nomirror b-basis a-basis) c-basis d-basis))
	   (other-pseudo (null-cartesian-to-pseudo *other-vector* (nomirror b-basis a-basis) c-basis d-basis))
	   ;(velocity-correction-b 0)
	   (velocity-correction-a (- (* (/ F-b zeta) (4vector-c observer-pseudo) (4vector-c other-pseudo))))
	   (velocity-correction-c (* F-b (4vector-c observer-pseudo) (4vector-b other-pseudo))))
      (prog1 (4vector+ (4vector-scale a-basis velocity-correction-a) ;(4vector-scale b-basis velocity-correction-b)
		       (4vector-scale c-basis velocity-correction-c))
	(deallocate 4vectors c-basis d-basis)
	;;These are cached: a-basis b-basis v-basis u-basis
	))))

;; In the event that the source diamond's left tip is the observer diamond's right tip,
;; we can use what we call the "direct method" to construct the numerator and denominator
;; of the logarithm terms. This is done by moving directly along diamond surfaces.
(defun b-crossing-direct (source-diamond lo-x hi-x)
  (declare (optimize speed)
	   (double-float lo-x hi-x))
  (let* ((L-b (diamond-b-t source-diamond))
	 (a-basis (diamond-a-wrap-normalized-half source-diamond)) ;e_(u), regardless of mirroring
	 (b-basis (diamond-b-wrap-normalized-half source-diamond)) ;e_(v).
	 (c-basis)
	 (d-basis)
	 (zeta (4dot-null-same-t a-basis b-basis)))
    (multiple-value-setq (c-basis d-basis) (find-pseudo-basis (nomirror b-basis a-basis) *observer-vector*))
    (let* ((observer-part-b (/ (4dot-null *observer-vector* a-basis) zeta))
	   (b-correction-q (4vector-parallel-p (diamond-a *observer-diamond*) *observer-vector* fudge-coordinates))
	   ;(hi-y (/ (funcall (if b-correction-q #'+ #'-) hi-x *observer-L*) 2)) ;these apparently broke the optimizer? complained of type ambiguities
	   ;(lo-y (/ (funcall (if b-correction-q #'+ #'-) lo-x *observer-L*) 2))
	   (hi-y (if b-correction-q (/ (+ hi-x *observer-L*) 2) (/ (- hi-x *observer-L*) 2)))
	   (lo-y (if b-correction-q (/ (+ lo-x *observer-L*) 2) (/ (- lo-x *observer-L*) 2)))
	   (start-to-edge-b
	    (if b-correction-q
		(- (/ (4dot-null (diamond-b *observer-diamond*) a-basis) (* 2 zeta)))
	      (/ (4dot-null (diamond-a *observer-diamond*) a-basis) (* 2 zeta))))
	   (end-to-edge-b (- start-to-edge-b (* 2 L-b)))
	   (F-b (* 2 (+ (log-small-ep end-to-edge-b hi-y lo-y observer-part-b)
			(log-small-ep start-to-edge-b lo-y hi-y observer-part-b))))
	   (observer-pseudo (null-cartesian-to-pseudo *observer-vector* (nomirror b-basis a-basis) c-basis d-basis))
	   (other-pseudo (null-cartesian-to-pseudo *other-vector* (nomirror b-basis a-basis) c-basis d-basis))
	   ;(velocity-correction-b 0)
	   (velocity-correction-a (- (* (/ F-b zeta) (4vector-c observer-pseudo) (4vector-c other-pseudo))))
	   (velocity-correction-c (* F-b (4vector-c observer-pseudo) (4vector-b other-pseudo))))
      (prog1 (4vector+ (4vector-scale a-basis velocity-correction-a) ;(4vector-scale b-basis velocity-correction-b)
		       (4vector-scale c-basis velocity-correction-c))
	(deallocate 4vectors c-basis d-basis)
	;;These are cached: a-basis b-basis v-basis u-basis
	))))

(defun b-crossing (source-diamond lo-x hi-x)
  (if (eq (diamond-left source-diamond) (diamond-right *observer-diamond*))
      (b-crossing-direct source-diamond lo-x hi-x)
    (b-crossing-general source-diamond lo-x hi-x)))
)

;; For a given observation and source diamond, and for the minimum value of some null parameter x of the observation
;; diamond's midline which sees the source diamond, we find the correction to the appropriate null vector of the
;; observation diamond due to the source diamond.
(defun find-correction (source-diamond min-x)
  (declare (optimize speed)
	   (double-float min-x))
  ;; If the direction of integration is common to one or the other of the source null vectors,
  ;; the acceleration vanishes. (This relies on the total derivative vanishing!)
  (when (or (< (abs (4dot-null-same-t (diamond-a-wrap-normalized source-diamond) *observer-vector*))
	       fudge-coordinates)
            (< (abs (4dot-null-same-t (diamond-b-wrap-normalized source-diamond) *observer-vector*))
	       fudge-coordinates))
    (return-from find-correction zero-4vector))
  (flet ((hi-x (x &optional y)	;find minimum of L, x, and y if they are not NIL
	       (if y
		   (without-compiler-notes ;avoid unused code warnings when y not supplied
		    (if x (min *observer-L* x y) (min *observer-L* y)))
		 (if x (min *observer-L* x) *observer-L*))))
    (declare (inline hi-x))
    (if (diamond-sw source-diamond)
	(if (diamond-se source-diamond)	;SW&SE: past crossing ends with right or left tip
	    (past-crossing source-diamond min-x (hi-x (diamond-left-x source-diamond) (diamond-right-x source-diamond)))
	  ;;SW&NE: b-crossing ends with left tip
	  (b-crossing source-diamond min-x (hi-x (diamond-left-x source-diamond))))
      (if (diamond-se source-diamond)	;SE&NW: a-crossing ends with right tip
	  (a-crossing source-diamond min-x (hi-x (diamond-right-x source-diamond)))
	;;NW&NE: future-crossing ends with diamond end
	(future-crossing source-diamond min-x (hi-x (diamond-end-x source-diamond))))
      )))

;;Add source diamond to calendar if it has both north neighbors, neither is the observer diamond, and the X
;;value at which the source diamond would be advanced as before the limiting X of the observer diamond.
(defun maybe-null-calendar-add (source-diamond)
  ;;If this diamond has both its successors already and is not adjacent to the observer diamond
  (when (and (diamond-ne source-diamond) (not (eq (diamond-ne source-diamond) *observer-diamond*))
	     (diamond-nw source-diamond) (not (eq (diamond-nw source-diamond) *observer-diamond*)))
   (let ((x-cross (diamond-end-x source-diamond)))
     (when (and x-cross			;If observer crosses into future light cone of end
		(< x-cross *observer-L*)) ;before end of observer diamond
       (calendar-add *null-calendar* x-cross source-diamond)))))

;; Walks around an intersection ring, accumulating diamonds' contributions to change to A'/B'.  Additionally finds
;; those diamonds with two north neighbors, and adds them to the calendar if their futuremost point's future lightcone
;; intersects the observation point's null future on the observation diamond.
(defun walk-ring ()
  (let ((source-diamond *most-se*)
        (delta-null-vector (make-zero-4vector))) ;The correction to the null vector is initially zero.
    (loop do (4vector-incf delta-null-vector (find-correction source-diamond (- *observer-L*)))
	  do (maybe-null-calendar-add source-diamond) ;Set up to advance this diamond if needed
	  until (eq source-diamond *most-sw*) ;Stops once we have processed *most-sw*
	  do (setf source-diamond (diamond-e source-diamond)))
    delta-null-vector))                 ;Returns the accumulated correction to the null vector thus far

;; Once we have walked an intersection ring, we evolve the observation point along the null midline, changing the
;; nature of the intersection ring as we go. As we add each diamond to the worldsheet, we must do two things:
;;   1) Find the totality of its change to the null vector, and
;;   2) Check to see if the northern neighbors of the diamond we just removed now need to be in the calendar.
(defun evolve-observer ()
  (loop with d-null-vector = (make-zero-4vector) ;The correction to the null vector is initially zero.
	until (calendar-empty-p *null-calendar*)
	do
	(multiple-value-bind (add-x old-diamond) (calendar-next *null-calendar*) ;Diamond to advance
	  (let* ((w-diamond (diamond-nw old-diamond))
		 (e-diamond (diamond-ne old-diamond))
		 (start (diamond-end old-diamond))
		 (left (diamond-end w-diamond))
		 (right (diamond-end e-diamond))
		 (end (4vector+ left right (4vector- start)))
		 (new-diamond (make-diamond ;Create the new diamond
			       :start start
			       :left left
			       :right right
			       :end end)))
	    ;;If the source diamond was *most-sw* or *most-se*, set the variable to old-diamond's ne or nw neighbor
	    (mirror-images (when (eq old-diamond *most-sw*) (setf *most-sw* (diamond-ne old-diamond))))
	    ;; We are done with the source diamond. Return it.
	    (return-diamond old-diamond)
	    ;; Link the new diamond to the source's northern neighbors, then remove links between the source and its
	    ;; northern neighbors.
	    (mirror-images 
	     (setf (diamond-ne w-diamond) new-diamond)
	     (setf (diamond-sw new-diamond) w-diamond)
	     (setf (diamond-sw e-diamond) nil))
	    ;; Compute the corrections to the null vector due to the newly added diamond, as well as the new crossing
	    ;; types of its west and east neighbors.
	    (4vector-incf d-null-vector
			  (4vector+ 
			   (find-correction new-diamond add-x)
			   (find-correction w-diamond add-x)
			   (find-correction e-diamond add-x)))
	    ;; Since our diamond has been advanced, neighboring diamonds may belong in calendar
	    (mirror-images (maybe-null-calendar-add e-diamond))))
	finally (return d-null-vector)))

;; Finds the most-sw or most-se diamond in a direct line from the given diamond
(mirror-images

(defun find-most-sw (diamond)
  (let ((test-diamond (diamond-sw diamond)))
    (loop while (diamond-sw test-diamond)
          do (setf test-diamond (diamond-sw test-diamond)))
    test-diamond))

;; Finds the correction to a single segment of a. Returns a value in units of NGmu.
;; For b-segment-change, b-hats and b-dsigmas come first in argument list
(defun a-segment-change (a-hats a-dsigmas b-hats b-dsigmas iter-a
				&key ((:observer-a *observer-diamond-center-a*) *observer-diamond-center-a*))
  (mirror-images (check-dsigmas a-dsigmas))
  (let ((d-a-prime (make-zero-4vector))
	(*other-vector* (3to4vector (aref a-hats iter-a) 1.0))) ;Unit null vector in a direction.  Doesn't depend on b.
    (dotimes (iter-b (length b-hats))	;Loop over b's that combine with this a.
      (let* ((*observer-vector* (3to4vector (aref b-hats iter-b) 1.0)) ;The direction which the observer point moves
	     (*observer-L* (/ (aref b-dsigmas iter-b) 2.0))	      ;The length of the observer line segment
	     (*segment-direction* :a)
	     *observer-diamond*)				      ;Must bind here
	(create-initial-ring :a-hats a-hats :a-dsigmas a-dsigmas ;Sets *observer-diamond*
			     :b-hats b-hats :b-dsigmas b-dsigmas
			     :start-a iter-a :start-b iter-b :start-x (- *observer-L*))
	(make-null-calendar *observer-L* (length a-hats))
        ;; Find and set the *most-sw* and *most-se* diamonds in the ring. These will be the ending and starting
        ;; points, respectively, of the walk-ring function.
	(mirror-image-let ((*most-sw* (find-most-sw *observer-diamond*)))
	  (4vector-incf d-a-prime (walk-ring))	    ;Find the total correction from all diamonds existing now.
	  (4vector-incf d-a-prime (evolve-observer))		;Evolve for one cycle
	  (loop for d = (diamond-e *observer-diamond*) then next ;Return all diamonds remaining in ring
		until (eq d *observer-diamond*)
		for next = (diamond-e d) ;Must get next diamond before returning D
		do (return-diamond d)))))
    d-a-prime))
)                                       ;mirror-images

;; Calculates the effect of gravitational backreaction on a loop specified by the input a and b hats and dsigmas.
;; We have let the effect accumulate over many oscillations, and so the scaling parameter is NGmu.
;;
;; The first loop steps over all segments of a (b), finding the correction to each by evolving through all combinations
;; of this segment with all segments of b (a) (handled in the second loop).
;;
;; Returns a set of a and b hats and sigmas which represent the modified loop.
(defun gravitational-backreaction (a-hats a-dsigmas b-hats b-dsigmas NGmu)
  (mirror-images (check-dsigmas a-dsigmas))
  (mirror-image-let ((new-a-hats (make-array (length a-hats)))
		     (new-a-dsigmas (make-array (length a-hats) :element-type 'double-float)))
    ;; We want to modify each segment of a, so loop over them
    (mirror-images
     (thread-dispatch (length a-hats) *threads*
       #'(lambda (iter-a)
	   (with-local-resources
	       ;; As we loop over a, we must evolve each a segment over one oscillation, accumulating changes as we go
	    (let ((d-a-prime (4vector-scale (a-segment-change a-hats a-dsigmas b-hats b-dsigmas iter-a) NGmu)))
	      ;;(format t "~S~%" (4dot (3to4vector (aref a-hats iter-a) 1) d-a-prime))
	      ;;(if (< 0 (4vector-t d-a-prime)) (format t "~S-~D: ~f~%" :a iter-a (4vector-t d-a-prime)))
	      ;;(format t "starting ~S-~D~%" :a iter-a)
		 ;; With the correction to this segment of a-prime, we may find the new a-prime, as well as the new
		 ;; a-hat and a-d-sigma values.
		 ;; The a-hats must be explicitly set to unit length.  Just dividing by 1+Delta as we did previously  
		 ;; introduces errors in the length of order (NGmu)^2, because adding the spatial parts of A and delta-A
		 ;; changes the length instead of just rotating the vector.
		 (setf (aref new-a-hats iter-a) (3vector-normalize
						 (3vector+ (aref a-hats iter-a) (4to3vector d-a-prime)))
		       (aref new-a-dsigmas iter-a) (* (aref a-dsigmas iter-a) (1+ (4vector-t d-a-prime))))
		 ;;(format t "Old dsigma: ~S | New dsigma: ~S~%" (aref a-dsigmas iter-a) (aref new-a-dsigmas iter-a))
		 ;;(format t "New a' ~S, a-hat ~S~%" (aref new-a-primes iter-a) (aref new-a-hats iter-a))
		 )))))
    (values new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas)))

;;Acceleration of loop according to backreaction calculation
;;The loop velocity is v = sum(da_i)/L.  In the initial conditions it is zero.  Thus a = dv/dt = 1/L sum(dda_i/dt)
;;and dL/dt does not contribute.
;;We return separate values for a's and b's
;;INTERPOLATE splits each segment into at least the given number of subsegments and uses the center of each one
;;MIN-LENGTH splits until points are at most this far apart in u or v
(defun backreaction-acceleration (a-hats a-dsigmas b-hats b-dsigmas &key (interpolate 1) min-length exponential)
  (mirror-images (check-dsigmas a-dsigmas))
  (assert (not (and min-length exponential)))
  (mirror-image-let ((a-accel (make-zero-4vector))
		     (a-parts (make-array (length a-hats)))
		     (a-counts (make-array (length a-hats) :initial-element 0)))
    (mirror-images
     (thread-dispatch (length a-hats) *threads*
       #'(lambda (iter-a)
	   (with-local-resources
	     (let ((this-accel (make-zero-4vector)) ;Accumulate just this segment in this thread
		   (count interpolate))		    ;at least this many segments
	       (when min-length
		 (setq count (max count (ceiling (/ (aref a-dsigmas iter-a) min-length)))))
	       (dotimes (i count)
		 (let ((pos (if exponential
				(- 1.0 (/ 1.0 (expt 2 (1+ i)))) ;1/2, 3/4, 7/8, ...
			      (/ (1+ (* 2.0 i)) 2 count))) ;1/N, 2/N, ...
		       (weight (if exponential 
				   (cond ((= count 1) 1.0)
					 ((zerop i) (/ 5.0 8)) ;center of segment gets 5/8
					 ;;Then 3/16, 3/32, 3/64, but last is duplicated
					 (t (/ 3.0 (expt 2 (+ 3 (min i (- count 2)))))))
				 (/ 1.0 count) ;1/N
				 )))
		   ;;dda_i/da is da'_i/dt scaled by dsigma_i
		   (4vector-incf this-accel
				 (4vector-scale
				  (a-segment-change a-hats a-dsigmas b-hats b-dsigmas iter-a ;rate of change of a'
						    :observer-a pos)
				  (* (aref a-dsigmas iter-a) weight)))))
	       (setf (aref a-parts iter-a) this-accel) ;Store for later
	       (incf (aref a-counts iter-a) count))
	     (format t ".") 
	     (force-output)
	     ))))
    (let ((total-count (+ (reduce #'+ a-counts) (reduce #'+ b-counts))))
      (format t "~&~D evaluations, ~$ per segment on average~%"  total-count (/ total-count (+ (length a-hats) (length b-hats)))))
    (mirror-images			;Add up accelerations after multithreading is done to avoid races
     (dotimes (i (length a-hats))
       (4vector-incf a-accel (aref a-parts i))))
    (values a-accel b-accel (4vector-euclidian-length (4vector- b-accel a-accel))
	    (/ (+ (4vector-t a-accel) (4vector-t b-accel)) -2 ;Average loss of length
	       (/ (total-sigma a-dsigmas) 2)))		      ;Divide by oscillation time to get Gamma
    ))
		  

;;Split hats into pieces so they can be backreacted separately.  The function is called with
;;this hat, the next hat, and the amount of sigmas spent at this hat.  It returns the number of 
;;pieces to use
(defun split-hats (hats dsigmas function)
  (check-dsigmas dsigmas)
  (loop for index below (length hats)
	for hat = (aref hats index)
	for next-hat = (aref hats (next-array-index-wrapping index hats))
	for dsigma = (aref dsigmas index)
	for count = (funcall function hat next-hat dsigma) ;Find number of pieces
	when (> count 1)
	do (format t "Splitting a segment of length ~S before a kink of angle ~S into ~D pieces~%"
		   dsigma (spherical-angle hat next-hat) count)
	nconc (make-list count :initial-element hat) into new-hats ;Copies of hat
	nconc (make-list count :initial-element (/ dsigma count)) into new-dsigmas ;divide sigma
	finally
	(let ((new-hats-1 (coerce new-hats 'vector))
	      (new-dsigmas-1 (coerce new-dsigmas '(vector double-float))))
	  (format t "Split ~D hats into ~D~%" (length hats) (length new-hats-1))
	  (return (values new-hats-1 new-dsigmas-1)))))

;;Split each into N
(defun split-hats-n (hats dsigmas n)
  (split-hats hats dsigmas #'(lambda (&rest ignore) (declare (ignore ignore)) n)))

;;Split up segments so they are no longer than threshold
(defun split-long-hats (hats dsigmas min-length)
  (split-hats hats dsigmas 
	      #'(lambda (hat next-hat dsigma)
		  (declare (ignore hat next-hat))
		  (ceiling (/ dsigma min-length)))))

;;Split hat that comes before (in increasing u or v) a large-angle kink
(defun split-large-kinky-hats (hats dsigmas n angle-threshold)
    (split-hats hats dsigmas 
	      #'(lambda (hat next-hat dsigma)
		  (declare (ignore dsigma)) ;Split into N if over threshold
		  (if (> (spherical-angle next-hat hat) angle-threshold) n 1))))

;;Split hats that come before kinks so that the number of hats is enough to distribute over the kink
;;to give angles less than the threshold
(defun split-kinky-hats (hats dsigmas angle-threshold)
  (split-hats hats dsigmas 
	      #'(lambda (hat next-hat dsigma)
		  (declare (ignore dsigma))
		  ;;Split into enough pieces that they could be distributed over the angle in kinks
		  ;;less than the threshold
		  (ceiling (spherical-angle next-hat hat) angle-threshold))))

;;Split so that segment is not too large based on all subsequent kinks.  Each segment should be no larger than
;;sqrt{6 d a} + a where d is the distance from this segment to a kink of angle K and a = 2 c L/K.
;;Unfortunately, this does not fit into the paradigm above.
(defun split-hats-c (hats dsigmas c)
  (check-dsigmas dsigmas)
  (let ((new-hats nil)
	(new-dsigmas nil)
	(kink-data nil)			;A list of (kink-angle . distances-to-kink) in order of increasing angle
	(loop-length (total-sigma dsigmas)))
    (flet ((add-kink-data (angle)	;Account new kink at zero distance
	     (loop while (and kink-data (<= (caar kink-data) ;old angle
					    angle))	     ;loop until we find an old angle that is larger
		   do (pop kink-data))
	     (setq kink-data (cons (cons angle 0.0) kink-data))) ;Put new entry on front
	   (advance-kink-data (dsigma)				    ;Increase all distances by length of new segment
	     (loop for item in kink-data			    ;For those that remain
		   do (incf (cdr item) dsigma)))		    ;account additional distance
	   (max-dsigma ()			   ;Maximum segment size that can come next
;;	      (format t "~&data ~S~%" kink-data)
	     (loop for (old-angle . d) in kink-data ;Scan old data.  There always is at least one kink.
		   for coefficient = (/ (* 2 c loop-length) old-angle)
;;		   do (format t "~S@~S " coefficient d)
		   minimize (+ (sqrt (* d coefficient)) coefficient))))
      ;;First pass: arrive at the beginning of the list with kink-data set up for the last segment
      (loop for index from (1- (length hats)) downto 0
	    for hat = (aref hats index)
	    for next-hat = (aref hats (next-array-index-wrapping index hats))
	    do (add-kink-data (spherical-angle hat next-hat)) (advance-kink-data (aref dsigmas index)))
      ;;Now actually do the work
      (loop for index from (1- (length hats)) downto 0
	    for hat = (aref hats index)
	    for next-hat = (aref hats (next-array-index-wrapping index hats))
	    for angle = (spherical-angle hat next-hat)
	    for dsigma = (aref dsigmas index)
	    ;;Put new angle on front of list.  This only happens once becaus eafter this new segments have angle 0.
	    do (add-kink-data angle)
	    do
	    (loop for max = (max-dsigma) ;maximum allowable size
		  until (<= dsigma max) ;return if small enough now
		  do (format t " dsigma = ~S, max = ~S ~% " dsigma max)
		  do (setq max (min max (/ dsigma 2))) ;don't split off more than half
		  do (format t " splitting index ~D old size ~S new ~S~%" index dsigma max)
		  do (push hat new-hats) (push max new-dsigmas)
		  do (advance-kink-data max) ;Advance lengths by length of new segment
		  do (decf dsigma max))	;decrease the size of segment being processed.  Continue splitting
	    ;;Here when segment does not need (more) splitting
	    (push hat new-hats) (push dsigma new-dsigmas)
	    (advance-kink-data dsigma))) ;Account length of final segment added
	(setq new-hats (coerce new-hats 'vector)
	      new-dsigmas (coerce new-dsigmas '(vector double-float)))
	(format t "Split ~D hats into ~D~%" (length hats) (length new-hats))
	(describe-split-hats new-hats new-dsigmas)
	(values new-hats new-dsigmas)))

;;Describe situation after splitting
(defun describe-split-hats (hats dsigmas &optional (threshold (/ pi 2)) (n 3))
  (loop for index from (1- (length hats)) downto 0
	    for hat = (aref hats index)
	    for next-hat = (aref hats (next-array-index-wrapping index hats))
	    for angle = (spherical-angle hat next-hat)
	    when (>= angle threshold)
	    do
	    (format t "Kink of angle ~S; segment sizes" angle)
	    (describe-split-hats-segments dsigmas index n)))

(defun describe-split-hats-segments (dsigmas start-index n)
  (loop repeat n
	for index = start-index then (previous-array-index-wrapping index dsigmas)
	do (format t " ~S" (aref dsigmas index)))
  (terpri))


;; Given a set of a and b hats and dsigmas, performs ITERATIONS instances of backreaction of strength NGMU.
(defun do-backreaction (a-hats a-dsigmas b-hats b-dsigmas &key (NGmu 0.0001) (iterations 100) (print-at 20)
                               (out "backreaction") (track-kinks nil) (offset-count 0))
  (loop for old-L = nil then L
	for L = (/ (+ (total-sigma a-dsigmas) (total-sigma b-dsigmas)) 2)
	for count from 0 do
	;; At the first, final, and PRINT-AT iterations, we dump all of the HATS and DSIGMAS data to file
	(when (and out (or (zerop (mod count print-at)) (= count iterations)))
	  (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas
			      (format nil "~A-hs-~4,'0D.dat" out (+ count offset-count))))
	(when out
	  (with-open-file  ;at each step, we dump the loop's Gamma, length, total angular momentum, and kinks(?) to file
	      (s (format nil "~A-gljk.dat" out) :direction :output :if-does-not-exist :create :if-exists :append)
	    (format s "~3$ ~8$ ~8$" (* (/ 2 NGmu) (- 1 (/ L (or old-L L)))) L
		    (3vector-length (loop-angular-momentum-ab a-hats (dsigma-to-sigma a-dsigmas)
							      b-hats (dsigma-to-sigma b-dsigmas))))
	    (if track-kinks		;only dump the N largest kinks if TRACK-KINKS is set to N
		(mirror-image-let
		    ((a-kinks (hat-angles a-hats :largest track-kinks)))
		  (mirror-images (dotimes (i track-kinks) (format s " ~8$" (nth i a-kinks))))))
	    (format s "~%")))
	until (= count iterations)
	do (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
	     (gravitational-backreaction a-hats a-dsigmas b-hats b-dsigmas NGmu))
	do (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
	     (close-4vectors a-hats a-dsigmas b-hats b-dsigmas))))


;; We compute the part of the angular momentum associated with a set of hats and sigmas (used to be called
;; COMPUTE-A-ANGULAR-MOMENTUM)
(defun compute-hats-angular-momentum (hats sigmas)
  (let*((vectors (integrate-hats hats sigmas))
	(count (length hats))
	(loop-length (aref sigmas (1- count)))
	(angular-momentum zero-3vector))
    ;;(print (aref a-vectors 0))
    (loop for index below count
	  for previous = (previous-index-wrapping index count)
	  for next = (next-index-wrapping index count)
	  for delta-sigma = (aref sigmas index) then (- (aref sigmas index) (aref sigmas previous))
	  do (setf angular-momentum (3vector+ angular-momentum 
					      (3vector-scale (3vector-cross-product (aref vectors index)
										    (aref hats index))
							     delta-sigma))))
    ;;We divide by the length of the loop squared to make it a dimensionless quantity
    (3vector-scale angular-momentum (/ 1.0 (expt loop-length 2.0)))))

;; Obtains the total angular momentum of a loop from the a and b hats and sigmas.
;; If there is some maximum angular momentum we wish to measure relative to, divide by that.
(defun loop-angular-momentum-ab (a-hats a-sigmas b-hats b-sigmas &optional max-J)
  (3vector-scale (3vector+ (compute-hats-angular-momentum a-hats a-sigmas)
			   (compute-hats-angular-momentum b-hats b-sigmas))
		 (if max-J (/ (* 4 max-J)) (/ 4))))

;;Call first close-sigmas than close-hats.
(defun close-4vectors (a-hats a-dsigmas b-hats b-dsigmas	
		      &optional (tolerance close-hats-default-tolerance) (max-steps 10))
  (multiple-value-bind (new-a-dsigmas new-b-dsigmas) (close-sigmas a-dsigmas b-dsigmas)
    (mirror-image-let ((new-a-hats (close-hats a-hats :dsigmas new-a-dsigmas
					       :tolerance tolerance :max-tries max-steps)))
      (values new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas))))

;;Rescale dsigmas so that the total a and total b are the same.
(defun close-sigmas (a-dsigmas b-dsigmas)
  (mirror-images (check-dsigmas a-dsigmas))
  (mirror-image-let* ((a-length (total-sigma a-dsigmas))
		      (period (nomirror (/ (+ a-length b-length) 2)))
		      (new-a-dsigmas (make-array (length a-dsigmas) :element-type 'double-float)))
    (mirror-images
     (dotimes (index (length a-dsigmas)) ;Rescale to make lengths match
       (setf (aref new-a-dsigmas index) (* (aref a-dsigmas index) (/ period a-length)))))
    (values new-a-dsigmas new-b-dsigmas)))

;;Solution using differential equation techniques

;;How elements are stored in flat arrays

;;Number of segments of a.  This must be bound so that the derivatives function knows how to split its argument
(defvar *na*)				

;;AB is :A or :B.  *NA* must be set.
(defun flatten-hats-index (ab index component)
  (+ (* (+ (ecase ab (:a 0) (:b *na*)) index) 3)
     component))

;;Convert to array with the components of da and db.  *NA* must be set already.
(defun flatten-hats (a-hats a-dsigmas b-hats b-dsigmas)
  (mirror-image-let* ((n-a (length a-hats)))
    (let* ((ny (* 3 (+ n-a n-b)))
	   (y (make-array ny :element-type 'double-float)))
      (mirror-images
       (dotimes (index n-a)
	 (dotimes (component 3)
	   (setf (aref y (flatten-hats-index :a index component))
		 (* (aref a-dsigmas index) (3vector-component (aref a-hats index) component))))))
      y)))

;;Convert flat array into hats and dsigmas and go to rest frame.  *NA* must be set
(defun unflatten-hats (y)
  (let* ((ny (length y))
	 (nab (/ ny 3))
	 (n-a *na*)
	 (n-b (- nab n-a)))
    (mirror-image-let ((a-hats (make-array n-a))
		       (a-dsigmas (make-array n-a :element-type 'double-float)))
      (mirror-images
       (loop for index below n-a
	     for d-a = (make-3vector)
	     do (dotimes (component 3)
		  (setf (3vector-component d-a component) (aref y (flatten-hats-index :a index component))))
	     do (setf (aref a-hats index) (3vector-normalize d-a))
	     do (setf (aref a-dsigmas index) (3vector-length d-a))))
;;      (mirror-images (format t "~A: ~S~%" :a (total-hats-dsigmas a-hats a-dsigmas)))
      (rest-frame-close-hats a-hats a-dsigmas b-hats b-dsigmas))))

;;First close-sigmas, then boost to rest frame, then remove any residual with close-4vectors.
;;We call close-sigmas twice, once to get a well-closed loop to boost, and a second time after boosting
;;to clean up the results
(defun rest-frame-close-hats (a-hats a-dsigmas b-hats b-dsigmas	
				     &optional (tolerance close-hats-default-tolerance) (max-steps 10))
  (multiple-value-setq (a-dsigmas b-dsigmas) (close-sigmas a-dsigmas b-dsigmas))
  (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas) (rest-frame-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas))
  (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
    (close-4vectors a-hats a-dsigmas b-hats b-dsigmas tolerance max-steps))
  (values a-hats a-dsigmas b-hats b-dsigmas))

;;Compare hats and sigmas recovered after evolution with those from before.  The new ones start at a randomly
;;different place in the loop.
(defun compare-hats-sigmas (hats dsigmas hats-1 dsigmas-1)
  (assert (= (length hats) (length dsigmas) (length hats-1) (length dsigmas-1)))
  (let ((length (length hats))
	(best nil)
	best-quality)
    (loop for index below length
	  for this = (abs (- (aref dsigmas 0) (aref dsigmas-1 index)))
	  when (or (null best) (< this best-quality))
	  do (setq best index best-quality this))
    (values (loop for index below length
		  for index-1 = (mod (+ index best) length)
		  maximize (abs (- (aref dsigmas index) (aref dsigmas-1 index-1))))
	    (loop for index below length
		  for index-1 = (mod (+ index best) length)
		  maximize (3vector-distance (aref hats index) (aref hats-1 index-1))))))
				   
;;Find average rotation of a' and b' due to backreaction
;;We're trying to find a rotation matrix R that will rotate each a-hat into the corresponding new a-hat, and the
;;b-hats also.  Of course it cannot be done exactly, but we're trying to find the closest approximation.
;;To find it, we write R as a 9-vector (Rxx Rxy Rxz Ryx ... Rzz).  The x component of R acting on a' is then
;;(a'x a'y a'z 0 0 0 0 0 0) . R, the y component is (0 0 0 a'x a'y a'z 0 0 0) . R, and so on.  We assemble these
;;row 9-vectors into a matrix M with 3 rows coming from each a-hat.  We then want to solve as well as possible the
;;overdetermined set of equations M.R = C where C is a column vector with a-hats and b-hats concatenated.
;;We can do this by singular value decomposition.
;;I think there is a better way to do this.
(defun find-backreaction-rotation (a-hats a-dsigmas new-a-hats b-hats b-dsigmas new-b-hats)
  (mirror-images (assert (= (length a-hats) (length new-a-hats))))
  (let* ((equations (* 3 (+ (length a-hats) (length b-hats)))) ;Total number of elements of the hats
	 (m (make-array (list equations 9) :element-type 'double-float :initial-element 0.0))
	 (c (make-array equations :element-type 'double-float))
	 (r (make-array '(3 3) :element-type 'double-float))
	 (row 0))
    (mirror-images
     (loop for index below (length a-hats)
	   for a = (3vector-scale (aref a-hats index) (aref a-dsigmas index)) ;Weight by dsigma
	   for new-a = (3vector-scale (aref new-a-hats index) (aref a-dsigmas index)) ;Same length
	   do (loop for new-axis below 3 ;coordinate of new a' that we're trying to match
		    do (loop for axis below 3 ;coordinate of old a' contributing
			     do (setf (aref m row (+ (* new-axis 3) axis)) (aref a axis)))
		    do (setf (aref c row) (3vector-component new-a new-axis))
		    do (incf row))))
    (handler-case 
	(progn
	  (multiple-value-bind (u w v) (svdcmp m) ;u is the same array as m.  w is just the diagonal elements.
	    (dotimes (i1 3)
	      (loop for i2 below 3
		    for i from (* 3 i1)	;index into flat R vector
		    do (setf (aref r i1 i2)
			     (loop for j below 9
				   sum (* (/ (aref v i j) (aref w j))
					  (loop for k below equations
						sum (* (aref u k j) (aref c k)))))))))
	  (multiple-value-bind (u w v) (svdcmp r) ;Now coerce it to the closest rotation matrix
	    (declare (ignore w))		    ;by making eigenvalues 1.
	    (dotmm u (transpose v))))
      (svdcmp-not-converging-error () nil)))) ;NIL on failure

(defun rotation-matrix-info (r)
  (multiple-value-bind (eigenvalues eigenvectors) (eigensystem-nonsymmetric r)
    (loop with axis = (make-3vector)
	  for index below 3
	  when (zerop (imagpart (aref eigenvalues index))) ;Find real eigenvalue
	  do
	  (dotimes (i 3)
	    (setf (3vector-component axis i) (realpart (aref eigenvectors i index)))) ;Copy column
	  (return
	   (values (phase (aref eigenvalues (mod (1+ index) 3)))		   ;Amount of rotation
		   axis)))))		;I don't know how to find sign of rotation
	  

(defun find-gravitational-rotation (a-hats a-dsigmas b-hats b-dsigmas)
  (mirror-image-let* ((length (nomirror (total-sigma a-dsigmas)))
		      (n-a (length a-hats))
		      (d-a (make-array n-a))
		      (d-a-rot (make-array n-a))
		      (total-a-rot (make-zero-3vector))
		      (m-a (make-array '(3 3) :element-type 'double-float :initial-element 0.0))
		      a-theta
		      theta)
    (mirror-images			;Compute all segment changes
     (thread-dispatch n-a *threads*
       #'(lambda (iter-a)
	   (with-local-resources
	     (let ((d-a-prime (a-segment-change a-hats a-dsigmas b-hats b-dsigmas iter-a))) ;rate of change of A'
	       (setf (aref d-a iter-a) d-a-prime)
	       (setf (aref d-a-rot iter-a) (3vector- d-a-prime ;Spatial part
						     (3vector-scale (aref a-hats iter-a) (3vector-t d-a-prime))))))))
     (dotimes (i 3)
       (setf (aref m-a i i) length))
     (loop for hat across a-hats
	   for dsigma across a-dsigmas
	   for rot across d-a-rot
	   do (dotimes (i 3)
		(dotimes (j 3)
		  (decf (aref m-a i j) (* (3vector-component hat i) (3vector-component hat j)
					  dsigma))))						   ;Weighted sum
	   do (3vector-incf total-a-rot (3vector-scale (3vector-cross-product hat rot) dsigma)))   ;Total A'_rot
     (setq a-theta (linear-solve m-a total-a-rot)))
    (setq theta (linear-solve (matrix+ m-a m-b) (3vector+ total-a-rot total-b-rot))) ;Best Theta w/ both
    (mirror-image-let (error-a len-a pos-len-a)
      (mirror-images
       (loop for hat across a-hats
	     for dsigma across a-dsigmas
	     for rot across d-a-rot
	     for d-a-prime across d-a
	     ;;total squared length of longitudinal changes
	     for len = (* dsigma (expt (4vector-t d-a-prime) 2))
	     sum len into len-change
	     when (plusp (4vector-t d-a-prime)) sum len into pos-len-change
	     ;;total squared length of spatial changes
	     sum (* dsigma (3vector-squared-length rot)) into space-change
	     ;;weighted average of distance from prediction
	     sum (* dsigma (3vector-squared-length (3vector- rot (3vector-cross-product theta hat)))) into error
	     finally (setq error-a (/ error space-change) ;Fraction of weighted total not predicted by rotation
			   len-a (/ len-change (+ len-change space-change))
			   pos-len-a (/ pos-len-change (+ len-change space-change)))
	     	     (format t "pos-len ~S len ~S space ~S~%" pos-len-change len-change space-change )))
      (values a-theta b-theta theta error-a error-b len-a len-b pos-len-a pos-len-b))))
			 
;;Compare results of different calculation techniques
(defun compare-hats-dsigmas (a-hats-1 a-dsigmas-1 b-hats-1 b-dsigmas-1 a-hats-2 a-dsigmas-2 b-hats-2 b-dsigmas-2)
  (compare-hats-dsigmas-1 (concatenate 'vector a-hats-1 b-hats-1) (concatenate 'vector a-dsigmas-1 b-dsigmas-1)
			  (concatenate 'vector a-hats-2 b-hats-2) (concatenate 'vector a-dsigmas-2 b-dsigmas-2)))

(defun compare-hats-dsigmas-1 (hats-1 dsigmas-1 hats-2 dsigmas-2)
  (let ((angle (loop for hat-1 across hats-1
		     for hat-2 across hats-2
		     maximize (spherical-angle hat-1 hat-2))))
    (loop for dsigma-1 across dsigmas-1
	  for dsigma-2 across dsigmas-2
	  maximize (abs (- dsigma-1 dsigma-2)) into absolute
	  maximize (/ (abs (- dsigma-1 dsigma-2)) dsigma-1) into relative
	  finally (format t "Maximum angle ~S, maximum segment length difference ~S abs, ~S rel"
			  angle absolute relative))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; CREATE HATS AND SIGMAS FOR DIFFERENT LOOP CLASSES
;;;; BACKREACTION ON DIFFERENT LOOP CLASSES

;; Does backreaction on a class of loop specified by CLASS-MAKER. OUT is the output file, while the arguments
;; MAKER-ARGUMENTS is a list of arguments that is to be given to the MAKE-*-HATS. To change any of the other keywords
;; of DO-BACKREACTION, specify them after the OUT argument.
(defun backreaction-on-class (class-maker maker-arguments out &rest backreaction-arguments)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (apply class-maker maker-arguments)
		       (apply #'do-backreaction a-hats a-dsigmas b-hats b-dsigmas :out out backreaction-arguments)))

;; Does backreaction on the a and b hats and sigmas loaded from a .DAT file of the sort produced by WRITE-HATS-DSIGMAS.
;; To change any of the keywords of DO-BACKREACTION, specify them after the FILE-NAME argument.
(defun backreaction-on-data (file-name &rest backreaction-arguments)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas file-name)
   (apply #'do-backreaction a-hats a-dsigmas b-hats b-dsigmas backreaction-arguments)))

;; Does backreaction on a Burden loop. See MAKE-BURDEN-HATS for the meanings of the arguments.
(defun do-burden-backreaction (m n phi &key (out (format nil "b~D~D" m n)) (L (* 2 pi))
				 (points 32) (NGmu 0.0001) (iterations 100) (print-at 20))
  (backreaction-on-class #'make-burden-hats (list m n phi :L L :points points)
			 out :NGmu NGmu :iterations iterations :print-at print-at))

;; Does backreaction on a broken Burden/smeared GV loop. See MAKE-BROKEN-BURDEN-HATS for the meanings of the arguments.
(defun do-broken-burden-backreaction (m n theta &key (out (format nil "bb~D~D" m n)) (L (* 2 pi)) (psi (/ pi 2))
					(psi-a psi) (psi-b psi) (points 32) (NGmu 0.0001) (iterations 100)
					(print-at 20) (track-kinks 4))
  (backreaction-on-class #'make-broken-burden-hats (list m n theta :L L :psi psi :psi-a psi-a :psi-b psi-b :points points)
			 out :NGmu NGmu :iterations iterations :print-at print-at :track-kinks track-kinks))

;; Does backreaction on an ACO loop, either smooth or polygonal. See MAKE-ACO-HATS for the meanings of the arguments.
(defun do-aco-backreaction (m n &key (out (format nil "aco~D~D" m n)) (L (* 2 pi)) (polygon nil) (print-at 20)
			      (points 32) (NGmu 0.0001) (iterations 100) (track-kinks 4))
  (backreaction-on-class #'make-aco-hats (list m n :L L :polygon polygon :points points)
			 out :NGmu NGmu :iterations iterations :print-at print-at :track-kinks track-kinks))

;; Does backreaction on a KT loop. See MAKE-KT-HATS for the meanings of the arguments.
(defun do-kt-backreaction (alpha phi &key (out (format nil "kt~3$" alpha)) (L (* 2 pi)) (print-at 20)
				 (points 32) (NGmu 0.0001) (iterations 100))
  (backreaction-on-class #'make-kt-hats (list alpha phi :L L :points points)
			 out :NGmu NGmu :iterations iterations :print-at print-at))

;; Does backreaction on a broken KT loop. See MAKE-BROKEN-KT-HATS for the meanings of the arguments.
(defun do-broken-kt-backreaction (alpha phi &key (psi (/ pi 2)) (out (format nil "bkt~3$" alpha)) (L (* 2 pi))
					(points 32) (NGmu 0.0001) (iterations 100) (print-at 20) (track-kinks 4))
  (backreaction-on-class #'make-broken-kt-hats (list alpha phi :psi psi :L L :points points)
			 out :NGmu NGmu :iterations iterations :print-at print-at :out out :track-kinks track-kinks))

;; Does backreaction on a twice-broken KT loop. See MAKE-2BROKEN-KT-HATS for the meanings of the arguments.
(defun do-2broken-kt-backreaction (alpha phi &key (out (format nil "bb~3$" alpha)) (L (* 2 pi)) (psi (/ pi 2))
					(psi-a psi) (psi-b psi) (points 32) (NGmu 0.0001) (iterations 100)
					(print-at 20) (track-kinks 4))
  (backreaction-on-class #'make-2broken-kt-hats (list alpha phi :L L :psi psi :psi-a psi-a :psi-b psi-b :points points)
			 out :NGmu NGmu :iterations iterations :print-at print-at :track-kinks track-kinks))


;; Creates :POINTS hats of a loop worldsheet function of length :L according to some HAT-FUNCTION, as well as an
;; evenly-spaced sigma list.
;; set :START to begin at some location on the worldsheet function other than what HAT-FUNCTION thinks "zero" is
;; set :RANGE to only create part of the worldsheet function. The resulting HATS and SIGMAS should only be used in
;;     CREATE-WORLDSHEET with :OPEN-LOOP T.
;; The HAT-FUNCTION must take, in order: the X (where 0<X<POINTS) of the hat to be created,
;;                                       the POINTS in total in the worldsheet function, and
;;                                       a list PARAMETERS of keywords.
(defun make-hats-dsigmas (hat-function &key (L (* 2 pi)) (points 32) (start 0) (range points))
  (let ((hats (make-array range))
	(dsigmas (make-array range :element-type 'double-float)))
    (loop for i below range
	  do (setf (aref hats i) (funcall hat-function (/ (mod (+ i start (/ 2)) points) points));+1/2 avoids intersexns
		   (aref dsigmas i) (/ L points)))
    (values hats dsigmas)))

;; Provides a function that creates 3vectors for use in generating KT- and Burden-type hats lists.
;; set :M to set the winding number.
;; set :PHI to rotate out of the xy plane.
;; set :ALPHA to introduce wigglyness, per the alpha parameter in the KT model.
(defun circle-hat-function (&key (m 1) (phi 0) (alpha 0))
  (lambda (x) (make-3vector (+ (* (- 1 alpha) (cos (* 2 pi m x))) (* alpha (cos (* 6 pi m x))))
			    (- (* (cos phi) (+ (* (- 1 alpha) (sin (* 2 pi m x))) (* alpha (sin (* 6 pi m x)))))
			       (* (sin phi) 2 (sqrt (* alpha (- 1 alpha))) (sin (* 2 pi m x))))
			    (+ (* (cos phi) 2 (sqrt (* alpha (- 1 alpha))) (sin (* 2 pi m x)))
			       (* (sin phi) (+ (* (- 1 alpha) (sin (* 2 pi m x))) (* alpha (sin (* 6 pi m x)))))))))

(defun make-kt-hats-1 (alpha phi &key (L (* 2 pi)) (points 32) (points-a points) (points-b points)
			     (start-a 0) (start-b 0) (range points) (range-a range) (range-b range))
  (multiple-value-bind (a-hats a-dsigmas) (make-hats-dsigmas (circle-hat-function :alpha alpha)
							     :L L :points points-a :start start-a :range range-a)
    (multiple-value-bind (b-hats b-dsigmas) (make-hats-dsigmas (circle-hat-function :phi phi)
							     :L L :points points-b :start start-b :range range-b)
      (values a-hats a-dsigmas b-hats b-dsigmas))))

;; Creates hats and dsigmas for the Burden loop. M is the number of times a' goes around, N is the number of times b'
;; goes around, and PHI is the angle between a' and b'
;; Note that here, as in the broken hats case, :POINTS is the number of points per cycle in a or b.
(defun make-burden-hats (m n phi &key (L (* 2 pi)) (points 32))
  (flet
   ((check-ratio (m name)
		 (let ((ratio (/ points (* 2 m))))
		   (when (and (oddp (numerator ratio)) (oddp (denominator ratio)))
		     (warn "POINTS/2~A = ~S, so there will be diamonds moving at the speed of light" name ratio)))))
   (check-ratio m :m) (check-ratio n :n))
  (let ((a-dsigmas (make-array (* m points) :element-type 'double-float))
	(b-dsigmas (make-array (* n points) :element-type 'double-float))
        (a-hats (make-array (* m points)))
        (b-hats (make-array (* n points))))
    (dotimes (i (* m points))
      (setf (aref a-hats i) ((lambda (x) (make-3vector (cos x) (sin x) 0.0)) (/ (+ (* 2 pi i) pi) points))
	    (aref a-dsigmas i) (/ L (* m points))))
    (dotimes (j (* n points))
      (setf (aref b-hats j) ((lambda (x) (make-3vector (* (cos x) (cos phi)) (sin x) (* (cos x) (sin phi))))
			     (/ (+ (* 2 pi j) pi) points))
            (aref b-dsigmas j) (/ L (* n points))))
    (values a-hats a-dsigmas b-hats b-dsigmas)))

;; Creates hats and sigmas for the worldsheet functions of the "broken burden loop"/"smeared GV loop". This loop is the
;; Burden loop, but with jumps of angle psi-a and psi-b in the a and b functions centered on the cusps. That is, we
;; insert kinks by hand to prevent cusps from appearing. This is equivalent to taking the planar GV loop and smearing
;; its a' and b' in the out-of-plane direction on the KT sphere.
;; m, n, and theta are the usual parameters passed to a Burden loop. psi is the total angular amount to remove from a'
;; or b' per cusp per cycle; psi-a and psi-b allow for different amounts to be removed from a' or b' if so desired.
;; Note that here, :POINTS is the number of points per winding of a or b! E.G., a has M*POINTS segments total.
(defun make-broken-burden-hats (m n theta &key (L (* 2 pi)) (psi (/ pi 2)) (psi-a psi) (psi-b psi) (points 32))
  (let ((a-dsigmas (make-array (* m points) :element-type 'double-float))
	(b-dsigmas (make-array (* n points) :element-type 'double-float))
        (a-hats (make-array (* m points)))
        (b-hats (make-array (* n points)))
	(per-side (/ points 2))
        (phi (/ theta 2))
        (a-polar 0)
        (b-polar 0))
    (unless (integerp per-side) ;if this is true, the loop would be asymmetric
      (error "POINTS/2 must be integer-valued for MAKE-BROKEN-BURDEN-HATS to function correctly."))
    (dotimes (i (* m points))
      (setf a-polar (- (* (/ (rem i per-side) (1- per-side)) (- pi psi-a)) (- (/ pi 2) (/ psi-a 2)))
            (aref a-hats i) (3vector-scale (make-3vector (* (cos phi) (cos a-polar))
							 (- (* (sin phi) (cos a-polar))) (sin a-polar))
					   (if (evenp (floor i per-side)) 1 -1))
            (aref a-dsigmas i) (/ L (* m points))))
    (dotimes (j (* n points))
      (setf b-polar (- (* (/ (rem j per-side) (1- per-side)) (- pi psi-b)) (- (/ pi 2) (/ psi-b 2)))
            (aref b-hats j) (3vector-scale (make-3vector (* (cos phi) (cos b-polar))
							 (* (sin phi) (cos b-polar)) (sin b-polar))
					   (if (evenp (floor j per-side)) 1 -1))
	    (aref b-dsigmas j) (/ L (* n points))))
    (values a-hats a-dsigmas b-hats b-dsigmas)))

;; Creates hats and dsigmas for the smooth and polygonal ACO loops. The parameter M tells how many times the polar
;; worldsheet function goes back and forth, while N tells how many times to go around for the equatorial worldsheet
;; function.
;; If the key :POLYGON is set, the hats created are for a polygonal equatorial worldsheet function which goes around N
;; times. Thus, POINTS/(POLYGON * N) must be an integer, or else we'll have an uneven distribution of points!
(defun make-aco-hats (m n &key (L (* 2 pi)) (polygon nil) (points 32))
  (let ((sigmas (make-array points :element-type 'double-float))
        (a-hats (make-array points))
        (b-hats (make-array points)))
    (if polygon (unless (integerp (/ points (* n polygon)))
		  (error "The ratio POINTS/(N*POLYGON) must be an integer for MAKE-ACO-HATS to function correctly.")))
    (unless (integerp (/ points (* 2 m)))
      (error "The ratio POINTS/2M must be an integer for MAKE-ACO-HATS to function correctly."))
    (dotimes (i points)
      (setf (aref a-hats i) (make-3vector 0.0 0.0 (if (evenp (floor (/ (* 2 i m) points))) 1.0 -1.0))
            (aref b-hats i) (if polygon
				(make-3vector (cos (/ (* 2 pi (truncate (* i polygon n) points)) polygon))
					      (sin (/ (* 2 pi (truncate (* i polygon n) points)) polygon))
					      0.0)
                              (make-3vector (cos (/ (* 2 pi i n) points)) (sin (/ (* 2 pi i n) points)) 0.0))
            (aref sigmas i) (/ L points)))
    (values a-hats sigmas b-hats sigmas)))

;; Creates hats and sigmas for the Kibble-Turok loop. ALPHA sets the perturbation of a' away from circular, while PHI
;; sets the angle between a' and b'.
(defun make-kt-hats (alpha phi &key (L (* 2 pi)) (points 32))
  (let ((sigmas (make-array points :element-type 'double-float))
	(a-hats (make-array points))
	(b-hats (make-array points)))
    (dotimes (i points)
      (setf (aref a-hats i) ((lambda (x) (make-3vector (+ (* (- 1 alpha) (cos x)) (* alpha (cos (* 3 x))))
						       (+ (* (- 1 alpha) (sin x)) (* alpha (sin (* 3 x))))
						       (* 2 (sqrt (* alpha (- 1 alpha))) (sin x))))
			     (/ (+ (* 2 pi i) pi) points))
	    (aref b-hats i) ((lambda (x) (make-3vector (cos x) (* (cos phi) (sin x)) (* (sin phi) (sin x))))
			     (/ (+ (* 2 pi i) pi) points))
            (aref sigmas i) (/ L points)))
    (values a-hats sigmas b-hats sigmas)))

;; Creates hats and sigmas for the broken Kibble-Turok loop. ALPHA and PHI are as in MAKE-KT-HATS, while the keyword
;; :PSI sets the angular size of the wedge taken out from b' around the cusps.
(defun make-broken-kt-hats (alpha phi &key (L (* 2 pi)) (points 32) (psi (/ pi 2)))
  (let ((sigmas (make-array points :element-type 'double-float))
        (a-hats (make-array points))
        (b-hats (make-array points))
        (decider (/ points 2)))
    (unless (integerp decider) (error "POINTS/2 must be an integer for MAKE-BROKEN-KT-HATS to function correctly"))
    (dotimes (i points)
      (setf (aref a-hats i) ((lambda (x) (make-3vector (+ (* (- 1 alpha) (cos x)) (* alpha (cos (* 3 x))))
						       (+ (* (- 1 alpha) (sin x)) (* alpha (sin (* 3 x))))
						       (* 2 (sqrt (* alpha (- 1 alpha))) (sin x))))
			     (/ (+ (* 2 pi i) pi) points))
	    (aref b-hats i) ((lambda (x) (make-3vector (cos x) (* (cos phi) (sin x)) (* (sin phi) (sin x))))
			     (if (< i decider) (+ (/ psi 2) (/ (* i (- pi psi)) (1- decider)))
			       (+ pi (/ psi 2) (/ (* (- i decider) (- pi psi)) (1- decider)))))
	    (aref sigmas i) (/ L points)))
    (values a-hats sigmas b-hats sigmas)))

;; replace the NEAREST segments in B' below the kinks with SPLIT segments of equal dsigma, equally spread out across
;; the arc from the original segment's position to just before the position of the segment prior to it.
;; note that the final number of points in B' is POINTS+2*NEAREST*(SPLIT-1)
(defun make-broken-kt-hats-half-dense (alpha phi nearest split &key (L (* 2 pi)) (points 32) (psi (/ pi 2)))
  (mirror-image-let
   (a-hats a-dsigmas)
   (multiple-value-setq
    (a-hats a-dsigmas b-hats b-dsigmas)
    (make-broken-kt-hats alpha phi :L L :points points :psi psi))
   (setf b-hats (coerce b-hats 'list)
	 b-dsigmas (coerce b-dsigmas 'list))
   (let* ((split-dsigmas (make-list (* nearest split) :initial-element (/ L (* points split))))
	  (dangle (/ (* 2 (- pi psi)) (* (- points 2) split)))
	  (split-angles (loop for i from 1 to (* nearest split) collect (* (- i (* split nearest)) dangle)))
	  new-b-hats new-b-dsigmas)
     (setf new-b-dsigmas (nconc (subseq b-dsigmas 0 (- (/ points 2) nearest))
				(copy-seq split-dsigmas)
				(subseq b-dsigmas (/ points 2) (- points nearest))
				(copy-seq split-dsigmas))
	   new-b-hats (nconc (subseq b-hats 0 (- (/ points 2) nearest))
			     (loop for da in split-angles
				   for x = (+ pi da (/ (- psi) 2))
				   collect (make-3vector (cos x) (* (cos phi) (sin x)) (* (sin phi) (sin x))))
			     (subseq b-hats (/ points 2) (- points nearest))
			     (loop for da in split-angles
				   for x = (+ (* 2 pi) da (/ (- psi) 2))
				   collect (make-3vector (cos x) (* (cos phi) (sin x)) (* (sin phi) (sin x))))))
     (values a-hats a-dsigmas (coerce new-b-hats 'vector) (coerce new-b-dsigmas '(vector double-float))))))

;; Creates hats and sigmas for the broken Kibble-Turok loop with dense kinks. All parameters are as in MAKE-KT-HATS,
;; and the additional parameters set how the dense kinks are formed. The NEAREST segments on both sides of the kinks
;; are divided into SPLIT segments each, with evenly-distributed dsigmas.
;;
;; N.B.! The final length of the b hats and sigmas will be :POINTS+4*NEAREST*(SPLIT-1)
(defun make-broken-kt-hats-dense (alpha phi nearest split &key (L (* 2 pi)) (points 32) (psi (/ pi 2)))
  (let*((a-dsigmas (make-array points :element-type 'double-float))
	(a-hats (make-array points))
	(new-points (+ points (* 4 nearest (1- split))))
	(b-dsigmas (make-array new-points :element-type 'double-float))
	(b-hats (make-array new-points))
	(dense-dsigma (/ L (* points split)))
	(sparse-dsigma (/ L points))
	(dsigma 0)
	(dangle (/ (- pi psi) (1- (/ points 2))))
	(next-angle 0) (j 0) (side 0))
    (unless (< nearest (/ points 4)) (error "NEAREST is too large. It cannot be greater than POINTS/4"))
    (unless (integerp (/ points 2)) (error "POINTS/2 must be an integer for MAKE-BROKEN-KT-HATS to function correctly"))
    (dotimes (i points)
      (setf (aref a-hats i) ((lambda (x) (make-3vector (+ (* (- 1 alpha) (cos x)) (* alpha (cos (* 3 x))))
						       (+ (* (- 1 alpha) (sin x)) (* alpha (sin (* 3 x))))
						       (* 2 (sqrt (* alpha (- 1 alpha))) (sin x))))
			     (/ (+ (* 2 pi i) pi) points))
	    (aref a-dsigmas i) (/ L points)))
    (dotimes (i new-points)
      (multiple-value-setq (side j) (truncate i (/ new-points 2)))
      (cond ((= j 0)
	     (setf next-angle (+ (/ psi 2) (* side pi))
		   dsigma dense-dsigma))
	    ((< j (* nearest split))
	     (setf dsigma dense-dsigma)
	     (when (integerp (/ j split))
	       (setf next-angle (+ next-angle dangle))))
	    ((< j (+ (/ points 2) (* nearest (- split 2))))
	     (setf next-angle (+ next-angle dangle)
		   dsigma sparse-dsigma))
	    ((= j (+ (/ points 2) (* nearest (- split 2))))
	     (setf next-angle (+ next-angle dangle)
		   dsigma dense-dsigma))
	    ((< j (+ (/ points 2) (* 2 nearest (1- split))))
	     (setf dsigma dense-dsigma)
	     (when (integerp (/ (- j (+ (/ points 2) (* nearest (- split 2)))) split))
	       (setf next-angle (+ next-angle dangle))))
	    (t
	     (setf next-angle (+ next-angle dangle)
		   dsigma dense-dsigma)))
      (setf (aref b-hats i) ((lambda (x) (make-3vector (cos x) (* (cos phi) (sin x)) (* (sin phi) (sin x)))) next-angle)
	    (aref b-dsigmas i) dsigma))
    (values a-hats a-dsigmas b-hats b-dsigmas)))

;; Creates hats and sigmas for the doubly-broken Kibble-Turok loop. ALPHA and PHI are as in MAKE-KT-HATS, while the
;; keywords :PSI-A and :PSI-B set the angular size of the wedges taken from a' and b', respectively, around the cusps.
(defun make-2broken-kt-hats (alpha phi &key (L (* 2 pi)) (points 32) (psi (/ pi 2)) (psi-a psi) (psi-b psi))
  (let ((dsigmas (make-array points :element-type 'double-float))
        (a-hats (make-array points))
        (b-hats (make-array points))
        (decider (/ points 2)))
        (unless (integerp decider) (error "POINTS/2 must be an integer for MAKE-2BROKEN-KT-HATS to function correctly"))
        (dotimes (i points)
          (setf (aref a-hats i) ((lambda (x) (make-3vector (+ (* (- 1 alpha) (cos x)) (* alpha (cos (* 3 x))))
                                                           (+ (* (- 1 alpha) (sin x)) (* alpha (sin (* 3 x))))
                                                           (* 2 (sqrt (* alpha (- 1 alpha))) (sin x))))
                                 (if (< i decider) (+ (/ psi-a 2) (/ (* i (- pi psi-a)) (1- decider)))
                                                   (+ pi (/ psi-a 2) (/ (* (- i decider) (- pi psi-a)) (1- decider)))))
                (aref b-hats i) ((lambda (x) (make-3vector (cos x) (* (cos phi) (sin x)) (* (sin phi) (sin x))))
                                 (if (< i decider) (+ (/ psi-b 2) (/ (* i (- pi psi-b)) (1- decider)))
                                                   (+ pi (/ psi-b 2) (/ (* (- i decider) (- pi psi-b)) (1- decider)))))
                (aref dsigmas i) (/ L points)))
    (values a-hats dsigmas b-hats dsigmas)))

;; Makes broken KT hats, but the kink jumps in the b' no longer avoid the crossing. Thus, there are two cusps and two
;; kinks. All arguments are otherwise the same as in MAKE-BROKEN-KT-HATS
(defun make-mixed-kt-hats (alpha phi &key (L (* 2 pi)) (points 32) (psi (/ pi 2)))
  (let ((dsigmas (make-array points :element-type 'double-float))
        (a-hats (make-array points))
        (b-hats (make-array points))
        (decider (/ points 2)))
        (unless (integerp decider) (error "POINTS/2 must be an integer for MAKE-MIXED-KT-HATS to function correctly"))
        (dotimes (i points)
          (setf (aref a-hats i) ((lambda (x) (make-3vector (+ (* (- 1 alpha) (cos x)) (* alpha (cos (* 3 x))))
                                                           (+ (* (- 1 alpha) (sin x)) (* alpha (sin (* 3 x))))
                                                           (* 2 (sqrt (* alpha (- 1 alpha))) (sin x))))
                                 (/ (+ (* 2 pi i) pi) points))
                (aref b-hats i) ((lambda (x) (make-3vector (- (sin x)) (* (cos phi) (cos x)) (* (sin phi) (cos x))))
                                 (if (< i decider) (+ (/ psi 2) (/ (* i (- pi psi)) (1- decider)))
                                                   (+ pi (/ psi 2) (/ (* (- i decider) (- pi psi)) (1- decider)))))
                (aref dsigmas i) (/ L points)))
    (values a-hats dsigmas b-hats dsigmas)))

;; Makes twice-broken burden hats, but the kink jumps in the tangent vectors no longer avoid the crossing. Thus, there
;; are both cusps and kinks. All arguments are otherwise the same as in MAKE-2BROKEN-BURDEN-HATS
(defun make-mixed-burden-hats (m n theta &key (L (* 2 pi)) (psi (/ pi 2)) (psi-a psi) (psi-b psi) (points 32))
  (let ((a-dsigmas (make-array (* m points) :element-type 'double-float))
	(b-dsigmas (make-array (* n points) :element-type 'double-float))
        (a-hats (make-array (* m points)))
        (b-hats (make-array (* n points)))
	(per-side (/ points 2))
        (phi (/ theta 2))
        (a-polar 0)
        (b-polar 0))
    (unless (integerp per-side) ;if this is true, the loop would be asymmetric
      (error "POINTS/2 must be integer-valued for MAKE-BROKEN-BURDEN-HATS to function correctly."))
    (dotimes (i (* m points))
      (setf a-polar (- (* (/ (rem i per-side) (1- per-side)) (- pi psi-a)) (- (/ pi 2) (/ psi-a 2)))
            (aref a-hats i) (3vector-scale (make-3vector (* (cos phi) (sin a-polar))
							 (- (* (sin phi) (sin a-polar))) (- (cos a-polar)))
					   (if (evenp (floor i per-side)) 1 -1))
            (aref a-dsigmas i) (/ L (* m points))))
    (dotimes (j (* n points))
      (setf b-polar (- (* (/ (rem j per-side) (1- per-side)) (- pi psi-b)) (- (/ pi 2) (/ psi-b 2)))
            (aref b-hats j) (3vector-scale (make-3vector (* (cos phi) (sin b-polar))
							 (* (sin phi) (sin b-polar)) (- (cos b-polar)))
					   (if (evenp (floor j per-side)) 1 -1))
	    (aref b-dsigmas j) (/ L (* n points))))
    (values a-hats a-dsigmas b-hats b-dsigmas)))

;;Return hats and sigmas for a and b for the Garfinkle-Vachaspati loop.  Theta is the angle between a' and b',
;;pi/2 for the square loop.  Count is the is the number of kinks in a or b.  In the simplest case it is 2.
(defun make-gv-hats (theta count &key (L 4.0) perturbation)
  (assert (evenp count))
  (setq L (double-float L))
  (let* ((phi (/ theta 2))
     (h (sin phi))
     (w (cos phi))
     (a-hat (make-3vector w (- h) 0.0))
     (b-hat (make-3vector w h 0.0)))
    (mirror-image-let ((a-hats (make-array count))
		       (a-dsigmas (make-array count :element-type 'double-float)))    
      (mirror-images
       (loop for index below count
         for this-hat = (if (< index (floor count 2)) a-hat (3vector- a-hat))
         when perturbation
         do (setq this-hat (3vector-normalize
			    (3vector+ this-hat (make-3vector 0.0 (random perturbation) 0.0))))
         do (setf (aref a-dsigmas index) (/ l count)
		  (aref a-hats index) this-hat))
       (when perturbation
	 ;;Because a' consists of vectors pointing almost in the exact same direction,
	 ;;perpendicular changes don't change the length of any vector, making it difficult
	 ;;for close-hats to do its job, so we need looser tolerance
	 (setq a-hats (close-hats a-hats :tolerance 1e-12)))) ;Account for changes so loop is closed
      (values a-hats a-dsigmas b-hats b-dsigmas))))

;;Make a helical standing wave wrapping one of the periodicity directions, with n turns of the helix.
;;Slope is the slope of the a or b function relative to  the axis of the helix.
(defun make-helical-standing-wave-hats (turns slope points)
  (let* ((periodicity (first (periodicity-vectors))) ;Arbitrary periodicity direction
	 (spatial-length (3vector-length periodicity))
	 (z (3vector-scale periodicity (/ 1 spatial-length))))
    (multiple-value-bind (x y) (two-perpendicular-vectors z) ;arbitrary basis for perpendicular plane
      (let* ((apz (sqrt (- 1.0 (expt slope 2))))	     ;Constant a'_z
	     (max-sigma (/ spatial-length apz))
	     (dsigma (/ max-sigma points))
	     (k (/ (* 2 pi turns) max-sigma))) ;k sigma = 0...2piN
	(flet ((ab (sigma sigma-sign rotation-direction) ;make A' or B'
		 (setq sigma (* sigma sigma-sign))	 ;The argument to A is really -sigma
		 (3vector+ (3vector-scale x (* slope sigma-sign (cos (* k sigma))))
			   (3vector-scale y (* slope sigma-sign rotation-direction (sin (* k sigma))))
			   (3vector-scale z (* apz sigma-sign)))))
	  (mirror-image-let ((a-hats (make-array points))
			     (a-dsigmas (make-array points :element-type 'double-float :initial-element dsigma)))
	    (loop for point below points
		  for sigma = (* dsigma point)
		  do (setf (aref a-hats point) (ab sigma -1.0 1.0)
			   (aref b-hats point) (ab sigma 1.0 1.0)))
	    (values a-hats a-dsigmas b-hats b-dsigmas)))))))

;;;; FILE INPUT/OUTPUT
;;;; PRINTING/ANALYZING RESULTS
;;;; PRODUCING PLOTS

;; Given a string of data delimited by :SEPARATOR, converts the string into an array of data.
(defun line-to-list (line &key (separator #\space))
  (unless (stringp line) (return-from line-to-list nil))
  (loop for i = 0 then (1+ j)
        for j = (position separator line :start i)
        collect (subseq line i j)
        while j))

;; Given a string which we know is secretly a float-valued bit of data, turns it into a float.
(defun parse-double (input)
  (with-input-from-string (in input) (read in)))

;; Writes the a and b hats and dsigmas to a file
(defun write-hats-dsigmas (a-hats a-dsigmas b-hats b-dsigmas file-name)
  (mirror-images (check-dsigmas a-dsigmas))
  (with-open-file (s file-name :element-type '(unsigned-byte 64) :direction :output
		     :if-does-not-exist :create :if-exists :supersede)
    (mirror-images
     (write-double (length a-hats) s)
     (loop for i across a-hats do (write-4vector i s))
     (write-double (length a-dsigmas) s)
     (loop for i across a-dsigmas do (write-double i s)))))

;; Reads the a and b hats and dsigmas from a file whose format we know to be the one of WRITE-HATS-DSIGMAS
(defun read-hats-dsigmas (file-name)
  (mirror-image-let
   ((a-hats (make-array 0 :fill-pointer 0 :adjustable t))
    (a-dsigmas (make-array 0 :fill-pointer 0 :adjustable t :element-type 'double-float)))
   (with-open-file
    (s file-name :element-type '(unsigned-byte 64) :direction :input :if-does-not-exist :error)
    (let ((x 0))
      (mirror-images
       (setf x (truncate (read-double s)))
       (dotimes (i x) (vector-push-extend (read-4vector s) a-hats))
       (setf x (truncate (read-double s)))
       (dotimes (i x) (vector-push-extend (read-double s) a-dsigmas)))))
   (mirror-images (setf a-hats (copy-seq a-hats)
			a-dsigmas (copy-seq a-dsigmas)))
   (values a-hats a-dsigmas b-hats b-dsigmas)))

;; Reads the Gamma, L, total J, and (possibly) kinks from a file whose format we know to be derived from DO-BACKREACTION
(defun read-gljk (file-name)
  (let ((gammas (list))
	(Ls (list))
	(Js (list))
	(kinks (list)))
    (with-open-file
     (s file-name :direction :input :if-does-not-exist nil)
     (loop for x = (map 'list #'parse-double (line-to-list (read-line s nil nil)))
	   while x
	   do (push (nth 0 x) gammas)
	   do (push (nth 1 x) Ls)
	   do (push (nth 2 x) Js)
	   do (if (< 3 (length x)) (push (subseq x 3) kinks))))
    (setf gammas (cdr (reverse gammas))
	  Js (reverse Js)
	  Ls (reverse Ls))
    (if kinks (setf kinks (reverse kinks)))
    (values gammas Ls Js kinks)))

;;Read dump and distribute loops as hats and dsigmas.  We don't go to the rest frame, but we do close,
;;because that is the default in get-a-data-dsigmas.
;;NIL for input-directory means use existing data instead of reading
(defun read-distribute-loops (input-directory output-directory
					      &key (dump-step 0) (read-verbose nil) (smallest-xi 0.01)
					      ;;Number of this run for data directory.  We assume that the output
					      ;;is going into something like "backreaction/17"
					      (run-number (car (last (pathname-directory
								      (format nil "~A/" output-directory))))))
  (when input-directory			;input-directory = NIL means don't read
    (setq input-directory (merge-pathnames input-directory batch-root-directory))
    (format t "Reading dump from ~A" input-directory)
    (read-dumps input-directory :step dump-step :verbose read-verbose))
  (setq output-directory (merge-pathnames output-directory batch-root-directory))
  (ensure-directories-exist (format nil "~A/" output-directory))
  (with-group-write-access
   (with-open-file (list-file (format nil "~A/loop-list.text" output-directory) :direction :output)
     (loop for (diamond xi) in (loops-by-xi) ;Go through loops in decreasing order of xi
	   for count from 0		     ;This was 1 in an older version
	   ;;The loop we want are the ones preserved because of their xi, but we also have loops that were made recently
	   ;;and have not been removed yet, including those that have not oscillated once, whose xi is 0.0 meaning unknown
	   until (< xi smallest-xi)
	   when (loop-ok-for-backreaction diamond)
	   do
	   (format t " ~D" count) (force-output)
	   (format list-file "~D	~D	~F	~F~%"
		   run-number count (4vector-t (tag-created-position (diamond-tag diamond))) xi)
	   (multiple-value-call #'write-hats-dsigmas (get-ab-data-dsigmas diamond)
				(format nil "~A/loop-hats-~D.dat" output-directory count)
				)))))

;;We need to know that this loop is non-self-intersecting and has not already intersected with something
;;At the first engulfment, we determine xi and set countup = 1.  This countup then propagates around the loop
;;and when it engulfs again, we know the loop is non-self-intersecting.  If it rejoins something, then 
;;places there will have countup = 0.  So if all countups are 1 or more it is a good loop.
(defun loop-ok-for-backreaction (diamond)
  (loop	with start-diamond = diamond
	do (setq diamond (diamond-e diamond))
	when (eq diamond start-diamond) return t ;successfully looped, so OK
	always (plusp (diamond-countup diamond)))) ;NIL if countup=0
  
(defun diagnose-loop-tags (diamond loop-number)
  (loop with data = (make-hash-table)
	with start-diamond = diamond
	for tag-table = (gethash (diamond-tag diamond) data nil)
	unless tag-table do (setf (gethash (diamond-tag diamond) data) (setq tag-table (make-hash-table)))
	do (incf (gethash (diamond-countup diamond) tag-table 0))
	do (setq diamond (diamond-e diamond))
	until (eq diamond start-diamond)
	finally (cond ((= (hash-table-count data) 1) ;Normal case: there is only one tag.
		       )
		      (t
		       (format t "~&Loop ~D. There are ~D tags~%" loop-number (hash-table-count data))
		       (cond ((loop for tag being the hash-keys of data always (zerop (tag-xi tag)))
			      (format t "All xi = 0, so this is a loop that has never engulfed"))
			     (t
			      (maphash #'(lambda (tag tag-table)
					   (format t "Created ~F, xi = ~F countups = ~S~%"
						   (4vector-t (tag-created-position tag)) (tag-xi tag)
						   (loop for countup being the hash-keys of tag-table
							 using (hash-value n)
							 collect (list countup n))))
					
			       data)))))))

(defun test-loop-corpus ()
  (loop for (diamond xi) in (loops-by-xi) ;Go through loops in decreasing order of xi
	for count from 1
	do (diagnose-loop-tags diamond count)))

;;Loops in decreasing order of their xi.  Returns a list of (DIAMOND xi)
(defun loops-by-xi ()
  (let ((result nil))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore)))
     #'(lambda (d)
	 (push (list d (tag-xi (diamond-tag d)))
 result)))
    (sort result #'> :key #'second)))

;; Calls a plotting function to operate on .DAT files produced by WRITE-HATS-DSIGMAS, and so the file names run in
;; increments of STEP from MIN-STEP up to MAX-STEP. Any additional arguments to the plot command can be passed in
;; PLOTTER-ARGUMENTS.
(defun plot-sequence (plotter file-name step min-step max-step &rest plotter-arguments)
  (unless (integerp (/ (- max-step 0) step)) (error "The ratio MAX-STEP/STEP must be an integer"))
  (loop for i from min-step to max-step by step
	do (apply plotter (format nil "~A-~4,'0D" file-name i) plotter-arguments)))

;; Plots the a and b hats and sigmas taken from a single .DAT file produced by WRITE-HATS-DSIGMAS.
(defun plot-hats-dsigmas (file-name)
  (multiple-value-bind
   (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (mirror-image-let
    ((a-J (compute-hats-angular-momentum a-hats (dsigma-to-sigma a-dsigmas)))
     (a-max (* 1.05 (apply 'max (loop for i across a-dsigmas collect i)))))
    (multiplot
     (:prelude (format nil "set terminal postscript enhanced color~%set output \"~A.eps\"~%" file-name)
	       :columns 2 :rows 2
	       :command-text (format nil "title 'Iteration ~D'"
				     (parse-integer (subseq file-name (- (length file-name) 4)))))
     (plot-tangent-vectors
      (list (list a-hats "a" "linespoints lw 1 lc rgb '#cf0f03'")
	    (list b-hats "b" "linespoints lw 1 lc rgb '#11cc62'"))
      :title "Unit sphere, A' equatorial" :axes (list a-J b-J))
     (plot-tangent-vectors
      (list (list a-hats "a" "linespoints lw 1 lc rgb '#cf0f03'")
	    (list b-hats "b" "linespoints lw 1 lc rgb '#11cc62'"))
      :title "Unit sphere, B' equatorial" :axes (list b-J a-J))
     (gnuplot 1 (length a-hats) #'(lambda (plot point) (declare (ignore plot))
				    (if (not (eq point :title)) (values point (aref a-dsigmas point))))
	      :xrange (cons 0 (length a-dsigmas)) :yrange (cons 0 a-max)
	      :styles '("linespoints pt 5 lc rgb '#cf0f03'") :title "Energy in segments of a"
	      :xlabel "Segment \#")
     (gnuplot 1 (length b-hats) #'(lambda (plot point) (declare (ignore plot))
				    (if (not (eq point :title)) (values point (aref b-dsigmas point))))
	      :xrange (cons 0 (length b-dsigmas)) :yrange (cons 0 b-max)
	      :styles '("linespoints pt 5 lc rgb '#11cc62'") :title "Energy in segments of b"
	      :xlabel "Segment \#")))))

;; Plots the a and b hats, at fixed axes, taken from a single .DAT file produced by WRITE-HATS-SIGMAS.
(defun plot-hats-fixed (file-name &key (axes nil) (extension "-fixed") (titled nil))
  (multiple-value-bind
   (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (declare (ignore a-dsigmas) (ignore b-dsigmas))
   (plot-tangent-vectors
    (list (list a-hats "a" "points lw 1 lc rgb '#cf0f03'")
	  (list b-hats "b" "points lw 1 lc rgb '#11cc62'"))
    :title (if titled
	       (format nil "Iteration ~D" (parse-integer (subseq file-name (- (length file-name) 4))))
	     nil)
    :prelude (concatenate 'string
			  (format nil "set terminal postscript epsf enhanced color font \"Helvetica,20\"~%")
			  (format nil "set output \"~A~A.eps\"~%" file-name extension)
			  (format nil "unset border~%unset xtics~%unset ytics~%unset arrow~%"))
    :axes axes)))

(defun plot-dsigmas-by-segment (file-name) 
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (declare (ignore a-hats) (ignore b-hats))
   (mirror-image-let
    ((a-max (* 1.05 (apply 'max (loop for i across a-dsigmas collect i)))))
    (multiplot
     (:prelude (format nil "set terminal postscript enhanced color size 10,4~%set output \"~A-seg.eps\"~%" file-name)
	       :columns 2 :rows 1
	       :command-text (format nil "title 'Iteration ~D'"
				     (parse-integer (subseq file-name (- (length file-name) 4)))))
     (gnuplot 1 (length a-dsigmas) #'(lambda (plot point) (declare (ignore plot))
				       (if (not (eq point :title)) (values point (aref a-dsigmas point))))
	      :xrange (cons 0 (length a-dsigmas)) :yrange (cons 0 a-max)
	      :styles '("linespoints pt 5 lc rgb '#cf0f03'") :title "Energy in segments of a"
	      :xlabel "Segment \#")
     (gnuplot 1 (length b-dsigmas) #'(lambda (plot point) (declare (ignore plot))
				       (if (not (eq point :title)) (values point (aref b-dsigmas point))))
	      :xrange (cons 0 (length b-dsigmas)) :yrange (cons 0 b-max)
	      :styles '("linespoints pt 5 lc rgb '#11cc62'") :title "Energy in segments of b"
	      :xlabel "Segment \#")))))
  
(defun plot-dsigmas-by-length (file-name) 
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (declare (ignore a-hats) (ignore b-hats))
   (mirror-image-let
    ((a-max (* 1.05 (apply 'max (loop for i across a-dsigmas collect i))))
     (a-sigmas (dsigma-to-sigma a-dsigmas)))
    (multiplot
     (:prelude (format nil "set terminal postscript enhanced color size 10,4~%set output \"~A-len.eps\"~%" file-name)
	       :columns 2 :rows 1
	       :command-text (format nil "title 'Iteration ~D'"
				     (parse-integer (subseq file-name (- (length file-name) 4)))))
     (gnuplot 1 (length a-dsigmas) #'(lambda (plot point) (declare (ignore plot))
				       (if (not (eq point :title))
					   (values (aref a-sigmas point) (aref a-dsigmas point))))
	      :xrange (cons 0 (aref a-sigmas (1- (length a-sigmas)))) :yrange (cons 0 a-max)
	      :styles '("lines lc rgb '#cf0f03'") :title "Energy by length in a"
	      :xlabel "L_a")
     (gnuplot 1 (length b-dsigmas) #'(lambda (plot point) (declare (ignore plot))
				       (if (not (eq point :title))
					   (values (aref b-sigmas point) (aref b-dsigmas point))))
	      :xrange (cons 0 (aref b-sigmas (1- (length b-sigmas)))) :yrange (cons 0 b-max)
	      :styles '("lines lc rgb '#11cc62'") :title "Energy by length in b"
	      :xlabel "L_b")))))

;; Plots the a and b hats in a series taken from multiple .DAT files produced by WRITE-HATS-DSIGMAS.
;; For each STEP from :START to MAX-STEP, the hats of that iteration of backreaction are plotted beside one
;; another. The entire series is read from the prefixes (everything before the -hs) given by FILE-NAME.
;;    :*-I set the index in a' and b' which we label with an arrow.
;;    :*-BASE-L set the length of the first (non-capped) segment of the respective arrows.
;;    :*-HEAD-L set the length of the second (capped) segment of the respective arrows.
;;    :OFFSET must be a list of lists, each of which is like (dx,dy), which is the distance to displace the
;;        second (capped) segment of one of the arrows. :OFFSET must always be length 2*2*rows. To change
;;        the a' arrow of row r, column c, use offset's (4*r+2*c)th index. For b', use the (4*r+2*c+1)th index.
(defvar zerol (list 0.0 0.0))
(defvar offsets (list (list "kt0100" (list (list 0.1 -0.35) zerol (list -0.25 -0.2) (list 0.1 0.0) (list 0.1 -0.1) zerol
					   (list 0.1 -0.1) zerol zerol zerol zerol zerol))
		      (list "bkt0100" nil)
		      (list "2bkt0100" (list zerol zerol zerol zerol (list -0.3 0.0) zerol
					     (list -0.35 0.0) zerol (list -0.3 0.2) zerol (list -0.3 0.2) zerol))
		      (list "mkt0100" (list (list 0.1 -0.35) zerol (list -0.25 -0.2) zerol zerol (list -0.3 0.1)
					    zerol (list -0.25 -0.15) zerol (list -0.3 0.15) zerol (list -0.25 -0.15)))))
(defun plot-hats-column (file-name step max-step &key (start 0) (a-color "rgb '#cf0f03'") (b-color "rgb '#11cc62'")
				   (a-i 0) (b-i 0) (offset nil) (a-base-l 0.25) (b-base-l 0.25) (a-head-l 0.4) (b-head-l 0.4))
  (unless (integerp (/ (- max-step start) step)) (error "The ratio MAX-STEP/STEP must be an integer"))
  (mirror-image-let
   ((a-hats) (a-dsigmas) (a-J) (a-arrow-1) (a-arrow-2) (a-arrow-3) (a-arrow-base) (a-arrow-mid) (a-arrow-head) (dd-a) (a-offset))
   (let* ((ax) (ay) (az) (bx) (by) (bz) ;these are the vectors for Mollweide projections for the a' and b' columns
	  (a-inc 0) (b-inc 1) ;these are for calculating which index of :OFFSET to use
	  (rows (1+ (/ (- max-step start) step)))
	  (plot-prelude (format nil "unset border~%unset xtics~%unset ytics~%unset arrow~%")) ;we want only the hats+ellipse
	  (arrow-tip (format nil "size 0.08,35 filled front"))
	  (arrow-width 4)
	  (multiplot-prelude
	   (concatenate 'string
			(format nil "set terminal postscript enhanced color size 10,~D font 18~%" (* 3.5 rows))
			(format nil "set output \"~A-h.eps\"~%" file-name)
			;(format nil "set style line 1 lc rgbcolor \"orange-red\"~%") ;style for the a'
			;(format nil "set style line 2 lc rgbcolor \"forest-green\"~%") ;style for the b'
			(format nil "set label at screen 0.275,0.975 center \"initial A' near equator\"~%")
			(format nil "set label at screen 0.775,0.975 center \"initial B' near equator\"~%")
			(apply #'concatenate 'string
			       (loop for i from start to max-step by step collect
				     (format nil "set label at screen 0.025,~3$ center rotate \"step ~D\"~%"
					     (* 1.0 (- 1.0 0.025 (/ (+ 0.5 (/ i step)) (1+ (/ max-step step))))) i))))))
     (when offset
       (unless (= (length offset) (* 2 2 rows)) (error (format nil ":OFFSET must be length ~D" (* 2 2 rows)))))
     (flet
      ((list+ (v1 v2) (map 'list #'(lambda (e1 e2) (+ e1 e2)) v1 v2)))
      (multiplot
       (:prelude multiplot-prelude :columns 2 :rows rows :command-text "offset 0.025,-0.025")
       (loop for i from start to max-step by step
	     do (multiple-value-setq
		 (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name i)))
	     do (when (= i start)
		  (mirror-images (setf a-J (compute-hats-angular-momentum a-hats (dsigma-to-sigma a-dsigmas))))
		  (multiple-value-setq (ax ay az) (choose-Mollweide-axes (list a-J (3vector- b-J))))
		  (multiple-value-setq (bx by bz) (choose-Mollweide-axes (list b-J a-J))))
	     do (mirror-images
		 (setf a-arrow-1 (aref a-hats a-i)
		       dd-a (3vector-normalize (3vector- (aref a-hats (+ 20 a-i)) a-arrow-1))
		       a-arrow-2 (3vector-normalize
				  (3vector+ a-arrow-1
					    (3vector-normalize 
					     (3vector-cross-product a-arrow-1 dd-a)
					     a-base-l)))
		       a-arrow-3 (3vector+ a-arrow-2 (3vector-normalize dd-a a-head-l))))
	     do (mirror-images
		 (setf a-arrow-base (multiple-value-list (mollweide-vector a-arrow-1 ax ay az))
		       a-offset (if offset (nth (+ (* 4 (/ i step)) a-inc) offset) (list 0.0 0.0))
		       a-arrow-mid (list+
				    (multiple-value-list (mollweide-vector a-arrow-2 ax ay az))
				    a-offset)
		       a-arrow-head (list+
				     (multiple-value-list (mollweide-vector a-arrow-3 ax ay az))
				     a-offset)))
	     do (plot-tangent-vectors
		 (list (list a-hats "a" (format nil "linespoints lc ~A" a-color))
		       (list b-hats "b" (format nil "linespoints lc ~A" b-color)))
		 :key nil :loop t :font-size 18 :axes (list a-J (3vector- b-J)) :solid t
		 :prelude (concatenate
			   'string plot-prelude
			   (format nil "set arrow 1 from ~{~S~^,~} to ~{~S~^,~} nohead lc ~A lw ~D~%"
				   a-arrow-base a-arrow-mid a-color arrow-width)
			   (format nil "set arrow 2 from ~{~S~^,~} to ~{~S~^,~} head ~A lc ~A lw ~D~%"
				   a-arrow-mid a-arrow-head arrow-tip a-color arrow-width)
			   (format nil "set arrow 3 from ~{~S~^,~} to ~{~S~^,~} nohead lc ~A lw ~D~%"
				   b-arrow-base b-arrow-mid b-color arrow-width)
			   (format nil "set arrow 4 from ~{~S~^,~} to ~{~S~^,~} head ~A lc ~A lw ~D~%"
				   b-arrow-mid b-arrow-head arrow-tip b-color arrow-width)))
	     do (mirror-images
		 (setf a-arrow-base (multiple-value-list (mollweide-vector a-arrow-1 bx by bz))
		       a-offset (if offset (nth (+ (* 4 (/ i step)) 2 a-inc) offset) (list 0.0 0.0))
		       a-arrow-mid (list+ 
				    (multiple-value-list (mollweide-vector a-arrow-2 bx by bz))
				    a-offset)
		       a-arrow-head (list+
				     (multiple-value-list (mollweide-vector a-arrow-3 bx by bz))
				     a-offset)))
	     do (plot-tangent-vectors
		 (list (list a-hats "a" (format nil "linespoints lc ~A" a-color))
		       (list b-hats "b" (format nil "linespoints lc ~A" b-color)))
		 :key nil :loop t :font-size 18 :axes (list b-J a-J) :solid t
		 :prelude (concatenate
			   'string plot-prelude
			   (format nil "set arrow 1 from ~{~S~^,~} to ~{~S~^,~} nohead lc ~A lw ~D~%"
				   a-arrow-base a-arrow-mid a-color arrow-width)
			   (format nil "set arrow 2 from ~{~S~^,~} to ~{~S~^,~} head ~A lc ~A lw ~D~%"
				   a-arrow-mid a-arrow-head arrow-tip a-color arrow-width)
			   (format nil "set arrow 3 from ~{~S~^,~} to ~{~S~^,~} nohead lc ~A lw ~D~%"
				   b-arrow-base b-arrow-mid b-color arrow-width)
			   (format nil "set arrow 4 from ~{~S~^,~} to ~{~S~^,~} head ~A lc ~A lw ~D~%"
				   b-arrow-mid b-arrow-head arrow-tip b-color arrow-width)))))))))

;; Plots the a and b dsigmas in a series taken from multiple .DAT files produced by WRITE-HATS-DSIGMAS.
;; For each STEP from START to MAX-STEP, the dsigmas of that iteration of backreaction are plotted beside
;; another. The entire series is read from the prefixes (everything before the -hs) given by FILE-NAME
(defun plot-dsigmas-column (file-name step max-step &optional (start 0))
  (unless (integerp (/ (- max-step start) step)) (error "The ratio MAX-STEP/STEP must be an integer"))
  (mirror-image-let
   ((a-hats) (a-dsigmas) (a-max))
   (multiplot
    (:prelude (format nil "set terminal postscript enhanced color size 10,~D~%set output \"~A-s.eps\"~%"
		      (* 3.5 (1+ (/ (- max-step start) step))) file-name)
	      :columns 2 :rows (1+ (/ (- max-step start) step)))
    (loop for i from start to max-step by step
	  do (multiple-value-setq
	      (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name i)))
	  do (mirror-images (setf a-max (* 1.05 (vector-mm a-dsigmas))))
	  do (gnuplot 1 (length a-hats)
		      #'(lambda (plot point) (declare (ignore plot))
			  (if (not (eq point :title)) (values point (aref a-dsigmas point))))
		      :xrange (cons 0 (length a-hats)) :yrange (cons 0 a-max)
		      :styles '("linespoints pt 5 lc rgb '#cf0f03'") :xlabel (if (= i max-step) "Segment \#" "")
		      :title (if (= i start) "Energy in each segment of a" ""))
	  do (gnuplot 1 (length b-hats)
		      #'(lambda (plot point) (declare (ignore plot))
			  (if (not (eq point :title)) (values point (aref b-dsigmas point))))
		      :xrange (cons 0 (length b-hats)) :yrange (cons 0 b-max)
		      :styles '("linespoints pt 5 lc rgb '#11cc62'") :xlabel (if (= i max-step) "Segment \#" "")
		      :title (if (= i start) "Energy in each segment of b"))))))

;; Plots the a and b dsigmas in a series taken from multiple .DAT files produced by WRITE-HATS-DSIGMAS.
;; For each STEP from START to MAX-STEP, the dsigmas of that iteration of backreaction are plotted atop one
;; another. The entire series is read from the prefixes (everything before the -hs) given by FILE-NAME.
(defun plot-dsigmas-stack (file-name step max-step &optional (start 0))
  (unless (integerp (/ (- max-step start) step)) (error "The ratio MAX-STEP/STEP must be an integer"))
  (mirror-image-let
   ((a-hats) (a-dsigmas-l (make-array (1+ (/ (- max-step start) step)))) (a-dsigmas) (a-max))
   (loop for i from start to max-step by step
	 do (multiple-value-setq
	     (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name i)))
	 do (mirror-images (setf (aref a-dsigmas-l (/ (- i start) step)) a-dsigmas)))
   (mirror-images (setf a-max (* 1.05 (vector-mm (aref a-dsigmas-l 0)))))
   (multiplot
    (:prelude (format nil "set ~A size ~D,~D font ~D~%set output \"~A-ss.eps\"~%"
		      "terminal postscript enhanced color" 7.5 2.5 10 file-name)
	      :columns 2)
    (gnuplot (length a-dsigmas-l) (length a-hats)
	     #'(lambda (plot point)
		 (if (not (eq point :title)) (values point (aref (aref a-dsigmas-l plot) point))))
	     :xrange (cons 0 (length a-hats)) :yrange (cons 0 a-max)
	     :styles (subseq (styles-shades 207 15 3 (1+ (length a-dsigmas-l)) "lines lt 1 lw 2")
			     0 (length a-dsigmas-l))
	     :xlabel "Segment \#" :title "Energy in each segment of A" :font-size "10")
    (gnuplot (length b-dsigmas-l) (length b-hats)
	     #'(lambda (plot point)
		 (if (not (eq point :title)) (values point (aref (aref b-dsigmas-l plot) point))))
	     :xrange (cons 0 (length b-hats)) :yrange (cons 0 b-max)
	     :styles (subseq (styles-shades 17 204 98 (1+ (length b-dsigmas-l)) "lines lt 1 lw 2")
			     0 (length b-dsigmas-l))
	     :xlabel "Segment \#" :title "Energy in each segment of B" :font-size "10"))))

;; Precisely as PLOT-DSIGMAS-STACK, but now the energy of each segment is normalized to the value of the energy
;; of that segment for the file selected by START.
(defun plot-dsigmas-stack-normed (file-name step max-step &optional (start 0))
  (unless (integerp (/ (- max-step start) step)) (error "The ratio MAX-STEP/STEP must be an integer"))
  (mirror-image-let
   ((a-hats) (a-dsigmas-l (make-array (1+ (/ (- max-step start) step)))) (a-dsigmas) (a-min) (a-norm))
   (loop for i from start to max-step by step
	 do (multiple-value-setq
	     (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name i)))
	 do (when (= i start) (mirror-images (setf a-norm (copy-seq a-dsigmas))))
	 do (mirror-images
	     (dotimes (j (length a-dsigmas)) (setf (aref a-dsigmas j) (/ (aref a-dsigmas j) (aref a-norm j)))))
	 do (mirror-images (setf (aref a-dsigmas-l (/ (- i start) step)) a-dsigmas)))
   (mirror-images (setf a-min (* 0.95 (vector-mm (aref a-dsigmas-l (1- (length a-dsigmas-l))) :max nil))))
   (multiplot
    (:prelude (format nil "set ~A size ~D,~D font ~D~%set output \"~A-ssn.eps\"~%"
		      "terminal postscript enhanced color" 7.5 2.5 10 file-name)
	      :columns 2)
    (gnuplot (length a-dsigmas-l) (length a-hats)
	     #'(lambda (plot point)
		 (if (not (eq point :title)) (values point (aref (aref a-dsigmas-l plot) point))))
	     :xrange (cons 575 625) :yrange (cons a-min 1.05)
	     :styles (styles-shades 208 35 27 (length a-dsigmas-l) "points pt 1 ps 0.5")
	     :xlabel "Segment \#" :title "Energy in each segment of a" :font-size "10")
    (gnuplot (length b-dsigmas-l) (length b-hats)
	     #'(lambda (plot point)
		 (if (not (eq point :title)) (values point (aref (aref b-dsigmas-l plot) point))))
	     :xrange (cons 575 625) :yrange (cons b-min 1.05)
	     :styles (styles-shades 0 111 82 (length b-dsigmas-l) "points pt 2 ps 0.5")
	     :xlabel "Segment \#" :title "Energy in each segment of b" :font-size "10"))))

;; Plot the gamma, L, J, and potentially kinks for a loop from a .DAT file whose formatting is assumed to be of the sort
;; produced by DO-BACKREACTION. FILE-NAME should be without the -GLJK.DAT part, i.e. just the loop identifier.
(defun plot-gljk (file-name)
  (multiple-value-bind
   (gammas Ls Js kinks) (read-gljk (format nil "~A-gljk.dat" file-name))
   (gnuplot 1 (length gammas) #'(lambda (plot point) (declare (ignore plot))
				  (if (not (eq point :title)) (values (1+ point) (nth point gammas))))
	    :yrange (cons (/ (apply 'min gammas) 1.05) (* (apply 'max gammas) 1.05))
	    :xlabel "Iterations" :ylabel "{/Symbol G}" :styles '("lines lt 1 lw 8 lc rgb '#9bcfcc'")
	    :color t :epsf (format nil "~A-ag.eps" file-name) :font-size "22")
   (gnuplot 1 (length Ls) #'(lambda (plot point) (declare (ignore plot))
			      (if (not (eq point :title)) (values point (nth point Ls))))
	    :yrange (cons (/ (apply 'min Ls) 1.05) (* (apply 'max Ls) 1.05))
	    :xlabel "Iterations" :ylabel "L" :styles '("lines lt 1 lw 8 lc rgb '#8f4099'")
	    :color t :epsf (format nil "~A-al.eps" file-name) :font-size "22")
   (gnuplot 1 (length Js) #'(lambda (plot point) (declare (ignore plot))
			      (if (not (eq point :title)) (values point (nth point Js))))
	    :yrange (cons (/ (apply 'min Js) 1.05) (* (apply 'max Js) 1.05))
	    :xlabel "Iterations" :ylabel "J_{tot}" :styles '("lines lt 8 lw 5 lc rgb '#fa6700'")
	    :color t :epsf (format nil "~A-aj.eps" file-name) :font-size "22")
   (if kinks
       (let* ((kink-num (/ (length (car kinks)) 2))
	      (a-kinks (map 'list #'(lambda (x) (subseq x 0 kink-num)) kinks))
	      (b-kinks (map 'list #'(lambda (x) (subseq x kink-num)) kinks)))
	 (gnuplot kink-num (length a-kinks)
		  #'(lambda (plot point)
		      (if (not (eq point :title)) (values point (nth plot (nth point a-kinks)))))
		  :color t :ylabel "Angle (rads)" :epsf (format nil "~A-ka.eps" file-name)
		  :styles (loop repeat kink-num collect "lines lt 1 lw 8 lc rgb '#cf0f03'")
		  :xlabel "Iterations" :yrange (cons (/ (apply 'min (apply 'append a-kinks)) 1.05)
						     (* (apply 'max (apply 'append a-kinks)) 1.05))
		  :font-size "32" :title "Kinks in a of the broken loop")
	 (gnuplot kink-num (length b-kinks)
		  #'(lambda (plot point)
		      (if (not (eq point :title)) (values point (nth plot (nth point b-kinks)))))
		  :color t :ylabel "Angle (rads)" :epsf (format nil "~A-kb.eps" file-name)
		  :styles (loop repeat kink-num collect "lines lt 1 lw 8 lc rgb '#11cc62'")
		  :xlabel "Iterations" :yrange (cons (/ (apply 'min (apply 'append b-kinks)) 1.05)
						     (* (apply 'max (apply 'append b-kinks)) 1.05))
		  :font-size "32" :title "Kinks in b of the broken loop")))))
  
;; Plots the spatial extent of the a and b worldsheet functions for a single .DAT file of the type produced by 
;; WRITE-HATS-DSIGMAS.
(defun plot-spatial-ab (file-name)
  (multiple-value-bind
   (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (mirror-image-let
    ((a-spatial (make-array 1 :initial-element (3vector-scale (aref a-hats 0) (aref a-dsigmas 0))
			    :fill-pointer 1 :adjustable t)))
    (mirror-images
     (loop for i from 1 upto (1- (length a-hats))
	   do (vector-push-extend (3vector+ (3vector-scale (aref a-hats i) (aref a-dsigmas i)) (aref a-spatial (1- i)))
				  a-spatial))
     (vector-push-extend (aref a-spatial 0) a-spatial))
    (multiplot
     (:prelude
      (format nil "set terminal postscript enhanced color~%set output \"~A-ab.eps\"~%set xyplane 0~%set view equal xy~%"
	      file-name) :columns 2
	      :command-text (format nil "title 'Iteration ~D'"
				    (parse-integer (subseq file-name (- (length file-name) 4)))))
     (gnuplot 1 (length a-spatial) #'(lambda (plot point) (declare (ignore plot))
				       (if (not (eq point :title))
					   (values-list (3vector-list (aref a-spatial point)))))
	      :styles '("lines lt 1 lc rgb '#cf0f03'") :title "Spatial A" :color t :3D t)
     (gnuplot 1 (length b-spatial) #'(lambda (plot point) (declare (ignore plot))
				       (if (not (eq point :title))
					   (values-list (3vector-list (aref b-spatial point)))))
	      :styles '("lines lt 1 lc rgb '#11cc62'") :title "Spatial B" :color t :3D t)))))

;; Plots the gravitational power on the sphere from a filename
(defun plot-gravitational-power-from-file (file-name n steps &optional (cbmax nil))
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (plot-gravitational-power a-hats (dsigma-to-sigma a-dsigmas) b-hats (dsigma-to-sigma b-dsigmas) n steps
			       :epsf (format nil "~A-power-n10E~Dby~D.eps" file-name (floor (log n 10)) steps)
			       :color t :cbmax cbmax)))

;; Prints out the cusp parameters (unit vector in the direction, magnitude of a'', magnitude of b'') for a given data file.
(defun print-cusps (file-name)
  (multiple-value-bind
   (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (let* ((cusp-parameters (find-discrete-cusps a-hats a-dsigmas b-hats b-dsigmas))
	  (current nil) (dda nil) (ddb nil))
     (format t "~%Iteration ~D:" (parse-integer (subseq file-name (- (length file-name) 4))))
     (loop while cusp-parameters
	   do (setf current (car cusp-parameters)
		    cusp-parameters (cdr cusp-parameters)
		    dda (nth 2 current)
		    ddb (nth 3 current))
	   do (format t "~%~%Cusp at ~S:~%Direction = ~S~%Naive power = ~S~%a'' = ~S, |a''| = ~8$~%b'' = ~S, |b''| = ~8$"
		      (nth 0 current) (nth 1 current)
		      (expt (* (3vector-length dda) (3vector-length ddb)) (- (/ 3)))
		      dda (3vector-length dda) ddb (3vector-length ddb))))))

;; Prints out the naive and exact powers of all of the cusps on a given loop, as well as the segment number of A and B
;; just below the cusp and the loop's invariant length. :THETA sets the angular width of the cusp beam we count in the
;; exact power. :VERBOSE sets if the output contains words or not.
(defun print-cusp-powers (file-name &key (theta 0.1) (verbose nil))
  (multiple-value-bind
   (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A.dat" file-name))
   (let* ((cusp-parameters (find-discrete-cusps a-hats a-dsigmas b-hats b-dsigmas))
;;	  (observer nil)
	  (current nil)
	  (iter (parse-integer (subseq file-name (- (length file-name) 4))))
	  (L (total-sigma a-dsigmas))	;TOTAL-POWER-FROM-CUSP-DISCRETE doesn't include the 1/L factor
	  (pf (if verbose "~%Cusp at (~3$,~3$), loop length ~8$~%Naive power: ~8$~%Exact power: ~8$"
		"~%~3$ ~3$ ~8$ ~8$ ~8$")))
     (if verbose
	 (format t "~%Iteration ~D:" iter)
       (format t "~%~D" iter))
     (loop while cusp-parameters
	   do (setf current (car cusp-parameters)
;;		    observer (3vector-normalize (3vector+ (nth 1 current) (3vector-normalize (nth 2 current) 0.1)))
		    cusp-parameters (cdr cusp-parameters))
	   do (format t pf
		      (car (nth 0 current)) (cadr (nth 0 current)) L
		      (expt (* (3vector-length (nth 2 current)) (3vector-length (nth 3 current))) (- (/ 3)))
		      (total-power-from-cusp-discrete current L :theta theta))))))

;; Compares the change in a' and b' for the segment INDEX between iterations ITER-I and ITER-F to the a'' and b''
;; during ITER-I. Loads the data from FILE-NAME.
(defun compare-hat-changes (file-name a-index b-index iter-i iter-f &key (verbose nil))
  (multiple-value-bind
   (a-hats-i a-dsigmas-i b-hats-i b-dsigmas-i) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name iter-i))
   (declare (ignore a-dsigmas-i b-dsigmas-i))
   (multiple-value-bind
    (a-hats-f a-dsigmas-f b-hats-f b-dsigmas-f) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name iter-f))
    (declare (ignore a-dsigmas-f b-dsigmas-f))
    (mirror-image-let
     ((delta-a-p (3vector-normalize (3vector- (aref a-hats-f a-index) (aref a-hats-i a-index))))
      (a-pp (3vector-normalize (3vector- (aref a-hats-i (mod (1+ a-index) (length a-hats-i))) (aref a-hats-i a-index)))))
     (let ((pf (if verbose "da' = ~S~%db' = ~S~%a'' = ~S~%b'' = ~S~%" "~S~%~S~%~S~%~S~%")))
       (format t pf delta-a-p delta-b-p a-pp b-pp))))))

(defun track-cusp-motion (file-name step max-step &key (start 0) (cusp-number 0) (degrees nil))
  (let ((a-hats) (a-dsigmas) (b-hats) (b-dsigmas) (initial-direction) (current-direction))
    (loop for i from start to max-step by step do
	  (multiple-value-setq
	   (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name i)))
	  (if (= i start) (setf initial-direction
				(nth 1 (nth cusp-number (find-discrete-cusps a-hats a-dsigmas b-hats b-dsigmas)))))
	  (setf current-direction
		(nth 1 (nth cusp-number (find-discrete-cusps a-hats a-dsigmas b-hats b-dsigmas))))
	  (format t "~D ~5$~%" i (* (if degrees (/ 360 (* 2 pi)) 1.0)
				    (spherical-angle current-direction initial-direction))))))

;; Dumb function for comparing the dsigma and hat of the initial setup and some later step for the segments just before
;; and after a kink. "dumb" because it assumes that the kink happens between the last segment and the first segment.
(defun dumb-kink-comparison (file-name later-step &optional (verbose nil))
  (mirror-image-let
   ((a-hat-im) (a-hat-ip) (a-hat-fm) (a-hat-fp) (a-dsigma-im) (a-dsigma-ip) (a-dsigma-fm) (a-dsigma-fp)
    (a-hats) (a-dsigmas) (a-L))
   (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas) (read-hats-dsigmas (format nil "~A-hs-000.dat" file-name)))
   (mirror-images (setf a-L (1- (length a-hats))
			a-hat-im (aref a-hats a-L)
			a-hat-ip (aref a-hats 0)
			a-dsigma-im (aref a-dsigmas a-L)
			a-dsigma-ip (aref a-dsigmas 0)))
   (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
			(read-hats-dsigmas (format nil "~A-hs-~4,'0D.dat" file-name later-step)))
   (mirror-images (setf a-hat-fm (aref a-hats a-L)
			a-hat-fp (aref a-hats 0)
			a-dsigma-fm (aref a-dsigmas a-L)
			a-dsigma-fp (aref a-dsigmas 0)))
   (if verbose
       (mirror-images (format t "Changes in ~A:~%At ~8$ (-) ~8$ (+)~%dsigmas: ~8$ (-), ~8$ (+)~%angle: ~8$ (-), ~8$ (+)~%"
			      :a (/ a-dsigma-fm 2) (/ a-dsigma-fp 2)
			      (/ a-dsigma-fm a-dsigma-im) (/ a-dsigma-fp a-dsigma-ip)
			      (spherical-angle a-hat-fm a-hat-im) (spherical-angle a-hat-fp a-hat-ip)))
     (mirror-images (format t "~8$ ~8$ ~8$ ~8$ ~8$ ~8$~%"
			    (/ a-dsigma-fm) (/ a-dsigma-fp)
			    (/ a-dsigma-fm a-dsigma-im) (/ a-dsigma-fp a-dsigma-ip)
			    (spherical-angle a-hat-fm a-hat-im) (spherical-angle a-hat-fp a-hat-ip))))))

;; Accepts some set of hats (as a vector) and returns the spherical angles between successive hats (as a list)
;; If a value is given for LARGEST, then we return only the LARGEST largest angles in the list.
(defun hat-angles (hats &key (largest nil))
  (let ((len (length hats))
	(angles (list)))
    (dotimes (i (1- len))
      (push (spherical-angle (aref hats (1+ i)) (aref hats i)) angles))
    (push (spherical-angle (aref hats (1- len)) (aref hats 0)) angles)
    (when largest (if (< largest len)
		      (setf angles (subseq (sort angles #'>) 0 largest))
		    (warn "LARGEST > LEN in HAT-ANGLES, so just returning the full unsorted list.")))
    (reverse angles)))

;; Follows after FIND-SMOOTH-CUSPS, but presuming you enter in the a and b hats and dsigma data directly.
(defun find-discrete-cusps (a-hats a-dsigmas b-hats b-dsigmas)
  (let ((n-a (length a-hats))
	(n-b (length b-hats))
	(results nil)
	(cusp-point nil))
    (dotimes (i n-a)
      (loop with a1 = (aref a-hats i)
	    with a2 = (aref a-hats (mod (1+ i) n-a))
	    for j below n-b
	    for b1 = (aref b-hats j)
	    for b2 = (aref b-hats (mod (1+ j) n-b))
	    do (setf cusp-point (spherical-intersection a1 a2 b1 b2))
	    when cusp-point
	    do (pushnew (list (list (/ i n-a) (/ j n-b)) ;x-a and x-b
			      cusp-point ;cusp direction
			      (3vector-scale (3vector- a2 a1)
					     (/ (sqrt (* (aref a-dsigmas i) (aref a-dsigmas (mod (1+ i) n-a))))));a''
			      (3vector-scale (3vector- b2 b1)
					     (/ (sqrt (* (aref b-dsigmas j) (aref b-dsigmas (mod (1+ j) n-b)))))));b''
			results
			:test #'(lambda (t1 t2) (equal (car t1) (car t2))))))
    results))

;; Produces a list appropriate for use in a :STYLES command of a gnuplot function. Takes an initial color's RGB
;; values RI, GI, BI (in decimal), then produces evenly-spaced shades on a straight line of evenly-spaced steps
;; terminating at a final color's RGB values RF, GF, BF (also in decimal). STEPS sets the length of the list.
;; DATA-STYLING lets the user specify line/point types, etc.
(defun styles-spectrum (ri gi bi rf gf bf steps data-styling)
  (loop for k from 0 to (1- steps)
	collect (format nil "~A lc rgb '#~2,'0X~2,'0X~2,'0X'" data-styling
			(+ ri (floor (* k (/ (- rf ri) (1- steps)))))
			(+ gi (floor (* k (/ (- gf gi) (1- steps)))))
			(+ bi (floor (* k (/ (- bf bi) (1- steps))))))))
(defun styles-shades (r g b steps data-styling)
  (styles-spectrum r g b 0 0 0 steps data-styling))
(defun styles-tints (r g b steps data-styling)
  (styles-spectrum r g b 255 255 255 steps data-styling))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; TESTING FUNCTIONS, ETC.


;;Sets up an intersection ring for a particular choice of observer diamond.
;; the -A-RING treats B as the observer direction of motion, and v.v.
(mirror-images
(defun set-up-a-ring (a-hats a-dsigmas b-hats b-dsigmas iter-a iter-b)
  (mirror-images (check-dsigmas a-dsigmas))
  (setf *other-vector* (3to4vector (aref a-hats iter-a) 1.0)
	*observer-vector* (3to4vector (aref b-hats iter-b) 1.0)
	*observer-L* (/ (aref b-dsigmas iter-b) 2.0)
	*observer-diamond* nil)
  (create-initial-ring :a-hats a-hats :a-dsigmas a-dsigmas ;Sets *observer-diamond*
		       :b-hats b-hats :b-dsigmas b-dsigmas
		       :start-a iter-a :start-b iter-b :start-x (- *observer-L*))
  (make-null-calendar *observer-L* (length a-hats))
  (mirror-images (setf *most-sw* (find-most-sw *observer-diamond*))))
)

;;Apply random perturbations of magnitude LEVEL to a list of hats
(defun perturb-hats (hats &optional (level 1e-15))
  (let* ((n (length hats))
	 (new-hats (make-array n)))
    (dotimes (i n)
      (setf (aref new-hats i)
	    (3vector+ (aref hats i) (3vector-scale (random-direction) level))))
    new-hats))

;;Check the stability of a loop by computing the dA' and dB' of all segments for several random perturbations
;;of the worldsheet, then finding for each segment the maximum difference between any two outcomes. This max
;;difference is normalized with respect to that segment's change for the unperturbed worldsheet.
(defun test-stability (a-hats a-dsigmas b-hats b-dsigmas &key (level 1e-15) (count 4))
  (mirror-image-let* ((n-a (length a-hats))
		      (a-segment-changes (make-array n-a))
		      (a-base-changes (make-array n-a))
		      (max-a-difference (make-array n-a))
		      (a-results (make-array count))
		      (a-differences (make-array (/ (* count (1- count)) 2)))
		      new-a-hats
		      iter)
    (mirror-images
     (dotimes (index n-a)
       (setf (aref a-base-changes index)
	     (4vector-Euclidian-length (a-segment-change a-hats a-dsigmas b-hats b-dsigmas index)))))
    (dotimes (i count) ;make the perturbed changes to A and B
      (mirror-images (setf new-a-hats (perturb-hats a-hats level)))
      (mirror-images
       (dotimes (index n-a)
	 (setf (aref a-segment-changes index) (a-segment-change new-a-hats a-dsigmas new-b-hats b-dsigmas index)))
       (setf (aref a-results i) (copy-seq a-segment-changes))))
    (dotimes (j count) ;for all changes, find the differences by segment
      (dotimes (k j)
	(setf iter (+ (/ (* j (1- j)) 2) k))
	(mirror-images
	 (setf (aref a-differences iter)
	       (map 'vector #'(lambda (x y z) (/ (4vector-Euclidian-length (4vector- x y)) z))
		    (aref a-results j) (aref a-results k) a-base-changes)))))
    (mirror-images ;for all segments, find the max change
     (dotimes (i n-a)
       (setf (aref max-a-difference i)
	     (apply #'max (loop for j from 0 to (1- (length a-differences)) 
				collect (aref (aref a-differences j) i))))))
    (values max-a-difference max-b-difference)))

(defun test-gv-backreaction (theta count &key perturbation (NGmu 0.0001) (L 4.0))
  (let* ((phi (/ theta 2))
	 (h (sin phi))
	 (w (cos phi))
	 (wwhh (expt (/ w h) 2))
	 (r (1+ (* 2 wwhh)))
	 (rr (1+ (/ 2 wwhh))))
    (format t "8 NGmu ln RR' = ~S~%" (* 8 NGmu (log (* r rr))))
    (format t "dE/E = ~S" (* -16 NGmu (+ (/ (log w) (expt h 2)) (/ (log h) (expt w 2)))))
    (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	(make-gv-hats theta count :perturbation perturbation :L L)
     (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
			  (gravitational-backreaction a-hats a-dsigmas b-hats b-dsigmas NGmu))
     (format t "Length fraction lost: ~S~%~%" (- 1 (/ (total-sigma b-dsigmas) L)))
     )))
		     
;;Explict min-x, max-x
(defun find-correction-1 (source-diamond &key min-x max-x)
  (if (diamond-sw source-diamond)
      (if (diamond-se source-diamond)	;SW&SE: past-type.
	  (past-crossing source-diamond min-x max-x)
	(b-crossing source-diamond min-x max-x))	;SW&NE (u-type/b-type). 
    (if (diamond-se source-diamond)	;SE&NW (v-type/a-type).
	(a-crossing source-diamond min-x max-x)
      (future-crossing source-diamond min-x max-x)))) ;NW&NE (future-type)

;;Infinitesimal amount of observer line
(defun walk-ring-change-1 (min-x &optional (epsilon 1e-6))
  (let ((diamond *most-se*)
        (delta-null-vector (make-zero-4vector))) ;The correction to the null vector is initially zero.
    (loop
      do (4vector-incf delta-null-vector (find-correction-1 diamond :min-x min-x :max-x (+ min-x epsilon)))
      until (eq diamond *most-sw*) ;Stops once we have processed *most-sw*
      do (setf diamond (diamond-e diamond)))
    (4vector-scale delta-null-vector (/ epsilon))))

;;Find rate of change for a single observer position
;;start-x is is the x value on the observer line, between -L and L.
;;For b-segment-change-1, b-hats and b-dsigmas come first in argument list and iter-b comes before iter-a
(mirror-images

(defun a-segment-change-1 (a-hats a-dsigmas b-hats b-dsigmas iter-a iter-b start-x &key (delta 1e-6))
  (mirror-images (check-dsigmas a-dsigmas))
  (let* ((*observer-vector* (3to4vector (aref b-hats iter-b) 1.0)) ;The direction which the observer point moves
	 (*observer-L* (/ (aref b-dsigmas iter-b) 2.0))	      ;The length of the observer line segment
	 *observer-diamond*)
    (mirror-image-let ((*most-sw* (find-most-sw *observer-diamond*)))
      (create-initial-ring :a-hats a-hats :a-dsigmas a-dsigmas :b-hats b-hats :b-dsigmas b-dsigmas
			   :start-a iter-a :start-b iter-b :start-x start-x)
      (walk-ring-change-1 start-x delta))))

(defun test-a-segment-change (theta count &key perturbation (segment 0))
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (make-gv-hats theta count :perturbation perturbation)
    (format t "A hat being modified is ~S" (3to4vector (aref a-hats segment) 1.0))
    (initialize)
    (let ((d-a-prime (a-segment-change a-hats a-dsigmas b-hats b-dsigmas
				       segment)))
      d-a-prime)))

(defun test-a-segment-change-1 (theta count &key perturbation (a-segment 0) (b-segment 0)
				      (observer-b fudge-global-coordinates))
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (make-gv-hats theta count :perturbation perturbation)
    (format t "A hat being modified is ~S" (3to4vector (aref a-hats a-segment) 1.0))
    (initialize)
    (a-segment-change-1 a-hats a-dsigmas b-hats b-dsigmas
			a-segment b-segment observer-b)))

(defun a-segment-change-single-diamond (a-hats a-dsigmas b-hats b-dsigmas iter-a iter-b)
  (mirror-images (check-dsigmas a-dsigmas))
  (let* ((d-a-prime (make-zero-4vector))
	 (*other-vector* (3to4vector (aref a-hats iter-a) 1.0))
	 (*observer-vector* (3to4vector (aref b-hats iter-b) 1.0)) ;The direction which the observer point moves
	 (*observer-L* (/ (aref b-dsigmas iter-b) 2.0))	      ;The length of the observer line segment
	 *observer-diamond*)
    (create-initial-ring :a-hats a-hats :a-dsigmas a-dsigmas :b-hats b-hats :b-dsigmas b-dsigmas
			 :start-a iter-a :start-b iter-b :start-x (- *observer-L*))
    (make-null-calendar *observer-L* (length a-hats))
    (mirror-image-let ((*most-sw* (find-most-sw *observer-diamond*)))
      (4vector-incf d-a-prime (walk-ring))
      (4vector-incf d-a-prime (evolve-observer))
      (loop for d = (diamond-e *observer-diamond*) then next
	    until (eq d *observer-diamond*)
	    for next = (diamond-e d)
	    do (return-diamond d)))
    d-a-prime))
)


;; For debugging purposes. From START-DIAMOND, walks east through the worldsheet, indicating if the intersection line
;; steps N or S as appropriate.      
(defun print-worldsheet-connection (start-diamond)
  (format t "O - ~{~A~^ - ~}~%" (get-worldsheet-connection start-diamond)))
(defun get-worldsheet-connection (start-diamond)
  (let ((test-diamond start-diamond))
    (loop
     collect (if (diamond-ne test-diamond) "N" "S")
     do (setf test-diamond (diamond-e test-diamond))
     until (eq test-diamond start-diamond))))


(mirror-images
(defun number-e (from-diamond to-diamond)
  (let ((test-diamond from-diamond)
	(inc 0))
    (loop until (eq test-diamond to-diamond)
	  do (setf inc (1+ inc)
		   test-diamond (diamond-e test-diamond)))
    inc))
(defun number-a-e (from-diamond to-diamond)
  (let ((test-diamond from-diamond)
	(inc 0))
    (loop until (eq test-diamond to-diamond)
	  do (if (diamond-se test-diamond) (setf inc (1+ inc)))
	  do (setf test-diamond (diamond-e test-diamond)))
    inc))
(defun number-b-e (from-diamond to-diamond)
  (let ((test-diamond from-diamond)
	(inc 0))
    (loop until (eq test-diamond to-diamond)
	  do (if (diamond-ne test-diamond) (setf inc (1+ inc)))
	  do (setf test-diamond (diamond-e test-diamond)))
    inc))
)

;; For producing the uvcd acceleration components on the b segment just below the kink in a broken kt loop.
;; Also prints the c and d acceleration coefficients predicted by theory.
(defun kink-acceleration (&key (a-fraction 0) (a-location 0.5) (phi (/ pi 2)) (psi (/ pi 2)) (alpha 0.1) (L (* 2 pi))
			       (kink-offset 9) (start 400) (end 4000) (step-size 200)
			       (suppress-numerics nil) (suppress-theory nil) (theory-combined t)
			       (change-type #'b-segment-change-1) (change-keys nil))
  (loop for i from start to end by step-size
	do (multiple-value-bind (t-a-hats t-a-dsigmas t-b-hats t-b-dsigmas)
	       (make-broken-kt-hats alpha phi :points i :psi psi :L L)
	     (let* ((a-index (floor (* a-fraction i)))
		    (dv (4vector-scale
			 (apply change-type t-a-hats  t-a-dsigmas t-b-hats t-b-dsigmas
				a-index (- (/ i 2) (1+ kink-offset)) a-location change-keys)
			 0.5))		; divide by 2 to get acceleration from the change to B'
		    (u-basis (4vector-scale (3to4vector (aref t-a-hats a-index) 1.0) 0.5))
		    (v-basis (4vector-scale (3to4vector (aref t-b-hats (- (/ i 2) (1+ kink-offset))) 1.0) 0.5))
		    (c-basis)
		    (d-basis)
		    (pseudo-dv))
	       (multiple-value-setq (c-basis d-basis) (find-pseudo-basis u-basis v-basis u-basis))
	       ;;(format t "u: ~S~%v: ~S~%c: ~S~%d: ~S~%" u-basis v-basis c-basis d-basis)
	       (setf pseudo-dv (cartesian-to-pseudo dv u-basis v-basis c-basis d-basis))
	       ;;(format t "reg dv: ~S~%psd dv: ~S~%" dv pseudo-dv)
	       ;;(format t "~D ~8$~%" i (/ (4vector-c pseudo-dv) (expt i (/ 1.0 3))))
	       ;;(format t "~D ~S ~S ~$ ~$~%" i c-basis d-basis (4dot c-basis c-basis) (4dot d-basis d-basis))
	       (unless suppress-numerics
		 (format t "~D ~8$ ~8$ ~8$ ~8$~%" i (4vector-u pseudo-dv) (4vector-v pseudo-dv)
			 (4vector-c pseudo-dv) (4vector-d pseudo-dv)))
	       (if (and (= i end) suppress-numerics suppress-theory)
		   (format t "~8$ ~8$~%" (* (4vector-c pseudo-dv) (expt (/ (* L (+ kink-offset 0.5)) i) (/ 3.0)))
			   (* (4vector-d pseudo-dv) (expt (/ (* L (+ kink-offset 0.5)) i) (/ 3.0)))))
	       (if (and (= i end) (not suppress-theory))
		   (let* ((dda (4vector-scale (4vector- (aref t-a-hats (mod (1+ a-index) i)) (aref t-a-hats (mod (1- a-index) i)))
					      (/ i (* 2 L))))
			  (da (3to4vector (aref t-a-hats a-index) 1.0))
			  (dbm (3to4vector (aref t-b-hats (- (/ i 2) (1+ kink-offset))) 1.0))
			  (dbp (3to4vector (aref t-b-hats (/ i 2)) 1.0))
			  (Cm (4dot da dbm))
			  (Cp (4dot da dbp))
			  (ddasq (4vector-squared-length-spacelike dda))
			  (dda-c (4vector-c (cartesian-to-pseudo dda u-basis v-basis c-basis d-basis)))
			  (dbp-c (4vector-c (cartesian-to-pseudo dbp u-basis v-basis c-basis d-basis)))
			  (dda-d (4vector-d (cartesian-to-pseudo dda u-basis v-basis c-basis d-basis)))
			  (dbp-d (4vector-d (cartesian-to-pseudo dbp u-basis v-basis c-basis d-basis)))
			  (db-star-dda (+ (* dbp-c dda-c) (* dbp-d dda-d)))
			  (overall (* (- 1) (/ 2 Cp) (expt (/ (* ddasq Cm Cm) 3) (/ 3.0))))
			  (theory-c-db (* overall dbp-c))
			  (theory-d-db (* overall dbp-d))
			  (theory-c-dda (* overall (/ db-star-dda ddasq) dda-c))
			  (theory-d-dda (* overall (/ db-star-dda ddasq) dda-d)))
		     ;;(format t "~8$~%" (* (/ (* 8 Cm) delta) (+ (* (/ dda-c ddasq) (- (/ dCp Cp) (/ dCm Cm))))))
		     (if theory-combined
			 (format t "~8$ ~8$~%" (+ theory-c-db theory-c-dda) (+ theory-d-db theory-d-dda))
		       (format t "~8$ ~8$ ~8$ ~8$~%" theory-c-db theory-d-db theory-c-dda theory-d-dda)))
		 )))))

#| Not updated for changes to bind special variables before calling create-initial-ring

;; For a worldsheet with a kink in B, finds the correction to the observer diamond located at (iter-a,iter-b), which has
;; KINK-OFFSET diamonds between it and the kink (i.e. the diamond just below the kink has KINK-OFFSET of 0), but only considers
;; the N diamonds between the observer and the kink plus :BEYOND*N diamonds on the other side of the kink.
(defun b-segment-change-1-kink
  (a-hats a-dsigmas b-hats b-dsigmas iter-a iter-b observer-a &key (kink-offset 0) (delta 1e-10) (beyond 1))
  (let* ((first-diamond (create-first-worldsheet-diamond a-hats a-dsigmas b-hats b-dsigmas iter-a iter-b 0.0))
	 (diamond first-diamond)
	 (observer (diamond-position first-diamond :a observer-a :b 0.5))
	 (min-x (* (- (* 2 observer-a) 1.0) (diamond-a-t first-diamond)))
	 (max-x (+ min-x delta))
	 (midline-vector (4vector-normalize-time (diamond-a first-diamond)))
	 (L (diamond-a-t first-diamond))
	 (other-vector (4vector-normalize-time (diamond-b first-diamond)))
	 (a-index iter-a)
	 (b-index iter-b)
	 (a-steps 0)
	 (b-steps 0)
	 (goal-steps (length a-hats))
	 (dv (make-4vector)))
    (flet
     ((test-function (test-diamond)
		     (let ((offset (4vector- observer (diamond-right test-diamond))))
		       (> (4vector-t offset) (3vector-length offset))))
      (do-crossing (test-diamond)
		   (if (diamond-sw test-diamond)
		       (if (diamond-se test-diamond)
		   (past-crossing first-diamond test-diamond midline-vector L other-vector min-x max-x)
			 (b-crossing first-diamond test-diamond midline-vector L other-vector min-x max-x))
		     (if (diamond-se test-diamond)
			 (a-crossing first-diamond test-diamond midline-vector L other-vector min-x max-x)
		       (future-crossing first-diamond test-diamond midline-vector L other-vector min-x max-x)))))
     ;; First, find the kink crossing
     (loop
      (cond ((= b-steps (1+ kink-offset))
	     (4vector-incf dv (future-crossing first-diamond (diamond-w diamond) midline-vector L other-vector min-x max-x))
	     (return nil))
	    ((test-function diamond)
	     (setf b-index (previous-index-wrapping b-index (length b-hats))
		   diamond (build-north-diamond diamond b-hats b-dsigmas b-index)
		   b-steps (1+ b-steps)))
	    (t
	     (setf a-index (previous-index-wrapping a-index (length a-hats))
		   diamond (build-south-diamond diamond a-hats a-dsigmas a-index)
		   a-steps (1+ a-steps)))))
     ;; Then, loop until we have gone :BEYOND times as far beyond the kink as we did before the kink
     (setf goal-steps (floor (* (1+ beyond) (+ a-steps b-steps))))
     (loop
      do (cond ((test-function diamond)
		(setf b-index (next-index-wrapping b-index (length b-hats))
		      diamond (build-north-diamond diamond b-hats b-dsigmas b-index)
		      b-steps (1+ b-steps))
		(4vector-incf dv (do-crossing (diamond-w diamond))))
	       (t
		(setf a-index (previous-index-wrapping a-index (length a-hats))
		      diamond (build-south-diamond diamond a-hats a-dsigmas a-index)
		      a-steps (1+ a-steps))
		(4vector-incf dv (do-crossing (diamond-w diamond)))))
      ;; only terminate once we've included at least as many source diamonds above the kink as below, and we just stepped
      ;; across a line of constant u
      until (and (<= goal-steps (+ a-steps b-steps)) (diamond-sw diamond)))
     ;;(format t "~D ~D ~D~%" a-steps b-steps goal-steps)
     (4vector-scale dv (/ delta)))))

|#
     
;; Prints the TXYZ and WRCD acceleration components an observer located A-(B-)OFFSET a-(b-)segments away from a cusp
;; on the canonical KT loop.
(defun cusp-acceleration (a-offset b-offset &key (alpha 0.1) (phi (/ pi 2)) (L (* 2 pi)) (points 1000) (max-x 1e-9))
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (make-kt-hats alpha phi :points points :L L)
    (let* ((a-index (mod (- (/ points 2) a-offset) points))
	   (b-index (mod (- (/ points 2) b-offset) points))
	   (dv (4vector-scale (a-segment-change-1 a-hats a-dsigmas b-hats b-dsigmas
						  a-index b-index 0.0 :delta max-x)
			      0.5)) ;scale by 1/2 to get acceleration from Delta-A'
	   (w-basis (3to4vector (3vector-normalize (3vector+ (aref b-hats (/ points 2)) (aref b-hats (1+ (/ points 2))))
					       0.5) 0.5))
	   (r-basis (3to4vector (3vector- (4to3vector w-basis)) 0.5))
	   (c-basis) (d-basis))
      (multiple-value-setq (c-basis d-basis) (two-perpendicular-vectors w-basis))
      ;(format t "~S~%" w-basis)
      ;; put everything on one line, but split it in case we later want to enable/disable certain parts
      (format t "~D ~D ~D " points a-offset b-offset)
      (format t "~8$ ~8$ ~8$ " (* 2 (4dot r-basis dv)) (* 2 (4dot w-basis dv))
	      (sqrt (+ (expt (4dot c-basis dv) 2) (expt (4dot d-basis dv) 2))))
      (format t "~8$ ~8$ ~8$ ~8$~%" (4vector-t dv) (4vector-x dv) (4vector-y dv) (4vector-z dv)))))

#| Not updated for changes to bind special variables before calling create-initial-ring

;; Prints the WRCD components of the contributions of individual source diamonds on the u branch of an observer
;; located A-(B-)OFFSET a-(b-)steps away from a cusp on the canonical KT loop. The contributions are printed every
;; STEP-SIZE steps in a or b until we reach the contribution which is RANGE steps in total from the observer.
(defun cusp-acceleration-integrand (a-offset b-offset range &key (alpha 0.1) (phi (/ pi 2)) (L (* 2 pi))
					     (points 1000) (max-x 1e-9) debug)
  (mirror-image-let ((a-index (mod (- (/ points 2) a-offset) points)))
    (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	(make-kt-hats-1 alpha phi :points points :L L
			:start-a a-index :start-b (mod (- b-index range) points) :range range)
      (let* ((diamond (create-initial-ring :a-hats b-hats :a-dsigmas b-dsigmas :b-hats a-hats
					   :b-dsigmas a-dsigmas :start-a (1- range) :start-b 0
					   :observer-a 0.5 :observer-b 0.5 :open-loop t))
	     (w-basis (3to4vector (3vector-normalize (3vector+
						      (funcall (circle-hat-function :phi phi) (- (/ 2) (/ (* 2 points))))
						      (funcall (circle-hat-function :phi phi) (+ (/ 2) (/ (* 2 points)))))
						     0.5) 0.5))
	     (r-basis (3to4vector (3vector- (4to3vector w-basis)) 0.5))
	     (c-basis) (d-basis)
	     (dv (make-4vector)))
	(multiple-value-setq (c-basis d-basis) (two-perpendicular-vectors w-basis))
	(format t "~D ~D ~D~%" points a-offset b-offset)
	(if debug (return-from cusp-acceleration-integrand diamond))
	;; we will collect all corrections from diamonds in a block which always has a shape of the logical negation sign,
	;; \neg in TeX, and then assign that value to the position |b|=sqrt(b1*b2) 
	(loop with count = 0
	      with test-diamond = (find-most-se diamond)
	      with b1 = (number-a-w test-diamond diamond) ;starting b for the first diamond block
	      with b2 = b1				  ;just need some starting value
	      until (<= range count)
	      do (4vector-incf dv (4vector-scale (find-correction-1 diamond test-diamond
								    (4vector-normalize-time (diamond-a diamond))
								    (diamond-a-t diamond)
								    (4vector-normalize-time (diamond-b diamond))
								    :min-x 0.5 :max-x (+ 0.5 max-x))
						 (/ max-x)))
	      do (setf test-diamond (diamond-e test-diamond)
		       count (1+ count))
	      do (if (and (diamond-nw test-diamond) (diamond-ne test-diamond))
		     (progn
		       (setf b2 (number-a-w (diamond-w test-diamond) diamond)	   ;step back to get the ending b
			     dv (4vector-scale dv (/ points (* L (sqrt (* b1 b2)))))) ;u value of inverse of mean
		       (format t "~3$ " (sqrt (* b1 b2)))
		       (format t "~8$ ~8$ ~8$~%" (* 2 (4dot r-basis dv)) (* 2 (4dot w-basis dv)) ; factor 2 from \eta^{wr}
			       (sqrt (+ (expt (4dot c-basis dv) 2) (expt (4dot d-basis dv) 2))))
		       (setf dv (make-4vector) ;clear the dv
			     b1 (1+ b2)))))    ; the new starting b is the b of the current test diamond
	))))

|#


(defun general-acceleration-a (a-index b-index &key (alpha 0.1) (phi (/ pi 2)) (L (* 2 pi)) (points 1000) (max-x 1e-9))
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (make-broken-kt-hats alpha phi :points points :L L)
    (let  ((dv (4vector-scale (a-segment-change-1 a-hats a-dsigmas b-hats b-dsigmas
						 a-index b-index 0.5 :delta max-x) 0.5)))
      (format t "~8$ ~8$ ~8$ ~8$~%" (4vector-t dv) (4vector-x dv) (4vector-y dv) (4vector-z dv)))))
(defun general-acceleration-b (a-index b-index &key (alpha 0.1) (phi (/ pi 2)) (L (* 2 pi)) (points 1000) (max-x 1e-9))
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (make-broken-kt-hats alpha phi :points points :L L)
    (let  ((dv (4vector-scale (b-segment-change-1 a-hats a-dsigmas b-hats b-dsigmas
						 a-index b-index 0.5 :delta max-x) 0.5)))
      (format t "~8$ ~8$ ~8$ ~8$~%" (4vector-t dv) (4vector-x dv) (4vector-y dv) (4vector-z dv)))))


(defvar *saved-random-state*)
(defun save-random-state ()
  (setq *saved-random-state* (make-random-state *random-state*))
  t)
(defun restore-random-state ()
  (setq *random-state* (make-random-state *saved-random-state*))
  t)


