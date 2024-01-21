;;;Solve equations for diamond intersection
(in-package "CL-USER")

;;Index vectors
(defvar 4indices
 (make-resource :name "4index" :constructor #'(lambda () (make-array 4 :element-type '(unsigned-byte 16))) :max 10))

(defvar *check-projected-overlap* t)	;If set, check some 4vector projections before solving for intersection

(defun handle-possible-intersection (d1 d2)
  (when (and (not (eq (diamond-e d1) :BH))
	     (not (eq (diamond-w d1) :BH))
	     (not (eq (diamond-e d2) :BH))
	     (not (eq (diamond-w d2) :BH))
	     (not (diamond-bh d1))
	     (not (diamond-bh d2)) ;intersections in bh diamonds are excluded
	     (diamond-box-overlap d1 d2) ;Check quickly that the diamond boxes overlap
	     (or (not *check-projected-overlap*) (diamond-projected-overlap d1 d2)) ;maybe try some projections also
	     )
    (if (eq *era* :flat)
	(handle-possible-intersection-flat d1 d2)
      (handle-possible-intersection-curved d1 d2))))

;;If these diamonds intersect, their projections into any one-dimensional subspace must overlap.
(defun diamond-projected-overlap (d1 d2)
  (declare (optimize speed))
  (macrolet ((try (direction-form)
	       `(let ((direction ,direction-form))
		  (prog1 (diamond-projected-overlap-1 d1 d2 direction) (deallocate 4vectors direction)))))
    (and (try (diamond-a d1))
	 (try (diamond-b d1))
	 (try (diamond-a d2))
	 (try (diamond-b d2)))))


;;Minimum and maximum of spatial inner product between position of any point on the diamond and the given
;;direction vector.
;;This is the 4vector inner products, so it ignores displacements along the direction vector and considers only those
;;in other spatial directions to and the opposite spacetime direction.
;;These are macros instead of inline functions, because otherwise there seems to be a problem with consing
;;arguments to MIN and MAX.
(defmacro diamond-projected-minimum (diamond direction)
  `(min (the double-float (4vector-dot-spacelike (diamond-start ,diamond) ,direction))
	(the double-float (4vector-dot-spacelike (diamond-left ,diamond) ,direction))
	(the double-float (4vector-dot-spacelike (diamond-right ,diamond) ,direction))
	(the double-float (4vector-dot-spacelike (diamond-end ,diamond) ,direction))))

(defmacro diamond-projected-maximum (diamond direction)
  `(max (4vector-dot-spacelike (diamond-start ,diamond) ,direction)
	(4vector-dot-spacelike (diamond-left ,diamond) ,direction)
	(4vector-dot-spacelike (diamond-right ,diamond) ,direction)
	(4vector-dot-spacelike (diamond-end ,diamond) ,direction)))


(defun diamond-projected-overlap-1 (d1 d2 direction)
  (declare (optimize speed))
  (let ((min1 (diamond-projected-minimum d1 direction))
	(min2 (diamond-projected-minimum d2 direction))
	(max1 (diamond-projected-maximum d1 direction))
	(max2 (diamond-projected-maximum d2 direction)))
    (and (>= max2 min1) (>= max1 min2)))) ;overlap unless one min is > the other max.

;;See if diamonds intersect.  Return INTERSECTION structure or NIL.
(defun handle-possible-intersection-flat (d1 d2)
  (let* ((b-vector (4vector- (diamond-start d2) (diamond-start d1))) ;Vector from start d1 to start d2
	 (m (fill-up-coeff-matrix d1 d2)) ;Create the matrix with the p and q as columns
	 (lu (copy-matrix-4-4 m))  ;Make a copy of M to give to LUDCMP
	 (index (allocate 4indices))	;Permutation vector
	 (vv (allocate 4vectors)))	;Working storage for ludcmp
    (prog1 (and (ludcmp lu index vv) ;LU-decompose: result in LU, permutation in INDEX.  NIL if singular
		(find-the-intersection lu m index b-vector d1 d2)) ;Find the actual intersection.  Return it or NIL
      (deallocate 4vectors b-vector vv)
      (deallocate 4indices index)
      (deallocate 4x4 m lu))))

(declaim (inline fill-up-Jacobian))	;Avoid consing of floats

;;Set up the matrix M to be the Jacobian giving the derivatives of the position in d2 minus the position in d1
;;as functions of the variables a1, b1, a2, b2.
;;offset is 0 or 2 to determine which diamond we're doing.  Sign is +1 for d2, -1 for d1.
(without-compiler-notes ;Avoid notes about not optimizing check that type is correct
(defun fill-up-Jacobian (m offset sign d a b)
  (declare (double-float sign a b)
	   (optimize speed)
	   (type (simple-array double-float (4 4)) m))
  (let ((start (diamond-start d))
	(left (diamond-left d))
	(right (diamond-right d))
	(end (diamond-end d)))
    (declare (type 4vector start left right end))
    (dotimes (i 4)
      (let* ((p (- (aref left i) (aref start i))) ;Component of P
	     (q (- (aref right i) (aref start i))) ;Component of Q
	     (new-p  (- (aref end i) (aref right i))) ;Component of P'
	     (quadratic (- new-p p)))	;Coefficient of quadratic term
	(setf (aref m i offset) (* sign (+ p (* quadratic b)))
	      (aref m i (1+ offset)) (* sign (+ q (* quadratic a))))))))
)

(defvar max-Newton-iterations 6)
(declaim (fixnum max-Newton-iterations))

(defvar *intersections-checked* 0)
(defvar *intersections-found* 0)
(defvar *worse-intersections-found* 0)
(defvar *max-iterations-needed* 0)


(declaim (inline solution-error))
;;Find error of trial solution, store into 4vector error
(defun solution-error (d1 d2 solution error)
  ;;It would be better to have the results stored automatically as part of the vector manipulation system
  (let ((new (4vector- (diamond-position d2 :a (vector-parameters-a2 solution) :b (vector-parameters-b2 solution))
		       (diamond-position d1 :a (vector-parameters-a1 solution) :b (vector-parameters-b1 solution)))))
    (copy-4vector new error)
    (deallocate 4vectors new)))


;;Curved space version.  Does not polish with 4-d Newton's method, since that makes it worse about 1/3 of the time
(defun handle-possible-intersection-curved (d1 d2)
  (incf *intersections-checked*)
  (catch 'solution-failed
    (let ((solution (kite-solution d1 d2)))
      (when (and solution (check-solution-inside-diamond solution d1 d2)) ;If it really is inside
	(incf *intersections-found*)
	(build-intersection-structure solution d1 d2)))))


;; Old method, using 4-d Newon's method.
(defun old-handle-possible-intersection-curved (d1 d2)
  (incf *intersections-checked*)
  (using-resources ((solution 4vectors) ;Trial solution
		    (error 4vectors))	;coordinate error in solution
    (declare (type (or null 4vector) solution))
    (dotimes (i 4) (setf (aref solution i) 0.5))	;Initial guess is center of diamonds
    (solution-error d1 d2 solution error) ;Compute error of initial guess
    (loop for last-error double-float = 0.0 then this-error
	  for this-error double-float = (4vector-Euclidian-length error) ;Length of error in result
	  for count from 0
	  with worse = nil
;;	  do (format t "~&Guess ~D: ~S: Error ~S" count solution this-error)
	  when (< this-error solution-accuracy) ;Desired accuracy reached
	    do 
	    (case (one-newton-step d1 d2 solution error) ;Try to improve solution once more
		  (:singular (warn "Final error improvement failed.  Giving up.") ;singular on last go
			     (return nil))
		  ((t)			;Changed solution
		   ;;		(format t "~&Final solution ~S, error ~S" solution (4vector-Euclidian-length error))
		   (when (> (4vector-Euclidian-length error) solution-accuracy) ;If still a good solution
		     (warn "Error got worse in last step from ~S to ~S.  Giving up." this-error (4vector-Euclidian-length error))
		     (return nil))))	;No action when solution is unchanged or actually improved
	    (setq *max-iterations-needed* (max *max-iterations-needed* count))
	    (when worse (incf *worse-intersections-found*))
	    (return (and (check-solution-inside-diamond solution d1 d2) ;If it really is inside
			 (prog1 (build-intersection-structure solution d1 d2) ;return it.
			   (setq solution nil)))) ;Clear variable, so that vector is not recycled
	  when (and (plusp count) (> this-error last-error))
	    do (setq worse t)
	    and do (format t "worse")
	  until (> count max-Newton-iterations)	;Converging too slowly.
	  while (eq (one-Newton-step d1 d2 solution error) t) ;take one step and loop unless it fails or does not change the result
	  )))

;;This is only called by old-handle-possible-intersection-curved, and so is obsolete.
;;Take one step of Newton's method.  SOLUTION is the trial solution, which is updated in place to the new trial.
;;ERROR is the error of the original solution.  We update it to the error of the new solution.
;;We limit the solution to parameters 0...1.0.
;;Returns T if successful, :SINGULAR if the Jacobian is singular,
;;NIL if there was no change to the result because of being stuck against the edge of valid parameter space.
(defun one-Newton-step (d1 d2 solution error)
  (declare (type 4vector solution error)
	   (optimize speed))
  (using-resources ((j 4x4)		;Matrix of equations.
		    (lu 4x4)		;LU matrix
		    (index 4indices)	;Permutation vector
		    (vv 4vectors))	;Working storage for ludcmp
    (without-compiler-notes ;Avoid notes about not optimizing check that type is correct
     (fill-up-Jacobian j 0 -1.0 d1 (vector-parameters-a1 solution) (vector-parameters-b1 solution)) ;Compute Jacobian in J
     (fill-up-Jacobian j 2 1.0 d2 (vector-parameters-a2 solution) (vector-parameters-b2 solution)))
    (copy-matrix-4-4 j lu)		;Copy Jacobian to LU.
    (cond ((ludcmp lu index vv)	;Compute LU matrix in place. NIL if Jacobian is singular.
	   (lubksb lu index error) ;Convert ERROR from desired change in result to negative change in solution
	   (let ((factor 1.0)) ;Amount of proposed change to apply to stay within bounds
	     (dotimes (i 4) ;First make sure proposed changes not too large
	       (let* ((old (aref solution i))
		      (change (aref error i))
		      (new (- old (* factor change)))) ;New value with present factor
		 (cond ((minusp new)	;Oops: negative
			(setq factor (/ old change))) ;decrease factor so result is zero
		       ((> new 1.0)	;Oops: too big
			(setq factor (/ (- old 1.0) change))) ;decrease factor so result is 1.0
		       )))
;;	     (format t "factor ~S ~%" factor)
	     (cond ((plusp factor) ;If we are stuck against the edge, return NIL
		    (dotimes (i 4)	 ;Otherwise update solution to new try
		      (setf (aref solution i)
			    (max 0.0 (min 1.0	;We could still get out of range by floating point error
					  (- (aref solution i) (* factor (aref error i)))))))
		    (solution-error d1 d2 solution error) ;Compute new error
		    t)			;Indicate success
		   )))	
	  (t :singular)				;Singular Jacobian
	  )))

;; Create a matrix whose columns are the p and q vectors of the d1 and -p and -q of d2.
(defun fill-up-coeff-matrix (d1 d2)
  (declare (optimize speed)
	   (optimize (safety 0)))	;Do not check array rank and size
  (let ((m (allocate 4x4))
	(start1 (diamond-start d1))
	(left1 (diamond-left d1))
	(right1 (diamond-right d1))
	(start2 (diamond-start d2))
	(left2 (diamond-left d2))
	(right2 (diamond-right d2)))
    (declare (type 4vector start1 left1 right1 start2 left2 right2)
	     (type (simple-array double-float (4 4)) m))
    (dotimes (i 4)
      (setf (aref m i 0) (- (aref left1 i) (aref start1 i)) ;p1
	    (aref m i 1) (- (aref right1 i) (aref start1 i)) ;q1
	    (aref m i 2) (- (aref start2 i) (aref left2 i)) ;-p2
	    (aref m i 3) (- (aref start2 i) (aref right2 i)) ;-q2
	    ))
    m))

;;Computes the solution of the possible intersection and its position in spacetime.
;;Returns INTERSECTION structure.
(defun find-the-intersection (lu m index b d1 d2)
  (declare (type 4vector b)
	   (type (simple-array double-float (4 4)) m))
  (let ((solution (copy-4vector b))) ;Initialize SOLUTION with right hand side
    (declare (optimize speed)
	     (type 4vector solution))
    (lubksb lu index solution)	     ;Substitute.  Result in SOLUTION
    (when (check-solution-roughly-inside-diamond solution) ;Make sure roughly right.  Could fail if nearly singular or no intersection.
      (using-resource (error 4vectors)
	(declare (type 4vector error))
	(flet ((compute-error () ;Compute distance between points that supposedly touch
		 (dotimes (i 4)
		   (setf (aref error i) ;ERROR =  B - M . SOLUTION
			 (- (aref b i)
			    (loop for j below 4
				  sum (* (aref m i j) (aref solution j)) double-float))))))
	  (declare (inline compute-error))
	  (compute-error)	      ;Get error from initial solution
	  (lubksb lu index error) ;Substitute.  ERROR now holds error in SOLUTION
	  (dotimes (i 4)	 ;Adjust SOLUTION to improve error
	    (incf (aref solution i) (aref error i)))
	  (compute-error)		;Compute error again
	  (when (and (< (4vector-Euclidian-length error) solution-accuracy) ;If good enough, we found intersection.  Could be large if singular.
		     (check-solution-inside-diamond solution d1 d2))
	    (build-intersection-structure solution d1 d2)))))))

;; It looks at the values of the parameters to see if there is any hope
;; that the solution could be within the diamonds
(defun check-solution-roughly-inside-diamond (solution)
  (declare (optimize speed)
	   (type 4vector solution))
  (loop for index below 4
	always (< -1.0d0 (aref solution index) 2.0d0)))

;; It looks at the values of the parameters that parametrize the diamond
;; that are always between 0 and 1. We include a fudge factor to allow for
;; possible numerical error.
(defun check-solution-inside-diamond (solution d1 d2)
  ;; we don't need to be incredibly precise about when the intersection happens, so we use the globalized time of
  ;; the solution parameters for d1. d2 should work just as well.
  (let* ((intersection-time (global-time (diamond-position-time d1 :a (vector-parameters-a1 solution) :b (vector-parameters-b1 solution))))
	 (reference-size (* minimum-ending-diamond-width (scale-factor-ratio *overall-end-time* intersection-time)))
	 (physicala1 (/ reference-size (4vector-t (diamond-p d1))))
	 (physicalb1 (/ reference-size (4vector-t (diamond-q d1))))
	 (physicala2 (/ reference-size (4vector-t (diamond-p d2))))
	 (physicalb2 (/ reference-size (4vector-t (diamond-q d2))))
	 (physical-fudge-factor (make-4vector physicala1 physicalb1 physicala2 physicalb2)))
    (loop for index below 4
	  always (< (aref physical-fudge-factor index) (aref solution index) (- 1.0d0 (aref physical-fudge-factor index))))))

;;Creates an intersection structure for the obtained solution
(defun build-intersection-structure (solution-parameters d1 d2)
  (let ((spacetime (solution-to-spacetime-position solution-parameters d1 d2)))
    (make-intersection
     :diamond-1 d1
     :diamond-2 d2
     :vector-parameters solution-parameters
     :spacetime spacetime)))



;; Returns the relative sign of intersection for the diamonds
;; d1 and d2.  This uses the same sign conventino as the old-flat-intersection-number.
(defun intersection-number (d1 d2)
  (let* ((segment (segment-in-common d1 d2))
	 (lambdaA (first segment))
	 (lambdaB (second segment))
	 (dtA (and segment (delta-t-of-lambda d1 d2 lambdaA)))
	 (dtB (and segment (delta-t-of-lambda d1 d2 lambdaB))))
    (cond ((null segment)
	   0)
	  ((or (eq (plusp dta) (plusp dtb)) ;both non-positive or both non-negative
	       (< (abs dta) fudge-coordinates) ;neither too close to zero
	       (< (abs dtb) fudge-coordinates))
	   0)
	  ((> dta dtb)			;negative slope at root
	   1)
	  (t				;positive slope at root
	   -1))))


;; Returns the (a1 b1 a2 b2) location of an intersection if it exists, else nil
;; Uses Newton's method to make guesses, but when this falls outside 
;; of the interval, the midpoint is used, keeping the root always bracketed.
(without-compiler-notes
(defun kite-solution (d1 d2)
  (declare (optimize speed)
           (diamond d1 d2))
  (let* ((segment (or (segment-in-common d1 d2) (return-from kite-solution nil)))
	 (lambdaA (first segment))
	 (lambdaB (second segment))
	 (dtA (delta-t-of-lambda d1 d2 lambdaA))
	 (dtB (delta-t-of-lambda d1 d2 lambdaB))
	 (best-lambda (and (not (eq (plusp dtA) (plusp dtB))) ;require transversal intersection
			   (> (abs dtA) fudge-coordinates) ;don't consider intersections at edge of a diamond
			   (> (abs dtB) fudge-coordinates)
			   (loop for count below 40 ;almost all are found when count is below 3. 
				 with lambda-min = lambdaA
				 with lambda-max = lambdaB
				 with dt-min = dtA
				 with dt-max = dtB
				 with new-dt
				 for old-lambda = lambdaA then new-lambda
				 for old-dt = dtA then new-dt
				 for old-dt-prime = (delta-t-prime d1 d2 old-lambda)
				 for new-lambda = (- old-lambda (/ old-dt old-dt-prime))
				 unless (and (<= lambda-min new-lambda lambda-max) ;keep the root bracketed and
					     (eq (plusp old-dt-prime) (plusp (- dtB dtA))) ;if 3 intersections, bracket a good one.
					     (< count 10)) ;if Newton's method is failing, bracketing will certainly work.  
				   do (setf new-lambda (/ (+ lambda-min lambda-max) 2.0))
				 do (setq new-dt (delta-t-of-lambda d1 d2 new-lambda))
				 when (< (abs new-dt) solution-accuracy)                                ;when goal is reached,
				   return (- new-lambda (/ new-dt (delta-t-prime d1 d2 new-lambda))) ;do one last step
				 when (eq (plusp new-dt) (plusp dt-min))
				   do (setq dt-min new-dt)
				      (setq lambda-min new-lambda)
				 when (eq (plusp new-dt) (plusp dt-max))
				   do (setq dt-max new-dt)
				      (setq lambda-max new-lambda)
				 finally (return new-lambda) ;is it worth checking it anyway?
				 )))
	 (ab1 (if best-lambda (a1-b1-of-lambda d1 d2 best-lambda) (return-from kite-solution nil)))
	 (ab2 (a2-b2-of-lambda d1 d2 best-lambda))
	 (mismatch (4vector-euclidian-distance (diamond-position d2 :a (first ab2) :b (second ab2))
					       (diamond-position d1 :a (first ab1) :b (second ab1)))))
    (if (< mismatch solution-accuracy) (list-4vector (nconc ab1 ab2))))) ;only return solutions that are within solution-accuracy
)

;; time difference of string crossing point
;; lambda of segment-in-common.
(defun delta-t-of-lambda (d1 d2 lambda)
  (- (t2-of-lambda d1 d2 lambda)
     (t1-of-lambda d1 d2 lambda)))

;; Derivative of delta-t-of-lambda
(defun delta-t-prime (d1 d2 lambda)
  (- (t2-prime d1 d2 lambda)
     (t1-prime d1 d2 lambda)))



;; The derivative of t1-of-lambda
(without-compiler-notes
(defun t1-prime (d1 d2 lambda)
  (declare (optimize speed)
	   (diamond d1 d2)
	   (double-float lambda))
  (let* ((p-t (4vector-t (diamond-p d1)))
	 (q-t (4vector-t (diamond-q d1)))
	 (del-t (- (4vector-t (diamond-new-p d1)) p-t))
	 (abdadb (abdadb1-of-lambda d1 d2 lambda))
	 (a (first abdadb))
	 (b (second abdadb))
	 (da (third abdadb))
	 (db (fourth abdadb)))
    (+ (* p-t da) 
       (* q-t db) 
       (* del-t (+ (* da b) (* a db)))))))

(defun t2-prime (d1 d2 lambda)
  (- (t1-prime d2 d1 (- lambda))))


;; time of d1 at displacement lambda along segment-in-common
(defun t1-of-lambda (d1 d2 lambda)
  (4vector-t (x1-of-lambda d1 d2 lambda)))

(defun t2-of-lambda (d1 d2 lambda)
  (t1-of-lambda d2 d1 (- lambda)))

;; spacetime position in d1 of displacement lambda along segment in common
(defun x1-of-lambda (d1 d2 lambda)
 ; (declare (optimize speed)
	;   (diamond d1 d2)
	;   (double-float lambda))
  (let* ((ab (a1-b1-of-lambda d1 d2 lambda))
	 (a (nth 0 ab))
	 (b (nth 1 ab)))
    (diamond-position d1 :a a :b b)))

(defun x2-of-lambda (d1 d2 lambda)
  (x1-of-lambda d2 d1 (- lambda)))

(without-compiler-notes
;;(a,b) location in d1 of displacement lambda along segment in common
(defun a1-b1-of-lambda (d1 d2 lambda)
  (declare (optimize speed)
	   (diamond d1 d2)
	   (double-float lambda))
  (let* ((alpha (a-alphabeta1-of-lambda d1 d2 lambda))
	 (beta (b-alphabeta1-of-lambda d1 d2 lambda))
	 (p (diamond-p d1))
	 (q (diamond-q d1))
	 (delta (4vector- (diamond-new-p d1) p))
	 (pxq (3vector-cross-product p q))
	 (pperp (3vector-cross-product p pxq))
	 (qperp (3vector-cross-product pxq q))
	 (pqperp (3vector-dot p qperp))
	 (delpperp (3vector-dot delta pperp))
	 (delqperp (3vector-dot delta qperp))
	 (termb (+ pqperp (* alpha delpperp) (* beta (- delqperp))))
	 (terma (- (* 2.0 pqperp) termb))
	 (root2 (+ (* 4.0 pqperp delqperp beta)
		   (expt termb 2.0)))
	 (root (cond ((>= root2 0.0) (real-sqrt root2))
		     (t
		      (warn "Solution failed because (a,b) not even close to (alpha,beta).  Diamond sizes ~S, ~S, ~S, ~S"
			    (4vector-t (diamond-a d1)) (4vector-t (diamond-b d1))
			    (4vector-t (diamond-a d2)) (4vector-t (diamond-b d2)))
		      (throw 'solution-failed nil))))
	 (a (if (zerop (- terma root)) 0.0 ;a = alpha = 0
	      (/ (* 2.0 pqperp alpha)
		 (- terma root))))
	 (b (if (zerop (- termb root)) 0.0 ;b = beta = 0
	      (/ (* 2.0 pqperp beta)
		 (- termb root)))))
    (deallocate 4vectors p q delta)
    (deallocate 3vectors pxq pperp qperp)
    (if (< root2 0.0) (or (format t "~%a,b = ~D ~D" a b) (error "asdf")))
    (list a b)))
)

(defun a2-b2-of-lambda (d1 d2 lambda)
  (a1-b1-of-lambda d2 d1 (- lambda)))

;;In d1, the values (a,b) and their derivatives w.r.t. lambda
;;Output is '(a, b, a', b').  This makes a1-b1-of-lambda obsolete, unless speed is important
(without-compiler-notes
(defun abdadb1-of-lambda (d1 d2 lambda)
  (declare (optimize speed)
	   (diamond d1 d2)
	   (double-float lambda))
  (let* ((alpha (a-alphabeta1-of-lambda d1 d2 lambda))
	 (beta (b-alphabeta1-of-lambda d1 d2 lambda))
	 (dalpha (d-a-alphabeta1-of-lambda d1 d2))
	 (dbeta (d-b-alphabeta1-of-lambda d1 d2))
	 (p (diamond-p d1))
	 (q (diamond-q d1))
	 (delta (4vector- (diamond-new-p d1) p))
	 (pxq (3vector-cross-product p q))
	 (pperp (3vector-cross-product p pxq))
	 (qperp (3vector-cross-product pxq q))
	 (pqperp (3vector-dot p qperp))
	 (delpperp (3vector-dot delta pperp))
	 (delqperp (3vector-dot delta qperp))
	 (termb (+ pqperp (* alpha delpperp) (* beta (- delqperp))))
	 (terma (- (* 2.0 pqperp) termb))
	 (root (real-sqrt (+ (* 4.0 pqperp delqperp beta)
			(expt termb 2.0))))
	 (a (if (zerop (- terma root)) 0.0
	      (/ (* 2.0 pqperp alpha)
		 (- terma root))))
	 (b (if (zerop (- termb root)) 0.0
	      (/ (* 2.0 pqperp beta)
		 (- termb root))))
	 (dadalpha (if (zerop b) 1.0
		     (/ (- (* pqperp beta))
			b root)))
	 (dadbeta (/ (* a delqperp) root))
	 (dbdalpha (/ (* b delpperp) root))
	 (dbdbeta (if (zerop a) 1.0
		    (/ (- (* pqperp alpha))
		       a root)))
	 (da (+ (* dadalpha dalpha) (* dadbeta dbeta)))
	 (db (+ (* dbdalpha dalpha) (* dbdbeta dbeta))))
    (deallocate 4vectors p q delta)
    (deallocate 3vectors pxq pperp qperp) 
    (list a b da db))))

;; Project diamonds to 3-space and see if they overlap.
;; Kites are flat planes with straight edges in R^3
;; lambda is the direction parallel to both kites
;; Returns the interval of overlap '(lambdaA,lambdaB), or nil if it is empty.
(defun segment-in-common (d1 d2)
;  (declare (optimize speed))
  (let* ((Lambda1 (segment1 d1 d2))
	 (Lambda2 (segment2 d1 d2))
	 (lambdaA (and Lambda1 Lambda2 (max (first Lambda1) (first Lambda2))))
	 (lambdaB (and Lambda1 Lambda2 (min (second Lambda1) (second Lambda2)))))
    (and lambdaA lambdaB (< lambdaA lambdaB) (list lambdaA lambdaB))))

;; segment-in-common is the intersection of segment1 and segment2
(without-compiler-notes
(defun segment1 (d1 d2)
  (declare (optimize speed))
  (let* ((p1 (diamond-p d1))
	 (q1 (diamond-q d1))
	 (s1 (diamond-start d1))
	 (delta1 (4vector- (diamond-new-p d1) p1))
	 (pxq1 (3vector-cross-product p1 q1))
	 (p1perp (3vector-cross-product p1 pxq1))
	 (q1perp (3vector-cross-product pxq1 q1))
	 (alpha1end (/ (3vector-dot q1perp (3vector+ delta1 p1)) (3vector-dot q1perp p1)))
	 (beta1end (/ (3vector-dot p1perp (3vector+ delta1 q1)) (3vector-dot p1perp q1)))
	 (p2 (diamond-p d2))
	 (q2 (diamond-q d2))
	 (pxq2 (3vector-cross-product p2 q2))
	 (lambda-raw-vector (3vector-cross-product (3vector-normalize pxq1) (3vector-normalize pxq2))) ;normalize to get angle
	 (lambda-raw-length (3vector-length lambda-raw-vector))	;sin(angle) >= 0
	 (lambda-hat (if (> lambda-raw-length solution-accuracy) ;how parallel are d1,d2?
			 (3vector-scale lambda-raw-vector (/ 1.0 lambda-raw-length))
		       (return-from segment1 nil)))
	 (ds (3vector- (diamond-start d2) s1))
	 (dspxq2 (3vector-dot ds pxq2))
	 (q1pxq2 (3vector-dot q1 pxq2))
	 (p1pxq2 (3vector-dot p1 pxq2))
	 (lambda1-se-vec (3vector+ s1 (3vector-scale q1 (/ dspxq2 (if (zerop q1pxq2) (return-from segment1 nil)
								    q1pxq2)))))
	 (lambda1-se (3vector-dot lambda-hat lambda1-se-vec))
	 (lambda1-sw-vec (3vector+ s1 (3vector-scale p1 (/ dspxq2 (if (zerop p1pxq2) (return-from segment1 nil)
								    p1pxq2)))))
	 (lambda1-sw (3vector-dot lambda-hat lambda1-sw-vec))
	 (denom-nw-1 (+ (* (1- alpha1end) p1pxq2) (* beta1end q1pxq2)))
	 (lambda1-nw-vec (3vector+ s1 
				   (3vector-scale p1 (/ (+ (* beta1end q1pxq2) (* (1- alpha1end) dspxq2)) 
							(if (zerop denom-nw-1) (return-from segment1 nil)
							  denom-nw-1)))
				   (3vector-scale q1 (/ (* beta1end (- dspxq2 p1pxq2))
							denom-nw-1))))
	 (lambda1-nw (3vector-dot lambda-hat lambda1-nw-vec))
	 (denom-ne-1 (+ (* (1- beta1end) q1pxq2) (* alpha1end p1pxq2)))
	 (lambda1-ne-vec (3vector+ s1 
				   (3vector-scale q1 (/ (+ (* alpha1end p1pxq2) (* (1- beta1end) dspxq2)) 
							(if (zerop denom-ne-1) (return-from segment1 nil)
							  denom-ne-1)))
				   (3vector-scale p1 (/ (* alpha1end (- dspxq2 q1pxq2))
							denom-ne-1))))
	 (lambda1-ne (3vector-dot lambda-hat lambda1-ne-vec))
	 (sw-1-test (<= 0.0 (a-alphabeta1-of-lambda d1 d2 lambda1-sw) 1.0))
	 (ne-1-test (<= 0.0 (a-alphabeta1-of-lambda d1 d2 lambda1-ne) alpha1end))
	 (se-1-test (<= 0.0 (b-alphabeta1-of-lambda d1 d2 lambda1-se) 1.0))
	 (nw-1-test (<= 0.0 (b-alphabeta1-of-lambda d1 d2 lambda1-nw) beta1end))
	 (lambda1-list nil)
	 (Lambda1 nil)
	 )
    (mirror-images
    (if sw-1-test (push lambda1-sw lambda1-list))
    (if nw-1-test (push lambda1-nw lambda1-list))
    );mirror images
    (setq Lambda1 (sort lambda1-list #'<))
    (unless (= (length Lambda1) 2) 	;might miss possibility of segment entering diamond precisely at a corner.  rewrite?
      (return-from segment1 nil))
    (unless (> (- (second Lambda1) (first Lambda1)) fudge-coordinates) 
      (return-from segment1 nil)) ;segment merely cuts off a tiny a corner
    (deallocate 4vectors p1 q1 delta1 p2 q2)	;Why can't I put s1 and s2 in this list?
    (deallocate 3vectors pxq1 p1perp pxq2 lambda-raw-vector lambda-hat ds lambda1-se-vec lambda1-sw-vec lambda1-nw-vec lambda1-ne-vec)
    Lambda1
    )))

(defun segment2 (d1 d2)
  (let ((segment (segment1 d2 d1)))
    (and segment (list (- (second segment))
		       (- (first segment))))))



(mirror-images
;; Calculates alpha1 (or beta1) as a function of
;; lambda for two intersecting kites.
;; Caution: lambda-hat flips for d1 <-> d2
;; so this function can be used for alpha2 provided
;; the arguments are given as (d2 d1 -lambda) 
(without-compiler-notes
(defun a-alphabeta1-of-lambda (d1 d2 lambda) 
  (declare (optimize speed)
	   (diamond d1 d2)
	   (double-float lambda))
  (let* ((p1 (diamond-a d1))
	 (q1 (diamond-b d1))
	 (s1 (diamond-start d1))
	 (pxq1 (3vector-cross-product p1 q1))
	 (pxq2 (3vector-cross-product (diamond-a d2) (diamond-b d2)))
	 (lambda-hat (3vector-normalize (3vector-cross-product pxq1 pxq2)))
	 (ds (3vector- (diamond-start d2) s1))
	 (num-vec (3vector- (3vector-scale q1 lambda)
			    (3vector-scale q1 (3vector-dot s1 lambda-hat))
			    (3vector-scale ds (3vector-dot q1 lambda-hat))))
	 (denom-vec (3vector- (3vector-scale q1 (3vector-dot p1 lambda-hat))
			      (3vector-scale p1 (3vector-dot q1 lambda-hat)))))
    (deallocate 4vectors p1 q1)	;Why can't I put s1 in this list?
    (deallocate 3vectors pxq1 lambda-hat ds)
    (/ (3vector-dot num-vec pxq2)
       (3vector-dot denom-vec pxq2)))))

;;The derivative of a-alphabeta1-of-lambda
(defun d-a-alphabeta1-of-lambda (d1 d2) 
;  (declare (optimize speed)
;	   (diamond d1 d2))
  (let* ((p1 (diamond-a d1))
	 (q1 (diamond-b d1))
	 (pxq1 (3vector-cross-product p1 q1))
	 (pxq2 (3vector-cross-product (diamond-a d2) (diamond-b d2)))
	 (lambda-hat (3vector-normalize (3vector-cross-product pxq1 pxq2)))
	 (num-vec q1)
	 (denom-vec (3vector- (3vector-scale q1 (3vector-dot p1 lambda-hat))
			      (3vector-scale p1 (3vector-dot q1 lambda-hat))))
	 )
    (deallocate 4vectors p1 q1)
    (deallocate 3vectors pxq1 lambda-hat)
    (/ (3vector-dot num-vec pxq2)
       (3vector-dot denom-vec pxq2))))

(without-compiler-notes
(defun a-alphabeta2-of-lambda (d1 d2 lambda)
  (declare (optimize speed)
	   (diamond d1 d2)
	   (double-float lambda))
  (a-alphabeta1-of-lambda d2 d1 (- lambda))))

(defun d-a-alphabeta2-of-lambda (d1 d2)
  (- (d-a-alphabeta1-of-lambda d2 d1)))
) ;mirror-images
