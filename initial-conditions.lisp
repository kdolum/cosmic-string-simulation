;;;Vachaspati-Vilenkin initial conditions for string simulation
;;Developed from field.lisp of the lattice simulation code by Jose and Ken

(in-package "CL-USER")

#| Data structures:

SITE -- a fixnum with 3 bytes giving X, Y, Z coordinates

AXIS -- a number giving a coordinate axis.  1 means X, 2 means Y, and so on.

DIRECTION -- a direction along the links of the lattice.  Positive
numbers mean to move along that axis in the positive direction.
Negative numbers mean to move along the same axis in the opposite
direction.  0 is not a valid direction.

LATTICE -- a simple array of lattice sites with a value on
each site.

FACES -- a quantity stored on each face of the lattice.  An array of arrays
in the order (X Y) (X Z) (Y Z)....  Thus element (X Y Z ...) of the first 
array represents the value on the face centered around
(X+1/2, Y+1/2, Z ...).

|#

(defvar *vv-p* nil)			;Set if using VV initialization

;;Constant vv-index-size and types vv-index, site, axis in direction are in in definitions.lisp

;;Extract coordinate from SITE object.  AXIS = 1..3.
(defmacro site-ref (axis site)
  `(ldb (byte vv-index-size (* (1- ,axis) vv-index-size)) ,site))

(defvar *vv-size*)			;Size of the (square) local lattice
(declaim (type vv-index *vv-size*))

;;Here are the lattices that we will work with
(defvar *vv-psi*)			;VV phase array
(defvar *vv-faces*)			;data about faces of the array
(defvar *vv-status*)			;Information about who needs this cube

(deftype vv-status-type () '(unsigned-byte 8))
(declaim (type (simple-array vv-status-type) *vv-status*))

(declaim (inline make-site site-index))

(declaim (ftype (function (site t &optional t) (values site t)) standardize-vv))

;;Return an uninitialized site object
(defun make-site () 0)

(defun site-list (site)
  (loop for axis from 1 to 3 collect (site-ref axis site)))

(defun list-site (list)
  (let ((site (make-site)))
    (loop for axis from 1 to 3
	  for coordinate in list
	  do (setf (site-ref axis site) coordinate))
    site))

;;Compute row-major index into a lattice array
(defun site-index (site)
  (declare (type site site)
	   (optimize speed))
  (+ (* *vv-size* (+ (* *vv-size* (site-ref 1 site)) (site-ref 2 site))) (site-ref 3 site)))

(defmacro lattice-ref (array site)
  `(aref ,array (site-index ,site)))

;;Array of face arrays.
(deftype faces ()
  `(simple-array t (3)))

;;Create a FACES object.
(defun make-faces (total-size &rest options)
  (make-array 3 :initial-contents
	      (loop repeat 3
		    collect (apply #'make-array total-size options))))

;;Iterate over lattice sites, binding variable SITE
;;LEFT-AXES is a list of directions in which we don't use the right boundary point
(defmacro do-lattice ((site &rest left-axes) &body body)
  (let ((maxima (gensym)))
    `(let ((,site (make-site))
	   (,maxima (without-compiler-notes ;Suppress coercion warning from internal byte-assembling code
		     (make-array 3 :element-type 'vv-index :initial-element *vv-size*)))) ;Iteration limits
       (declare (type site site))
       ,@(loop for axis in left-axes	;left-axes is list of forms for axes 1..3.
	       collect `(decf (aref ,maxima (1- ,axis)))) ;If axis in list, reduce maximum
       (loop for x from 0 below (aref ,maxima 0)
	     do (setf (site-ref 1 ,site) x)
	     do (loop for y below (aref ,maxima 1)
		      do (setf (site-ref 2 ,site) y)
		      do (loop for z below (aref ,maxima 2)
			       do (setf (site-ref 3 ,site) z)
			       do (progn ,@body)))))))

(defmacro do-faces ((site axis1 axis2) &body body)
  `(loop for (,axis1 ,axis2) in '((1 2) (2 3) (3 1))
	 do (do-lattice (,site ,axis1 ,axis2) ,@body)))

;;Access the phase lattice.  SITE should already be standardized
(defmacro vv-lattice-ref (site)
  `(lattice-ref *vv-psi* ,site))

;;VV-status system.  Get status of site, which should already be standardized.
(defmacro vv-status (site)
  `(lattice-ref *vv-status* ,site))

;;The lowest order 4 bits tell whether the cube intersects each successor, including its future
(define-constant vv-status-all-successors-byte (byte 4 0))

(defconstant vv-status-ours-bit 5)	   ;Are we doing it?
;;Was cube done by predecessor?  If so, no other bits matter.  DONE means that the predecessor created
;;phases, face-points, and either diamonds or pseudo-diamonds as necessary.
(defconstant vv-status-predecessor-bit 6)

;;Does this cell intersect successor (including future)?  First 4 bits
(defmacro vv-status-successor-p (site successor)
  `(logbitp ,successor (vv-status ,site)))
(defun vv-status-all-successors (site)
  (ldb vv-status-all-successors-byte (vv-status site)))

;;Are we doing this cell?
(defmacro vv-status-ours-p (site)
  `(logbitp ,vv-status-ours-bit (vv-status ,site)))

;;Did predecessor do it?
(defmacro vv-status-predecessor-p (site)
  `(logbitp ,vv-status-predecessor-bit (vv-status ,site)))


(defun describe-vv-status (site)
  (format t "~&~S = ~S: predecessor: ~A, ours: ~A, successors ~4,'0B"
	  site (site-list site)
	  (vv-status-predecessor-p site) (vv-status-ours-p site) (vv-status-all-successors site))
  (force-output))

;;Return third direction such that direction1, direction2, third direction
;;are a righthanded triplet.
(defun right-hand-rule (direction1 direction2)
  (let ((axis1 (abs direction1))
	(axis2 (abs direction2)))
    (* (signum direction1) (signum direction2)
       (ecase (mod (- axis2 axis1) 3)	;see if in order
	 (1 1)
	 (2 -1))
       (- 6 axis1 axis2))))

(defun left-hand-rule (direction1 direction2)
  (right-hand-rule direction2 direction1))

;;Workhorse for face-ref and face-set.  Given faces list and
;;coordinates, returns array, new site, and relative sign.
(defun get-face-array (faces site direction1 direction2)
  (multiple-value-bind (site axis1 axis2 sign)
      (canonicalize-face site direction1 direction2)
    (declare (type site site)
	     (type axis axis1 axis2))
    (values (aref faces
		    (+ (- (* (1- axis1) 3)
			  (the integer
			    (/ (* axis1 (1- axis1)) 2)
			    )
		       )
		       (- axis2 axis1 1)))
	      site sign)))

;;Put face in standard form.  Returns SITE AXIS1 AXIS2.  Axis variables are positive
;;and in order.  SIGN is +/-1 depending on whether the new face is now the inverse of the original
(defun canonicalize-face (site direction1 direction2)
  (declare (type site site)
	   (type direction direction1 direction2))
  (let ((sign 1))
    (declare (type (integer -1 1) sign))
    ;;If direction is negative, then adjust position and change sign
    (when (minusp direction1)
      (setq site (site-step site direction1))
      (setq direction1 (- direction1)
	    sign (- sign)))
    (when (minusp direction2)
      (setq site (site-step site direction2))
      (setq direction2 (- direction2)
	    sign (- sign)))
    ;;If directions in wrong order, reverse sign.
    (when (> direction1 direction2)
      (rotatef direction1 direction2)
      (setq sign (- sign)))
    (values site direction1 direction2 sign)))

(mirror-images
;;Find cube on the side of the face given by the right hand rule
(defun face-right-cube (site direction1 direction2)
  (multiple-value-bind (site axis1 axis2 sign)
      (canonicalize-face site direction1 direction2)
    (let ((direction3 (* sign (right-hand-rule axis1 axis2)))) ;direction into our cube by RH rule
      ;;If negative, must step to correct site
      (when (minusp direction3) (setq site (site-step site direction3))))
    site)))

;;See if face has string through it according to Vachaspati-Vilenkin condition.
(defun face-string (site direction1 direction2)
  (let* ((site2 (site-step site direction1))
	 (site3 (site-step site2 direction2))
	 (site4 (site-step site direction2)))
    (signum
     (+ (vv-winding site site2)	;Sum of windings is -3, 0, or +3
	(vv-winding site2 site3)
	(vv-winding site3 site4)
	(vv-winding site4 site)))))

;;Returns winding according to vv phases: -1, 0, or +1
(defun vv-winding (from-site to-site)
  (let* ((from (vv-lattice-ref from-site))
	 (to (vv-lattice-ref to-site)))
    (1- (mod (1+ (- to from)) 3))))

;;Call function with (site direction1 direction2) for each face of cube.
;;Sites are standardized.  Orientation is always outward. 
(defun map-cube-faces (fun site)
  (funcall fun site 2 1)
  (funcall fun site 3 2)
  (funcall fun site 1 3)
  (funcall fun (site-step site 1) 2 3)
  (funcall fun (site-step site 2) 3 1)
  (funcall fun (site-step site 3) 1 2))

(defstruct face-point
  (location nil)			;4vector at center
  (left	nil)				;Diamonds to left and right of this face.
  (right nil)				;The vector given by the right hand rule points from left to right
  )


;;Return the face from the *vv-faces* found by starting at SITE and moving first in the
;;direction given by DIRECTION1 and then in the direction given by DIRECTION2.
;;Second value NIL if face has canonical orientation, T if reverse.
;;SITE should already be standardized
(defun vv-face-ref (site direction1 direction2)
  (multiple-value-bind (array site sign)
       (get-face-array *vv-faces* site direction1 direction2)
    (values (lattice-ref array site) (minusp sign))))

;;Set value of given face.  This is insensitive to orientation.
;;SITE should already be standardized
(defsetf vv-face-ref vv-face-set)
(defun vv-face-set (site direction1 direction2 new)
  (multiple-value-bind (array site)
      (get-face-array *vv-faces* site direction1 direction2)
    (setf (lattice-ref array site) new)))

;;Get at slots of face-point structure
(defun face-ref-location (site direction1 direction2)
  (let ((structure (vv-face-ref site direction1 direction2)))
    (and structure (face-point-location structure))))

(defsetf face-ref-location face-set-location)
(defun face-set-location (site direction1 direction2 new)
  (let ((structure (vv-face-ref site direction1 direction2)))
    (cond (structure			;structure exists already
	   (setf (face-point-location structure) new)) ;fill slot
	  (t
	   (setf (vv-face-ref site direction1 direction2)
		 (make-face-point :location new))
	   new))))

;;LEFT means this side of the face from which it is directed by the right hand rule.
;;IT doesn't have to do with the sense of the string that might go through this face.
(mirror-images
(defun face-ref-right (site direction1 direction2)
  (multiple-value-bind (structure reversed) (vv-face-ref site direction1 direction2)
    (and structure
	 (if reversed (face-point-left structure) (face-point-right structure)))))
(defsetf face-ref-right face-set-right)
(defun face-set-right (site direction1 direction2 new)
  (multiple-value-bind (structure reversed)  (vv-face-ref site direction1 direction2)
    (cond (structure			;structure exists already
	   (if reversed
	       (setf (face-point-left structure) new) ;fill opposite slot
	     (setf (face-point-right structure) new)))
	  (t
	   (setf (vv-face-ref site direction1 direction2)
		 (if reversed (make-face-point :left new) (make-face-point :right new)))
	   new))))
)					;mirror-images


(defvar *global-vv-size*)		;Number of cells in each direction of the global array
(defvar *vv-cell-size*)			;Size of each cell
(defvar *vv-offsets*)			;Site of our (0,0,0) cell in the global array.

(declaim (double-float *vv-cell-size*)
	 (fixnum *global-vv-size*)
	 (site *vv-offsets*))

;;Step in given DIRECTION, standardize resulting site.
(defun site-step (site direction &optional (amount 1) (error-p t))
  (declare (type site site)
	   (type direction direction)
	   (type fixnum amount))
  (let ((axis (abs direction)))
    (when (minusp direction)
      (setq amount (- amount)))
    (let ((new (+ (site-ref axis site) amount)))
      (when (minusp new)
	(incf new *global-vv-size*))	;Find equivalent site with positive coordinates
      (setf (site-ref axis site) new))
    (standardize-vv site t error-p)))

;;The magnitude of the largest perturbation to the initial points.
(define-simulate-variable *max-position-perturbation* 0.1 double-float)

;;The maximum (transverse) speed of an initial string segment.
(define-simulate-variable *max-initial-speed* 0.5 double-float)

;;Return a position perturbation vector, chosen randomly from square of half-side *max-position-perturbation*
;;lying on the face
(defun position-perturbation (axis1 axis2)
  (let ((result (make-zero-3vector)))
    (incf (aref result (1- axis1)) (* (- (random 2.0) 1.0) *max-position-perturbation*))
    (incf (aref result (1- axis2)) (* (- (random 2.0) 1.0) *max-position-perturbation*))
    result))

;;Return a random velocity for the string whose tangent vector is given.  It will lie in a random direction
;;perpendicular to the tangent, and its squared length is uniformly distributed between 0 and *max-initial-speed*^2
;;This is similar to the velocity distribution that goes with uncorrelated a and b.
(defun random-velocity (tangent)
  (random-transverse-vector tangent *max-initial-speed*))

;;Return 2 unit vectors perpendicular to each other and the given vector
;;The first value cross the second value gives the original vector (if it is normalized)
(defun two-perpendicular-vectors (vector)
  (let* ((x (3vector-x vector))
	 (y (3vector-y vector))
	 (z (3vector-z vector))
	 (r2 (+ (expt x 2) (expt y 2))))
    (if (plusp r2)			;If original vector in z direction, use x and y.
	(values (3vector-normalize (make-3vector y (- x) 0.0)) ;otherwise get 2 perpendicular and normalize
		(3vector-normalize (make-3vector (* z x) (* z y) (- r2))))
	(values (make-3vector 1.0 0.0 0.0) (make-3vector 0.0 1.0 0.0))
	)))

(defun random-transverse-vector (tangent max-magnitude)
  (if (zerop max-magnitude) zero-3vector
    (multiple-value-bind (v1 v2)
	(two-perpendicular-vectors tangent)
      (let* ((angle (random (* 2 pi)))
	     (speed2 (random (expt max-magnitude 2)))
	     (speed (sqrt speed2)))
	(3vector+ (3vector-scale v1 (* speed (sin angle))) (3vector-scale v2 (* speed (cos angle))))))))
	 
;;Maximum size of the VV cell.  It will be reduced somewhat to be compatible with the periodicity constraints.
(defvar maximum-vv-cell-size 1.0)

;;Determine the size of the largest diamond that can be result from the initial conditions,
;;which is twice the maximum length of a p or q vector.  That in turn is half the maximum length of
;;a perturbed vector across a cell times the Lorentz boost.
(defun setup-diamond-span ()
  (setq diamond-span 
	;;Position perturbations could be in opposite directions, but they are transverse to the unperturbed vector
	(/ (sqrt (+ (expt maximum-vv-cell-size 2) (expt (* 2 *max-position-perturbation*) 2)))
	   (sqrt (- 1.0 (expt *max-initial-speed* 2))))))
		 
;;The coordinate of the global origin in the VV lattice.  Because the lattice has an even number of
;;elements numbered from 0 to size-1, this is (size-1)/2, not an integer
(defun global-vv-origin ()
  (/ (1- *global-vv-size*) 2.0))

;;Standardize a set of VV coordinates, by choosing the image that is closest to the center of the lattice.
;;We can add or subtract N = *global-vv-size*/2 from one coordinate if we add or subtract it from another.
;;Global lattice coordinates range from 0 to 2N-1, so the center of the lattice is at N-1/2.  
;;To avoid ambiguities, we treat this as N-1/2+epsilon, so that N is closer than N-1.  N-epsilon is equivalent.
;;In the case of the local lattice (which sometimes still has multiple images of a point), let S be
;;*vv-size*.  If S is even we use (S-1)/2+epilson, or S/2-epsilon as the center.
;;If S is odd, we use (S-1)/2-epsilon.  In either case, we have some integer -epsilon.
;;We first add or subtract a multiple of N to each coordinate to make it as close as possible to the center.
;;If N is odd, there may be 2 possible such values, and we take the smaller.
;;That if we have added an odd number of N's, we find the coordinate furthest from the center and swap it
;;to the other side.  If there is a tie, we use the first one.
;;If we are working locally, there might be no image of the given site in the local lattice.  In that case
;;we signal an error if ERROR-P is set, and otherwise return a site not in the local lattice and second value T
(defun standardize-vv (site local-p &optional (error-p t))
  (declare (optimize speed))
  (let* ((n (/ *global-vv-size* 2))	;Only the global lattice has actual periodicity
	 (goal (if local-p (truncate *vv-size* 2) n)) ;Goal is this -epsilon
	 (min (- goal (truncate n 2))) ;The smallest value that is not improved by adding N.
	 (site site)			;Avoid warnings about changing arguments defined in function type proclamation
	 (total 0))
    (declare (type site site)
	     (fixnum n goal min total)
	     (inline floor))		;It is declared maybe-inline in the sbcl source
    (loop for axis from 1 to 3
	  do (multiple-value-bind (quotient remainder) ;Get number of N's to subtract to get into range goal...goal+n-1
		 (floor (- (site-ref axis site) min) n)
	       (declare (fixnum quotient remainder)
			(optimize (safety 0))) ;Don't check that they really are fixnums
	       (when (minusp (+ remainder min)) ;Might have min<0
		 (incf remainder n)
		 (decf quotient 1))
	       (incf total quotient)
	       (setf (site-ref axis site) (+ remainder min)) ;map into desired range
	       ))
    (when (oddp total)			;If we added even number of N's, we're done.  Otherwise must do one more
      (let* ((best-axis 1)
	     (best-value (site-ref best-axis site)))
	(flet ((check-axis (axis)	;Go over other axes looking for best slot to change
		 (let ((value (site-ref axis site)))
		   ;;If this value is further from its goal than the best value is from its, then we should modify
		   ;;this one.  If the distances are the same we apply epsilon: if this one is below, while
		   ;;the previous best is above, then the current one is really further.
		   (when (or (> (abs (- value goal)) (abs (- best-value goal))) ;New choice further
			     (and (= (abs (- value goal)) (abs (- best-value goal))) ;Equally good.  Apply epsilon
				  (>= value goal) (< best-value goal)) ;New one on high side, so really further
			     )		;If new and old values are the same, we keep the old one
		     (setq best-axis axis best-value value)))))
	  (declare (inline check-axis))
	  (check-axis 2)		;Check remaining 2 axes for improvements
	  (check-axis 3)
	  (let* ((old (site-ref best-axis site)) ;Adjust once more, so total number of changes is even
		 (new (if (< old goal)	;Smaller than goal now
			  (+ old n)	;Add
			(- old n))))	;Bigger: subtract
	    (declare (fixnum new old)
		     (optimize (safety 0))) ;Don't check that they really are fixnums
	    ;;Now we have the desired value, but if it is negative, we can't actually install it
	    (when (minusp new)
	      (unless local-p
		(error "Unexpected trouble in standardizing a global site"))
	      (incf new *global-vv-size*)) ;Get a valid site.  It is not in the local lattice 
	    (setf (site-ref best-axis site) new) ;Install new value
	    ))))
    ;;Now check that it is valid
    (if (and local-p
	     (loop for axis from 1 to 3
		   thereis (>= (site-ref axis site) *vv-size*))) ;Too large?
      (if error-p
	  (error "The given site has no image in the local lattice")
	(values site t)) ;Return outside site and error indicator
      site)))

;;Convert VV location indices to global
(defun globalize-vv (site)
  (loop for axis from 1 to 3
	do (setf (site-ref axis site) (mod (+ (site-ref axis site) (site-ref axis *vv-offsets*)) *global-vv-size*)))
  (standardize-vv site nil))

;;Convert VV location indices to local
(defun localize-vv (site)
  (loop for axis from 1 to 3
	do (setf (site-ref axis site) (mod (- (site-ref axis site) (site-ref axis *vv-offsets*)) *global-vv-size*)))
  (standardize-vv site t))

;;Convert global VV location into unstandardized global position 4vector
;;AXES are directions in to offset the result by half a cell, allowing us to
;;get the positions of edge, face, or cube centers
(defun unstandardized-vv-position (site &optional axes)
  (3to4vector
   (list-3vector
    (loop for axis from 1 to 3
	  collect (* (+ (if (member axis axes) 0.5 0.0)
			(- (site-ref axis site) (global-vv-origin)))
		     *vv-cell-size*)))
   (global-time *initial-time*)		;Put it at the initial time because that's when we are doing initial conditions
   ))

;;Convert global VV location into standard global position
(defun global-vv-position (location &rest axes)
  (standardize-position (unstandardized-vv-position location axes)))

;;Convert local VV location into local position
(defun vv-position (location &rest axes)
  (localize-position (unstandardized-vv-position (globalize-vv location) axes)))

;;Convert global position into into global VV location.  We mod the result by the global size, because
;;negative numbers can't be stored in sites, but don't standardize otherwise.
(defun unstandardized-position-vv (position &optional (round :down))
  (loop with site = (make-site)
	for axis from 1 to 3
	as raw = (+ (global-vv-origin)
		    (/ (+ (3vector-component position (1- axis)) fudge-global-coordinates) *vv-cell-size*)) ;floating point result
	do (setf (site-ref axis site)
		 (mod (ecase round
			(:down (floor raw))
			(:up (ceiling raw))
			(:nearest (round raw)))
		      *global-vv-size*))
	finally (return site)))

;;Convert global position into into global VV location
(defun global-position-vv (position &optional (round :down))
  (standardize-vv (unstandardized-position-vv position round) nil))

;;Convert local local position into VV location.
(defun position-vv (position)
  (localize-vv (unstandardized-position-vv (globalize-position position))))


;;;Iterators

;;Iterate over the sites in the lattice, omitting those that are non-standard representations of the same points
;;If ANY-STATUS is given (at compile time), do them all.
;;Otherwise, do those that adjoin cubes handled by us.
;;AXES is a list of forms yielding directions from the site that we have already gone half a step.
;;We move only perpendicular to these when looking for looking for adjoining cubes.  So if 2 are given, we check
;;the two cubes adjoining this face, and so on.  The sites of the lattice have the full range from 0
;;to *vv-size*.  Everything else is constrained to be within that volume, so
;;we don't consider, for example, faces that stick outward from the last lattice point.
;;Do not change the value of the SITE variable in the body.
(defmacro do-vv-lattice ((site &key axes any-status)  &body body)
  (unless (listp axes) (error "AXES should be a list, not ~S" axes))
  (unless (member any-status '(t nil))
    (warn "ANY-STATUS needs to be known at compile time.  Probably you're making a mistake by supplying ~S" any-status))
  (let* ((axes-var (and axes (gensym)))	;Generate variable to store cached axes, or use NIL if none.
	 (code
	  `(do-lattice (,site ,@axes)	;Go over all cells in our lattice but not external bits
	     (when (= ,site (standardize-vv ,site t)) ;Not if non-standard
	       ,(if any-status
		    `(progn ,@body)	;No extra check
		  `(when (vv-site-interesting-p site :axes ,axes-var)
		     ,@body))))))
    (if (and axes				;If axes supplied, we should make list only once
	     (not any-status))			;But we won't use it if this was set
	`(let ((,axes-var (list ,@axes)))
	   ,code)
      code)))

;;Iterate over faces that adjoin our cubes.
(defmacro do-vv-faces ((site axis1 axis2 &key any-status) &body body)
  `(loop for (,axis1 ,axis2) in '((1 2) (2 3) (3 1))
	 do (do-vv-lattice (,site :axes (,axis1 ,axis2) :any-status ,any-status)
	      ,@body)))

;;Iterate over our cubes.
(defmacro do-vv-cells ((site &key any-status) &body body)
  `(do-vv-lattice (,site :axes (1 2 3) :any-status ,any-status)
     ,@body))

;;Call function with cubes that adjoining SITE, offset by half a cell in the given directions,
(defun map-adjoining-cells (function site axes)
  (dotimes (x (if (member 1 axes) ;Consider stepping backward in x, unless X among axes
		  1 2))
    (dotimes (y (if (member 2 axes)
		    1 2))
      (dotimes (z (if (member 3 axes)
		      1 2))
	(let ((new site)
	      (error nil))
	  ;;Optionally take steps in 3 directions.  ERROR is set if final result could not be standardized.
	  (when (plusp x) (multiple-value-setq (new error) (site-step new 1 -1 nil)))
	  (when (plusp y) (multiple-value-setq (new error) (site-step new 2 -1 nil)))
	  (when (plusp z) (multiple-value-setq (new error) (site-step new 3 -1 nil)))
	  (unless error			;Don't consider sites that are not in the lattice at all
	    (funcall function new)))))))

(declaim (inline last-common-predecessor)) ;This allows calls to member with constant axes to be optimized out

;;Find the last common predecessor of an object of those cubes that adjoin the given site, offset half a cell
;;along the given directions.
;;We work by looking at the corners of the brick formed by taking steps 0 and 1 along the directions given
;;and -1 and 1 along the directions not given.
;;We return the LOCAL IJKL of the last common predecessor, which are all 0 if that is us.
(defun last-common-predecessor (site axes)
  (let ((min-ijkl (make-4vector most-positive-double-float most-positive-double-float
				most-positive-double-float most-positive-double-float))
	;;Convert VV to position once and for all.  All other corners will be offset from this one, so that 
	;;we work in a single image
	(position (vv-position site)))
    (declare (type 4vector position min-ijkl)
	     (optimize speed))
    (macrolet ((offset (variable axis) ;Variable is 0 or 1.  Change 0 to -1 if axis not listed, then mult by cell size
		`(without-compiler-notes ;No duplicate code warnings in case of constant axes
		  (* (if (member ,axis axes) ,variable (- (* ,variable 2) 1)) *vv-cell-size*))))
      (dotimes (x 2)			;Loop over corners of cube
	(dotimes (y 2)
	  (dotimes (z 2)
	    (let* ((x (4vector+ position (make-4vector (offset x 1) (offset y 2) (offset z 3) 0.0))) ;Location of this corner
		   (ijkl (xtoi x)))	;in local IJKL coordinates
	      (dotimes (index 4)	;Keep minimum coordinate values
		(setf (aref min-ijkl index)
		      (min (aref min-ijkl index) ;Keep minimum value in each coordinate
			   ;;In case of uncertainty, we consider the point to be in the earlier of the
			   ;;2 possibilities.  This may still not fix problems where, for example,
			   ;;a point is not ours even though an adjoining face is because of floating point error
			   ;;about the corners
			   (fast-ffloor (- (aref ijkl index) fudge-ijkl))))
		)
	      (deallocate 4vectors x ijkl))))))
    ;;MIN-IJKL is now the last common predecessor in the covering space.  Convert to standard representative
    ;;of this volume and return
    (prog1 (standardize-ijkl min-ijkl)
      (deallocate 4vectors min-ijkl)
      )))

;;Tell if we need to deal with the object given by SITE offset in the given directions.
;;If a predecessor did any adjoining cell, we don't do this object.
;;Otherwise, we see if we are the LCP
(defun vv-site-interesting-p (site &key axes)
  (map-adjoining-cells
   #'(lambda (site)
       (when (vv-status-predecessor-p site) ;If any cell already done by predecessor, don't do
	 (return-from vv-site-interesting-p nil)))
     site axes)
  ;;If predecessor did one, we already returned NIL.  Otherwise see if we are LCP.
  (4vector= (last-common-predecessor site axes) zero-4vector fudge-ijkl))

(defun vv-successor-needs-site-p (site successor &key axes)
  (map-adjoining-cells
   #'(lambda (site)
       (when (vv-status-successor-p site successor)		      ;Successor wants any adjacent cell?
	 (return-from vv-successor-needs-site-p t))		      ;return T
       )
     site axes)
  nil)					;No interesting cells

;;Setup status array.  Initialized to 0 on creation.  Then predecessor's data read.
;;If predecessor didn't do a cell, we consider doing it.
(defun setup-vv-status ()
  (declare (optimize speed))
  (do-vv-lattice (site :any-status t)		       ;All sites
    (loop for successor below 4
	  do (when (vv-cube-intersects-volume-p site :successor successor :future t) ;Does it intersect at all?
	       (setf (vv-status-successor-p site successor) t)))			     ;Say so
    (unless (vv-status-predecessor-p site)		   ;If already done by predecessor, not ours to do
      (setf (vv-status-ours-p site)			   ;See if If LCP is us
	    (4vector= (last-common-predecessor site '(1 2 3)) zero-4vector fudge-ijkl)))))

;;Tell whether any point in a VV cube is in our simulation volume.
;;If SUCCESSOR, consider instead his simulation cube and all future cubes.
(defun vv-cube-intersects-volume-p (site &key successor future)
  (let ((position (vv-position site)))
    (when successor			;Wants successor's volume instead?
      (setq position			;Switch to successor's coordinates
	    (4vector- position (destination-offset successor))))
    (vv-cube-position-intersects-volume-p position future)))

;;See if VV cube starting from the given 4vector position intersects the cube IJKL=0...1.
;;If FUTURE-P is set, check instead for IJKL>=0
;;It appears to be sufficient to check whether either one of the cube vertices is inside the volume
;;or whether an edge of the cube intersects a face of the volume.
(defun vv-cube-position-intersects-volume-p (position future-p)
  (let ((corners (make-array '(2 2 2) :initial-element nil)))
    (prog1
	(block done
	  (dotimes (x 2)		;Set up ijkl locations of corners
	    (dotimes (y 2)
	      (dotimes (z 2)
		(let* ((xyzt (4vector+ position (make-4vector ;Location of this corner
						 (* x *vv-cell-size*) (* y *vv-cell-size*) (* z *vv-cell-size*) 0.0)))
		       (ijkl (xtoi xyzt)))
		  (deallocate 4vectors xyzt)
		  (setf (aref corners x y z) ijkl) ;Install it
		  (when (loop for index below 4 ;Check if corner is inside region
			      always (and (< (- fudge-ijkl) (aref ijkl index))
					  (or future-p (< (aref ijkl index) (+ 1.0 fudge-ijkl)))))
		    (return-from done t)))))) ;Return true if any corner inside
	  ;;Now check for edges intersecting volume boundaries
	  (dotimes (direction 3)	;Direction to go in
	    (dotimes (x (if (= direction 0) 1 2)) ;Go over corners, but if moving in X, only start from x = 0
	      (dotimes (y (if (= direction 1) 1 2))
		(dotimes (z (if (= direction 2) 1 2))
		  (let ((first (aref corners x y z))
			(second (aref corners (if (= direction 0) 1 x) ;Corner of cube given by step in direction direction
				      (if (= direction 1) 1 y)
				      (if (= direction 2) 1 z))))
		    (when (line-intersects-simulation-face-p first second future-p)
		      (return-from done t)))))))
	  nil)					;Didn't find any intersection, give NIL
      (dotimes (x 2)			;Deallocate corner vectors
	(dotimes (y 2)
	  (dotimes (z 2)
	    (let ((corner (aref corners x y z)))
	      (when corner (deallocate 4vectors corner))))))
      )))
    
;;See if the line between two given ijkl 4vectors crosses a face of our simulation volume
;;If FUTURE is set, check for intersections not with our volume, but with the edge of the volume that includes
;;our cube and all cubes in the future of our cube.
(defun line-intersects-simulation-face-p (from to &optional future)
  (declare (optimize speed)
	   (type 4vector from to))
  (dotimes (direction 4)		;Check faces perpendicular to this direction
    (let ((denominator (- (aref to direction) (aref from direction))))
      (unless (< (abs denominator) fudge-ijkl) ;Parallel to edge: say no intersection
	(loop for value below (if future 1 2)  ;If future, only consider past edges
	      ;;Find parameter (0 = from; 1 = to) at which line intersects 3-plane of face
	      for parameter double-float = (/ (- (float value 0.0) (aref from direction)) denominator)
	      when (and (< 0.0 parameter 1.0) ;Valid parameter range
			;;Find position of intersection
			(let ((position (4vector+ (4vector-scale from (- 1.0 parameter)) (4vector-scale to parameter))))
			  (loop for index below 4
				always (or (= index direction) ;Don't check direction we just solved
					   (and (> (+ (aref position index) fudge-ijkl) 0.0) ;must be inside face
						(or future				     ;Future: anything > 0 OK.
						    (< (- (aref position index) fudge-ijkl) 1.0)))))))
	      do (return-from line-intersects-simulation-face-p t)))))
  nil					;No intersection
  )


;;See if the line between two given 4vectors contains any point in our volume
(defun line-overlaps-simulation-p (from to)
  (or (point-mine from)			;If one in inside, yes.  This handles all-inside case
      (point-mine to)			;Quickly check other end as well
      (line-intersects-simulation-face-p from to) ;Might cross a corner of our volume.
      ))

;;Create a junction for diamonds that meet at a face.  The junction can have an orientation.
(defun create-vv-junction (site axis1 axis2)
  (make-vv-junction :site site :axis1 axis1 :axis2 axis2))


;;Here are the functions called by SETUP-VV

(defun create-vv-lattice ()
  (let ((total (expt *vv-size* 3)))
    (setq *vv-psi* (make-array total :initial-element nil) ;phase, or NIL if not initialized
	  *vv-status* (make-array total :element-type 'vv-status-type :initial-element 0)
	  *vv-faces* (make-faces total :initial-element nil))))

;;Make sure we are far from places where the volume changes shape
(defun check-initial-vv-time ()
  (let* ((time (- *initial-time* (job-start-t))) ;Time since start of volume
	 (side-time (/ (job-duration) 4)) ;job-duration/4 is the time between cube vertices
	 (offset (/ (mod time side-time) side-time)) ;Fractional distance from last vertex plane
	 (edge (* *size* (min offset (- side-time offset)))) ;edge length of smallest tetrahedron
	 (radius (/ edge 2 (sqrt 6.0)))) ;radius of inscribed sphere
    (unless (> radius (+ (* (sqrt 3.0) *vv-cell-size*) fudge-coordinates)) ;must be larger than main diagonal of cell
      (error "Initial time too close to cube corner time or size too small"))))

;;Set up our part of the lattice.  We create a rectangular lattice large enough to contain anything we could care
;;about.  Probably this could be improved.
(defun setup-vv-lattice ()
  ;;Choose an even number of cells in the cardinal directions such that cells are no larger than desired
  (let ((coordinate-range (* *total-size* (sqrt 2.0)))) ;Total range in each coordinate
    (setq *global-vv-size* (* 2 (ceiling coordinate-range (* maximum-vv-cell-size 2))) ;Round up to even
	  *vv-cell-size* (/ coordinate-range *global-vv-size*))) ;Size of each cell
  (let ((minima (copy-3vector *global-location*)) ;Construct minima and maxima in each coordinate
	(maxima (copy-3vector *global-location*)))
    (dotimes (index 3)
      (incf (aref minima index)
	    (- (job-slice-coordinate-minimum *initial-time*)
	       diamond-span ;Any point within this range is potentiallly of interest
	       *vv-cell-size*	;This provides room for an uninteresting boundary cell.
	       *vv-cell-size*	;We also need an extra cell side to allow for common successors.  See the notes.
	       ))
      (incf (aref maxima index) (+ (job-slice-coordinate-maximum *initial-time*)
				   *vv-cell-size* *vv-cell-size* diamond-span)))
;;    (format t "~& Maxima = ~S, Minima = ~S" maxima minima)
    (setq *vv-offsets* (unstandardized-position-vv minima :down)) ;Convert coordinates separately to VV position
    (let ((max-vv (unstandardized-position-vv maxima :up))) ;Round so as to include slot outside our region
;;      (format t "~& Maxima = ~S, Minima = ~S" (site-list max-vv) (site-list *vv-offsets*))
      (setq *vv-size*			;Big enough for any axis (anisotropy possible only by rounding)
	    (loop for axis from 1 to 3
		  for size = (mod (- (site-ref axis max-vv) (site-ref axis *vv-offsets*)) *global-vv-size*)
		  when (zerop size) do (setq size *global-vv-size*) ;This much overlap OK
		  when (< size 3)
		  do (error "VV size came out to be ~D, meaning that there isn't enough room between the job volume and the overall lattice" size)
		  maximize size)))
    (create-vv-lattice))
  nil)

;;Install random phases 0, 1 or 2, in the interior cells of the lattice and boundaries that were not filled already
(defun create-vv-phases ()
  (do-vv-lattice (site)				 ;Examine all corners of cells that we're doing
    (unless (vv-lattice-ref site)		 ;Already set by predecessor?
      (setf (vv-lattice-ref site) (random 3)))	 ;No: install random phase
    ))

;;Create a point in the center of every face that has a string going through it, unless done by predecessor
(defun install-face-points () 
  (do-vv-faces (site axis1 axis2)		   ;Examine all faces of cubes that we need to do
    (unless (zerop (face-string site axis1 axis2)) ;If there's any string through this face
      (unless (vv-face-ref site axis1 axis2) ;Done already
	(install-face-point site axis1 axis2)))))

(defun install-face-point (site axis1 axis2)
  (let ((position (vv-position site axis1 axis2)))
    (setq position (3vector+ position (position-perturbation axis1 axis2))) ;Perturb position
    (setq position (3to4vector position *initial-time*)) ;Convert to 4vector
    (setf (face-ref-location site axis1 axis2) position)))

(defun describe-cube-faces (site)
  (map-cube-faces
   #'(lambda (site direction1 direction2)
       (format t "~&~S+~D+~D: ~D" site direction1 direction2 (face-string site direction1 direction2))
       (force-output))
   site))

;;Make diamonds for the string going through each cube of the lattice and put them in the faces structures
(defun make-vv-diamonds ()
  (do-vv-cells (site)			;Examine the cells that we're doing
    (let ((entries nil)
	  (exits nil))
      (map-cube-faces #'(lambda (site direction1 direction2)
			  (ecase (face-string site direction1 direction2)
			    (-1 (push (list site direction1 direction2) entries)) ;Entering cube
			    (1 (push (list site direction1 direction2) exits)) ;Exiting
			    (0)		;No string
			    ))
		      site)
      (ecase (length entries)
	((0 1))				;OK, nothing special
	(2 (if (plusp (random 2));Must randomize relationship between entries and exits
	       (setq entries (nreverse entries)))))
      (loop for entry in entries
	    for exit in exits
	    do (make-vv-diamond entry exit)
	    ))))

;;Data about a diamond that we are not making but merely passing on
(defstruct pseudo-diamond
  left					;left corner
  right					;right corner
  start)				;starting position

;;Handle 1 diamond for make-vv-diamonds.  We create a diamond if it has not already been created.  If the start point
;;does not lie in our volume, we make a pseudo-diamond instead.
(defun make-vv-diamond (entry exit)
  ;;Both face-refs are LEFT, because faces are directed outward from the cube.
  (destructuring-bind (entry entry-1 entry-2) entry
    (destructuring-bind (exit exit-1 exit-2) exit
      (if (face-ref-left entry entry-1 entry-2) ;Already have string?
	  (unless (face-ref-left exit exit-1 exit-2) ;Must have this end also
	    (error "Inconsistent data structures in MAKE-VV-DIAMOND"))
	(multiple-value-bind (left-diamond right-diamond)
	    (create-vv-segment entry entry-1 entry-2 exit exit-1 exit-2) ;Create string going through diamond
	  (setf (face-ref-left entry entry-1 entry-2) left-diamond) ;Install ends in the face array
	  (setf (face-ref-left exit exit-1 exit-2) right-diamond))))))

(define-simulate-variable *split-initial-diamonds* nil fixnum)	;Number of pieces into which to split initial diamond

;;Function to split diamonds 
;;Returns list of splitting positions not including left but including right, or NIL if splitting failed
(define-simulate-variable *split-initial-diamonds-function* nil)

;;Just introduce extra points on the same segments without perturbing anything.
(defun split-initial-diamond-trivially (left right &rest ignore)
  (declare (ignore ignore))
  (append
   ;;Introduce N-1 interior points on a straight line from left to right
   (loop for j from 1 below *split-initial-diamonds*
	 for split = (3to4vector (3vector+ (3vector-scale left (/ (- *split-initial-diamonds* j) *split-initial-diamonds*))
					   (3vector-scale right (/ j *split-initial-diamonds*)))
				 *initial-time*)
	 collect split
	 unless (point-inside-past-p split) ;Generated point not inside our past, so can't advance here?
	 do (return-from split-initial-diamond-trivially nil))
   (list right)))

;;Interpolate smooth circles.  If you use this, you probably want max-initial-velocity small.
(defun split-initial-diamond-interpolate (left right entry entry-1 entry-2 exit exit-1 exit-2)
  (declare (ignore entry exit))
  (let* ((delta-x (standardize-position (3vector- right left))) ;Vector between corners of original diamond
	 (entry-index (- 5 (abs entry-1) (abs entry-2)))	;Coordinate index of entry direction. 0=x
	 (exit-index (- 5 (abs exit-1) (abs exit-2)))
	 (other-index (- 3 entry-index exit-index))) ;Index in unused direction
    (unless (= entry-index exit-index)	;If straight segment, nothing to do
      (append
       ;;Introduce N-1 interior points on a quarter circle from left to right
       (loop for j from 1 below *split-initial-diamonds*
	     for angle = (* (/ pi 2) (/ j *split-initial-diamonds*))
	     for split = (copy-4vector left) ;Copy entry point
	     do (incf (aref split entry-index) (* (aref delta-x entry-index) (sin angle)))
	     do (incf (aref split exit-index) (* (aref delta-x exit-index) (- 1 (cos angle))))
	     ;;Linear interpolation in other index
	     do (incf (aref split other-index) (* (aref delta-x other-index) (/ j *split-initial-diamonds*)))
	     do (setf (4vector-t split) *initial-time*)
	     collect split
	     unless (point-inside-past-p split) ;Generated point not inside our past, so can't advance here?
	       do (return-from split-initial-diamond-interpolate nil))
       (list right)))))
	
(defun split-initial-diamond-Fourier-amplitude (n) ;Size of perturbation in nth frequency
  (* 0.05 (expt n -0.75)))		;Fit phenomological power law P(k) = 2kA(k)^2 ~ k^-0.5
  
;;Split by random Fourier components
(defun split-initial-diamond-Fourier (left right &rest ignore)
  (declare (ignore ignore))
  (let* ((delta-x (standardize-position (3vector- right left))) ;Vector between corners of original diamond
	 (components (loop for n from 1 below *split-initial-diamonds* ;N-1 frequencies can be excited
			   collect (random-transverse-vector delta-x (split-initial-diamond-Fourier-amplitude n)))))
    (append
     (loop with split0 = (copy-3vector left) ;Advance along string
	   with split-fraction = (/ 1.0 *split-initial-diamonds*)
	   with split
	   with step = (3vector-scale delta-x split-fraction) ;amount to advance to next point
	   for j from 1 below *split-initial-diamonds* ;Introduce N-1 interior points
	   do (setq split0 (prog1 (3vector+ split0 step)
			     (deallocate 3vectors split0)))
	   do (loop with split1 = (copy-3vector split0)	;Adjust by offsets
		    for n from 1
		    for component in components
		    do (setq split1 (prog1 (3vector+ split1 (3vector-scale component
									   (sin (/ (* pi n j) *split-initial-diamonds*))))
				      (deallocate 3vectors split1)))
		    finally (setq split (3to4vector split1 *initial-time*)) ;Convert to 4vector
		    (deallocate 3vectors split1))
	   collect split
	   unless (point-inside-past-p split) ;Generated point not inside our past, so can't advance here?
	   do (return-from split-initial-diamond-Fourier nil)	;give up
	   finally (deallocate 3vectors split0 step))
     (list right))))

;;Make a string segment going through the diamond.  Returns the first and last diamonds.
;;KLUGE: At the moment, this makes a single diamond or pseudo-diamond if the diamonds
;;would otherwise start in the future.
(defun create-vv-segment (entry entry-1 entry-2 exit exit-1 exit-2)
  (let ((left (face-ref-location entry entry-1 entry-2)) ;First and last locations
	(right (face-ref-location exit exit-1 exit-2)))
    (block split			;Attempt to split up the string in this cell
      (unless *split-initial-diamonds* (return-from split nil)) ;Not requested
      (let* ((splits (funcall *split-initial-diamonds-function* left right entry entry-1 entry-2 exit exit-1 exit-2))
	     (starts (loop initially (unless splits (return-from split nil)) ;splitting failed
			   for previous = left then split
			   for split in splits
			   for start = (random-velocity-diamond-start previous split)
			   collect start
			   unless (point-inside-past-p start) ;Give up if not ours
			   do (return-from split nil))))
	;;Everything is OK.  Make new diamonds and link together
	(loop for start in starts	;Make diamond for each starting point
	      for (left-corner right-corner) on (cons left splits)
	      with first-diamond = nil
	      for previous-diamond = nil then diamond
	      for diamond = (create-vv-diamond-1 start left-corner right-corner) ;Make diamond on string
	      if previous-diamond	;First time?
	      do (advance-diamond-1 previous-diamond diamond nil) ;No, link to previous diamond
	      else do (setq first-diamond diamond) ;Yes, remember first diamond
	      finally (return-from create-vv-segment (values first-diamond diamond)))))	;Success: exit function
    ;;Here when not splitting or splitting failed.
    (let* ((start (random-velocity-diamond-start left right))
	   (diamond			;If start point is ours or in our past (see notes), actually make the diamond.
	    (if (point-inside-past-p start)
		(create-vv-diamond-1 start left right) ;Make real diamond
	      (make-pseudo-diamond :start start :left left :right right)))) ;Not ours
      (values diamond diamond))))

;;Find the starting point for a diamond with random velocity given left and right corners
;;The 2 locations might not be compatible with each other because of standardization.
;;We compute everything based on the left corner.  When we get to the starting position, we standardize that.
(defun random-velocity-diamond-start (left right)
  (let* ((delta-x (standardize-position (4vector- right left))) ;Vector reaching across diamond
	 (v (random-velocity delta-x))	;Random perpendicular velocity vector
	 (center (3vector+ left (3vector-scale delta-x 0.5))) ;Centerpoint
	 (speed  (3vector-length v))
	 (gamma  (/ 1 (sqrt (- 1 (* speed speed))))) 
	 (energy (* gamma (3vector-length delta-x)))
	 (delta-xt  (3vector-scale v energy) )
	 (delta-xt2 (3vector-scale delta-xt 0.5))
	 (p0  (/ energy 2)))
    (prog1
	(standardize-position		;New start.
	 (3to4vector (3vector- center delta-xt2) (- *initial-time* p0)))
      (deallocate 3vectors v center delta-xt delta-xt2)
      (deallocate 4vectors delta-x))))

;;When start is ours, we can make a real diamond.
(defun create-vv-diamond-1 (start left right)
  ;;Start is already standardized in attempt to get it into our volume.
  ;;It is possible that right or left is in one of our images, since it was just taken from the standardized
  ;;vv site.  We should restandardize relative to start.
  (let* ((global-start (globalize-position start))
	 (diamond (make-diamond :start start
				:left (standardize-position left start)
				:right (standardize-position right start)
					;Each diamond is considered to have been created by an event at its start
				:tag (create-loop-tag global-start)
				:a-kink-created global-start
				:b-kink-created global-start
				)))
    (setf (diamond-end diamond) (compute-diamond-end-flat diamond)) ;No expansion in initial conditions
    (handle-new-diamond diamond)
    diamond))

;;Make diamonds to the future of face points, so everything is linked together.  When we cannot do so,
;;install a junction instead.
(defun advance-vv-diamonds ()
  ;;Examine faces of cells for us and predecessors.
  ;;It could be that predecessor did diamonds on both sides of face, but only we can advance.
  (do-vv-faces (site direction1 direction2 :any-status t)
    (when (vv-face-ref site direction1 direction2)	     ;Something here?
      (let ((string (face-string site direction1 direction2))) ;Get string direction
	(when (zerop string)				       ;None?
	  (error "Face structure set, but no string"))
	(when (minusp string)				     ;String goes backward through face?
	  (rotatef direction1 direction2))		     ;Reverse orientation
	;;Now we have a string going forward through the face, meaning that as we go eastward through the diamonds
	;;we go forward through the face
	(advance-vv-diamond (face-ref-left site direction1 direction2)
			    (face-ref-right site direction1 direction2)
			    site direction1 direction2))
      )))

;;Make one future diamond connecting diamonds in different cells
;;Face is oriented to point eastward.
(defun advance-vv-diamond (west east site direction1 direction2)
  ;;If both sides of the face have been filled in and not just by pseudo-diamonds, we may be able to advance
  (cond ((and (diamondp west) (diamondp east)
	      (point-mine (diamond-right west))) ;If this point isn't ours, someone else has to do it
	 (mirror-images (setf (right-rejoining-junction west) nil)) ;Clear junctions from diamonds
	 (assert (4vector= (diamond-right west) (diamond-left east) fudge-coordinates))
	 (setf (diamond-left east) (diamond-right west))
	 (advance-diamond-1 west east nil))			    ;Make future diamond
	(t							    ;If not advancing, make junctions as needed
	 (mirror-images
	  (when (diamondp west)			    ;Real diamond to west: make junction to right of diamond
	    (unless (right-rejoining-junction west) ;Made already by predecessor
	      (setf (right-rejoining-junction west) (create-vv-junction site direction1 direction2))))))))

;;Do everything to set up a VV lattice.
(defun setup-vv ()
  (format t "~&VV setup: ") (force-output)
  (setq *vv-p* t)
  (let ((*initializing* t))
    (setup-vv-lattice)			;Figure out which part is ours, set up arrays
    (check-initial-vv-time)		;Check that we are far enough away from neighbors
    (read-initial-strings)		;Read files from predecessors, if any, with both VV data and strings
    (format t "status...") (force-output)
    (setup-vv-status)			;Figure out which cubes in the lattice are interesting to us
    (format t "phases...") (force-output)
    (create-vv-phases)			;Install phases on vertices
    (format t "faces...") (force-output)
    (install-face-points)		;Put points at center of faces
    (with-modification-group
      (format t "cell diamonds...") (force-output)
      (make-vv-diamonds)		 ;Install diamonds connecting points
      (format t "face diamonds...") (force-output)
      (advance-vv-diamonds)		 ;Put in our own next series of diamonds or junctions at ends of string
      ))
  (format t "~&VV setup done.~%") (force-output))


;;Receive and install a vertex phase from a predecessor.  Put it where he said.
(defun receive-vv-phase (phase site)
  (let* ((old (vv-lattice-ref site)))
    (when (and old (not (= old phase)))
      (error "Received phase ~D for ~{~D, ~D, ~D~}, but we already have ~D there"
	     phase (site-list site) old))
    (setf (vv-lattice-ref site) phase) ;Install phase
    ))

;;Mark that predecessor took care of a cell
(defun receive-vv-cell-done (site)
  (setf (vv-status-predecessor-p site) t))

(defun receive-vv-face-point (point site axis1 axis2)
  ;;The point might be outside our volume, so it is possible that we have the wrong image of this point
  ;;Standardizing here avoid conflicts where we receive the same point through different paths and so
  ;;get different coordinates for it.
  (setq point (standardize-position point))
  (let* ((old (vv-face-ref site axis1 axis2))
	 (old-location (and old (face-point-location old))))
    (if old-location				;Have it already?
	(unless (< (4vector-Euclidian-distance old-location point) fudge-coordinates) ;Should be same
	  (error "Received a face-point ~S for ~S, ~S, ~S, but we already have ~S"
		 point (site-list site) axis1 axis2 old))
      (let ((offset (4vector- point (vv-position site axis1 axis2)))) ;New point.  Check that it is reasonable.
	(setq offset (standardize-position offset))
	(when (> (4vector-Euclidian-length offset) (* *vv-cell-size* (sqrt 2)))
	  ;;Too far to be on face
	  (error "Point ~S is too far from the center of face ~D+~D+~D at ~S.  (Offset ~S)"
		 point (site-list site) axis1 axis2 (vv-position site axis1 axis2) offset))
	(setf (face-ref-location site axis1 axis2) point)) ;Install new point
      )))

;;Receive instructions from predecessor about where to create a diamond
;;Both faces point outward from the cell with the diamond.
(defun receive-pseudo-diamond (start face1 face2)
  (setq start (standardize-position start)) ;Standardize position.  Perhaps this should be done in communication
  (destructuring-bind (site1 axis11 axis12) face1
    (destructuring-bind (site2 axis21 axis22) face2
      (let ((old1 (face-ref-left site1 axis11 axis12))
	    (old2 (face-ref-left site2 axis21 axis22)))
	(cond ((or old1 old2)		;Perhaps this is old news
	       (unless (eq old1 old2)	;Previous version should be the same in both slots
		 (error "Inconsistent data structures in RECEIVE-PSEUDO-DIAMOND"))
	       (unless (< (4vector-Euclidian-distance
			   (if (pseudo-diamond-p old1) (pseudo-diamond-start old1) ;Previous start location
			     (diamond-start old1))				   ;Perhaps we upgraded to a diamond
			   start)
			  fudge-coordinates)
		 (error "Received a pseudo-diamond start ~S for ~S-~S, but we already have ~S"
			start face1 face2 old1)))
	      (t			;New news
	       (let ((left (face-ref-location site1 axis11 axis12))
		     (right (face-ref-location site2 axis21 axis22)))
		 (mirror-images
		  (unless left
		    (error "Predecessor supplied pseudo-diamond whose ~(~A~) face has no point" :west)))
		 (let ((diamond		       ;Create diamond if we own start point.  Otherwise pass on pseudo
			(if (point-inside-past-p start) ;Can we create it?  See notes.
			    (create-vv-diamond-1 start left right) ;Do it
			  (make-pseudo-diamond :start start :left left :right right))))
		   ;;Whatever we did, put it in structure to find later
		   ;;Code here is somewhat duplicative of make/create-vv-diamond.
		   (setf (face-ref-left site1 axis11 axis12) diamond
			 (face-ref-left site2 axis21 axis22) diamond)
		   ))))))))

;;When VV junction appears at start of string from predecessor, put diamond in lattice
;;so we can (perhaps) connect it later.  The face in the junction points forward along the string,
;;so the right of the faces points toward the string.
(defun handle-vv-start-junction (junction)
  (setf (face-ref-right (vv-junction-site junction) (vv-junction-axis1 junction) (vv-junction-axis2 junction))
	(junction-right-diamond junction)))

;;End of string.  Now the face points away, so the left of the face as the string.
(defun handle-vv-end-junction (junction)
  (setf (face-ref-left (vv-junction-site junction) (vv-junction-axis1 junction) (vv-junction-axis2 junction))
	(junction-left-diamond junction)))


;;Output of initial data to successors

(defun write-vv-phases ()
  (do-vv-lattice (site :any-status t) ;Successors might need data from predecessors
    (let ((phase (vv-lattice-ref site))) ;Have a phase?
      (when phase
	(dotimes (successor 4)		       ;Check all successors
	  ;;Phase is interesting to him or his successors?
	  ;;Really this is too generous.  If it is only of interest to our grandsuccessor, we send it to several
	  ;;immediate successors.
	  (when (vv-successor-needs-site-p site successor)
	    (transmit-vv-phase successor site phase)))))))

(defun write-vv-face-points ()
  (do-vv-faces (site axis1 axis2 :any-status t)
    (let* ((structure (vv-face-ref site axis1 axis2))
	   (position (and structure (face-point-location structure))))
      (when position			;Is there anything there to send?
	;;If we have actual diamonds on both sides, the successor will not need the face point for anything.
	(unless (and (diamondp (face-point-right structure))
		     (diamondp (face-point-left structure))) ;If both, we are done with this
	  (dotimes (successor 4)					;Check all successors
	    (when (vv-successor-needs-site-p site successor :axes (list axis1 axis2))
	      (transmit-vv-face-point successor site axis1 axis2 position))))))))

(defstruct pseudo-diamond-data
  diamond
  left-face				;(site axis1 axis2)
  right-face)

(defun write-pseudo-diamonds-and-status ()
  ;;Data storage: index: which of the 2 possible diamonds in this cube are we dealing with
  (let ((data (make-array 2)))
    (dotimes (i 2) (setf (aref data i) (make-pseudo-diamond-data))) ;make structures to reuse
    (do-vv-cells (site :any-status t)			    ;Go over all cubes where anything could be
      (when (or (vv-status-predecessor-p site)		    ;Predecessor did this cube?
		(vv-site-interesting-p site :axes '(1 2 3)))	    ;Or we did it
	(dotimes (i 2) (setf (pseudo-diamond-data-diamond (aref data i)) nil)) ;Clear structures
	(map-cube-faces
	 #'(lambda (site direction1 direction2)			       ;Outward-directed face
	     (let ((diamond (face-ref-left site direction1 direction2))) ;Diamond inside cube
	       (when (pseudo-diamond-p diamond)				 ;Ours to handle?
		 (loop for number below 2
		       until (eq diamond (pseudo-diamond-data-diamond (aref data number))) ;Seen already
		       finally
		       (when (= number 2) ;Didn't find it
			 (setq number
			       (if (pseudo-diamond-data-diamond (aref data 0))		;First slot used?
				   1							;Used second
				 0))							;Use first
			 (setf (pseudo-diamond-data-diamond (aref data number)) diamond)) ;Install diamond in slot
		       (let ((face (list site direction1 direction2)))
			 ;;Find which side of diamond goes to this face
			 (cond ((eq (pseudo-diamond-left diamond) (face-ref-location site direction1 direction2))
				(setf (pseudo-diamond-data-left-face (aref data number)) face))
			       ((eq (pseudo-diamond-right diamond) (face-ref-location site direction1 direction2))
				(setf (pseudo-diamond-data-right-face (aref data number)) face))))))))
	 site)
	(dotimes (successor 4)		;See who needs to know that this cube was done
	  (when (vv-status-successor-p site successor)
	    (transmit-vv-cell-done successor site)) ;Say that we took care of this cell
	  ;;Now we have up to 2 diamonds to transmit for this site
	  ;;Sometimes they can be transmitted to people who don't intersect the cube.
	  (dotimes (number 2)
	    (let* ((data (aref data number))
		   (diamond (pseudo-diamond-data-diamond data)))
	      (when diamond
		;;We have a pseudo-diamond here because the start point is not in our past.  We need to pass it forward
		;;in some directions so that it will be in the past of some successor eventually.
		(let ((ijkl (xtoi (pseudo-diamond-start diamond))))
		  (loop for successor below 4
			;;There's definitely some direction in which this coordinate is larger than 1.0
			when (> (aref ijkl successor) 1.0) ;Find one such direction
			;;Pass the diamond to that successor.  Eventually it will end up where it needs to go.
			do (return (transmit-pseudo-diamond successor (pseudo-diamond-start diamond)
							    (pseudo-diamond-data-left-face data)
							    (pseudo-diamond-data-right-face data)))))))))))))


(defun write-vv-output ()
  (when *vv-p*				;If using this system
					;If successors (they are all the same) start before initialization time
    (when (> *initial-time* (4vector-t (destination-offset 0)))
      (format t "~&VV output: phases...") (force-output)
      (write-vv-phases)
      (format t "face points...") (force-output)
      (write-vv-face-points)
      (format t "pseudodiamonds+status...") (force-output)
      (write-pseudo-diamonds-and-status)
      (format t "done~%")
      )))

(defun transmit-vv-phase (successor site phase)
  (send-vv-phase (aref *successor-streams* successor) phase site))

(defun transmit-vv-cell-done (successor site)
  (send-vv-cell-done (aref *successor-streams* successor) site))

(defun transmit-vv-face-point (successor site axis1 axis2 position)
  (let ((stream (aref *successor-streams* successor))
	(*send-offset* (destination-offset successor))) ;Offset position into his coordinates
    (send-vv-face-point stream position site axis1 axis2) ;Site will be globalized
    ))

;;Face lists are (site axis1 axis2)
(defun transmit-pseudo-diamond (successor start left-face right-face)
  (let ((stream (aref *successor-streams* successor))
	(*send-offset* (destination-offset successor))) ;Offset position into his coordinates
    ;;It is possible that the faces to which the pseudo-diamond connects have not been transmitted, because the successor only
    ;;needs the face because of perturbations to the start of this diamond.  In that case, we better transmit them anyway
    (mirror-images
     (destructuring-bind (site axis1 axis2) left-face
       (unless (vv-successor-needs-site-p (car left-face) successor :axes (list axis1 axis2))	;If not done before
	 (transmit-vv-face-point successor site axis1 axis2 (face-ref-location site axis1 axis2)))))
    (send-pseudo-diamond stream start left-face right-face)))


;;Hash table (direction string-number diamond-number)-> junction
;;The diamond-number is the number of the (start of the) past diamond
;;DIRECTION = :LEFT if future diamond is to the left of past diamond

(defvar *initial-junctions* nil)

(define-simulate-variable *explicit-initial-directory* nil) ;Where to read initial conditions

(defun explicit-initial-strings-file (&optional (directory (or *explicit-initial-directory* (input-directory))))
  (format nil "~A/initial-strings.dat" directory))

;;Explicit initial conditions.  Starting points of diamonds are given in a file.
(defun setup-explicit ()
  (format t "~&Explicit initial conditions...~%") (force-output)
  (let ((*initializing* t))
    (setq *initial-junctions* (make-hash-table :test #'equal))
    (read-initial-strings)		;Read files from predecessors.
    (format t "Initial strings...") (force-output)
    (read-explicit-initial-conditions))
  (when *check-data-structures* (check-data-structures))
  (format t "~&Explicit initial conditions done.~%") (force-output))

;;Read an integer expressed as a double-float
(defun read-double-fix (stream)
  (let* ((double (read-double stream))
	 (result (round double)))
    (unless (fudge= double (double-float result) fudge-coordinates)
      (error "~S was supposed to be an integer" double))
    result))

(defun read-explicit-initial-conditions ()
  (with-open-file (stream (explicit-initial-strings-file) :element-type '(unsigned-byte 64))
    (format t "Reading ~A..." (pathname stream))
    (loop for string-number from 0
	  for string = (read-explicit-string stream)
	  while string
	  do (process-explicit-string string string-number))))

;;Read one string from datafile, return array of starts in local coordinates.  NIL if EOF
(defun read-explicit-string (stream)
  (let ((count (handler-case (read-double-fix stream) (end-of-file () nil))))
    (and count				;NIL if end of file
	 (let ((starts (make-array count)))
	   (dotimes (index count)
	     (setf (aref starts index) (localize-position (read-4vector stream))))
	   starts))))

(defun process-explicit-string (starts string-number)
  (format t "(~D)" string-number)
  (loop for index below (length starts)
	unless (point-mine (aref starts index)) ;Find a point which is not in our volume, if any
	return (process-explicit-string-1 starts string-number index)
	finally (process-explicit-string-loop starts))) ;All points are in our volume

;;Handle case where some starting points are in our past.  INDEX is one that isn't.
(defun process-explicit-string-1 (starts string-number index)
  (loop with start-index = index
	do (setq index (next-array-index-wrapping index starts))
	until (= index start-index)		 ;Wrapped around, so all segments (if any) have been done.  Exit.
	when (point-mine (aref starts index)) ;Found a point that is ours.  Process segment.
	do (setq index				;Set index to last start that is ours, increment, go on.
		 (process-explicit-open-string starts string-number index))))


;;Install the diamonds of a loop entirely in our volume.  There's no issue of junctions
(defun process-explicit-string-loop (starts)
  (format t "L~D " (length starts))
  (loop with first-diamond = (make-diamond :start (aref starts 0))
	for previous-diamond = first-diamond then diamond
	for index from 1 below (length starts)
	for start = (aref starts index)
	for diamond = (make-diamond :start start)
	do (link-explicit-diamond previous-diamond diamond) ;Connect via directional links
	finally
	(link-explicit-diamond diamond first-diamond) ;Link last one
	(finish-explicit-string first-diamond)))
	
(defun link-explicit-diamond (previous-diamond diamond)
  (if (> (4vector-t (diamond-start diamond)) (4vector-t (diamond-start previous-diamond))) ;We start after previous?
      (setf (diamond-sw diamond) previous-diamond					   ;So he is SW of us
	    (diamond-ne previous-diamond) diamond)
    (setf (diamond-nw diamond) previous-diamond	;We start before previous, so she is NW of us
	  (diamond-se previous-diamond) diamond)))

;;A string which started by crossing from another job
;;START-INDEX is the index of the first diamond that we are to create.
;;Returns index of last start that is inside our volume.
(defun process-explicit-open-string (starts string-number start-index)
  (unless *initial-junctions* (error "Unexpectedly processing open initial string with *initial-junctions* not set"))
  (let* ((first-diamond (make-diamond :start (aref starts start-index)))
	(diamond first-diamond))
    (explicit-link-left-end first-diamond string-number starts start-index) ;Attach to predecessor or make junction
    (loop for count from 1
	  for previous-diamond = diamond
	  for index = start-index then next-index
	  for next-index =  (next-array-index-wrapping index starts) ;Next index to use
	  for start = (aref starts next-index)
	  while (point-mine start)	;Stop when we reach a point that isn't ours
	  do (setq diamond (make-diamond :start start))	      ;Ours.  Make diamond with given starting point
	  do (link-explicit-diamond previous-diamond diamond) ;Connect directional links
	  finally					      ;start is outside our volume.  INDEX is last of ours
	    (explicit-link-right-end diamond string-number starts index) ;Attach to right or make junction
	    (finish-explicit-string first-diamond)
	    (format t "~D " count)
	    (return index))))

(mirror-images
;;Set up end of chain of of explicit diamonds.  DIAMOND is the last diamond in our volume.  STRING-NUMBER and INDEX
;;referr to it.  LEFT is either the next one left from our predecessor or a junction to link to our successor.
(defun explicit-link-left-end (diamond string-number starts index)
  (let ((start (aref starts index))	;Starting point of our diamond
	(left-start (aref starts (west-array-index-wrapping index starts)))) ;Starting point of diamond in other job
    (if (< (4vector-t left-start) (4vector-t start))			     ;Other diamond earlier?
	;;Predecessor should have already sent us the junction, associated with his diamond
	(let ((junction (gethash (list :right string-number (west-array-index-wrapping index starts))
				 *initial-junctions*)))
	  (unless junction
	    (error "Couldn't find predecessor diamond at ~A of string ~D, index ~D" :left string-number index))
	  (let ((left-diamond (junction-left-diamond junction)))
	    (unless (4vector= (diamond-right left-diamond) start fudge-coordinates)
	      (error "~A of ~S from predecessor disagrees with our start ~S" :right left-diamond start))
	    ;;We need to use EQ 4vectors for shared corners.  Which should we use?  It's possible that our left,
	    ;;which is his end, is already shared with the next further left diamond, so use it.  But our start,
	    ;;which is his right, cannot be shared with any other diamond leftward.  It might be shared with another
	    ;;diamond to our right, so use it.	    
	    (setf (diamond-left diamond) (diamond-end left-diamond)
		  (diamond-right left-diamond) (diamond-start diamond)
		  (diamond-ne left-diamond) diamond ;We are NE of predecessor's diamond
		  (diamond-sw diamond) left-diamond))
	  )
      ;;Previous in future. We must make a junction to pass to our successor.
      ;;diamond-start is already filled in
      (setf (diamond-left diamond) left-start ;Fill in left corner since we will not be creating left diamond
	    (diamond-left-rejoining-junction diamond) ;make junction and store
	    (make-initial-junction :string string-number :diamond index :right-diamond diamond)
	    )))))

;;Here with the first diamond that we created in a string.  Add created and tag slots and then
;;add all diamonds to data structures.
(defun finish-explicit-string (diamond)
  (loop with first-diamond = diamond
	for global-start = (globalize-position (diamond-start diamond))
	do (explicit-diamond-end diamond) ;Make sure all points of diamond are filled in
	do (mirror-images (setf (diamond-a-kink-created diamond) global-start))
	do (setf (diamond-tag diamond) (create-loop-tag global-start))
	do (handle-new-diamond diamond)
	do (setq diamond (diamond-e diamond))
	until (or (eq diamond first-diamond) ;Loop
		  (null diamond)	     ;No diamond because it is in the future
		  (diamond-predecessorp diamond)))) ;diamond from processor
	
;;Recursively fill in slots of diamond
(mirror-images
(defun explicit-diamond-left (diamond)
  (or (diamond-left diamond)		;Already known?
      (prog1 (setf (diamond-left diamond)			   ;No.  Install from left predecessor
		   (if (diamond-nw diamond)			   ;Diamonds to left in future?
		       (diamond-start (diamond-nw diamond))	   ;Yes, so her start is our left
		     (explicit-diamond-end (diamond-sw diamond)))) ;No, use her end
	(assert (< (- (4vector-t (diamond-left diamond)) (4vector-t (diamond-start diamond))) 1.0))))))
	     
(defun explicit-diamond-end (diamond)
  (or (diamond-end diamond)		;Already known?
      (setf (diamond-end diamond)	;No, compute from other corners
	    (progn
	      (explicit-diamond-left diamond) ;Make sure other corners known
	      (explicit-diamond-right diamond)
	      (compute-diamond-end diamond)))))

;;Put received initial junction into hash table.  Direction is :left or :right, the direction leading toward the future
(defun handle-initial-junction (junction direction)
  (setf (gethash (list direction (initial-junction-string junction) (initial-junction-diamond junction))
		 *initial-junctions*)
	junction))

(defvar *test-string-explicit-output-overwrite* t)

;;Create initial-strings.dat file.  Set *test-string-explicit-output* to come here
(defun output-explicit-test-string (starts)
  (let ((file (explicit-initial-strings-file *output-directory*)))
    (with-open-file (stream file :element-type '(unsigned-byte 64) :direction :output
			    :if-exists (if *test-string-explicit-output-overwrite* :supersede :error))
      (write-double (double-float (length starts)) stream)
      (loop for start across starts
	    do (write-4vector starts stream)))
    file))


;;Accepts discrete functions a(t-sigma) and b(t+sigma).  These functions are lists of 3-vector values, and the index
;;in the list has no special meaning.  The implicit t+/sigma argument changes by the spatial distance at each step,
;;so the tangent vector automatically has magnitude 1.  The total argument change in the two functions must be the same.
;;We return a list of starting points for a diamond world sheet overlapping *initial-time*
(defun explicit-string-from-discrete-ab (a b)
  (loop with result = (make-array (+ (length a) (length b)))
	with result-index = 0 and a-index = 0 and b-index = 0
	with a-done = nil and b-done = nil
	with first-time = (- *initial-time* (* fudge-coordinates 10)) ;Start slightly earlier than initial surface
	with time = first-time
	do (setf (aref result result-index)
		 (3to4vector (3vector-scale (3vector+ (aref a a-index) (aref b b-index)) 0.5) time)) ;Add this point
	do (let* ((b-next (next-array-index-wrapping b-index b)) ;Tentatively advance b.  Where does that take us?
		  (next-b-time (+ time (/ (3vector-length (3vector- (aref b b-next) (aref b b-index))) 2))))
	     (cond ((< next-b-time *initial-time*) ;If this is in the past of the initial surface, we guessed correctly
		    (when b-done		  ;If b has wrapped, we shouldn't be moving this direction again
		      (error "Overwrapping b when a isn't done"))
		    (setq b-index b-next ;So make tentative advance permanent
			  time next-b-time)
		    (when (zerop b-index) ;b wrapped around?
		      (when a-done (loop-finish)) ;a also, so we're finished
		      (setq b-done t)))	     ;Set flag
		   (t			;No, instead we need to advance backwards in a and into the past
		    (let ((a-previous (previous-array-index-wrapping a-index a)))
		      (decf time (/ (3vector-length (3vector- (aref a a-previous) (aref a a-index))) 2))
		      (setq a-index a-previous)
		      (when a-done		  ;If a has wrapped, we shouldn't be moving this direction again
			(error "Overwrapping a when b isn't done"))
		      (when (zerop a-index)	    ;a wrapped around?
			(when b-done (loop-finish)) ;a also, so we're finished
			(setq a-done t)))))	    ;Set flag
	     (incf result-index))		    ;Fill next slot
	finally
	(unless (fudge= time first-time (* 10 fudge-coordinates))
	  (error "Loop did not have same amount of sigma in the two directions"))
	(unless (= result-index (1- (length result)))
	  (error "Result did not have right number of components"))
	(return result)))
	

;;Make a string from functions for a and b.  MAX-SIGMA is the length of the string.
;;Functions will be called with POINTS values from 0 below MAX-SIGMA.  
;;The argument of a is tau-sigma, meaning that increasing argument goes around the string the opposite way.
;;You should give functions whose tangent vector has magnitude 1.  Nevertheless,
;;the total distance of the piecewise linear approximation will not be max-sigma.
;;So we rescale a and b (keeping the average of each fixed) to make them the right length.
(defun explicit-string-from-ab (a-function b-function max-sigma points &key (offset zero-3vector))
  (mirror-image-let ((delta-sigma (/ max-sigma points))
		     (a-data (make-array points))
		     (a-length 0.0)	;Total length of piecewise-linear a
		     (a-total zero-3vector))	;Total of a positions, for averaging
    (mirror-images
     (loop for index below points	;Collect needed samples of a and b.
	   for previous = nil then a
	   for a = (funcall a-function (* delta-sigma index))
	   do (setf (aref a-data index) a)
	   do (setf a-total (3vector+ a-total a))
	   when previous
	   do (incf a-length (3vector-distance a previous))
	   finally (incf a-length (3vector-distance a (aref a-data 0)))) ;Account last segment
     (loop with scale = (/ max-sigma a-length) ;Some number slightly more than 1 by which the loop needs to be magnified
	   with a-average = (3vector-scale a-total (/ 1.0 points)) ;Average position around which to scale
	   initially (format t "Rescaling ~A by ~F~%" :a scale)
	   for index below points
	   do (setf (aref a-data index)
		    (3vector+ offset a-average (3vector-scale (3vector- (aref a-data index) a-average) scale)))))
     ;;Now we have discrete functions a(tau-sigma) and b(tau+sigma) with the total distance between samples
     ;;max-sigma in both cases.
     (explicit-string-from-discrete-ab a-data b-data)))

(defvar *test-string-explicit-output* nil) ;If set, output test strings instead of actually setting them up

;;Make and install string from functions a and b. See explicit-string-from-ab
;;STRING-NUMBER is used only if the you are initializing for multiple simulation volumes, and the default
;;is fine unless you have more than one string that you're trying to make.
(defun make-test-ab (a-function b-function max-sigma points &key (offset zero-3vector) (string-number 0))
  (let ((starts (explicit-string-from-ab a-function b-function max-sigma points :offset offset)))
    (if *test-string-explicit-output* (output-explicit-test-string starts)
      (process-explicit-string starts string-number))))
