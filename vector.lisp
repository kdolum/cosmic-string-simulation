;;Fast code for manipulation vectors and coordinates
(in-package "CL-USER")

;;In sbcl, arrays have a 2-word header and are allocated on 2-word boundaries.  Thus a
;;3-element array and a 4-element array both take up 6 words.  So there is no space saved
;;by distinguishing.  Meanwhile, sbcl does not know that it doesn't need to check the bounds
;;when it is referencing element 0-2 of an object of type
;;(or (simple-array double-float (3)) (simple-array double-float (4)))
;;Thus it seems to be more efficient to always use 4vectors
(eval-when (:compile-toplevel :load-toplevel :execute)
  (push :long-3vectors *features*))

;;A 3vector is a double-float array (x,y,z).  It could represent either a
;;vector or a position
(defstruct (3vector
	      (:type (vector coordinate))
	      (:constructor nil)	;see below
	      (:copier nil))
    x
    y
    z
    #+long-3vectors t)			;Never used

;;DEFSTRUCT with :TYPE doesn't automatically give a named type.
(deftype strict-3vector () `(simple-array coordinate (#-long-3vectors 3 #+long-3vectors 4)))

;;Needed at compile time for constant vectors
(eval-when (:compile-toplevel :load-toplevel :execute)

  (defun construct-3vector ()
    (make-array #-long-3vectors 3 #+long-3vectors 4 :element-type 'coordinate))

  (defvar 3vectors (make-resource :name "3vector" :constructor #'construct-3vector :max 100)))

;;Make 3vector.  Arguments given positionally.  Supply all or none.
;;If none, the slots are not initialized.
(defmacro make-3vector (&optional x y z)
  (when x (unless z (error "Must supply all slots or none")))
  (let ((result (gensym)))
    `(let ((,result (allocate 3vectors)))
       ,@(when x
	   `((set-3vector ,result ,x ,y ,z)))
       (the strict-3vector ,result))))

;;Set components
(defmacro set-3vector (vector x y z)
  (let ((v (gensym)))
    `(let ((,v ,vector))
       (setf (3vector-x ,v) ,x
	     (3vector-y ,v) ,y
	     (3vector-z ,v) ,z)
       ,v)))				;Return it, like SETF

;;Construct it out of a list
(defun list-3vector (list)
  (let ((result (make-3vector)))
    (loop for i from 0
	  for value in list
	  do (setf (aref result i) value))
    result))

;;Much like 3vector; see above.
(defstruct (4vector
	    (:type (vector coordinate))
	    (:include 3vector)
	    (:constructor nil)
	    (:copier nil))
  #-long-3vectors t)

(defmacro 4vector-time (v) `(4vector-t ,v))

(deftype 4vector () '(simple-array coordinate (4)))

;;Needed at compile time for constant vectors
(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun construct-4vector ()
    (make-array 4 :element-type 'coordinate))

  (defvar 4vectors (make-resource :name "4vector" :constructor #'construct-4vector :max 100)))

(defmacro make-4vector (&optional x y z time)
  (when x (unless time (error "Must supply all slots or none")))
  (let ((result (gensym)))
    `(let ((,result (allocate 4vectors)))
       ,@(when x
	   `((set-4vector ,result ,x ,y ,z ,time)))
       (the 4vector ,result))))

(defmacro set-4vector (vector x y z time)
  (let ((v (gensym)))
    `(let ((,v ,vector))
       (setf (3vector-x ,v) ,x
	     (3vector-y ,v) ,y
	     (3vector-z ,v) ,z
	     (3vector-t ,v) ,time
	     )
       ,v)))


(defun list-4vector (list)
  (let ((result (make-4vector)))
    (loop for i from 0
	  for value in list
	  do (setf (aref result i) value))
    result))

(defun copy-4vector (vector &optional (result (make-4vector)))
  (declare (optimize speed)
	   (type 4vector vector result))
  (loop for i from 0 below 4
	do (setf (aref result i) (aref vector i)))
  result)

;;We can consider a 4vector to be a 3vector
(deftype 3vector () '(or strict-3vector 4vector))

(defun copy-3vector (vector &optional (result (make-3vector)))
  (declare (optimize speed)
	   (type 3vector vector result))
  (loop for i from 0 below 3
	do (setf (aref result i) (aref vector i)))
  result)

(defmacro 3vector-component (vector component)
  `(aref ,vector ,component))
(defmacro 4vector-component (vector component)
  `(aref ,vector ,component))

(eval-when (:compile-toplevel :load-toplevel :execute)
;;Return zeroed vectors.  Needed at compile time for constants.
(defun make-zero-3vector ()
  (make-3vector zero-coordinate zero-coordinate zero-coordinate))
(defun make-zero-4vector ()
  (make-4vector zero-coordinate zero-coordinate zero-coordinate zero-coordinate)))

;;Constants.  Do not use these to initialize slots that could then be modified or deallocated to the resource
(define-constant zero-3vector (make-zero-3vector))
(define-constant zero-4vector (make-zero-4vector))

(define-constant unit-3vector (make-3vector 1.0 1.0 1.0))
(define-constant unit-4vector (make-4vector 1.0 1.0 1.0 1.0))

(defun 3vector-list (x)
  (loop for i below 3
	collect (aref x i)))

(defun 4vector-list (x)
  (loop for i below 4
	collect (aref x i)))

(defun print-3vector (v &optional (stream t))
  (format stream "(~{~$~^, ~})" (3vector-list v)))
  
(defun print-4vector (v &optional (stream t))
  (format stream "(~{~$~^, ~})" (4vector-list v)))

;;See if this is a list of operators of those opposed to a single operator
(defmacro operator-listp (operator)
  `(and (consp ,operator) (not (eq (car ,operator) 'lambda))))

;;General componentwise operations.  We bind a variable to each of the VECTORS and SCALARS.  We apply
;;the function OPERATOR (usually a lambda expression) to each component of each vector and the scalars.
;;We then call RESULT-FUNCTION on the resulting N numbers.
;;RESULT-FUNCTION can also be a list, in which case the N the numbers are appended.  This allows you
;;to supply static arguments to a macro before the slots of the vector are added.
;;OPERATOR can be a 2-element list was first element isn't LAMBDA, in which case the
;;first element is applied to the spatial slots and the last to the time slot
(defmacro nvector-operate-flat (n result-function operator vectors scalars)
  (when (atom result-function) (setq result-function (list result-function)))
  (unless (operator-listp operator) (setq operator (list operator operator)))
  (let ((type (ecase n (3 '3vector) (4 '4vector)))
	(vector-variables (loop repeat (length vectors) collect (gensym)))
	(scalar-variables (loop repeat (length scalars) collect (gensym))))
    `(let (,@(mapcar #'list vector-variables vectors) ;Bind arguments
	   ,@(mapcar #'list scalar-variables scalars)) ;Bind arguments
       (declare (type ,type ,@vector-variables)) ;declare type of bindings
       (,@result-function		;Call given function
	,@(loop for i below n
		for op = (if (= i 3) (second operator) (first operator))
		collect
		`(,op ,@(loop for variable in vector-variables collect `(aref ,variable ,i)) ;for vector, get slot
		      ,@scalar-variables ;Scalars are just themselves
		      ))))))

;;General recursive componentwise operations
;;(3vector+ (3vector-scale a s) (3vector- b c))
;; => (nvector-operate 3 + (3vector- a b) c)
;; => (nvector-operate 3 + (nvector-operate-flat 3 (lambda (g1 g2) (- g1 g2)) a b) c)
;; => (nvector-operate-flat 3 (lambda (g3 g4 g5) (+ ((lambda (g1 g2) (- g1 g2)) g3 g4) g5)) a b c)
(defmacro nvector-operate (n result-function operator vectors scalars &environment environment)
  (unless (operator-listp operator) (setq operator (list operator operator)))
  (loop with scalar-variables = (loop repeat (length scalars) collect (gensym))	;Variables for binding the original scalars
	for vector in vectors
	for (new-flat-vectors new-flat-scalars new-vector-lambdas new-scalar-lambdas argument-s argument-t)
	= (nvector-operate-one-argument n vector environment)
	append new-flat-vectors into flat-vectors
	append new-flat-scalars into flat-scalars
	append new-vector-lambdas into vector-lambda-list
	append new-scalar-lambdas into scalar-lambda-list
	collect argument-s into arguments-s
	collect argument-t into arguments-t
	finally
	(return 
	 `(nvector-operate-flat ,n ,result-function
				,(list `(lambda ,(append vector-lambda-list scalar-lambda-list scalar-variables)
						;;sub operator results in arguments.  Scalars passed on.
					  (,(first operator) ,@arguments-s ,@scalar-variables)) ;Spatial components
				       `(lambda ,(append vector-lambda-list scalar-lambda-list scalar-variables)
					  (,(second operator) ,@arguments-t ,@scalar-variables))) ;time component
				,flat-vectors
				,(append flat-scalars scalars)))))

;;Processing vector argument to nvector-operate.  We return 5 lists:
;;New vector objects to be bound at top level
;;New scalar objects to be bound at top level
;;New vector variables for the lambda list in the new operator
;;New scalar variables for the lambda list in the new operator
;;A new argument to the old operator for spatial parts of the vector
;;A new argument to the old operator for time component of the vector
(eval-when (:compile-toplevel :load-toplevel :execute) ;Needed for macro expansion
(defun nvector-operate-one-argument (n vector environment)
  (declare (ignore n))			;Was used for checking against sub-n, but seems unused now
  (loop
   (when (and (consp vector) (eq (first vector) 'nvector-operate-flat)) ;Recursive call, with sub-recursion expanded
     (destructuring-bind (sub-n sub-result-function sub-operator sub-vectors sub-scalars) (cdr vector)
       (declare (ignore sub-n))
       ;;If sub-n=4, inner code constructs a 4vector.  It's OK if only 3 elements are used,
       ;;in which case we won't construct the fourth at all.
       ;;If sub-n=3, but n=4, inner code constructs a 3vector.  Outer code will then reference fourth element of it,
       ;;which works, but only because 3vector are really 4vectors.  The operator in this case should ignore
       ;;the fourth element.
       (unless (operator-listp sub-operator) (setq sub-operator (list sub-operator sub-operator)))
       (when (member sub-result-function '(make-3vector make-4vector)) ;Sub-call returns a vector result?
	 ;;We can expand what is below us
	 (let ((vector-variables (loop repeat (length sub-vectors) collect (gensym)))
	       (scalar-variables (loop repeat (length sub-scalars) collect (gensym))))
	   (return
	    (list sub-vectors ;These need to be processed at top level
		  sub-scalars
		  vector-variables ;Outer operator takes all these variables now
		  scalar-variables
					;Our operator acts on result of sub-operation
		  `(,(first sub-operator) ,@vector-variables ,@scalar-variables) ;spatial case
		  `(,(second sub-operator) ,@vector-variables ,@scalar-variables) ;time component
		  ))))))
   ;;Doesn't look recursive so far.  Expand it
   (multiple-value-bind (new expanded-p) (macroexpand-1 vector environment)
     (if expanded-p		   ;Expanded, so we should check again
	 (setq vector new)
       ;;Done expanding with no recursion: simple case: vector just becomes argument to operator
       (let ((variable (gensym))) 
	 (return (list (list vector)
		       nil
		       (list variable)
		       nil
		       variable		;Operator just acts on variable
		       variable
		       ))))))))

(defmacro 3vector+ (&rest vectors)
  `(nvector-operate 3 make-3vector + ,vectors nil))

(defmacro 4vector+ (&rest vectors)
  `(nvector-operate 4 make-4vector + ,vectors nil))

;;Increase one vector by another.  The first vector is modified.
(defmacro 3vector-incf (v delta)
  `(nvector-operate 3 (set-3vector ,v) + (,v ,delta) nil))

(defmacro 4vector-incf (v delta)
  `(nvector-operate 4 (set-4vector ,v) + (,v ,delta) nil))

;;Subtract 3vectors, or negate if only one.
(defmacro 3vector- (&rest vectors)
  `(nvector-operate 3 make-3vector - ,vectors nil))

(defmacro 4vector- (&rest vectors)
  `(nvector-operate 4 make-4vector - ,vectors nil))

(defmacro 3vector-decf (v delta)
  `(nvector-operate 3 (set-3vector ,v) - (,v ,delta) nil))

(defmacro 4vector-decf (v delta)
  `(nvector-operate 4 (set-4vector ,v) - (,v ,delta) nil))

;;Multiply vector by scalar
(defmacro 3vector-scale (x scale)
  `(nvector-operate 3 make-3vector * (,x) (,scale)))

;;Multiply vector by scalar
(defmacro 4vector-scale (x scale)
  `(nvector-operate 4 make-4vector * (,x) (,scale)))

;;Euclidian length of vector
(defmacro 3vector-length (x)
  `(real-sqrt (3vector-squared-length ,x)))

(defmacro 3vector-squared-length (x)
  `(nvector-operate 3
		    +				   ;Add results
		    (lambda (part) (* part part))  ;Square each component
		    (,x) nil))

;;Euclidian length of 4vector
(defmacro 4vector-Euclidian-length (x)
  `(real-sqrt
    (nvector-operate 4
		     +			;Add results
		     (lambda (part) (* part part)) ;Square each component
		     (,x) nil)))

;;invariant length with signature (+ + + -)
(defmacro 4vector-squared-length-spacelike (x)
  `(nvector-operate 4
		    +			;Add results
		    ((lambda (part) (* part part)) ;Square each component and multiply by metric
		     (lambda (part) (- (* part part)))) ;Time component negative
		    (,x)
		    nil))

(declaim (inline 4vector-spacelike-length 4vector-timelike-length))
;;The vector must be spacelike.  We return its length
(defun 4vector-spacelike-length (x)
  (let ((length2 (4vector-squared-length-spacelike x))) ;Squared length.  Spacelike is positive.
    (when (minusp length2) (error "Vector ~S is not spacelike" x))
    (real-sqrt length2)))

;;The vector must be timelike.  We return its length
(defun 4vector-timelike-length (x)
  (let ((length2 (4vector-squared-length-spacelike x))) ;Squared length.  Spacelike is positive.
    (when (plusp length2) (error "Vector ~S is not timelike" x))
    (real-sqrt (- length2))))

;;Rescale length to 1 or given length
(defmacro 3vector-normalize (v &optional (length 1.0))
  (let ((vector (gensym)))
    `(let ((,vector ,v))
       (3vector-scale ,vector (/ ,length (3vector-length ,vector))))))

;;Rescale vector so time component is 1 or given value
(defmacro 4vector-normalize-time (v &optional (time 1.0))
  (once-only (v)
    `(4vector-scale ,v (/ ,time (4vector-t ,v)))))

;;Euclidian distance between vectors.
(defmacro 3vector-distance (x1 x2)
  `(3vector-length (3vector- ,x1 ,x2)))

;;Euclidian distance between 4vectors.
(defmacro 4vector-Euclidian-distance (x1 x2)
  `(4vector-Euclidian-length (4vector- ,x1 ,x2)))

(defmacro 4vector-squared-distance-spacelike (x1 x2)
  `(4vector-squared-length-spacelike (4vector- ,x1 ,x2)))

;;dot product
(defmacro 3vector-dot (v1 v2)
  `(nvector-operate 3
		    +			;Add results
		    *			;Multiply components
		    (,v1 ,v2) nil))


;;scalar product with signature (+ + + -)
(defmacro 4vector-dot-spacelike (v1 v2)
  `(nvector-operate 4 + (* (lambda (x y) (- (* x y)))) (,v1 ,v2) nil))
		     
;;3-d cross product
(defun 3vector-cross-product (v1 v2)
  (make-3vector
   (- (* (3vector-y v1) (3vector-z v2)) (* (3vector-z v1) (3vector-y v2)))
   (- (* (3vector-z v1) (3vector-x v2)) (* (3vector-x v1) (3vector-z v2)))
   (- (* (3vector-x v1) (3vector-y v2)) (* (3vector-y v1) (3vector-x v2)))))

;;See if vectors are the same within fudge factor
(defun 3vector= (v1 v2 fudge-factor)
  (loop for index below 3
	always (fudge= (aref v1 index) (aref v2 index) fudge-factor)))

;; decides if two 4-vectors are parallel by first converting them to 4-velocities and using 4-vector=
;; fudge-factor is rescaled so that tiny 4-vectors are parallel to everything.
(defun 4vector-parallel-p (v1 v2 fudge-factor)
  (let ((min-t (min (4vector-time v1) (4vector-time v2)))
	(vel1 (4vector-scale v1 (/ 1.0 (4vector-time v1))))
	(vel2 (4vector-scale v2 (/ 1.0 (4vector-time v2)))))
    ;(assert (plusp min-t))
    (4vector= vel1 vel2 (/ fudge-factor min-t))))

;;Test for exact equality
(declaim (inline 4vector-exactly=))
(defun 4vector-exactly= (v1 v2)
  (declare (type 4vector v1 v2))
  (loop for index below 4
	always (= (aref v1 index) (aref v2 index))))

;;Turn 3vector into 4vector by adding time
(defmacro 3to4vector (v time)
  `(nvector-operate 4 make-4vector ((lambda (slot time) (declare (ignore time)) slot) ;Use slot of 3vector
				    (lambda (slot time) (declare (ignore slot)) time))
		    (,v) (,time)))

;;3vector operations give 0.0 for the time component
(declaim (inline zero-time-component))
(defun zero-time-component (&rest arguments)
  (declare (ignore arguments))
  0.0)

;;Convert a 4vector to 3vector.  When applied to an input 4vector, it makes a new vector with time = 0.
;;When part of a chain of operations, the compiler will optimize away the computation of the fourth slot.
(defmacro 4to3vector (p)
  `(nvector-operate 3 make-3vector ((lambda (x) x) zero-time-component) (,p) nil))

;;See if vectors are the same up to periodicity, within fudge factor
(defun 3vector-equivalent (v1 v2 &optional (fudge-factor fudge-global-coordinates))
  (loop for dual-vector in (periodicity-dual-basis)
	always (fudge= (mod (+ (3vector-dot (3vector- v1 v2)
					    dual-vector) ;should be an integer
			       0.5)	;should be an integer plus 0.5
			    1.0)	;should be 0.5
		       0.5
		       fudge-factor)))


(defun 4vector= (v1 v2 fudge-factor)
  (loop for index below 4
	always (fudge= (aref v1 index) (aref v2 index) fudge-factor)))

(defun 4vector-equivalent (v1 v2 &optional (fudge-factor fudge-global-coordinates))
  (and (fudge= (4vector-t v1) (4vector-t v2) fudge-factor)
       (3vector-equivalent (4to3vector v1) (4to3vector v2) fudge-factor)))

;;Evenly distributed random unit vector
(defun random-direction ()
  (let* ((costheta (- (random 2.0) 1.0))
	 (sintheta (sqrt (- 1.0 (expt costheta 2))))
	 (phi (random (* 2 pi))))
    (make-3vector (* sintheta (sin phi)) (* sintheta (cos phi)) costheta)))

;;Copy a simple vector of 3vectors
(defun copy-3vector-array (array)
  (let ((copy (make-array (length array))))
    (dotimes (index (length array))
      (setf (aref copy index) (copy-3vector (aref array index))))
    copy))