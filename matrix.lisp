;;Some matrix manipulation routines
(in-package "CL-USER")
;;;Copyright (c) 1998 Ken Olum.

;;Not using speeded-up code
(defmacro float-aref (&rest args)
  (aref args))

(defmacro with-fast-float-arrays (specs &body body)
  (declare (ignore specs))
  `(progn ,@body))

#|

;;Optimize calculations with double-floats.  It does not work to declare
;;the type of the array unless it is a simple-vector
(defmacro float-aref (&rest args)
  `(the float (aref ,@args)))

(eval-when (compile load eval)
  (defvar *compiler-fast-float-info* nil
    "Alist of actual-var -> fast-vector-var, fast-dims-var"))

;;Fancy optimization for arrays of double-floats.  We get the underlying
;;array for multidimensional arrays and get it at it with array-row-major-index
;;This triggers lots of compiler optimizations.  You don't need float-aref
;;if you are using this.
;;Currently this works only with 1 and 2-d arrays
(defmacro with-fast-float-array (var &body body)
  (let ((fast-vec-var (make-symbol (format nil "fast-vec-~A" var)))
	(fast-dims-var (make-symbol (format nil "fast-dims-~A" var))))
    ;;dimensions are all the dimensions after the first
    `(compiler-let ((*compiler-fast-float-info*
		    (cons (list ',var ',fast-vec-var ',fast-dims-var)
			  *compiler-fast-float-info*)))
       (let ((,fast-vec-var (system:underlying-simple-vector ,var))
	     (,fast-dims-var (cdr (array-dimensions ,var))))
	 (declare (type system:simple-double-float-vector ,fast-vec-var))
	 ,@body))))

(defmacro fast-aref (array &rest indicies)
  (let ((spec (assq array *compiler-fast-float-info*)))
    (unless spec
      (error "~S has not be declared with ~S" array 'with-fast-float-array))
    (let ((vec-var (second spec))
	  (dims-var (third spec)))
      `(aref ,vec-var (fast-row-major-index ,dims-var ,indicies)))))

;;Does the computation for array-row-major-index doing as much as possible
;;in advance.  DIMS is a list of the array-dimensions not counting the
;;first.
;;The idea is that
;; (A) -> A
;; (A B) -> (+ (* (FIRST DIMS) A) B)
;; (A B C) -> (+ (* (SECOND DIMS) (+ (* (FIRST DIMS) A) B)) C)
;; etc.
(defmacro fast-row-major-index (dims indicies)
  (fast-row-major-index-1 dims (reverse indicies)))

(defun fast-row-major-index-1 (dims rindicies)
  (cond ((null (cdr rindicies))		;First index
	 (car rindicies))		;Just return it
	(t
	 `(+ (*				;Do preceding indicies
	      ,(fast-row-major-index-1 dims (cdr rindicies)) 
	      (nth ,(- (length rindicies) 2) ,dims)) ;Multiply by next dim
	     ,(car rindicies)
	     ))))

(defmacro with-fast-float-arrays (specs &body body)
  (if specs				;Anything to do?
      `(with-fast-float-array ,(car specs)
          (with-fast-float-arrays ,(cdr specs) ,@body))
    `(progn ,@body)))

|#


(defun print-matrix (mat &key mathematica (format "~7F")
			 (start-row 0) (end-row (array-dimension mat 0))
			 (start-col 0) (end-col (array-dimension mat 1)))
  (cond (mathematica
	 (format t "~&{")
	 (loop for row from start-row below end-row do
	   (format t "{")
	   (loop for col from start-col below end-col do
	     (format t format (aref mat row col))
	     (unless (= col (1- (array-dimension mat 1)))
	       (format t ", ")))
	   (format t "}")
	   (unless (= row (1- (array-dimension mat 0)))
	     (format t ",~%"))))
	(t
	 (loop for row from start-row below end-row do
	   (format t "~&")
	   (loop for col from start-col below end-col do
	     (format t format (aref mat row col))
	     (format t " "))))))

(defun make-identity (n)
  (let ((mat (make-array (list n n) :initial-element 0)))
    (dotimes (i n)
      (setf (aref mat i i) 1))
    mat))

(defun make-identity-float (n)
  (let ((mat (make-array (list n n) :element-type 'double-float :initial-element 0.0)))
    (dotimes (i n)
      (setf (aref mat i i) 1.0))
    mat))

(defun make-diagonal (elts)
  (let* ((n (length elts))
	 (mat (make-array (list n n) :initial-element 0)))
    (loop for i from 0
	  for elt in elts
	  do (setf (aref mat i i) elt))
    mat))

(defun transpose (mat)
  (let ((nmat (make-array (reverse (array-dimensions mat))
                          :element-type (array-element-type mat))))
    (dotimes (i (array-dimension mat 0))
      (dotimes (j (array-dimension mat 1))
        (setf (aref nmat j i) (aref mat i j))))
    nmat))

#|
(defun float-transpose (mat)
  (check-type mat (array float))
  (let ((nmat (make-array (reverse (array-dimensions mat))
                          :element-type 'float)))
    (with-fast-float-arrays (mat nmat)
      (dotimes (i (array-dimension mat 0))
        (dotimes (j (array-dimension mat 1))
          (setf (fast-aref nmat j i) (fast-aref mat i j))))
      nmat)))
|#

;;Multiply matricies or vectors
(defun dot (&rest mats)
  ;;This is rather a kluge.  We accept a list as a one-dimensional array
  (loop for cons on mats
	for mat = (car cons)
	when (listp mat)
	  do (setf (car cons) (make-array (length mat) :initial-contents mat)))
  (let ((mat1 (first mats))
        (mat2 (second mats))
        (more (cddr mats)))
    (let ((res 
	   (ecase (array-rank mat1)
	     (2 (ecase (array-rank mat2)
		  (2 (dotmm mat1 mat2))
		  (1 (dotmv mat1 mat2))))
	     (1 (ecase (array-rank mat2)
		  (2 (dotvm mat1 mat2))
		  (1 (dotvv mat1 mat2)))))))
      (if more                          ;Recurse to do rest of matricies
	  (apply #'dot res more)
	res))))

(defun dotmm (mat1 mat2)
  ;;We are multiplying an LxM matrix by an MxN matrix to get an LxN matrix
  (let ((l (array-dimension mat1 0))
	(m (array-dimension mat1 1))
	(n (array-dimension mat2 1)))
    (unless (= m (array-dimension mat2 0))
      (error "Matricies for dot product aren't compatible in size"))
    (let ((new-mat (make-array (list l n)
			       :element-type (array-element-type mat1))))
      (dotimes (i l)
	(dotimes (k n)
	  (setf (aref new-mat i k)
		(loop for j below m
		      sum (* (aref mat1 i j) (aref mat2 j k))))))
      new-mat)))

;;Matrix times vector 
(defun dotmv (mat vec)
  (let ((l (array-dimension mat 0))
	(m (array-dimension mat 1)))
    (unless (= m (length vec))
      (error "Matrix and vector for dot product aren't compatible in size"))
    (let ((new-vec (make-array l :element-type (array-element-type mat))))
      (dotimes (i l)
	(setf (aref new-vec i)
	      (loop for j below m
		    sum (* (aref mat i j) (aref vec j)))))
      new-vec)))

;;Vector times matrix.  Vector is taken as a row vector.
(defun dotvm (vec mat)
  (let ((l (array-dimension mat 0))
	(m (array-dimension mat 1)))
    (unless (= l (length vec))
      (error "Matrix and vector for dot product aren't compatible in size"))
    (let ((new-vec (make-array m :element-type (array-element-type vec))))
      (dotimes (j m)
	(setf (aref new-vec j)
	      (loop for i below l
		    sum (* (aref vec i) (aref mat i j)))))
      new-vec)))

;;Inner product.
(defun dotvv (vec1 vec2)
  (let ((l (length vec1)))
    (unless (= (length vec2) l)
      (error "Vectors to be dotted aren't equal in length"))
    (loop for i below l
	  sum (* (aref vec1 i) (aref vec2 i)))))

;;Multiple array or list by scalar
;;Kluge
(defun smult (scalar array)
  (etypecase array
    (list
     (mapcar #'(lambda (x) (* scalar x)) array))
    (array
     (let ((result (make-array (array-dimensions array)
			       :element-type (array-element-type array))))
       (ecase (array-rank array)
	 (1
	  (dotimes (i (array-dimension array 0))
	    (setf (aref result i) (* (aref array i) scalar))))
	 (2
	  (dotimes (i (array-dimension array 0))
	    (dotimes (j (array-dimension array 1))
	      (setf (aref result i j) (* (aref array i j) scalar))))))
       result))))

;; 3-d cross product
(defun vector-cross (v1 v2)
  (destructuring-bind (x1 y1 z1) v1
    (destructuring-bind (x2 y2 z2) v2
      (list (- (* y1 z2) (* z1 y2))
	    (- (* z1 x2) (* x1 z2))
	    (- (* x1 y2) (* y1 x2))))))

;; Returns a copy of a matrix with the same element type.
(defun copy-matrix (array)
  (let* ((rows (array-dimension array 0))
	 (cols (array-dimension array 1))
	 (new-ar (make-array (list rows cols)
			     :element-type (array-element-type array))))
    (dotimes (i rows)
      (dotimes (j cols)
	(setf (aref new-ar i j)
	      (aref array i j))))
    new-ar))

(defun matrix+ (&rest matricies)
  (let* ((mat1 (first matricies))
	 (shape (array-dimensions mat1))
	 (rows (first shape))
	 (cols (second shape)))
    (dolist (mat (cdr matricies))
      (unless (equal shape (array-dimensions mat))
	(error "Matricies don't have the same shape and size")))
    (let ((result (make-array shape :element-type (array-element-type mat1))))
      (dotimes (i rows)
	(dotimes (j cols)
	  (setf (aref result i j)
		(loop for mat in matricies
		      sum (aref mat i j)))))
      result)))

(defun matrix- (mat1 &rest more)
  (let* ((shape (array-dimensions mat1))
	 (rows (first shape))
	 (cols (second shape)))
    (dolist (mat more)
      (unless (equal shape (array-dimensions mat))
	(error "Matricies don't have the same shape and size")))
    (let ((result (make-array shape :element-type (array-element-type mat1))))
      (dotimes (i rows)
	(dotimes (j cols)
	  (setf (aref result i j)
		(- (aref mat1 i j)
		   (loop for mat in more
			 sum (aref mat i j))))))
      result)))

;;Random matricies for testing
(defun make-random-float-matrix (rows cols &optional symmetric-p (type 'float))
  (when symmetric-p
    (unless (= rows cols)
      (error "Symmetric matricies must be square")))
  (let ((res (make-array (list rows cols) :element-type type))
	(one (coerce 1 type)))
    (dotimes (i rows)
      (dotimes (j (if symmetric-p (1+ i) cols))
	(let ((val (random one)))
	  (setf (aref res i j) val)
	  (when symmetric-p
	    (setf (aref res j i) val)))))
    res))


;;Code specialized to simple double-float arrays

;;Resource of 4x4 double-float matrices
(defvar 4x4
 (make-resource :name "4x4" :constructor #'(lambda () (make-array '(4 4) :element-type 'double-float)) :max 10))

;;Copy a 4x4 matrix of double-floats
(defun copy-matrix-4-4 (array &optional (new-array (allocate 4x4)))
  (declare (optimize speed
		     (safety 0))	;Don't type-check arguments
	   (type (simple-array double-float (4 4)) array new-array))
  (dotimes (i 4)
    (dotimes (j 4)
      (setf (aref new-array i j) (aref array i j))))
  new-array)

;;Multiply double-float 4x4 matrix by vector
(defun dotmv44 (m v)
  (declare (optimize speed
		     (safety 0))	;Don't type-check arguments
	   (type (simple-array double-float (4 4)) m)
	   (type 4vector v))
  (let ((result (allocate 4vectors)))
    (declare (type 4vector result))
    (dotimes (i 4)
      (setf (aref result i)
	    (loop for j below 4 sum (* (aref m i j) (aref v j)) double-float)))
    result))

;;Multiply double-float 3x3 matrix by 3vector.  dotmv does not handle this case because 3vectors are actually 4
;;elements long.
(defun dotmv33 (m v)
  (declare (optimize speed)
	   (type (simple-array double-float (3 3)) m)
	   (type 4vector v))
  (let ((result (allocate 3vectors)))
    (declare (type 3vector result))
    (dotimes (i 3)
      (setf (aref result i)
	    (loop for j below 3 sum (* (aref m i j) (aref v j)) double-float)))
    result))

;;Solve m.x = b
(defun linear-solve (m b)
  (let ((lu (copy-matrix m))		;They will be clobbered
	(result (copy-seq b))
	(vv (make-array (length b) :element-type 'double-float))
	(index (make-array (length b) :element-type '(unsigned-byte 16))))
    (ludcmp lu index vv)		;LU decompose in place
    (lubksb lu index result)		;Substitute in place
    result))
