;;;Divide simulation volume into cells.  Keep track of which cells a diamond
;;;overlaps and which diamonds overlap each cell
(in-package "CL-USER")

;;Default size of cells if we're not decreasing them because of expansion.
(defvar cell-goal-size 0.5)

(defvar *sparse-cells* nil)		;Use sparse array for cells

(defvar *cells*)			;3-dimensional array of cell contents, or hash table
(defvar *cells-flat*)			;Actual data stored in this one-dimensional array

(defvar *cell-size*)			;Size of each cell
(defvar *cells-start-coordinate*)	;Minimum coordinate value covered by array
(defvar *cells-dimension*)		;Dimension in each direction of array

(declaim (double-float *cell-size* *cells-start-coordinate*)
	 (type (unsigned-byte 20) *cells-dimension*) ;Small enough that compiler can see overflows will never occur
	 (type (simple-array list (*)) *cells-flat*))

(declaim (inline cells-ref)
	 (ftype (function (fixnum) list) cells-ref))
(defun cells-ref (index)
  (if *sparse-cells*
      (gethash index *cells*)
    (aref *cells-flat* index)))

(defsetf cells-ref (index) (value)
  `(if *sparse-cells*
       (if ,value
	   (setf (gethash ,index *cells*) ,value)
	 (remhash ,index *cells*))
     (setf (aref *cells-flat* ,index) ,value)))

;;Initialize cells system.
(defun initialize-cells ()
  (let ((distance (and *size* (job-coordinate-range))) ;Distance that we need to cover
	(goal-size cell-goal-size))	;Start with specified goal size
    (when *longest-edge*	;If we know the set of diamonds that we have, use half of largest edge
      (setq goal-size (/ *longest-edge* 2)))
    (if *size*
	(setq *cells-start-coordinate* (job-coordinate-minimum)
	      *cells-dimension* (ceiling distance goal-size) ;Number of cells in each direction
	      *cell-size* (/ distance *cells-dimension*) ;Size of each cell.
	      *sparse-cells* (> *cells-dimension* 100)) ;If more than a million, store sparsely
      ;;For infinite volume, just use one giant cell
      (setq *cells-start-coordinate* -5d19 *cells-dimension* 1 *cell-size* 1d20 *sparse-cells* nil)))
  (format t "~D^3 cells of size ~S.  " *cells-dimension* *cell-size*)
  (if *sparse-cells*
      (setq *cells* (make-hash-table :test #'eq))
    (setq *cells-flat* (make-array (expt *cells-dimension* 3) :initial-element nil) ;1-dim array to store actual data
	  *cells* (make-array (list *cells-dimension* *cells-dimension* *cells-dimension*) :displaced-to *cells-flat*))))

;;Evaluate body with index bound successively to all cells that diamond intersects
(defmacro do-diamond-cells ((index diamond) &body body)
  `(if *size*
       (do-diamond-cells-1 (,index ,diamond) ,@body)
     (let ((,index 0)) ;If infinite volume, only one cell
       ,@body)))

;;Return first and last indices to cover given coordinate range in given direction
(declaim (inline diamond-cells-index-range))
(defun diamond-cells-index-range (minima maxima axis)
  (declare (type 3vector minima maxima))
  (values (max 0 			;If it extends outside our space, use starting box
	       (fixnum-floor (/ (- (aref minima axis) *cells-start-coordinate* fudge-coordinates)
				*cell-size*))) ;earliest possible slot
	  (min (1- *cells-dimension*)	;Largest possible coordinate
	       (fixnum-floor (/ (+ (- (aref maxima axis)  *cells-start-coordinate*) fudge-coordinates)
				*cell-size*))) ;latest possible slot
	  ))

(defmacro do-diamond-cells-1 ((index diamond) &body body)
  (let ((minima (gensym)) (maxima (gensym))
	(x (gensym)) (xmin (gensym)) (xmax (gensym))
	(y (gensym)) (ymin (gensym)) (ymax (gensym))
	(z (gensym)) (zmin (gensym)) (zmax (gensym)))
    `(let ((,minima (diamond-minima ,diamond))
	   (,maxima (diamond-maxima ,diamond))
	   (,index 0))
       (declare (fixnum index))
       (multiple-value-bind (,xmin ,xmax) (diamond-cells-index-range ,minima ,maxima 0)
	 (multiple-value-bind (,ymin ,ymax) (diamond-cells-index-range ,minima ,maxima 1)
	   (multiple-value-bind (,zmin ,zmax) (diamond-cells-index-range ,minima ,maxima 2)
	     (declare (type (unsigned-byte 20) ,xmin ,xmax ,ymin ,ymax ,zmin ,zmax))
	     (loop for ,x from ,xmin to ,xmax
		   do (loop for ,y from ,ymin to ,ymax
			    do (loop for ,z from ,zmin to ,zmax
				     do (locally (declare (optimize (safety 0))) ;Don't worry about fixnum overflow
					  (setq index (+ ,z (* *cells-dimension* (+ ,y (* *cells-dimension* ,x))))))
				     do (locally ,@body))))))))))
	       
;;Store diamond in cell structures
(defun add-diamond-cells (diamond)
  (do-diamond-cells (index diamond)
    (push diamond (cells-ref index))))

;;Delete diamond from cell structures
(defun delete-diamond-cells (diamond)
  (declare (optimize speed))
  (do-diamond-cells (index diamond)
    (let ((list (cells-ref index)))
      (if (eq (car list) diamond)	;First element?
	  (setf (cells-ref index) (cdr list))
	(loop for rest = (cdr list)
	      unless rest
	      do (error "~D not found in cell ~D" diamond index)
	      when (eq (car rest) diamond) ;Found it
	      do (setf (cdr list) (cdr rest)) ;splice out
	      and return nil
	      do (setq list rest))))))

	

(defun maybe-delete-diamond-cells (diamond)
  (unless (or *read-and-store*
	      (diamond-inertp diamond))	;If inert, should not be in cells
    (delete-diamond-cells diamond)))

;;Call function for each diamond in any cell overlapped by the given diamond.
;;If you set the processed bit of any diamonds before or during the execution of this function,
;;we will not call your function with that diamond.  We will clear all the bits of the diamonds that we find
;;but if you set any others, you are responsible for clearing them.
;;Do not change the cell structure during the execution of this function.
(defun map-diamond-cells-diamonds (function diamond add)
  (declare (optimize speed)
	   (function function)
	   (diamond diamond))
  (using-diamond-processed
   (declare (optimize (safety 0)))	;Don't check declarations
   (unwind-protect
       (do-diamond-cells (index diamond)
	 (dolist (d (cells-ref index))
	   (declare (diamond d))
	   (unless (diamond-processedp d)   ;Not seen yet
	     (setf (diamond-processedp d) t) ;Set flag
	     (funcall function d))	     ;and call function
	   )
	 (when add (push diamond (cells-ref index)))) ;Now add to cell
     (do-diamond-cells (index diamond)
       (dolist (d (cells-ref index))
	 (declare (diamond d))
	 (setf (diamond-processedp d) nil))) ;Clear flag for next time
     )))

;;Evaluate body with D bound successively to all diamonds that are in cells DIAMOND intersects
;;and add DIAMONDS to the cells
;;See MAP-DIAMOND-CELLS-DIAMONDS
(defmacro do-diamond-cells-diamonds ((d diamond &key add) &body body)
  (let ((function (gensym)))
    `(flet ((,function (,d) ,@body))
       (declare (dynamic-extent (function ,function))) ;stack-cons closure
       (map-diamond-cells-diamonds #',function ,diamond ,add))))

;;;Debugging
;;Decide if a diamond is ours to advance, but not in the calendar
(defun diamond-abandoned-p (diamond)
  (let ((end (diamond-end diamond)))
    (if (or (object-in-calendar-p *calendar* (4vector-t end) diamond);in calendar?
	    (not (point-mine end)))	;or not ours?
	nil				;then not abandoned
      t)))

;;Return a list of lines as 2-element lists of 3-element lists giving the outline of the a bounding box
(defun bounding-box-outline (minima maxima)
  (loop for index1 below 3
	nconc (loop for index2 below index1 ;All pairs of 2 directions
		    for index3 = (- 3 index1 index2)
		    nconc (loop for value1 below 2
				nconc (loop for value2 below 2
					    collect
					    (loop for value3 below 2
						  for result = (make-list 3)
						  do (setf (nth index1 result) (aref (if (plusp value1) maxima minima) index1)
							   (nth index2 result) (aref (if (plusp value2) maxima minima) index2)
							   (nth index3 result) (aref (if (plusp value3) maxima minima) index3))
						  collect result
						  ))))))

(defun plot-bounding-boxes (list &rest gnuplot-keys)
  (let ((data (mapcar #'(lambda (diamond)
			  (multiple-value-call #'bounding-box-outline (diamond-bounding-box diamond)))
		      list)))
    (apply #'gnuplot (length list) 48
	   #'(lambda (plot point)
	       (unless (eq point :title)
		 (multiple-value-bind (edge slot) (truncate point 4)
		   (values-list (nth slot (nth edge (nth plot data)))))))
	   :styles :lines
	   :3d t
	   gnuplot-keys)))


(defun diamond-cells-data ()
  (let ((table (make-hash-table :test #'eq))
	(maximum 0)
	(total 0)
	(total-square 0)
	max-index)
    (flet ((do-entry (index entry)
	     (let ((length (length entry)))
	       (dolist (diamond entry) (setf (gethash diamond table) t))
	       (incf total length)
	       (incf total-square (expt length 2))
	       (when (> length maximum)
		 (setq maximum length max-index index)))))
      (if *sparse-cells*
	  (maphash #'do-entry *cells*)
	(dotimes (index (length *cells-flat*))
	  (do-entry index (cells-ref index)))))
    (if (plusp total)
	(format t "There are ~D diamonds in ~D entries, the longest list (index ~D) has ~D, and the average entry is in a list of length ~$~%"
		(hash-table-count table) total max-index maximum (/ (float total-square) total))
      (format t "There are no diamonds ~%"))))
