;;;Multidimensional linear interpolating functions for caching slow-to-compute results
(in-package "CL-USER")

;;Data for a cube starting at START (a vector) and going by the elements of INTERVALS in each direction.
(defstruct (interpolation-data (:print-object print-interpolation-data))
  (start nil :type (vector double-float))
  (intervals nil :type (vector double-float))
  (table nil :type (array double-float))					;Array of values indexed by step number list
  )

(defun print-interpolation-data (data stream)
  (print-unreadable-object (data stream :type t)
    (format stream "start ~S, intervals ~S, counts ~S"
	    (interpolation-data-start data)
	    (interpolation-data-intervals data)
	    (array-dimensions (interpolation-data-table data)))))

;;Create interpolation data by calling the function over cube.  START is the starting location,
;;INTERVALS gives the deltas in the various directions and COUNTS the number of steps taken each direction
;;FUNCTION takes N arguments and returns the value.
;;Arguments can be any sequences.
(defun make-interpolation (function start intervals counts)
  (let ((table (make-array (coerce counts 'list) :element-type 'double-float)))
    (setq intervals (coerce intervals '(vector double-float)))
    (setq start (coerce start '(vector double-float)))
    (setq counts (coerce counts 'vector))
    (fill-interpolation-table function table start intervals counts 0 0)
    (make-interpolation-data :start start :intervals intervals :table table)))
	  
(defun fill-interpolation-table (function table start intervals counts row-major-start axis)
  (if (= axis (length start))
      (setf (row-major-aref table row-major-start) (apply function (coerce start 'list)))
    (let ((our-start (copy-seq start))
	  (stride (reduce #'* counts :start (1+ axis))))
      (dotimes (step (aref counts axis))
	(setf (aref our-start axis) (+ (aref start axis) (* (aref intervals axis) step)))
	(fill-interpolation-table function table our-start intervals counts
				  (+ row-major-start (* stride step)) (1+ axis))))))
    
;;Create and read an array of double-floats in row major order
(defun read-double-float-array (stream &rest dimensions)
  (let ((array (make-array dimensions :element-type 'double-float)))
    (dotimes (index (array-total-size array))
      (setf (row-major-aref array index) (read-double stream)))
    array))

;;Write an array of double-floats  in row major order without regard to shape.
(defun write-double-float-array (array stream)
  (dotimes (index (array-total-size array))
    (write-double (row-major-aref array index) stream)))

;;Read from file.  Format is DIMENSIONS, START, INTERVALS, COUNTS, data in row major order
;;Everything is a double float
(defun read-interpolation-data (filename)
  (with-open-file (stream filename :element-type '(unsigned-byte 64))
    (read-interpolation-data-stream stream)))

(defun read-interpolation-data-stream (stream)
  (let ((dimensions (round (read-double stream))))
    (let* ((start (read-double-float-array stream dimensions))
	   (intervals (read-double-float-array stream dimensions))
	   (counts (read-double-float-array stream dimensions))
	   (table (apply #'read-double-float-array stream (map 'list #'round counts))))
      (make-interpolation-data :start start :intervals intervals :table table))))
  
(defun write-interpolation-data (data filename)
  (with-open-file (stream filename :element-type '(unsigned-byte 64) :direction :output :if-exists :supersede)
    (let ((table (interpolation-data-table data)))
      (write-double (double-float (array-rank table)) stream)
      (write-double-float-array (interpolation-data-start data) stream)
      (write-double-float-array (interpolation-data-intervals data) stream)
      (write-double-float-array (map '(vector double-float) #'double-float (array-dimensions table)) stream)
      (write-double-float-array table stream))))

(defun interpolate (data &rest arguments)
  (do-interpolate data arguments 0 0))

(defun do-interpolate (data arguments row-major-location axis)
  (let ((table (interpolation-data-table data)))
    (if (= axis (array-rank table))								;No axes left?
	(row-major-aref table row-major-location)						;Nothing to interpolate
      (multiple-value-bind (index rest)
	  (floor (/ (- (nth axis arguments) (aref (interpolation-data-start data) axis)) ;Index in this direction
		    (aref (interpolation-data-intervals data) axis)))
	(cond ((>= (1+ index) (array-dimension table axis)) ;Must be able to get subsequent index
	       (error "Attempted to interpolate beyond end of data"))
	      ((< index 0)
	       (error "Attempted to interpolate before beginning of data")))
	(let* ((stride (reduce #'* (array-dimensions table) :start (1+ axis))))
	  (+ (* (do-interpolate data arguments (+ row-major-location (* stride index)) (1+ axis)) (- 1.0 rest))
	     (* (do-interpolate data arguments (+ row-major-location (* stride (1+ index))) (1+ axis)) rest)))))))

;;Simpler interface for one dimension
(defun make-interpolation-1 (function start interval count)
  (make-interpolation
   function
   (make-array 1 :element-type 'double-float :initial-element start)
   (make-array 1 :element-type 'double-float :initial-element interval)
   (make-array 1 :initial-element count)))
