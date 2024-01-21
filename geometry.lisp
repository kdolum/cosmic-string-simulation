(in-package "CL-USER")

(defvar *ijkl-origin*)			;The location in global xyzt coordinates of ijkl = 0.
(defvar *xtoi* nil)			;Matrix to convert ijkl to xyzt
(defvar *itox* nil)			;Matrix to convert xyzt to ijkl
(declaim (double-float *job-end*))
(defvar *job-end*)

;;Vectors pointing from our space to successors  The fifth successor is dumping, which we do in global coordinates
(defvar *successor-offsets* (make-array 5))

(declaim (ftype (function (4vector) 4vector) xtoi itox)
	 (inline itox xtoi))

(without-compiler-notes			;Avoid complaints about function type declaration to complex to check

;;Convert a 4vector in local xyzt coordinates to local ijkl
(defun xtoi (xyzt)
  (dotmv44 *xtoi* xyzt))

;;Convert a 4vector in local ijkl coordinates to local xyzt coordinates
(defun itox (ijkl)
  (dotmv44 *itox* ijkl))
)					;without-compiler-notes

;;Set up the coordinate transformations.  Each simulation job volume will be
;;a 4-cube whose main diagonal points in the time direction.
;;SIZE specifies the size of our cube as the periodicity distance.
;;TOTAL-SIZE is the periodicity distance of the overall lattice
;;IJKL gives the starting position of this job in units of the edge length of the job volumes.
;;JOB-NUMBER gives the number of the string, from which we compute IJL
;;In either case, we set *GLOBAL-IJKL* and *JOB-NUMBER*.
;;IJKL-ORIGIN specifies the global XYZT coordinates of IJKL = 0.
(defun setup-geometry (&key total-size (split-factor 1) ijkl (ijkl-origin zero-4vector) job-number)
  (setup-split-factor split-factor)	;Setup data structures for job numbering
  (setq *total-size* total-size
	*size* (/ total-size split-factor)
	*ijkl-origin* ijkl-origin)
  (unless ijkl				;Position not specified
    (setq ijkl
	  (if job-number		   ;Did caller give job number?
	      (job-coordinates job-number) ;Infer position
	    zero-4vector)))		  ;Nothing supplied: use zero
  ;;This symmetric matrix consists of 4 orthonormal vectors whose sum points in the t direction
  ;;The periodicity distance is the distance between two vertices of the cube with these edges
  ;;that are 2 edges apart, which is sqrt(2).
  (let ((orthonormal (make-array '(4 4) :element-type 'double-float
			   :initial-contents '(( 0.5 -0.5 -0.5  0.5)
					       (-0.5  0.5 -0.5  0.5)
					       (-0.5 -0.5  0.5  0.5)
					       ( 0.5  0.5  0.5  0.5)))))
    (let ((edges (smult (/ *size* (sqrt 2.0)) orthonormal)) ;Rescale so periodicity distance is SIZE
	  ;;Generate perpendicular vectors.  edge_i dot perpendicular_j = delta_ij
	  (perpendiculars (smult (/ (sqrt 2.0) *size*) orthonormal)))
      ;;Now squish cubes by factor sqrt{3} in the time direction.
      ;;Edges become more spacelike, perpendiculars more timelike
      ;;Edges are the columns of the edges matrix.  Perpendiculars are the rows of the perpendiculars matrix.
      ;;The two matrices are inverse.
      ;;We decrease the last row of the edges matrix and increase the last column of the perpendiculars
      (dotimes (i 4)
	(setf (aref edges 3 i) (/ (aref edges 3 i) (sqrt 3.0))
	      (aref perpendiculars i 3) (* (aref perpendiculars i 3) (sqrt 3.0))))
      (setq *itox* edges *xtoi* perpendiculars)
      (dotimes (j 4)
	(setf (aref *successor-offsets* j)
	      (list-4vector (loop for i below 4 collect (aref edges i j)))))))
  ;;Now that conversion matrices are set up, we can figure out where we are
  (setq *global-ijkl* ijkl
	*global-location* (4vector+ (itox ijkl) *ijkl-origin*)
	*job-number* (number-of-job ijkl)
	*job-end* (job-end-t))
  (when job-number (unless (= *job-number* job-number)
		     (error "Something went wrong.  We ended up with job number ~D instead of supplied ~D"
			    *job-number* job-number)))
  (setf (aref *successor-offsets* dump-destination) (4vector- *global-location*)) ;Offset dumps to global coords
  )

(declaim (inline ijkl-mine-p))
(defun ijkl-mine-p (ijkl &optional definitely)
  (declare (optimize speed)
	   (4vector ijkl))
  (loop with fudge double-float = (if definitely (- fudge-ijkl) ;invert usual fudging
				    fudge-ijkl)
	for index below 4 always (< (- fudge) (aref ijkl index) (+ 1.0 fudge))
	))

;;Tell if point is in our volume.  If DEFINITELY is set, uncertain points give NIL, otherwise T
(defun point-mine (point &optional definitely)
  (declare (optimize speed))
  (if *total-size*
      (let ((ijkl (xtoi point)))
	(prog1 (ijkl-mine-p ijkl definitely)
	  (deallocate 4vectors ijkl)))
    t					;infinite volume
    ))

;;Tell if point is inside us and our future, or the same for the successor
;;In the successor case, we attempt to deal with the fact that periodicity gives several ways to reach
;;a given point.
(defun point-inside-future-p (position &optional successor)
  (let ((ijkl (xtoi position)))
    (when successor
      (decf (aref ijkl successor) 1.0)	  ;Convert into his coordinates
      (setq ijkl (standardize-ijkl ijkl))) ;Get simplest path to this point
    (loop for index below 4 always
	  (< (- fudge-ijkl) (aref ijkl index)))))

;;Tell if point is inside our volume or in our past
(defun point-inside-past-p (position)
  (let ((ijkl (xtoi position)))
    (loop for index below 4 always
	  (< (aref ijkl index) (+ 1.0 fudge-ijkl)))))

;;Amount to add to output coordinates to convert them used by successor in given direction.
;;We need to subtract the edge, which is the column of the *itox* matrix.
(defun destination-offset (destination)
  (aref *successor-offsets* destination))

;;Local time of beginning corner of job volume
(defun job-start-t ()
  0.0)

;;Local time of ending corner.  Cube edge length is size over sqrt{2}.  Main diagonal on on squished cube is
;;2 x edge.  Squishing factor sqrt{3}.
(defun job-end-t (&optional (size *size*))
  (* size (sqrt (/ 2.0 3.0))))

;;Total range of time coordinate of this job.
(defun job-duration (&optional (size *size*))
  (job-end-t size))

;;Total 4-volume of this job
(defun job-volume ()
  (/ (expt *size* 4) 4 (sqrt 3.0)))

;;Total 3-volume of simulation.  Ignore squishing.  Then 4 cubes of edge length e, starting at different
;;times advance the simulation in time by the main diagonal of 1 cube, which is 2e.  The 3 volume is thus
;;4 e^4/(2 e) = 2 e^3.  The periodicity distance *total-size* is sqrt{2} e, because it is generated by 2
;;cube edges.
(defun total-3volume ()
  (/ (expt *total-size* 3) (sqrt 2.0)))

;;Volume of space between starting and current time
(defun job-volume-so-far (time)
  (let* ((ijkl-max (/ (* 4 time) (job-duration))) ;i+j+k+l < this
	 (ijkl-volume (/ (expt ijkl-max 4) 24)))  ;Volume in ijkl space
    (/ ijkl-volume (job-volume)))		  ;Divide by volume of unit ijkl-cube
  )

;;Maximum value of each spatial coordinate.
;;Minimum value is negative this
;;This is the distance from the center to a point of the octahedron.
(defun job-coordinate-maximum ()
  (/ *size* (sqrt 2.0)))

(defun job-coordinate-minimum ()
  (- (job-coordinate-maximum)))

;;Total range in each coordinate
;;This is used for cells and initial conditions.
;;It would be more efficient to reorient the axes to take the minimum bounding rectangular volume.
(defun job-coordinate-range ()
  (* *size* (sqrt 2.0)))

;;Maximum of each spatial coordinate in this job at the given time
(defun job-slice-coordinate-maximum (time)
  (setq time (min time (- (job-end-t) time))) ;Minimum distance from starting or ending corner
  (* (/ time (/ (job-end-t) 2))		      ;Interpolate between 0 and max at time goes from 0 to halfway point
     (job-coordinate-maximum)))

(defun job-slice-coordinate-minimum (time)
  (- (job-slice-coordinate-maximum time)))

(defun job-slice-coordinate-range (time)
  (* 2 (job-slice-coordinate-maximum time)))

;;;Job numbers

(defvar *layer-job-location*)		;Array: (jobnumber, index) -> pj, pk, or pl within a layer
(defvar *layer-location-job*)		;Array: (pj pk pl) -> job number within layer

(defun setup-split-factor (n)
  (setq *split-factor* n)
  (setq *layer-job-location* (make-array (list (expt n 3) 3))
	*layer-location-job* (make-array (list n n n)))
  (let ((job 0))
    ;;Go through jobs in order of j+k+l
    (dotimes (jkl (- (* n 3) 2))		;j+k+l = 0..3n-3
      (loop for jk from (max 0 (- jkl (1- n)))	;Minimum value for j+k.  Can't have l larger than n-1.
	    to (min jkl (* (1- n) 2))	;Maximum value for j+k.  Each can be as large as n-1.
	    as l = (- jkl jk)
	    do (loop for j from (max 0 (- jk (1- n))) ;Minimum value for j.  Can't have k larger than n-1.
		     to (min jk (1- n))	;Maximum value for j.
		     as k = (- jk j)
		     do (setf (aref *layer-job-location* job 0) j
			      (aref *layer-job-location* job 1) k
			      (aref *layer-job-location* job 2) l
			      (aref *layer-location-job* j k l) job)
		     do (incf job))))))

;;A job code is 4 values: s, p1, p2, p3
(defun job-code (job-number)
  (multiple-value-bind (layer sub) (truncate job-number (expt *split-factor* 3)) ;Get layer and location within layer
    (values layer (aref *layer-job-location* sub 0) (aref *layer-job-location* sub 1) (aref *layer-job-location* sub 2))))

;;Convert code to job number
(defun code-job-number (s p1 p2 p3)
  (+ (* s (expt *split-factor* 3)) (aref *layer-location-job* p1 p2 p3)))

;;Given job number, return canonical ijkl , which is standard version of
;;s(1,0,0,0) + p1(1,-1,0,0) + p2(1,0,-1,0) + p3(1,0,0,-1)
(defun job-coordinates (job-number)
  (multiple-value-bind (time p1 p2 p3) (job-code job-number)
    (standardize-ijkl (make-4vector (float (+ time p1 p2 p3) 0.0) (float (- p1) 0.0)
				    (float (- p2) 0.0) (float (- p3) 0.0)))))

;;Given integer ijkl coordinates of starting point, return job number
;;This is invariant under periodicity, meaning adding the split-factor N to one component and subtracting
;;it from another.  We write ijkl = s(1,0,0,0) + p1(1,-1,0,0) + p2(1,0,-1,0) + p3(1,0,0,-1)
(defun number-of-job (ijkl)
  (declare (optimize speed)		;Used in loops-graph
	   (4vector ijkl))
  (let* ((sum (fixnum-round (loop for value double-float across ijkl sum value double-float)))
	 (p3 (mod (- (fixnum-round (aref ijkl 3))) *split-factor*))
	 (p2 (mod (- (fixnum-round (aref ijkl 2))) *split-factor*))
	 (p1 (mod (- (fixnum-round (aref ijkl 1))) *split-factor*)))
    (declare (type (integer 0 10000) sum p1 p2 p3)
	     (optimize (safety 0)))	;Don't check type assertion
    (code-job-number sum p1 p2 p3)))

;;Give job number of successor.  This does not require initialization.
(defun successor-job-number (successor-index &optional (job-number *job-number*))
  (multiple-value-bind (s p1 p2 p3) (job-code job-number)
    (incf s)				;This moves in i direction
    (case successor-index
      (1 (when (minusp (decf p1)) (setq p1 (1- *split-factor*)))) ;Switch to j direction.  See notes.
      (2 (when (minusp (decf p2)) (setq p2 (1- *split-factor*))))			;Switch to k direction
      (3 (when (minusp (decf p3)) (setq p3 (1- *split-factor*)))))			;Switch to l direction
    (code-job-number s p1 p2 p3)))

;;Inverse of above
(defun predecessor-job-number (predecessor-index &optional (job-number *job-number*))
  (multiple-value-bind (s p1 p2 p3) (job-code job-number)
    (decf s)				;This moves in -i direction
    (case predecessor-index
      (1 (when (= (incf p1) *split-factor*) (setq p1 0))) ;Switch to -j direction.  See notes.
      (2 (when (= (incf p2) *split-factor*) (setq p2 0)))	;Switch to k direction
      (3 (when (= (incf p3) *split-factor*) (setq p3 0)))) ;Switch to l direction
    (code-job-number s p1 p2 p3)))

;;Return list of predecessor job numbers.  If split-factor is 1, we give the same number 4 times.
(defun predecessor-job-numbers (job)
  (unless (< job (expt *split-factor* 3))	;Early jobs have no predecessors
    (loop for predecessor below 4
	  collect (predecessor-job-number predecessor job))))

;;Return list of successor job numbers
(defun successor-job-numbers (job)
  (loop for successor below 4
	collect (successor-job-number successor job)))

;;Get all predecessors of job, including the job itself
;;Split factor must be set up
(defun all-ancestors (job)
  (sort
   (reduce #'union
	   (cons (list job)
		 (loop for predecessor in (predecessor-job-numbers job) ;Examine all predecessors
		       collect (all-ancestors predecessor))))
   #'<))


;;Starting time for given job number
;;The periodicity distance is SIZE.  The total length of the cube side before squishing is size/sqrt{2}.
;;The component in the time direction is size/(2sqrt{2}).  After squishing this becomes size/(2sqrt{6}).
;;Used to read dumps
(defun global-job-start (total-size time-offset job-number)
  (+ (* (reduce #'+ (job-coordinates job-number))
	(/ total-size *split-factor* 2 (sqrt 6.0)))
     time-offset))

(defun global-job-end (total-size time-offset job-number)
  (+ (* (+ (reduce #'+ (job-coordinates job-number)) 4) ;4 more edges to end
	(/ total-size *split-factor* 2 (sqrt 6.0)))
     time-offset))

(declaim (inline global-time local-time))
(defun global-time (local-time)
  (+ local-time (4vector-t *global-location*)))

(defun local-time (global-time)
  (- global-time (4vector-t *global-location*)))

;;End time of predecessor's volume.  It is the 3/4 point of ours.
(defun predecessor-global-end-t ()
  (global-time (/ (* (job-end-t) 3) 4)))

;;The total amount of time since the initial conditions completely covered by a given number of layers of a given
;;size.  Each layer covers one fourth the duration associated with that size, but after the first
;;4 layers we have only covered duration/8, because the first layer ends there.
(defun layers-duration (size layers)
  (* (job-duration size) (/ (- layers 3.5) 4.0)))

;;Give the number of layers of job-size SIZE needed to cover a given duration
(defun duration-layers (size duration)
  (values (ceiling (+ 3.5 (* (/ 4.0 (job-duration size)) duration)))))

(defun duration-jobs (total-size split-factor duration)
  (* (expt split-factor 3) (duration-layers (/ total-size split-factor) duration)))

(defun jobs-duration (total-size split-factor jobs)
  (layers-duration (/ total-size split-factor) (/ jobs (* (expt split-factor 3)))))

(defun job-layer (job-number &optional (split-factor *split-factor*))
  (floor job-number (expt split-factor 3)))

(defun successor-name (successor)
  (char "IJKL" successor))

;;A list of periodicity 3vectors
(defun periodicity-vectors ()
  (loop with offset = (/ *total-size* (sqrt 2.0)) ;Periodicity vectors have +/- this in exactly two slots
	for index below 3
	nconc (loop for index-2 from (1+ index) below 3
		    nconc (loop for sign in '(1 -1)
				nconc (loop for sign-2 in '(1 -1)
					    for result = (make-zero-3vector)
					    do (setf (aref result index) (* sign offset)
						     (aref result index-2) (* sign-2 offset))
					    collect result)))))

;; three linearly independent periodicity vectors with positive components
(defun periodicity-basis ()
  (loop with offset = (/ *total-size* (sqrt 2.0)) ;Periodicity vectors have +/- this in exactly two slots
	for index below 3
	nconc (loop for index-2 from (1+ index) below 3
		    for result = (make-zero-3vector)
		    do (setf (aref result index) offset
			     (aref result index-2) offset)
		    collect result)))


;; dual vectors to periodicity-basis
(defun periodicity-dual-basis ()
  (loop for index below 3
	with basis = (periodicity-basis)
	for vector = (nth index basis)
	for others = (remove vector basis)
	for orthogonal = (3vector-cross-product (nth 0 others) (nth 1 others))
	for ortho-normal = (3vector-scale orthogonal (/ 1.0 (3vector-dot orthogonal vector)))
	collect ortho-normal))


;;Return a list of lines as 2-element lists of 3vectors giving the outline of the rhombic dodecahedron
;;The order-3 vertices are just the space part of having one of i,j,k,l=+/-1, and the others zero.
;;The order-4 vertices (the points of the octahedron) are found by setting two of these coords to the same value
;;You must set up geometry first
(defun unit-cell-outline ()
  (mapcar #'itox			;Actually returns 4vectors, but that should be OK
	  ;;Get lists in ijkl coordinates
	  (loop for sign in '(1.0 -1.0)
		nconc (loop for index below 4
			    for corner3 = (make-zero-4vector)
			    do (setf (aref corner3 index) sign)
			    nconc (loop for index-2 below 4
					unless (= index index-2)
					collect (list corner3 (let ((corner4 (copy-4vector corner3)))
								(setf (aref corner4 index-2) sign)
								corner4)
						      ))))))

;;Similar but give faces to gnuplot as 2x2 grids
(defun unit-cell-faces ()
  (mapcar #'(lambda (point)
	      (and point (4to3vector (itox point))))
	  ;;Get lists in ijkl coordinates
	  (loop for index below 4
		for corner1 = (make-zero-4vector)
		do (setf (aref corner1 index) 1.0) ;Each positive face
		nconc (loop for index-2 below 4
			    unless (= index index-2) ;Each negative face adjoins
			    nconc (loop for index-3 below 4
					for index-4 = (- 6 index index-2 index-3)
					unless (or (= index-3 index) (= index-3 index-2)) ;other two indicies
					nconc (list corner1 ;First 3-vertex
						    (let ((corner (copy-4vector corner1)))
						      (setf (aref corner index-3) 1.0)
						      corner) ;First 4-vertex
						    nil ;Break line
						    (let ((corner (copy-4vector corner1)))
						      (setf (aref corner index-4) 1.0) ;Second 4-vertex
						      corner)
						    (let ((corner (make-zero-4vector)))
						      (setf (aref corner index-2) -1.0) ;Second 3-vertex
						      corner)
						    nil nil) ;New surface
					)))))

;;This generates the rhombic dodecahedron figure.
(defun plot-cell-faces (&rest gnuplot-keys)
  (let ((faces (unit-cell-faces)))
    (apply #'gnuplot 1 (length faces)
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (and (numberp point)
		    (values-list (coerce (nth point faces) 'list))))
	   :3d t
	   :styles :lines
	   gnuplot-keys)))

;;Convert ijkl into standard form by adding and subtracting the split factor N, which represents going to a different
;;image of the point in the covering space.  Suppose i is the smallest component and j the largest.  Then we should
;;consider adding N to i and subtracting N from j, i.e. moving the point by N(a_i - a_j).
;;This helps as long as j-i>N, so that |j-i| is decreased by this operation.
;;So define the excess in a component to be the multiple of N that it is larger than the smallest component.
;;We should redistribute the total excess in a uniform way.
;;Often i,j,k,l are integers.  In that case we standardize on distributing the excess into the first components first.
;;See standardize-vv for a probably superior algorithm.
;;This always returns a new 4vector, even if the original one is standard
(defun standardize-ijkl (ijkl &optional (split-factor *split-factor*))
  (declare;; (optimize speed)
	   (type 4vector ijkl)
	   (type (unsigned-byte 20) split-factor)) ;Let compiler know that this will not cause overflow
  (let ((min (reduce #'min ijkl))	;Smallest component
	(total-excess 0)		;Count total number of multiples of N above smallest component
	(result (copy-4vector ijkl)))	;Don't clobber input vector
    (declare (double-float min)
	     (type (unsigned-byte 20) total-excess)
	     (type 4vector result))
    (dotimes (index 4)
      (let ((excess (fixnum-floor (/ (+ (- (aref result index) min) fudge-ijkl) ;Round up exact excesses
				     split-factor))))
	(declare (type (unsigned-byte 20) excess))
	(decf (aref result index) (* excess split-factor)) ;For now, remove excesses
	(incf total-excess excess)))			   ;Keep track of total
    ;;It would now be sufficient to go through the total excess and add it back in, choosing the smallest coordinate
    ;;to advance in each case.  But in case the excess is large, we can distribute the multiple of 4
    ;;right away first.
    (multiple-value-bind (each extra) (truncate total-excess 4) ;Split up excess across coordinates
      (dotimes (index 4)
	(incf (aref result index) (* each split-factor))) ;Add in uniform part
      ;;Find the smallest ijkl coordinates to give the extra excess.
      (dotimes (i extra)
	(loop with least = 0		;Start with first slot
	      for index from 1 below 4		;Try others
	      ;;Look for the smallest slot to increment, but if there is a near tie, use the earlier slot
	      when (< (aref result index) (- (aref result least) fudge-ijkl)) ;This one significantly better
	      do (setq least index)
	      finally (incf (aref result least) split-factor))))
    result))

;;Map global position into the fundamental rhombic dodecahedron around the origin or a given point
;;Works for local positions too
(defun standardize-position (position &optional (center zero-4vector))
  (cond (*total-size*
	 ;;Since 3vectors and 4vectors are indistinguishable, this doesn't work
	 ;;(unless (typep position '4vector) ;If not a 4vector already, make one
	 ;;(setq position (3to4vector position 0.0)))
	 ;;Get offset from center, convert to ijkl, standardize that, convert back, and add center
	 (let ((result (4vector+ (itox (standardize-ijkl (xtoi (4vector- position center)))) center)))
	   (if (4vector= result position fudge-coordinates) ;Didn't change?
	       position					    ;Return original 4vector
	     result)))
	(t				;Infinite volume
	 position)))

;;Convert a local position to a standardized global position
(defun globalize-position (position)
;;Since 3vectors and 4vectors are indistinguishable, this doesn't work
;;  (unless (typep position '4vector)
;;    (setq position (3to4vector position 0.0)))
  (standardize-position (4vector+ position *global-location*)))

;;Convert a global position into a standardized local position
(defun localize-position (position)
;;Since 3vectors and 4vectors are indistinguishable, this doesn't work
;;  (unless (typep position '4vector)
;;    (setq position (3to4vector position 0.0)))
  (standardize-position (4vector- position *global-location*)))

;;Split a line into points in the unit cell.
(defun split-line (start end)
  (let ((offset (standardize-position (3vector- end start))))
    (unless (< (3vector-length offset) (+ diamond-span fudge-coordinates))
      (error "Line to be split is too long"))
    (let* ((splits (loop for p in (periodicity-vectors)
			 ;;Find where line intersects face perp to p
			 for a = (/ (- (3vector-dot p p) 2) (3vector-dot offset p))
			 when (<= 0.0 a 1.0)
			 collect (list a p)))
	   (result (loop for (a p) in splits
			 for point = (3vector+ start (3vector-scale offset a))
			 collect point
			 collect (3vector+ point p))))
      (nconc (list start) result (list end)))))

;;Tell if point is in box.  Dimension is 3 or 4.
;;If DEFINITELY is set, we return T only if the point is further than the fudge
;;factor from the boundary.  If not, we return T even if it is a little outside
(defun point-in-box (location minima maxima dimensions &optional definitely)
  (declare (optimize speed)
	   (type 4vector location minima maxima)
	   (type fixnum dimensions))
  (loop with fudge double-float = (if definitely (- fudge-coordinates) ;invert usual fudging
				    fudge-coordinates)
	for i below dimensions
	as this double-float = (aref location i)
	as min double-float = (aref minima i)
	as max double-float = (aref maxima i)
	always (and (< (- min fudge) this)
		    (< this (+ max fudge)))))


;;Gives vector given a spherical coordinates
(defun spherical-coordinates (theta phi &optional (r 1.0))
  (make-3vector (* r (sin theta) (cos phi)) (* r (sin theta) (sin phi)) (* r (cos theta))))

(declaim (inline spherical-angle))

;;Angle between 2 unit vectors x and y in range 0..pi
;;The simple formula would be acos(x . y), but numerical stability for small angles is better by using this:
;;(x-y)^2 = 2(1-cos theta) = 4 sin^2 theta/2, so |x-y| = 2 sin theta/2, and theta = 2 asin(|x-y|/2)
(defun spherical-angle (x-hat y-hat)
  (declare (optimize speed)
	   (muffle-conditions compiler-note)) ;Avoid warning about consing return value
  (* 2 (asin (the (double-float 0.0 1.0)
					;avoid >1 from numerical error when almost opposite
	       (min (/ (3vector-length (3vector- x-hat y-hat)) 2) 1.0)
	     ))))

;;Azimuthal angle in spherical polar coordinates defined by unit vectors Z X Y.  The angle is the one from
;;X toward Y of the direction projected into the X-Y plane.  If you don't supply y, we'll figure it out from Z and X
(defun azimuthal-angle (direction z x &optional (y (3vector-cross-product z x)))
  (atan (3vector-dot direction y) (3vector-dot direction x)))
				  
;;Given two unit vectors, return a vector lying on the great circle between the other two the given fraction away
(defun spherical-interpolation (a b fraction)
  (if (3vector= a b 1e-15) a
    (let* ((angle (spherical-angle a b))
	   (cos (cos angle))
	   (left (cos (* angle fraction)))	;Cos of angle from a
	   (right (cos (* angle (1- fraction))))) ;Cos of angle from b
      ;;Is this right?
      (3vector-normalize (3vector+ (3vector-scale a (- left (* right cos)))
				   (3vector-scale b (- right (* left cos))))))))

;;Centroid of a triangle.  We use the planar formula and then renormalize.
(defun triangle-center (triangle)
  (destructuring-bind (a b c) triangle
    (3vector-normalize (3vector+ a b c))))

;;The same but with the arguments spread
(defun spherical-triangle-center (a b c)
  (3vector-normalize (3vector+ a b c)))

(defun spherical-midpoint (a b)
  (3vector-normalize (3vector+ a b)))

;;A spherical triangle is represented by a list of 3 unit 3-vectors, which we take to go around the triangle 
;;counterclockwise, as seen from outside the sphere.


;;Return an array of triangles giving the faces of an icosahedron
(defun icosahedron ()
  (loop with z = (/ 1 (sqrt 5.0))	;z coordinate of 5 vertices
	with r = (* 2 z)		;their r coordinate
	with top = (make-3vector 0.0 0.0 1.0)
	with bottom = (make-3vector 0.0 0.0 -1.0)
	with faces = (make-array 20 :fill-pointer 0)
	for i below 10
	for phi = (/ (* i 2 pi) 10)	;ten points evenly spaced around circle
	if (evenp i)
	collect (make-3vector (* r (cos phi)) (* r (sin phi)) z) into north ;north 5 vertices
	else collect (make-3vector (* r (cos phi)) (* r (sin phi)) (- z)) into south	;south 5 vertices
	finally
	(loop for this below 5
	      for next = (mod (1+ this) 5)
	      do (vector-push (list top (nth this north) (nth next north)) faces)
	      do (vector-push (list bottom (nth next south) (nth this south)) faces)
	      do (vector-push (list (nth next north) (nth this north) (nth this south)) faces)
	      do (vector-push (list (nth next north) (nth this south) (nth next south)) faces)
	      )
	(return faces)))
	      
(defvar *triangulate-sphere-rotation* nil)

;;Return set of triangular faces that cover the sphere.  We start with an icosahedron and then
;;split each triangle into 4 LEVELS times by taking the center of each edge, as recommended by
;;Atkinson, "Numerical integration on the sphere", J. Austral. Math. Soc. (Series B) 23, 332 (1982)
;;ROTATION is a rotation matrix to apply to the positions for testing
(defun triangulate-sphere (levels &key (rotation *triangulate-sphere-rotation*))
  (let ((triangulation (split-triangulation (icosahedron) levels)))
    (cond (rotation
	   (when (= (array-dimension rotation 1) 3) (setq rotation (3to4matrix rotation)))
	   (map 'vector
		#'(lambda (triangle) (mapcar #'(lambda (position) (dotmv rotation position)) triangle))
		triangulation))
	  (t triangulation))))
      
;;Rotation matrix for an active intrinsic rotation around z by alpha, then x by beta, then z' by gamma
(defun euler-rotation-matrix (alpha beta gamma)
  (let ((c1 (cos alpha))
	(c2 (cos beta))
	(c3 (cos gamma))
	(s1 (sin alpha))
	(s2 (sin beta))
	(s3 (sin gamma)))
    (make-array '(3 3) :element-type 'double-float
		:initial-contents (list (list (- (* c1 c3) (* c2 s1 s3)) (- (+ (* c1 s3) (* c2 c3 s1))) (* s1 s2))
					(list (+ (* c3 s1) (* c1 c2 s3)) (- (* c1 c2 c3) (* s1 s3)) (- (* c1 s2)))
					(list (* s2 s3) (* c3 s2) c2)))))

;;Convert 3x3 to 4x4 matrix for DOTMV, etc.
(defun 3to4matrix (m)
  (let ((result (make-array '(4 4) :element-type 'double-float :initial-element 0.0)))
    (dotimes (i 3)
      (dotimes (j 3)
	(setf (aref result i j) (aref m i j))))
  result))


;;Triangulate sphere, separating those triangles that have any point near a cusp.
;;Returns an array of (triangle . cusp-info-list).  The cusp-directions are the EQ values from the list.
;; Array of lists of triangles near each cusp, corresponding to the cusp directions given.
;;A triangle which is not in the first list appears in each of the other lists when the distance to the corresponding
;;cusp is less than threshold2.
(defun triangulate-sphere-separate-cusps (levels cusp-data threshold1 threshold2)
  (assert (<= threshold1 threshold2))
  (map 'vector
       #'(lambda (triangle)
	   (cons triangle
		 (and (loop for info in cusp-data
			    thereis (< (spherical-triangle-point-distance triangle (cusp-info-direction info))
				       threshold1))	     ;Close to any cusp?
		      (loop for info in cusp-data ;Close to this cusp with weaker threshold?
			    when (< (spherical-triangle-point-distance triangle (cusp-info-direction info))
				       threshold2)
			    collect info))))
       (triangulate-sphere levels)))

;;Tell whether the point is inside the spherical triangle.  Each side of the triangle is part of a great circle,
;;and each great circle divides the sphere into two parts.  Since triangle points go around counterclockwise,
;;If the point is inside the triangle it should be to the left of each great circle, defined by going forward
;;through the list of points.
(defun point-inside-spherical-triangle (point triangle)
  (loop for (a b) on triangle
	unless b do (setq b (first triangle)) ;Loop around
	always (plusp (3vector-dot (3vector-cross-product a b) point))))

;;Minimum distance of given spherical triangle from given point on sphere(
(defun spherical-triangle-point-distance (triangle point)
  (if (point-inside-spherical-triangle point triangle) 0.0
    (min (arc-point-distance (first triangle) (second triangle) point) ;Outside: use closest edge
	 (arc-point-distance (second triangle) (third triangle) point)
	 (arc-point-distance (third triangle) (first triangle) point))))

;;Find minimum spherical angle between POINT and any point on the (shorter) great circle segment from a to b.
;;These points can be written
;;x(phi) = (a sin(theta-phi) + b sin phi)/sin theta, where theta is the angle between a and b and phi = 0...theta
;;The closest approach to a point c is the place where c.x is largest, which is either a, b, or a point where
;;(d/dphi)(c.x) = 0, so (c.b cos phi - a cos(theta-phi)) = 0 or tan phi = (c.b - c.a cos theta)/(c.a sin theta)
;;There are two such points, with phi differing by pi.  The closest approach is the one where x.c > 0.
;;If c.a > 0, then this is the one where x is within a quarter circle of alpha, i.e. cos phi > 0,
;;whereas if c.a < 0, it is the one where x is in the other half circle where cos phi < 0.
(defun arc-point-distance (a b point)
  (let* ((theta (spherical-angle a b))
	 (pa (3vector-dot a point))
	 (pb (3vector-dot b point))
	 (phi (if (zerop pa)					     ;Avoid errors if perpendicular
		  (* (/ pi 2))					     ;extrema in direction of b/-b then
		(atan (/ (- (/ pb pa) (cos theta)) (sin theta)))))) ;Angle of extremum between -pi/2 and pi/2
    (when (minusp pa)						 ;If c.a<0, minimum is other extremum
      (setq phi (+ phi pi)))
;;    (format t "theta = ~S, pa = ~S, pb = ~S, phi = ~S~%" theta pa pb phi)
    (if (< 0.0 phi theta)		;Point in segment?  If so, it is closest.
	(spherical-angle point (3vector+ (3vector-scale a (/ (sin (- theta phi)) (sin theta)))
					 (3vector-scale b (/ (sin phi) (sin theta)))))
      (min (spherical-angle point a) (spherical-angle point b))))) ;Otherwise closer end

;;Split at least min-levels, and then only if split-p-function returns true
(defun split-triangulation (faces levels &optional split-p-function min-levels)
  (dotimes (i levels)
    (let ((result (make-array (* (length faces) 4) :fill-pointer 0)))
      (loop for (a b c) across faces
	    if (or (null split-p-function) (< i min-levels)
		   (funcall split-p-function a b c)) ;Split it?
	    do (let ((ab (3vector-normalize (3vector+ a b))) ;Center of each edge
		     (ac (3vector-normalize (3vector+ a c)))
		     (bc (3vector-normalize (3vector+ b c))))
		 (vector-push (list a ab ac) result)
		 (vector-push (list b bc ab) result)
		 (vector-push (list c ac bc) result)
		 (vector-push (list ab bc ac) result))
	    else do (vector-push (list a b c) result) ;Keep original unsplit
	    )
      (setq faces result)))		;Replace input with next split
  faces)

(defparameter spherical-triangle-degeneracy-threshold 1e-13)

;;See if spherical triangle is degenerate, in the sense that two sides lie nearly on top of each other.
(defun spherical-triangle-degenerate-p (a b c)
  (let* ((ab (spherical-angle a b))	;Angles of sides, with care for small angles
	 (ac (spherical-angle a c))
	 (bc (spherical-angle b c)))
    (< (min (- (+ ab ac) bc) (- (+ bc ab) ac) (- (+ ac bc) ab))
       (* (max ab ac bc) spherical-triangle-degeneracy-threshold)))) ;Triangle inequality only barely satisfied

#|
;;Find the unsigned interior angle at a of a spherical triangle using spherical law of cosines
;;Old version
(defun spherical-triangle-angle (a b c)
  (let* ((ab (spherical-angle a b))	;Angles of sides, with care for small angles
	 (ac (spherical-angle a c))
	 (bc (spherical-angle b c))
	 (min (min ab ac bc)))
    (when (< min 1e-14)
      (error "Points in spherical-triangle-angle are distant only by ~S" min))
    (unless (> (min (- (+ ab ac) bc) (- (+ bc ab) ac) (- (+ ac bc) ab)) -1e-14) ;allow tiny violation
      (error "Spherical triangle does not satisfy triangle inequality with side lengths ~S, ~S, ~S"
	     ab ac bc))
    (let ((cos (/ (- (cos bc) (* (cos ab) (cos ac))) (sin ab) (sin ac))))
      (if (> cos 1)			;This can happen by numerical error when the triangle is almost degenerate
	  0.0
	(acos cos)))))
|#
	   
;;Find the unsigned interior angle at a of a spherical triangle.
(defun spherical-triangle-angle (a b c)
  (spherical-angle (3vector-normalize (3vector-cross-product a b))
		   (3vector-normalize (3vector-cross-product a c))))

;;Signed angle at a.  Positive if looking down on the sphere at a, we go counterclockwise from c to b.
;;Thus a normal spherical triangle with the three corners in counterclockwise order has negative angles
(defun signed-spherical-triangle-angle (a b c)
  (* (spherical-triangle-angle a b c) (signum (3vector-dot (3vector-cross-product a c) b))))

;;Find the area of the spherical triangle given by the 3 vertices
(defun spherical-triangle-area (a b c)
  (+ (spherical-triangle-angle a b c)
     (spherical-triangle-angle b c a)
     (spherical-triangle-angle c a b)
     (- pi)))

;;Find the point, if any, where the great circle from a1 to a2 intersects that from b1 to b2.
;;We do it by defining the plane spanned by each pair of vectors, finding the line on which the two
;;planes intersect, and normalizing the vector that generates this line to 1.
(defun spherical-intersection (a1 a2 b1 b2)
  (when (fast-no-spherical-intersection a1 a2 b1 b2)
    (return-from spherical-intersection nil))
  (let ((a (make-array '(4 4) :element-type 'double-float :initial-element 0.0)))
    ;;There are 3 equations (for x, y and z) in 4 unknown and coefficients given by
    ;;c1 a1 + c2 a2 - d1 b1 - d2 b2 = 0
    ;;The columns of A are c1, c2, d1, d2.  The rows are x,y,z and an extra row to make it square
    (dotimes (coordinate 3)
      (setf (aref a coordinate 0) (3vector-component a1 coordinate)
	    (aref a coordinate 1) (3vector-component a2 coordinate)
	    (aref a coordinate 2) (3vector-component b1 coordinate)
	    (aref a coordinate 3) (3vector-component b2 coordinate)))
    (multiple-value-bind (u w v) (svdcmp a)
      (declare (ignore u))
      (let ((zero-index nil))
	(loop for index below 4		;Find zero element
	      when (< (abs (aref w index)) 1e-15)
	      do (if zero-index (error "Two zero indices finding sphere intersection")
		   (setq zero-index index)))
	(unless zero-index (error "Couldn't find 0 in ~S" w))
	;;Now the given column of v is the null space
	(let ((x (3vector-normalize
		  (3vector+ (3vector-scale a1 (aref v 0 zero-index)) ;x = c1 a1 + c2 a2 = d1 b1 + d2 b2
			    (3vector-scale a2 (aref v 1 zero-index)))))) ;so it lies in both planes
	  ;;We need to know if x is on the short segment of the great circle including a1 and a2,
	  ;;rather than the long segment.  First the minimum of the angles between x and a1 and x and a2 should be less
	  ;;than the corresponding angles with -x, so we are in the right hemisphere
	  ;;The we check that angles between x and a1 and x and a2 are each less than the angle between a1 and a2
	  (let ((aa (3vector-dot a1 a2))
		(xa1 (3vector-dot x a1))
		(xa2 (3vector-dot x a2)))
	    (when (< (max xa1 xa2) (max (- xa1) (- xa2))) ;Wrong direction?
	      (setq xa1 (- xa1) xa2 (- xa2) x (3vector- x))) ;Switch
	    ;;We need to fudge here, so we don't miss cases where the two endpoints are the same or very close
	    (when (and (> (+ xa1 1e-14) aa)	;x closer to a1 than a2 is?
		       (> (+ xa2 1e-14) aa))	;x closer to a2 than a1 is?
	      (let ((bb (3vector-dot b1 b2))
		    (xb1 (3vector-dot x b1))
		    (xb2 (3vector-dot x b2)))
		(and (>= (max xb1 xb2) (max (- xb1) (- xb2))) ;Right hemisphere?
		     (> (+ xb1 1e-14) bb)		;x closer to b1 than b2 is?
		     (> (+ xb2 1e-14) bb)		;x closer to b2 than b1 is?
		     x)			;return x
		))))))))

;;Quickly find if the great circle segments are short and far away
;;If there is intersection, the distance from any point a to a point b must be no more than the
;;length of a plus the length of b
(defun fast-no-spherical-intersection (a1 a2 b1 b2)
  (< (+ (spherical-angle a1 a2) (spherical-angle b1 b2))
     (spherical-angle a1 b1)))

(defun test-spherical-intersection ()
  (flet ((random-direction () (3vector-normalize ;Random unit vector, not evenly distributed
			       (make-3vector (- (random 2.0) 1.0) (- (random 2.0) 1.0) (- (random 2.0) 1.0)))))
    (let* ((a1 (random-direction))
	   (a2 (random-direction))
	   (b1 (random-direction))
	   (b2 (random-direction))
	   (intersection (spherical-intersection a1 a2 b1 b2)))
      (format t "~&~S -- ~S~%~S -- ~S~%~S" a1 a2 b1 b2 intersection)
      (gnuplot 3 100
	       #'(lambda (plot point)
		   (and (numberp point)
			(if (= plot 2)
			    (and (= point 0)
				 intersection
				 (values-list (3vector-list intersection)))
			  (let ((from (if (zerop plot) a1 b1))
				(to (if (zerop plot) a2 b2))
				(interpolate (/ point 99)))
			    (values-list (3vector-list (3vector-normalize
							(3vector+ (3vector-scale from (- 1 interpolate))
								  (3vector-scale to interpolate)))))))))
	       :styles '(:lines :lines :points)
	       :3d t :xrange '(-1.0 . 1.0) :yrange '(-1.0 . 1.0) :zrange '(-1.0 . 1.0)
	       :reuse t
	       ))))

(defun rotate-vectors (vectors rotation)
  (map 'vector #'(lambda (v) (dotmv33 rotation v)) vectors))
