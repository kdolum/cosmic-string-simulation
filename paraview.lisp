;;; Paraview VTK output
(in-package "CL-USER")

;;;VTK output format

;;Write a single-float in big-endian format
(defun vtk-write-float (float stream)
  (setq float (single-float float))
  (setf (sap-ref-single conversion-data 0) float)
  (dotimes (byte 4)
    ;;If we are big-endian then layout in order.  Otherwise reverse.
    (write-byte (sap-ref-8 conversion-data #+big-endian byte #+little-endian (- 3 byte)) stream)))

(defun vtk-write-3vector (vector stream)
  (vtk-write-float (3vector-x vector) stream)
  (vtk-write-float (3vector-y vector) stream)
  (vtk-write-float (3vector-z vector) stream))

(defun vtk-write-32 (int stream)
  (setf (sap-ref-32 conversion-data 0) int)
  (dotimes (byte 4)
    ;;If we are big-endian then layout in order.  Otherwise reverse.
    (write-byte (sap-ref-8 conversion-data #+big-endian byte #+little-endian (- 3 byte)) stream)))

(defmacro with-open-vtk-file ((stream file &rest arguments) &body body)
  `(with-group-write-access
    (with-open-file (,stream ,file :element-type '(unsigned-byte 8) ,@arguments)
      ,@body)))

;;Write text to binary VTK file
(defun vtk-format (stream &rest format-args)
  (loop for char across (apply #'format nil format-args)
	do (write-byte (char-code char) stream)))

;;Write time-slice for paraview

;;Break plot-data into individual strings.  Returns a list of lists of tag point v point v ... point
(defun flatten-plot-data (paths range)
  (loop for path in paths
	for tag = (path-bh path)
	nconc (loop
	       ;;This returns pos v pos v ... pos with NIL instead of v when there is a break
	       with string = (coerce (get-plot-positions path (first range) (second range) t)
				     'list)
	       while string	;Nothing to do, or complete
	       collect (cons tag
			     (loop for thing = (pop string) ;Point or velocity
				   while thing ;Stop if internal NIL or end of string
				   collect thing)))))

;;Write time slice for VTK
;;See VTK User's Guide section 19.3

(defun paraview-time-slice (file &key (time (current-time)) range)
  (let ((data (flatten-plot-data (get-paths time) range)))	;List of (original-length . segment)
    (loop with starts = (make-array (length data)) ;Starting indices of strings
	  for (tag . string) in data
	  for index from 0
	  for length = (/ (1+ (length string)) 2) ;Number of points.  Velocities one less
	  unless (integerp length) do (error "Position and velocity list should have had odd length")
	  do (setf (aref starts index) total-length) ;Total length of previous strings is starting index for this one
	  sum length into total-length
	  finally
      (format t "~&~D string~:P with ~D segments in all~%" (length data) total-length)
      (with-open-vtk-file (output file :direction :output :if-exists :supersede :if-does-not-exist :create)
	(vtk-format output "# vtk DataFile Version 3.0~%")
	(vtk-format output "String dump time ~$~@[, range ~S~]~%" time range)
	(vtk-format output "BINARY~%DATASET POLYDATA~%")
	(vtk-format output "POINTS ~D float~%" total-length) ;Give total number of points
	(loop for (tag . string) in data		  ;Output points
	      do (loop for point in string by #'cddr
		       do (vtk-write-3vector point output)))
	;;Connect strings. 
	;;Parameters are number of strings and total length of connections including initial count elements
	(vtk-format output "LINES ~D ~D~%" (length data) (+ (length data) total-length))
	(loop for (tag . string) in data
	      for length = (/ (1+ (length string)) 2)
	      for index from 0
	      do (vtk-write-32 length output) ;Number of points in this line
	      do (loop for point from (aref starts index)
		       repeat length
		       do (vtk-write-32 point output) ;List of point numbers
		       ))
	(vtk-format output "POINT_DATA ~D~%" total-length)	 ;Give colors to points
	(vtk-format output "SCALARS loop-lengths float 1~%") ;Number of segments in loop
	(vtk-format output "LOOKUP_TABLE default~%")
	(loop for (tag . string) in data
	      do (loop for point on string by #'cddr ;For each point
		       do (vtk-write-float tag output)))
	;;Write velocity of segment ending at this point, because that is what Paraview will use for the tube
	;;if we turn off all interpolation
	(vtk-format output "VECTORS velocity float~%")
	(loop for (tag . string) in data
	      do (loop for v in (cons zero-4vector string) by #'cddr ;First zero, which isn't used, then velocities
		       do (vtk-write-3vector v output)))
	))))

;;Average of the function a(sigma) or b(sigma)
(defun hats-dsigmas-com (hats dsigmas)
  (loop with position = (make-zero-3vector)
	with total = (make-zero-3vector)
	for index below (length hats)
	for hat = (aref hats index)
	for dsigma = (aref dsigmas index)
	for last-dsigma = (aref dsigmas (previous-array-index-wrapping index dsigmas))
	;;Weight is average of previous length and this
	do (3vector-incf total (3vector-scale position (/ (+ dsigma last-dsigma) 2)))
	do (3vector-incf position (3vector-scale hat dsigma)) ;On to next point
	sum dsigma into length
	finally (return (3vector-scale total (/ 1.0 length)))))


;;Write a (not a') or b.  If zoom, rescale by string length.
(defun paraview-string-half (file hats dsigmas &key time zoom)
  (let ((count (length hats)))
    (with-open-vtk-file (output file :direction :output :if-exists :supersede :if-does-not-exist :create)
      (vtk-format output "# vtk DataFile Version 3.0~%")
      (vtk-format output "Dump of a/b at time ~S~%" time)
      (vtk-format output "BINARY~%DATASET POLYDATA~%")
      (vtk-format output "POINTS ~D float~%" count) ;number of segments
      (loop with scale = (if zoom (/ 1.0 (total-sigma dsigmas)))
	    with com = (hats-dsigmas-com hats dsigmas)
	    with position = (make-zero-3vector)
	    for index below count
	    for hat = (aref hats index)
	    for dsigma = (aref dsigmas index)
	    do (vtk-write-3vector (3vector-scale (3vector- position com) scale) output)
	    do (3vector-incf position (3vector-scale hat dsigma)))
      ;;Connect curve
      ;;Parameters are number of strings and total length of connections including initial count elements
      (vtk-format output "LINES 1 ~D~%" (+ count 2))
      (vtk-write-32 (1+ count) output)	;Number of points in connections.
      (loop for number below (length hats)
	    do (vtk-write-32 number output)) ;List of point numbers
      (vtk-write-32 0 output)			;Back to first point
      )))

(defun paraview-unit-cell (file)
  (let ((data (unit-cell-outline)))
    (with-open-file (output file :direction :output :if-exists :supersede :if-does-not-exist :create)
      (format output "# vtk DataFile Version 3.0~%")
      (format output "Unit cell~%")
      (format output "ASCII~%DATASET POLYDATA~%")
      (format output "POINTS 48 double~%")
      (dolist (line data)
	(dolist (point line)
	  (format output "~F ~F ~F~%" (3vector-x point) (3vector-y point) (3vector-z point))))
      (format output "LINES 24 72~%")
      (loop for (x1 x2) in data
	    for index from 0 by 2
	    do (format output "2 ~D ~D~%" index (1+ index))))))



(defun paraview-ab-movie (input-filename-prefix output-directory start end step &key zoom rotate)
  (when rotate (unless (= step 1) (error "Rotate generally doesn't work unless step = 1, because of reconnections")))
  (loop with r = (and rotate (make-identity-float 3))	;Rotation matrix
	for time from start to end by step
	for index from 0
	do (format t "~D " time) (force-output)
	do (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	       (read-hats-dsigmas (format nil "~A-~D.dat" input-filename-prefix time))
	     (mirror-images
	      (paraview-string-half (format nil "~A/~A~:[~;-zoom~]~:[~;-rotate~].~D.vtk" output-directory :a zoom
					    rotate index)
				    (if rotate
					(rotate-vectors a-hats r) ;Rotate to match original
				      a-hats)
				    a-dsigmas :zoom zoom))
	     (when rotate						;Compute rotation for later use
	       (multiple-value-bind (new-a-hats new-a-dsigmas new-b-hats) ;Get backreacted but not evolved
		   (read-hats-dsigmas (format nil "~A-~D-0.dat" input-filename-prefix (+ time step)))
		 (declare (ignore new-a-dsigmas))
		 ;;Accumulate total rotation to transform new vectors into original directions
		 (setq r (dotmm r (transpose (find-backreaction-rotation
					      a-hats a-dsigmas new-a-hats b-hats b-dsigmas new-b-hats)))))))))

;;Make paraview movie
(defun make-paraview-movie (directory start end step)
  (setq directory (merge-pathnames directory batch-root-directory))
  (loop for time from start to end by step
	for index from 0
	do (make-paraview-frame directory time index)))

(defun make-paraview-movie-submit (directory start end step &key (jobs 10))
  (load "load")
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((frames (round (- end start) step))) ;Number of frames
    (format t "~&Making ~D frame~:P in ~A~%" frames directory)
    (loop for job below jobs
	  do (do-submit directory (format nil "paraview-~D" job) (format nil "paraview-~D-" job) batch-flags
			`(make-paraview-frames ,directory ,start ,step 
					       ,@(loop for frame from job by jobs below frames collect frame))
			:load-file (format nil "~A/load" (sb-unix:posix-getcwd))))))

;;Output one frame for paraview
(defun make-paraview-frame (directory time frame)
  (format t "~&Making paraview frame ~D at time ~4$~%" frame time)
  (let ((file (format nil "~A/data.~D.vtk" directory frame)))
    (if (probe-file file)
	(warn "File ~A already exists" file)
      (let ((*allow-missing-predecessor-files* t))
	(read-dumps directory :time time)
	(paraview-time-slice file :time time)))))

;;Output many frames in batch job
(defun make-paraview-frames (directory start step &rest frames)
  (loop for frame in frames
	do (make-paraview-frame directory (+ start (* step frame)) frame)))

;;Write loop-creation positions for paraview
;;Write some paraview objects to his stream.  DATA should be a list of (points polygons &rest data)
;;In each entry, POINTS is a list of points, POLYGONS is a list of polygons using those points and
;;each element pair of DATA is some POINT_DATA associated with each point of the object
(defun write-paraview-objects (stream data title &rest data-spec)
  (vtk-format stream "# vtk DataFile Version 3.0~%~A~%" title)			      
  (vtk-format stream "BINARY~%DATASET POLYDATA~%")
  (let ((count-points 0)		;Total number of points
	(count-polygons 0)		;Total number of polygons
	(total-polygon-length 0)	;Total number of points in all polygons
	(table (make-hash-table :test #'eq)) ;Point -> point-number
	(point-number 0))
    (loop for (points polygons) in data
	  do (incf count-points (length points))
	  do (dolist (polygon polygons)
	       (incf count-polygons)
	       (incf total-polygon-length (length polygon))))
    (vtk-format stream "POINTS ~D float~%" count-points)
    (dolist (entry data)
      (dolist (point (first entry))
	(vtk-write-3vector point stream)
	(setf (gethash point table) point-number)
	(incf point-number)))
    (vtk-format stream "POLYGONS ~D ~D~%" count-polygons ;Number of polygons
		(+ count-polygons total-polygon-length)) ;Total length including length spec for each one
    (dolist (entry data)
      (dolist (polygon (second entry))
	(vtk-write-32 (length polygon) stream)
	(dolist (point polygon)
	  (vtk-write-32 (gethash point table) stream))))
    (vtk-format stream "POINT_DATA ~D~%" count-points)
    (loop for (type name) on data-spec by #'cddr
	  for slot from 2		;Slot in list of (points polygons ... data ...)
	  do (ecase type
	       (:scalar (vtk-format stream "SCALARS ~A float 1~%" name)
			(vtk-format stream "LOOKUP_TABLE default~%"))
	       (:vector (vtk-format stream "VECTORS ~A float~%" name)))
	  do (dolist (entry data)
	       (let ((value (nth slot entry)))
		  (dolist (point (first entry))	;Same value for each point
		    (declare (ignore point))
		    (funcall (case type (:scalar #'vtk-write-float) (:vector #'vtk-write-3vector))
			     value stream)))))))

(defun paraview-loop-positions (directory start end &key velocities max-worker)
  (with-open-vtk-file (output (format nil "~A/loop-positions.vtk" directory)
			  :direction :output :if-exists :supersede :if-does-not-exist :create)
    (let* ((data (read-loop-positions directory start end :velocities velocities :max-worker max-worker))
	   (points (length data)))
      (vtk-format output "# vtk DataFile Version 3.0~%")
      (vtk-format output "Loop positions, times ~$--~$~%" start end)
      (vtk-format output "BINARY~%DATASET POLYDATA~%")
      (vtk-format output "POINTS ~D float~%" points) ;Give total number of points
      (dotimes (index points)		;Point positions
	(when (zerop (mod index 1000)) (format t "~D " index) (force-output))
	(vtk-write-3vector (loop-data-position (aref data index)) output))
      (vtk-format output "VERTICES ~D ~D~%" points (* points 2))
      (dotimes (index points)
	(vtk-write-32 1 output)
	(vtk-write-32 index output))
      (vtk-format output "POINT_DATA ~D~%" points)
      (vtk-format output "SCALARS time float 1~%")
      (vtk-format output "LOOKUP_TABLE default~%")
      (dotimes (index points)		;Times
	(vtk-write-float (loop-data-time (aref data index)) output))
      (vtk-format output "SCALARS length float 1~%")
      (vtk-format output "LOOKUP_TABLE default~%")
      (dotimes (index points)		;Times
	(vtk-write-float (loop-data-length (aref data index)) output))
      (when velocities
	(vtk-format output "VECTORS velocity float~%")
	(dotimes (index points)
	  (vtk-write-3vector (loop-data-velocity (aref data index)) output))))))

;;Return an arrow for paraview.  The length of VECTOR is the overall length.  SHAFT-RADIUS is the shaft radius.
;;TIP-LENGTH is the fraction to devote to the tip, and TIP-RADIUS is the radius of the tip relative to the shaft.
(defun paraview-arrow (start vector shaft-radius &key (sides 6) (tip-length 0.3) (tip-radius 2.0))
  (format t "Length ~$, radius ~$~%" (3vector-length vector) shaft-radius)
  (let* ((tip (3vector+ start vector))	;Tip of arrow
	 (shaft-end (3vector+ start (3vector-scale vector (- 1.0 tip-length)))) ;center of end of shaft
	 (x (3vector-x vector))
	 (y (3vector-y vector))
	 (z (3vector-z vector))
	 (r (sqrt (+ (expt x 2) (expt y 2))))
	 (v1 (if (plusp r)		;Get 2 unit vectors perpendicular to each other and VECTOR
		 (3vector-normalize (make-3vector (- y) x 0.0))
	       (make-3vector 1.0 0.0 0.0))) ;If r = 0, above vector could not be normalized, so use this
	 (v2 (if (plusp r)
		 (3vector-normalize (make-3vector (* z x) (* z y) (- (expt r 2))))
	       (make-3vector 0.0 1.0 0.0)))
	 (starts (make-array sides))	;Points at base of arrow
	 (middles (make-array sides))	;points at end of shaft
	 (outers (make-array sides))) 	;points at beginning of tip.
    (values
     (nconc				;List of points
      (loop for side below sides
	    for angle = (/ (* 2 pi side) sides)
	    for radius = (3vector+ (3vector-scale v1 (* shaft-radius (sin angle)))
				   (3vector-scale v2 (* shaft-radius (cos angle))))
	    for start-point = (3vector+ start radius)
	    for middle-point = (3vector+ shaft-end radius)
	    for outer-point = (3vector+ shaft-end (3vector-scale radius tip-radius))
	    collect start-point
	    collect middle-point
	    collect outer-point
	    do (setf (aref starts side) start-point
		     (aref middles side) middle-point
		     (aref outers side) outer-point))
      (list tip))
     (cons
      (loop for side below sides collect (aref starts side)) ;Base
      (loop for side below sides
	    for next-side = (mod (1+ side) sides)
	    collect (list (aref starts side) (aref middles side) (aref middles next-side) (aref starts next-side))
	    collect (list (aref middles side) (aref outers side) (aref outers next-side) (aref middles next-side))
	    collect (list (aref outers side) tip (aref outers next-side)))))))

;;Construct arrows for each loop in a given range
(defun paraview-loop-arrows (directory start end minima maxima
				       &key velocities max-worker
				       (length-scale 3) ;Multiply velocity by this to get arrow length
				       (radius-scale 0.3) ;Multiply length by this to get radius
				       (min-loop-length 0.1) ;Smaller loops are treated as having this length
				       (max-loop-length 10) ;Larger loops are treated as having this length
				       )
  (with-open-vtk-file (output (format nil "~A/loop-arrows.vtk" directory)
			      :direction :output :if-exists :supersede :if-does-not-exist :create)
    (let* ((loops (read-loop-positions directory start end :velocities velocities :max-worker max-worker
				      :range (and minima (list minima maxima))))
	   (data nil))
      (format t "~D loops ~%" (length loops))
      (loop for loop across loops
	    for position = (loop-data-position loop)
	    for velocity = (loop-data-velocity loop)
	    for length = (loop-data-length loop)
	    for arrow-length = (* (3vector-length velocity) length-scale)
	    for arrow-radius = (* radius-scale (max min-loop-length (min max-loop-length length)))
	    ;;Offset position by 3 oscillations to go backward from the position of deletion
	    ;;to the position of emission
	    do (setq position (3vector- position (3vector-scale velocity (* length 1.5))))
	    do (format t "Time ~F at ~S: " (loop-data-time loop) position)
	    do (multiple-value-bind (points polygons)
		   (paraview-arrow position (3vector-scale (loop-data-velocity loop) arrow-length) arrow-radius)
		 (push (list points polygons (loop-data-time loop) (loop-data-length loop) (loop-data-velocity loop)) data)))
      (write-paraview-objects output data (format nil "Loop positions, times ~$--~$" start end)
			      :scalar "time" :scalar "length" :vector "velocity"))))

;;Can't use dotmv because 3vectors are really 4vectors
(defun 3vector-rotate (v x)
  (list-3vector
   (loop for i below 3
	 collect (loop for j below 3
		       sum (* (aref v i j) (3vector-component x j))))))

;;Make some still pictures (could be easily modified to make a movie) from backreaction data
(defun make-paraview-loop-frames (directory steps)
    (block scan
      (scan-backreaction-dumps
       #'(lambda (this-step a-hats a-dsigmas b-hats b-dsigmas this-r total-r)
	   (declare (ignore this-r))
	   (when (= this-step (car steps))
	     (pop steps)
	     (initialize)
	     (let ((r (transpose total-r))) ;Matrix to rotate new positions back to old frame
	       (mirror-images
		(setq a-hats (map 'vector #'(lambda (a) (3vector-rotate r a)) a-hats))))
	     (let ((com (3vector-scale (3vector+ (hats-dsigmas-com a-hats a-dsigmas) ;center of mass position
						 (hats-dsigmas-com b-hats b-dsigmas))
				       0.5)))
	       (create-worldsheet-at-time a-hats a-dsigmas b-hats b-dsigmas
					  (3to4vector (3vector-scale com -1) 0.0)
					  :install t))
	     (paraview-time-slice (format nil "~A/loop-frame-~D.vtk" directory this-step))
	     (unless steps		;Nothing left to do
	       (return-from scan nil))))
       directory)))
