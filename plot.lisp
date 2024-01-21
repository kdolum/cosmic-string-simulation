;;;Plotting
(in-package "CL-USER")

;;;Plot a few diamonds, showing edges of processor volumes, etc.

(defvar *loop-positions*)

;;Plot the world sheet of a single diamond, annotated with the edges of our simulation volume.
;;We project into the plane of the diamond by ignoring the t coordinate.  We then rotate coordinates so that the y
;;coordinate on the plot is in the direction of phat + qhat and the x coordinate in the direction of qhat - phat
(defun plot-diamond (diamond &rest gnuplot-keys
			     &key ((:dump-time *dump-time*) *dump-time*) (steps 20) &allow-other-keys)
  (let* ((p (diamond-p diamond))
	 (p5 (make-pq5 p))
	 (q (diamond-q diamond))
	 (q5 (make-pq5 q))
	 (phat (3vector-normalize p))
	 (qhat (3vector-normalize q))
	 (x (3vector-normalize (3vector- qhat phat)))
	 (y (3vector-normalize (3vector+ phat qhat)))
	 (delta (make-pq5 (4vector- (compute-diamond-end-flat diamond) (diamond-end diamond))))	;end - flat-end
	 (togo (make-togo (diamond-start diamond))) ;offsets to ijkl = 1, plus dump line
	 (tostart (4vector- (xtoi (diamond-start diamond)))) ;offsets to ijkl = 0
	 (lines				;Make 5 curves
	  (loop for end in '(nil t)	;Start and end of our volume
		nconc (loop for index below (if (and *dump-time* end) 5 4)
			    as pj = (aref p5 index)
			    as qj = (aref q5 index)
			    as deltaj = (aref delta index)
			    as offset = (aref (if end togo tostart) index)
			    collect (loop for step below steps
					  for value = (- (* 1.2 (/ step steps)) 0.1) ;From -0.1 to 1.1
					  when (if (> (abs qj) (abs pj)) ;QJ bigger so it's better to let A vary
						      (let ((b (/ (- offset (* pj value)) (- qj (* value deltaj)))))
							(and (< -0.1 b 1.1) (list value b)))
						      (let ((a (/ (- offset (* qj value)) (- pj (* value deltaj))))) ;Other way: let B vary
							(and (< -0.2 a 1.2) (list a value))))
					  collect it))))
	 (points '((0.0 0.0) (0.0 1.0) (1.0 1.0) (1.0 0.0) (0.0 0.0)))	;(a b) coordinates of corners
	 (styles '("linespoints" "lines lt 1 lw 1" "lines lt 2 lw 1" "lines lt 3 lw 1" "lines lt 4 lw 1"
		   "lines lt 1 lw 2" "lines lt 2 lw 2" "lines lt 3 lw 2" "lines lt 4 lw 2")))
	(if *dump-time* (setq styles (append styles '("lines lt -1")))) ;Style for dumping
	(flet ((plot-coordinates (&optional a b)	;Give position for plotting from (a b)
		 (and a			;Pass on NIL
		      (let ((position (diamond-position diamond :a a :b b)))
			(values (3vector-dot position x) (3vector-dot position y) (4vector-t position))))))
	  (apply
	   #'gnuplot
	   (1+ (length lines))
	   (max 5 steps)
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (if (= plot 9)
		       (format nil "T = ~A" (format-float-reasonably *dump-time*))
		     (if (zerop plot)
			 (format nil "~S" diamond)
			 (nth (1- plot) '("I=0" "J=0" "K=0" "L=0" "I=1" "J=1" "K=1" "L=1"))))
		 (apply #'plot-coordinates
			(cond ((zerop plot) ;The diamond itself
			       (nth point points))
			      (t (nth point (nth (1- plot) lines)))))))
	   :prelude (format nil "~:{set label '~C' at ~{~F, ~F, ~F~} offset character 0.3,0.3~%~}"
			    (loop for char across "srel" ;Same order as in POINTS
				  for ab in points
				  collect (list char (multiple-value-list (apply #'plot-coordinates ab)))))
	   :styles styles
	   gnuplot-keys))))


;;This plots in space-time and doesn't know about curved diamonds
(defun old-plot-diamond (diamond &rest gnuplot-keys &key ((:dump-time *dump-time*) *dump-time*) &allow-other-keys)
  (let* ((p (make-pq5 (diamond-p diamond)))
	 (q (make-pq5 (diamond-q diamond)))
	 (togo (make-togo (diamond-start diamond))) ;offsets to ijkl = 1, plus dump line
	 (tostart (4vector- (xtoi (diamond-start diamond)))) ;offsets to ijkl = 0
	 (lines				;Make 5 lines of 2 points each
	  (loop for end in '(nil t)	;Start and end of our volume
		 nconc (loop for index below (if (and *dump-time* end) 5 4)
			     as pj = (aref p index)
			     as qj = (aref q index)
			     as offset = (aref (if end togo tostart) index)
			     collect (loop for value in '(-1d0 2d0)
					   collect (if (> (abs qj) (abs pj))
						       (diamond-position
							diamond :a value :b (/ (- offset (* pj value)) qj))
						     (diamond-position
							diamond :a (/ (- offset (* qj value)) pj) :b value))))))
	 (points (list (diamond-start diamond) (diamond-right diamond) (diamond-end diamond)
		       (diamond-left diamond) (diamond-start diamond))))
    (multiple-value-bind (minima maxima)
	(diamond-bounding-box diamond)
      (let ((ranges (loop for index below 3
			  as min = (aref minima index)
			  as max = (aref maxima index)
			  as range = (- max min)
			  collect (cons (- min range fudge-coordinates)
					(+ max range fudge-coordinates))))
	    (styles '("linespoints" "lines lt 1 lw 1" "lines lt 2 lw 1" "lines lt 3 lw 1" "lines lt 4 lw 1"
		      "lines lt 1 lw 2" "lines lt 2 lw 2" "lines lt 3 lw 2" "lines lt 4 lw 2")))
	(if *dump-time* (setq styles (append styles '("lines lt -1")))) ;Style for dumping
	(apply
	 #'gnuplot (1+ (length lines))
	 5
	 #'(lambda (plot point)
	     (if (eq point :title)
		 (if (= plot 9)
		     (format nil "T = ~A" (format-float-reasonably *dump-time*))
		   (nth plot '("diamond" "I=0" "J=0" "K=0" "L=0" "I=1" "J=1" "K=1" "L=1")))
	       (cond ((zerop plot)	;The diamond itself
		      (values-list (3vector-list (nth point points))
		      )
		      )
		     ((> point 1) nil)
		     (t  (values-list (3vector-list
					(nth point (nth (1- plot) lines))))))))
	 :prelude (format nil "~:{set label '~C' at ~{~F, ~F, ~F~} offset character 0.3,0.3~%~}"
			  (list (list #\s (3vector-list (diamond-start diamond)))
				(list #\l (3vector-list (diamond-left diamond)))
				(list #\r (3vector-list (diamond-right diamond)))
				(list #\e (3vector-list (diamond-end diamond)))))
	 :styles styles
	 :xrange (first ranges)
	 :yrange (second ranges)
	 :zrange (third ranges)
	 gnuplot-keys)))))

	  

;;Plot a set of diamonds
(defun plot-diamonds (diamonds &rest gnuplot-keys)
  (let* ((points-list			;List of lists of boundary points
	  (mapcar
	   #'(lambda (diamond)		;Get list of points around diamond
	       (let ((points (list (diamond-left diamond)
				   (diamond-start diamond)
				   (diamond-right diamond)
				   (diamond-end diamond))))
		 (push (car (last points)) points)))	;Make loop
	   diamonds)))
    (apply #'gnuplot (length points-list) (loop for points in points-list maximize (length points))
	   #'(lambda (plot point)
	       (if (eq point :title) (format nil "~S" (nth plot diamonds))
		 (let ((location (nth point (nth plot points-list))))
		   (and location (values-list (3vector-list location))))))
	   :styles :linespoints
	   gnuplot-keys)))


;;Plots intersection of a given time and the strings that we know about.
(defun plot-time-slice (&rest gnuplot-keys &key (time (current-time)) cell (style "lines lt 1")
			      marks 	;Mark where loops were deleted
			      range
			      title
			      &allow-other-keys)
  (let ((data (get-plot-data :time time :range range))
	(outline (and cell (unit-cell-outline))))
    (format t "~&~D string~:P with ~D segments in all~%" (length data)
	    (reduce #'+ (map 'list #'length data)))
    (setq data (coerce data 'simple-vector))
    (let* ((length (length data))	;Number of segments to plot
	   (mark-plot (and marks length)) ;Plot # for marks
	   (outline-plot (and outline (if marks (1+ length) length))) ;Plot # for outline
	   (styles (make-list length :initial-element style))) ;All same color
      (when mark-plot
	(setq styles (append styles '("points lt 3"))))
      (when outline
	(setq styles (append styles '("lines lt 2"))))
      (apply #'gnuplot
	     (+ (if outline 1 0) (if marks 1 0) length)
	     (max (* 4 (length outline)) (* 3 (length marks)) (loop for string across data maximize (length string)))
	     #'(lambda (plot point)
		 (unless (eq point :title)
		   (let ((vector
			  (cond ((eq plot outline-plot)	;Outline
				 (multiple-value-bind (line part) (truncate point 4) ;Two blank lines between segments
				   (nth part (nth line outline))))
				((eq plot mark-plot)
				 (multiple-value-bind (line part) (truncate point 3) ;Two blank lines between segments
				   (and (zerop part)
					(nth line marks))))
				(t		;Regular data
				 (let* ((string (aref data plot))
					(length (length string)))
				   (and (< point length) (aref string point)))))))
		     (and vector
			  (values-list (3vector-list vector))))))
	     :title (or title (format nil "Time ~A" (format-float-reasonably time)))
	     :styles styles
	     :3d t
	     (append (and range
			  (list :xrange (cons (3vector-x (first range)) (3vector-x (second range)))
				:yrange (cons (3vector-y (first range)) (3vector-y (second range)))
				:zrange (cons (3vector-z (first range)) (3vector-z (second range)))))
		     gnuplot-keys)))))


;;Get data for plotting: a list of vectors of positions
(defun get-plot-data (&key (time (current-time)) range velocities-p)
  (let ((data nil))
    (map-string-paths
     #'(lambda (diamond &rest rest)
	 (let* ((path (apply #'get-path time diamond rest)) ;Get path of string
		(positions (get-plot-positions path (first range) (second range) velocities-p))) ;Get positions
	   (unless (zerop (length positions)) ;Empty string or all outside range
	     (push positions data)))))
    data))


;;Make a cubical range for plotting
(defun make-range (min max)
  (list (make-3vector min min min) (make-3vector max max max)))

(defun plot-unit-cell ()
  (let ((data (unit-cell-outline)))
    (gnuplot 1 (* 4 (length data))
	     #'(lambda (plot point)
		 (declare (ignore plot))
		 (unless (eq point :title)
		   (multiple-value-bind (line part) (truncate point 4) ;Two blank lines between segments
		     (let ((result (nth part (nth line data))))
		       (and result (values-list (3vector-list result)))))))
	     :3d t
	     :styles :lines)))


;;Accepts a list of (DIAMOND A1 B1 A2 B2)
;;Returns a vector of global positions to be plotted, with 2 consecutive NIL's where there is a break.
;;This function works whether or not the string is a loop
;;If velocities set, put velocity between each pair of positions.
(defun get-plot-positions (string &optional minima maxima velocities-p)
  (let ((data (make-array 100 :fill-pointer 0)))
    (loop for (diamond a1 b1 a2 b2) in string
	  with previous = nil
	  for start = (diamond-position-wrap diamond a1 b1) ;Start of the segment
	  for end = (diamond-position-wrap diamond a2 b2) ;End of the segment
	  do
	  (when minima			;only including range?
	    (multiple-value-setq (start end) (segment-inside-range start end minima maxima)))
	  (when start			;Nothing to do if string was entirely outside
	    (cond ((null previous)	;First time: start with start point
		   (vector-push-extend start data))
		  ((3vector= start previous fudge-global-coordinates) ;Previous end same as start?
		   )			;Don't need start point
		  (t			;There is a previous, but it is not the same as this point,
		   (vector-push-extend nil data) ;Break line here and start again
		   ;(vector-push-extend nil data) ;gnuplot splot requires two blank lines or else you get cross-lines
		   (vector-push-extend start data)))
	    ;;Always use end point.  If there's a break between the first start and the last end, we need both.
	    ;;If not, we still need both, so that the loop will close in the plot
	    (when velocities-p
	      (vector-push-extend (diamond-velocity diamond) data))
	    (vector-push-extend end data)
	    (setq previous end)))
    data))

;;Return the portion of the segment that lies inside the given box, if any
(defun segment-inside-range (start end minima maxima)
  (if (point-in-box start minima maxima 3) ;Start is inside?
      (if (point-in-box end minima maxima 3) ;End too?
	  (values start end)		;Yes, nothing to do
	(values start (line-box-intersection start end minima maxima)))
    (if (point-in-box end minima maxima 3) ;Start is outside.  What about end?
	(values (line-box-intersection end start minima maxima) end) ;End inside
      (let ((results (line-box-intersections end start minima maxima)))
	(when results			;NIL if doesn't intersect
	  (unless (= (length results) 2)
	    (error "We expected 2 intersections or none, but there were ~D" (length results)))
	  (values-list results))))))

(defun line-box-intersection (from to minima maxima)
  (let ((results (line-box-intersections from to minima maxima)))
    (unless (= (length results) 1)
      (error "We expected 1 intersection, but there were ~D" (length results)))
    (first results)))		 

;;Return place or places where line intersects box
(defun line-box-intersections (from to minima maxima)
  (let ((results nil))
    (dotimes (direction 3)		;Check faces perpendicular to this direction
      (let ((denominator (- (aref to direction) (aref from direction))))
	(unless (< (abs denominator) fudge-global-coordinates) ;Parallel to edge: say no intersection
	  (flet ((try (value)
		   (let ((parameter (/ (- value (aref from direction)) denominator)))
		     (when
			 (and (< 0.0 parameter 1.0) ;Valid parameter range
			      ;;Find position of intersection
			      (let ((position (4vector+ (4vector-scale from (- 1.0 parameter))
							(4vector-scale to parameter))))
				(when (loop for index below 3
					    always (or (= index direction) ;Don't check direction we just solved
						       (and (> (+ (aref position index) fudge-global-coordinates)
							       (aref minima index))
							    (< (- (aref position index) fudge-global-coordinates)
							       (aref maxima index)))))
				  (unless (loop for previous in results ;Already there?
						thereis (4vector= previous position fudge-global-coordinates))
				    (push position results)))))))))
	    (try (aref minima direction))
	    (try (aref maxima direction))))))
    results))

;;Get diamond position, taking into account the possibility that the diamond may have wrappings between corners
;;The position is standardized relative to the start
;;All diamond corners should have been standardized on reading
(defun diamond-position-wrap (diamond a b)
  (let* ((start (diamond-start diamond)) ;Wrap all corners so that they are near start
	 (end (standardize-position (diamond-end diamond) start))
	 (left (standardize-position (diamond-left diamond) start))
	 (right (standardize-position (diamond-right diamond) start)))
    ;;Now the diamond does not have any wrappings between its corners.  Get position the usual way
    (cond ((and (= a 0.0) (= b 0.0)) start)	;Special cases are important to overlap-diamond
	  ((and (= a 0.0) (= b 1.0)) right)
	  ((and (= a 1.0) (= b 1.0)) end)
	  ((and (= a 1.0) (= b 0.0) left))
	  (t (let ((a1 (- 1.0 a))
		   (b1 (- 1.0 b)))
	       (4vector+ (4vector-scale start (* a1 b1))
			 (4vector-scale end (* a b))
			 (4vector-scale left (* a b1))
			 (4vector-scale right (* b a1))))))))

;;Diamond position, but wrap if we are reading dumps
(defun diamond-position-wrap-dumps (diamond &key a b)
  (if *reading-dumps* (diamond-position-wrap diamond a b)
    (diamond-position diamond :a a :b b)))

;;Edges of diamond with wrapping
(mirror-images
(defun diamond-a-wrap (diamond)
  (let ((start (diamond-start diamond)))
    (4vector- (standardize-position (diamond-left diamond) start) start)))

(defun diamond-new-a-wrap (diamond)
  (let ((right (diamond-right diamond)))
    (4vector- (standardize-position (diamond-end diamond) right) right))))

;;Give the amount of wrapping as we go over the loop
(defun loop-wrapping (start-diamond)
  (loop with start = (diamond-right start-diamond)
	with position = start
	for diamond = (diamond-e start-diamond) then (diamond-e diamond) ;Walk over loop to east
	do (setq position (standardize-position (diamond-right diamond) position)) ;Standardize each w.r.t. last
	until (eq diamond start-diamond)
	finally (return (3vector- position start))))

;;Get velocity vector of this diamond (at the start, in case of expansion), taking into account wrapping
(defun diamond-velocity (diamond)
  (let* ((start (standardize-position (diamond-start diamond))) ;Wrap corners so that they are near start
	 (left (standardize-position (diamond-left diamond) start))
	 (right (standardize-position (diamond-right diamond) start)))
    ;; v = (unit p + unit q)/2.  Normalize each to 1/2, then add
    (if (and (zerop (3vector-length (3vector- left start)))
	     (zerop (3vector-length (3vector- right start))))
	0.0
      (if (zerop (3vector-length (3vector- left start)))
	  (3vector-normalize (3vector- right start) 0.5)
	(if (zerop (3vector-length (3vector- right start)))
	    (3vector-normalize (3vector- left start) 0.5)
	  (3vector+ (3vector-normalize (3vector- left start) 0.5) (3vector-normalize (3vector- right start) 0.5)))))))

;;Read and plot a snapshot
;;LOOP-POSITIONS is an array of loop positions or T to read them
(defun plot-read-dumps (directory time &rest gnuplot-keys &key loop-positions &allow-other-keys)
  (when loop-positions
    (unless (or (integerp time) (= time (float (round time) time)))
      (error "Time ~F is not an integer" time))
    (setq time (round time))
    (when (eq loop-positions t)
      (setq loop-positions (read-loop-positions directory (1- time) time)))
    (setq gnuplot-keys (append `(:marks ,(aref loop-positions time)) gnuplot-keys)))
  (read-dumps directory :time time)
  (apply #'plot-time-slice :time time gnuplot-keys))

(defun make-movie (directory start end step &rest gnuplot-keys)
  (loop for time from start by step below end
	for count from 0
	for title = (format nil "Time ~A" (format-float-reasonably time))
	do (format t "~A~%" title)
	do (apply #'plot-read-dumps directory time
		  :prelude (format nil "set output '~A/plot-~D.jpg'~%set terminal jpeg~%"
				   directory count)
		   :title title :key nil
		  :exit t
		  gnuplot-keys)))

(defun make-movie-submit (directory start end step &rest gnuplot-keys)
  (setq directory (merge-pathnames directory batch-root-directory))
  (loop for time from start by step below end
	for count from 0
	for title = (format nil "Time ~A" (format-float-reasonably time))
	for file = (format nil "~A/plot-command-~D.lisp" directory count)
	do (with-group-write-access
	    (with-open-file (stream file :direction :output :if-exists :supersede :if-does-not-exist :create)
	      (write
	       `(plot-read-dumps ,directory ,time
				 :prelude ,(format nil "set output '~A/plot-~D.jpg'~%set terminal jpeg~%"
						   directory count)
				 :title ,title :key nil
				 ;;				:loop-positions t
				 :exit t
				 ,@(loop for item in gnuplot-keys
					 collect `',item))
	       :stream stream)))		;Put command into file
	(do-submit directory (format nil "plot-~D" count) (format nil "plot-~D-" count) batch-flags
		   `(load ,file))))

;;Histogram loop sizes read from dump
;;We plot dL/d(ln l), where dL and dn are the total length and the number of loops whose lengths are between l and l+dl
;;so dL = l dn.
;;According to Alex and Tanmay, dn/dl goes as l^{-5/2}, so dn/d(ln l) goes as l^{-3/2}, and dL/dl as l^{-1/2}
(defun histogram-loop-sizes (&rest gnuplot-keys &key (lengths (path-lengths)) (bins 20) (min 1.0) (max *total-size*)
				   &allow-other-keys)
  (setq max (float max 0.0))
  (let* ((bin-size (/ (log (/ max min)) bins)) ;Logarithmic size of a bin
	 (data (make-array bins :element-type 'double-float :initial-element 0.0)) ;Total length in this bin
	 (infinite-total 0.0)
	 (finite-total 0.0)
	 (short-total 0.0)
	 )
    (dolist (length lengths)
      (declare (double-float length))
      (if (> length min)
	  (let ((index (fixnum-floor (/ (real-log (/ length min)) bin-size))))
	    (cond ((>= index bins)	;Too large?
		   (incf infinite-total length))
		  (t
		   (incf finite-total length)
		   (incf (aref data index) (/ length bin-size (total-3volume)) ;Convert to dL/(dV d(ln l))
			 ))))
	(incf short-total length)))
    (format t "~&Length in short loops: ~$, finite: ~$ infinite: ~$, total ~$"
	    short-total finite-total infinite-total (+ short-total finite-total infinite-total))
    (flet ((fit (l)			
		(* 0.4 (expt l -0.5) ;Fit from paper, but 0.2 seems to be closer to what we have
		   ;;In VV paper, all cells have unit length of string.  But we have only about 1/sqrt{2} if the
		   ;;string turns a corner.  This computes our average length. 0.13 and 0.22 should be rationals
		   (/ (+ (* (sqrt 0.5) 0.22) 0.13) (+ 0.22 0.13))
		   )))
		   
      (apply #'gnuplot
	     2 bins
	     #'(lambda (plot point)
		 (if (eq point :title) (nth plot '("data" "VV paper result"))
		   (ecase plot
		     (0 (let ((length (* min (exp (* bin-size (+ point 0.5)))))	;Central length value
			      (total (aref data point))) ;Total length in loops in this bin
			  (and (plusp total) ;no point, rather than 0 if no data
			       (values-list (list length total)))))
		     (1 (case point
			  (0 (values min (fit min)))
			  (1 (values max (fit max)))
			  (t nil))))))
	     :styles :linespoints
	     :logscale '(:x :y)
	     gnuplot-keys))))

;;New spectrum plotting

;;Parameter m in nonuniform FT.  See Fourier.tex.
(defparameter *nufft-m* 8)


(mirror-images
;;Advance (west for A, east for B) to the next diamond with a different A or different B
;;If you do this for both A and B from the same point, you get places that will see each other in the future.
;;This returns a diamond that has a SE rather than a NE.  Thus if there are several diamonds that have the same
;;A sides, this returns the latest in time.  If reconnections
;;decreased the length of the side, this is the shortest one.
;;If we reach something that is not a diamond, return it
;;The meanings of PREVIOUS and NEXT were reversed in November 2015.
(defun diamond-next-a (diamond)
  (if (diamondp diamond)
      (or (diamond-nw diamond)		;New A?  Return it
	  (diamond-next-a (diamond-sw diamond)))
    diamond))

;;Reverse direction: east for A, west for B.  This gives places that have seen each other in the past.
;;if there are several diamonds that have the same A sides, this returns the earliest in time.  
(defun diamond-previous-a (diamond)
    (or (diamond-se diamond)		;New A?  Return it
	(diamond-previous-a (diamond-ne diamond)))) ;Same A: recurse

;;Returns values:
;; A-HATS: SIMPLE-VECTOR with the values of the unit tangent 3-vector a'
;; DSIGMAS: (SIMPLE-ARRAY DOUBLE-FLOAT (*)) with length of corresponding a'
;; CREATED: (SIMPLE-ARRAY DOUBLE-FLOAT (*)) when kink between this a' and next was generated.
;;The A-HAT vectors are in order going to the west, meaning that each vector is applied after the previous one.
;;In the mirror image case, the B-HAT vectors go around to the east, with the same meaning.
;;We work along the future edge of the existing world sheet, so the total spatial change in a and the total change in b
;;are equal for a closed, non-wrapping loop, and for such a loop in the rest frame, they are both zero.
;;This also means that if there are intersections in the past, we will not be confused by them.
;;If SKIP is 0 (the default), the first tangent vector will be the one for the current diamond, or its furthest
;;NE neighbor (which will be the same unless there are intersections).
;;If SKIP > 0, we will skip forward the given number of a' values relative to that.  If SKIP < 0 we skip backward
;;If this is called for A and B with the same parameters, and SKIP > 0, you get segments that have not seen each
;;other.  If SKIP < -COUNT, you get segments that have seen each other.
;;The order of returned vectors and thus the meaning of positive and negative SKIP were reversed in November 2015.
;;If CLOSE is set, we will adjust the vectors to make sure that the total of hat*sigma is the same as the
;;original total of the diamond-a's, which might not otherwise happen because of numerical errors.
(defun get-a-data-dsigmas (diamond &key count skip (close t))
  (unless skip (setq skip 0))
  (when count (setq close nil))	;Don't bother if only using part of the loop
  ;;We have to make sure here to get a diamond that begins a new a going westward, so we will find it again if loop
  (cond ((plusp skip)			;Skipping forward (to left)
	 (loop repeat skip do (setq diamond (diamond-next-a diamond)))) ;Skip this many.
	(t			;Backward or no skip.  First go back one more than skip requested.
	 (loop repeat (- 1 skip) do (setq diamond (diamond-previous-a diamond)))
	 (setq diamond (diamond-next-a diamond)))) ;Now forward so we are at the beginning of a stretch of a
  (loop with this = diamond			   ;First count a's in loop
	do (setq this (diamond-next-a this))	   ;Move to next diamond
	count t into loop-count
	when (eq this diamond)			  ;Looped?  Each A in loop counted once
	return (setq count loop-count)		  ;Found loop before count ran out
	until (and count (= count loop-count)))	  ;Reached max without looping?  Use given count
  (loop with a-hats = (make-array count) ;Now we'll store exactly COUNT data
	with dsigmas = (make-array count :element-type 'double-float)
	with created = (make-array count :element-type 'double-float)
	with total-a = (make-zero-3vector)
	with this = diamond
	for index from 0 below count
	for a = (diamond-new-a-wrap this) ;Future edge
	for dsigma = (* 2 (3vector-t a))  ;Use time as sigma, so loop closes in time properly
	unless (plusp dsigma)
	do (warn "Zero length segment in ~S" this)
	unless (fudge= (3vector-length a) (3vector-t a) fudge-coordinates)
	do (warn "~D has non-null side, length=~F, t=~F, difference ~S" this (3vector-length a)  (3vector-t a)
		 (/ (- (3vector-t a) (3vector-length a)) (3vector-t a)))
	do
	(3vector-incf total-a a)			 ;Get 3vector total of sides
	(setf (aref a-hats index) (3vector-normalize a)) ;Store present a'.  Normalize so it is always unit
	(setf (aref dsigmas index) dsigma) ;Store length
	(setf (aref created index) (4vector-t (diamond-a-kink-created this))) ;Kink from this to new diamond
	(setq this (diamond-next-a this)) ;Move to next diamond
	finally	(return
		 (values (if close
			     (close-hats a-hats :dsigmas dsigmas ;Fix up hats to give right total a
					 :goal (3vector-scale total-a 2.0)) ;double because a side is half hat*dsigma
			   a-hats)
			 dsigmas created))))

(defun get-a-data-sigmas (diamond &key count skip smooth-range smooth-interval)
  (multiple-value-bind (hats dsigmas created) (get-a-data-dsigmas diamond :count count :skip skip)
    (if smooth-range			;Smooth out data?
	(multiple-value-bind (new-hats new-sigmas)
	    (smooth-data hats (dsigma-to-sigma dsigmas) smooth-range smooth-interval)
	  (values new-hats new-sigmas created))
      (values hats (dsigma-to-sigma dsigmas) created))))
  
;;Get just the a-hats, checking that the a-dsigmas are equal
(defun get-a-hats (diamond &key (tolerance fudge-coordinates))
  (multiple-value-bind (a-hats a-dsigmas) (get-a-data-dsigmas diamond)
    (loop with average = (/ (total-sigma a-dsigmas) (length a-dsigmas))
	  for dsigma across a-dsigmas
	  minimize dsigma into min
	  maximize dsigma into max
	  finally
	  (let ((spread (/ (- max min) average)))
	    (when (> spread tolerance)
	      (error "Fractional spread of ~A delta-sigma values is ~S" :a spread))))
    a-hats))

;;Check a loop for diamonds whose sides aren't null
(defun check-non-null-a (diamond &key (tolerance fudge-coordinates))
  (setq diamond (diamond-next-a diamond)) ;forward so we are at the beginning of a stretch of a
  (loop with worst = nil
	and worst-error
	with this = diamond
	for a = (diamond-new-a-wrap this) ;Future edge
	for time = (3vector-t a)
	for length = (3vector-length a)
	for error = (abs (- time length))
	count t into segments
	unless (fudge= time length tolerance)
	count t into bad
	and sum error into total-error
	and when (or (null worst) (> error worst-error))
    	    do (setq worst this worst-error error)
	do (setq this (diamond-next-a this)) ;Move to next diamond
	until (eq this diamond)
	finally (format t "Out of ~D ~A segments, ~[all are null~:;~:*~D are non-null with average error ~S, worst ~S~]~%"
			segments :a bad (and (plusp bad) (/ total-error bad)) worst-error)
	(return worst)))

;;Return length sigma and count of segments in A
(defun loop-length-and-count-a (start-diamond)
  (setq start-diamond (diamond-next-a start-diamond)) ;Advance to diamond that starts a new A
  (loop with diamond = start-diamond
	do (setq diamond (diamond-next-a diamond)) ;Get next A segment
	count t into count
	;;The vector q points, say, from sigma=t=0 to sigma=t=q_t, since sigma-t is fixed along the right edge
	;;Thus the argument of b(t+sigma) is 2 q_t.  Check: x=(a+b)/2 so delta-x = delta-b/2
	;;= (b(sigma_1) - b(sigma_0))/2 = b-hat (sigma_1-sigma_0)/2 so sigma_1-sigma_0 = 2 |delta-x|
	sum (* 2 (3vector-length (diamond-a-wrap diamond))) into length
	when (eq diamond start-diamond)	;Looped around?  This last segment has been done
	return (values length count)))

;;Return total of <(a - prev-a)^2> and count of segments
(defun loop-length-and-angle-a (start-diamond)
  (setq start-diamond (diamond-next-a start-diamond)) ;Advance to diamond that starts a new A
  (loop with diamond = start-diamond
	with a with dsigma with a-hat
	for previous-a = (3vector-normalize (diamond-a-wrap start-diamond)) then a-hat
	do (setq diamond (diamond-next-a diamond) ;Advance to new diamond
		 a (diamond-a-wrap diamond)	 ;This a' vector
		 dsigma (* 2 (3vector-length a))	   ;Amount of sigma here.  See loop-length-and-count-a
		 a-hat (3vector-scale a (/ 2 dsigma))) ;unit tangent vector
	sum dsigma into length
	sum (* 2 (- 1.0 (3vector-dot previous-a a-hat))) into total-angle ;(a1 - a2)^2 = 2(1 - a1.a2)
	when (eq diamond start-diamond)	;Looped around?  This last segment has been done
	return (values length total-angle)))

;;Find power spectrum of a single string loop, as defined by Fourier.tex.
;;Returns array with P(k) where element i corresponds to frequency i/l, where l is the length of the loop.
;;We always return 0 in slot 0, because P(0)=0 by definition
;;If EXTEND-K is given, we return frequencies through that value
(defun loop-power-a (start-diamond &optional extend-k)
  (setq start-diamond (diamond-next-a start-diamond))			      ;Advance to diamond that starts a new A
  (multiple-value-bind (l nn) (loop-length-and-count-a start-diamond) ;Number of segments N in loop and total length L
    (declare (double-float l)
	     (fixnum nn))
    (let* ((needed (if extend-k (max nn (ceiling (* l extend-k) (* 2 pi)))
		     nn))		;Normally give N frequencies, but can request more
					;Size of FFT: At least twice needed frequencies, rounded up to power of 2
	   (n (expt 2 (1+ (ceiling (log needed 2)))))
	   (xi (make-array nn :element-type 'double-float)) ;Point locations for Fourier transform
	   (fi (make-array nn :element-type 'double-float)) ;Function values for Fourier transform
	   ;;Our calculations will yield n frequencies ranging from -(n-1)/(2L) to n/(2L).  But negative ones are just
	   ;;the conjugates of positive ones and those in the larger half are not well approximated by our procedure.
	   ;;So we return n/4 powers.
	   (power (make-array (/ n 4) :element-type 'double-float :initial-element 0.0)))
      (locally
	  (declare (optimize speed)
	       (type (simple-array double-float (*)) xi fi power)
	       (fixnum n))
      (dotimes (component 3)
	(loop with sigma double-float = 0d0		;sigma = 0 represents the end of this segment
	      with diamond = start-diamond
	      with a
	      and dsigma double-float
	      and a-hat
	      with index = 0
	      with previous-a of-type 3vector = (3vector-normalize (diamond-a-wrap start-diamond))
	      do (setq diamond (diamond-next-a diamond) ;Advance to new diamond
		       a (diamond-a-wrap diamond) ;This a' vector
		       dsigma (* 2 (3vector-length a)) ;Amount of sigma here.  See loop-length-and-count-a
		       a-hat (and (plusp dsigma) (3vector-scale a (/ 2 dsigma)))) ;unit tangent vector
	      when a-hat
	        ;;We come here with each pair of diamonds.  The last time, diamond is the starting diamond
	        ;;sigma is the coordinate at which one diamond gives way to the next
	        do (setf (aref xi index) (/ sigma l) ;Scale into 0..1
			 (aref fi index) (- (aref a-hat component) (aref previous-a component))) ;FT of this: see notes
		and do (incf index)
		and do (incf sigma dsigma)	;Advance over present diamond
		and do (setq previous-a a-hat)
	      else do (warn "Zero length segment")
	      until (eq diamond start-diamond) ;Stop when start diamond has just been processed
	      )
	(let ((spectrum (nufft xi fi *nufft-m* n))) ;Get Fourier components
	  (declare (type (simple-array double-float (*)) spectrum))
	  ;;Frequency j in SPECTRUM (stored in slots 2j, 2j+1) corresponds to wave number k = 2 pi j/L.
	  ;;The value is sum_i [a'_i - a'_{i-1}] e^{-i k sigma_i}.
	  ;;Thus A(k) = 1/sqrt{2piL} (1/(ik)) spectrum(j).
	  ;;We are trying to find P = 2 k |A|^2 =  |spectrum(j)|^2 / (pi L k) = |spectrum(j)|^2 / (2 pi^2 j)
	  (loop for j from 1 below (the fixnum (/ n 4))	;Skip freq 0
		do (incf (aref power j)
					;squared magnitude of complex number stored as 2 reals
			 (/ (+ (expt (aref spectrum (* 2 j)) 2) (expt (aref spectrum (1+ (* 2 j))) 2))
			    (* 2 pi pi j))))))
    (values power l)))))

)					;mirror-images

(defun get-ab-data-sigmas (diamond)
  (multiple-value-bind (a-hats a-sigmas) (get-a-data-sigmas diamond)
    (multiple-value-bind (b-hats b-sigmas) (get-b-data-sigmas diamond)
      (values a-hats a-sigmas b-hats b-sigmas))))

(defun get-ab-data-dsigmas (diamond)
  (multiple-value-bind (a-hats a-sigmas) (get-a-data-dsigmas diamond)
    (multiple-value-bind (b-hats b-sigmas) (get-b-data-dsigmas diamond)
      (values a-hats a-sigmas b-hats b-sigmas))))

;;Get hats and sigmas from loop and call function with given arguments
(defun call-with-hats (diamond function &rest args)
  (multiple-value-bind (a-hats a-sigmas b-hats b-sigmas) (get-ab-data-sigmas diamond)  
    (apply function a-hats a-sigmas b-hats b-sigmas args)))

;;Here when the function is supposed to be applied to the arguments
(defun apply-with-hats (diamond function &rest args)
  (multiple-value-bind (a-hats a-sigmas b-hats b-sigmas) (get-ab-data-sigmas diamond)  
    ;;args is a list of arguments, the last one of which is supposed to be spread, so spread ARGS and call APPLY to spread last arg
    (apply #'apply function a-hats a-sigmas b-hats b-sigmas args)))

;;Add one power spectrum to bins.  Add P(k) to DATA, k values to KS.
;;K here is always k, not kt.
(defun bin-power-1 (data counts ks k-min k-max k-bins power l)
  (loop with bin-size = (/ (log (/ k-max k-min)) k-bins)
	for j from 1 below (length power) ;Skip freq 0.  Index j corresponds and wave number k = 2 pi j/L.
	for p = (aref power j)
	for k = (/ (* 2 pi j) l)
	for bin = (floor (/ (log (/ k k-min)) bin-size))
	when (and (>= bin 0) (< bin k-bins))
	do
	(incf (aref data bin) p)
	(incf (aref counts bin))
	(incf (aref ks bin) k)))

;;Spectrum in high-K limit is <(a_i - a_{i+1})^2> / (pi k <seg size>) = (total (a-a)^2)/(pi total-length)
;;We return the coefficient
(defun high-k-power (&key (min-length 0.0) (scaling t))
  (let ((total-length 0.0)
	(total-angle 0))		;Total <a . a>
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Ignore open strings
     #'(lambda (diamond)		;Loops
	 (multiple-value-bind (l angle) (loop-length-and-angle-a diamond)
	   (when (> l min-length)
	     (incf total-length l)
	     (incf total-angle angle)
	     (multiple-value-bind (lb angle) (loop-length-and-angle-b diamond)
	       (unless (< (/ (- l lb) l) 0.01)
		 (warn "A length is ~F, but B length is ~F" l lb)
		 (incf total-length l)
		 (incf total-angle angle)))))))
    (let ((result (/ total-angle pi total-length)))
      (if scaling (* result (current-time)) result)))) ;In scaling case, user will want to divide by kt
	  
(mirror-images
;; Fourier series of a'(sigma)
(defun slow-right-moving-k-spectrum (path k)
  (loop for component below 3
	for sigma = 0d0
	collect (/ (loop for segment in path
			 for dL = (apply #'path-segment-length segment)
			 for p-hat = (3vector-normalize (4to3vector (diamond-a (car segment))))
			 do (incf sigma dL)
			 sum (* (aref p-hat component)  
				(/ (exp (* #C(0.0 -1.0)
					   k
					   sigma))
				   (* #C(0.0 1.0) k))
				(- (exp (* #C(0.0 1.0)
					   k
					   dL))
				   1.0)))
		   (sqrt (* 2 pi sigma)))))


;; power spectrum of path for wavenumber k
(defun slow-right-moving-power (path k)
  (let ((spectrum (slow-right-moving-k-spectrum path k))
	)
    (loop for component below 3
	  sum (abs (* 2
		      k
		      (nth component spectrum)
		      (conjugate (nth component spectrum)))))))

)

;;Get power for single wavenumber by direct integration
(defun slow-power (k &optional (min-length 3.0))
  (let ((total-power 0)
	(total-length 0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Ignore open strings
     #'(lambda (diamond)
	 (let* ((path (get-path (current-time) diamond))
		(length (path-length path)))
	   (when (> length min-length)
	     (incf total-length length)
	     (mirror-images		;Add right and left power
	      (incf total-power (* (slow-right-moving-power path k) length))) ;Weight by length
	     ))))
    (/ total-power total-length 2)	;Extra 2 for A and B
    ))


(defstruct power
  k-min
  k-max
  extend-k
  data					;Average power in each bin
  ks)					;Average k value in each bin

;;Read dumps, process power and two-point-function as requested
;;This works by calling start-process-power and start-process-two-point-function to get some closures which are
;;then called twice for each loop in the dump, with the loop, its length, the a- or b-power, and the total a or b,
;;and once at the end with no arguments, to say were finished
;;If EXTEND-K, we Fourier transform enough points to give that value of k.  Otherwise we use only as many
;;frequencies as there are points, rounded up to a power of 2.
(defun process-dumps (directory time
				&key (power t) (two-point-function t) ;What to do
				(delta-min 1d-3) (delta-max 1d5) (two-point-bins 80)
				(k-min 1d-3) (k-max 1d2) (power-bins 80)
				extend-k ;This applies to both operations
				min-length)
  (setq directory (merge-pathnames directory batch-root-directory))
  (cond ((or power two-point-function)
	 (read-dumps directory :time time)
	 (unless min-length (setq min-length *total-size*)) ;Now that we have read the dump, we can default this
	 (with-group-write-access
	  (let ((handlers nil))
	    (when power
	      (push (start-process-power directory time k-min k-max power-bins extend-k) handlers))
	    (when two-point-function
	      (push (start-process-two-point-function directory time delta-min delta-max two-point-bins extend-k)
		    handlers))
	    (map-string-paths
	     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Ignore open strings
	     #'(lambda (diamond)		;Loops
		 (let ((l (loop-length diamond)))
		   (when (> l min-length)
		     (format t "Length ~$: " l)
		     (mirror-images	    ;A and B
		      (format t "~A " :a) (force-output)
		      (let ((power (loop-power-a diamond extend-k))
			    (total-a (loop-total-a diamond)))
			(mapc #'(lambda (handler) (funcall handler diamond l power total-a)) handlers)))))))
	    (mapcar #'funcall handlers)
	    (terpri))))
	 (t (warn "PROCESS-DUMPS called with nothing to do"))))


;;Process all times available in this directory
(defun process-dumps-all (directory &rest keys)
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((times (dump-times directory)))
    (format t "Times: ~S~%" times)
    (dolist (time (reverse times))	;Latest to earliest, so memory problems occur later
      (format t "Time ~S:~%" time)
      (apply #'process-dumps directory time keys))))

(defun start-process-power (directory time k-min k-max k-bins extend-k)
  (let ((data (make-array k-bins :element-type 'double-float :initial-element 0d0))
	(counts (make-array k-bins :element-type 'fixnum :initial-element 0))	   ;Number of samples in each bin
	(ks (make-array k-bins :element-type 'double-float :initial-element 0d0))) ;Total of k values in each bin
    #'(lambda (&optional diamond l power total-a)				   ;Loop to process, or no arg when done
	(declare (ignore total-a))						   ;Only used by two-point-function
	(if diamond
	    (bin-power-1 data counts ks k-min k-max k-bins power l)
	  ;;Finished: write file
	  (with-open-file (stream (power-file time directory) :direction :output
				   :if-does-not-exist :create :if-exists :supersede)
	    (loop for index below k-bins	;Average power
		  for count = (aref counts index)
		  when (plusp count)
		  do (setf (aref data index) (/ (aref data index) count)
			   (aref ks index) (/ (aref ks index) count)))
	    (prin1 (make-power :k-min k-min :k-max k-max :extend-k extend-k
			       :data data :ks ks)
		   stream)
	    (pathname stream))))))

;;Read processed data from file
(defun read-power (directory time)
  (with-open-file (stream (power-file time directory)) (read stream)))
    
;;Accept list of (DIRECTORY TIME &OPTIONAL TITLE K-UNIT) or just DIRECTORY and plot functions
;;In the latter case we use all times.
;;If no TITLE is given, we use the directory name and time
;;If LENGTH-SCALE is given, k values in data are multiplied by it.
(defun plot-power-1 (list &rest gnuplot-keys &key scaling bins &allow-other-keys)
  (setq list (loop for entry in list
		   append (if (consp entry) (list entry)
			    (let* ((directory (merge-pathnames entry batch-root-directory))
				   (times (dump-times directory)))
			      (unless times
				(error "No dumps were made in ~A" directory))
			      (mapcar #'(lambda (time) (list directory time)) times)))))
  (loop with arbitrary-scaling = nil
	for (directory time title length-scale) in list
	for power = (read-power (merge-pathnames directory batch-root-directory) time)
	when bins
	do (setq power (rebin-power power bins))
	collect (list power time length-scale (or title (format nil "~A: ~D" directory time))) into data-list
	when length-scale do (setq arbitrary-scaling t)
	maximize (length (power-data power)) into max-points
	finally
	(apply #'gnuplot (length data-list) max-points
	       #'(lambda (plot point)
		   (destructuring-bind (power time length-scale title) (nth plot data-list)
		     (if (eq point :title) title
		       (when (< point (length (power-data power)))
			 (let* ((p (aref (power-data power) point))
				(k (aref (power-ks power) point)))
			   (unless (zerop p)			   ;don't plot points with no data
			     (if scaling (setq k (* k time))) ;Convert to kt for horizontal axis
			     (if length-scale (setq k (* k length-scale))) ;Multiply by arbitrary length unit
			     (values k p)))))))
	       :xlabel (cond (arbitrary-scaling "given units") (scaling "kt") (t "k"))
	       :ylabel "correlator"
	       :styles :linespoints
	       :logscale :x
	       ;;Put lambda = 2pi/k on opposite axis.
	       :prelude (format nil "~@{~A~%~}" "set xtics nomirror" ;Get rid of mirrored tics on top.
				"set link x2 via log10(2*pi)-x inverse log10(2*pi)-x" ;Map logarithm of x
				"set x2tics"					  ;Label tics on top
				(format nil "set x2label '~A'"
					(cond (arbitrary-scaling "2pi/given") (scaling "lambda t") (t "lambda"))))
	       gnuplot-keys)))

;;Merge data into fewer bins to reduce fluctuations in graphs
;;Return new power struct
(defun rebin-power (power bins)
  (let ((old-bins (length (power-data power))))
    (if (= old-bins bins) power		;nothing to do
      (let ((factor (floor old-bins bins))
	    (new (copy-power power)))
	(setf (power-data new) (make-array bins) ;New arrays to use
	      (power-ks new) (make-array bins))
	(unless (= (* bins factor) old-bins)
	  (error "New bin count ~D does not go evenly into old count ~D" bins old-bins))
	 (loop for bin below bins
	       do (loop repeat factor for old-index from (* bin factor)
			for k = (aref (power-ks power) old-index)
			sum (aref (power-data power) old-index) into total-data
			when (plusp k)
			  sum k into total-k
			  and count t into k-count
			do (setf (aref (power-data new) bin) (/ total-data factor)
				 (aref (power-ks new) bin) (if (zerop k-count) 0.0
							     (/ total-k k-count)))))
	 new))))
      
(defun plot-power-all (directory &rest keys)
  (let ((title directory))
    (setq directory (merge-pathnames directory batch-root-directory))
    (apply #'plot-power directory :times (dump-times directory) :title title keys)))

(defun plot-power (directory
		   &rest keys 
		   &key
		   time	    ;A single time
		   times	    ;A list of times
		   &allow-other-keys
		   )
  (when time (push time times))	;Convert single time to list
  (apply #'plot-power-1
	 (loop for time in times
	       collect (list directory time (format-float-reasonably time)))
	 keys))

;;Two-point function

(mirror-images
;;Return two-point function <a'(sigma+delta) . a'(sigma)> in an array
;;Slot j stores distance jL/2N.
(defun two-point-function-a (diamond &optional extend-k)
  (let ((l (loop-length diamond))
	(power (loop-power-a diamond extend-k)))	;Array of P = 2k|A(k)|^2
    (two-point-function-1 power l (loop-total-a diamond))))
)
 
;;Compute two-point-function from power and total a' or b'
(defun two-point-function-1 (power l total-a)
  (let ((data (make-array (length power) :element-type 'double-float))) ;Working array
    (loop for n from 1 below (length power)
	  do (setf (aref data n) (* (aref power n) (/ l 4 pi n)))) ;Divide by 2k to get |A^2|
    ;;Do cosine transform. 
    ;;The two-point function at distance x is (2pi/L) sum_{n=-infinity}^infinity |A(k_n)|^2 exp(i k_n x),
    ;;with k_n = 2 pi n/L.  We can write it sum_{n=1}^infinity |A(k_n)|^2 2 cos(i k_n x) + A_0^2 and compute it
    ;;with cosft, if we first divide A_0^2 by 2.
    (setf (aref data 0) (/ (expt (3vector-length total-a) 2) ;Zero mode
			   (* 4 pi) l))	;Normalized half of others because there's only one, not +k/-k
    ;;Slot n of the input data stores P(k_n).  So the desired phase is 2 pi n x/L,
    ;;while the phase in the cosine transform is pi k j/N..  Thus slot j of the output data stores distance x_j = jL/2N,
    ;;and distances range from 0 to L/2 as j goes from 0 to N.
    (cosft data)
    ;;Now multiply by prefactor (2pi/L) and an extra factor 2 from 2 cos(..) above.
    (let ((factor (/ (* 4 pi) l)))
      (dotimes (j (length data))
	(setf (aref data j) (* (aref data j) factor))))
    data))


;;Compute two-point function for a single spacing by brute force without smoothing
(defun slow-two-point-function (hats sigmas delta)
  (loop with n = (length hats)
	with l = (aref sigmas (1- n))
	with result = 0.0
	with sigma1 = 0.0		;Two positions separated by delta
	with sigma2 = delta
	with index1 = 0			;Ranges of sigma where these positions are found
	with index2 = (loop for i from 0
			    while (> sigma2 (aref sigmas i))
			    finally (return i))
	for max1 = (aref sigmas index1)	;End of this range
	for max2 = (aref sigmas index2)
	for range1 = (- max1 sigma1)	;Remaining distance to end of range
	for range2 = (- max2 sigma2)
	when (< range1 range2)		;Range 1 will run out first
	do
	(incf result (* (3vector-dot (aref hats index1) (aref hats index2)) range1))
	(incf index1)			;Advance left segment.
	(when (= index1 n)		;Left index reached end of array, so finished
	  (return (/ result l)))
	(setq sigma1 max1		;Advance left position to end of range 1.
	      sigma2 (mod (+ max1 delta) l)) ;Advance right by same amount.  Wrap if needed.
	else do				;Range 2 will run out first
	(incf result (* (3vector-dot (aref hats index1) (aref hats index2)) range2))
	(setq index2 (mod (1+ index2) n)) ;Advance right segment.  Wrap if needed.
	(setq sigma2 (if (zerop index2)	  ;Advance right position to end of range 2
			 0.0		  ;But if wrapping, beginning of wrapped segment = 0
		       max2)
	      sigma1 (mod (- max2 delta) l)) ;Advance left position.  Wrap backward if needed.
	))

(defun bin-two-point-function (sigma-min sigma-max bins &key extend-k (min-length *total-size*))
  (let ((data (make-array bins :element-type 'double-float :initial-element 0d0))
	(counts (make-array bins :element-type 'fixnum :initial-element 0)) ;Number of samples in each bin
	(lengths (make-array bins :element-type 'double-float :initial-element 0d0))) ;Total of sample lengths
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Ignore open strings
     #'(lambda (diamond)				 ;Loops
	 (let ((l (loop-length-and-count-a diamond)))
	   (when (> l min-length)	;Don't do loops that are too small
	     (format t "Length ~$: " l)
	     (mirror-images
	      (format t "~A " :a) (force-output)
	      (bin-two-point-function-1 data counts lengths sigma-min sigma-max bins
					(two-point-function-a diamond extend-k)
					l))))))
    (loop for index below bins		;Average correlators over different strings
	  for count = (aref counts index)
	  when (plusp count)
	  do (setf (aref data index) (/ (aref data index) count)
		   (aref lengths index) (/ (aref lengths index) count))) ;Average length actually used in this bin
    (terpri)
    (values data lengths)))



;;Add one string to bins
(defun bin-two-point-function-1 (data counts lengths sigma-min sigma-max bins l correlation)
  (loop with bin-size = (/ (log (/ sigma-max sigma-min)) bins)
	with n = (length correlation)
	for j from 1 below (length correlation) ;Don't do sigma=0.  Index j corresponds to dsigma = jL/(2N)
	for c = (aref correlation j)
	for sigma = (/ (* j l) 2 n)
	for bin = (floor (/ (log (/ sigma sigma-min)) bin-size))
	when (and (>= bin 0) (< bin bins))
	do
	(incf (aref data bin) c)
	(incf (aref counts bin))
	(incf (aref lengths bin) sigma)))

(defstruct two-point-function
  sigma-min
  sigma-max
  extend-k
  data
  lengths)

(defun start-process-two-point-function (directory time sigma-min sigma-max bins extend-k)
  (let ((data (make-array bins :element-type 'double-float :initial-element 0d0))
	(counts (make-array bins :element-type 'fixnum :initial-element 0)) ;Number of samples in each bin
	(lengths (make-array bins :element-type 'double-float :initial-element 0d0))) ;Total of sample lengths
    #'(lambda (&optional diamond l power total-a) ;Loop to process, or no arg when done
	(if diamond
	    (bin-two-point-function-1 data counts lengths sigma-min sigma-max bins l
				      (two-point-function-1 power l total-a))
	  (with-open-file (stream (two-point-file time directory) :direction :output
				  :if-does-not-exist :create :if-exists :supersede)
	    (loop for index below bins		;Average correlators over different strings
		  for count = (aref counts index)
		  when (plusp count)
		  do (setf (aref data index) (/ (aref data index) count)
					;Average length actually used in this bin
			   (aref lengths index) (/ (aref lengths index) count)))
	    (prin1 (make-two-point-function :sigma-min sigma-min :sigma-max sigma-max :extend-k extend-k
					    :data data :lengths lengths)
		   stream)
	    (pathname stream)
	    )))))

;;Read processed data from file
(defun read-two-point-function (directory time)
  (with-open-file (stream (two-point-file time directory)) (read stream)))
    
;;Accept list of (directory time &optional title delta-unit) and plot functions
;;If no TITLE is given, we use the directory name and time
;;If DELTA-UNIT is given, scales in data are divided by it.
(defun plot-two-point-functions-1 (list &rest gnuplot-keys &key scaling &allow-other-keys)
  (loop with arbitrary-scaling = nil
	for (directory time title delta-unit) in list
	for data = (read-two-point-function (merge-pathnames directory batch-root-directory) time)
	collect (list data time delta-unit (or title (format nil "~A: ~D" directory time))) into data-list
	when delta-unit do (setq arbitrary-scaling t)
	maximize (length (two-point-function-data data)) into max-points
	finally
	(apply #'gnuplot (length data-list) max-points
	       #'(lambda (plot point)
		   (destructuring-bind (data time delta-unit title) (nth plot data-list)
		     (if (eq point :title) title
		       (when (< point (length (two-point-function-data data)))
			 (let* ((c (aref (two-point-function-data data) point))
				(delta (aref (two-point-function-lengths data) point)))
			   (unless (zerop c)			   ;don't plot points with no data
			     (if scaling (setq delta (/ delta time))) ;Convert to l/t for horizontal axis
			     (if delta-unit (setq delta (/ delta delta-unit))) ;Scale by given delta unit
			     (values delta c)))))))
	       :xlabel (cond (arbitrary-scaling "given units") (scaling "delta/t") (t "delta"))
	       :ylabel "correlator"
	       :styles :linespoints
	       gnuplot-keys)))

;;Plot data processed by process-two-point-function
(defun plot-two-point-function (directory
				&rest keys 
				&key
				time	    ;A single time
				times	    ;A list of times
				&allow-other-keys
				)
  (when time (push time times))	;Convert single time to list
  (apply #'plot-two-point-functions-1
	 (loop for time in times
	       collect (list directory time (format-float-reasonably time)))
	 keys))

(defun plot-two-point-function-all (directory &rest keys)
  (let ((title directory))
    (setq directory (merge-pathnames directory batch-root-directory))
    (apply #'plot-two-point-function directory :times (dump-times directory) :title title keys)))

(defun plot-two-point-functions (directories &rest keys)
  (apply #'plot-two-point-functions-1
	 (loop for directory in directories
	       nconc (loop for time in (dump-times (merge-pathnames directory batch-root-directory))
		      collect (list directory
				    time
				    (format nil "~A, T=~A" directory (format-float-reasonably time)))))
	 keys))

;;Plot time evolution of correlator at fixed distance.
(defun plot-two-point-evolution (directory delta &rest keys)
  (setq directory (merge-pathnames directory batch-root-directory))
  (loop for time in (dump-times directory)
	for data = (read-two-point-function directory time)
	for lengths = (two-point-function-lengths data)
	for values = (two-point-function-data data)
	for value = (loop for index below (1- (length lengths))
			  when (<= (aref lengths index) delta (aref lengths (1+ index)))
			  return (+ (aref values index) (* (/ (- (aref values (1+ index)) (aref values index))
							      (- (aref lengths (1+ index))) (aref lengths index))
							   (- delta (aref lengths index))))) ;Interpolate
	collect (list time value) into result
	finally (apply #'gnuplot 1 (length result)
		       #'(lambda (plot point)
			   (declare (ignore plot))
			   (if (eq point :title) (format nil "d = ~$" delta)
			     (values-list (nth point result))))
		       keys)))
		


;;Create Lorentz-transformation matrix.  The time coordinate is element 3.  V is a 3-vector velocity.
;;(dotmv (Lorentz-transform-matrix v) x) will give translate the 4-vector x to a coordinate system moving in the v
;;direction.
(defun Lorentz-transform-matrix (v)
  (let* ((transform (make-array '(4 4) :element-type 'double-float :initial-element 0.0))
	 (v2 (expt (3vector-length v) 2))
	 (gamma (sqrt (/ 1.0 (- 1.0 v2)))))
    (dotimes (i 4)
      (dotimes (j 4)
	(setf (aref transform i j)
	      (if (= i 3)
		  (if (= j 3)		;time-time component
		      gamma
		    (* -1 gamma (aref v j))) ;time-space component
		(if (= j 3)
		    (* -1 gamma (aref v i)) ;space-time component
		  (+ (if (= i j) 1.0 0.0)
		     (if (plusp v2)	;In case v=0, just return identity without division by zero error
			 (/ (* (aref v i) (aref v j) (- gamma 1)) v2)
		       0.0)))))))
    transform))

(mirror-images
;;Give \int a' dsigma = a(L)-a(0).  Since a'_t = 1, the time component is L
;;We work along the future edge of the existing world sheet.
(defun loop-total-a (start-diamond)
  (setq start-diamond (diamond-next-a start-diamond)) ;Advance to diamond that starts a new A
  (loop with total = (make-zero-4vector)
	with diamond = start-diamond
	for a = (diamond-new-a-wrap diamond)   ;This is delta-x = delta-a/2
	do (4vector-incf total a)	       ;Add up all a' vectors
	do (setq diamond (diamond-next-a diamond))
	until (eq diamond start-diamond)
	finally (return (4vector-scale total 2)))) ;cancel factor 1/2 above

;;Give (1/L)\int a dsigma, where a = a(0) + \int a' dsigma
(defun loop-com-a (start-diamond)
  ;;Skip to diamond that starts new A and new B both, i.e., it is the futuremost diamond with this A and B
  (nomirror				;Give consistent results when this is called for A and B
   (loop until (and (diamond-se start-diamond) (diamond-sw start-diamond))
	 do (setq start-diamond (diamond-e start-diamond))))
  (loop	with integral = (make-zero-3vector)	      ;Accumulate \int a dsigma
	;;We define a(0) = b(0) = this x.  Time component measures length
	with a = (3to4vector (diamond-end start-diamond) 0.0)
	with diamond = start-diamond
	do (setq diamond (diamond-next-a diamond)) ;Since we started with diamond-end, first step is in next diamond
	do (let* ((dx (diamond-new-a-wrap diamond))   ;This is delta-x = delta-a/2
		  (da (4vector-scale dx 2)))	      ;actual change in a 
	     (3vector-incf integral (4vector-scale (4vector+ a dx)  ;center of segment of a
						   (4vector-t da))) ;weight
	     (4vector-incf a da))				    ;Go to next a
	until (eq diamond start-diamond)			    ;Exit when start-diamond done
	finally
	(unless (3vector= (diamond-end start-diamond) a fudge-global-coordinates)
	  (warn "~S failed to close by ~S.  Maybe your loop is not in the rest frame"
		:a (3vector-distance (diamond-end start-diamond) a)))
	;;Total weight is length of loop, which is time component of a
	(return (3vector-scale integral (/ 1 (4vector-t a))))))
)

;;Find center of mass velocity of loop (1/L) \int x-dot dsigma = (1/(2L)) \int (b'+a') dsigma
;;The distance from the start to the end of the loop is \int x' dsigma = \int (b'-a') dsigma
;;For a closed, non-wrapping loop, this is zero, so \int a' dsigma = \int b' dsigma.  But to handle
;;wrapping we keep both terms.
;;We work along the future edge of the existing world sheet.
(defun loop-velocity (start-diamond)
  (mirror-image-let* ((total-a (loop-total-a start-diamond))
		      (velocity-a (3vector-scale total-a (/ 1 (4vector-t total-a))))) ;1/L \int a' dsigma
    (3vector-scale (3vector+ velocity-a velocity-b) 0.5)))

;;Loop velocity from hats and dsigmas
(defun loop-velocity-hats (a-hats a-dsigmas b-hats b-dsigmas)
  (mirror-image-let* ((total-a (total-hats-dsigmas a-hats a-dsigmas))
		      (velocity-a (3vector-scale total-a (/ 1 (total-sigma a-dsigmas))))) ;1/L \int a' dsigma
    (3vector-scale (3vector+ velocity-a velocity-b) 0.5)))

;;3vector center of mass of loop.  Not meaningful unless in rest frame.
(defun loop-com (start-diamond)
  (3vector-scale (3vector+ (loop-com-a start-diamond) (loop-com-b start-diamond)) 0.5))

;;Transform the coordinates of the diamonds of a loop by calling the function TRANSFORM on each
;;The global positions in the kink-created data are just copied.  The tag is lost.
;;EQ coordinates get EQ transforms.
(defun transform-loop-coordinates (start-diamond transform)
  (let ((table (make-hash-table :test #'eq)))
    (flet ((do-transform (x)
	     (or (gethash x table)
		 (setf (gethash x table) (funcall transform x)))))
      (loop for previous = nil then diamond
	    for diamond = start-diamond then (diamond-e diamond)
	    for finishing = (and (eq diamond start-diamond) previous) ;diamond is start-diamond for the second time
	    for previous-new = nil then new
	    ;;The without-compiler-notes in the next line is to suppress an unused-code warning.  Somehow some copy
	    ;;of this iteration code must appear a context where finishing is known to be false, but I don't know how
	    ;;that happens
	    for new = (if finishing (without-compiler-notes first-new) ;Already have boosted version of first diamond
			(make-diamond :start (do-transform (diamond-start diamond)) ;Diamond in new coordinates
				      :left (do-transform (diamond-left diamond))
				      :right (do-transform (diamond-right diamond))
				      :end (do-transform (diamond-end diamond))
				      :a-kink-created (diamond-a-kink-created diamond)
				      :b-kink-created (diamond-b-kink-created diamond)))
	    for first-new = new then first-new ;Keep first new diamond for later
	    when previous		       ;Except first time
	    do (cond ((eq (diamond-ne previous) diamond)
		      (setf (diamond-ne previous-new) new ;Match link direction
			    (diamond-sw new) previous-new))
		     ((eq (diamond-se previous) diamond)
		      (setf (diamond-se previous-new) new
			    (diamond-nw new) previous-new))
		     (t (error "Couldn't understand diamond linkage")))
	    when finishing		;Exit after setting up last linkage
	    return first-new))))

;;Transform coordinates of the points of the loop by going to a frame moving with velocity V
;;We return the transformed first diamond.  Otherwise the diamonds do not go anywhere.
(defun Lorentz-transform-loop (start-diamond v)
  (let ((m (Lorentz-transform-matrix v)))
    (transform-loop-coordinates
     start-diamond 
    #'(lambda (x) (dotmv44 m x)))))

;;Return (one diamond of) the loop in its rest frame
(defun rest-frame-loop (diamond)
  (Lorentz-transform-loop diamond (loop-velocity diamond)))

;;Boost hats and dsigmas to rest frame
(defun rest-frame-hats-dsigmas (a-hats a-dsigmas b-hats b-dsigmas)
  (let ((m (lorentz-transform-matrix (loop-velocity-hats a-hats a-dsigmas b-hats b-dsigmas))))
    (mirror-image-let ((new-a-hats (make-array (length a-hats)))
		       (new-a-dsigmas (make-array (length a-hats) :element-type 'double-float)))
      (mirror-images
       (loop for index below (length a-hats)
	     for delta-a = (4vector-scale (3to4vector (aref a-hats index) 1.0) (aref a-dsigmas index))
	     for new-delta-a = (dotmv44 m delta-a)
	     do (setf (aref new-a-hats index) (3vector-normalize (4to3vector new-delta-a))
		      (aref new-a-dsigmas index) (4vector-t new-delta-a))))
      (values new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas))))


;;Put the loop centered at the spatial origin
(defun center-loop (diamond)
  (let ((com (3to4vector (loop-com diamond) 0.0)))
    (transform-loop-coordinates
     diamond #'(lambda (x) (4vector- x com)))))

;;Standardize loop positions relative to the first diamond, so loop does not wrap identification surfaces
(defun standardize-loop (diamond)
  (transform-loop-coordinates
   diamond #'(lambda (x) (standardize-position x (diamond-start diamond)))))

;;Loops graphing

(defmacro read-single-double (stream)
  `(double-float (read-single ,stream)))

(defvar *loops-file-double* nil)	;If set, read double-floats

(declaim (inline read-loop-datum))
(defun read-loop-datum (stream)
  (if *loops-file-double* (read-double stream)
    (read-single-double stream)))

(defun loops-file-element-type ()
  (list 'unsigned-byte (if *loops-file-double* 64 32)))
  
;; print out the first n loops in a file
(defun loop-raw-data (directory n)
  (with-open-file (stream (loop-spectrum-file directory 0)
			  :element-type (loops-file-element-type))
    (loop for i below n
	  for length = (read-loop-datum stream)
	  while length
	  for time = (read-loop-datum stream)
	  with x and y and z
	  when *log-loop-positions*
	  do (setq x (read-loop-datum stream) y (read-loop-datum stream) z (read-loop-datum stream))
	  do (format t "~&(length , time) = ( ~D , ~D )" length time)
	  when *log-loop-positions*
	  do (format t "pos = ( ~D , ~D, ~D )" x y z))))

(declaim (inline split-by-position))	;Avoid float consing

;;Return number from 0 to errorbar-splits based on position of point in already-setup geometry
(defun split-by-position (x y z)
  (declare (optimize speed))
  (let ((ijkl (xtoi (make-4vector x y z 0.0)))) ;Convert to coordinates with split-factor = errorbar-split
    ;;Now we want the number of the level-0 job that contains this spatial position
    (let ((result (make-4vector))
	  (remainder (make-4vector))
	  excess)
      (dotimes (i 4)
	(let ((q (fast-ffloor (aref ijkl i)))) ;Round down
	  (setf (aref result i) q)	;New candidate ijkl
	  (setf (aref remainder i) (- (aref ijkl i) q))))
      ;;If sum is zero, we're done.  Otherwise, rounding got us below level 0
      (setq excess (loop with sum fixnum = 0
			 for value double-float across result
			 do (locally (declare (optimize (safety 0))) ;Don't check that new value still fixnum
			      (incf sum (- (fixnum-round value))))
			 finally (return sum)))
      (dotimes (i excess)	;Must add this many
	(loop with least = 0		;Start with first slot
	      for index from 1 below 4		;Try others
	      ;;Look for largest remainder to increment.
	      when (> (aref remainder index) (aref remainder least)) ;This one better
	      do (setq least index)
	      finally (incf (aref result least)) ;Increment this slot
	      (decf (aref remainder least) 1.0))) ;Adjust remainder so we don't find it again
      (prog1 (number-of-job result)
	(deallocate 4vectors ijkl result remainder)))))


;;Information about prebin files in this run
(defstruct prebin-info
  era
  scaling				;Normal case if set: bin x = l/t, both comoving.  If NIL, bin l
  physical				;scaling = NIL, physical = T: bin physical l.  If scaling set, does not matter
  (p-bins 1 :type fixnum)	       
  (p-min 0.0 :type double-float)	
  (p-max 0.0 :type double-float)
  (x-bins 0 :type fixnum)
  (t-bins 0 :type fixnum)
  (errorbar-split 1 :type fixnum)
  (x-min 0.0 :type double-float)
  (x-bin-size 0.0 :type double-float)
  ;;Comment here said "this is physical time", but I think it is conformal -- kdo, 11/20/18
  (t-min 0.0 :type double-float)
  (t-bin-size 0.0 :type double-float)
  (total-size 0.0 :type double-float)	;Size of simulation
  )

(defun write-prebin-info (info directory integrated)
  (with-group-write-access
   (with-open-file (stream (prebin-info-file directory :integrated integrated)
			   :direction :output :if-does-not-exist :create :if-exists :supersede)
     (prin1 info stream))))

(defun read-prebin-info (directory integrated)
  (with-open-file (stream (prebin-info-file directory :integrated integrated))
    (let ((info (read stream)))
      (unless (prebin-info-scaling info)(error "Only scaling prebinning is currently implemented"))
      (unless (= (prebin-info-p-bins info) 1) (error "Only files with p-bins = 1 are currently supported.  Sorry, but you must prebin again"))
      info)))

;;Read all loops files for runs and bin data into files.
(defun bin-loops-files (directory
			&key (t-bins 24) (x-bins 120) (p-bins 1) 
			(t-min 50.0) t-max
			(x-min 1d-8) x-max
			(p-min 1d-6) ;Set very low since we're not using it
			(p-max 1d3)
			(scaling t)	;If NIL, use l instead of x, if :mass use x = m/t instead of E/t
			(errorbar-split 1) ;Compute error bars by splitting
			(integrated nil) ;if NIL compute just f(x)dxdt binning, if T compute n(x,t)dx
			submit)		;Submit batch jobs to do work.  # to submit or T
  (unless scaling (error "Only scaling prebinning is currently implemented"))
  (setq directory (merge-pathnames directory batch-root-directory))
  (unless (probe-file directory)
    (error "No such directory ~A" directory))
  (format t "~&Binning loops...") (force-output)
  (let* ((run-info (read-run-info-file directory))
	 (*era* (run-info-era run-info)) ;Used for cosmological data
	 (total-size (double-float (run-info-total-size run-info)))
	 (start (run-info-start-time run-info))
	 (initial-time (run-info-start-time run-info))
	 ;;The following may be wrong if you didn't actually run until the default ending time
	 (t-end (default-simulation-end total-size initial-time nil)) ;FIX 
	 (t-min (or t-min start))	
	 (t-max (cond (t-max t-max)	;last time when a loop can be produced which is relevant
		      (t (+ start (light-crossing-time total-size)))))
	 (x-max (or x-max (/ (loop-length-i t-max t-max t-end) t-max 3))) ;FIX
	 (prebin-info (make-prebin-info :era *era* :physical nil :scaling scaling
					:p-bins p-bins
					:p-min p-min
					:p-max p-max
					:x-bins x-bins :t-bins t-bins :errorbar-split errorbar-split
					:x-min x-min :x-bin-size (/ (log (/ x-max x-min)) x-bins)
					:t-min t-min 
					:t-bin-size (/ (log (/ t-max t-min)) t-bins)
					:total-size total-size)))
    (write-prebin-info prebin-info directory integrated)
    (let ((end (loop for worker from 0						;Find how many files there are
		     unless (probe-file (worker-subdirectory directory worker))	;don't care if worker has a file or not
		     return worker)))
      (cond (submit			;Do with batch job
	     (when (eq submit t)
	       (setq submit (if integrated 100 10))) ;default job number
	     (load "load")			    ;Must recompile or we get races
	     (format t "Deleting old success files...") (force-output)
	     (let ((old-files (directory (prebin-success-files directory))))
	       (when old-files
		 (dolist (file old-files) (delete-file file))
		 (format t "done.~%")))
	     (loop for job below submit
		   do (do-submit directory
				 (format nil "prebin-~D" job) (format nil "prebin-~D-" job) batch-flags
				 `(prebin-file-job
				   ,job
				   ,integrated
				   ,@(loop for worker from job by submit below end ;0, 10, 20 ...
					   collect worker))
				 :load-file (format nil "~A/load" (sb-unix:posix-getcwd))))
	     (format t "Waiting for jobs...") (force-output)
	     (loop for job below submit
		   do (loop until (probe-file (prebin-success-file directory job)) do (sleep 1))
		   do (format t "~D " job) (force-output)))
	    (t				;Do right now
	     (dotimes (worker end)
	       (format t "~D " worker) (force-output)
	       (prebin-file directory worker prebin-info integrated))))
      (combine-prebin-files directory :integrated integrated))))

;;Run in batch job in top-level directory
(defun prebin-file-job (job integrated &rest workers)
  (let ((info (read-prebin-info "." integrated)))
    (format t "Prebin job ~D~%" job)
    (dolist (worker workers)
      (let ((file (loop-spectrum-file "." worker)))
	(format t "Processing file ~A~%" file)
	(prebin-file "." worker info integrated)
	)))
  (with-open-file (stream (prebin-success-file nil job) :direction :output))) ;Create empty file to say success

;;Success file location.  If directory is NIL, use current directory
(defun prebin-success-file (directory job)
  (format nil "~@[~A/~]prebin-~D-success" directory job))
(defun prebin-success-files (directory)
  (format nil "~A/prebin-*-success" directory))

;;Read one loops.dat file, write one loops.prebin file, consisting of an array of
;; p-bins*x-bins*t-bins*errorbar-split^3 single floats.  p-bins varies fastest
(defun prebin-file (directory worker prebin-info integrated)
  (write-prebin-data directory
		     worker
		     (bin-loops-file (loop-spectrum-file directory worker) prebin-info integrated)
		     prebin-info
		     integrated))

(defun write-prebin-data (directory worker data prebin-info integrated)
  (with-group-write-access
   (with-open-file (stream (prebin-data-file directory :integrated integrated :worker worker)
			   :direction :output :element-type '(unsigned-byte 32) :if-does-not-exist :create
			   :if-exists :supersede)
     (let* ((errorbar-split (prebin-info-errorbar-split prebin-info))
	    (p-bins (prebin-info-p-bins prebin-info))
	    (x-bins (prebin-info-x-bins prebin-info))
	    (t-bins (+ (if integrated 1 0) (prebin-info-t-bins prebin-info))))
       (dotimes (split (expt errorbar-split 3))
	 (dotimes (t-bin t-bins)
	   (dotimes (x-bin x-bins)
	     (dotimes (p-bin p-bins)
	       (write-single (aref data split t-bin x-bin p-bin) stream)))))))))

(defvar *offset-loop-times* 0.0)	;Add this to conformal times read from file

;;Read one loops.dat file, return array DATA indexed by split, t-bins, x-bins, p-bins
;;Times are conformal.
;;X is l/t or m/t if :mass in PREBIN-INFO
;;Quantity added to bin is x t^3/(Vc x-bin-size p-bin-size t-bin-size), 
;; i.e. we estimate t^3 x f dt dx dp where f is comoving loop production function;
;;See notes.text
(defun bin-loops-file (file prebin-info integrated)
  (let* ((*era* (prebin-info-era prebin-info))
	 (p-bins (prebin-info-p-bins prebin-info))
	 (p-min (prebin-info-p-min prebin-info))
	 (p-max (prebin-info-p-max prebin-info))
	 (p-bin-size (/ (real-log (/ p-max p-min)) p-bins))
	 (x-bins (prebin-info-x-bins prebin-info))
	 (t-bins (prebin-info-t-bins prebin-info))
	 (errorbar-split (prebin-info-errorbar-split prebin-info))
	 (x-min (prebin-info-x-min prebin-info))
	 (x-bin-size (prebin-info-x-bin-size prebin-info))
	 (t-min (prebin-info-t-min prebin-info))
	 (t-bin-size (prebin-info-t-bin-size prebin-info))
	 (total-size (prebin-info-total-size prebin-info))
	 (comoving-volume (/ (expt total-size 3) (sqrt 2.0)))
	 (scaling (prebin-info-scaling prebin-info))
	 (physical (prebin-info-physical prebin-info))
	 (data (make-array (list (expt errorbar-split 3) 
				 (if integrated (1+ t-bins) t-bins) 
				 x-bins p-bins) :element-type 'double-float :initial-element 0.0)))
    (declare (type (simple-array double-float 4) data))
    (when (> errorbar-split 1)		;Splitting?
      (unless *log-loop-positions*
	(error "Can't compute splits unless ~S is set" '*log-loop-positions*))
      (setup-geometry :total-size total-size :split-factor errorbar-split)) ;Use geometry system to split positions
    (with-open-file (stream file :element-type (loops-file-element-type) :if-does-not-exist nil) ;ignore those never created
      (unless stream (return-from bin-loops-file data)) ;and just return nothing
;;      (declare (optimize speed))		;Declaration here so that type verification on entry does not cause notes
      (handler-case			;Exit on EOF.  Using NIL return spoils type optimizations.
          (loop for x-energy double-float = (read-loop-datum stream) ;scaling length (energy) in units of horizon distance
		for p-momentum double-float = (read-loop-datum stream) ;dimensionless momentum p = v gamma
		for time double-float = (+ (read-loop-datum stream) *offset-loop-times*)
		for l-energy = (* x-energy time)
		for physical-l = (* l-energy (relative-scale-factor time))
		for x-mass = (/ x-energy (sqrt (1+ (expt p-momentum 2.0))))
		for this double-float = (cond ((eq scaling :mass) x-mass)
					      (scaling x-energy)
					      (physical physical-l)
					      (t l-energy)) ;not scaling and not physical means plot vs. comoving l
		with split = 0		;Default split if not splitting
		when (zerop p-momentum) ;in the bh case it is always zero
		do (setf p-momentum p-min) ;set it to minimum value to not create problems
		when *log-loop-positions*
		do (let ((x (read-loop-datum stream)) ;Read position
			 (y (read-loop-datum stream))
			 (z (read-loop-datum stream)))
		     (when (> errorbar-split 1) ;Splitting?
		       (setq split (split-by-position x y z))))
		when *log-loop-velocities*
		do (read-loop-datum stream) (read-loop-datum stream) (read-loop-datum stream) ;Skip velocity
;		if (> x-energy largest-x) do (setf largest-x x-energy) and do (setf largest-x-t time)
;		if (< p-momentum smallest-p) do (setf smallest-p p-momentum) and do (setf smallest-p-t time)
 		if (minusp x-energy) do (warn "ignoring negative energy ~F" x-energy)
		else
		when (or integrated	 ;Decide whether to keep loop.  If integrating, always
			 (> time t-min))		   ;Not integrating: keep if after t-min
		do (let ((t-bin (fixnum-floor (/ (real-log (/ time t-min)) t-bin-size)))
			 (x-bin (max 0
				     (min (1- x-bins)	
				  (fixnum-floor (/ (real-log (/ this x-min)) x-bin-size)))))
			 (p-bin (max 0                ;don't ignore any loops, overfill first bin instead
				     (min (1- p-bins) ;don't ignore any loops, overfill last bin instead
					  (fixnum-floor (/ (real-log (/ p-momentum
									p-min)) 
							   p-bin-size)))))) ;bin log pf
		     (when (< t-bin t-bins) ;Below maximum time
		       (cond (integrated ;here, t-bin could be negative
			      (loop for t-slice from (max 0 (1+ t-bin)) to t-bins ;scan over all later time slices (up to t-bins + 1 of them)
				    for t-prime = (* t-min (exp (* t-slice t-bin-size))) ;time of slice = future edge of t-bin
				    for x-prime = (if (eq scaling :mass) (dynamic-loop-x-mass x-mass time t-prime)  
						    (dynamic-loop-x x-mass p-momentum time t-prime)) ;value of x at time of slice
				    while (<= x-min x-prime) ;give up if loop shrinks out of range
				    for x-bin-prime = (fixnum-floor (/ (real-log (/ x-prime x-min)) x-bin-size))
				    when (< x-bin-prime x-bins)
				    ;; number of loops divided by log interval in x gives x p n(x,p,t) a^3 V_c
				    ;; we want d_h^3 x p n(x,p,t) so extra factor of t^3 / V_c
				    do (incf (aref data split t-slice x-bin-prime p-bin)
					     (/ (expt t-prime 3.0)
						comoving-volume x-bin-size p-bin-size))))	;t^3 x p n(x,p,t) dlogx dlogp
			     (t									;Not integrating
			      (when (< x-bin x-bins)
				(incf (aref data split t-bin x-bin p-bin)
				      ;;If we just take the number of loops and divide by the
				      ;;logarithmic intervals in t, p, and l, we get t l p g(l,p,t) a^3 V_c
				      ;;We want x^2 p f(x,p) = d_h^3 l^2 p g(l,p,t) so we need an extra
				      ;;factor of t^3 x / V_c.  We work in conformal time, so d_h = a t with no 2 or 3.
				      (/ (* this (expt time 3.0))
					 comoving-volume t-bin-size x-bin-size p-bin-size)
				      )))))))
	(end-of-file () nil)))
    data))

(defun trim-loops-files (directory &rest args)
  (let ((end (loop for worker from 0	;Find how many workers there are
		   unless (probe-file (worker-subdirectory directory worker))	;don't care if worker has a file or not
		   return worker)))
  (loop for worker from 0 below end
	for file = (loop-spectrum-file directory worker)
	when (probe-file file)
	do (format t "~&~D: " worker) (force-output)
	(with-simple-restart (skip "Skip this worker")
	    (apply #'trim-loops-file file args)))))

;;Remove early loops from file.  Everything before KEEP-TIME is discarded.  IN-N is the number of floats in each entry
;;and OUT-N of the number to keep.  Using 5 and 2 will delete positions that come after time and size
(defun trim-loops-file (file keep-time in-n out-n
			     &key (max-time 4000.0)) ;Last reasonable time.  Watch out for conformal vs. physical times
  (assert (typep keep-time 'float))
  (assert (<= 2 in-n 10))
  (assert (<= 2 out-n in-n))
  (with-group-write-access
   (let ((save-file (format nil "~A.old" (truename file)))
	 (data (make-array (- in-n 2) :element-type 'single-float))
	 (in-entries 0)
	 (out-entries 0))
     (when (probe-file save-file)
       (error "File ~A already exists" save-file))
     (rename-file file save-file)
     (with-open-file (input save-file  :element-type '(unsigned-byte 32))(- in-n 2)
		     (with-open-file (output file :direction :output :element-type '(unsigned-byte 32) :if-does-not-exist :create)
		       (handler-case	;Exit when file runs out
			   (loop for length of-type single-float = (read-single input)
				 for time of-type single-float = (read-single input)
				 unless (and (< 0 time max-time) (< 0 length (* time 2.0))) ;For some reason there are a few loops w/ x>1
				 do (error "Unreasonable time ~G, length ~G.  Is stride right?" time length)
				 do (incf in-entries)
				 do (loop for i below (- in-n 2)
					  do (setf (aref data i) (read-single input)))
				 unless (< time keep-time)
				 do (write-single length output)
				 (write-single time output)
				 (loop for i below (- out-n 2)
				       do (write-single (aref data i) output))
				 (incf out-entries))
			 (end-of-file () nil))))
     (format t "Copied ~D of ~D entries~%" out-entries in-entries))))

;; returns the scaling energy x of a loop of scaling energy x-i & momentum p-i at conformal time t-i measured at later conformal time time
(defun dynamic-loop-x (xm-i p-i t-i time)
  (let* ((p (* p-i (scale-factor-ratio t-i time)))
	 (gamma (sqrt (1+ (expt p 2.0))))
	 (xm (dynamic-loop-x-mass xm-i t-i time)))
    (* xm gamma)))

;; returns future scaling mass of loop of initial mass xm-i at time t-i
(defun dynamic-loop-x-mass (xm-i t-i time)
  (let ((ar (scale-factor-ratio t-i time))
	(tr (/ t-i time)))
    (* xm-i ar tr)))

;;Combine all prebin files under this directory into one
(defun combine-prebin-files-1 (directory info integrated)
  (let* ((*era* (prebin-info-era info))
	 (p-bins (prebin-info-p-bins info))
	 (x-bins (prebin-info-x-bins info))
	 (t-bins (+ (if integrated 1 0) (prebin-info-t-bins info))) ;One more "bin" if integrated
	 (errorbar-split (prebin-info-errorbar-split info))
	 (result (make-array (list (expt errorbar-split 3) t-bins x-bins p-bins)
			     :element-type 'double-float :initial-element 0.0)))
    (declare (fixnum p-bins x-bins t-bins errorbar-split)
	     (type (simple-array double-float 4) result)
;;	     (optimize (speed 3))
	     )
    (format t "~&Combining prebin files...") (force-output)
    (loop for worker from 0
	  as file = (prebin-data-file directory :worker worker :integrated integrated)
	  while (probe-file file)	;Stop when file doesn't exist.  Even never-started workers should have one.
	  do (format t "~D " worker) (force-output)
	  do (with-open-file (stream file :element-type '(unsigned-byte 32))
	       (dotimes (split (expt errorbar-split 3))
		 (dotimes (t-bin t-bins)
		   (dotimes (x-bin x-bins)
		     (dotimes (p-bin p-bins)
		       (incf (aref result split t-bin x-bin p-bin) (read-single-double stream))))))))
    (write-prebin-data directory nil result info integrated)))
   

;;Combines all worker prebin files into a single one.
(defun combine-prebin-files (directory &key integrated)
  (setq directory (merge-pathnames directory batch-root-directory))
  (combine-prebin-files-1 directory
			  (read-prebin-info directory integrated)
			  integrated))

(defmacro check-unimplemented-keywords (rest &rest features)
  (let ((keywords (mapcar #'(lambda (symbol) (intern (symbol-name symbol) "KEYWORD")) features)))
    `(loop for keyword in ,rest by #'cddr	;Supplied keywords
	   when (member keyword ',keywords)
	   collect keyword into trouble
	   finally (when trouble (error "Unimplemented features ~S" trouble)))))

;;Graph loop production function from prebinned data
;;Many options to this code have been deleted, including all plotting of momenta.
;;The older code is in old.lisp
(defun loops-graph (directories &rest gnuplot-keys
		    &key (t-bins 12) (x-bins 40) (x-min 1e-8) (x-max 1.0)
		    (integrated nil)	;Plot n(x,t) instead of f(x)
		    (scaling t)		;If NIL, use l instead of t.
					; :physical means l_phys/t_phys instead of l_com/t_conf = l/d_h
		    errorbars		;Compute error bars using multiple runs
		    title styles
		    loglog
		    xlabel
		    x-power	;plot f(x) or n(x) times x^x-power.  The default is 2 for f, 1 for n.
		    (bh nil) 
		    &allow-other-keys)
  ;;Error when user tries something that we removed.  Must check because otherwise keywords will go to gnuplot 
  ;;which will ignore them.
  (check-unimplemented-keywords gnuplot-keys 
				normalize x-shift subtract-coef
				t-min t-max p-min p-max 3d-p-bins 2d-p-bins p-power errorbar-split
				physical x-spread early-loops split-power)
  (unless x-power
    (setq x-power (if integrated 1 2)))
  (unless (listp directories) 
    (setq directories (list directories)))
  (setq directories (loop for directory in directories
			  collect (merge-pathnames directory batch-root-directory))) 
  (let* ((info (read-prebin-info (first directories) integrated))
	 (*era* (prebin-info-era info))
	 (pre-x-min (prebin-info-x-min info))
	 (x-prebins (prebin-info-x-bins info))
	 (t-prebins (prebin-info-t-bins info))
	 (x-prebin-size (prebin-info-x-bin-size info))
	 (x-merge (/ x-prebins x-bins)) ;Number of bins to combine into one
	 (x-bin-size (* x-prebin-size x-merge))
	 (t-min (prebin-info-t-min info))
	 (t-prebin-size (prebin-info-t-bin-size info))
	 (t-merge (/ t-prebins t-bins))	;Number of filenames to combine into one graph bin
	 (t-bin-size (* t-prebin-size t-merge))
	 (pre-p-min (prebin-info-p-min info))
	 (pre-p-max (prebin-info-p-max info))
	 (p-prebin-size (log (/ pre-p-max pre-p-min))) ;There's only one bin now, checked in read-prebin-info
	 (errorbar-split (prebin-info-errorbar-split info))
	 (splits (cond (errorbars	;Inter-run variation?
			(unless (= errorbar-split 1)
			  (error "Can't compute error bars using both methods"))
			(length directories))
		       (errorbar-split
			(expt errorbar-split 3)) ;No, use splits if any
		       (t 1)))
	 (data (make-array				   ;Array for storing data
		(list splits				   ;Subsets for errorbars.  1 if none.
		      (if integrated (1+ t-bins) t-bins)   ;Number of time steps to plot
		      x-bins)				   ;Number of bins in x
		:element-type 'double-float :initial-element 0.0)))
    (loop for directory in (cdr directories)
	  unless (equalp info (read-prebin-info directory integrated))
	  do (error "~%All directories must have been prebinned the same way"))
    (unless (integerp x-merge)
      (error "x-bins = ~D, which does not evenly divide the ~D bins in the prebin file" x-bins x-prebins))
    (unless (integerp t-merge)
      (error "t-bins = ~D, which does not evenly divide the ~D bins in the prebin file" t-bins t-prebins))
    ;;Read and combine prebinned data
    (format t "~&Reading prebinned data...") (force-output)
    (loop for directory in directories	;Read all data
	  for directory-number from 0
	  do
	  (with-open-file (stream (if bh
				      (prebin-bh-data-file directory :integrated integrated)
				    (prebin-data-file directory :integrated integrated))
				  :element-type '(unsigned-byte 32))
	    (dotimes (split (expt errorbar-split 3))
	      (dotimes (t-bin (if integrated (1+ t-prebins) t-prebins)) ;Number of time steps in file
		(dotimes (x-bin x-prebins)
		  (let* ((new-t-bin (floor t-bin t-merge))
			 (new-x-bin (floor x-bin x-merge))
			 (pre-x (* pre-x-min (exp (* x-bin x-prebin-size)))) ;x=l/d_h as stored in file
			 (x (if (eq scaling :physical) ;x that we use.  If x = l_phys/t_phys, convert to that
				(* pre-x (horizon-size-factor))
			      pre-x))
			 (sample (if errorbars directory-number split))
			 (integrand (read-single-double stream)))
		    (unless (and integrated (plusp (mod t-bin t-merge))) ;For n(x,t), just ignore other times
		      (let ((increment (* integrand 
					  (/ (expt x x-power) ;What we will use
					     (expt pre-x (if integrated 1 2))) ;divide by what is in file already
					  (if integrated 1.0
					    (/ 1.0 t-merge))		       ;average when combining for f(x),
					  (if (eq scaling :physical)	       ;x=l_phys/t_phys?
					      ;;3 powers to convert # in d_h^3 to # in t^3, one for range of
					      ;;loop size and one for time interval if not integrated
					      (expt (horizon-size-factor) (if integrated -4 -5))
					    1.0)
					  p-prebin-size)		       ;Undo normalization from binning
				       ))
			(incf (aref data sample new-t-bin new-x-bin) ;Accumulate average value
			      (/ increment x-merge (length directories)))
			))))))))
    (terpri)
    (multiple-value-bind (values sigma) (compute-loops-sigma data)
      (apply				; x f(x) 2d plot
       #'gnuplot
       (if integrated (1+ t-bins) t-bins)
       x-bins
       #'(lambda (plot point)
	   (if (eq point :title)
	       (if integrated
		   (format nil "~5F" (* t-min (exp (* plot t-bin-size))))
		 (format nil "~5F to ~5F" 
			 (* t-min (exp (* plot t-bin-size)))
			 (* t-min (exp (* (1+ plot) t-bin-size)))))
	     (let* ((x (* pre-x-min (exp (* (+ point 0.5) ;Central x value
					    x-bin-size))))
		    (scaled-f (aref values plot point))
		    (error (and sigma (aref sigma plot point))))
	       (and (plusp scaled-f)	;no point, rather than 0 if no data
		    (< x-min x x-max)
		    (let ((xval (ecase scaling
				  (t x) ;Normal case.  File always has scaling now.
				  ;;Undo scaling by multiplying by central time
				  ((nil) (* x t-min (exp (* (+ plot 0.5) t-bin-size))))
				  ;;x=l/t instead of x=l/d_h.
				  (:physical (* x (horizon-size-factor))))))
		      (if error (values xval scaled-f error)
			(values xval scaled-f)))))))
       :styles (or styles (if (> splits 1) :errorlines :linespoints))
       :title (or title
		  (format nil "x^{~A} ~:[f~;n~](x)" x-power integrated))
       :logscale (if loglog (list :x :y) :x)
       :xlabel (or xlabel
		   (ecase scaling
		     (t "x = l/d_h")
		     ((nil) "l")
		     (:physical "x = l/t")))
       gnuplot-keys))))
	   

;;Accept array indexed by split, t, x and return value array and standard-deviation array indexed by t, x
;;If no split, second arg NIL.
(defun compute-loops-sigma (data)
  (destructuring-bind (splits t-bins x-bins) (array-dimensions data)
    (let ((values (make-array (list t-bins x-bins)))
	  (sigma (and (> splits 1) (make-array (list t-bins x-bins)))))
      (dotimes (t-bin t-bins)
	(dotimes (x-bin x-bins)
	  (loop for split below splits
		for this = (aref data split t-bin x-bin)
		sum this into total
		sum (expt this 2) into total-square
		finally
		(setf (aref values t-bin x-bin) total)
		(when (> splits 1)
		  ;;Let Xi = value in this bin in ith split, X = Sum Xi.
		  ;;The mean of the Xi is Xm = X/N.
		  ;;The variance of the Xi can be estimated by Sum(Xi - Xm)^2/(N-1)
		  ;; = (Sum(Xi^2) - N Xm^2)/(N-1) = (Sum(Xi^2) - X^2/N)/(N-1).
		  ;;The expected variance of X is then N times larger, so it is
		  ;;(N Sum(Xi^2) - X^2)/(N-1).  The -1 is Bessel's correction.
		  (setf (aref sigma t-bin x-bin) (sqrt (/ (- (* total-square splits) 
							     (expt total 2)) 
							  (1- splits)  ;Bessel's correction
							  )))))))
      (values values sigma))))

(defun compute-mean-sigma (data)
  (let ((n (if (listp data) (length data) (array-dimension data 0)))
	mean 
	sigma)
    (loop for i below n
	  for this = (if (listp data) (nth i data) (aref data i))
	  sum this into total
	  sum (expt this 2) into total-square
	  finally
	  (setf mean (/ total n))
	  (setf sigma (sqrt (/ (- (* total-square n)
				  (expt total 2))
			       n (if (< 1 n) (1- n) n))))) ;-1 is Bessel's correction
    (values mean sigma)))


;; Returns a list of the form '("radiation/200.1" ... "radiation/200.5") from the input "radiation/200" 1 5 or "radiation/200" '(1 2 3 4 5)
(defun directory-sequence (base a &optional b)
  (if a
      (let ((chosen (if (listp a) a
		      (loop for n from a to (if b b a)
			    collect n))))
	(loop for i in chosen
	      collect (format nil "~A.~D" base i)))
    (format nil "~A" base)))

;;Draw line demonstrating a power law in gnuplot.  The command goes in :prelude with a new line after
(defun gnuplot-power-law (from to coefficient power)
  (format nil "set arrow from ~A,~A to ~A,~A nohead"
	  from (* coefficient (expt from power))
	  to (* coefficient (expt to power))))

(defun gnuplot-power-law-pivot (from to power pivot-point pivot-value)
  (gnuplot-power-law from to (/  pivot-value (expt pivot-point power)) power))

(defvar *loop-spectrum*)

;;Info about a loop
(defstruct loop-data
  (length 0.0 :type double-float)	;physical length of loop
  (time 0.0 :type double-float)		;conformal time of formation
  (position zero-3vector :type 3vector)
  (velocity zero-3vector :type 3vector))

;;Read loop files.  Return array of loop-data structures.  If RANGE given, include only those created in that range.
(defun read-loop-positions (directory start end &key max-worker velocities range)
  (let ((data (make-array 1000 :adjustable t :fill-pointer 0)))
    (loop for worker from 0
	  while (or (not max-worker) (< worker max-worker))
	  as file = (loop-spectrum-file directory worker)
	  while (probe-file file)	;Stop when worker output file doesn't exist
	  do (format t "~D " worker) (force-output)
	  do (with-open-file (stream file :element-type (loops-file-element-type))
	       (flet ((read-3vector ()
			(make-3vector (read-loop-datum stream) (read-loop-datum stream) (read-loop-datum stream))))
		 (handler-case		;Exit when file runs out
		     (loop for xi of-type double-float = (read-loop-datum stream)
			   for p-i of-type double-float = (read-loop-datum stream)
			   for ti double-float = (read-loop-datum stream)
			   for position = (read-3vector)
			   for velocity = (and velocities (read-3vector))
			   when (and (<= start ti) (< ti end)
				     (or (null range)
					 (point-in-box position (first range) (second range) 3)))
			   do (vector-push-extend
			       (make-loop-data :length (* xi ti (relative-scale-factor ti)) :time ti
					       :position position :velocity (or velocity zero-3vector))
			       data))
		   (end-of-file () nil)))))
    data))

;; Convert binary file "loops.dat" in given directory from double-floats to single-floats
(defun convert-double-to-single (directory)
  (setq directory (truename directory))	;Avoid problems in rename
  (let ((file2 (format nil "~A/loops.dat2" directory))
	(file (format nil "~A/loops.dat" directory)))
    (when (probe-file file2)
      (error "File ~A already exists" file2))
    (rename-file file file2)
    (format t "~A " file) (force-output)
    (with-open-file (input-stream file2 :element-type '(unsigned-byte 64))
      (with-group-write-access
       (with-open-file (output-stream file :direction :output :if-does-not-exist :create :element-type '(unsigned-byte 32))
	 (handler-case
	     (loop (write-single (read-double input-stream) output-stream))
	   (end-of-file () nil)))))))

(defun convert-double-to-single-directory (directory)
  (dolist (file (directory (format nil "~A/worker-*/loops.dat" directory)))
    (convert-double-to-single (directory-namestring file))))

;;Convert a run with batch jobs
(defun convert-double-to-single-submit (directory)
  (setq directory (merge-pathnames directory batch-root-directory))
  (loop for worker from 0
	for worker-directory = (worker-subdirectory directory worker)
	while (probe-file worker-directory) ;Stop when no more workers
	;;Run job in worker subdirectory, with output to "convert-output".
	do (do-submit worker-directory (format nil "convert-~D" worker) "convert-" batch-flags
		      '(convert-double-to-single ".") :load-file "/cluster/home/k/o/kolum/strings/parallel/load")))

;;Write time slice for Mathematica
(defun write-plot-data (directory start end step min max &key parallel)
  (loop for time from start by step below end
	for count from 0
	do (if parallel
	       (do-submit nil (format nil "plot-~D" count) (format nil "~A/plot-~D-" directory count)
			  (if (eq parallel t) "normal_public" parallel)
			  `(write-plot-data-1 ,directory ,time ,count (make-range ,min ,max))
			  :load-file "load")
	     (write-plot-data-1 directory time count (make-range min max)))))

(defun write-plot-data-1 (directory time count range)
  (warn "On 12 December 2018, write-plot-data-1 was modified by kdo to account changes to get-plot-data, and not tested")
  (with-open-file (stream (format nil "~A/plot-data-~D.m" directory count)
			  :direction :output :if-exists :supersede :if-does-not-exist :create)
    (read-dumps directory :time time)
    (format stream "time = ~A~&data = {" (format-float-reasonably time))
    (let ((first-string t))
      (dolist (string (get-plot-data :time time :range range :velocities-p t))
	(unless (zerop (length string))
	  (format stream "~:[},~&~;~]{" first-string) ;close previous string if any, open next
	  (setq first-string nil)
	  (loop with first = t		;First element in string
		for vector across string
		if vector
		do (format stream "~:[,~%~;~]{~F,~F,~F}" ;new line except first time, then data
			   first
			   (3vector-x vector) (3vector-y vector) (3vector-z vector))
		(setq first nil)
		else			;Break in string: start new
		do (unless first	;Repeated NIL?  Don't do anything.
		     (format stream "},~%{") (setq first t))))))
    (format stream "}};~%")		;close last string and list of strings, and add semicolon 
    ))


(defun histogram-angles (&rest gnuplot-keys &key (min 1e-10) (bins 10) &allow-other-keys)
  (let ((data (make-array bins :initial-element 0))
	(bin-size (/ (- (log min)) bins))
	(zero 0))
    (dolist (string (get-all-strings))
      (loop for (diamond next) on string
	    while next
	    unless (eq next (diamond-ne diamond)) ;Ignore if it has the same a
	    do (let* ((p1 (diamond-a-wrap diamond))
		      (p2 (diamond-a-wrap next))
		      (cos (/ (3vector-dot p1 p2) (3vector-length p1) (3vector-length p2))))
		 (when (> cos 1.0)
		   (when (> cos (+ 1.0 fudge-coordinates))
		     (error "Unreasonable cosine ~S" cos))
		   (setq cos 1.0))
		 (let ((angle (acos cos)))
		   (if (< angle min) (incf zero)
		     (incf (aref data (floor (log (/ angle pi min)) bin-size))))))))
    (apply #'gnuplot 1 bins
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (unless (eq point :title)
		 (values (* min (exp (* bin-size (+ point 0.5)))) ;Central value
			 (aref data point))))
	   :logscale :x
	   :styles :linespoints
	   gnuplot-keys)
    (values zero data)))

(defun histogram-edges (directory start end step
				  &rest gnuplot-keys
				  &key (min 1e-16) (max diamond-span) (bins 30) &allow-other-keys)
  (let* ((steps (1+ (floor (- (+ end fudge-global-coordinates) start) step)))
	 (data (make-array (list steps bins) :initial-element 0))
	 (bin-size (/ (log (/ max min)) bins)))
    (loop for index from 0 below steps
	  for time = (+ start (* index step))
	  do
      (read-dumps directory :time time)
      (map-all-diamonds
       #'(lambda (diamond)
	   (let ((start (diamond-start diamond)))
	     (mirror-images
	      ;;Standardization is not really right for dumps
	      (let ((length (3vector-distance (standardize-position (diamond-left diamond) start) start)))
		(incf (aref data index (floor (log (/ length min)) bin-size)))))))))
    (apply #'gnuplot steps bins
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (format nil "time ~A" (format-float-reasonably (+ start (* plot step))))
		 (values (* min (exp (* bin-size (+ point 0.5)))) ;Central value
			 (aref data plot point))))
	   :logscale :x
	   gnuplot-keys)))

(defun plot-interstring-distance (directories &rest gnuplot-keys &key (titles t) (styles :lines)
					      (offset-time 0.0) (over-time t) velocities t-rescale
					      &allow-other-keys)
  (unless (listp directories) (setq directories (list directories)))
  (let* ((data (map 'vector #'(lambda (directory)
				(read-interstring-distance directory :offset-time offset-time :over-time over-time
							   :velocities velocities))
		    directories))
	 (sample-values nil))
    (apply #'gnuplot (length data) (loop for table across data maximize (length table))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (and titles (nth plot directories))
		 (let* ((plotdata (aref data plot))
			(tpoints (length plotdata))
			(pointdata (if (< point tpoints) (aref plotdata point)))
			(scale (cond ((eq t-rescale :tf)
				      (first (aref plotdata (1- tpoints))))
				     ((eq t-rescale :ti)
				      (first (aref plotdata 0)))
				     (t
				      1.0)))
			(rescaled-data (and pointdata (second pointdata) ;Have this time and a length for it?
					    (list (/ (first pointdata) scale) (second pointdata)))))
		   (if (= point (round (1- tpoints))) (push (second pointdata) sample-values))
		   (values-list rescaled-data))))
	   :styles styles
	   gnuplot-keys)
    (format t "~%mean and standard deviation of last time value of ~%~D ~%are~%" sample-values)
    (compute-mean-sigma sample-values)))

;;Plot distance to which information has traveled.  Nakib Protik.
(defun plot-information-flow-distance (directories &rest gnuplot-keys &key (titles t) (styles :lines)
					      (offset-time 0.0) t-rescale
					      &allow-other-keys)
  (unless (listp directories) (setq directories (list directories)))
  (let* ((data (map 'vector #'(lambda (directory)
				(read-information-flow-distance directory :offset-time offset-time))
		    directories))
	 (sample-values nil))
    (apply #'gnuplot (length data) (loop for table across data maximize (length table))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (and titles (nth plot directories))
		 (let* ((plotdata (aref data plot))
			(tpoints (length plotdata))
			(pointdata (if (< point tpoints) (aref plotdata point)))
			(scale (cond ((eq t-rescale :tf)
				      (first (aref plotdata (1- tpoints))))
				     ((eq t-rescale :ti)
				      (first (aref plotdata 0)))
				     (t
				      1.0)))
			(rescaled-data (and pointdata (second pointdata) ;Have this time and a length for it?
					    (list (/ (first pointdata) scale) (second pointdata)))))
		   (if (= point (round (1- tpoints))) (push (second pointdata) sample-values))
		   (values-list rescaled-data))))
	   :styles styles
	   gnuplot-keys)
    (format t "~%mean and standard deviation of last time value of ~%~D ~%are~%" sample-values)
    (compute-mean-sigma sample-values)))

;;Read information-flow data file.  Nakib Protik.
(defun read-information-flow-distance (directory &key (offset-time 0.0))
  (setq directory (merge-pathnames directory batch-root-directory))
  (with-open-file (stream (information-flow-distance-file directory) :element-type '(unsigned-byte 64))
    (coerce (loop for time = (handler-case (read-double stream) (end-of-file () nil))
		  while time
		  for information-flow-distance = (read-double stream)
		  when offset-time do (incf time offset-time) ;Adjust sample time
		  collect (list* time (list information-flow-distance))) 'vector)))

(defun plot-average-velocity (directories &rest gnuplot-keys &key (titles t) (styles :lines) &allow-other-keys)
  (unless (listp directories) (setq directories (list directories)))
  (let ((data (map 'vector #'(lambda (directory) (read-interstring-distance directory :velocities t))
		   directories)))
    (apply #'gnuplot (length data) (loop for table across data maximize (length table))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (and titles (nth plot directories))
		 (let ((subset (aref data plot)))
		   (and (< point (length subset))
			(let ((entry (aref subset point)))
			  (values (first entry) (third entry)))))))
	   :title "<v^2>"
	   :styles styles
	   gnuplot-keys)))
	  
;;Read processed interstring distance.  Returns an array of (time distance)
;;If VELOCITIES, returns (time distance average-v^2)
;;If RAW-DATA, give length instead of distance
;;If OFFSET-TIME is set, add it to times read in
;;If OVER-TIME is NIL, don't divide by time.
(defun read-interstring-distance (directory &key (over-time t) (offset-time 0.0) (raw-data nil) (velocities nil))
  (setq directory (merge-pathnames directory batch-root-directory))
  (let* ((info (read-run-info-file directory))
	 (comoving-volume (/ (expt (run-info-total-size info) 3) (sqrt 2.0))))
    (with-open-file (stream (length-file directory) :element-type '(unsigned-byte 64))
      (coerce (loop for time = (handler-case (read-double stream) (end-of-file () nil))
		    while time
		    for length = (read-double stream)
		    for velocity = (and velocities (read-double stream))
		    when offset-time do (incf time offset-time) ;Adjust sample time
		    collect (list* time
				   (if raw-data length
				     (let ((density (/ length comoving-volume))) ;Density at that time
				       (and (plusp density) (plusp time)
					    (/ (expt density -1/2) (if over-time time 1.0)))))
				   (and velocities
					(list (cond (raw-data velocity)
						    ((zerop length) 0.0) ;Don't crash
						    (t (/ velocity length)))))))
	      'vector))))

(defun parse-output-starts (file)
  (with-open-file (stream file)
    (loop for line = (read-line stream nil)
	  while line
	  for rest = (or (line-starts-p line "Initializing...Starting at global time ")
			 (line-starts-p line "Starting at global time ")) ;If GC bewteen
	  when (and rest (read-from-string rest t nil :end (position #\. rest :end 7 :from-end t)))
	  collect it)))

(defvar *read-line-skip-gc-message-next* nil)
(defvar *read-line-skip-gc-message-line-number*) ;Only we know the line number
;;Repeated calls to this function will give successive lines and skip GC messages
(defun read-line-skip-gc-message (stream)
  (let ((line (or *read-line-skip-gc-message-next*
		  (progn (incf *read-line-skip-gc-message-line-number*)
			 (read-line stream nil)))))
    (when line				;Return NIL if nothing more
      (setq *read-line-skip-gc-message-next* nil) ;If this was set, we used it
      (let ((next (read-line stream nil))) ;Get subsequent line
	(incf *read-line-skip-gc-message-line-number*)
	(cond ((line-starts-p next "Next GC in region 2 at") ;GC message.
	       (incf *read-line-skip-gc-message-line-number*)
	       (format nil "~A~A" line (read-line stream))) ;Concatenate next line
	      (t (setq *read-line-skip-gc-message-next* next) ;Save next line for next time
		 line))))))

;;Figure out how long the lengths.dat file should be very reading the output log
;;SIZE is the size of the individual runs.
(defun parse-output-length (file)
  (with-open-file (stream file)
    (setq *read-line-skip-gc-message-next* nil *read-line-skip-gc-message-line-number* 0)
    (loop with start = nil
	  for line = (read-line-skip-gc-message stream)
	  while line
	  finally (when start (warn "End of file in the middle of a run"))
        	  (return (* entries 16)) ;Each entry should be double floats for time and length
	  for start-text = (line-starts-p line "Initializing...Starting at global time ")
	  for new-start = (and start-text
			       (read-from-string start-text t nil :end (position #\. start-text :end 7 :from-end t)))
	  for end-text = (line-starts-p line "Finished with global time ")
	  for new-end = (and end-text (read-from-string end-text))
;;	  for restart = (line-starts-p line "This is SBCL")
	  when new-start
	  do (when start
	       (warn "Two starts in a row in line ~S" *read-line-skip-gc-message-line-number*))
	  (setq start (round-to-job-start new-start))
	  when new-end
	  do (unless start (error "End without a start in line ~S" *read-line-skip-gc-message-line-number*))
	  and sum (1+ (- (floor (round-to-job-start new-end))
			 (ceiling start))) into entries ;number of dumps for this run
	  and do (setq start nil)
;;	  when (and restart entries) do (format t "Lisp restart after ~D entries~%" entries)
	  )))

;;Convert time read for output file into actual job starting to.  Ware kluges!
(defun round-to-job-start (time)
  (if (fudge= time *initial-time* 1e-14)
      time				;Initial time
    (let* ((step (/ (job-end-t *size*) 4)) ;Time between layers
	   (offset (+ *initial-time* (/ step 2)))) ;First job end
      (+ (* step (round (- time offset) step))
	 offset))))

(defun check-output-length (directory worker *size* *initial-time*)
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((length (with-open-file (stream (length-file directory worker)) (file-length stream)))
	(expect (parse-output-length (format nil "~A/output" (worker-subdirectory directory worker)))))
    (unless (= length expect)
      (format t "~&~A: File length is ~D instead of ~D, so ~D bytes were lost~%"
	      (worker-subdirectory directory worker) length expect (- expect length)))))
      
(defun check-output-lengths (directory size initial-time max-worker)
  (loop for worker from 0 below max-worker
	do (format t "~D " worker) (force-output)
	do (check-output-length directory worker size initial-time)))

(defun check-output-lengths-submit (directory size initial-time max-worker)
  (do-submit directory "check-lengths" "check-" batch-flags
	     `(check-output-lengths ,directory ,size ,initial-time ,max-worker)
	     :load-file "/cluster/home/k/o/kolum/strings/parallel/load"))
      

(defun debug-lengths-file (directory worker &key (initial-time 6.0) (size 30.0))
  (with-open-file (stream (length-file directory worker) :element-type '(unsigned-byte 64))
    (loop with print = 0
	  with history = nil
	  with job-duration = (job-end-t size) ;Duration of each job
	  with first-end = (+ initial-time (/ job-duration 8)) ;End of first box
	  with start
	  with current
	  for previous = nil then time
	  for previous-length = nil then length
	  for time = (read-double stream)
	  for length = (read-double stream)
	  unless start do (setq start time) ;first time
	  do (push time current)
	  when (and (not (zerop length)) (fudge= length (round length) 1e-14))
	  do (warn "Integer length ~F at time ~F at position ~D" length time (file-position stream))
	  when (and previous		;If this extends the previous run, just keep going.  Otherwise...
		    (not (fudge= time (1+ previous) 1e-14)))
	  do				;End of run
	  (cond ((zerop time)		;Bad data?
		 (format t "~&Bad time ~S.  Prev. Len ~S. Current list: ~S" time previous-length current)
		 (format t "~&Next: ")
		 (format t "~F " length) ;Length for bad time
		 (dotimes (i 47)
		   (format t "~F " (read-double stream)))
;;		 (return (nreverse history))
		 (setq print 30 time nil start nil))
		(t			;New run.  Record all boxes handled.
		 (when (plusp print)
		   (decf print)
		   (format t "~S -- ~S ~%" start previous))
		 ;;Convert last dump in run into a ending time of region
		 (let ((run-end (+ (* (ceiling (- previous first-end) (/ job-duration 4))
				      (/ job-duration 4))
				   first-end))
		       (run-start (+ (* (floor (- start first-end) (/ job-duration 4)) ;Starting time of first region
					(/ job-duration 4))
				     first-end)))
		   (when (< run-start (- job-duration))
		     (error "Unreasonable run start ~S" run-start))
		   (unless (= previous (floor run-end))
		     (error "Last dump ~D not the last in a region" previous))
		   (cond ((fudge= start initial-time 1e-14) ;Started in middle at initial time
			  (push initial-time history)
			  (setq run-start
				(- run-end
				   (* (floor (- run-end initial-time) job-duration) ;number of jobs after initial
				      job-duration))))
			 (t
			  (unless (= start (ceiling run-start))
			    (error "First dump ~D not the first in a region" previous))))
		   (let ((jobs (/ (- run-end run-start) job-duration))) ;How many jobs in this run
		     (unless (fudge= jobs (round jobs) 1e-10)
		       (error "Not an even number of jobs: ~S" jobs))
		     (loop repeat (round jobs)
			   for this from run-start by job-duration
			   do (push this history))))
		 (setq start time current (list time)) ;Start new run
		 )))))


;;Read output file to determine sizes of layers.
(defun parse-layer-sizes-1 (file data split-factor)
  (with-open-file (stream file)
    (setq *read-line-skip-gc-message-next* nil *read-line-skip-gc-message-line-number* 0)
    (loop with job
	  for line = (read-line-skip-gc-message stream)
	  while line
	  for job-text = (line-starts-p line "********** JOB ")	;Job starting (or completed)
	  for reading-text = (line-starts-p line "Reading /cluster") ;Reading from predecessor
	  when job-text do (setq job (parse-integer job-text :junk-allowed t)) ;Remember job number
	  when reading-text
	  do (let* ((position (search ".dat: " reading-text))
		    (count (parse-integer reading-text :start (+ 6 position) :junk-allowed t)))
	       (incf (aref data (floor job (expt split-factor 3))) count)))))

(defun parse-layer-sizes (directory split-factor)
  (let ((data (make-array (* split-factor 10) ;Should be plenty for number of layers
			  :element-type 'fixnum :initial-element 0)))
    (loop for worker from 0
	  for file = (format nil "~A/output" (worker-subdirectory directory worker))
	  while (probe-file file)
	  do (format t "~D " worker) (force-output)
	  do (parse-layer-sizes-1 file data split-factor))
    (loop for index downfrom (1- (length data))	;Trim off trailing zeros
	  while (zerop (aref data index))
	  finally (return (subseq data 1 index)) ;No diamonds enter first layer
	  )))
  
(defun plot-layer-sizes (directory split-factor &rest gnuplot-keys)
  (let ((data (parse-layer-sizes directory split-factor)))
    (apply #'plot-data-list
	   (loop for index from	1
		 for value across data
		 collect (list (log index) (log value)))
	   gnuplot-keys)))



(defun process-interstring-distance-submit (directory)
  (setq directory (merge-pathnames directory batch-root-directory))
  (do-submit directory "interstring" "process-interstring-" batch-flags
	     `(process-interstring-distance ,directory)
	     :load-file "/cluster/home/k/o/kolum/strings/parallel/load"))

;;Read string length from worker files, combine into one file in top directory
(defun process-interstring-distance (directory &key velocities)
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((times (length-times directory)))
    (unless times
      (error "Lengths were not dumped in ~A" directory))
    (let ((data (make-hash-table))	;time -> total length
	  (velocity-data (and velocities (make-hash-table))))
      (dolist (time times)
	(setf (gethash time data) 0.0)	;Initialize tables
	(when velocities (setf (gethash time velocity-data) 0.0))) 
      (with-group-write-access
       (with-open-file (out-stream (length-file directory) :direction :output :element-type '(unsigned-byte 64)
				   :if-does-not-exist :create :if-exists :supersede)
	 (loop for worker from 0
	       as file = (length-file directory worker)
	       while (probe-file (worker-subdirectory directory worker)) ;Stop if directory was never created
	       do (with-open-file (stream file :element-type '(unsigned-byte 64) :if-does-not-exist nil)
		    (cond (stream
			   (format t "~D " worker) (force-output)
			   (process-interstring-distance-1 stream data times velocity-data))
			  (t (format t "No such file: ~A~%" file))))) ;Some workers might never have started
	 (loop for time in times
	       do (write-double time out-stream)
	       do (write-double (gethash time data) out-stream)
	       when velocities do (write-double (gethash time velocity-data) out-stream)
	       ))))))

;;Read one lengths file and increment lengths in hash table
;;2/1/16: When I changed this to use hash tables, I removed the code for dealing with corrupted files.  If it is
;;necessary to read an old, corrupted file, see old code in version 1.171 or earlier.
;;This code depends on having the time to which the dump was done written in run-info.lisp, and communicated to the
;;workers to do the dumps, in such a way that there's no corruption of lower bits in the floating point representation
;;of the time, so that it remains EQL.
(defun process-interstring-distance-1 (stream data times &optional velocity-data)
  (let ((reading-time nil))		;Flag to detect file ending in the middle of an entry
    (handler-case
         (loop
	  (setq reading-time t)		     ;EOF OK now
	  (let* ((time (read-double stream)) ;Time of report
		 (old (gethash time data)))  ;Look up in table
	    (unless (or old
			;;Extra times are OK if they are significantly after the ending time.  They result when we
			;;are finishing a job that covers the next dump time which is not covered completely.
			(> time (+ (car (last times)) 1e-12)))
	      (error "Found time ~F in file ~A, but that is not a time at which lengths were dumped."
		     time (pathname stream)))
	    (setq reading-time nil)	;EOF would be an error now
	    (let ((length (read-double stream))
		  (velocity (and velocity-data (read-double stream))))
	      (when old
		(setf (gethash time data) (+ old length))
		(when velocity
		  (incf (gethash time velocity-data) velocity))))))
      (end-of-file ()
        (unless reading-time			;Ran out after reading half entry?
	  (error "File ~S ended in the middle of a record" (pathname stream)) ;Corrupted file
	  )))))

;;Plotting of tangent vectors on the unit sphere

;;Mollweide projection.  Returns x and y values.  See http://en.wikipedia.org/wiki/Mollweide_projection
;;Results in range -2<x<2, -1<y<1
(defun Mollweide (x y z)
  (let* ((lambda (atan y x))		;longitude
	 (phi (atan z (sqrt (+ (expt x 2) (expt y 2))))) ;latitude
	 (rhs (* pi (sin phi)))
	 (alpha (* phi 2)))		;Initial guess for alpha = 2 theta
    ;;Solve alpha + sin alpha = pi sin phi
    (loop for count from 1
	  for result = (+ alpha (sin alpha))
	  until (fudge= result rhs 1e-10)
	  do (decf alpha  (/ (- result rhs) (+ 1 (cos alpha))))	;Newton's method
	  when (>= count 100) do (error "Newton's method failed to converge in Mollweide iteration"))
    (values (* (/ 2 pi) lambda (cos (/ alpha 2))) (sin (/ alpha 2)))))

(defun inverse-Mollweide (x y)
  (let* ((alpha (* 2 (asin y)))
	 (lat (asin (/ (+ alpha (sin alpha)) pi)))
	 (long (/ (* pi x) 2 (cos (/ alpha 2)))))
    (values (* (cos lat) (cos long)) (* (cos lat) (sin long)) (sin lat))))

;;List of (X Y Z) axes (unit 3vectors)
(defvar *Mollweide-last-axes* nil)

;;AXES can be
;; :LAST -- do whatever we did last
;;  NIL -- cardinal axes
;; :CENTER -- center values to be plotted, with arbitrary rotation
;;            in this case POINTS should be a sequence of arrays of 3vector, (3vector ...) or NIL
;; a 3vector -- use this as center, arbitrary rotation
;; a list of 2 3vectors -- these are Z and X, so the cut is from Z to -X to -Z.
(defun choose-Mollweide-axes (axes &optional points)
  (let (x y z)
    (cond ((null axes)
	   (setq x (make-3vector 1.0 0.0 0.0)
		 y (make-3vector 0.0 1.0 0.0)
		 z (make-3vector 0.0 0.0 1.0)))
	  ((eq axes :last)
	   (setq x (first *Mollweide-last-axes*) y (second *Mollweide-last-axes*) z (third *Mollweide-last-axes*)))
	  ((eq axes :center)
	   (let ((center zero-3vector))
	     (map nil #'(lambda (vectors)
			  (loop for p across vectors
				when (consp p) do (setq p (car p))
				when p do (setf center (3vector+ center p))))
		  points)
	     (setq x (3vector-normalize center))
	     (multiple-value-setq (y z) (two-perpendicular-vectors x))))
	  ((typep axes '3vector)
	   (setq x (3vector-normalize axes))
	   (multiple-value-setq (y z) (two-perpendicular-vectors x)))
	  ((typep axes 'cons)
	   (setq z (3vector-normalize (first axes))
		 ;;Component perpendicular to z in case user did not supply them perpendicular
		 x (3vector-normalize (3vector- (second axes) (3vector-scale z (3vector-dot (second axes) z))))
		 y (3vector-cross-product z x)))
	  (t (error "Don't understand Mollweide AXES specification ~S" axes)))
    (setq *Mollweide-last-axes* (list x y z))
    (values x y z)))

;;Project the 3vector V using given axes
(defun Mollweide-vector (v x y z)
  (Mollweide (3vector-dot v x) (3vector-dot v y) (3vector-dot v z)))

;;Inverse transformation
(defun inverse-Mollweide-vector (x y)
  (multiple-value-bind (vx vy vz) (inverse-Mollweide x y)
    (3vector+ (3vector-scale (first *Mollweide-last-axes*) vx)
	      (3vector-scale (second *Mollweide-last-axes*) vy)
	      (3vector-scale (third *Mollweide-last-axes*) vz))))

;;Give ellipse bounding Mollweide projection.  Returns array of (x y)
(defun Mollweide-ellipse (&optional (points 200))
  (let ((data (make-array (1+ points))))
    (loop for step to points
	  for angle = (* 2 pi (/ step points))
	  do (setf (aref data step) (list (* 2 (cos angle)) (sin angle))))
    data))

;;See if segment goes through identification surface
(defun Mollweide-segment-wraps-identification-p (x y v1 v2)
  (let ((phi1 (atan (3vector-dot y v1) (3vector-dot x v1))) ;-pi...pi from x axis
	(phi2 (atan (3vector-dot y v2) (3vector-dot x v2))))
    (and (/= (signum phi1) (signum phi2)) ;Different hemispheres?
	 (> (abs (- phi1 phi2)) pi))))	  ;If distance > pi, goes around back

;;Break between vectors when lines wrap identification
(defun break-Mollweide-wrapping (x y vectors)
  (loop with new-vectors = (make-array (* 2 (length vectors)) :adjustable t :fill-pointer 0)
	for index below (length vectors)
	for last = (aref vectors (1- (length vectors))) then this
	for this = (aref vectors index)
	when (Mollweide-segment-wraps-identification-p
	      x y
	      (if (consp last) (car last) last) ;each vector may actually be (vector . extra-data)
	      (if (consp this) (car this) this))
	  do (vector-push-extend nil new-vectors) ;break line
	do (vector-push-extend this new-vectors)
	finally (return new-vectors)))

;;Plot some lists of tangent vectors on the unit sphere
;;DATA is a sequence of (VECTORS TITLE STYLE).  VECTORS is an array of 3vectors or lists (3vector . extra-data)
;;The extra data will be returned as extra values in gnuplot, so you can use it for circle radii, etc.
;;A null element of VECTORS causes a break in the line
;;TITLE and STYLE are passed to gnuplot.
;;Each VECTORS array is plotted separately.  If LOOP is set, we plot line from the last point to the first of each
;;If you specify STYLES, it overrides the styles in the data.  Your vectors come first, then the ellipse.
;;If NO-WRAP, break lines between points that would have gone through the Mollweide projection identification surface
;;If 3D T, plot without projecting.
;;MARK-LOCATIONS is a list of (command 3vector) to mark the specified locations.
;;The 3vector will be projected and the command will go in the prelude followed by " at " and the projected coordinates
(defun plot-tangent-vectors (data &rest gnuplot-keys &key 3d styles (axes :center) loop prelude show-axes
				  (no-wrap t) mark-locations
				  &allow-other-keys)
  (setq data (coerce data 'vector))	;Make into vector for fast processing
  (unless styles
    (setq gnuplot-keys (list* :styles
			      (append (loop for (vectors title style) across data collect style)
				      (and (not 3d) '(:lines))) ;Ellipse
			      gnuplot-keys)))
  (when show-axes
    (setq mark-locations
	  (append mark-locations
		  `(("set label 'X'" ,(make-3vector 1.0 0.0 0.0))
		    ("set label 'Y'" ,(make-3vector 0.0 1.0 0.0))
		    ("set label 'Z'" ,(make-3vector 0.0 0.0 1.0))))))
  (unless 3d
    ;;In case of :style :circle, avoid circles forcing larger ranges
    (setq gnuplot-keys (append gnuplot-keys (list :xrange '(-2.0 2.0) :yrange '(-1.0 1.0)))))
  (let ((ellipse (and (not 3d) (Mollweide-ellipse)))
	x y z)
    (unless 3d
      (multiple-value-setq (x y z) (choose-Mollweide-axes axes (loop for (vectors) across data collect vectors)))
      (when mark-locations
	(setf (getf gnuplot-keys :prelude)
	      (format nil "~@[~A~%~]~:{~A at ~F, ~F~%~}"
		      prelude
		      (loop for (command direction) in mark-locations
			    collect (cons command (multiple-value-list (Mollweide-vector direction x y z)))))))
      (when no-wrap			;Break lines where they would wrap identification
	(setq data
	      (map 'vector
		   #'(lambda (entry)
		       (destructuring-bind (vectors . rest) entry
			 (cons (break-Mollweide-wrapping x y vectors) rest)))
		   data))))
    (apply #'gnuplot (+ (length data) (if 3d 0 1))
	   (max (length ellipse) (loop for (vectors) across data maximize (1+ (length vectors))))
	   #'(lambda (plot point)
	       (if (= plot (length data)) ;Ellipse
		   (unless (eq point :title)
		     (and (< point (length ellipse))
			  (destructuring-bind (sx sy) (aref ellipse point)
			    (values sx sy))))
		 (destructuring-bind (vector title style) (aref data plot) ;Regular point
		   (declare (ignore style))
		   (cond ((eq point :title) title)
			 ((null vector) nil)
			 (t
			  (when (and loop (= point (length vector)))
			    (setq point 0))
			  (let ((p (and (< point (length vector)) (aref vector point)))
				(data nil))
			    (when (consp p)
			      (psetq p (car p)
				     data (cdr p)))
			    (and p
				 (if 3d
				     (values-list (append (3vector-list p) data))
				   (multiple-value-bind (sx sy) (Mollweide-vector p x y z)
				     (values-list (list* sx sy data)))
				   ))))))))
	   gnuplot-keys)))

;;Accepts DATA as in PLOT-TANGENT-VECTORS.  Returns a new sequence whose last element connects
;;last element of each set of vectors in DATA to the first element of the next set.  If LOOP, the last
;;set is connected to the first 
(defun add-connecting-segments (data title style &key loop)
  (let ((connections (make-array 100 :adjustable t :fill-pointer 0)))
    (flet ((connect (from-index to-index) ;Make a connecting dotted line from end of from to start of to.
	     (let* ((from (car (aref data from-index)))
		    (to (car (aref data to-index))))
	     (vector-push-extend (aref from (1- (length from))) connections) ;Last of these
	     (vector-push-extend (aref to 0) connections) ;To first of these
	     (vector-push-extend nil connections)))) ;break line
      (loop for index from 0 below (1- (length data))
	    do (connect index (1+ index)) ;All except last: connect each to next
	    finally (when loop (connect index 0)))) ;If loop, connect last to first, otherwise nowhere
    (concatenate 'simple-vector data (list (list connections title style)))))

#|
;;Is this useful anymore, or has it been subsumed by plot-a?

(defun plot-steps (diamond &rest keys &key count (skip 0) &allow-other-keys)
  (setq diamond (diamond-next-a diamond)) ;Advance to diamond that starts a new A
  (loop repeat skip do (setq diamond (diamond-next-a diamond))) ;Skip this many
  (loop with sigma = 0.0
	for first-time = t then nil
	for this = diamond then (diamond-next-a this)
	for p = (diamond-a-wrap this)
	for dsigma = (* 2 (3vector-length p))
	for p-hat = (3vector-scale p (/ 2 dsigma))
	;;Make list of (p-hat sigma)
	collect (list p-hat sigma) into data
	do (incf sigma dsigma)
	until (and (not first-time) (eq this diamond)) ;exit when we have done the starting diamond for a second time.
	until (and count (zerop (decf count))) ;Limit to given number
	finally
	(apply #'plot-steps-1
	       (make-array 1 :initial-element (coerce data 'vector)) ;Need array of arrays
	       keys)))

;;Do the plotting.  DATA is an array of arrays of (p-hat sigma), or a list two such arrays, one for a and one for b
(defun plot-steps-1 (data &rest args &key (projection t) &allow-other-keys)
  (apply (if projection #'plot-steps-projection #'plot-steps-no-projection) data args))

(defun plot-steps-no-projection (data &rest keys &key 3d &allow-other-keys)
  (cond (3d
	 (apply #'gnuplot (length data) (loop for cluster across data maximize (length cluster))
		#'(lambda (plot point)
		    (unless (eq point :title)
		      (destructuring-bind (p-hat sigma) (aref (aref data plot) point)
			(declare (ignore sigma))
			(values-list (3vector-list p-hat)))))
		:styles :linespoints keys))
	(t (unless (= (length data) 1)
	     (error "xyz plotting of multiple segments is not implemented"))
	   (setq data (aref data 0))
	   (apply #'gnuplot 3 (length data)
		  #'(lambda (plot point)
		      (if (eq point :title) (char "xyz" plot)
			(destructuring-bind (p-hat sigma) (aref data point)
			  (values sigma (3vector-component p-hat plot)))))
		  :styles :steps keys))))

(defun plot-steps-projection (data &rest keys &key loop &allow-other-keys)
  (let ((divider nil))
    (when (listp data)			;Plotting both a and b
      (setq divider (length (first data)) ;Remember what is what
	    data (concatenate 'simple-vector (first data) (second data)))) ;Combine arrays
 	(cond (divider			;A and B separate
	       (connect-range 0 divider)
	       (connect-range divider (length data)))
	      (t			;Just one
	       (connect-range 0 (length data)))))
      (apply #'plot-tangent-vectors
	     (nconc
	      (loop for lt from 1	;Line types
		    for vector across data 
		    collect (list vector (format nil "Length ~D" (length vector))
				  (format nil "lp lt ~D" lt)))
	      (list (list connections nil "lines lt 0"))) ;Connections dotted
	     :prelude (format nil "unset border~%unset xtics~%unset ytics~%unset ztics~%set view 0,0,1.2,1.2~%")
	     keys))))

|#

;;Smoothing

;;Smoothing method dispatch
(defun smooth-data (hats sigmas standard-deviation spacing)
  (Gaussian-smooth-data hats sigmas standard-deviation spacing))
	 
;;Return new HATS, SIGMAS with even spacing and smoothing using Fourier transforms.
;;See Fourier.tex
;;Smoothing is done with a Gaussian with given standard deviation
;;We don't keep the zero-mode, so if the a' values don't integrate to zero, strange things will result.
;;Spacing is the minimum spacing, but if the standard deviation is small we will use smaller spacing also
(defun Gaussian-smooth-data (hats sigmas standard-deviation spacing)
  (let* ((points (length sigmas))	   ;Number of segments
	 (l (aref sigmas (1- points)))	   ;Length of loop
	 (needed (/ l standard-deviation)) ;Number of frequencies needed for accuracy.  See Fourier.tex.
					;Size of FFT: At least twice needed frequencies, rounded up to power of 2
	 (nft (max 8 (expt 2 (1+ (ceiling (log needed 2)))))) ;and never less than 8
	 (xi (make-array points :element-type 'double-float)) ;Point locations for Fourier transform
	 (fi (make-array points :element-type 'double-float)) ;Function values for Fourier transform
	 ;;Number of new points, round up to power of 2, but no less than frequencies to be used.
	 (result-count (max nft (expt 2 (ceiling (log (/ l spacing) 2)))))
	 (result1 (make-array result-count :element-type 'double-float))    ;Place to do inverse transform
	 (new-sigmas (make-array result-count :element-type 'double-float)) ;New sigmas to return
	 (new-hats (make-array result-count)))				    ;Array of 3vectors to return.
    (format t "NFT = ~D, result-count = ~D~%" nft result-count)
    (loop for index below result-count
	  do (setf (aref new-sigmas index) (/ (* l (1+ index)) result-count)) ;New sigmas evenly spaced
	  do (setf (aref new-hats index) (make-3vector)))		      ;Create 3vectors to return
    (dotimes (index points)						      ;Build xi
      (setf (aref xi index) (/ (aref sigmas index) l)))			      ;Scale into 0..1
    (dotimes (coordinate 3)						      ;Smooth cardinal directions separately
      (dotimes (index points)						      ;Build fi for this coordinate
	(setf (aref fi index)						      ;fi = (ai - a(i-1))
	      (- (3vector-component (aref hats (mod (1+ index) points)) coordinate)
		 (3vector-component (aref hats index) coordinate))))
      (let* ((transform (nufft xi fi *nufft-m* nft))) ;Do non-uniform Fourier transform
	;;TRANSFORM is the sum of e^{-i omega_j sigma_{i+1} (f_{i+1} - f_i)}
	;;What we need is really the complex conjugate of this, multiplied by i and divided by omega_j
	;;We scale that by the Fourier transform of the Gaussian, exp(-2pi^2 j^2 s^2/L^2), and transform back.
	;;Frequency j in TRANSFORM is stored in slots 2j, 2j+1.
	;;Real values for frequencies 0 and nft/2 are stored in slots 0 and 1
	;;Frequency 0 should not be present in the NUFFT result because of cancellation in f_{i+1} - f_i.
 	;;Frequency nft/2 should be scaled by an infinitesimal number on the tail of the Gaussian, so just ignore it.
	;;Scale and copy coefficients from TRANSFORM into RESULT1, which we will then FFT in place
	(fill result1 0.0)				      ;Clear unused frequencies
	(loop with expt1 = (* -2 (expt (* pi (/ standard-deviation l)) 2)) ;exponent except for j^2
	      with factor1 = (/ l 2 pi)				       ;prefactor except for 1/j
	      for index from 1 below (/ nft 2)			       ;Number of valid coefficients
	      for index2 = (* 2 index)
	      for factor = (* (/ factor1 index) (exp (* expt1 (expt index 2)))) ;Total scaling = (1/omega_j)exp(...j^2)
	      ;;Exchanging real and complex parts implements taking complex conjugate and multiplying by i
	      do (setf (aref result1 index2) (* (aref transform (1+ index2)) factor)
		       (aref result1 (1+ index2)) (* (aref transform index2) factor)))
	;;F_0 = \int a' dsigma.  The Gaussian does not affect F_0.
	(setf (aref result1 0) (loop for index below points
				     sum (* (3vector-component (aref hats index) coordinate)
					    (- (aref sigmas index) (if (zerop index) 0 (aref sigmas (1- index)))))))
	(realft result1 (/ result-count 2) :isign -1) ;Inverse transformation in place
	;;Now RESULT contains RESULT-COUNT real numbers with slot k giving string position at location kL/result-count
	;;REALFT is defined to return the inverse Fourier transform multiplied by N = result-count/2.  The
	;;Numerical Recipes inverse Fourier transform is the sum divided by 2N.  Thus the sum itself is the
	;;the returned result times 2.  We want the sum divided by L, so we should scale by 2/L, but we're not going to
	;;bother because we have to normalize the 3vectors anyway.
	(dotimes (index result-count)
	  (setf (3vector-component (aref new-hats index) coordinate) (aref result1 index)))))
    (loop for index below result-count	;Now that all components have been installed, normalize resulting vectors
	  ;;	  do (setf (aref new-hats index) (3vector-scale (aref new-hats index) (/ 2 l))) ;unrenormalized
	  for old-length = (* (3vector-length (aref new-hats index)) (/ 2 l)) ;Length from actual calculation
	  minimize old-length into smallest-old-length
	  do (setf (aref new-hats index) (3vector-normalize (aref new-hats index))) ;renormalize
	  finally (format t "Largest renormalization factor ~F~%" (/ 1 smallest-old-length))
	  )
    (values new-hats new-sigmas)))

;;;Smoothing by small Lorentzians

;;Smooth evenly-spaced hats using Lorentzians.  The smoothing fraction is the scale parameter w
;;in units of loop length L.  The Lorentzian is w/pi/(w^2+x^2).  The Fourier transform of it is exp(-w |k|)
;;where k = 2 pi n/L
;;The individual steps are Lorentzians with width epsilon times the spacing between points.
;;If SIGMAS is not supplied, then HATS represents a piecewise-constant function with pieces all the same width.
;;If SIGMAS is supplied, it gives the sigma values for the end of each piece of the function.
;;In any case, we return the new HATS representing a piecewise-constant function with pieces all the same width.
;;If sigmas is given, N is number of data points to use in the calculation, so we will get n/2 frequencies
(defun Lorentzian-smooth (hats smoothing-fraction &key (epsilon 0.05) sigmas (n (length hats)))
  (let* ((steps (ceiling (* smoothing-fraction n) epsilon)) ;Number of smoothing steps
	 (step-width (/ smoothing-fraction steps)))	     ;Width each step as a fraction of L
    (check-rest-frame hats sigmas)
;;    (format t "~D steps:" steps) (force-output)
    (loop for step below steps
	  do (setq hats (Lorentzian-smooth-1 hats step-width n sigmas)
		   sigmas nil)			;Pass sigmas only first time
;;	  when (zerop (mod step 100))
;;	  do (format t "~D " step) (force-output)
	  finally (return hats))))

;;Smooth once by Lorentzian convolution.  The step-width is the scale parameter for the Lorentzian
;;in units of loop length, w1 = f1 L.  The mode of harmonic j has wavenumber k = 2 pi j/L, so it is damped
;;by the Fourier transform on the Lorentzian, exp(-w k) = exp(- 2 pi j f1)
;;You might worry that if f1 < 1/N, since j is at most N/2, the Fourier transform will never be small.
;;This represents the worry that the Lorentzian is too narrow to be accurately represented discretely, and so
;;we are actually convolving with something else.  But if we did not renormalize, repeated convolution would
;;yield a much wider Lorentzian, so it should be OK.
(defun Lorentzian-smooth-1 (hats step-width n sigmas)
  (let* ((nhats (length hats))
	 (data (make-array n :element-type 'double-float)) ;Array for Fourier transform
	 (new-hats (make-array n)))			   ;Array for result of this step
    (unless (= n (expt 2 (1- (integer-length n))))
      (error "The number of points must be a power of two instead of ~D" n))
    (dotimes (index n)
      (setf (aref new-hats index) (make-3vector))) ;Create 3vectors to return
    (dotimes (coordinate 3)			   ;Do cardinal directions separately
      (cond (sigmas				   ;Uneven spacing?
	     (Fourier-components-as-realft hats sigmas coordinate data)) ;nufft and install data
	    (t				;Even spacing
	     (dotimes (index nhats)	;Only use first part of array if fewer points
	       (setf (aref data index) (3vector-component (aref hats index) coordinate))) ;Copy data to transform
	     (when (> n nhats)
	       (fill data 0.0 :start  nhats)) ;Clear rest, if any
	     (realft data (/ nhats 2))))	  ;Transform nhats real numbers
      ;;Now DATA has the Fourier transform with frequency j stored in slots 2j, 2j+1.
      ;;It is not actually the Fourier transform of the piecewise-constant function, but the discrete Fourier transform of
      ;;the hat component as discrete data.  The right coefficients would be smaller by sinc(j/N) where sinc x = sin(pi x)/(pi x).
      ;;As a result when we do the inverse transform as though we were generating discrete data, we are in fact generating
      ;;the correct values for the smoothed piecewise-constant function.
      ;;Real values for frequencies 0 and n/2 are stored in slots 0 and 1
      ;;Rescale by Fourier transform of Lorentzian
      (loop for j from 1 to (/ n 2)
	    for coefficient = (exp (* -2 pi j step-width))
	    if (= j (/ n 2)) do (setf (aref data 1) (* coefficient (aref data 1)))
	    else do (setf (aref data (* 2 j)) (* coefficient (aref data (* 2 j)))
			  (aref data (1+ (* 2 j))) (* coefficient (aref data (1+ (* 2 j))))))
      (realft data (/ n 2) :isign -1)	;Inverse Fourier transform
      ;;The numbers in DATA should now be scaled by 2/n, but we aren't going to bother
      ;;because we have to normalize anyway.
      (dotimes (index n)
	(setf (3vector-component (aref new-hats index) coordinate) (aref data index))))
    (dotimes (index n)
      (setf (aref new-hats index) (3vector-normalize (aref new-hats index)))
;;Do not normalize: (setf (aref new-hats index) (3vector-scale (aref new-hats index) (/ 2.0 n)))
      )
    (close-hats new-hats)		;Renormalization might make sum not zero, so fix it.
    ))

;;Given hats and sigmas, do the non-uniform FFT for one component and install the result in the RESULT as though it were
;;the discrete Fourier transform of evenly-spaced a' values, in the format used by REALFT.  This is scaled differently
;;from the actual Fourier transform of the piecewise-constant function.  See Fourier.tex
;;We use the normalization used by REALFT, which is that transform of 1 has 0-frequency component N
;;We return 0 for frequencies 0 and N/2, the former because we should have a closed loop, and the latter because
;;you should have made N large enough that this component does not matter.
(defun Fourier-components-as-realft (hats sigmas coordinate result)
  (let* ((n (length result))		;Total number of frequencies
	 (points (length sigmas))	;Number of segments
	 (l (aref sigmas (1- points)))	;Length of loop
	 ;;To get n/2 frequencies accurately, we need to do the nufft with at least n positive frequencies
	 (nft (* n 2))			;so 2n frequencies in all.
	 (xi (make-array points :element-type 'double-float))	;Point locations for Fourier transform
	 (fi (make-array points :element-type 'double-float)))	;Function values for Fourier transform
    (dotimes (index points)				  ;Build xi
      (setf (aref xi index) (/ (aref sigmas index) l)))	  ;Scale into 0..1.  (aref xi 0) = sigma_1/L
    (dotimes (index points)
      (setf (aref fi index)				  ;fi = (ai - a_(i-1)).  See Fourier.tex
	    (- (3vector-component (aref hats (mod (1+ index) points)) coordinate)
	       (3vector-component (aref hats index) coordinate))))
    (let* ((transform (nufft xi fi *nufft-m* nft))) ;Do non-uniform Fourier transform
      ;;TRANSFORM is the sum of e^{-i omega_j sigma_{i+1} (f_{i+1} - f_i)}
      ;;Now if we divide by 2 i sin(pi j/N), we get something analagous to the sum of e^{-2 pi i l j / N}.  See
      ;;Fourier.tex.  To get something corresponding to the REALFT convention, we should now take the complex conjugate
      ;;Dividing by i and taking the complex conjugate is just exchanging the real and imaginary parts.
      ;;Frequency j in TRANSFORM is stored in slots 2j, 2j+1.
      ;;Real values for frequencies 0 and nft/2 are stored in slots 0 and 1
      ;;Slot 0 should contain (essentially) 0, because the data we passed to nufft is f_{i+1}-f_i, and
      ;;the sum of these quantities is zero, but we do not look.
      (setf (aref result 0) 0.0		;frequency 0
	    (aref result 1) 0.0)	;frequency N/2.  See header comment.
      (loop for j from 1 below (floor n 2)
	    for index2 = (* 2 j)
	    for divisor = (* 2 (sin (/ (* pi j) n)))
	    do (setf (aref result index2)     ;Real return
		     (/ (aref transform (1+ index2)) divisor)) ;from imaginary part
	    do (setf (aref result (1+ index2)) ;Imag return.
		     (/ (aref transform index2) divisor)))))) ;from real part

;;Signal an error unless the hats or hats and sigmas represent a loop in the rest frame
(defun check-rest-frame (hats &optional sigmas)
  (let ((total (make-zero-3vector))
	(length 0.0))
    (loop for index below (length hats)
	  for hat = (aref hats index)
	  for dsigma = (if sigmas (- (aref sigmas index) (if (zerop index) 0.0 (aref sigmas (1- index))))
			 1.0)
	  do (incf length dsigma)
	  do (setq total (prog1 (3vector+ total (3vector-scale hat dsigma))
			   (deallocate 3vectors total))))
    (let ((result (/ (3vector-length total) length)))
      (unless (< result 1e-4)
	(error "Your string is not in the rest frame.  Length of offset ~S" result)))))

;;Returns evenly spaced sigmas for functions that want that
(defun evenly-spaced-sigmas (length points)
  (let ((sigmas (make-array points :element-type 'double-float)))
    (dotimes (j points)
      (setf (aref sigmas j) (* (1+ j) (/ length points))))
    sigmas))


;;Make sure function works
(defun check-Fourier-components-as-realft ()
  (let ((hats (make-array 64))
	(data (make-array 64 :element-type 'double-float))
	(result (make-array 64 :element-type 'double-float))
	(sigmas (evenly-spaced-sigmas 3.0 64)))
    (loop for i below 64
	  for x = (make-3vector (- (random 2.0) 1.0) 0.0 0.0)
	  do (setf (aref hats i) x)
	  sum (3vector-x x) into total
	  finally (loop for i below 64
			do (setf (aref data i)
				 (decf (3vector-x (aref hats i)) (/ total 64)))) ;"rest frame"
	  )
    (Fourier-components-as-realft hats sigmas 0 result)
    (realft data 32)			;Fourier transform in place
    (setf (aref data 0) 0.0		;0, N/2 are zero in non-uniform case.
	  (aref data 1) 0.0)
    ;;The string in the DFT starts with a_0 at position 0, but really it should be at the central position sigma_1/2 = L/(2N)
    ;;To account for this, we change the phase of component j by pi j/N
    (print (loop for j from 1 below 32
		 collect (+ (expt (aref data (* j 2)) 2) (expt (aref data (1+ (* j 2))) 2))))
    (loop for j from 1 below 32
	  for index = (* j 2)
	  for angle = (/ (* j pi) 64)
	  for s = (sin angle)
	  for c = (cos angle)
	  do (psetf (aref data index) (- (* (aref data index) c) (* (aref data (1+ index)) s))
		    (aref data (1+ index)) (+ (* (aref data (1+ index)) c) (* (aref data index) s))))
    (print (loop for j from 1 below 32
		 collect (+ (expt (aref data (* j 2)) 2) (expt (aref data (1+ (* j 2))) 2))))
    (values result data)))
    

;;;Smoothing 3 points at a time

;;Total up an array of 3vectors
(defun 3vector-total (vectors)
  (loop with total = (make-zero-3vector)
	for vector across vectors
	do (3vector-incf total vector)
	finally (return total)))

;;Total up an array of 3vectors
(defun 3vector-total-length (vectors)
  (loop for vector across vectors
	sum (3vector-length vector)))

;;Total up an array of 4vectors
(defun 4vector-total (vectors)
  (loop with total = (make-zero-4vector)
	for vector across vectors
	do (4vector-incf total vector)
	finally (return total)))

(defun total-hats-dsigmas (hats dsigmas)
  (loop with total = (make-zero-3vector)
	for vector across hats
	for dsigma across dsigmas
	do (3vector-incf total (3vector-scale vector dsigma))
	finally (return total)))

;;Default tolerance for close-hats.  The total offset must be less than this multiplied by sqrt{N}
(defvar close-hats-default-tolerance 1e-15)

;;Offset some hats so that they total to 0 or some desired goal
;;Tolerance is the allowable total divided by the square root of the number of elements
(defun close-hats (hats &key dsigmas (tolerance close-hats-default-tolerance) (max-tries 12) (goal zero-3vector))
  (setq hats (copy-seq hats))
  (loop with n = (length hats)
	with last-hats = (make-array n) ;Save previous in case new is worse
	with length = (if dsigmas (total-sigma dsigmas) (double-float n))
	with Jacobian = (make-array (list 3 3) :element-type 'double-float) ;dtotal_i/doffset_j
	with indx = (make-array 3 :element-type '(unsigned-byte 16))	    ;Permutation array for LUDCMP
	with vv = (make-array 3 :element-type 'double-float)		    ;Working array for LUDCMP
	for count from 0
	for last-error = nil then error
	for total = (3vector- (if dsigmas (total-hats-dsigmas hats dsigmas) (3vector-total hats)) goal) ;to make 0
	for error = (/ (3vector-length total) (sqrt (double-float n))) ;error rescaled by sqrt{N} for adding in quadrature
;;	do (format t "Total ~S, length ~S rescaled ~D~%" total (3vector-length total) error)
	when (zerop error) return hats		 ;If perfect, exit now.
	when last-error				 ;Otherwise never exit immediately.  Always try once.
	do (cond ((>= error last-error)		 ;No improvement? Stop.
		  (when (> last-error tolerance) ;Last error was less.  Is it good enough?
		    (error "Failed to achieve tolerance ~S because error ~S did not improve in try ~D" tolerance last-error count))
		  (return last-hats))	;Previous try was better
		 ((>= count max-tries)	;Error improved, but we are out of tries
		  (when (> error tolerance)	
		    (error "Failed to achieve tolerance ~S.  Error is ~S after ~D tries." tolerance error count))
		  (return hats)))	;new values are better, but we won't try to improve them further
	;;If we offset each element by dx, the resulting change to a_i will be dx_i - sum(dx_j a_j) a_i, because we
	;;can only move in the perpendicular direction.  The change to the total will be the sum of these,
	;;each multiplied by its sigma. ;so we set Jacobian to the sum of sigma (delta_ij - a_i a_j).
	;;Now we solve Jacobian*x = total.  If we offset each vector by -x, then the total will become zero at linear order.
	do (dotimes (i 3)		;Reset to diagonal matrix of length.
	     (dotimes (j 3)
	       (setf (aref Jacobian i j) (if (= i j) length 0.0))))
	do (loop for index below n	;Subtract sum of a_i a_j
		 for hat = (aref hats index)
		 for sigma = (if dsigmas (aref dsigmas index) 1.0)
		 do (dotimes (i 3)
		      (dotimes (j 3)
			(decf (aref Jacobian i j) (* sigma (3vector-component hat i) (3vector-component hat j))))))
	unless (ludcmp Jacobian indx vv) ;Solve Jacobian*x = total
	  do (error "LUDCMP failed")
	do (lubksb Jacobian indx total) ;Now TOTAL becomes X
;;	do (format t "Subtracting ~S~%" total)
	do (dotimes (index n) ;Subtract X from each
	     (setf (aref last-hats index) (aref hats index)) ;Save old in case new worse
	     (setf (aref hats index) (3vector-normalize (3vector- (aref hats index) total))))))


;;Handling a string loop (in the rest frame) in terms of its Fourier components.  We handle components from 0
;;to N, but component 0 is always 0 in the rest frame.  Thus there are 6N numbers representing the 3 components
;;of the real and imaginary part of the Fourier amplitudes for N nonzero frequencies.
;;We store these in one vector in order a1_x, a1_y, a1_z, b1_x, b1_y, b1_z, a2_x,...
;;Conventions as in DeLaney and Brown, http://dx.doi.org/10.1103/PhysRevLett.63.474
;;Note unusual sign of B, factor of 2.

;;Give index into 6N-element vector of given component of A_j or B_j.
(defun Fourier-ab-index (ab i coordinate)
  (+ (* (1- i) 6) (if ab 3 0) coordinate))

;;Get N from an array of components passed as an argument
(defun Fourier-components-n (components)
  (let ((n (/ (length components) 6)))
    (unless (integerp n)
      (error "~S was not an array of 6N numbers" components))
    n))

;;Get N Fourier components for a tangent vector and return as a 6N-element vector
(defun get-Fourier-components (hats sigmas n)
  (let* ((points (length sigmas))	 ;Number of segments
	 (l (aref sigmas (1- points)))	 ;Length of loop
	 ;;To get n frequencies accurately, we need to do the nufft with at least 2n positive frequencies
	 (nft (expt 2 (ceiling (log (* n 4) 2))))	      ;so 4n frequencies in all.  Round up to power of 2
	 (xi (make-array points :element-type 'double-float))	;Point locations for Fourier transform
	 (fi (make-array points :element-type 'double-float))	;Function 1values for Fourier transform
	 (result (make-array (* n 6))))
    (dotimes (index points)				  ;Build xi
      (setf (aref xi index) (/ (aref sigmas index) l)))	  ;Scale into 0..1
    (dotimes (coordinate 3)				  ;Do cardinal directions separately
      (dotimes (index points)				  ;Build fi for this coordinate
	(setf (aref fi index)				  ;fi = (ai - a(i-1))
	      (- (3vector-component (aref hats (mod (1+ index) points)) coordinate)
		 (3vector-component (aref hats index) coordinate))))
      (let* ((transform (nufft xi fi *nufft-m* nft))) ;Do non-uniform Fourier transform
	;;TRANSFORM is the sum of e^{-i omega_j sigma_{i+1} (f_{i+1} - f_i)}
	;;C_n is this, multiplied by -i and divided by omega_j L
	;;A_n = 2 Re(C_n) = imaginary return value * 2/(omega_j L)
	;;B_n = -2 Im(C_n) = real return value * 2/(omega_j L)
	;;Frequency j in TRANSFORM is stored in slots 2j, 2j+1.
	;;Real values for frequencies 0 and nft/2 are stored in slots 0 and 1
	;;Slot 0 should contain (essentially) 0, because the data we passed to nufft is f_{i+1}-f_i, and
	;;the sum of these quantities is zero, but we do not look.
	(loop for j from 1 to n
	      for index2 = (* 2 j)
	      for divisor = (* pi j) ;Omega_j L/2 = pi j
	      do (setf (aref result (Fourier-ab-index nil j coordinate)) ;Aj
		       (/ (aref transform (1+ index2)) divisor))
	      do (setf (aref result (Fourier-ab-index t j coordinate)) ;Bj
		       (/ (aref transform index2) divisor)))))
    result))

;;Generate a string loop from the given Fourier component vector
;;Returns hats and sigmas, without normalizing the hats.
;;POINTS gives the number of points to use, but this is increased to at least the number of frequencies in use,
;;and rounded up to power of 2, so we can do the Fourier transform.
(defun string-from-Fourier-components (components l points)
  (let* ((nfreq (Fourier-components-n components))		     ;Number of frequencies (not including 0)
	 (count (expt 2 (ceiling (log (max points (* nfreq 2)) 2)))) ;Round up to power of 2
	 (result1 (make-array count :element-type 'double-float))    ;Place to do inverse transform
	 (sigmas (make-array count :element-type 'double-float))     ;New sigmas to return
	 (hats (make-array count)))				     ;Array of 3vectors to return.
    (loop for index below count
	  do (setf (aref sigmas index) (/ (* l (1+ index)) count)) ;New sigmas evenly spaced
	  do (setf (aref hats index) (make-3vector)))		   ;Create 3vectors to return
    (dotimes (coordinate 3)					       ;Do cardinal directions separately
      (fill result1 0.0)					       ;Clear unused frequencies
      (loop for j from 1 to nfreq
	    do (setf (aref result1 (* j 2)) ;Real part.
		     (/ (aref components (Fourier-ab-index nil j coordinate)) 2)) ;Fix norm.
					;Imaginary part.  Opposite sign convention cancels difference in FT convention
	    do (setf (aref result1 (1+ (* j 2)))
		     (/ (aref components (Fourier-ab-index t j coordinate)) 2)))
      (realft result1 (/ count 2) :isign -1)			      ;Inverse FFT in place
      ;;Now RESULT contains COUNT real numbers with slot k giving string position at location (k/count)L.
      ;;REALFT is defined to return the inverse Fourier transform multiplied by N = count/2, but the FT convention
      ;;of numerical recipes is sum_1^count(C_n e^{i omega_n x})/count, so we must multiply by 2.
      (dotimes (index count)
	(setf (3vector-component (aref hats index) coordinate) (* (aref result1 index) 2))))
    (values hats sigmas)
    ))

;;More code from here is in old.lisp

(defun plot-Fourier-magnitudes (magnitude-lists names)
  (gnuplot (length magnitude-lists) (length (first magnitude-lists))
	   #'(lambda (plot point)
	       (if (eq point :title) (nth plot names)
		 (values (1+ point) (nth point (nth plot magnitude-lists)))))
	   :logscale :y
	   :xrange (cons 0  (1+ (length (first magnitude-lists))))
	   :styles :linespoints))

(defun plot-Fourier-magnitudes-fit (magnitude-lists names)
  (let ((n (length magnitude-lists))
	(fits (loop for magnitudes in magnitude-lists
		    collect (multiple-value-list (fit-Fourier-magnitudes magnitudes)))))
    (gnuplot (* 2 n) (length (first magnitude-lists))
	     #'(lambda (plot point)
		 (cond ((eq point :title)
			(nth plot names)) ;Name points, not fits
		       ((< plot n)	  ;Point
			(values (1+ point) (nth point (nth plot magnitude-lists))))
		       (t (destructuring-bind (b a) (nth (- plot n) fits)
			    (values (1+ point) (exp (- a (* b (1+ point)))))))))
	     :logscale :y
	     :xrange (cons 0  (1+ (length (first magnitude-lists))))
	     :styles (append (make-list n :initial-element :points)
			     (loop for index from 1 to n
				   collect (format nil "lines lt ~D" index))))))

;;Compute the "wiggliness parameter": sum_{j=2}^infinity j^2 a_j^2
;;include j = 1, because we don't want to call a circular string wiggly.
(defun Fourier-wiggliness (magnitudes)
  (loop for j from 2
	for magnitude in (cdr magnitudes)
	sum (expt (* j magnitude) 2)))

;;If we took these magnitudes and smoothed them by exp(- alpha j), what wiggliness would we have?
(defun Fourier-wiggliness-goal (magnitudes alpha)
  (Fourier-wiggliness (loop for j from 1 for magnitude in magnitudes
			    collect (/ magnitude (exp (* alpha j))))))


;;A least-squares fit to the decline of Fourier component magnitudes with frequency
;;Maybe the component magnitudes look like A exp(-bj), so their logarithms like a - b j
;;Fits and returns B and A as values.
(defun fit-Fourier-magnitudes (magnitudes)
  (loop for magnitude in magnitudes
	for x from 1
	for y = (log magnitude)
	sum 1 into s
	sum x into sx
	sum y into sy
	sum (* x x) into sxx
	sum (* x y) into sxy	
	finally
	(return (values (- (/ (- (* s sxy) (* sx sy)) (- (* s sxx) (* sx sx))))
			(/ (- (* sxx sy) (* sx sxy)) (- (* s sxx) (* sx sx)))))))

;;Return a list of total magnitudes (Euclidean length of 6-vector) of components for each harmonic
(defun Fourier-component-magnitudes (components)
  (loop for j below (fourier-components-n components)
	collect (sqrt (loop for i from (* j 6) repeat 6 sum (expt (aref components i) 2)))))


;;Smoothing by direct calculation.  This does not do Gaussian convolution, but rather
;;uses a cheap system of considering each segment to be located at the sigma of its midpoint.
(defun slow-smooth-data (hats sigmas standard-deviation spacing)
  (let* ((loop-length (aref sigmas (1- (length sigmas)))) ;Total sigma
	 (count (round loop-length spacing)) ;Number of new points
	 (step (/ loop-length count))	;Actual spacing
	 (smooth-hats (make-array count))
	 (smooth-sigmas (make-array count :element-type 'double-float)))
    (loop for index below count
	  for sigma = (* step index)
	  do (setf (aref smooth-sigmas index) sigma
		   (aref smooth-hats index) (slow-smooth-data-1 hats sigmas sigma standard-deviation)))
    (values smooth-hats smooth-sigmas)))

;;Return a'(sigma) computed with Gaussian smoothing, sort of
(defun slow-smooth-data-1 (hats sigmas sigma standard-deviation)
  (let* ((count (length sigmas)) 
	 (loop-length (aref sigmas (1- count)))) ;Total sigma
    (unless (< (* standard-deviation 2) loop-length)
      (error "Loop length is only ~F, so I can't smooth with range ~F" loop-length standard-deviation))
    (flet ((do-mod (x)			;Map x into -L/2...L/2
	     (- (mod (+ x (/ loop-length 2)) loop-length) (/ loop-length 2))))
      (loop with result = zero-3vector
	    for index from 0 below count
	    for length = (- (aref sigmas index) ;Length of this segment
			    (if (zerop index) 0.0 (aref sigmas (1- index))))
	    for center = (- (aref sigmas index) (/ length 2)) ;sigma of center of this segment
	    for distance = (do-mod (- center sigma))	      ;signed distance away from sigma of evaluation
	    for weight = (* length
			    (loop for wrap from -5 to 5 ;Wrap Gaussian
				  as wrapped-distance = (+ distance (* wrap loop-length))
				  sum (exp (- (/ (expt wrapped-distance 2) 2 (expt standard-deviation 2))))))
	    do (setq result (3vector+ result (3vector-scale (aref hats index) weight)))
	    finally (return (3vector-normalize result))))))

;;Length of loop or NIL if not a loop
(defun loop-length (start-diamond)
  (setq start-diamond (diamond-next-a start-diamond)) ;Advance to diamond that starts a new A
  (and start-diamond
       (loop with diamond = start-diamond
	     do (setq diamond (diamond-next-a diamond)) ;Get next A segment
	     while (diamondp diamond)
	     ;;The vector q points, say, from sigma=t=0 to sigma=t=q_t, since sigma-t is fixed along the right edge
	     ;;Thus the argument of b(t+sigma) is 2 q_t
	     sum (* 2 (3vector-length (diamond-a-wrap diamond))) into length
	     when (eq diamond start-diamond) ;Looped around?  This last segment has been done
	     return length)))

;;Loops in decreasing order of length.  Returns a list of (DIAMOND LENGTH)
(defun loops-by-length ()
  (let ((result nil))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore)))
     #'(lambda (d)
	 (push (list d (loop-length d)) result)))
    (sort result #'> :key #'second)))

(defun loop-lengths ()
  (mapcar #'second (loops-by-length)))

;;Return a diamond in the longest loop.  If an argument is given,
;;return the nth longest.  The longest is number 0.
(defun longest-loop (&optional (nth 0))
  (car (nth nth (loops-by-length))))
	   
;;Loop with a given number of diamonds
(defun loop-of-length (count &optional (nth 1))
  (let ((serial 0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore)))
     #'(lambda (d)
	 (let ((this (loop-count d)))
	   (when (and (= this count)
		      (= (incf serial) nth))
	     (return-from loop-of-length d)))))
    (if (zerop serial)
	    (error "There is no loop of length ~D" count)
      (error "There are only ~D loop~:P of length ~D" serial count))))
  
;;Compute average angle between changes in a' at successive kinks
(defun average-kink-kink-angle (diamond &optional (n 1000))
  (let* ((center diamond)
	 (center-a (3vector-normalize (diamond-a-wrap center)))
	 (right (diamond-next-a center))
	 (right-a (3vector-normalize (diamond-a-wrap right)))
	 (length-2 (3vector-dot center-a right-a)) ;cosine of segment from center to right
	 length-1
	 left-a
	 far-length
	 angle)
    (loop repeat n
	  do (setq length-1 length-2	;Move quantities one step to the left
		   left-a center-a
		   center-a right-a
		   right (diamond-next-a right) ;Get new diamond
		   right-a (3vector-normalize (diamond-a-wrap right))
		   length-2 (3vector-dot center-a right-a)
		   far-length (3vector-dot left-a right-a) ;Length of far side
		   angle (acos (- (min (/ (- far-length (* length-1 length-2)) ;Spherical law of cosines gives angle at center
					  (sqrt (- 1 (expt length-1 2)))
					  (sqrt (- 1 (expt length-2 2))))
				    1.0))))
;;	  do (format t "l: ~F~%c: ~F~%r ~F~%" left-a center-a right-a)
;;	  do (format t "1: ~F, 2: ~F, Far ~F~%" (acos length-1) (acos length-2) (acos far-length))
;;	  do (format t "~F " angle)
	  sum angle into total
	  finally (return (/ total n)))))
	  

#|
      (defun random-unit-disk-point ()
(loop for x = (- (random 2.0) 1)
for y = (- (random 2.0) 1)
when (< (+ (expt x 2) (expt y 2)) 1)
return (list x y)))

      (defun unit-disk-average-angle (&optional (n 1000))
(loop repeat n
for (x y) = (random-unit-disk-point)
for (x1 y1) = (random-unit-disk-point)
for (x2 y2) = (random-unit-disk-point)
for d1 = (make-3vector (- x1 x) (- y1 y) 0.0)					 ;Don't have 2vectors
for d2 = (make-3vector (- x2 x) (- y2 y) 0.0)					 ;Don't have 2vectors
sum (3vector-dot (3vector-normalize d1) (3vector-normalize d2))	into sum	 ;cos theta
finally (return (/ sum n))))
      |#

;;Cusps, half-cusps and the like.

(mirror-images
(defun a-kinks ()
  (let ((kinks nil))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only loops
     #'(lambda (diamond)
	 (setq diamond (diamond-next-a diamond)) ;Start with a new a
	 (loop with d = diamond
	       do (push (diamond-a-kink-created d) kinks)
	       do (setq d (diamond-next-a d))
	       until (eq d diamond))))
    kinks)))

;;Return longest sequence of old kink creation positions
(defun longest-old-a-sequence (diamond &optional (threshold-time 4.5))
  (setq diamond (diamond-next-a diamond)) ;Skip to beginning of a segment
  (let (best-list
	(best-steps 0)
	best-start 
	(list nil)
	(steps 0)			;Number of steps in sequence
	(startup t)			;Skipping excursions found at the beginning
	(d diamond)
	(looped nil))
    (loop for kink = (diamond-a-kink-created d)
	  for index from 0
	  do
	  (cond ((< (4vector-t kink) threshold-time) ;old kink
		 (unless startup	;unless still starting
		   (push kink list)
		   (incf steps)))
		(t			;new kink
		 (setq startup nil)	;not still starting
		 (when (> steps best-steps) ;this is longer than previous
		   (setq best-list list best-steps steps
			 best-start (- index steps))) ;Starting index of best list
		 (setq steps 0 list nil) ;Forget old list
		 (when looped		;looped around?
		   (return (values best-list best-start)))))
	  (setq d (diamond-next-a d))	;Advance to next a
	  (when (eq d diamond)		;looped around?
	    (setq looped t)))))		;Exit when we can


;;Count places where a' goes back and forth between 2 values
;;If FIND, returns diamonds
(defun count-excursions (&key (min-length (current-time)) (excursion-angle 0.1) (return-angle 0.01) find)
  (let ((result nil)
	(sigma 0.0)
	(count 0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
     #'(lambda (d)
	 (multiple-value-bind  (length as) (loop-length-and-count-a d)
	   (when (> length min-length)	;Not short loops
	     (incf sigma length)
	     (incf count as)
	     (setq result (nconc (find-excursions d nil excursion-angle return-angle) result))))))
    (let ((found (length result)))
      (format t "In ~D a values covering sigma ~6F we found ~D excursion~:P.  Average separation com ~6F scaling ~7E.~%Counts: "
	      count sigma found (/ sigma found) (/ sigma found (current-time))))
    (let ((counts (make-array (1+ (loop for (steps) in result maximize steps)) :initial-element 0))
	  (found nil))
      (loop for (steps d) in result
	    do (incf (aref counts steps))
	    when (and find (>= steps find))
	    do (push (list steps d) found))
      (loop for count across counts
	    for steps from 0
	    when (plusp count)
	    do (format t "~D(~D) " count steps))
      (terpri)
      found)))


;;Find excursions and return list of (steps diamond p len p len p len p ...)
(defun find-excursions (diamond count excursion-angle return-angle)
  (setq diamond (diamond-next-a diamond)) ;Skip to beginning of a segment
  (let ((results nil)
	(steps 0)			;Number of steps in excursion.  Out and back is 2 steps.
	(start nil)			;Where excursion started
	start-diamond
	excursion			;Where it went
	set
	(sigma 0.0)
	(startup t)			;Skipping excursions found at the beginning
	(d diamond)
	(looped nil))
    (loop for a = (diamond-a-wrap d)
	  for length = (* (3vector-length a) 2)
	  for this = (3vector-normalize a)
	  for index from 1
	  do
	  (incf sigma length)
	  (cond ((null start)		;Very first time
		 (setq start this	;Copy first a and that's it
		       start-diamond d))
		((zerop steps)		;No excursion in progress
		 (cond ((< (3vector-dot start this) (cos excursion-angle)) ;Kink is larger than minimum angle?
			(setq excursion this steps 1) ;Excursion potentially beginning here.
			(setq set (list this)))
		       (t (setq start this start-diamond d
				startup nil)))) ;Not starting excursion: keep scanning
		((<= (spherical-angle this (if (evenp steps) excursion start))
		     return-angle)	;Does excursion continue to alternate between two places?
		 (incf steps)		;Count one more step
		 (push this set)
		 (push length set))
		(t			;Excursion over: record and reset
		 (when (> steps 1)	;Don't record 1-step excursions
		   (unless startup	;Don't record an excursion that we find that at the very beginning
		     (push (list steps start-diamond (nreverse set)) results)))
		 (setq start this steps 0 startup nil start-diamond d))) ;Ready to look for a new excursion
	  (setq d (diamond-next-a d))	;Advance to next a
	  (when (eq d diamond)		;looped around?
	    (setq looped t))		;Exit when we can
	  (when (and (zerop steps)		;Not in the middle of an excursion
		     (or looped (and count (>= index count))))
	    (return results)))))

;;Give exterior (bending) angle between two lines on the unit sphere
(defun kink-kink-angle (left center right)
  (- pi (spherical-triangle-angle center right left)))

;;Sign is positive if path on the unit sphere, seen from the inside, turns to the left
(defun signed-kink-kink-angle (one two three)
  (let ((angle (signed-spherical-triangle-angle two three one)))
    (if (minusp angle)			;Same sign, supplement of angle
	(- (+ pi angle))
      (- pi angle))))

;;Find segments which have no more than maximum-angle kink-kink angles
;;Returns list of (INDEX STEPS SIGMA SPREAD)
;; STEPS is the number of kinks that have small bends between them
;; SIGMA is the total sigma including the segments which are at the centers of the large bends
;; SPREAD is the total of the kink angles
(defun find-straight-pieces (diamond &key (max-bend (/ pi 2)) min-sigma (min-steps 2) count skip (loop t)
				     (big-kink-fraction 0.75)) ;If a single kink is more than this, split piece
  (assert min-sigma)
  (multiple-value-bind (hats sigmas) (get-a-data-sigmas diamond :count count :skip skip)
    (let ((n (length hats))
	  (result nil)
	  (index 0)			;Index of center of bend
	  (looped nil)
	  start-index
	  (startup t))			;Skipping segments found at the beginning
      ;;Process a piece of string.  At positions FROM and TO there are large bends, and small between.
      ;;This is called for every large bend
      (labels ((do-piece (from to)
	         (let ((steps (mod (- to from) n)))
		   (when (>= steps min-steps) ;If enough steps
		     ;;Length in sigma including both ends
		     (let ((sigma (- (aref sigmas to) (aref sigmas (mod (1- from) n)))))
		       (when (minusp sigma) (incf sigma (aref sigmas (1- (length sigmas))))) ;Wrap
		       (when (>= sigma min-sigma) ;If long enough in sigma
			 (loop for i = from then (mod (1+ i) n) until (= i to) ;Loop over kinks
			       as angle = (spherical-angle (aref hats i) (aref hats (mod (1+ i) n)))
			       with max-angle = 0
			       with max-angle-index
			       sum angle into spread
			       when (>= angle max-angle) ;Accumulate largest angle
			       do (setq max-angle angle max-angle-index i)
			       finally
			       (cond ((> (/ max-angle spread) big-kink-fraction) ;One kink accounts for too much
				      (do-piece from max-angle-index) ;Split pieces around big kink, try again
				      (do-piece (mod (1+ max-angle-index) n) to))
				     (t
				      (push (list from steps sigma spread) result))))))))))
	(loop
	 (when (> (kink-kink-angle (aref hats (mod (1- index) n)) (aref hats index) (aref hats (mod (1+ index) n)))
		  max-bend)		;Large bend here?
	   (unless startup (do-piece start-index index)) ;See if piece ending here is interesting
	   (when looped			;Not accumulating, so see if time to exit
	     (return result))
	   (setq startup nil start-index index)) ;New straight piece might start here
	 ;;If small angle, just keep scanning
	 (incf index)
	 (when (= index (length hats))
	   (when startup (return result)) ;Never found any large bend, so return instead of looping forever
	   (if loop
	       (setq index 0 looped t)	;Exit when we can
	     (return result))))))))

(defun find-straight-pieces-hats-sigmas (hats sigmas &key (max-bend (/ pi 2)) (min-sigma 0.001) (min-steps 2) (loop t)
				     (big-kink-fraction 0.75)) ;If a single kink is more than this, split piece
  (assert min-sigma)
  (let ((n (length hats))
	(result nil)
	(index 0)			;Index of center of bend
	(looped nil)
	start-index
	(startup t))			;Skipping segments found at the beginning
    ;;Process a piece of string.  At positions FROM and TO there are large bends, and small between.
    ;;This is called for every large bend
    (labels ((do-piece (from to)
		       (let ((steps (mod (- to from) n)))
			 (when (>= steps min-steps) ;If enough steps
			   ;;Length in sigma including both ends
			   (let ((sigma (- (aref sigmas to) (aref sigmas (mod (1- from) n)))))
			     (when (minusp sigma) (incf sigma (aref sigmas (1- (length sigmas))))) ;Wrap
			     (when (>= sigma min-sigma) ;If long enough in sigma
			       (loop for i = from then (mod (1+ i) n) until (= i to) ;Loop over kinks
				     as angle = (spherical-angle (aref hats i) (aref hats (mod (1+ i) n)))
				     with max-angle = 0
				     with max-angle-index
				     sum angle into spread
				     when (>= angle max-angle) ;Accumulate largest angle
				     do (setq max-angle angle max-angle-index i)
				     finally
				     (cond ((> (/ max-angle spread) big-kink-fraction) ;One kink accounts for too much
					    (do-piece from max-angle-index) ;Split pieces around big kink, try again
					    (do-piece (mod (1+ max-angle-index) n) to))
					   (t
					    (push (list from steps sigma spread) result))))))))))
	    (loop
	     (when (> (kink-kink-angle (aref hats (mod (1- index) n)) (aref hats index) (aref hats (mod (1+ index) n)))
		      max-bend)		;Large bend here?
	       (unless startup (do-piece start-index index)) ;See if piece ending here is interesting
	       (when looped			;Not accumulating, so see if time to exit
		 (return result))
	       (setq startup nil start-index index)) ;New straight piece might start here
	     ;;If small angle, just keep scanning
	     (incf index)
	     (when (= index (length hats))
	       (when startup (return result)) ;Never found any large bend, so return instead of looping forever
	       (if loop
		   (setq index 0 looped t)	;Exit when we can
		 (return result)))))))


;;Returns number of half-cusps, average steps, average sigma, average spread, total number of segments
;;and total length
(defun count-straight-pieces (&rest args)
  (let ((count 0)
	(total-steps 0)
	(total-sigma 0)
	(total-spread 0.0)
	(total-length 0.0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
     #'(lambda (diamond)
	 (loop for (index steps sigma spread) in (apply #'find-straight-pieces diamond args)
	       do (incf count)
	       do (incf total-steps steps)
	       do (incf total-sigma sigma)
	       do (incf total-spread spread))
	 (multiple-value-bind (len steps) (loop-length-and-count-a diamond)
	   (incf total-length len)
	   (incf total-steps steps))))
    (values count (/ (double-float total-steps) count) (/ total-sigma count) (/ total-spread count)
	    total-steps total-length)))

(defun cusps-from-straight-pieces (&key (max-bend (/ pi 2)) min-sigma (min-steps 2)
					(alpha 0.05) ;length of loops as a fraction of current time
					(time (current-time)))
  (multiple-value-bind (cusps steps sigma spread total-steps total-length)
      (count-straight-pieces :max-bend max-bend :min-sigma min-sigma :min-steps min-steps)
    (declare (ignore steps sigma))
    (format t "~&~D total steps in length ~S~%" total-steps total-length)
    (cusps-from-straight-pieces-1 cusps spread total-length alpha time)))

(defun cusps-from-straight-pieces-1 (cusps spread total-length alpha time)
  (let* ((density (/ cusps total-length)) ;qualifying half-cusps per unit length
	 (per-loop (* density alpha (current-time))) ;Number of half-cusps for loop
	 (pairs (expt per-loop 2))	;Number of half-cusp pairs per loop
	 (probability (* (expt spread 2) (/ 2 pi pi))) ;Probability of 2 half-cusps intersecting
	 (result (* pairs probability))) ;Average number of cusps per loop
    (format t "~&~D half-cusps in total length ~F, so density is ~D~%" cusps total-length density)
    (format t "Loop of size ~F has typically ~D half-cusps with intersection probability ~F"
	    (* alpha time) per-loop probability)
    result))
					
;;Supposing every old kink to be a half-cusp, how many cusps would there be?
(defun cusps-from-old-kinks (&key (creation-threshold 4.5))
  (let ((total-spread 0.0)
	(old-count 0)
	(total-count 0)
	(total-length 0.0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
     #'(lambda (diamond)
	 (multiple-value-bind (a-hats a-sigmas a-created) (get-a-data-sigmas diamond)
	   (loop for i below (length a-hats)
		 do (incf total-count)
		 when (< (aref a-created i) creation-threshold) ;Old kink in a?
		 do (incf old-count)
		 and do (incf total-spread (spherical-angle (aref a-hats i) (aref a-hats (mod (1+ i) (length a-hats))))))
	   (incf total-length (aref a-sigmas (1- (length a-sigmas)))))))
    (let* ((spread (/ total-spread old-count)) ;average spread
	   (density (/ old-count total-length)) ;density of old kinks
	   (probability (* (expt spread 2) (/ 2 pi pi))) ;Probability of old kinks crossing
	   (total-expected-crossings 0.0)
	   (total-actual-crossings 0)
	   (loop-count 0)
	   (result nil))
      (format t "~&Total of ~D old kinks (out of ~D total)" old-count total-count)
      (format t "~&Old kink density ~S, spread ~S, intersection probability ~S~%" density spread probability)
      (map-string-paths
       #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
       #'(lambda (diamond)
	   (when (>= (diamond-countup diamond) 2) ;Non-self-intersecting loops
	     (let ((length (loop-length-and-count-a diamond)))
	       (when (< length *total-size*) ;No wrapping loops
		 (let* ((expected (* density length))
			(pairs (expt expected 2)) ;expected number of pairs on this loop
			(crossings (* pairs probability))) ;expected number of crossing pairs
		   (format t "Loop segs ~D len ~8F made ~$.  Expct. count ~8F, crossings ~8F.  Actual crossings"
			   (loop-count diamond) length (4vector-t (tag-created-position (diamond-tag diamond)))
			   expected crossings)
		   (force-output)
		   (incf total-expected-crossings crossings)
		   (incf loop-count)
		   (let ((actual (length (find-old-crossings diamond :creation-threshold creation-threshold))))
		     (format t " ~D~%" actual)
		     (push (list (4vector-t (tag-created-position (diamond-tag diamond))) actual)
			   result)			       
		     (incf total-actual-crossings actual))))))))
      (format t "Total actual ~D, total expected ~$, ratio ~8F~%" total-actual-crossings total-expected-crossings
	      (/ total-actual-crossings total-expected-crossings))
      (format t "~D crossings in ~D loops, ~$ average ~%" total-actual-crossings loop-count
	      (/ total-actual-crossings loop-count))
      result)))


(defun find-half-cusps (&rest args)
  (let ((coverage 0.0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
     #'(lambda (diamond)
	 (incf coverage (apply #'find-half-cusps-1 diamond args))))
    coverage))

;;Look for half-cusps meeting giving criteria.  Returns fraction of unit sphere covered
;;by circles with half-cusp segments as diameters.
(defun find-half-cusps-1 (diamond &key min-sigma min-steps min-angle
				  (print t) (excursion-threshold 0.01))
  (multiple-value-bind (hats sigmas) (get-a-data-sigmas diamond)
    (loop with length = (loop-count diamond)
	  with last-last = nil		;Last index of previously reported segment
	  with sum-angle2 = 0.0
	  with count = 0
	  for start below (length hats)
	  do (multiple-value-bind (last excursions excursion-length)
		 (find-half-cusp hats sigmas start :excursion-threshold excursion-threshold)
	       (let ((steps (- last start))
		     (sigma (- (aref sigmas last) (aref sigmas start))) ;Sigma in good points
		     (angle (acos (min 1.0 (3vector-dot (aref hats last) (aref hats start))))))
		 (when (< last start)		;Wrapped?
		   (incf steps (length hats))
		   (incf sigma (aref sigmas (1- (length hats)))))
		 (when (and (or (not min-steps) (>= (- steps (* 2 excursions)) min-steps))
			    (or (not min-sigma) (>= sigma min-sigma))
			    (or (not min-angle) (>= angle min-angle))
			    ;;Don't report subsegment (starting later and ending same or earlier)
			    (not (and last-last (<= last last-last)))
			    (loop for i from start below last
				  never (> (acos (min 1.0 (3vector-dot (aref hats (1+ i)) (aref hats i))))
					   (* angle 2)))) ;Reject if one angle more than half the total
		   (when print
		     (format t "Len ~D idx ~D: ~D steps with length ~G and angle ~G~[~:; and ~:*~D excursion~:P of length ~G~]~%" length start steps sigma angle excursions excursion-length))
		   (incf sum-angle2 (* angle angle))
		   (incf count)
		   (setq last-last last))))
	  finally (return (/ sum-angle2 16)))))

;;See if there is a half-cusp here, looking only forward
;;Returns index of last a' that was further than its predecessor from the starting point and count and sigma of excursions.
;;Single-point excursions shorter than EXCURSION-THRESHOLD relative to accumulated sigma so far will be ignored
(defun find-half-cusp (hats sigmas index &key (excursion-threshold 0.01))
  (loop with total-sigma = (aref sigmas (1- (length sigmas)))
	with excursion-sigma = 0.0
	with excursions = 0
	with excursions-sigma = 0.0
	with start = (aref hats index)	;Starting direction
	with start-sigma = (aref sigmas index) ;End of first a'
	for last-j = index then j
	for j = (mod (1+ last-j) (length hats)) ;Advance and wrap
	for last-sigma = nil then sigma	;Start of this segment
	for sigma = (aref sigmas j)	;End of this a'
	with last-angle = nil
	for angle = (acos (min 1.0 (3vector-dot start (aref hats j))))
	do (cond ((or (null last-angle) (> angle last-angle)) ;New point is further than original.  Good.
		  (when (plusp excursion-sigma) ;On excursion?
		    (incf excursions)	;Count it
		    (incf excursions-sigma excursion-sigma) ;And its length
		    (setq excursion-sigma 0.0)) ;No longer on excursion
		  (setq last-angle angle)) ;Save new largest angle
		 ((plusp excursion-sigma) ;Already on excursion?
		  (return (values last-j excursions excursions-sigma)))	;Done
		 (t (setq excursion-sigma (mod (- sigma last-sigma) total-sigma)) ;Length of this segment
		    (when (> excursion-sigma
			     (* (mod (- sigma start-sigma) total-sigma) excursion-threshold)) ;Too long for excursion
		      (return (values last-j excursions excursions-sigma))) ;Done
		    ;;Otherwise continue scanning in case excursion over.  Don't reset last-angle.
		    ))))

;;Find old kinks that cross each other on a loop, forming "cusps"
(defun find-old-crossings (diamond &key count skip (creation-threshold 4.5))
  (multiple-value-bind (a-hats a-sigmas a-created) (get-a-data-sigmas diamond :count count :skip skip)
    (declare (ignore a-sigmas))
    (multiple-value-bind (b-hats b-sigmas b-created) (get-b-data-sigmas diamond :count count :skip skip)
      (declare (ignore b-sigmas))
      ;;Brute force
      (loop for i below (length a-hats)
	    when (< (aref a-created i) creation-threshold) ;Old kink in a?
	    nconc (loop with a1 = (aref a-hats i)
			with a2 = (aref a-hats (mod (1+ i) (length a-hats)))
			for j below (length b-hats)
			when (and (< (aref b-created j) creation-threshold) ;Old kink in b?
				  (spherical-intersection a1 a2 (aref b-hats j)
							  (aref b-hats (mod (1+ j) (length b-hats)))))
			collect it)))))

;;Plot tangent vectors

(defun plot-a (diamond &rest gnuplot-keys
		       &key count skip (style :lines)
		       creation-threshold ;Flag new vs. old kinks at this threshold
		       loop
		       smooth-range smooth-interval
		       &allow-other-keys)
  (multiple-value-bind (a-hats a-sigmas a-created)
      (get-a-data-sigmas diamond :count count :skip skip :smooth-range smooth-range :smooth-interval smooth-interval)
    (declare (ignore a-sigmas))
    (apply #'plot-tangent-vectors
	   (if creation-threshold
	       (multiple-value-bind (a-old a-new)
		   (split-by-creation-time a-hats a-created creation-threshold :loop loop)
		 (list (list a-old "A old" (plot-ab-style style nil nil))
		       (list a-new "A new" (plot-ab-style style nil t))))
	     (list (list a-hats "A" style)))
	   gnuplot-keys)))

;;Take an array of a-hats or b-hats and interpolate points so that great circles come out curved in projection
;;Either make all angles less than max-angle or split each segment into a split-count pieces
(defun interpolate-hats (hats &key max-angle split-count (loop t))
  (assert (or max-angle split-count))
  (assert (not (and max-angle split-count)))
  (loop with new = (make-array (* 2 (length hats)) :adjustable t :fill-pointer 0)
	for index below (length hats)
	for this = (aref hats index)
	for next = (aref hats (next-array-index-wrapping index hats))
	do (vector-push-extend this new)
	when (or loop (< index (1- (length hats)))) ;Don't interpolate if after last if not looping
	do (let* ((angle (spherical-angle next this)))
	     (when max-angle
	       (setq split-count (ceiling angle max-angle))) ;Number of segments to use
	     (loop for j from 1 repeat (1- split-count)
		   do (vector-push-extend (spherical-interpolation this next (/ (double-float j) split-count))
					  new)))
	finally (return new)))

;;Return a vector of (hat radius) for :style :circle
(defun hats-circle-data (hats dsigmas)
  (map 'vector #'(lambda (hat dsigma)
		   ;;1.0 -> radius 0.06.  <= 1e-5 -> radius 0.01
		   (list hat (/ (max (+ (log dsigma 10.0) 6.0) 1.0) 100.0)))
       hats dsigmas))

;;See old.lisp for version that splits by creation time
(defun plot-ab (diamond &rest keys &key count skip
			smooth		;Not implemented any more
			&allow-other-keys)
  (when smooth
    (error "Smoothing is no longer implemented in PLOT-AB"))
  (multiple-value-bind (a-hats a-dsigmas)
      (get-a-data-dsigmas diamond :count count :skip skip)
    (multiple-value-bind (b-hats b-dsigmas)
	(get-b-data-dsigmas diamond :count count :skip skip)
      (apply #'plot-ab-1 a-hats a-dsigmas b-hats b-dsigmas keys))))

(defun plot-ab-1 (a-hats a-dsigmas b-hats b-dsigmas
			 &rest keys
			 &key
			 (style :points)	;Style for a' and b'
			 circles		;Show segments sizes with circles
			 creation-threshold ;Not implemented anymore
			 interpolate
			 show-strongest-cusp
			 &allow-other-keys)
  (when creation-threshold
    (error "Splitting by creation time is no longer implemented in PLOT-AB"))
  (mirror-image-let ((smooth-a-hats (and interpolate (interpolate-hats a-hats :split-count interpolate))))
    (when show-strongest-cusp
      (multiple-value-bind (total-coefficient direction max-coefficient)
	  (total-cusp-gravitational-power-coefficient (find-cusp-parameters a-hats a-dsigmas b-hats b-dsigmas))
	(format t "Strongest cusp ~S, total ~S~%" max-coefficient total-coefficient)
	(setq keys (append (list :mark-locations (list (list "set object circle radius 0.05 lw 2 fc lt 7" ;Red
							     direction)))
			   keys))))
    (apply #'plot-tangent-vectors
	   (append (list (list a-hats "A" style)
			 (list b-hats "B" style))
		   (and interpolate
			(list (list smooth-a-hats nil "lines lt 1")
			      (list smooth-b-hats nil "lines lt 2")))
		   (and circles
			(list (list (hats-circle-data a-hats a-dsigmas) nil "circles lt 1")
			      (list (hats-circle-data b-hats b-dsigmas) nil "circles lt 2"))))
	   :loop t
	   keys))
						       )

;can plot various kink angle trends in a directory over time based on the function its given
;fn options: largest-kink-angle, average-kink-angle, and all-kink-angles
(defun plot-kink-angle-trends (directory fn &rest gnuplot-keys)
  (let* ((kink-angles (kink-angle-trends directory fn))
	 (data nil)
	 (a-angles nil)
	 (b-angles nil))
    (loop for kink in kink-angles
	  do (destructuring-bind
		 (kink-angles timestamp) kink
	       (push (list (first kink-angles) timestamp) a-angles)
	       (push (list (second kink-angles) timestamp) b-angles)))
    (push (coerce (reverse b-angles) 'vector) data)
    (push (coerce (reverse a-angles) 'vector) data)
    (apply #'gnuplot (length data) (max (length a-angles) (length b-angles))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (ecase plot
			  (0 "a-kink-angles")
			  (1 "b-kink-angles"))
		 (let* ((angles (nth plot data)))
		   (when (< point (length angles))
		     (values (nth 1 (aref angles point)) (nth 0 (aref angles point)))))))
	   :style :linespoints
	   gnuplot-keys)))

;Circularity determines whether the speed metric calculated using 2nd derivatives will be scaled relative to that of a
;circle
(defun plot-cusp-2nd-der-trends (directory starting-timestamp Circularity &rest gnuplot-keys)
  (let* ((data (all-cusps-2nd-der-trends directory starting-timestamp Circularity)))
    (setf data (coerce data 'vector))
    (apply #'gnuplot (length data) (loop for cusp across data maximize (length (second cusp)))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (values "")
		 (let* ((cusp-2nd-ders (second (aref data plot))))
		   (when (< point (length cusp-2nd-ders))
		     (values (+ starting-timestamp point) (aref cusp-2nd-ders point))))))
	   :style :linespoints
	   gnuplot-keys)))

(defun plot-b-int-2nd-der-sq-trends (directory &rest gnuplot-keys)
  (let* ((data (b-int-2nd-derivative-squared-trends directory :vector t)))
    (setf data (coerce data 'vector))
    (apply #'gnuplot (loop for timestamp across data maximize (length (first timestamp))) (length data)
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (write-to-string plot)
		 (when (< plot (length (nth 0 (aref data point))))
		   (values (nth 1 (aref data point)) (aref (nth 0 (aref data point)) plot)))))
	   :style :linespoints
	   gnuplot-keys)))

(defun plot-smoothing-by-percent-slow-segments (directory threshold &rest gnuplot-keys)
  (let* ((data (percent-slow-segments directory threshold)))
    (apply #'gnuplot 2 (max (length (aref data 0)) (length (aref data 1)))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (ecase plot
			  (0 "a slow segment %")
			  (1 "b slow segment %"))
		 (let* ((percents (nth plot data)))
		   (when (< point (length percents))
		     (values (nth 1 (aref percents point)) (nth 0 (aref percents point)))))))
	   :style :linespoints
	   gnuplot-keys)))

(defun plot-smoothing-by-int-2nd-der-sq (directory &rest gnuplot-keys)
  (let* ((data (2nd-der-int-data directory)))
    (apply #'gnuplot 2 (max (length (aref data 0)) (length (aref data 1)))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (ecase plot
			  (0 "a integral 2nd derivative squared")
			  (1 "b integral 2nd derivative squared"))
		 (let* ((percents (nth plot data )))
		   (when (< point (length percents))
		     (values (nth 1 (aref percents point)) (nth 0 (aref percents point)))))))
	   :style :linespoints
	   gnuplot-keys)))

(defun plot-smoothing-by-percent-straight-segments (directory &rest gnuplot-keys)
  (let* ((data (percent-straight-segments directory)))
    (apply #'gnuplot 2 (max (length (aref data 0)) (length (aref data 1)))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (ecase plot
			  (0 "a num straight segments")
			  (1 "b num straight segments"))
		 (let* ((percents (nth plot data)))
		   (when (< point (length percents))
		     (values (nth 1 (aref percents point)) (aref (nth 0 (aref percents point)) 0))))))
	   :style :linespoints
	   gnuplot-keys)
    (apply #'gnuplot 2 (max (length (aref data 0)) (length (aref data 1)))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (ecase plot
			  (0 "a percent pieces in straight segments")
			  (1 "b percent pieces in straight segments"))
		 (let* ((percents (nth plot data)))
		   (when (< point (length percents))
		     (values (nth 1 (aref percents point)) (aref (nth 0 (aref percents point)) 1))))))
	   :style :linespoints
	   gnuplot-keys)
    (apply #'gnuplot 2 (max (length (aref data 0)) (length (aref data 1)))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (ecase plot
			  (0 "a sigma in straight segments")
			  (1 "b sigma in straight segments"))
		 (let* ((percents (nth plot data)))
		   (when (< point (length percents))
		     (values (nth 1 (aref percents point)) (aref (nth 0 (aref percents point)) 2))))))
	   :style :linespoints
	   gnuplot-keys)
    ))

(defun plot-all-kink-angles (directory a-or-b &rest gnuplot-keys)
  (let* ((kink-angles (kink-angle-trends directory 'all-kink-angles))
	 (data nil))
    (loop for kinks in kink-angles
	  do (destructuring-bind
	      (kink-angles timestamp) kinks
	      (cond ((string= a-or-b "a")
		     (push (list (coerce (first kink-angles) 'vector) timestamp) data))
		    ((string= a-or-b "b")
		     (push (list (coerce (second kink-angles) 'vector) timestamp) data))
		    (t
		     (break "not 'a' or 'b'")))))
    (setf data (coerce (reverse data) 'vector))
    (apply #'gnuplot (loop for timestamp across data maximize (length (first timestamp))) (length data)
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (values "")
		 (when (< plot (length (first (aref data point))))
		   (values
		    (second (aref data point))
		    (aref (first (aref data point)) plot)))))
	   :style :linespoints
	   gnuplot-keys)))

(defun plot-ab-style (style b-p new-p)
  (ecase style
    (:lines (if b-p (if new-p "lines lt 3" "lines lt 5") ; ;new b blue, old b cyan
	      (if new-p "lines lt 6" "lines lt 8"))) ;new a brown, old a coral
    (:linespoints	;As above, but different point types.  They overlap each other at transitions
     (if b-p (if new-p "linespoints lt 3 pt 2" "linespoints lt 5 pt 4")
       (if new-p "linespoints lt 6 pt 2" "linespoints lt 8 pt 4")))
    (:points '(if b-p "points lt 6 pt 1" "points lt 5 pt 2"))))

;;Returns to arrays with old and new kinks separately.  Segments which connect old and new kinks
;;appear in both.
(defun split-by-creation-time (hats created threshold &key loop)
  (loop with old = (make-array 100 :adjustable t :fill-pointer 0)
	with new = (make-array 100 :adjustable t :fill-pointer 0)
	with this = nil and last = nil
	for index below (length hats)
	when this			;Except first time
	do (vector-push-extend (aref hats index) this) ;Point that comes after previous kink
	do (setq this (if (<= (aref created index) threshold) ;Is next kink early?
			  old new))
	;;If next kink is a different color than previous one, must put this point in old list too
	unless (eq this last)
	do (when last (vector-push-extend nil last)) ;Break present line, except first time
	and do (vector-push-extend (aref hats index) this) ;Switched colors: add present point again
	do (setq last this)
	finally (when loop (vector-push-extend (aref hats 0) this)) ;First point is end of last kink
	(return (values old new))))

;;Split up a' at kinks larger than threshold.
;;If EXCURSION-THRESHOLD is set, 1-point excursions that return to within this threshold will be ignored
(defun plot-a-break (diamond threshold &rest gnuplot-keys &key count skip loop
			     (style :lines)
			     (excursion-threshold (/ threshold 10)) &allow-other-keys)
  (multiple-value-bind (hats sigmas) (get-a-data-sigmas diamond :count count :skip skip)
    (let ((n (length hats))
	  (data nil)
	  (joins (make-array 100 :adjustable t :fill-pointer 0))) ;One array for joining segments
      (flet ((add-piece (start end)	;Add piece from start through end inclusive
	       (let ((piece (make-array (1+ (mod (- end start) n)))))
		 (loop for j = start then (mod (1+ j) n) ;Copy pieces from start through end
		       for k from 0
		       do (setf (aref piece k) (aref hats j))
		       until (= j end))
		 (push (list piece (format nil "Length ~7F, count ~D"
					   (mod (- (aref sigmas end) (aref sigmas start)) (aref sigmas (1- n)))
					   (1+ (mod (- end start) n)))
				   style)
		       data))))
	(loop with looped = nil
	      with excursion = nil
	      with start = (unless loop -1) ;Position where last big kink began.  If loop, do first segment last
	      with index = 0
	      do
	      (cond (excursion		;Continue through excursion
		     (setq excursion nil))
		    ((< (spherical-angle (aref hats index) (aref hats (mod (1+ index) n))) threshold) ;Small kink?
		     )			;Just go on
		    ((< (spherical-angle (aref hats index) (aref hats (mod (+ 2 index) n))) excursion-threshold)
		     (setq excursion t)) ;Large kink, but excursion
		    (t			;Large kink: break here
		     (when start
		       (add-piece (1+ start) index)
		       (vector-push-extend (aref hats index) joins)
		       (vector-push-extend (aref hats (mod (1+ index) n)) joins)
		       (vector-push-extend nil joins))
		     (when looped (return nil)) ;Already looped: return now
		     (setq start index))) ;Start new segment
	      when (= (incf index) n)	;Onto next.  Check for end of array
	      do (cond (looped		;Looped twice, so there are no big kinks
			(push (list hats (format nil "Length ~7F, count ~D" (aref sigmas (1- n)) n) style)
			      data)
			(return nil))
		       (loop (setq index 0 looped t)) ;Back to finish last segment
		       (t (unless (= start (1- n))
			    (add-piece (1+ start) (1- n)))	;Do last piece
			  (return nil)))))
      (apply #'plot-tangent-vectors (cons (list joins nil "lines lt 0") (nreverse data)) gnuplot-keys))))
	       
(defun plot-all-ab (&key (style :points) rest-frame (rotate (not rest-frame)) (b t) (initial-time 4.5))
  (let ((loops nil))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
     #'(lambda (diamond)
	 (let ((length (loop-length-and-count-a diamond)))
	   (when (and (>= (diamond-countup diamond) 2) ;Non-self-intersecting loop
		      (< length *total-size*)) ;No wrapping loops
	     (push (list diamond length) loops)))))
    (setq loops (coerce loops 'vector))
    (setq loops (sort loops #'> :key #'second))
    (loop with index = 0
	  for (diamond len) = (aref loops index)
	  do (if b (plot-ab (if rest-frame (rest-frame-loop diamond) diamond)
			    :title (format nil "Length ~7F, count ~D" len (loop-count diamond))
			    :reuse t :style style :rotate rotate
			    :initial-time initial-time)
	       (plot-a (if rest-frame (rest-frame-loop diamond) diamond)
			 :title (format nil "Length ~7F, count ~D" len (loop-count diamond))
			 :reuse t :style style :rotate rotate
			 :initial-time initial-time))
	  do (let ((command (read-char)))
	       (cond ((char-equal command #\Newline) (incf index))
		     ((char= command #\b) (decf index))
		     ((char= command #\l) (setq style :lines))
		     ((char-equal command #\p) (setq style :points))
		     ((char= command #\L) (setq style :linespoints))
		     ((char-equal command #\d) (setq style :dots))
		     ((char-equal command #\A) (setq b nil))
		     ((char-equal command #\B) (setq b t))
		     ((char-equal command #\q) (return nil))
		     ((char= command #\() ;Evaluate lisp form
		      (unread-char command)
		      (return nil)))
	       (unless (char-equal command #\Newline) (read-line)))))) ;Discard rest of line

(mirror-images
;;Find segment near a given position.  Returns a-hat and length
(defun find-a-hat (diamond vector &key count skip)
  (multiple-value-bind (hats sigmas) (get-a-data-sigmas diamond :count count :skip skip)
    (loop with best-index = nil
	  with best-angle
	  for index from 0 below (length hats)
	  for this = (aref hats index)
	  for this-angle = (spherical-angle vector this)
	  when (or (null best-index)	;See if this one is better
		   (< this-angle best-angle))
	  do (setq best-index index best-angle this-angle)
	  finally (format t "~&Best angle ~8E at index ~D.~%" best-angle best-index)
	  (return (values (aref hats best-index)
			  (- (aref sigmas best-index)
			     (if (zerop best-index) 0.0 (aref sigmas (1- best-index)))))))))
)
		   

;;Clustering analysis

(defstruct (cluster
	    (:print-object print-cluster))
  (left nil)				;Left sub-cluster or NIL
  (right nil)
  sigma					;Amount of sigma covered
  vector				;Unit vector giving average direction of elements in this cluster
  (max-angle 0.0 :type double-float))	;Maximum angle between left and right any node in our subtree

(defun print-cluster (cluster stream)
  (print-unreadable-object (cluster stream :identity t :type t)
    (format stream "sigma ~F direction ~S" (cluster-sigma cluster) (cluster-vector cluster))))

;;Make a cluster out of two previous clusters.  The direction is the weighted average direction.
;;ANGLE is the angle between the vectors of the two clusters that are being merged
(defun merge-clusters (left right angle)
  (make-cluster :left left :right right
		:sigma (+ (cluster-sigma left) (cluster-sigma right))
		:vector (3vector-normalize (3vector+ (3vector-scale (cluster-vector left) (cluster-sigma left))
						     (3vector-scale (cluster-vector right) (cluster-sigma right))))
		:max-angle (max (cluster-max-angle left) (cluster-max-angle right) angle)))

;;To do the clustering, we keep a priority queue of the angles between clusters and always merge the
;;smallest angle.  Each angle is represented by the following data structure.
;;For a loop, each angle has an element in the queue.  For an open segment, there is also 
;;one of these that is not queued, for storing the leftmost cluster.
(defstruct cluster-queue-element
  cluster				;Cluster object at right of angle
  previous				;Previous cluster-queue-element
  next)					;Next cluster-queue-element

;;Return angle associated with this element
(defun cluster-queue-angle (cluster)
  (let ((cos (3vector-dot (cluster-vector (cluster-queue-element-cluster cluster))	;This is the object at the right
			  (cluster-vector (cluster-queue-element-cluster (cluster-queue-element-previous cluster)))))) ;And this the left
    (if (> cos 1) 0.0			;Avoid complex numbers due to numerical error
      (acos cos))))

(defvar *cluster-queue*)

(defun cluster-queue-add (element)
  (calendar-add  *cluster-queue* (cluster-queue-angle element) element))
(defun cluster-queue-delete (element)
  (calendar-delete *cluster-queue* (cluster-queue-angle element) element))
;;Return element to be merged and its angle
(defun cluster-queue-next ()
  (multiple-value-bind (angle element) (calendar-next *cluster-queue*)
    (values element angle)))

;;Make clusters of segments
;;We take the smallest distance and combine the elements
(mirror-images
(defun cluster-segments-a (diamond &key count skip loop)
  (multiple-value-bind (hats sigmas) (get-a-data-sigmas diamond :count count :skip skip)
    (let ((*cluster-queue* (make-calendar 1.0 :name "clustering"
					  :backups-allowed t))) ;Sometimes merging makes a cluster closer to another
      (loop with leftmost
	    for index from 0 below (length hats)
	    for cluster = (make-cluster :vector (aref hats index)
					:sigma (- (aref sigmas index)
						  ;;Subtract last sigma or zero for first.  Avoid warning that
						  ;;in the initialization copy of this code, index is known to be 0.
						  (if (zerop index) 0 (without-compiler-notes
								       (aref sigmas (1- index))))))
	    for previous = nil then queue-element
	    for queue-element = (make-cluster-queue-element :cluster cluster :previous previous)
	    if previous			;Except first time
	      do (setf (cluster-queue-element-next previous) queue-element) ;Connect previous element to this one
	         (cluster-queue-add queue-element) ;Put into queue
	    else do (setq leftmost queue-element) ;Otherwise remember it in case of loop
	    until (and count (zerop (decf count))) ;Limit to given number
	    finally
	    (when loop
	      (setf (cluster-queue-element-previous leftmost) queue-element ;Connect first element to last one
		    (cluster-queue-element-next queue-element) leftmost)
	      (cluster-queue-add leftmost))) ;Put first element in data structure
      (make-clusters)))))

(defun make-clusters ()
  (loop for (element angle) = (multiple-value-list (cluster-queue-next)) ;Get smallest angle, delete from queue
	for previous = (cluster-queue-element-previous element) ;Previous angle or first-cluster placeholder
	for next = (cluster-queue-element-next element) ;Next angle or NIL
	for right-cluster = (cluster-queue-element-cluster element) ;We store cluster to our right
	for left-cluster = (cluster-queue-element-cluster previous) ;He has cluster to our left
	for new-cluster = (merge-clusters left-cluster right-cluster angle) ;Make a new cluster of the two old
	do
	;;Remove angles on both sides of us from the queue, because we are modifying them.
	(when (cluster-queue-element-previous previous) ;If previous is a real angle
	  (cluster-queue-delete previous))	;Remove from queue
	(when (and next			;Have next?
		   (not (eq next previous))) ;Last time these are equal
	  (cluster-queue-delete next))	;Remove from queue
	(setf (cluster-queue-element-cluster previous) new-cluster) ;New cluster goes on right of previous element
	(setf (cluster-queue-element-next previous) next) ;Splice removed element out of thread
	(when next
	  (setf (cluster-queue-element-previous next) previous))
	(when (cluster-queue-element-previous previous) ;If previous is a real angle rather than placeholder
	  (cluster-queue-add previous))	;Put back in queue
	(when (eq next previous)	;Current element removed, next = previous, so done
	  (return (cluster-queue-element-cluster previous)))
	(when next			;Have next?
	  (cluster-queue-add next))
	;;Exit in case of loop above.  If not loop, exit when queue empty and only placeholder left.
	until (calendar-empty-p *cluster-queue*)
	finally (return (cluster-queue-element-cluster previous))))	;Return overall cluster

;;Split clusters whose elements are further away than the threshold and plot them separately
(defun plot-clusters (top threshold &rest keys &key loop &allow-other-keys)
  (let ((data (make-array 0 :fill-pointer 0 :adjustable t))
	(plot-type 0))			;First type is 1
    (walk-clusters
     #'(lambda (vectors sigma)
	 (vector-push-extend (list vectors (format nil "~D pt~:P, sigma ~,3,1G" (length vectors) sigma)
				   (format nil "lp lt ~D pt ~:*~D" (incf plot-type)))
			     data))
     top threshold)
    (setq data (add-connecting-segments data nil "lines lt 0" :loop loop)) ;Dotted lines between clusters
    (apply #'plot-tangent-vectors data :title (format nil "Total length ~,3,1G" (cluster-sigma top))
	   keys)))

(defun histogram-clusters (top threshold &rest gnuplot-keys
			       &key (bins 20) (min 1.0) (max *total-size*) &allow-other-keys)
  (setq max (float max 0.0))
  (let ((bin-size (/ (log (/ max min)) bins)) ;Logarithmic size of a bin
	(counts (make-array bins :element-type 'fixnum :initial-element 0)) ;Number of segments
	(total-l2 0.0)
	(total-l 0.0)
	)
    (walk-clusters
     #'(lambda (vectors sigma)
	 (declare (ignore vectors))
	 (incf total-l2 (* sigma sigma)) ;Length-weighted average length
	 (incf total-l sigma)
	 (when (> sigma min)
	   (let ((index (fixnum-floor (/ (real-log (/ sigma min)) bin-size))))
	     (when (< index bins)
	       (incf (aref counts index))))))
     top threshold)
    (apply #'gnuplot 1 bins
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (unless (eq point :title)
		 (let ((length (* min (exp (* bin-size (+ point 0.5))))) ;Central length value
		       (count (aref counts point))) ;Number of segments in this bin
		   (and (plusp count) ;no point, rather than 0 if no data
			(values-list (list length count))))))
	   :styles :linespoints
	   :logscale '(:x :y)
	   gnuplot-keys)
    (format t "~&Length-weighted average length ~S" (/ total-l2 total-l))
    counts))

;;Compute the length-weighted average of the "inverse wiggliness parameter",
;;|a(sigma_f) - a(sigma_i)|/(sigma_f - sigma_i)
(defun cluster-wiggliness (top threshold)
  (let ((total-l 0.0)			;Total of sigma
	(total-distance 0.0)		;Total of a(sigma_f) - a(sigma_i)
	(total-l2d 0.0))		;Total of sigma^2/d
    (walk-clusters
     #'(lambda (vectors length)
	 (unless (zerop length)
	   (incf total-l length)
	   (loop with position = zero-3vector
		 for index below (length vectors)
		 for p-hat = (first (aref vectors index))
		 for dsigma = (- (if (= (1+ index) (length vectors)) ;Values returned for SIGMA here are awkward
				     (+ (second (aref vectors 0)) length)
				   (second (aref vectors (1+ index)))) ;sigma at end of this segment
				 (second (aref vectors index))) ;sigma at beginning of this segment
		 do (setq position (3vector+ position (3vector-scale p-hat dsigma)))
		 finally (let ((distance (3vector-length position)))
			   (incf total-distance distance)
			   (incf total-l2d (* length (/ length distance)))))))
     top threshold)
    (values (/ total-distance total-l)
	    (/ total-l2d total-l))))
	 
#|
      (defun plot-clusters-2 (top-a top-b threshold &rest keys)
(mirror-image-let ((data-a (make-array 0 :fill-pointer 0 :adjustable t)))
(mirror-images (walk-clusters top-a threshold data-a))
(apply #'plot-steps-1 (list data-a data-b) keys)))
      |#

(defvar *walk-cluster-sigma*)

;;Recursively walk over the high-level clusters to find ones that are close enough to keep together
;;Call FUNCTION with array of (P-HAT SIGMA) and the total SIGMA in the cluster
(defun walk-clusters (function cluster threshold)
  (let ((*walk-cluster-sigma* 0.0))
    (walk-clusters-1 function cluster threshold)))

(defun walk-clusters-1 (function cluster threshold)
  (mirror-image-let ((left (cluster-left cluster)))
    (cond ((or (and (null left) (null right)) ;Terminal cluster encountered when splitting: 1 point in this cluster
	       (< (cluster-max-angle cluster) threshold)) ;Or all subclusters close together
	   ;;Found a cluster to output.  Make an array to store it
	   (let ((results (make-array 0 :fill-pointer 0 :adjustable t))) ;Make array for results
	     (walk-clusters-2 cluster results)	;fill it up
	     (funcall function results (cluster-sigma cluster)))) ;Call function
	  (t				;Far apart.  They will be separate clusters
	   (walk-clusters-1 function left threshold)
	   (walk-clusters-1 function right threshold)))))

;;Walk over cluster that has been determined to belong together and collect leaves.
(defun walk-clusters-2 (cluster results)
  (mirror-image-let ((left (cluster-left cluster)))
    (cond ((and left right)		;Children?
	   (walk-clusters-2 left results) ;Walk them
	   (walk-clusters-2 right results))
	  ((or left right)		;Trouble
	   (error "~S should have had either two children or none" cluster))
	  (t				;Terminal cluster: collect it
	   (vector-push-extend (list (cluster-vector cluster) *walk-cluster-sigma*) results)
	   (incf *walk-cluster-sigma* (cluster-sigma cluster))))))

;;Split at given angle and count number of segments
(defun count-clusters (cluster threshold)
  (mirror-image-let ((left (cluster-left cluster)))
    (cond ((and (null left) (null right)) ;Terminal cluster encountered when splitting: 1 point in this cluster
	   0)				;Don't count at all
	  ((< (acos (3vector-dot (cluster-vector left) (cluster-vector right))) threshold) ;close together
	   1)
	  (t
	   (+ (count-clusters left threshold)
	      (count-clusters right threshold))))))

(defun cluster-distance (cluster threshold)
  (/ (cluster-sigma cluster) (count-clusters cluster threshold)))
 
(defun long-string-cluster-distance (threshold)
  (let ((total-length 0.0)
	(count 0))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Only consider loops
     #'(lambda (diamond)
	 (let ((length (loop-length-and-count-a diamond)))
	   (when (> length (* *total-size* 3)) ;Only long loops
	     (incf total-length length)
	     (incf count (count-clusters (cluster-segments-a diamond) threshold))))))
    (values (/ total-length count) total-length)))

(defun make-range-near (position distance)
  (let ((displacement (make-3vector distance distance distance)))
    (list (3vector- position displacement) (3vector+ position displacement))))
