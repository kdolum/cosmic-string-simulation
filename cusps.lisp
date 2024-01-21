;;;Code for analyzing cusps (and some for kinks)
;;;Mostly by Vishnu Gade
(in-package "CL-USER")

;;Find kinks that cross each other on a loop, forming "cusps"
;;Returns info about the cusp.
(defun find-crossings-for-cusps (diamond &key count skip )
  (multiple-value-bind (a-hats a-dsigmas) (get-a-data-dsigmas diamond :count count :skip skip)
    (multiple-value-bind (b-hats b-dsigmas) (get-b-data-dsigmas diamond :count count :skip skip)
      (find-crossings-for-cusps-1 a-hats a-dsigmas b-hats b-dsigmas))))

;;Find crossings that may be cusps from hats and dsigmas
(defun find-crossings-for-cusps-1 (a-hats a-dsigmas b-hats b-dsigmas)
  (mirror-image-let ((n-a (length a-hats)))
    ;;Brute force
    (loop for i-a below n-a
	  for a-1 = (aref a-hats i-a)
	  for a-2 = (aref a-hats (mod (1+ i-a) n-a))
	  nconc (loop for i-b below n-b
		      for b-1 = (aref b-hats i-b)
		      for b-2 = (aref b-hats (mod (1+ i-b) n-b))
		      for intersection = (spherical-intersection a-1 a-2 b-1 b-2)
		      when intersection
		      collect (mirror-images-make-cusp-info
			       :direction intersection
			       :i-a i-a :a-1 a-1 :a-2 a-2 :a-3 (aref a-hats (mod (+ i-a 2) n-a))
			       :a-dsigma-1 (aref a-dsigmas i-a)
			       :a-dsigma-2 (aref a-dsigmas (mod (1+ i-a) n-a))
			       :a-dsigma-3 (aref a-dsigmas (mod (+ i-a 2) n-a)))))))
				   
;;Make a loop kind of like a 1,1 Burden loop, but a' goes around its circle unevenly, so a''' at the
;;cusp is not trivial
(defun setup-cusp-form-paper-check (factor points)
    (make-test-ab-prime
     #'(lambda (sigma)
	 (let ((alpha (- sigma (* factor (cos (* 2 sigma))))))
	   (make-3vector (cos alpha) 0.0 (sin alpha))))
     #'(lambda (sigma)
	 (make-3vector 0.0 (cos sigma) (sin sigma)))
    (* 2 pi)
    points))

;;Returns lists of a'', b'', a''', b'''
(defun check-cusp-form-paper (diamond)
  (declare (optimize debug))
  (let*((cusps (find-crossings-for-cusps diamond)) ;We find the crossings for a and b
	(result nil))
    (flet ((d3 (a1 a2 a3 ds12 ds23)
	       (3vector-scale
		(3vector-
		 (3vector-scale (3vector- a3 a2) (/ ds23))  ;first deriv 2->3
		 (3vector-scale (3vector- a2 a1) (/ ds12))) ;first deriv 1->2
		(/ 2 (+ ds12 ds23)))))			    ;divide by dist from 1&2 to 2&3
      (loop for (a1 a2 b1 b2 a-delta-sigma b-delta-sigma cusp-position i j a0 a3 b0 b3
		    a-dsigma-prev a-dsigma-next b-dsigma-prev b-dsigma-next)
	    in cusps
	    do
	    (let*((total-a-angle (spherical-angle a1 a2)) ;We compute the total spherical angle betwen a1 and a2
		  (total-b-angle (spherical-angle b1 b2))		     		   
		  (angle-cusp-a1 (spherical-angle a1 cusp-position)) ;Compute the angle between the cusp and a1     
		  (angle-cusp-b1 (spherical-angle b1 cusp-position))
		  (sigma-a1-to-cusp (* a-delta-sigma (/ angle-cusp-a1 total-a-angle))) ;Compute the amount of sigma proportial to the angular distance
		  (sigma-cusp-to-a2 (- a-delta-sigma sigma-a1-to-cusp))
		  (sigma-b1-to-cusp (* b-delta-sigma (/ angle-cusp-b1 total-b-angle)))
		  (sigma-cusp-to-b2 (- b-delta-sigma sigma-b1-to-cusp))
		  ;;We evaluate the second derivative by taking the average before and after the cusp
		  (a-second-derivative	;For a we go the wrong way, but for second derivative it's OK
		   (3vector-scale (3vector+ 
				   (3vector-scale (3vector- a2 cusp-position) (/ 1.0 sigma-cusp-to-a2))
				   (3vector-scale (3vector- cusp-position a1) (/ 1.0 sigma-a1-to-cusp)))
				  0.5))
		  (b-second-derivative
		   (3vector-scale (3vector+ 
				   (3vector-scale (3vector- b2 cusp-position) (/ 1.0 sigma-cusp-to-b2))
				   (3vector-scale (3vector- cusp-position b1) (/ 1.0 sigma-b1-to-cusp)))
				  0.5))						  ; The same thing for b
		  (a-third-derivative (d3 a1 a2 a3 a-delta-sigma a-dsigma-next))  ;not centered
		  (b-third-derivative (d3 b1 b2 b3 b-delta-sigma b-dsigma-next))) ;not centered
	      (push (list a-second-derivative b-second-derivative a-third-derivative b-third-derivative) result)))
      (reverse result))))		;Same order as find-crossings-for-cusps


;;Gade
(defun find-cusp-position-in-hats-with-hats (a-hats b-hats)
  (let*((cusp-positions nil))
    (loop for i below (length a-hats)
	  nconc (loop with a1 = (aref a-hats i)
		      with a2 = (aref a-hats (mod (1+ i) (length a-hats)))
		      for j below (length b-hats)
		      when (spherical-intersection a1 a2 (aref b-hats j)
						   (aref b-hats (mod (1+ j) (length b-hats))))
		      do (push (list i j) cusp-positions)))
    cusp-positions))

#| I don't understand what this plotted.  It could be revived by using find-cusp-parameters instead of 
   recalculating things.  -- kdo
(defun plot-second-derivatives (diamond &rest gnuplot-keys)
  (let ((cusp-positions (find-cusp-position-in-hats diamond))
	(a-hats (get-a-data-sigmas diamond))
        (b-hats (get-b-data-sigmas diamond)))
      (apply #'gnuplot 4 (length a-hats)
             #'(lambda (plot point)
                 (if (eq point :title) (nth plot '("A''" "B''" "A-CUSP" "B-CUSP"))
		   (if (= plot 2)
		       (if (< point (length cusp-positions))
			   (values (first (nth point cusp-positions)) 
				   (spherical-angle (aref a-hats (first (nth point cusp-positions)))
						    (aref a-hats (mod (1+ (first (nth point cusp-positions))) (length a-hats)))))
			 nil)
		   (if (= plot 3)
		       (if (< point (length cusp-positions))
			   (values  (second (nth point cusp-positions)) 
				    (spherical-angle (aref b-hats (second (nth point cusp-positions)))
						     (aref b-hats (mod (1+ (second (nth point cusp-positions))) (length b-hats)))))
			 nil)
		     (let* ((hats (if (zerop plot) a-hats b-hats))
			    (first (aref hats point))
			    (second (aref hats (mod (1+ point) (length hats))))
			    (angle (spherical-angle first second)))
		       (when (> angle 0.1)
			 (let ((*print-pretty* nil))
			   (format t "~&Large angle at position ~D in ~A: ~S to ~S is ~S~%"
				   point (nth plot '("A'" "B'"))
				   first second angle)))
                     (values point angle))))))
             gnuplot-keys))) |#

;; Returns the integral of hats using the sigmas, i.e., get A from A'
(defun integrate-hats (hats sigmas)
  (let*((count (length sigmas))
	(a (make-array count)))
    (setf (aref a 0) zero-3vector)
    (loop for index below count
	  for previous = (previous-index-wrapping index count)
	  for next = (next-index-wrapping index count)
	  for delta-sigma = (aref sigmas index) then (- (aref sigmas index) (aref sigmas previous))
	  do (setf (aref a next) (3vector+ (aref a index) (3vector-scale (aref hats index) delta-sigma))))
    (values a)))



;;Locate cusps and and annotate them with additional information.
(defun find-cusp-parameters (a-hats a-sigmas b-hats b-sigmas)
  (let ((cusps (find-crossings-for-cusps-1 a-hats a-sigmas b-hats b-sigmas))
	(l (average-total-sigma a-sigmas b-sigmas))) ;length of loop
    ;;Fill in more slots
    (loop for cusp in cusps
	  do (mirror-images
	      ;;We imagine that the time a' takes to travel over the distance between a1 and a2 is the
	      ;;average of the dwell time at the two ends.
	      (let* ((avg-a-dsigma (/ (+ (cusp-info-a-dsigma-1 cusp) (cusp-info-a-dsigma-2 cusp)) 2))
		     (total-a-angle (spherical-angle (cusp-info-a-1 cusp) (cusp-info-a-2 cusp)))
		     ;;Compute the amount of sigma proportial to the angular distance
		     (angle-cusp-a-1 (spherical-angle (cusp-info-a-1 cusp) (cusp-info-direction cusp)))
		     (sigma-a-1-to-cusp (* avg-a-dsigma (/ angle-cusp-a-1 total-a-angle)))
		     ;;Cross product gives the axis around which a' is rotating
		     (axis (3vector-normalize (3vector-cross-product (cusp-info-a-1 cusp) (cusp-info-a-2 cusp))))
		     ;;Second derivative is the tangent vector to the path of rotation at the cusp
		     (a-pp (3vector-scale (3vector-cross-product axis (cusp-info-direction cusp))
					  (* l (/ total-a-angle avg-a-dsigma))))) ;L a''
		(setf (cusp-info-sigma-a-1-to-cusp cusp) sigma-a-1-to-cusp
		      (cusp-info-sigma-cusp-to-a-2 cusp) (- avg-a-dsigma sigma-a-1-to-cusp)
		      (cusp-info-a-pp cusp) a-pp))))
    cusps))

;A wrapper function that takes a .dat filename, finds the a and b hats and sigmas and then locates the cusps using
;find-cusp-parameters-in-sigmas
(defun read-and-find-cusps (filename)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)(read-hats-dsigmas filename)
    (find-cusp-parameters a-hats a-dsigmas b-hats b-dsigmas)))

;Calculates and outputs the distance between 2 cusps in a-sigma and b-sigma
(defun distance-between-cusps-in-sigma (cusp-1 cusp-2)
  (let* ((output nil)
	 (a-dif 0)
	 (b-dif 0))
    (destructuring-bind
	(cusp-1-a cusp-1-a-fraction cusp-1-b cusp-1-b-fraction cusp-1-a1-to-cusp cusp-1-cusp-to-a2
		  cusp-1-b1-to-cusp cusp-1-cusp-to-b2 &rest ignore)
	cusp-1
      (declare (ignore ignore))
      (destructuring-bind
	  (cusp-2-a cusp-2-a-fraction cusp-2-b cusp-2-b-fraction cusp-2-a1-to-cusp cusp-2-cusp-to-a2
		   cusp-2-b1-to-cusp cusp-2-cusp-to-b2 &rest ignore)
	  cusp-2
	(declare (ignore ignore))
	(mirror-images
	 (cond ((eql cusp-1-a cusp-2-a)
					;if the two cusps are located in the same segment of string just take the
					;difference in their location in the segment
		(setf a-dif
		      (* (abs (- cusp-2-a-fraction cusp-1-a-fraction))
			 (+ cusp-2-a1-to-cusp cusp-2-cusp-to-a2))))
	       (t
					;if not sum the sigma from the first cusp to the hat and then from the hat to the
					;second cusp
		(cond ((> cusp-2-a cusp-1-a)
		       (setf a-dif
			     (+ cusp-2-a1-to-cusp cusp-1-cusp-to-a2)))
		      ((< cusp-2-a cusp-1-a)
		       (setf a-dif
			     (+ cusp-2-cusp-to-a2 cusp-1-a1-to-cusp))))))
	 (push a-dif output))))
    output))

;removes from the list of cusps any cusps with kink angles greater than 45 degrees and 2nd derivatives of a or b greater
;than 10
(defun remove-fake-cusps (list-of-cusps)
  (let* ((output nil))
    (loop for cusp in list-of-cusps
	  do (destructuring-bind
		 (cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
			 cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle
			 cusp-a-second-derivative cusp-b-second-derivative intersection)
		 cusp
	       (declare (ignore cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
				cusp-b1-to-cusp cusp-cusp-to-b2 intersection))
	       (cond ((and
		       (< cusp-a-second-derivative 10) (< cusp-b-second-derivative 10)
		       (< cusp-a-angle (/ pi 4)) (< cusp-b-angle (/ pi 4)))
		      (push cusp output)))))
    output))

;loops through all timestamps of a string in a directory and looks for future possible repeat cusps for each cusp that exists at
;the inputed
;initial timestamp. If there are multiple possible repeat cusps within the same segment the code will pick the closer one to the
;initial cusp.
(defun num-repeat-cusps (initial-timestamp directory)
  (let*((num-of-repeat-cusps 0)
	(repeat-cusps nil)
	(*backreaction-filename-digits* (find-backreaction-digits directory))
	(cusps-first-timestamp
	 (remove-fake-cusps (read-and-find-cusps (backreaction-dump-filename directory initial-timestamp)))))
    (multiple-value-bind
     (a-hats) (read-hats-dsigmas (backreaction-dump-filename directory 1))
     (loop for timestamp from (1+ initial-timestamp)
					;loop through all following timestamps starting from the inputed starting timestamp
	   while (probe-file (backreaction-dump-filename directory timestamp))
	   do (let* ((cusps-current-timestamp
		      (remove-fake-cusps (read-and-find-cusps (backreaction-dump-filename directory timestamp)))))
		(loop for cusp in cusps-first-timestamp
		      do (let* ((possible-repeats nil)
					;create list of possible repeat cusps for the cusp from the starting timestamp that the
					;loop is on 
				(best-cusp nil))
			   (destructuring-bind
			       (cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
				       cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle
				       cusp-a-second-derivative cusp-b-second-derivative)
			       cusp
			     (declare (ignore cusp-a1-to-cusp cusp-cusp-to-a2 cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle))
			     (loop for curr-cusp in cusps-current-timestamp
				   do (destructuring-bind
					  (curr-cusp-a curr-cusp-a-fraction curr-cusp-b curr-cusp-b-fraction &rest ignore)
					  curr-cusp 
					(declare (ignore curr-cusp-a-fraction curr-cusp-b-fraction ignore))
					(cond ((and (<= (abs (mod (- cusp-a curr-cusp-a) (length a-hats))) 1)
						    (<= (abs (mod (- cusp-b curr-cusp-b) (length a-hats))) 1))
					       (push curr-cusp possible-repeats)))))
					;if a cusp from the timestamp the code is on is within 1 a and b hat of the cusp from
					;the starting
			                ;timestamp add it to the list of possible repeats
			    (cond ((not (eql possible-repeats nil));given possible repeats is not empty
				   (setf best-cusp (first possible-repeats))
				   (loop for curr-cusp in possible-repeats
					 do
					 #| (destructuring-bind
					 (curr-cusp-a curr-cusp-a-fraction curr-cusp-b curr-cusp-b-fraction ;
					 curr-cusp-a1-to-cusp curr-cusp-cusp-to-a2 ;
					 curr-cusp-b1-to-cusp curr-cusp-cusp-to-b2 ;
					 curr-cusp-a-angle curr-cusp-b-angle ;
					 curr-cusp-a-second-derivative curr-cusp-b-second-derivative) ;
					 curr-cusp ;
					 (declare (ignore curr-cusp-a curr-cusp-a-fraction ;
					 curr-cusp-b curr-cusp-b-fraction ;
					;
					 (destructuring-bind ;
					 (best-cusp-a best-cusp-a-fraction best-cusp-b best-cusp-b-fraction ;
					 best-cusp-a1-to-cusp best-cusp-cusp-to-a2 ;
					 best-cusp-b1-to-cusp best-cusp-cusp-to-b2 ;
					 best-cusp-a-angle best-cusp-b-angle ;
					 best-cusp-a-second-derivative best-cusp-b-second-derivative) ;
					 best-cusp |#
					 (let* ((best-cusp-a-dif
						 (first (distance-between-cusps-in-sigma cusp best-cusp)))
						(best-cusp-b-dif
						 (second (distance-between-cusps-in-sigma cusp best-cusp)))
						(curr-cusp-a-dif
						 (first (distance-between-cusps-in-sigma cusp curr-cusp)))
						(curr-cusp-b-dif
						 (second (distance-between-cusps-in-sigma cusp curr-cusp))))
					   (cond ((< (sqrt (+ (expt curr-cusp-a-dif 2) (expt curr-cusp-b-dif 2)))
						     (sqrt (+ (expt best-cusp-a-dif 2) (expt best-cusp-b-dif 2))))
						  (setf best-cusp curr-cusp)))))
					;if the current cusp is closer to the cusp from the starting timestamp set it as the new
					;best cusp
				    (destructuring-bind
					(best-cusp-a best-cusp-a-fraction best-cusp-b best-cusp-b-fraction
						     best-cusp-a1-to-cusp best-cusp-cusp-to-a2
						     best-cusp-b1-to-cusp best-cusp-cusp-to-b2
						     best-cusp-a-angle best-cusp-b-angle
						     best-cusp-a-second-derivative best-cusp-b-second-derivative)
					best-cusp
				      (declare (ignore best-cusp-a1-to-cusp best-cusp-cusp-to-a2
						       best-cusp-b1-to-cusp best-cusp-cusp-to-b2
						       best-cusp-a-angle best-cusp-b-angle))
				      (push (list (list cusp-a cusp-a-fraction cusp-b cusp-b-fraction
							cusp-a-second-derivative cusp-b-second-derivative)
						  (list best-cusp-a best-cusp-a-fraction best-cusp-b best-cusp-b-fraction
							best-cusp-a-second-derivative best-cusp-b-second-derivative)
						  timestamp)
					    ;pushes the cusp from the starting timestamp 
					    repeat-cusps))
				    (incf num-of-repeat-cusps)
				    ))))))))
    (list num-of-repeat-cusps (reverse repeat-cusps) cusps-first-timestamp)))

;gets the data for all possible cusps at a certain timestamp
(defun possible-cusps-data (directory timestamp threshold)
  (let* ((cusps nil)
	 (*backreaction-filename-digits* (find-backreaction-digits directory))
	 (possible-cusps (read-and-find-cusps (backreaction-dump-filename directory timestamp)))
	 (num-2nd-der-under-thresh 0))
    (loop for cusp in possible-cusps
	  do (destructuring-bind
		 (cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
			 cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle
			 cusp-a-second-derivative cusp-b-second-derivative intersection)
		 cusp
	       (declare (ignore cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
			 cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle intersection))
	       (cond ((and
		       (< cusp-a-second-derivative threshold)
		       (< cusp-b-second-derivative threshold))
		      (setf num-2nd-der-under-thresh (1+ num-2nd-der-under-thresh))))
	       (push cusp cusps)))
    (list num-2nd-der-under-thresh cusps)
    ))

;finds a cusp in cusps-curr-timestamp that is most likely to be the future of the cusp given by the first input
(defun matching-cusp (cusp cusps-curr-timestamp length-a-hats length-b-hats)
  (let* ((possible-matches nil)
	 (best-cusp nil))
    (destructuring-bind
	(cusp-a cusp-a-fraction cusp-b &rest ignore)
	cusp
      (declare (ignore cusp-a-fraction ignore))
      (loop for curr-cusp in cusps-curr-timestamp
	    do (destructuring-bind
		   (curr-cusp-a curr-cusp-a-fraction curr-cusp-b &rest ignore)
		   curr-cusp
		 (declare (ignore curr-cusp-a-fraction ignore))
		 (cond ((and (<= (abs (mod (- cusp-a curr-cusp-a) length-a-hats)) 1)
			     (<= (abs (mod (- cusp-b curr-cusp-b) length-b-hats)) 1))
			(push curr-cusp possible-matches)))))
					;if a cusp from the timestamp the code is on is within 1 a and b hat of the cusp
					;from the starting timestamp add it to the list of possible repeats
      (setf possible-matches (reverse possible-matches))
      (setf best-cusp (first possible-matches))
      (loop for curr-cusp in possible-matches
	    do #| (destructuring-bind
		   (curr-cusp-a curr-cusp-a-fraction curr-cusp-b curr-cusp-b-fraction
				curr-cusp-a1-to-cusp curr-cusp-cusp-to-a2 curr-cusp-b1-to-cusp curr-cusp-cusp-to-b2
				curr-cusp-a-angle curr-cusp-b-angle
				curr-cusp-a-second-derivative curr-cusp-b-second-derivative intersection)
		   curr-cusp
		 (destructuring-bind
		     (best-cusp-a best-cusp-a-fraction best-cusp-b best-cusp-b-fraction
				  best-cusp-a1-to-cusp best-cusp-cusp-to-a2 best-cusp-b1-to-cusp best-cusp-cusp-to-b2
				  best-cusp-a-angle best-cusp-b-angle
				  best-cusp-a-second-derivative best-cusp-b-second-derivative intersection)
		     best-cusp |#
		   (let* ((best-cusp-difs (distance-between-cusps-in-sigma cusp best-cusp))
			  (best-cusp-a-dif (first best-cusp-difs))
			  (best-cusp-b-dif (second best-cusp-difs))
			  (curr-cusp-difs (distance-between-cusps-in-sigma cusp curr-cusp))
			  (curr-cusp-a-dif (first curr-cusp-difs))
			  (curr-cusp-b-dif (second curr-cusp-difs)))
		     (cond ((< (sqrt (+ (expt curr-cusp-a-dif 2) (expt curr-cusp-b-dif 2)))
			       (sqrt (+ (expt best-cusp-a-dif 2) (expt best-cusp-b-dif 2))))
					;if the current cusp is closer to the cusp from the starting timestamp set it as
					;the new best cusp
			    (setf best-cusp curr-cusp))))))
    best-cusp
    ))

(defun all-cusps-2nd-der-trends (directory
				 &key (circularity nil)
				 (starting-timestamp 1))
					;circularity is a multiplier that returns the cusp speed metric as a fraction of
					;that of a perfect circle where those of a circle would be equal to 1
  (princ "function started")
  (terpri)
  (let* ((output (make-hash-table))
					;a hash table with the initial cusp as the key and list of its speed metric
					;(based on 2nd derivatives) values over time
	 (data nil)
	 (*backreaction-filename-digits* (find-backreaction-digits directory))
;;	 (cusps-1st-timestamp (read-and-find-cusps (backreaction-dump-filename directory starting-timestamp)))
	 )
    (loop for timestamp from starting-timestamp
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  for curr-cusps = (read-and-find-cusps (backreaction-dump-filename directory timestamp))
	  do (let* ((in-table nil)
		    (circularity-metric 1))
	       (multiple-value-bind
		   (a-hats a-dsigmas b-hats b-dsigmas)
		   (read-hats-dsigmas (backreaction-dump-filename directory timestamp))
		 (cond (circularity
			(setf circularity-metric (/ (* (* 2 (sqrt 2)) pi)
						    (average-total-sigma a-dsigmas b-dsigmas)))))
		 (maphash #'(lambda (k v)
			      (let* ((next-cusp (matching-cusp k curr-cusps (length a-hats) (length b-hats))))
				(cond (next-cusp
					;if there is a next cusp calculate its speed metric and add it to the list of
					;2nd der values
				       (destructuring-bind
					   (cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
						   cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle
						   a-2nd-der b-2nd-der intersection)
					   next-cusp
					 (declare (ignore cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
						   cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle intersection))
					 (let* ((cusp-speed-metric (* (/ 1 (sqrt (+ (expt a-2nd-der 2)
										    (expt b-2nd-der 2))))
								      circularity-metric)))
					   (setf (gethash k output) (push cusp-speed-metric v)))
				       (push next-cusp in-table)))
				      (t
					;if there is no next cusp simply put 0 into the speed metric list
				       (setf (gethash k output) (cons 0 v))))))
			  output)
		 (let* ((0-list nil))
		   (loop for i from starting-timestamp
			 to (1- timestamp)
			 do (push 0 0-list))
					;a 0 list that will be used as the beginning of the data for any cusps at the
					;current timestamp that didnt exist before	       
		   (loop for cusp in curr-cusps ;any cusps that didnt exist before will be added to the table here
			 do (cond ((not (find cusp in-table))
				   (destructuring-bind
				       (cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
					       cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle
					       a-2nd-der b-2nd-der intersection)
				       cusp
				     (declare (ignore cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
					       cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle intersection))
				     (let* ((cusp-speed-metric (* (/ 1 (sqrt (+ (expt a-2nd-der 2) (expt b-2nd-der 2))))
								  circularity-metric))
					    (curr-cusp-0-list 0-list)
					    (2nd-ders-list (push cusp-speed-metric curr-cusp-0-list)))
				       (setf (gethash cusp output) 2nd-ders-list)
				     )))))))))
    (princ "formating data 1")
    (terpri)
    (maphash #'(lambda (k v) ;reverses the speed metric lists since the code uses push
		 (let* ((reversed-list (reverse v)))
		   (setf (gethash k output) reversed-list)))		  
	     output) ;turns each list into a vector to more easily access the data
    (princ "formating data 2")
    (terpri)
    (maphash #'(lambda (k v)
		 (let* ((vector (coerce v 'vector)))
		   (push (list k vector) data)))
	     output)
					;(maphash #'(lambda (k v) (format t "~a => ~a~%" k v)) output)
    (princ "done formating and output")
    (terpri)
    data))

;wrapper function that takes the ouput from all-cusps-2nd-der-trends and stores a list of the strongest cusp at each
;timestamp
;strongest cusp is defined as the one with the largest speed metric
(defun strongest-cusp-per-timestamp (directory)
  (let* ((output nil)
	 (*backreaction-filename-digits* (find-backreaction-digits directory))
	 (cusp-2nd-der-trends (all-cusps-2nd-der-trends directory :circularity t)))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (let* ((strongest-cusp (nth 0 (nth 0 cusp-2nd-der-trends)))
		    (strongest-cusp-strength (aref (nth 1 (nth 0 cusp-2nd-der-trends))
						   (1- timestamp))))
	       (loop for cusp in cusp-2nd-der-trends
		     do (let* ((cusp-strength (aref (nth 1 cusp) (1- timestamp))))
			  (cond ((> cusp-strength strongest-cusp-strength)
				 (setf strongest-cusp (nth 0 cusp))
				 (setf strongest-cusp-strength (aref (nth 1 cusp) (1- timestamp)))))))
	       (push (list timestamp strongest-cusp strongest-cusp-strength) output)))
    (setf output (reverse output))))

;lists all possible cusps and removes any cusps that are in the same a' or b' and have a distance in sigma between them
;less than the threshold the function takes as an input
(defun find-possible-cusps-remove-conseq-in-sigma (directory threshold)
  (let* ((cusps nil)
	 (*backreaction-filename-digits* (find-backreaction-digits directory)))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (progn
	       (print timestamp)
	       (let* ((possible-cusps (remove-fake-cusps
				       (read-and-find-cusps (backreaction-dump-filename directory timestamp))))
		      (prev-cusp (first possible-cusps)))
		 (loop for cusp in possible-cusps
		       for i from 0
		       do (if (not (= i 0))
			      (destructuring-bind
				  (prev-cusp-a prev-cusp-a-fraction prev-cusp-b prev-cusp-b-fraction
					       prev-cusp-a1-to-cusp prev-cusp-cusp-to-a2
					       prev-cusp-b1-to-cusp prev-cusp-cusp-to-b2
					       prev-cusp-a-angle prev-cusp-b-angle
					       prev-cusp-a-second-derivative prev-cusp-b-second-derivative intersection)
				  prev-cusp
				(declare (ignore prev-cusp-a-fraction prev-cusp-b-fraction prev-cusp-cusp-to-a2
						 prev-cusp-cusp-to-b2 prev-cusp-a-angle prev-cusp-b-angle
						 prev-cusp-a-second-derivative prev-cusp-b-second-derivative intersection))
				(destructuring-bind
				    (cusp-a cusp-a-fraction cusp-b cusp-b-fraction cusp-a1-to-cusp cusp-cusp-to-a2
					    cusp-b1-to-cusp cusp-cusp-to-b2 cusp-a-angle cusp-b-angle
					    cusp-a-second-derivative cusp-b-second-derivative intersection)
				    cusp
				  (declare (ignore cusp-a-fraction cusp-b-fraction cusp-cusp-to-a2 cusp-cusp-to-b2
						   cusp-a-angle cusp-b-angle cusp-a-second-derivative cusp-b-second-derivative intersection))
				  (cond ((and
					  (and
					   (= prev-cusp-a cusp-a)
					   (= prev-cusp-b cusp-b))
					  (or
					   (< (abs (- prev-cusp-a1-to-cusp cusp-a1-to-cusp)) threshold)
					   (< (abs (- prev-cusp-b1-to-cusp cusp-b1-to-cusp)) threshold)))
					 (remove prev-cusp possible-cusps)
					 (remove cusp possible-cusps)))))))
		 (push (list (length possible-cusps) possible-cusps timestamp) cusps))))
    (reverse cusps)
    ))

;Runs num-repeat-cusps for all timestamps and dumps the data
(defun repeat-cusps-all-timestamps (directory)
  (let* ((output nil)
	 (*backreaction-filename-digits* (find-backreaction-digits directory)))
    (loop for timestamp
	  from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (progn
	       (print timestamp)
	       (let* ((curr-repeat-cusps (num-repeat-cusps timestamp directory)))
	       (destructuring-bind (num-repeats repeat-cusps initial-cusp) curr-repeat-cusps
		 (declare (ignore num-repeats initial-cusp))
		 (push (list (first repeat-cusps) timestamp) output)))))
    (reverse output)))

;outputs the largest kink angle in a and b for a set of a and b hats
(defun largest-kink-angle (a-hats b-hats)
  (mirror-image-let* ((n-a (length a-hats))
		      (largest-a-angle (loop for i-a below n-a
					      maximize (spherical-angle (aref a-hats i-a)
									(aref a-hats (mod (1+ i-a) n-a))))))
    (list largest-a-angle largest-b-angle)))
     

;outputs the average kink angle in a and b for a set of a and b hats
(defun avg-kink-angle (a-hats b-hats)
  (mirror-image-let ((output nil)
		     (n-a (length a-hats))
		     (total-a-angle 0)
		     (avg-a-angle 0))
    (mirror-images
     (loop for i-a below n-a
	   for curr-a-angle = (spherical-angle (aref a-hats i-a) (aref a-hats (mod (1+ i-a) n-a)))
	   do (setf total-a-angle (+ total-a-angle curr-a-angle)))
     (setf avg-a-angle (/ total-a-angle n-a)))
    (setf output (list avg-a-angle avg-b-angle))))

;gives a list of all the kinks in and b to visualize overall trends
(defun all-kink-angles (a-hats b-hats)
  (mirror-image-let ((output nil)
		     (n-a (length a-hats))
		     (a-angles nil))
    (mirror-images
     (loop for i-a below n-a
	   for curr-a-angle = (spherical-angle (aref a-hats i-a) (aref a-hats (mod (1+ i-a) n-a)))
	   do (push curr-a-angle a-angles)))
    (setf output (list (reverse a-angles) (reverse b-angles)))))

;takes one of the above functions as an input and applies it to the data in a certain directory for all timestamps
(defun kink-angle-trends (directory fn)
  (let*((output nil)
	(*backreaction-filename-digits* (find-backreaction-digits directory)))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
		 (read-hats-dsigmas (backreaction-dump-filename directory timestamp))
	       (declare (ignore a-dsigmas b-dsigmas))
	       (push (list (funcall fn a-hats b-hats) timestamp) output)))
    (reverse output)))

;finds the number of segments with 2nd derivatives below a threshold at a specific time
(defun num-small-segments (directory timestamp threshold)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory)))
    (mirror-image-let* ((total-small-a-seg-angle 0)
					;total sigma held in segments with individual lengths under the threshold
			(total-large-a-seg-angle 0)
					;total sigma held in segments with individual lengths above the threshold
			(percent-small-a-segments 0))
      (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	  (read-hats-dsigmas (backreaction-dump-filename directory timestamp))
	(mirror-images
	 (let ((n-a (length a-hats)))
	   (loop for i-a below n-a
		 for a-1 = (aref a-hats i-a)
		 for a-2 = (aref a-hats (mod (1+ i-a) n-a))
		 do (let* ((a-dsigma (aref a-dsigmas i-a))
			   (a-next-dsigma (aref a-dsigmas (mod (1+ i-a) n-a)))
			   (avg-a-dsigma (/ (+ a-dsigma a-next-dsigma) 2))
			   (total-a-angle (spherical-angle a-1 a-2))
			   (a-second-derivative (* (/ (average-total-sigma a-dsigmas b-dsigmas) (* 2.0 pi))
						   (/ total-a-angle avg-a-dsigma))))
		      (cond ((< a-second-derivative threshold)
			     (incf total-small-a-seg-angle total-a-angle))
			    (t
			     (incf total-large-a-seg-angle total-a-angle))))))
	 (setf percent-small-a-segments (/ total-small-a-seg-angle
					   (+ total-small-a-seg-angle total-large-a-seg-angle)))))
      (values percent-small-a-segments percent-small-b-segments total-small-a-seg-angle total-large-a-seg-angle
	      total-small-b-seg-angle total-large-b-seg-angle))))

;returns the percentage of slow segments below a threshold over time using num-small-segments because segment length 
;is proportional with speed
(defun percent-slow-segments (directory threshold)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (data nil)
	 (a-smoothing nil) ;called smoothing because the percentage of slow segments is a measure of smoothing
	 (b-smoothing nil))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (multiple-value-bind (a-small-percent b-small-percent) (num-small-segments directory timestamp threshold)
	       (push (list a-small-percent timestamp) a-smoothing)
	       (push (list b-small-percent timestamp) b-smoothing)))
    (push (coerce (reverse b-smoothing) 'vector) data)
    (push (coerce (reverse a-smoothing) 'vector) data)))

;returns the number and percentage of hats and sigma located in straight segments of string
(defun num-straight-pieces (directory timestamp)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory)))
    (mirror-image-let* ((num-a-straight-segments 0)
			(num-a-hats-in-straight-segments 0)
			(percent-a-hats-in-straight-segments 0)
			(sigma-in-straight-a-segments 0))
      (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	  (read-hats-dsigmas (backreaction-dump-filename directory timestamp))
	(mirror-images
	 (let ((a-straight-pieces (find-straight-pieces-hats-sigmas a-hats a-dsigmas)))
	   (loop for straight-piece in a-straight-pieces
		 do (destructuring-bind
			(starting-index num-hats sigma sum-kink-angles) straight-piece
		      (declare (ignore starting-index sum-kink-angles))
		      (incf num-a-straight-segments)
		      (incf num-a-hats-in-straight-segments num-hats)
		      (incf sigma-in-straight-a-segments sigma)))
	   (setf percent-a-hats-in-straight-segments (float (/ num-a-hats-in-straight-segments (length a-hats)))))))
      (values (coerce (list num-a-straight-segments percent-a-hats-in-straight-segments sigma-in-straight-a-segments)
		      'vector)
	      (coerce (list num-b-straight-segments percent-b-hats-in-straight-segments sigma-in-straight-b-segments)
		      'vector)))))

;returns the percentage of stright segments over time using num-straight-segments
(defun percent-straight-segments (directory)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (data nil)
	 (a-smoothing nil)
	 (b-smoothing nil))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (multiple-value-bind (a-straight-pieces b-straight-pieces) (num-straight-pieces directory timestamp)
	       (push (list a-straight-pieces timestamp) a-smoothing)
	       (push (list b-straight-pieces timestamp) b-smoothing)))
    (push (coerce (reverse b-smoothing) 'vector) data)
    (push (coerce (reverse a-smoothing) 'vector) data)))

;calculates the integrals of the 2nd derivatives of a and b squared as a measure of smoothing
;should converge over time
(defun int-2nd-derivative-squared (directory timestamp)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (total-int-a-2nd-der-sq 0)
	 (total-int-b-2nd-der-sq 0))
    (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	(read-hats-dsigmas (backreaction-dump-filename directory timestamp))
      (mirror-images
       (let ((n-a (length a-hats)))
	 (loop for i-a below n-a
	       for a-1 = (aref a-hats i-a)
	       for a-2 = (aref a-hats (mod (1+ i-a) n-a))
	       do (let* ((a-dsigma (aref a-dsigmas i-a))
			 (a-next-dsigma (aref a-dsigmas (mod (1+ i-a) n-a)))
			 (avg-a-dsigma (/ (+ a-dsigma a-next-dsigma) 2))
			 (total-sigma (average-total-sigma a-dsigmas b-dsigmas))
			 (total-a-angle (spherical-angle a-1 a-2))
			 (a-second-derivative (* (/ total-sigma (* 2.0 pi)) (/ total-a-angle avg-a-dsigma)))
			 (a-2nd-der-sq (* a-second-derivative a-second-derivative))
			 (int-a-2nd-der-sq (* a-2nd-der-sq avg-a-dsigma)))
		    (incf total-int-a-2nd-der-sq int-a-2nd-der-sq))))))
    (values total-int-a-2nd-der-sq total-int-b-2nd-der-sq)))

;collates the data from int-2nd-derivative-squared
(defun 2nd-der-int-data (directory)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (data nil)
	 (a-2nd-der-int nil)
	 (b-2nd-der-int nil))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (multiple-value-bind (a-total-int-2nd-der-sq b-total-int-2nd-der-sq)
		 (int-2nd-derivative-squared directory timestamp)
					 (push (list a-total-int-2nd-der-sq timestamp) a-2nd-der-int)
					 (push (list b-total-int-2nd-der-sq timestamp) b-2nd-der-int)))
    (push (coerce (reverse b-2nd-der-int) 'vector) data)
    (push (coerce (reverse a-2nd-der-int) 'vector) data)))

;special function that calculates only the data for the integral of the 2nd derivative of b squared in order to study the
;observed biased smoothing
(defun b-int-2nd-derivative-squared-trends (directory
					    &key (vector nil))
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (data nil))
    (loop for timestamp from 1
	  while (probe-file (backreaction-dump-filename directory timestamp))
	  do (let* ((int-b-2nd-der-sq-list nil))
	       (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
		   (read-hats-dsigmas (backreaction-dump-filename directory timestamp))
		 (declare (ignore a-hats))
		 (let ((n-b (length b-hats)))
		   (loop for i-b below n-b
			 for b-1 = (aref b-hats i-b)
			 for b-2 = (aref b-hats (mod (1+ i-b) n-b))
		         do (let* ((b-dsigma (aref b-dsigmas i-b))
				   (b-next-dsigma (aref b-dsigmas (mod (1+ i-b) n-b)))
				   (avg-b-dsigma (/ (+ b-dsigma b-next-dsigma) 2))
				   (total-sigma (average-total-sigma a-dsigmas b-dsigmas))
				   (total-b-angle (spherical-angle b-1 b-2))
				   (b-second-derivative (* (/ total-sigma (* 2.0 pi)) (/ total-b-angle avg-b-dsigma)))
				   (b-2nd-der-sq (* b-second-derivative b-second-derivative))
				   (int-b-2nd-der-sq (* b-2nd-der-sq avg-b-dsigma)))
			      (push int-b-2nd-der-sq int-b-2nd-der-sq-list)))
		   (cond ((eql vector t)
			  (setf int-b-2nd-der-sq-list (coerce (reverse int-b-2nd-der-sq-list) 'vector)))
			 (t
			  (setf int-b-2nd-der-sq-list (reverse int-b-2nd-der-sq-list))))
		   (push (list int-b-2nd-der-sq-list timestamp) data)))))
    (setf data (reverse data))))

;locates any 2nd derivative b-hats above a threshold to be graphed, plotted, and studied individually
(defun locate-large-2nd-der-b-hats (directory threshold)
  (let* ((output (make-hash-table))
	 (trends (b-int-2nd-derivative-squared-trends directory)))
    (loop for timestamp in trends
	  for time-index from 0
	  do (let* ((b-int-2nd-der-sq-list (first timestamp)))
	       (loop for int-2nd-der-sq in b-int-2nd-der-sq-list
		     for hat-index from 0
		     do (cond ((and (> int-2nd-der-sq threshold)
				    (not (gethash hat-index output)))
					;if a b 2nd der is above the threshold and is not already in the table add it to
					;the table
			       (setf (gethash hat-index output) (list time-index int-2nd-der-sq)))))))
    (loop for value being the hash-values of output
	  using (hash-key key)
	  do (format t "~&~A -> ~A" key value))))

;locates any cusps with strengths (calculated using 2nd derivatives in all-cusps-2nd-der-trends) to be studied
;individually
(defun locate-large-cusp-strengths (directory start end threshold)
  (let* ((output (make-hash-table))
	 (trends (all-cusps-2nd-der-trends directory :circularity t)))
    (loop for cusp in trends
	  do (destructuring-bind
		 (cusp-loc trends) cusp
	       (loop for timestamp from start
		     while (< timestamp end)
		     do (let* ((strength (aref trends timestamp)))
			 (cond ((and (> strength threshold)
				    (not (gethash cusp-loc output)))
			       (setf (gethash cusp-loc output) (list strength timestamp))))))))
    (loop for value being the hash-values of output
	  using (hash-key key)
	  do (format t "~&~A -> ~A" key value))))

;;Find coefficient in P_n propto n^{-4/3}.  See cusp-notes.tex
(defun cusp-gravitational-power-coefficient (cusp)
  (mirror-image-let* ((a-pp (cusp-info-a-pp cusp)) ;L a'', L b''
		      (alpha-a (3vector-length a-pp)))
    (let ((delta (spherical-angle (3vector-normalize a-pp) (3vector-normalize b-pp))))
      ;;This is an unsigned version of delta, but it doesn't matter because of invariance under reflection.
      (* (/ (* 16 (expt 1.5 1/3)
	       2)				;Because we integrate only to pi instead of 2pi
	    (expt pi 7/3)
	    (expt (* alpha-a alpha-b) 1/3))
	 (integrate
	  #'(lambda (phi-a)
	      (cusp-coefficient-integrand alpha-a alpha-b phi-a (+ phi-a delta)))
	  0.0 pi
	  :relative-error 1e-4)		;Otherwise get "cannot reach tolerance because of roundoff error"
	 ))))

;;Integral of power over a wedge divided by wedge angle.  alpha's are the lengths
;;of a'' and b'', phi-a is the angle from the observer direction to a'' and delta the angle from a'' to b''
(defun cusp-gravitational-power-coefficient-phi (alpha-a alpha-b phi-a delta)
      (* (/ (* 16 (expt 1.5 1/3))
	    (expt pi 7/3)
	    (expt (* alpha-a alpha-b) 1/3))
	 (cusp-coefficient-integrand alpha-a alpha-b phi-a (+ phi-a delta))))

;;Testing Burden 1,1 loop
(defun burden-cusp-gravitational-power-coefficient ()
  (cusp-gravitational-power-coefficient (make-cusp-info :a-pp (make-3vector (* 2 pi) 0.0 0.0)
							:b-pp (make-3vector 0.0 (* 2 pi) 0.0))))

(defun burden-cusp-gravitational-power-coefficient-phi (phi-a)
  (cusp-gravitational-power-coefficient-phi (* 2 pi) (* 2 pi) phi-a (/ pi 2)))

;;Total coefficient, direction of strongest cusp, and coefficient of that cusp alone
(defun total-cusp-gravitational-power-coefficient (cusps)
  (loop with max-power = 0
	with strongest-cusp-direction
	for cusp in cusps
	for power = (cusp-gravitational-power-coefficient cusp)
;;	do (format t "power ~S " power)
	sum power into total-power
	when (> power max-power)
	do (setq max-power power
		 strongest-cusp-direction (cusp-info-direction cusp))
	finally (return (values total-power strongest-cusp-direction max-power))))

;;Compute h3pm(a)/|sin phi-a sin phi-b| where a = alpha-a/alpha-b |sin phi-b / sin phi-a|^3
(defun cusp-coefficient-integrand (alpha-a alpha-b phi-a phi-b)
  (mirror-image-let ((sin-phi-a (sin phi-a)))
    (mirror-images
     ;;If sin phi-a or sin phi-b is tiny (mostly important when it is strictly zero), the result approaches a constant
     ;;Return that to avoid divergence in calculation.
     (when (< (abs sin-phi-a) 1e-10)
       (return-from cusp-coefficient-integrand
	 (/ (* 2/5 (expt Pi 3/2) (Gamma 11/6) ;coefficient of leading term a^-1/3 in h3
	       (expt (abs (/ alpha-a alpha-b)) -1/3)) ;Part from first factor in a.  Tiny sin phi-a cancels
	    (expt sin-phi-b 2)))))		      ;One from second factor in a, one from denominator
    (let ((sins (* sin-phi-a sin-phi-b))	      ;Normal calculation
	  (aparam (* (/ alpha-a alpha-b) (expt (abs (/ sin-phi-b sin-phi-a)) 3)))) ;See header comment
      (/ (h3pm aparam (plusp sins))
	 (abs sins)))))
