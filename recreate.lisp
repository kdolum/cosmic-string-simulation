;;;Create data structures from information about a loop
(in-package "CL-USER")

;; worldsheet constructor functions will warn (w/ appropriate argument set) if they build
;; a diamond whose edge length is less than this value.
(defvar *create-worldsheet-minimum-edge* minimum-diamond-width)
(declaim (double-float *create-worldsheet-minimum-edge*))

;;Construct a piece of the worldsheet from a-hats, b-hats, and their corresponding delta-sigmas.
;;We make a loop on the worldsheet with no diamond being the successor of any other.  It will have have Na+Nb diamonds.
;;We make the initial diamond using indices START-A and START-B.  Then we make half the others going to the east from
;;the initial diamond than half going to the west.  We call the test-function with the latest diamond and the
;;direction we're going, :EAST or :WEST.  If it returns true, the next diamond will be north (and east or west)
;;of the previous one.  If NIL, it will be south.
;;We return the first time we create.
;;If TAGS is set, we put in loop and kink-created tags.
;;If INSTALL is set, we call handle-new-diamond on each diamond.
;;If OPEN-LOOP is set, we do not try to put the two parts of the loop together.
;;If SMALL-L is set, we print a warning when we construct a diamond with an edge length below *create-worldsheet-minimum-edge*.
;;INITIAL-POSITION becomes the left corner of the first diamond.  If not set, start at the origin with INITIAL-TIME or 0.
;;We work toward the east.  Since a-hats is in order toward the west, we have to go backwards through it.
(defun create-worldsheet (a-hats a-dsigmas b-hats b-dsigmas test-function
				 &key (start-a 0) (start-b 0) (initial-time nil) (initial-position nil)
				 tags install open-loop small-l)
  (declare (optimize speed)
	   (fixnum start-a start-b)
	   (simple-vector a-hats b-hats)
	   (function test-function))
  (when open-loop (warn "CREATE-WORLDSHEET :OPEN-LOOP T has not been tested in the rewritten version"))
  (mirror-image-let ((n-a (length a-hats)))
    (let* ((initial-position (or initial-position
				 (make-4vector 0.0 0.0 0.0 (or initial-time 0.0))))
	   (first-diamond (create-first-worldsheet-diamond
			   a-hats a-dsigmas b-hats b-dsigmas start-a start-b initial-position :small-l small-l))
	   (total-count (+ n-a n-b)) ;Number of diamonds to create including initial
	   (east-count (floor (1- total-count) 2))	     ;Number of new diamonds to create going eastward
	   (west-count (- total-count 1 east-count)))	     ;Number to make going westward
      (declare (fixnum total-count east-count west-count))
      ;;Variables for current diamond in both directions.  They will be left pointing to the last diamonds.
      (mirror-image-let ((east-diamond first-diamond)
			 (east-a-index start-a)
			 (east-b-index start-b)
			 (east-a-done nil) ;Flag set when we have wrapped around and returned to the original A index
			 (east-b-done nil))
	;;Need smaller data types than fixnum here, so that they can be added and subtracted.
	(declare (type (unsigned-byte 32) east-a-index east-b-index west-a-index west-b-index))
	(mirror-images				     ;Go both ways.  The code written below is for eastward.
	 (dotimes (i east-count)		     ;Make desired number of diamonds
	   (cond ((funcall test-function east-diamond :east) ;See which way to go next
		  ;;Hats lists are in order of increasing time, so whether we go east or west, we always
		  ;;want the next index
		  (when east-b-done	;We're back to the original segment and we're trying to advance again
		    (error "Ran out of ~A segments going ~A" :b :east))
		  (setf east-b-index (next-index-wrapping east-b-index n-b)) ;NE
		  (when (= east-b-index start-b)					 ;Back to original b'?
		    (setq east-b-done t))						 ;That's OK, but we can't advance again
		  (setf east-diamond (build-north-east-diamond east-diamond b-hats b-dsigmas east-b-index :small-l small-l)))
		 (t			;SE step
		  ;;Since we're going backward in time, we want the previous index
		  (when east-a-done
		    (error "Ran out of ~A segments going ~A" :a :east))
		  (setf east-a-index (previous-index-wrapping east-a-index n-a))
		  (when (= east-a-index start-a) ;Back to original a'?
		    (setq east-a-done t))
		  (setf east-diamond (build-south-east-diamond east-diamond a-hats a-dsigmas east-a-index :small-l small-l))))))
	(unless open-loop
	  ;;If connecting the diamonds requires a new b going from easternmost to the westernmost, the westernmost
	  ;;should be one larger.  If it requires a new a, the easternmost should be one larger
	  (mirror-image-let ((b-difference (- (mod (+ (- west-b-index east-b-index) (floor n-b 2)) n-b)
					      (floor n-b 2)))) ;map into -Nb/2 ... Nb/2
	    (block connected
	      (mirror-images
	       ;;To get a new B, we go NE from easternmost to westernmost.  For a new A, we go NW from westernmost
	       (when (and (zerop a-difference) (= b-difference 1))
		 (unless (funcall test-function east-diamond :east)
		   (error "To close worldsheet, we must go ~A from ~A diamond, but test function says ~A"
			  :ne :east :se))
		 (return-from connected (connect-created-worldsheet-ne east-diamond west-diamond total-count))))
	      (error "Halves of created worldsheet don't connect properly. A-difference = ~D, B-difference = ~D"
		     a-difference b-difference))))
	(when tags
	  (loop for this = west-diamond then (diamond-e this)
		do (add-diamond-tags this)
		until (eq this east-diamond)))
	(when install
	  (loop for this = west-diamond then (diamond-e this)
		do (handle-new-diamond this)
		until (eq this east-diamond))) ;Stop after doing last diamond
	first-diamond))))

;;This creates the first diamond of the wordlsheet from the a-hats and b-hats at the given time
(defun create-first-worldsheet-diamond (a-hats a-dsigmas b-hats b-dsigmas i-a i-b initial-position &key small-l)
  (let* ((left initial-position)
	 (a-prime (4vector-scale (3to4vector (aref a-hats i-a) 1.0) (/ (aref a-dsigmas i-a) 2)))
	 (b-prime (4vector-scale (3to4vector (aref b-hats i-b) 1.0) (/ (aref b-dsigmas i-b) 2)))
	 (start (standardize-position (4vector- left a-prime)))
	 (right (standardize-position (4vector+ start b-prime)))
	 (end (standardize-position (4vector+ left b-prime))))
    (mirror-images
     (when (and small-l (< (3vector-length a-prime) *create-worldsheet-minimum-edge*))      ;We check that there are not very small diamonds
       (warn "First diamond very small ~A: ~S" :a a-prime)))
    (make-diamond :start start
		  :left left
		  :right right
		  :end end
		  )))

(mirror-images
;;Builds a south diamond and links it to the rest of the worldsheet already constructed
(defun build-south-east-diamond (diamond a-hats a-sigmas a-index &key small-l)
  (declare (optimize speed)
	   (type simple-vector a-hats)
	   (type (simple-array double-float (*)) a-sigmas)
	   (type diamond diamond)
	   (type fixnum a-index))
  (let* ((left (diamond-start diamond))
	 (end (diamond-right diamond))
	 (aprime (4vector-scale (3to4vector (aref a-hats a-index) 1.0) (/ (aref a-sigmas a-index) 2)))
	 (start (standardize-position (4vector- left aprime)))
	 (right (standardize-position (4vector- end aprime))))
    (declare (type 4vector left start right end aprime))
    (when (and small-l (< (3vector-length aprime) *create-worldsheet-minimum-edge*))    ;We check that there are not very small diamonds
      (locally (declare (optimize (speed 0)))  ;Avoid compiler notes about WARN arguments
	(warn "South diamond very small A: ~S, sigma ~D index ~D" aprime (aref a-sigmas a-index) a-index)))
    (let ((new-diamond (make-diamond 
			:start start
			:left left
			:right right
			:end end
			:nw diamond
			)))
      (setf (diamond-se diamond) new-diamond)
      new-diamond)))

;;Builds a north diamond and links it to the rest of the worldsheet already constructed
(defun build-north-east-diamond (diamond b-hats b-sigmas b-index &key small-l)
  (declare (optimize speed)
	   (type simple-vector b-hats)
	   (type (simple-array double-float (*)) b-sigmas)
	   (type diamond diamond)
	   (type fixnum b-index))
  (let* ((left (diamond-end diamond))
	 (start (diamond-right diamond))
	 (bprime (4vector-scale (3to4vector (aref b-hats b-index) 1.0) (/ (aref b-sigmas b-index) 2)))
	 (right (standardize-position (4vector+ start bprime)))
	 (end (standardize-position (4vector+ left bprime))))
    (declare (type 4vector left start right end bprime))
    (when (and small-l (< (3vector-length bprime) *create-worldsheet-minimum-edge*))
      (locally (declare (optimize (speed 0))) ;Avoid compiler warnings about WARN arguments
	(warn "North diamond very small B: ~S, sigma ~D index ~D" bprime (aref b-sigmas b-index) b-index)))
    (let ((new-diamond (make-diamond 
		       :start start
		       :left left
		       :right right
		       :end end
		       :sw diamond
		       )))
      (setf (diamond-ne diamond) new-diamond)
      new-diamond)))

;;Connect the ends of the string made by create-worldsheet.  The second diamond goes NE of the first
(defun connect-created-worldsheet-ne (diamond1 diamond2 total-n)
  (create-worldsheet-check-closure (diamond-end diamond1) (diamond-left diamond2)
				   (diamond-right diamond1) (diamond-start diamond2) total-n)
  (setf (diamond-end diamond1) (diamond-left diamond2) ;Use points from diamond2 in both
	(diamond-right diamond1) (diamond-start diamond2))
  (when (diamond-nw diamond1)	;If penultimate diamond to NW, it was sharing end of diamond1, so change that too
    (setf (diamond-right (diamond-nw diamond1)) (diamond-end diamond1)))
  (setf (diamond-ne diamond1) diamond2
	(diamond-sw diamond2) diamond1))

)					;mirror-images


;;Accuracy required for the loop to close (i.e., the sum of a' and b' to vanish)
;;It is multiplied by the square root of the number of elements
(defvar create-worldsheet-closure-tolerance close-hats-default-tolerance)

;;Instead of points in the last diamond, we're going to use saved points from the first diamond.  Check they are
;;not too far.
(defun create-worldsheet-check-closure (last1 first1 last2 first2 total-n)
  (let ((distance (max (4vector-Euclidian-distance last1 first1) (4vector-Euclidian-distance last2 first2)))
	(tolerance (* create-worldsheet-closure-tolerance (sqrt (double-float total-n)))))
    (when (> distance tolerance)
      (warn "Adjusting coordinates of supposedly matching points by ~S.  Tolerance ~S" distance tolerance))))

;;Put in kink-created and loop tags
(defun add-diamond-tags (diamond)
  (let ((global-position (globalize-position (diamond-start diamond))))
    (mirror-images (setf (diamond-a-kink-created diamond) global-position))
    (setf (diamond-tag diamond) (create-loop-tag global-position))))

;;Create worldsheet with all diamonds overlapping a given time.  To avoid floating point comparison
;;problems we create the first diamond with its left corner slightly after the given time.
;;You can also give 4vector, in which case we keep the spatial position
(defun create-worldsheet-at-time (a-hats a-dsigmas b-hats b-dsigmas time-or-position &rest keys)
  ;;This fudge factor now matters only when there are other diamonds in the loop that have the same
  ;;starting time as this one
  (let* ((time (if (numberp time-or-position) time-or-position (4vector-t time-or-position)))
	 (initial-time (- time (* 10 fudge-global-coordinates)))) ;Make diamond definitely overlap given time
    (apply #'create-worldsheet a-hats a-dsigmas b-hats b-dsigmas
	   #'(lambda (diamond direction)	;If new corner is earlier than time, go north, otherwise south
	       (< (4vector-t (choose-mirror-image direction (diamond-right diamond))) time))
	   :initial-position (3to4vector (if (numberp time-or-position) zero-3vector time-or-position) initial-time)
	   keys)))

;;;Read in dump, smooth, and simulate again (in one processor).  Ken Olum, 2/2016
;;This does not work if there any strings linking the boundary conditions, because the single-processor
;;simulation does not have periodic boundary conditions.

(defun get-dump-strings (directory time)
  (setq directory (merge-pathnames directory batch-root-directory))
  (read-dumps directory :time time)
  (let ((result nil))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore))) ;Ignore open strings
     #'(lambda (d)					 ;"Loops", including everything but half-deleted small loops
	 (push d result)))
    result))

(defun resimulate (directory new-directory time initial-time end-time &rest simulate-keys)
  ;;Get one diamond from each loop in the dump.  This preserves them against reuse
  (setq new-directory (merge-pathnames new-directory batch-root-directory))
;;  (check-old-directory new-directory nil nil nil)
  (let ((strings (get-dump-strings directory time)))
    #| (dolist (string strings)
      (multiple-value-bind (a-hats a-sigmas) (get-a-data-sigmas string)
	(multiple-value-bind (b-hats b-sigmas) (get-b-data-sigmas string)
	  (mirror-image-let ((b-length (aref b-sigmas (1- (length b-sigmas)))))
	    (unless (fudge= a-length b-length fudge-global-coordinates)
	      (warn "A length ~S B length ~S.  Giving up." a-length b-length)
	      (return-from resimulate string)))))) |#
    (apply #'simulate #'(lambda () (recreate-old-strings strings))
	   :output-directory new-directory
	   :size *total-size*		;Must have same periodicity as before
	   :start initial-time :end end-time 
	   simulate-keys)))

;;Re-create strings from dump
(defun recreate-old-strings (strings &optional (time 0.0))
  (dolist (string strings)
    (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) (get-ab-data-dsigmas string)
      (create-worldsheet-at-time a-hats a-dsigmas b-hats b-dsigmas time :install t :tags t)
      )))

