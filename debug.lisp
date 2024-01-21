;;;Debugging
(in-package "CL-USER")

;;Tell whether object has given handle to decide whether to trace something
(defun object-has-handle (object code &optional job)
  (let ((handle (object-handle object nil))) ;handle if any
    (and handle (= (handle-code handle) code)
	 (or (not job)
	     (= job (handle-creator handle))))))

;;Find all copies of the object with a given handle on the creator.
(defun objects-with-handle (code job)
  (let ((results nil))
    (maphash
     #'(lambda (object handle)
	 (when (and (= code (handle-code handle))
		    (= job (handle-creator handle)))
	   (push object results)))
     *object-handle*)
    results))

;;Return the object with the given handle, if there's only one
(defun object-with-handle (code rank)
  (let ((objects (objects-with-handle code rank)))
    (when (cdr objects)
      (error "Multiple objects with handle code ~D and rank ~D: ~S" code rank objects))
    (car objects)))


(defun describe-calendar (calendar &optional start end)
  (map-calendar
   #'(lambda (time object)
       (when (and (or (null start) (<= start time))
		  (or (null end) (>= end time)))
	 (format t "~F: ~S~%" time object)))
   calendar))



(defvar *debug-trace-handles* nil)	;(CODE RANK CODE RANK ...) for handles to trace
(defvar *debug-trace-break* nil)	;Evaluate this to decide whether to break

;;Trace operations that involve the object with the specified handle.
;;Returns the object
(defun debug-trace-object (object format-string &rest arguments)
  (when (and *debug-trace-handles*
	     (loop for (code rank) on *debug-trace-handles* by #'cddr
		   thereis (object-has-handle object code rank)))
    (debug-trace-object-print object format-string arguments))
  object)

;;Trace operations that involve the handle itself
;;Returns the handle
(defun debug-trace-handle (handle format-string &rest arguments)
  (when (and *debug-trace-handles*
	     (loop for (code job) on *debug-trace-handles* by #'cddr
		   thereis (and (= code (handle-code handle))
				(= job (handle-creator handle)))))
     (debug-trace-object-print handle format-string arguments))
  handle)

(defun debug-trace-object-print (object format-string arguments)
  (declare (special *receive-source* *send-destination*))	;Not known yet
  (format t "~&Debug-trace: ~S: ~?" object format-string arguments)
  (cond (*debug-send-communicator*	;First check if in sending.  Send can be nested inside receive.
	 (format t " in ~S" (communicator-name *debug-send-communicator*))
	 (when *debug-send-argument-number*
	   (format t ", argument number ~D" *debug-send-argument-number*))
	 (format t " to ~S" *send-destination*))
	(*debug-receive-communicator*	;Now check for receive
	 (format t " in ~S" (communicator-name *debug-receive-communicator*))
	 (when *debug-receive-argument-number*
	   (format t ", argument number ~D" *debug-receive-argument-number*))
	 (format t " from ~S" *receive-source*)))
  (force-output)
  (when (eval *debug-trace-break*)
    (break "debug trace object")))
    


;;Return lists of a and b segment lengths
(defun segment-lengths ()
  (let ((as nil)
	(bs nil))
    (map-string-paths
     #'(lambda (diamond &optional start-junction end-diamond end-junction)
	 (declare (ignorable start-junction end-junction))
	 (format t "~:[L~;S~]" end-junction) (force-output)
	 (when end-diamond		;If not loop, do first diamond
	   (push (3vector-length (diamond-p diamond)) as)
	   (push (3vector-length (diamond-q diamond)) bs))
	 (let ((d diamond))
	   (loop
	    (when (eq d end-diamond)	;End diamond reached
	      (return nil))
	    (let ((next (diamond-ne d))) ;Try going northeast
	      (cond (next		;Something there
		     (push (3vector-length (diamond-q d)) bs)) ;New B
		    (t
		     (setq next (diamond-se d)) ;Nothing: try southeast
		     (unless next (error "no diamond to E"))
		     (push (3vector-length (diamond-p d)) as) ;New A
		     ))
	      (when (eq next diamond)	;Looped back to original diamond
		(return nil))
	      (setq d next)
	      )))
	 (format t ".") (force-output)))
    (values as bs)))

(defun histogram-segment-lengths (&rest gnuplot-keys &key (min 1e-10) (max 1.0) (bins 20) &allow-other-keys)
  (let* ((data (multiple-value-list (segment-lengths)))
	 (count 0)
	 (tables (loop repeat 2
		       collect (make-array bins :element-type 'fixnum :initial-element 0)))
	 (bin-size (/ (log (/ max min)) bins)))
    (loop for ab below 2
	  do (loop for length in (nth ab data)
		   with table = (nth ab tables)
		   when (> length min)
		   do (let ((slot (floor (log (/ length min)) bin-size)))
			(incf count)
			(unless (> slot bins)
			  (incf (aref table slot))))))
    (apply
     #'gnuplot 2 bins
     #'(lambda (plot point)
	 (if (eq point :title) (nth plot '("a" "b"))
	   (values (* min (exp (* (+ point 0.5) bin-size))) (aref (nth plot tables) point))))
    :logscale :x
    :title (format nil "Job number ~D, time ~$: ~D segments" *job-number* (current-time) count)
    gnuplot-keys)))
    
    
(defun plot-segment-lengths (&rest gnuplot-keys)
  (let ((data (multiple-value-list (segment-lengths))))
    (apply
     #'gnuplot 2 (loop for list in data maximize (length list))
     #'(lambda (plot point)
	 (if (eq point :title) (nth plot '("a" "b"))
	   (and (< point (length (nth plot data)))
		(values point (nth point (nth plot data))))))
    :logscale :y
    gnuplot-keys)))
	
(defvar *small-segment-count* 0)

;If count is above this, find-small-segments will signal an error
(defvar small-segment-threshold 1000)


(mirror-images
;;Go east to find a diamond with a new A value, or west to find one with new B.
;;If we encounter end-diamond before we find a new A, return NIL
;;If we advance to start-diamond, return NIL.
(defun diamond-debug-new-a (diamond start-diamond end-diamond)
  (loop until (eq diamond end-diamond)	;If we reach end without a new a, return NIL
	with found = nil
	for next = (diamond-ne diamond) ;Move NE: same A
	if next do (setq diamond next)	;keep going
	else do (setq diamond (diamond-se diamond) ;Go SE
		      found t)		;Found new A
	until (eq diamond start-diamond) ;If we looped, return NIL.  Don't return start-diamond, even if new a
	when found return diamond
	))

(defun find-small-segments-a (start-diamond end-diamond max-length min-count)
  (let ((diamond start-diamond))
    (block done
      (unless end-diamond		;If loop, see if all are small.  Otherwise, advance diamond to a large one
	(loop for count from 0
	      for length = (4vector-t (diamond-a diamond))
	      while (< length max-length) ;Skip qualifying diamonds
	      maximize length into max
	      unless (setq diamond (diamond-debug-new-a diamond start-diamond nil)) ;Advance to new A.  NIL if loop
	      do (when (>= count min-count)					    ;All segments are short.  Are there enough?
		   (format t "Loop with all ~D ~As shorter than ~6E~%" count :a max) ;Yes
		   (setq *small-segment-count* (max *small-segment-count* count))) ;Keep track for monster alerting
	         (return-from done nil))					   ;Return  in any case, because loop has been scanned
	(setq end-diamond diamond	;Now treat loop like a regular segment ending here and starting with next
	      diamond (diamond-e diamond)))
      (loop with count = 0
	    with max = 0.0
	    for length = (and diamond (4vector-t (diamond-a diamond)))
	    do (cond ((and diamond (< length max-length)) ;Qualifying segment?
		      (incf count)
		      (setq max (max max length)) ;Maximum length in this run
		      )
		     (t			;Large segment or end of string
		      (when (>= count min-count)
			(format t "~&~D consecutive ~As of max length ~8E at ~S~%" count :a max diamond)
			(setq *small-segment-count* (max *small-segment-count* count))) ;Keep track for monster alerting
		      (setq count 0 max 0.0)
		      (unless diamond (return nil))))
	    do (setq diamond (diamond-debug-new-a diamond nil end-diamond))))))
)

(defun find-small-segments (max-length min-count)
  (let ((*small-segment-count* 0))
    (map-string-paths
     #'(lambda (start-diamond &optional start-junction end-diamond end-junction)
	 (declare (ignore start-junction end-junction))
	 (find-small-segments-a start-diamond end-diamond max-length min-count)
	 ;;For B, go in reverse order so that mirror-images works
	 (if end-diamond
	     (find-small-segments-b end-diamond start-diamond max-length min-count)
	   (find-small-segments-b start-diamond nil max-length min-count))
	 ))
    (when (> *small-segment-count* small-segment-threshold)
      (error "~D small segments.  Check for monsters" *small-segment-count*))))


;;2^N complex numbers
(defun time-fft (n)
  (let* ((nn (expt 2 n))		;Count of complex numbers
	 (data (make-array (* 2 nn) :element-type 'double-float :initial-element 0.0)))
    (four1-time data nn)))

(defun test-fft (from to)
  (loop for n from from to to
	for time = (time-fft n)
	do (format t "~&~A: ~A: N=~D: ~Dms, ~Dns/op~%"
		   (machine-instance) (lisp-implementation-version)
		   n time (round (* 1000000.0 time) (* (expt 2 n) n)))))

(defun four1-time (data nn &key (isign 1))
  (declare (type (simple-array double-float (*)) data)
	   (type fixnum nn)
	   (type (integer -1 1) isign))
  (let (start
	(ops 0))
    (declare (fixnum ops))
    (locally (declare (optimize (safety 0) ;Avoid checking that various numbers still are fixnums
				speed))
      (prog ((wr 0d0) (wi 0d0) (wpr 0d0) (wpi 0d0) (wtemp 0d0) 
	     (theta 0d0) (tempr 0d0) (tempi 0d0) (j 0) (n 0) (m 0) 
	     (mmax 0) (istep 0))
	    (declare (type double-float wr wi wpr wpi wtemp theta tempr tempi)) 
	    (declare (type fixnum j n m mmax istep))

	    (setf n (* 2 nn)) 
	    (setf j 1) 
	    (do ((i 1 (+ i 2)))
		((> i n) t)
	      (declare (type fixnum i))
	      (when (> j i) 
		(setf tempr (aref data (1- j)))
		(setf tempi (aref data j)) 
		(setf (aref data (1- j)) (aref data (1- i)))
		(setf (aref data j) (aref data i)) 
		(setf (aref data (1- i)) tempr)
		(setf (aref data i) tempi))
	      (setf m (floor n 2))
	      label1
	      (when (and (>= m 2) (> j m))
		(setf j (- j m)) (setf m (floor m 2))
		(go label1))
	      (setf j (+ j m))) 

	    (setf mmax 2) 
	    (setq start (get-internal-run-time))
	    label2 
	    (when (> n mmax)
	      (setf istep (* 2 mmax))
	      (setf theta (/ 6.28318530717959d0 (* isign mmax)))
	      (setf wpr (* -2.0d0 (expt (sin (* 0.5d0 theta)) 2)))
	      (setf wpi (sin theta)) (setf wr 1.0d0) (setf wi 0.0d0)
	      (do ((m 1 (+ m 2)))
		  ((> m mmax) t)
		(declare (type fixnum m))
		(do ((i m (+ i istep)))
		    ((> i n) t)
		  (declare (type fixnum i))
		  (setf j (+ i mmax))
		  (setf tempr (+ (* wr (aref data (1- j)))
				 (* (* -1d0 wi) (aref data j))))
		  (setf tempi (+ (* wr (aref data j))
				 (* wi (aref data (1- j)))))
		  (setf (aref data (1- j)) (+ (aref data (1- i)) (- tempr)))
		  (setf (aref data j) (+ (aref data i) (* -1d0 tempi)))
		  (setf (aref data (1- i)) (+ (aref data (1- i)) tempr))
		  (setf (aref data i) (+ (aref data i) tempi))
		  ;;		  (incf ops 2)		;Two complex operations done here
		  )
		(setf wtemp wr)
		(setf wr (+ (+ (* wr wpr) (* (- wi) wpi)) wr))
		(setf wi (+ (+ (* wi wpr) (* wtemp wpi)) wi)))
	      (setf mmax istep)
	      (go label2))))
    (values (- (get-internal-run-time) start) ops)))

;;Find any diamond on a loop which does not have the normal relationship with its neighbors
(defun find-unusual-diamonds (start-diamond)
  (loop with diamond = start-diamond
	with *print-pretty* = nil
	for collect = nil
	do (setq diamond (diamond-e diamond))
	do (check-points-eq diamond)
	do (let ((nw (diamond-nw diamond))
		 (sw (diamond-sw diamond)))
	     (when (diamondp nw)
	       (unless (eq (diamond-end diamond) (diamond-right nw))
		 (format t "~S has end ~S but NW diamond ~D has right ~S~%"
			 diamond (diamond-end diamond) nw (diamond-right nw))
		 (setq collect t))
	       (unless (eq (diamond-left diamond) (diamond-start nw))
		 (format t "~S has left ~S but NW diamond ~D has start ~S~%"
			 diamond (diamond-left diamond) nw (diamond-start nw))
		 (setq collect t)))
	     (when (diamondp sw)
	       (unless (eq (diamond-left diamond) (diamond-end sw))
		 (format t "~S has left ~S but SW diamond ~D has end ~S~%"
			 diamond (diamond-left diamond) sw (diamond-end sw))
		 (setq collect t))
	       (unless (eq (diamond-start diamond) (diamond-right sw))
		 (format t "~S has start ~S but SW diamond ~D has right ~S~%"
			 diamond (diamond-start diamond) sw (diamond-right sw))
		 (setq collect t))))
	when collect collect diamond
	until (eq diamond start-diamond)))

;;Check that edges of diamonds created by intersection are not horribly non-null
(mirror-images
 (defun check-diamond-edge-a (diamond b)
   (declare (optimize debug))
   (let* ((edge (4vector- (diamond-position diamond :a 1.0 :b b)
			 (diamond-position diamond :a 0.0 :b b))) ;Line crossing through intersection point
	  (v (/ (3vector-length edge) (3vector-t edge))))		 ;Velocity of edge: should be 1
     (cond ((< v 0.6)
	    (error "~S being split with new ~A edge velocity ~F~%" diamond :a v))
	   ((< v 0.9)
	    (warn "~S being split with new ~A edge velocity ~F~%" diamond :a v))))))

(defun check-diamond-edges (diamond a b)
  (check-diamond-edge-a diamond b)
  (check-diamond-edge-b diamond a))

(defun check-intersection-edges (intersection)
  (let* ((d1 (intersection-diamond-1 intersection))
	 (d2 (intersection-diamond-2 intersection))
	 (vector-a-b (intersection-vector-parameters intersection))
	 (a1 (vector-parameters-a1 vector-a-b))
	 (b1 (vector-parameters-b1 vector-a-b))
	 (a2 (vector-parameters-a2 vector-a-b))
	 (b2 (vector-parameters-b2 vector-a-b)))
    (check-diamond-edges d1 a1 b1)
    (check-diamond-edges d2 a2 b2)))

(defun plot-nearby-segment-lengths (diamond)
  (loop initially (loop repeat 100 do (setq diamond (diamond-w diamond)))
	repeat 200 for index from -100
	collect (list index (4vector-t (diamond-a diamond)) (4vector-t (diamond-b diamond))) into data
	do (setq diamond (diamond-e diamond))
	finally (gnuplot 2 (length data)
			 #'(lambda (plot point)
			     (if (eq point :title) (nth plot '("A" "B"))
			       (let ((entry (nth point data)))
				 (values (first entry) (nth (1+ plot) entry))))))))

;;Load previously dumped loop and simulate again with same conditions
;;It might be better to use resimulate
(defun redo-one-loop (file end)
  (read-one-loop file)
  (let ((loop (transform-loop-coordinates (longest-loop) #'(lambda (x) (4vector+ x (3to4vector zero-3vector 50))))))
    (simulate #'(lambda() 
		  (initialize-more)
		  (reset-loop-tags loop)
		  (loop for diamond = (diamond-e loop) then (diamond-e diamond)
			do (handle-new-diamond diamond)
			until (eq diamond loop)))
	      :end end :overwrite t :era :radiation :start 6.0 :print-progress t :size 200 :time-offset -50.0)))

;;Reset tags and kink-creation data as though this loop had just been created
(defun reset-loop-tags (start-diamond)
  (loop for diamond = start-diamond then (diamond-e diamond)
	for start = (globalize-position (diamond-start diamond))
	do (setf (diamond-a-kink-created diamond) start
		 (diamond-b-kink-created diamond) start
		 (diamond-tag diamond) (create-loop-tag start))
  	until (eq (diamond-e diamond) start-diamond)))
