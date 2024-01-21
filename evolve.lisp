(in-package "CL-USER")

(defun main-calendar-next (&optional peek)
  (calendar-next *calendar* peek))

(defun main-calendar-add (time event)
  (calendar-add *calendar* time event))

(defun main-calendar-delete (time event)
  (calendar-delete *calendar* time event))

(defun main-calendar-empty-p ()
  (calendar-empty-p *calendar*))

(defun map-main-calendar (function)
  (map-calendar function *calendar*))

;;The current time is the time at which the last event was removed from
;;the main calendar
(defun current-time ()
  (or *reading-dumps*
      (calendar-current-time *calendar*)))

(defvar *advance-diamond-count*)
(defvar *last-advance-diamond-count* 0)

;;Set up to start a simulation
;;TOTAL-SIZE is the periodicity distance for the entire lattice.
;;SPLIT-FACTOR is the number of pieces that the lattice is split into along each of the 4 dimensions.
;;See SETUP-GEOMETRY
;;If READING-DUMPS is set, we initialize for that, rather than for simulation
;;If we have predecessors, we do only the first phase of initialization here, and finish after reading inputs
(defun initialize (&key total-size (split-factor 1) ijkl job-number reading-dumps (ijkl-origin zero-4vector)
			start	;Global starting time
			bh-size ; size of bh
			bh-number ;number of bhs
			bh-start) ;first time of the bh intersections
  (format t "~&Initializing...") (force-output)
  (when total-size
    (setq total-size (double-float total-size)))
  (setup-diamond-span)
  (setq *reading-dumps* reading-dumps)	;Say what we're doing
  (setq *vv-p* nil)			;Say not using VV unless that is initialized
  (setq *bh-size* bh-size)
  (setq *bh-number* bh-number)
  (setq *bh-start* bh-start)
  (cond (total-size				;Unless infinite volume
	 (setup-geometry :total-size total-size :split-factor split-factor :ijkl ijkl :job-number job-number
			 :ijkl-origin ijkl-origin))
	(t
	 (setq *job-number* (or job-number 0)
	       *global-location* zero-4vector
	       *total-size* nil)))	;Say infinite volume
  ;;Job starting time.  No earlier than beginning of cube.
  (let ((job-start (global-time (job-start-t))))
    (unless (and start (>= start job-start))
      (setq start job-start)))
  (format t "Starting at global time ~$.  " start)
  (setq *initial-time* (local-time start)) ;Convert to local time and install
  (setq *longest-edge* nil	       
	monster-length-scale (* monster-length-multiplier diamond-span)	;Reset in setup-longest-edge
	*calendar* nil			;Already read from predecessors first
	*cells* nil
	*final-diamonds* (make-hash-table :test #'eq)
	*read-diamonds* (make-hash-table :test #'eq)
	*advance-diamond-count* 0
	*last-advance-diamond-count* 0
	*compute-diamond-adjust-count* 0
	*intersections-performed* 0
	*suppressed-rejoinings* 0
	*rejoinings-performed* 0
	*intersections-unlucky* 0)
  (initialize-handles)
  (initialize-junctions)
  (cond ((or (predecessors-p)		;If we have predecessors, read from them before finishing init
	     *reading-dumps*)		;If reading dumps, just store everything
	 (setq *read-and-store* t))
	(t (initialize-more)			;If none, finish here
	   (setq *read-and-store* nil)))
  (gc :full t)				;Move static data to last generation, get rid of previous runs
  (format t "Initialization done.~%") (force-output)
  )

;;Second half of initialization, after reading from predecessors
(defun initialize-more ()
  (when (predecessors-p) (setup-longest-edge))
  (initialize-cells)			;Now that we know what we have to process, initialize cells system
  (setq *calendar* (make-main-calendar)) ;And calendar
  (setq *read-and-store* nil)		;Now go to normal processing
  (map-read-diamonds #'handle-new-diamond) ;Install diamonds in calendar and cells
  (clrhash *read-diamonds*)		;Get rid of previous copies
  )

(defun make-main-calendar ()
  (make-calendar (if *longest-edge* (/ *longest-edge* 2) diamond-span) :name "main"
		 ;;Start calendar at initial time.  This prevents intersections from being found
		 ;;an earlier times and detects bugs
		 :current-time (if *total-size* *initial-time* 0.0)))

;;Find the location (a 4vector) for any kind of event.
(defun event-location (event)
  (etypecase event
    (diamond (diamond-end event))
    (intersection (intersection-spacetime event))
    ))

(defun external-diamond-p (diamond)
  (handle-p diamond))

;;Process diamond found in calendar.  If diamond-finalp, all we do is remove from cells.  Otherwise
;;we advance the diamond
(defun evolve-diamond (diamond)
  (cond ((diamond-finalp diamond)
	 (maybe-delete-diamond-cells diamond)
	 (delete-from-calendar diamond))
	(t				;Normal case
	 (advance-diamond diamond))))



;;Move this diamond forward in time.
(defun advance-diamond (diamond)
  (with-modification-group
   (let* ((east (diamond-ne diamond))
	  (west (diamond-nw diamond)))
     (unless (and east west)
       ;; each future neighbors exist as a diamond or as :deleted
       (error "Can't advance diamond ~S unless future neighbors known E: ~S and W ~S South E: ~S and W ~S"
	      diamond east west (diamond-se diamond) (diamond-sw diamond)))
     (when (or (external-diamond-p east) (external-diamond-p west))
       (error "Can't advance diamond ~S which has external links"
	      diamond))
     (cond ((and (or (eq east :deleted) (eq east :BH) (eq east :BHeatit) (eq east :BHpropdel) (eq east :BHdeleted))
		 (or (eq west :deleted) (eq west :BH) (eq west :BHeatit) (eq west :BHpropdel) (eq west :BHdeleted)))
	    (report-progress "x"))
	   ((or (eq east :deleted) (eq west :deleted)) ;We're just propagating deletion
	    (report-progress "x")
	    (when (null (diamond-inertp diamond))
	      (format t "bh: ~S~%" (tag-bh (diamond-tag diamond)))
	      (format t "c: ~S~%" (diamond-countup diamond)))
	    (assert (diamond-inertp diamond)) ;only inert diamonds can border :deleted or :BH diamonds
	    (mirror-images		;Tell non-deleted side about deleted future
	     (unless (eq east :deleted)
	       (setf (diamond-nw east) :deleted))))
	   ((or (eq east :BHdeleted) (eq west :BHdeleted)) ;Propagating deletion of a BH string
	    (report-progress "x")
	    (mirror-images
	     (unless (eq east :BHdeleted)
	       (setf (diamond-nw east) :BHpropdel))))
	   ((or (eq east :BHpropdel) (eq west :BHpropdel)) ;Propagating deletion of a BH string
	    (unless (check-point-inbhhalo (diamond-start diamond)) ;check that the deleted bh diamond is inside a halo from the bh
	      (warn "Deleting diamond which is outside BH"))
	    (report-progress "x")
            (mirror-images
             (unless (eq east :BHpropdel)
               (setf (diamond-nw east) :BHpropdel))))
	   ((and (null (eq east :BH)) (null (eq west :BH)) (null (eq east :BHdeleted)) (null (eq west :BHdeleted))
		 (null (eq east :BHeatit)) (null (eq west :BHeatit))
		 (or (eq (diamond-e east) :BH) (eq (diamond-e east) :BHdeleted) (eq (diamond-e east) :BHpropdel) (eq (diamond-e east) :deleted))
		 (or (eq (diamond-w west) :BH) (eq (diamond-w west) :BHdeleted) (eq (diamond-w west) :BHpropdel) (eq (diamond-w west) :deleted))
		 (< (3vector-distance (diamond-end east) (diamond-end west)) 1E-8))
	    (format t "Deleting small loop from BH")
	    (setf (diamond-nw east) :BHpropdel)
	    (setf (diamond-ne west) :BHpropdel))
	   ((eq west :BH) (advance-diamond-BH-left  diamond))  ;added by SJW
           ((eq east :BH) (advance-diamond-BH-right diamond))  ;added by SJW
	   ((eq west :BHeatit) (advance-diamond-BH-left  diamond))
	   ((eq east :BHeatit) (advance-diamond-BH-right diamond)) 
	   ;;see BH-String.lisp for function advance-diamond-BH-left/right
	   (t
	    (advance-diamond-1 west east diamond))))
   (discard-object diamond)   ;Get rid of old diamond
   (incf *advance-diamond-count*)
   ))



;;Create a new diamond given the predecessors on the west and east and add it to the calendar.
;;PREDECESSOR is the diamond from which this one was made, or NIL if this diamond is being created in initialization.  
;;Intersections with that diamond are not permitted
(defun advance-diamond-1 (west east predecessor)
  (multiple-value-bind (tag count engulfp inertp deletep) (handle-loop-tags west east predecessor)
   (when deletep			;Delete the string?  No new diamond
     (when (tag-bh tag)
       (error "Deleting bh from middle"))
     (report-progress "x")
      (record-loop predecessor east) ;write that this loop has been deleted
      (incf *loop-counter*)
      (setf (diamond-nw east) :deleted ;"Install" deleted diamond
	    (diamond-ne west) :deleted))
    ;;Common code: Update tag on engulfment.  In deletion case, must have recorded loop first.
    (when engulfp		;Engulfment?  Update tag
      (report-progress "E")
      (setf (tag-last-position tag) (globalize-position (diamond-left east)))) ;Position of engulfment
    (unless deletep			;Normal advancement
      (report-progress ".")		;Normal advancement
      (let ((new (make-diamond :start (diamond-left east)
			       :left  (diamond-end west)
			       :right (diamond-end east)
			       :sw west 
			       :se east
			       :tag tag
			       :a-kink-created (diamond-a-kink-created west)
			       :b-kink-created (diamond-b-kink-created east)
			       :countup count
			       :inertp inertp)))
	(setf (diamond-nw east) new
	      (diamond-ne west) new
	      (diamond-se west) nil
	      (diamond-sw east) nil)
	(setf (diamond-end new) (compute-diamond-end new))
	(handle-new-diamond new :predecessor predecessor)
	new))))				;return new diamond

;;See if WEST-TAG was generated after EAST-TAG.  If by coincidence the creation times are the same, break the tie
;;in a consistent fashion
(defun tag-later-p (west-tag east-tag)
  (mirror-image-let ((west-time (4vector-t (tag-created-position west-tag))))
    (cond ((> west-time east-time) t)	;West later
	  ((> east-time west-time) nil)	;East later
	  (t
	   ;;Create handles if needed for future use and compare those to break tie
	   (compare-handles (get-tag-handle west-tag) (get-tag-handle east-tag))))))

;;Handle processing of loop tags for advance-diamond-1.  Returns new tag and count
;;and flags engulfp, inertp, and deletep.
(defun handle-loop-tags (west east predecessor)
  (mirror-image-let ((west-tag (diamond-tag west)))
    (cond ((eq west-tag east-tag)		;Same tag to both sides?
	   (cond ((eq west-tag (diamond-tag predecessor)) ;Predecessor same also?
		  (mirror-image-let ((west-count (diamond-countup west)))
		    ;;If one count larger, pass on
		    (cond ((> west-count east-count) (values west-tag west-count nil (diamond-inertp west)))
			  ((> east-count west-count) (values east-tag east-count nil (diamond-inertp east)))
			  ((= west-count (diamond-countup predecessor)) ;All 3 same
			   (values east-tag east-count nil (diamond-inertp east))) ;Pass on everything
			  (t		;East and west count are the same, but predecessor not
			   (assert (= (1- east-count) (diamond-countup predecessor))) ;Sanity check
			   (handle-loop-tags-engulf west east east-tag west-count)))))
		 (t			;East and west have same tag, but predecessor is different: first engulfment
		  (assert (and (zerop (diamond-countup west)) (zerop (diamond-countup east)))) ;Sanity check
		  (handle-first-engulfment east east-tag))))
	  ((tag-later-p west-tag east-tag) ;Different tags to east and west.  Return later time one.  Not inert.
	   (values west-tag 0))
	  (t (values east-tag 0)))))


;;First engulfment: compute loop length and store in tag
(defun handle-first-engulfment (east tag)
  (let ((ic-time (4vector-t (tag-created-position tag))) ;Global time of loop creation
	(now (global-time (4vector-t (diamond-left east))))) ;Global time of engulfment
    (setf (tag-xi tag) (/ (loop-length-i ic-time ic-time now) now))
    (values tag 1 t)))			;Return engulfment with count = 1

;;See if loop should be preserved because of *loop-preservation-threshold*
(defun preserve-loop-p (now tag)
  (and *loop-preservation-threshold*	;Feature turned on?
       (let* ((last-time (4vector-t (tag-last-position tag))) ;last global conformal engulfment time
	      (length (* 2 (- now last-time)))) ;Approximate present length is twice last period.
	 ;;Compare with specified fraction of current conformal time, but if that is too small, use 1.0
	 (>= length (* (max now 1.0) *loop-preservation-threshold*))))) ;Too big to delete

;;See if loop should be preserved for dumping
(defun preserve-loop-dump-p (now tag)
  (and *loop-preservation-dump-x*	;Feature turned on?
       (>= now *loop-preservation-dump-start*) ;Time to begin saving loops?
       (not (and *last-dump-time*		;But not after last dump?
		 (> now (+ *last-dump-time* fudge-global-coordinates))))
       (and (> (tag-xi tag) *loop-preservation-dump-x*) ;Keep if x_i larger than threshold
	    (not (tag-dumped tag)))))	;and not already dumped

;;Handle engulfment after the first.  Unless preserving, we increment the count and
;;sometimes make the diamond inert or deleted
(defun handle-loop-tags-engulf (west east tag count)
  (check-tag-count (incf count))	;New COUNT value
  (let ((now (global-time (4vector-t (diamond-left east))))) ;Global time of engulfment
    (cond ((or (preserve-loop-p now tag)
	       (preserve-loop-dump-p now tag) ;Wants it kept?
	       (< count (tag-minimum-inert-count tag))) ;Suppressed by previous unlucky intersections
	   (values tag count t))	;Engulfed, but not inert yet.
	  ((diamond-inertp west)	;Loop ready to be deleted.  Already inert?
	   (assert (diamond-inertp east)) ;Better be inert this way too
	   (values tag count t t t))	;Delete it.
	  (t				;Make inert now, delete next time
	   (values tag count t t)))))
	   
;;Minimum global time of the creation for a loop to be output, or NIL for all
(define-simulate-variable *loop-record-start* 25.0 double-float)

;; writes the deleted loop's x, p, and time of creation to a file.
;; Because the event-time of the last intercommutation is stored in predecessor
;; and the begin-inert time is stored in east, we can calculate both length and time with good accuracy. 
(defun record-loop (predecessor east)
  (when *loops-found* 			;Keeping these?
    (let* ((delete-event (globalize-position (diamond-end predecessor)))
	   (delete-time (4vector-t delete-event))
	   (tag (diamond-tag east))
	   (inert-event (tag-last-position tag)) ;Position of last engulfment, where it was made inert
	   (inert-time (4vector-t inert-event))	;Time it was made inert.
	   (ic-time (4vector-t (tag-created-position (diamond-tag predecessor)))) ;Time of event creating the loop
	   (4-momentum (standardize-position (4vector- delete-event inert-event))) ;comoving 4-momentum averaged over 1 period
	   (velocity (3vector-scale 4-momentum (/ 1 (4vector-t 4-momentum))))
	   (speed (3vector-length velocity))
	   (pf (/ speed 		;this is the momentum between inert and delete times
		  (sqrt (- 1.0 (expt speed 2.0)))))
	   (energy (loop-length-i ic-time inert-time delete-time))  ;energy is the comoving energy of the loop at time of formation.
	   (xi (/ energy ic-time))
	   (p-i (/ pf (scale-factor-ratio ic-time (/ (+ delete-time inert-time) 2.0)))))	;correct for redshifting since formation
      (when (or (null *loop-record-start*) ;Don't record loops at early times
		(> ic-time *loop-record-start*))
	(vector-push-extend xi *loops-found*) 
	(vector-push-extend p-i *loops-found*)
	(vector-push-extend ic-time *loops-found*)
	(when *log-loop-positions*
	  (let ((position (globalize-position (diamond-left east)))) ;The point at which the first :deleted diamond began
	    (vector-push-extend (4vector-x position) *loops-found*)
	    (vector-push-extend (4vector-y position) *loops-found*)
	    (vector-push-extend (4vector-z position) *loops-found*)))
	(when *log-loop-velocities*
	  (vector-push-extend (3vector-x velocity) *loops-found*)
	  (vector-push-extend (3vector-y velocity) *loops-found*)
	  (vector-push-extend (3vector-z velocity) *loops-found*))
	))))

(define-timer :loops)

;;Write loops from array into file
(defun write-recorded-loops ()
  (account-time :loops
    (format t "Writing ~D loops..." (/ (length *loops-found*) (+ 3
								 (if *log-loop-positions* 3 0))))
    (loop for x across *loops-found*
	  do (write-single x *loops-output*))
    (force-output *loops-output*)
    (format t "done.~%")))

;;Do various things that must be done when a new diamond has been created.
;;If the diamond has instead been reshaped, predecessor is :meiosis
;;Every diamond must come here after its end has been set or changed
(defun handle-new-diamond (new &key predecessor)
  (when *read-and-store*		;Reading dumps or pre-reading from predecessor
    (add-read-diamond new)		;Just store in table
    (return-from handle-new-diamond nil))
  (unless (point-mine (diamond-end new)) ;We don't own final point?  Then this diamond will only be sent to successors
    (setf (diamond-finalp new) t)
    (add-final-diamond new))		;Put in table to write to successors
  (maybe-add-to-calendar new)		;Put in calendars to evolve or to remove from cells if final
  (unless (diamond-inertp new)
    (install-diamond-bounding-box new)	;Set up bounding box
    (if (eq predecessor :meiosis)	;If shrinking, we moved intersections and don't need to check again
        (add-diamond-cells new)		;Just add to cells
        (check-for-intersections new predecessor) 
    ) ;Add and check intersection of new diamond with all other diamonds.
    ))

;;Diamond goes in calendar except if the end of the diamond would be after the end of our run, 
;;in which case it would never come up anyway.
(defun diamond-belongs-in-calendar (diamond)
  (and (not *read-and-store*)
       (or (null *total-size*)		;Infinite volume
	   (<= (4vector-t (diamond-end diamond)) *job-end*))))

;;Put diamond in calendar if it belongs there
(defun maybe-add-to-calendar (diamond)
  (when (diamond-belongs-in-calendar diamond) (add-to-calendar diamond)))

;;Delete diamond from calendars if it should be there
(defun maybe-delete-from-calendar (diamond)
  (when (diamond-belongs-in-calendar diamond) (delete-from-calendar diamond)))

(defun add-to-calendar (event)
  (main-calendar-add (4vector-t (event-location event)) event))

(defun delete-from-calendar (event)
  (main-calendar-delete (4vector-t (event-location event)) event))
 
;;Timestamps
(defstruct (timestamp (:include timed-event (function 'do-timestamp))
		      (:constructor make-timestamp (time interval)))
  (interval 0.1))

;;Print timestamp.  This is not controlled by *print-progress*.
(defun do-timestamp (structure)
  (let ((time (timestamp-time structure))
	(interval (timestamp-interval structure)))
    ;;Print necessary number of digits to distinguish one timestamp from another
    (format t " ~V$" (max 1 (ceiling (- (- (log interval 10.0)) 1e-12))) time)
    (force-output)
    (main-calendar-add (incf (timestamp-time structure) interval) ;Move to next time
		       structure)))

;;Memory usage reporting
(defstruct (usage-report (:include timed-event (function 'do-usage-report))
		      (:constructor make-usage-report (time interval)))
  (interval 1.0))

(define-simulate-variable *external-usage* nil) ;Whether to call the memory-usage program.

(defvar *start-real-time* 0)	      ;Set by SIMULATE at start of run
(defvar *start-run-time* 0)
(defvar *last-run-time* 0)
(defvar *last-real-time* 0)

;;Print memory usage
(defun do-usage-report (structure)
  (let* ((time (usage-report-time structure))
	 (interval (usage-report-interval structure))
	 (realtime (- (get-internal-real-time) *start-real-time*))
	 (runtime (- (get-internal-run-time) *start-run-time*))
	 (new-run-time (- runtime *last-run-time*))
	 (new-real-time (- realtime *last-real-time*))
	 (real-seconds (round realtime internal-time-units-per-second))
	 (run-seconds (round runtime internal-time-units-per-second))
	 (4-volume (- (job-volume-so-far time) (job-volume-so-far (- time interval))))
	 (speed (and (plusp new-run-time) (/ 4-volume (/ (float new-run-time) internal-time-units-per-second))))
	 (real-speed (and (plusp new-real-time) (/ 4-volume (/ (float new-real-time) internal-time-units-per-second))))
	 (advanced (- *advance-diamond-count* *last-advance-diamond-count*))
	 (diamond-speed (and (plusp new-run-time) (/ advanced (/ (float new-run-time) internal-time-units-per-second)))))
;;    (when (zerop (mod (round time 0.0001) 10000))
;;      (gc :full t)
;;      (format t "(Full GC)"))
    (format t "~&Time ~3$: Real ~D:~2,'0D:~2,'0D, CPU ~D:~2,'0D:~2,'0D~@[ (~D%)~]~
                ~@[, incremental CPU ~D%~]~@[, speed ~$~]~@[, real speed ~$~], calendar: ~D, advanced: ~D~@[, d/s: ~$~], ~D MB of dynamic space in use~:[.~%~;: ~]"
	    time
	    (truncate real-seconds 3600) (truncate (mod real-seconds 3600) 60) (mod real-seconds 60)
	    (truncate run-seconds 3600) (truncate (mod run-seconds 3600) 60) (mod run-seconds 60)
	    (and (plusp realtime) (round (* (/ runtime realtime) 100)))
	    (and (plusp new-real-time) (round (* (/ new-run-time new-real-time) 100)))
	    speed
	    real-speed
	    (calendar-event-count *calendar*)
	    advanced
	    diamond-speed
	    (round (sb-kernel:dynamic-usage) (expt 2 20))
	    *external-usage*)
    (force-output)
    (setq *last-run-time* runtime
	  *last-real-time* realtime
	  *last-advance-diamond-count* *advance-diamond-count*)
    (when *external-usage*
      (external-memory-usage))
    (main-calendar-add (incf (usage-report-time structure) interval) ;Move to next time
		       structure)))
  
;;Outside view of memory usage
(defun external-memory-usage ()
  (run-program "/cluster/home/k/o/kolum/strings/parallel/memory-usage" (list (format nil "~D" (sb-unix:unix-getpid))) :output t :error t))

(defun report-progress (format-control &rest format-arguments)
  (when *print-progress*
    (apply #'format t format-control format-arguments)
    (force-output)))

;;Do one thing and maybe print a character saying what we did.
(defun evolve-1 ()
  (when (main-calendar-empty-p)	;Nothing to do
    (error "There's nothing left to do"))
  (multiple-value-bind (time thing)	;Peek at next thing to do
      (main-calendar-next t)
    (etypecase thing
      (diamond
       (evolve-diamond thing))
      (intersection
       (report-progress "!")
       (perform-intersection thing))
      (timed-event ;Any kind of thing that happens at a pre-scheduled time
       (report-progress (timed-event-report-string thing))
       (main-calendar-delete time thing) ;All such things get removed first
       (funcall (timed-event-function thing) thing)))))

;;Stop at given time
(defstruct (stop-evolving (:include timed-event (function 'stop-evolving))
			  (:constructor make-stop-evolving (time))
			  ))

(defvar *stop-event* nil)

;;Perform all events that are before the given time.
(defun evolve-until (time)
  (evolve-until-1 time))

(defun evolve-until-1 (time)
  (if (> (current-time) time)
      "Time to evolve until already reached"
    (let ((*stop-event* (make-stop-evolving time)))
      (unwind-protect			;Don't leave stop event in calendar
	  (progn (main-calendar-add time *stop-event*) ;Stop at given time
		 (catch 'stop-evolving	;Throw here when time reached
		   (loop (evolve-1)))
		 (setq *stop-event* nil)) ;Normal exit -- don't need to delete
	(when *stop-event* (main-calendar-delete time *stop-event*)) ;On abort, remove event
	))))

;;Called when we have reached the time to stop.
(defun stop-evolving (structure)
  (declare (ignore structure))
  (throw 'stop-evolving nil)
  )

