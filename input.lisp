;;;Input from predecessors, reading of snapshots and spectra
(in-package "CL-USER")

(defvar *read-string-path*)		;Positions of string being read
(defvar *read-string-rank*)		;Rank of processor whose output is currently being read

(defun predecessors-p (&optional (job-number *job-number*) (split-factor *split-factor*))
  (>= job-number (expt split-factor 3)))

;;Tell whether this job intersects initial condition surface
(defun initialization-surface-p (&optional (job-number *job-number*) (split-factor *split-factor*))
  (< job-number (* 4 (expt split-factor 3))) ;4 layers of it
  )

;;Read predecessor or dump files.
;;The caller should call CONNECT-JUNCTIONS
(defun read-files (files &key (verbose t) report-n)
  (let (*receive-start-junction* *receive-previous-diamond* *receive-first-diamond*)
    (start-receive-string)
    (loop with count = 0
	  for file in files
	  for *receive-count-diamonds* = 0
	  do (when verbose (format t "~&Reading ~A: " file) (force-output))
	  do (if *read-dump-diamonds*
		 (read-file file :verbose verbose)
	       (read-file-segments file :verbose verbose))
	  do (when verbose
	       (format t "~D diamonds.~%" *receive-count-diamonds*)
	       (force-output))
	  do (incf count)
	  do (when (and report-n (zerop (mod count report-n)))
	       (format t "~D " count)
	       (force-output)))
    (when report-n (terpri))))

(defun read-initial-strings ()
  (when (predecessors-p)		;predecessors?
    (maybe-time :input
       (if *predecessor-hosts*
	   ;;If multiple sets of predecessors, offset successive ones to prevent handle conflicts
	   (loop for host in *predecessor-hosts*
		 for predecessor from 0
		 for *receive-handle-job-offset* = (* multiple-predecessor-job-offset-multiplier (floor predecessor 4))
		 when host
		 do (read-files (list (predecessor-file predecessor))))
	 (read-files
	  (loop for direction below 4
		collect (predecessor-file direction)))))
    (connect-junctions))
  (when *read-and-store*		;Not initialized yet
    (initialize-more))
  )

#|

;; Should output be a path? i.e. diamond a1 b1 a2 b2 
;; Reads in dump-files, averaging over p's and q's, giving a coarse grained version of the strings;
;; Output is an array of p's and q's, indexed by accumulated length (energy) sigma.
;; NOT YET WRITTEN
(defun read-coarse-grain-dumps  (*output-directory* &key time step wait)
  (multiple-value-bind (files time size split-factor time-offset) (dump-files *output-directory* :time time :step step)
    (when wait				;Wait for all files to be present if requested
      (dolist (file files)
	(loop until (probe-file file)
	      do (sleep 0.01))))
    ;;Initialize enough data structures to be able to read dumps
    (initialize :total-size size :split-factor split-factor :ijkl-origin (3to4vector zero-3vector time-offset)
		:reading-dumps time)
    (let ((*dump-time* time)
	  (*dump-length* 0.0)) 		;Accumulate total length here
      (read-files files :verbose nil)
      (connect-junctions)
      (values time *dump-length*))))
|#

;;Return a list of the times which dumps were made in the given directory
(defun dump-times (directory)
  (run-info-dump-times (read-run-info-file directory)))

(defun length-times (directory)
  (run-info-length-times (read-run-info-file directory)))

;;Returns a list of old dump files for the given time or step.
;;Returns run-info as second value, time used as third.
(defun dump-files (*output-directory* &key time step)
  (when (and time step)
    (error "Shouldn't specify both time and step"))
  (unless (or time step)
    (error "Must specify either time or step"))
  (let* ((info (read-run-info-file))
	 (total-size (run-info-total-size info))
	 (dump-times (run-info-dump-times info))
	 (time-offset (run-info-time-offset info)))
      (setup-split-factor (run-info-split-factor info))	;Install globally
      (unless dump-times
	(if (run-info-length-times info)
	    (error "Only lengths were dumped"))
	(error "According to ~A, no dumps were made" (run-info-file)))
      (cond (time
	     (setq time (coerce time 'double-float))
	     (setq step (position time dump-times :test #'fudge=global))
	     (unless step
	       (error "Could not find a dump at time ~F.  Dumps were made at ~S" time dump-times)))
	    (t
	     (setq time (nth step dump-times))))
      (values
       (if total-size
	   (loop for job from 0
		 for start = (global-job-start total-size time-offset job)
		 for end = (global-job-end total-size time-offset job)
		 until (> start time)	;Stop when we reach a job that started after the dump time
		 when (< time end)	;Did job run through the time?
		 collect (dump-file job
				    (if (run-info-dump-interval info) ;Old format.  Step numbers are per-job
					;;Decreased by any dumps that job didn't do at all
					(- step (loop for time in dump-times
						      count (<= time (- start fudge-global-coordinates))))
				      step)))

	 (list (dump-file 0 step)))		;Infinite volume: only one job
       info time)))
  
;;Read and reconnect all dump files for the given time or dump step
;;If WAIT, wait for files to appear if they are not present.
;;Returns time, even if step supplied.
;;Returns total length as second value
(defun read-dumps (*output-directory* &key time step wait (verbose nil))
  (multiple-value-bind (files info time) (dump-files *output-directory* :time time :step step)
    (when wait				;Wait for all files to be present if requested
      (dolist (file files)
	(loop until (probe-file file)
	      do (sleep 0.01))))
    ;;Initialize enough data structures to be able to read dumps
    (initialize :total-size (run-info-total-size info) :split-factor (run-info-split-factor info)
		:ijkl-origin zero-4vector
		:reading-dumps time)
    (setf *dumped-segments* nil)
    (clrhash *start-dumped-segments*)
    (clrhash *end-dumped-segments*)
    (let ((*dump-time* time)
	  (*dump-length* 0.0) 		;Accumulate total length here
	  (*allow-missing-predecessor-files* t) ;Dumps usually have skipped jobs which don't have files
	  (*count-missing-predecessor-files* 0)
	  (file-count (length files)))
      (format t "~&~D files: " file-count) (force-output)
      (read-files files :verbose verbose :report-n 1000)
      (when (plusp *count-missing-predecessor-files*)
	(if (= file-count *count-missing-predecessor-files*)
	    (error "All ~D files were missing" file-count)
	  (format t "~D files were missing out of ~D~%" *count-missing-predecessor-files* file-count)))
      (if *read-dump-diamonds*
	  (connect-junctions)
	(connect-junctions-segments))
      (when (null *read-dump-diamonds*)
	*dumped-segments*)
      (values time *dump-length*))))

;;Read a single loop written by dump-one-loop
(defun read-one-loop (file)
  ;;Initialize for infinite volume
  (initialize :total-size nil :reading-dumps :unknown)
  (let ((*dump-time* :unknown)
	(*dump-length* 0.0))
    (read-files (list file) :verbose t)
    *dump-length*))

;;Return number of bytes of disk space used by a dump
(defun dump-disk-usage (directory time &key (block-size 4096))
  (let ((blocks 0)
	(existing 0)
	(small 0)
	(all 0))
    (dolist (file (dump-files directory :time time))
      (incf all)
      (with-open-file (stream file :element-type '(unsigned-byte 8) :if-does-not-exist nil)
	(when stream
	  (incf existing)
	  (let ((length (file-length stream)))
	    (incf blocks (ceiling length block-size))
	    (when (< length 100) (incf small))))))
    (values (* blocks block-size) existing small all)))
      
(defun delete-dump (directory time)
  (let ((files (dump-files directory :time time)))
    (loop while files
	  do (do-run-program "rm" :args
			     (cons "-f"	;Ignore missing
				   (loop repeat 100 ;Delete in batches of 100
					 for file = (pop files)
					 while file
					 collect file))))))
    
(defun describe-dump-files (*output-directory* &key time step (terse t))
  (dolist (file (dump-files *output-directory* :time time :step step))
    (format t "~&~A: " file)
    (describe-file file terse)))

(defun describe-successor-files (job-number &optional (terse t))
  (loop for index below 4
	for file = (successor-file index job-number)
	do (format t "~&~A: " file)
	(describe-file file terse)))

;;Wait for dumps to appear and plot them
(defun watch-run (directory &rest plot-args &key sleep &allow-other-keys)
  (loop for step from 0
	for time = (read-dumps directory :step step :wait t)
	do (apply #'plot-time-slice :time time :reuse t plot-args)
	when sleep do (sleep sleep)))

;;Receive junction at start of string.
(defun receive-start-junction (junction)
  (setf (object-handle junction) (rejoining-junction-handle junction)) ;Be able to find it again
  (setq *receive-start-junction* junction))	;Save for later

(defun receive-start-vv-junction (junction)
  (setq *receive-start-junction* junction))	;Save for later

(defun receive-start-initial-junction (junction)
  (if *read-dump-diamonds*
      (progn
	(setq *receive-start-junction* junction)	;Save for later when we know diamond
	(handle-initial-junction junction :left))	;String in predecessor started here, so his future is to L into us
    (progn
      (setq *receive-start-junction* junction)))) 
    
(defun receive-diamond (start left right end loop-tag a-kink-created b-kink-created countup inertp)
   (let ((diamond (make-diamond
		  :start start :left left :right right :end end
		  :tag loop-tag
		  :countup countup
		  :a-kink-created a-kink-created
		  :b-kink-created b-kink-created
		  :inertp inertp
		  :predecessorp t)))
    (when *receive-start-junction*	;If string began with junction, set it up now
      (when (or *receive-previous-diamond* *receive-first-diamond*)
	(error "Should not have a junction and a previous or first diamond at the same time"))
      (setf (junction-right-diamond *receive-start-junction*) diamond) ;Fill in its diamond, now that we know it
      (when (vv-junction-p *receive-start-junction*)		       ;If it is a VV junction, put it in the lattice
	(handle-vv-start-junction *receive-start-junction*))
      (setf (left-rejoining-junction diamond) *receive-start-junction*) ;In either case, attach it to diamond
      (setq *receive-start-junction* nil))
    (when *receive-previous-diamond*			       ;Except first time, when not :deleted
      (receive-link-previous diamond *receive-previous-diamond*)) ;Make links to previous diamond
    (setq *receive-previous-diamond* diamond) ;Remember this for next time
    (unless *receive-first-diamond*
      (setq *receive-first-diamond* diamond) ;Remember first diamond of string
      )
    (incf *receive-count-diamonds*)
    (handle-new-diamond diamond)))


;;Figure out how we are related to the previous diamond and set up links.
(defun receive-link-previous (diamond previous)
  (cond ((eq previous :deleted)
	 (setf (diamond-nw diamond) :deleted))
	((eq previous :BH) ;take into account that now we have BH diamonds
	 (setf (diamond-nw diamond) :BH))
	((eq previous :BHdeleted) ;and BHdeleted ones too
	 (setf (diamond-nw diamond) :BHdeleted))
	((eq previous :BHeatit)
	 (setf (diamond-nw diamond) :BHeatit))
	((eq previous :BHpropdel)
         (setf (diamond-nw diamond) :BHpropdel))
	((eq (diamond-left diamond) (diamond-end previous)) ;Our left is his end, so we are NE of him
	 (setf (diamond-ne previous) diamond
	       (diamond-sw diamond) previous)
	 ;;Propagate a-kink from previous to this diamond.  In the case of a loop, DIAMOND can be the first
	 ;;transmitted diamond in the loop, and so we might need to propagate this kink to further diamonds
	 ;;that were set up before the last diamond in the loop was available.
	 (loop with kink  = (diamond-a-kink-created previous)
	       for d = diamond then (diamond-ne d)
	       while d
	       do (setf (diamond-a-kink-created d) kink)))
	((eq (diamond-end diamond) (diamond-right previous)) ;Our left is his start, so we are SE of him
	 (setf (diamond-se previous) diamond
	       (diamond-nw diamond) previous)
	 (loop with kink  = (diamond-b-kink-created previous) ;See NE case above.
	       for d = diamond then (diamond-se d)
	       while d
	       do (setf (diamond-b-kink-created d) kink)))
	(t
	 (error "Can't figure out how ~S is related to previous ~S" diamond previous))))

;;The first diamond goes after the one we just received to make a loop.
(defun receive-loop ()
  (if *read-dump-diamonds*
      (progn
	(receive-link-previous *receive-first-diamond* *receive-previous-diamond*)
	(start-receive-string))		;Get ready for next
    (progn
      (let ((segment (make-dump-segment)))
	(setf (dump-segment-dump-time segment) *dump-time*)
	(setf (dump-segment-start-junction segment) :loop)
	(setf (dump-segment-end-junction segment) :loop)
	(setf (dump-segment-length segment) *dump-length*)
	(push segment *dumped-segments*)
	(start-receive-segment))
      ))
  )

;; The string starts with a :deleted diamond
(defun receive-start-deleted ()
  (if *read-dump-diamonds*
      (progn
	(start-receive-string)
	(setq *receive-previous-diamond* :deleted))
    (progn
      (setq *receive-start-junction* :deleted))))

;added by SJW
(defun receive-start-BH ()
  (if *read-dump-diamonds*
      (progn
	(start-receive-string)
	(setq *receive-previous-diamond* :BH))
    (progn
      (setq *receive-start-junction* :BH))))


;string that starts with a :BHdeleted diamond 
(defun receive-start-BHdeleted ()
  (if *read-dump-diamonds*
      (progn
	(start-receive-string)
	(setq *receive-previous-diamond* :BHdeleted))
    (progn
      (setq *receive-start-junction* :BHdeleted))))

(defun receive-start-BHeatit ()
  (if *read-dump-diamonds*
      (progn
	(start-receive-string)
	(setq *receive-previous-diamond* :BHeatit))
    (progn
      (setq *receive-start-junction* :BHeatit))))

(defun receive-start-BHpropdel ()
  (if *read-dump-diamonds*
      (progn
	(start-receive-string)
	(setq *receive-previous-diamond* :BHpropdel))
    (progn
      (setq *receive-start-junction* :BHpropdel))))

;; The string ends with a :deleted diamond
(defun receive-end-deleted ()
  (if *read-dump-diamonds*
      (progn
	(setf (diamond-ne *receive-previous-diamond*) :deleted)
	(start-receive-string))
    (progn
      (let ((segment (make-dump-segment)))
        (setf (dump-segment-dump-time segment) *dump-time*)
        (setf (dump-segment-start-junction segment) *receive-start-junction*)
	(setf (dump-segment-start-bh segment) *receive-start-bh*)
        (setf (dump-segment-end-junction segment) :deleted)
        (setf (dump-segment-length segment) *dump-length*)
	(if (rejoining-junction-p *receive-start-junction*)
	    (setf (gethash (rejoining-junction-handle *receive-start-junction*) *start-dumped-segments*) segment)
	  (push segment *dumped-segments*))
	(start-receive-segment)))))

;added by SJW
(defun receive-end-BH ()
  (if *read-dump-diamonds*
      (progn
	(setf (diamond-ne *receive-previous-diamond*) :BH)
	(start-receive-string))
    (progn
      (let ((segment (make-dump-segment)))
        (setf (dump-segment-dump-time segment) *dump-time*)
        (setf (dump-segment-start-junction segment) *receive-start-junction*)
	(if (and *receive-start-bh* *receive-end-bh*)
	    (progn
	      (setf (dump-segment-start-bh segment) *receive-start-bh*)
	      (setf (dump-segment-end-bh segment) *receive-end-bh*))
	  (if (eq *receive-start-junction* :BHeatit)
	      (progn
                (setf (dump-segment-start-bh segment) *receive-start-bh*)
                (setf (dump-segment-end-bh segment) *receive-start-bh*))
	    (setf (dump-segment-end-bh segment) *receive-start-bh*)))
        (setf (dump-segment-end-junction segment) :BH)
        (setf (dump-segment-length segment) *dump-length*)
	(if (rejoining-junction-p *receive-start-junction*)
	    (setf (gethash (rejoining-junction-handle *receive-start-junction*) *start-dumped-segments*) segment)
	  (push segment *dumped-segments*))
	(start-receive-segment)))))


;string that ends with a :BHdeleted diamond
(defun receive-end-BHdeleted ()
  (if *read-dump-diamonds*
      (progn
	(setf (diamond-ne *receive-previous-diamond*) :BHdeleted)
	(start-receive-string))
    (progn
      (let ((segment (make-dump-segment)))
        (setf (dump-segment-dump-time segment) *dump-time*)
        (setf (dump-segment-start-junction segment) *receive-start-junction*)
	(if (and *receive-start-bh* *receive-end-bh*)
            (progn
              (setf (dump-segment-start-bh segment) *receive-start-bh*)
              (setf (dump-segment-end-bh segment) *receive-end-bh*))
          (setf (dump-segment-end-bh segment) *receive-start-bh*))
        (setf (dump-segment-end-junction segment) :BHdeleted)
        (setf (dump-segment-length segment) *dump-length*)
	(if (rejoining-junction-p *receive-start-junction*)
	    (setf (gethash (rejoining-junction-handle *receive-start-junction*) *start-dumped-segments*) segment)
	  (push segment *dumped-segments*))
	(start-receive-segment)))))

(defun receive-end-BHeatit ()
  (if *read-dump-diamonds*
      (progn
	(setf (diamond-ne *receive-previous-diamond*) :BHeatit)
	(start-receive-string))
    (progn
      (let ((segment (make-dump-segment)))
        (setf (dump-segment-dump-time segment) *dump-time*)
        (setf (dump-segment-start-junction segment) *receive-start-junction*)
	(if (and *receive-start-bh* *receive-end-bh*)
            (progn
              (setf (dump-segment-start-bh segment) *receive-start-bh*)
              (setf (dump-segment-end-bh segment) *receive-end-bh*))
	  (if (or (eq *receive-start-junction* :BHeatit)
		  (eq *receive-start-junction* :BH))
	      (progn
		(setf (dump-segment-start-bh segment) *receive-start-bh*)
		(setf (dump-segment-end-bh segment) *receive-start-bh*))
	    (setf (dump-segment-end-bh segment) *receive-start-bh*)))
	(setf (dump-segment-end-junction segment) :BHeatit)
        (setf (dump-segment-length segment) *dump-length*)
	(if (rejoining-junction-p *receive-start-junction*)
	    (setf (gethash (rejoining-junction-handle *receive-start-junction*) *start-dumped-segments*) segment)
	  (push segment *dumped-segments*))
	(start-receive-segment)))))

(defun receive-end-BHpropdel ()
  (if *read-dump-diamonds*
      (progn
	(setf (diamond-ne *receive-previous-diamond*) :BHpropdel)
	(start-receive-string))
    (progn
      (let ((segment (make-dump-segment)))
        (setf (dump-segment-dump-time segment) *dump-time*)
        (setf (dump-segment-start-junction segment) *receive-start-junction*)
	(setf (dump-segment-start-bh segment) *receive-start-bh*)
        (setf (dump-segment-end-junction segment) :BHpropdel)
        (setf (dump-segment-length segment) *dump-length*)
	(if (rejoining-junction-p *receive-start-junction*)
	    (setf (gethash (rejoining-junction-handle *receive-start-junction*) *start-dumped-segments*) segment)
	  (push segment *dumped-segments*))
	(start-receive-segment)))))

;;String ends with junction
(defun receive-end-junction (junction)
  (if *read-dump-diamonds*
      (progn
	(setf (object-handle junction) (rejoining-junction-handle junction)) ;Be able to find it again
	(receive-end-junction-1 junction))
    (progn
      (setf (object-handle junction) (rejoining-junction-handle junction)) ;Be able to find it again  
      (let ((segment (make-dump-segment)))
	(setf (dump-segment-dump-time segment) *dump-time*)
        (setf (dump-segment-start-junction segment) *receive-start-junction*)
	(setf (dump-segment-start-bh segment) *receive-start-bh*)
        (setf (dump-segment-end-junction segment) junction)
        (setf (dump-segment-length segment) *dump-length*)
	(when (rejoining-junction-p *receive-start-junction*)
	  (setf (gethash (rejoining-junction-handle *receive-start-junction*) *start-dumped-segments*) segment))
	(setf (gethash (rejoining-junction-handle junction) *end-dumped-segments*) segment)
	(start-receive-segment))
      )))

    
(defun receive-end-vv-junction (junction)
  (receive-end-junction-1 junction)	;Common code, including putting diamond in junction
  (handle-vv-end-junction junction)	;Put in lattice
  )

(defun receive-end-initial-junction (junction)
  (receive-end-junction-1 junction)	;Common code, including putting diamond in junction
  (handle-initial-junction junction :right)) ;Moving rightward, string in predecessor ends here, to future is to R

;;Common code for end-junctions
(defun receive-end-junction-1 (junction)
   ;;Junction slots are filled in by receive-argument-junction, except for diamond
  (setf (junction-left-diamond junction) *receive-previous-diamond*)
  ;;Say that the string ends in this junction, either VV or regular
  (setf (right-rejoining-junction *receive-previous-diamond*) junction)
  (start-receive-string)		;Get ready for next  
  )

(defun receive-mybh ()
  (setq *receive-mybh* t))


;; Handle the bh information of each diamond
(defun receive-bh (blackhole)
  (if *receive-mybh*
      (progn
	(when (null (gethash (blackhole-center blackhole) *my-bhs*))
	  (setf (gethash (blackhole-center blackhole) *my-bhs*) blackhole))
	(setq *receive-mybh* nil))
    (if *read-dump-diamonds*
	(setf (diamond-bh *receive-previous-diamond*) blackhole)
      (if *receive-start-bh*
	  (setf *receive-end-bh* blackhole)
	(setf *receive-start-bh* blackhole)))))

;;Handle the end of the string: get ready for the next
(defun start-receive-string ()
  (setq *receive-first-diamond* nil
	*receive-previous-diamond* nil
	*receive-start-junction* nil))

(defun start-receive-segment ()
  (setq *receive-start-junction* nil
	*receive-start-bh* nil
	*receive-end-bh* nil
	*dump-length* 0.0))

  
(mirror-images
(defun receive-note-left-junction (junction)
  ;;We don't set the junction-right-diamond for this dump junction.
  ;;I think this is right, because this type of junction note is never reassembled in our process
  ;;Instead the junction is dumped again as a start or end junction in a dump file and reassembled on reading that
  (push junction (received-left-dump-junctions *receive-previous-diamond*)))
)

(defun receive-dump-time (time)
  (unless *dump-time* (error "Received dump-time, but not reading dumps"))
  (if (eq *dump-time* :unknown)		;This allows reading a dump without knowing what time it was made
      (setq *dump-time* time *reading-dumps* time)
    (unless (fudge= time *dump-time* fudge-global-coordinates)
      (error "This dump is for time ~F, but we expected ~F" time *dump-time*))))

(defun receive-dump-length (length)
  (if *read-dump-diamonds*
      (progn
	(unless *dump-time* (error "Received dump-length, but not reading dumps"))
	(when (numberp *dump-length*)
	  (incf *dump-length* length)))	;Accumulate total length read
    (when (numberp *dump-length*)
      (setf *dump-length* length)))
  )

(defun receive-tag (handle created-position last-position xi minimum-inert-count dumped bh)
  (let ((old (handle-object handle nil))) ;See if we already received the same tag from a different predecessor
    (cond (old				;Yes.  Potentially update old tag, then use it.
	   (when bh
	     (when (null (tag-bh old))
	       (error "We have to update the bh in tag")))
	   (unless (4vector-exactly= created-position (tag-created-position old)) ;Never updated.  Should be exact.
	     (error "Received same tag with a different created position"))
	   (unless (zerop xi)		;If xi set by current predecessor, might need propagation to old tag
	     (if (zerop (tag-xi old))	;but not previous one
		 (setf (tag-xi old) xi)	;install
	       (unless (= xi (tag-xi old)) ;both set: must match
		 (error "Tag x_i different from different predecessors"))))
	   ;;Minimum your account to delete his larger of minima from different sides
	   (setf (tag-minimum-inert-count old) (max (tag-minimum-inert-count old) minimum-inert-count))
	   ;;If either predecessor dumped the loop (i.e. a part with count >= 2), then is done.
	   (setf (tag-dumped old) (or (tag-dumped old) dumped))
	   (when last-position		;Have position of last engulfment?
	     (if (tag-last-position old) ;Position of last engulfment from other side also?
		 (when (> (4vector-t last-position) (4vector-t (tag-last-position old))) ;Use later
		   (setf (tag-last-position old) last-position))
	       (setf (tag-last-position old) last-position)))) ;Set new position in old tag
	  (t				;New to us.  Create object.
	   (let ((tag (make-loop-tag :handle handle :created-position created-position
				     :last-position last-position
				     :xi xi
				     :minimum-inert-count minimum-inert-count
				     :dumped dumped
				     :bh bh)))
	     (setf (object-handle tag) handle)))	;Install in hash table
	  )))
	  
;;If we have both the left and right junctions for a diamond, put them together
(defun connect-junctions ()
  (let ((unconnected 0))
    (map-junctions
     #'(lambda (junctions)			      ;List of junctions with the same handle
	 (cond ((null (cdr junctions))		      ;Only one junction with this handle
;		(format t "Jucs: ~S~%" junctions)
		(incf unconnected)		      ;Nothing to do.  Pass it on later.
		(when (> (job-layer *job-number*) ;If we aren't rejoining, it should be in the previous layer
			 (1+ (job-layer (handle-creator (rejoining-junction-handle (car junctions))))))
		  (error "Stale unrejoined junction ~S" (car junctions))))
	       ((null (cddr junctions))		      ;Exactly 2 junctions
		(let ((left-junction (first junctions)) ;Guess that the first one is the left end of the string
		      (right-junction (second junctions)))
		  (unless (junction-right-diamond left-junction)              ;Wrong guess
		    (rotatef left-junction right-junction))		      ;Exchange them
		  (let ((left-diamond (junction-left-diamond right-junction)) ;Diamond to left of junction
			(right-diamond (junction-right-diamond left-junction))) ;Diamond to right of junction
		    (mirror-images (assert left-diamond)			;Make sure there is no confusion
				   (assert (null (junction-right-diamond right-junction))))
		    (if ;;These are either the same diamond or they are images of the same diamond.
			;;Only in the former case do we merge them, except when reading dumps
			(or *reading-dumps*
			    (< (3vector-distance (diamond-start left-diamond) (diamond-start right-diamond))
			       (+ diamond-span fudge-coordinates))) ;Once part of same diamond, so this close
			(merge-diamonds right-junction left-junction left-diamond right-diamond)
		      (incf unconnected)))))
	       (t
		(error "Too many junctions with one handle: ~S" junctions)))))
    (when (plusp unconnected)
      (format t "~&~D unconnected junction~:P~%" unconnected))))

(defun connect-junctions-segments ()
  (let ((unconnected 0)
	(tot 0.0d0))
    (map-junctions
     #'(lambda (junctions)
	 (cond ((null (cdr junctions))
		(incf unconnected)
		(when (> (job-layer *job-number*)
			 (1+ (job-layer (handle-creator (rejoining-junction-handle (car junctions))))))
		  (error "Stale unrejoined junction ~S" (car junctions))))
	       ((null (cddr junctions))
		(merge-segments (rejoining-junction-handle (first junctions)) (rejoining-junction-handle (second junctions)))))))
    (merge-loops)
    (loop for k being each hash-key of *start-dumped-segments*
	  do (setf tot (+ tot (dump-segment-length (gethash k *start-dumped-segments*)))) 
	  do (format t "Seg: ~S~%" tot))
    (when (or (> (hash-table-count *start-dumped-segments*) 0)
	      (> (hash-table-count *end-dumped-segments*) 0))
      (error "NOT empty hash~%"))))

  

(defun merge-segments (handle1 handle2)
  (let ((segSt nil)
	(segEn nil)
	(segment (make-dump-segment))
        (dump-time nil)
        (start-junction nil)
        (start-bh nil)
        (end-junction nil)
        (end-bh nil)
	(handle11 nil)
	(handle21 nil)
        (length 0.0d0))
    (when (null (equalp handle1 handle2))
      (error "Different handles in merging junctions"))
    (setf segSt (gethash handle1 *start-dumped-segments*))
    (when (rejoining-junction-p (dump-segment-end-junction segSt))
      (setf handle11 (rejoining-junction-handle (dump-segment-end-junction segSt))))
    (setf segEn (gethash handle2 *end-dumped-segments*))
    (when (rejoining-junction-p (dump-segment-start-junction segEn))
      (setf handle21 (rejoining-junction-handle (dump-segment-start-junction segEn))))
    
    (if (eq (dump-segment-dump-time segSt) (dump-segment-dump-time segEn))
        (setf dump-time (dump-segment-dump-time segSt))
      (error "Junctions with different Dump-time"))
    (if (equalp handle1 (rejoining-junction-handle (dump-segment-start-junction segSt)))
        (progn
          (setf end-junction (dump-segment-end-junction segSt))
          (setf end-bh (dump-segment-end-bh segSt)))
      (error "This is am end junction"))
    (if (equalp handle2 (rejoining-junction-handle (dump-segment-end-junction segEn)))
        (progn
          (setf start-junction (dump-segment-start-junction segEn))
          (setf start-bh (dump-segment-start-bh segEn)))
      (error "This is a start junction"))

    (if (and (rejoining-junction-p (dump-segment-end-junction segSt)) (equalp handle1 (rejoining-junction-handle (dump-segment-end-junction segSt))))
	(setf length (dump-segment-length segSt))
      (setf length (+ (dump-segment-length segSt) (dump-segment-length segEn))))
    
    ;(format t "~S + ~S = ~S~%" (dump-segment-length segSt) (dump-segment-length segEn) length)

        
    (setf (dump-segment-dump-time segment) dump-time)
    (setf (dump-segment-start-junction segment) start-junction)
    (setf (dump-segment-start-bh segment) start-bh)
    (setf (dump-segment-end-junction segment) end-junction)
    (setf (dump-segment-end-bh segment) end-bh)
    (setf (dump-segment-length segment) length)

    
    (remhash handle1 *start-dumped-segments*)
    (when handle11
      (remhash handle11 *end-dumped-segments*))
    (remhash handle2 *end-dumped-segments*)
    (when handle21
      (remhash handle21 *start-dumped-segments*))

    (when (rejoining-junction-p start-junction)
      (setf (gethash (rejoining-junction-handle start-junction) *start-dumped-segments*) segment))
    (when (rejoining-junction-p end-junction)
      (setf (gethash (rejoining-junction-handle end-junction) *end-dumped-segments*) segment))

    (when (and (null (rejoining-junction-p start-junction))
	       (null (rejoining-junction-p end-junction)))
      (push segment *dumped-segments*))
   ))

(defun merge-loops ()
  (let ((seg nil)
	(start-junction nil)
	(end-junction nil)
	(start-handle nil)
	(end-handle nil)
	(tot 0.0d0))
  (loop for k being each hash-key of *start-dumped-segments*
	do (setf seg (gethash k *start-dumped-segments*))
	do (setf start-junction (dump-segment-start-junction seg))
	do (setf start-handle (rejoining-junction-handle start-junction))
	do (setf end-junction (dump-segment-end-junction seg))
	do (setf end-handle (rejoining-junction-handle end-junction))
	do (when (equalp start-handle end-handle)
	     (setf (dump-segment-start-junction seg) :loop)
	     (setf (dump-segment-end-junction seg) :loop)
	     (setf tot (+ tot (dump-segment-length seg)))
	     (remhash k *start-dumped-segments*)
	     (remhash k *end-dumped-segments*)
	     (push seg *dumped-segments*)))
  (format t "TOT: ~S~%" tot)))

(defun merge-segments-old (seg1 seg2 junc1 junc2)
  (let ((segment (make-dump-segment))
	(dump-time nil)
	(start-junction nil)
	(start-bh nil)
	(end-junction nil)
	(end-bh nil)
	(length 0.0d0))
    (if (eq (dump-segment-dump-time seg1) (dump-segment-dump-time seg2))
	(setf dump-time (dump-segment-dump-time seg1))
      (error "Junctions with different Dump-time"))
    (if (eq junc1 (dump-segment-start-junction seg1))
	(progn
	  (setf start-junction (dump-segment-end-junction seg1))
	  (setf start-bh (dump-segment-end-bh seg1)))
      (progn
	(setf start-junction (dump-segment-start-junction seg1))
	(setf start-bh (dump-segment-start-bh seg1))))
    (if (eq junc2 (dump-segment-start-junction seg2))
	(progn
	  (setf end-junction (dump-segment-end-junction seg2))
	  (setf end-bh (dump-segment-end-bh seg2)))
      (progn
	(setf end-junction (dump-segment-start-junction seg2))
	(setf end-bh (dump-segment-start-bh seg2))))

    (setf length (+ (dump-segment-length seg1) (dump-segment-length seg2)))

    (setf (dump-segment-dump-time segment) dump-time)
    (setf (dump-segment-start-junction segment) start-junction)
    (setf (dump-segment-start-bh segment) start-bh)
    (setf (dump-segment-end-junction segment) end-junction)
    (setf (dump-segment-end-bh segment) end-bh)
    (setf (dump-segment-length segment) length)

    segment))
    
;;These two diamonds are really the same diamond, but they may have been cut by predecessors,
;;so we have two subsets of the same diamond.  We want the subset that is in both.
;;The left junction is the one that comes at the left end of the string.
;;The left diamond is the one that is to the left of the right-junction.  That means its left corner might have been
;;cut off by the predecessor.  Similarly, the right corner of the right diamond might have been cut off.
;;In this case, the left diamond extends further to the right and the right diamond further to the left.
;;Cuts will need to be propagated.
(defun merge-diamonds (right-junction left-junction left-diamond right-diamond)
   (let* ((left-a-j (rejoining-junction-a right-junction)) ;left diamond's a-coordinate for the junction
	 (left-b-j (rejoining-junction-b right-junction)) ;
	 (right-a-j (rejoining-junction-a left-junction)) ;right diamond's a-coordinate for the junction
	 (right-b-j (rejoining-junction-b left-junction))
	 (original-a-j (rejoining-junction-a0 left-junction)) ;original diamond's a coordinate for the junction
	 (original-b-j (rejoining-junction-b0 left-junction))
	 (left-a-right (if (= right-a-j original-a-j) 0.0 ;certain to give precisely 0.0 if aR-j = a0-j
		      (/ (* left-a-j (- right-a-j original-a-j)) ;a of right corner of overlap as measured by left-diamond
			 (* original-a-j (- right-a-j 1.0)))))
	 (right-a-left (if (= left-a-j original-a-j) 1.0	;certain to give precisely 1.0 if appropriate
		     (/ (- (* original-a-j (+ left-a-j right-a-j -1.0))
			   (* left-a-j right-a-j))
			(* left-a-j (- original-a-j 1.0)))))
	 (right-b-left (if (= left-b-j original-b-j) 0.0 ;certain to give precisely 0.0 if appropriate
		     (/ (* right-b-j (- left-b-j original-b-j)) ;left corner of overlap as measured by Right-diamond
			(* original-b-j (- left-b-j 1.0)))))
	 (left-b-right (if (= right-b-j original-b-j) 1.0	;certain to give precisely 1.0 if appropriate
		      (/ (- (* original-b-j (+ right-b-j left-b-j -1.0))
			    (* right-b-j left-b-j))
			 (* right-b-j (- original-b-j 1.0)))))
	 (overlap (overlap-diamond right-junction left-junction left-diamond right-diamond)))
     (clear-intersection left-diamond right-diamond) ;don't let any bad intersection through
    (mirror-images
     (divide-pending-intersections left-diamond :a left-a-right :b left-b-right :east nil :west overlap)
     (discard-rejoining-junction right-junction) ;Get rid of rejoining junction since we are joining it now.  Don't rescale it.
     (rescale-right-junctions left-diamond overlap :a left-a-right :b left-b-right)
     (rescale-left-junctions left-diamond overlap :a left-a-right :b left-b-right)
     (clear-links left-diamond)
     (discard-object left-diamond)
     )					;mirror images
    (handle-new-diamond overlap :predecessor :meiosis) ;don't check it for intersections
    ))

;; This creates the overlap region of the two diamonds, and propagates the cuts
;; to make it compatible with all neighbors.  Nothing is discarded. Links and points are put in.
;; It always makes a new diamond, even if neither original diamond had been cut.
;; Note that right-diamond extends to the left of left-diamond 
(defun overlap-diamond (right-junction left-junction left-diamond right-diamond)
  (mirror-images (assert (diamond-predecessorp left-diamond))) ;Sanity checking
  (assert (= (diamond-countup left-diamond) (diamond-countup right-diamond)))
  (assert (eq (diamond-tag left-diamond) (diamond-tag right-diamond)))
  (let* ((aL-j (rejoining-junction-a right-junction)) ;left diamond's (a,b) for rejoining junction location
	 (bL-j (rejoining-junction-b right-junction))
	 (aR-j (rejoining-junction-a left-junction))
	 (bR-j (rejoining-junction-b left-junction))
	 (a0-j (rejoining-junction-a0 left-junction))
	 (b0-j (rejoining-junction-b0 left-junction))
	 (left-a-start (if (= aR-j a0-j) 0.0 ;certain to give precisely 0.0 if aR-j = a0-j
			 (/ (* aL-j (- aR-j a0-j))
			    (* a0-j (- aR-j 1.0)))))
	 (right-a-end (if (= aL-j a0-j) 1.0	;certain to give precisely 1.0 if aL-j = a0-j
			(/ (- (* a0-j (+ aL-j aR-j -1.0))
			      (* aL-j aR-j))
			   (* aL-j (- a0-j 1.0)))))
	 (left-b-end (if (= bR-j b0-j) 1.0	;certain to give precisely 1.0 if bR-j = b0-j
		       (/ (- (* b0-j (+ bR-j bL-j -1.0))
			     (* bR-j bL-j))
			  (* bR-j (- b0-j 1.0)))))
	 (left (diamond-left left-diamond))
	 (right (diamond-right right-diamond))
					;diamond-position knows to check corners as special case
	 (start (diamond-position-wrap-dumps left-diamond :a left-a-start :b 0.0))
	 (end (diamond-position-wrap-dumps left-diamond :a 1.0 :b left-b-end)) 
	 (sw (diamond-sw left-diamond))
	 (nw (diamond-nw left-diamond))
	 (se (diamond-se right-diamond))
	 (ne (diamond-ne right-diamond))
	 (overlap (make-diamond :start start
				:end end
				:left left
				:right right
				:sw sw
				:nw nw
				:se se
				:ne ne
				;;Propagate kinks.  I'm not sure the choice really matters, because if there
				;;was an intersection, we already have a successor for the given kink, and
				;;if not it is the same from the two sides.
				:a-kink-created (diamond-a-kink-created right-diamond)
				:b-kink-created (diamond-b-kink-created left-diamond)
				:predecessorp t
				:tag (diamond-tag left-diamond)
				:countup (diamond-countup left-diamond)
				:inertp (diamond-inertp left-diamond))))
    ;;Check if diamonds have a bh linked and if it is the same. At the moment if the bh are not the same error it
    (when (or (diamond-bh left-diamond) (diamond-bh right-diamond))
      (if (and (diamond-bh left-diamond) (diamond-bh right-diamond))
	  (if (4vector= (blackhole-center (diamond-bh left-diamond)) (blackhole-center (diamond-bh right-diamond)) fudge-global-coordinates)
	      (setf (diamond-bh overlap) (diamond-bh left-diamond))
	    (format t "Left bh: ~S~% and Right bh: ~S~%" (diamond-bh left-diamond) (diamond-bh right-diamond)))
	(if (diamond-bh left-diamond)
	    (setf (diamond-bh overlap) (diamond-bh left-diamond))
	  (setf (diamond-bh overlap) (diamond-bh right-diamond)))))
	    
    ;;Sometimes junctions are not as we expect.  See email from Ken on "Rejoining dump junctions"
    (unless (<= 0.0 left-a-start 1.0) (warn "left-a-start = ~S out of range" left-a-start))
    (unless (<= 0.0 right-a-end 1.0) (warn "right-a-end = ~S out of range" right-a-end))
    (unless (<= 0.0 left-b-end 1.0) (warn "left-b-end = ~S out of range" left-b-end))
    (let ((left-position (diamond-position-wrap-dumps left-diamond :a aL-j :b bL-j))
	  (right-position (diamond-position-wrap-dumps right-diamond :a aR-j :b bR-j)))
      (unless (4vector= left-position right-position fudge-global-coordinates)
	(warn "Identical junction positions differ by ~D" (4vector- left-position right-position))))
    (mirror-images
     (cond ((diamondp nw) 
	    (setf (diamond-se nw) overlap) ;nw may be :deleted
	    (when (4vector= (diamond-start nw) ;nw may have been cut
			    left
			    fudge-coordinates)
	      (setf (diamond-left overlap) (diamond-start nw))))      ;make left point eq to left neighbor's start (not vice versa)
	   (sw 
	    (setf (diamond-ne sw) overlap)
	    (assert (eq (diamond-end sw) left)) ;check left neighbor's start eq to left point
	    (when (4vector= (diamond-start overlap) (diamond-right sw) fudge-coordinates) 
	      (setf (diamond-right sw) (diamond-start overlap)))) ;make it eq to start (not vice versa)
	   (t 
	    (assert (or (eq nw :deleted) (eq nw :BH) (eq nw :BHdeleted) (eq nw :BHeatit) (eq nw :BHpropdel) ;BH added
			(junction-p (left-rejoining-junction left-diamond))))))   ;there must be a proper left-junction in left-diamond
     (if (diamondp nw) 
	 (if (= left-b-end 1.0)
	     (setf (diamond-right nw) end)
	   (propagate-cut-nw overlap left-b-end)))
     )					;mirror-images
    overlap)				;output the diamond
  )


;; removes the links to neighbor diamonds, making it safe to discard this one
;; without affecting others.
(defun clear-links (diamond)
  (mirror-images
   (setf (diamond-nw diamond) nil)
   (setf (diamond-sw diamond) nil)
   )					;mirror images
  )

;; Removes merging diamonds from each other's pending intersection list, since they
;; may have been found to intersect
(defun clear-intersection (left-diamond right-diamond)
  (let ((pending (diamond-pending-intersections left-diamond))
	)
    (dolist (intersection pending)
      (when (or (eq (intersection-diamond-1 intersection) right-diamond)
		(eq (intersection-diamond-2 intersection) right-diamond))
	(discard-object intersection)))))

;;Read the information about the run
;;Handle old versions of the file format that used dump-start, dump-interval, dump-only-length
(defun read-run-info-file (&optional (directory (input-directory)))
  (let* ((info (with-open-file (stream (run-info-file directory))
		 (let ((*read-default-float-format* 'double-float))
		   (read stream))))
	 (start (run-info-dump-start info))
	 (end (run-info-end info))
	 (interval (run-info-dump-interval info)))
    (when interval
      (unless end
	(warn "No end-time in run-info file, assuming end=size")
	(setq end (run-info-total-size info)))
      (let ((times (loop for time from start by interval below (+ end fudge-global-coordinates)
			 unless (< time (run-info-start-time info)) ;Too early: was not made
			 collect time)))
	(setf (run-info-length-times info) times)
	(unless (run-info-dump-only-length info)
	  (setf (run-info-dump-times info) times))))
    info))

; interval
;; Gather final energy from output file (use conformal units)
(defun read-final-energy (directory)
  (loop for file in (directory (format nil "~A/rank-*/output" directory))
	sum (with-open-file (stream file)
	 		    (or (loop for line = (read-line stream nil)
				      while line
				      when (and (> (length line) 38)
						(string-equal line " Total energy in non-inert remains is " :end1 38))
				      return (read-from-string line nil nil :start 38 :end (- (length line) 1)))
				(error "unable to read energy from ~A" file))))) 

;;Find longest edge from diamonds read in and install in *longest-edge*
(defun setup-longest-edge ()
  (let ((longest nil))
    (map-read-diamonds
     #'(lambda (diamond)
	 (let ((this (max (3vector-length (diamond-p diamond)) (3vector-length (diamond-q diamond)))))
	   (setq longest (if longest (max longest this) this)))))
    (if longest
	(setq *longest-edge* longest
	      monster-length-scale (* monster-length-multiplier *longest-edge*))
      (warn "Simulation volume has no string"))))
