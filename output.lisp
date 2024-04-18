;;;Output for successors and snapshots
(in-package "CL-USER")

;;See notes.text for file system structure

;;Top-level directory for output
(define-simulate-variable *output-directory* "test")

;;Directory to read input files.  Defaults to *output-directory*
(define-simulate-variable *input-directory* nil)

(defun input-directory ()
  (or *input-directory* *output-directory*))

;;If set, this is the output directory on the server rather than the local scratch directory
(defvar *global-output-directory* nil)
(defun global-output-directory ()
  (or *global-output-directory* *output-directory*))

;;Construct filename for job given top-level directory and format controls for part after hyphen
(defun job-filename (directory job-number format-string &rest format-args)
  (multiple-value-bind (quotient remainder) (truncate job-number 1000)
    (format nil "~A/~D/~D-~?" directory quotient remainder format-string format-args)))

(defun worker-subdirectory (directory worker-number)
  (format nil "~A/worker-~A" directory worker-number))

;;Default place to store random seed for this processor
(defun random-seed-file (directory &optional (job-number *job-number*))
  (job-filename directory job-number "random-seed.lisp"))

(defun successor-file (successor &optional (job-number *job-number*))
  (job-filename *output-directory* job-number "output-~A.dat" (successor-name successor)))

(defun successor-files ()
  (loop for successor below 4 collect (successor-file successor)))

(defvar *predecessor-files-copied* nil)	;T if files have already been copied into input-....

;;Location of predecessor file.  If COPIED, input-...
;;If LOCAL, but not COPIED, local scratch directory of originating process
;;Otherwise global directory.
(defun predecessor-file (predecessor &key (copied *predecessor-files-copied*) (of-job-number *job-number*)
				     (local (and *worker-directory* ;If not worker, always on server
						 *local-data-files*))
				     user
				     )
  (multiple-value-bind (split direction) (floor predecessor 4) ;In case multiple sets of predecessors
    (if copied				;Already copied to input-...
	(job-filename (input-directory) of-job-number "input-~A.dat" (successor-name direction))
      (let ((predecessor-number (predecessor-job-number direction of-job-number))
	    (file (format nil "output-~A.dat" (successor-name direction)))
	    (split-dir (if (> (length *predecessor-hosts*) 4) (format nil "/splits/~D" split) "")))
	(cond (local			;Some other node local scratch
	       (unless (= of-job-number *job-number*)
		 (error "local data files by :nfs incompatible with :of-job-number"))
	       (job-filename (format nil "~A~A~A"
				     (remote-local-root-directory (nth predecessor *predecessor-hosts*) :user user)
				     *worker-directory* split-dir)
			     predecessor-number file))
	      (t (job-filename (format nil "~A~A" (input-directory) split-dir)
			       predecessor-number file)))))))

;;Files for dump output
(defun dump-directory (job-number)
  (job-filename (global-output-directory) job-number "dump"))

(defun dump-file (job-number step)
  (format nil "~A/step-~A.dat" (dump-directory job-number) step))

;;Temporary file to write dump output.  Renamed to above later.
(defun temporary-dump-file (&optional job-number step)
  (format nil "~A/step-~A.tmp" (dump-directory job-number) step))

;;All dump files, including temporary
(defun all-dump-files (job-number)
  (job-filename *output-directory* job-number "step-*.*"))

;;Regular place to store loop log file.
(defun loop-spectrum-file (directory worker-number)
  (format nil "~A/loops.dat" (worker-subdirectory directory worker-number)))

;;Regular place to store bh loop log file.
(defun bh-loop-spectrum-file (directory worker-number)
  (format nil "~A/bhloops.dat" (worker-subdirectory directory worker-number)))


;;Prebin data or info file.
;;Setting mass give "alpha" as part of filename, but
;;that is not implemented anywhere else.
(defun prebin-file-name (base-filename directory &key mass integrated worker)
  (format nil "~A/~:[~;alpha-~]~:[~;sum-~]~A"
	  (if worker (worker-subdirectory directory worker)
	    directory)
	  mass integrated
	  base-filename))

(defun prebin-data-file (directory &rest keys)
  (apply #'prebin-file-name "loops.prebin" directory keys))

(defun prebin-bh-data-file (directory &rest keys)
  (apply #'prebin-file-name "bhloops.prebin" directory keys))

(defun prebin-info-file (directory &rest keys)
  (apply #'prebin-file-name "prebin-info.lisp" directory keys))


;;Place to store log of lengths
(defun length-file (directory &optional worker-number)
  (format nil "~A/lengths.dat" (if worker-number (worker-subdirectory directory worker-number)
				 directory)))

;;Place to store log of information propagation distance
(defun information-flow-distance-file (directory &optional worker-number)
  (format nil "~A/info-flow-distances.dat" (if worker-number (worker-subdirectory directory worker-number)
				 directory)))

;;Per-job place for testing
(defun job-loop-spectrum-file (&optional (job-number *job-number*))
  (job-filename *output-directory* job-number "loops.dat"))

;;Per-job place for testing bh loop output
(defun job-bh-loop-spectrum-file (&optional (job-number *job-number*))
  (job-filename *output-directory* job-number "bhloops.dat"))

(defun run-info-file (&optional (directory (global-output-directory)))
  (format nil "~A/run-info.lisp" directory))

(defun two-point-file (time &optional (directory (global-output-directory)))
  (format nil "~A/two-point-~A.lisp" directory (format-float-reasonably time)))

(defun power-file (time &optional (directory (global-output-directory)))
  (format nil "~A/power-~A.lisp" directory (format-float-reasonably time)))

;;This attempts to cure the problem that floating-point numbers which are supposed to be integers or
;;numbers with a few decimal digits acquire different low-order bits, which lisp then attempts to print
(defun format-float-reasonably (x)
  (let ((string (format nil "~15F" x)))
    (subseq string
	    (loop for index from 0
		  while (char= (char string index) #\space)
		  finally (return index)))))

;;Format seconds as h:mm:ss.                                                                                                                                                                             
(defun format-seconds (seconds)
  (format nil "~D:~2,'0D:~2,'0D" (truncate seconds 3600) (truncate (mod seconds 3600) 60) (mod seconds 60)))

;;Find A and B where line from (a1, b1) to (a2, b2) intersects coordinate = constant surface.
;;Treats diamond as flat.
;;Arguments as keywords to help with mirror-images
;;If ERROR is NIL, we return NIL if the given edge does not overlap the given time.  Otherwise
;;we signal an error
(defun interpolate-find-coordinate (diamond coordinate value &key a1 b1 a2 b2 (error t))
  (let* ((value1 (aref (diamond-position diamond :a a1 :b b1) coordinate))
	 (value2 (aref (diamond-position diamond :a a2 :b b2) coordinate))
	 (parameter (/ (- value value1) (- value2 value1)))) ;interpolation parameter 0=(a1, b1), 1=(a2, b2)
    (cond ((<= 0 parameter 1)	      ;Normal case.
	   (values (+ (* a1  (- 1 parameter)) (* a2 parameter)) ;interpolate A and B
		   (+ (* b1  (- 1 parameter)) (* b2 parameter))))
	  ((and (minusp parameter)
		(< (- value1 value) fudge-coordinates)) ;Smaller than value1, but within fudge factor
	   (values a1 b1))
	  ((and (> parameter 1)
		(< (- value value2) fudge-coordinates)) ;Larger than value2, but within fudge factor
	   (values a2 b2))
	  (error
	   (error "The value ~F of coordinate ~D is not included in the given range" value coordinate))
	  (t nil))))

  
(mirror-images
;;Return A and B in DIAMOND giving the point where the right edge of the diamond has the given TIME
;;Both mirror image versions of this function return A, B.
;;If ERROR is NIL, we return NIL if the time does not lie on the edge.  Otherwise give an error.
 (defun find-right-edge-position (diamond time &optional (error t))
   (if (or (eq (diamond-e diamond) :BH)
	   (eq (diamond-e diamond) :BHdeleted))
       (if (> (4vector-t (diamond-position diamond :a 0.5 :b 0.5)) time) ;Before right corner?
	   (interpolate-find-coordinate diamond 3 time :a1 0.0 :b1 0.0 :a2 0.5 :b2 0.5 :error error) ;From (0,0) to (0,1)
	 (interpolate-find-coordinate diamond 3 time :a1 0.5 :b1 0.5 :a2 1.0 :b2 1.0 :error error)) ;From (0,1) to (1,1) 
     (if (> (4vector-t (diamond-right diamond)) time) ;Before right corner?
	 (interpolate-find-coordinate diamond 3 time :a1 0.0 :b1 0.0 :a2 0.0 :b2 1.0 :error error) ;From (0,0) to (0,1)
       ;;After right corner
       (interpolate-find-coordinate diamond 3 time :a1 0.0 :b1 1.0 :a2 1.0 :b2 1.0 :error error) ;From (0,1) to (1,1)
       )))
)

;;If T, get-path does not complain if the string comes to an end, but instead gives those segments that are valid
(defvar *allow-get-path-end* t)
(defvar *allow-get-path-deleted* nil)

;;Accept arguments from MAP-STRING-PATHS specifying a string segment.
;;Returns list of (DIAMOND A1 B1 A2 B2).
;;Return NIL if segment has :deleted at either end unless *get-path-deleted* set
(defun get-path (time start &optional start-junction end-diamond end-junction)
  (unless (and (or (eq start-junction :deleted) (eq end-junction :deleted))
	       (not *allow-get-path-deleted*))
    (loop with diamond = start
	  for left = (if (and start-junction (rejoining-junction-p start-junction))
			 (list (rejoining-junction-a start-junction) (rejoining-junction-b start-junction))
		       (multiple-value-list
			(find-left-edge-position diamond time (not *allow-get-path-end*)))) ;starting position
	  for right = (if (and (eq diamond end-diamond) end-junction (rejoining-junction-p end-junction))      ;If this one is last
			  (list (rejoining-junction-a end-junction) ;know ending position
				(rejoining-junction-b end-junction))
			(multiple-value-list
			 (find-right-edge-position diamond time (not *allow-get-path-end*)))) ;otherwise right edge
	  when (and (car right) (car left))
	  collect (cons diamond (nconc left right))
	  until (eq diamond end-diamond) ;If we have just collected end diamond, then all done
	  do (setq diamond (diamond-e diamond)
		   start-junction nil)	;Forget this if we started with it
	  until (and (eq diamond start) (not end-diamond)) ;If a loop, don't do start diamond again
	  )))

;;Get list of diamonds without worrying about time.
(defun get-string-diamonds (start &optional start-junction end-diamond end-junction)
  (declare (ignore start-junction end-junction))
  (loop with diamond = start
	collect diamond
	until (eq diamond end-diamond)
	do (setq diamond (diamond-e diamond))
	until (and (eq diamond start) (not end-diamond))))

  
(declaim (inline check-diamond-step))
;;Check the diamond that we got in scanning for the result of data structure problems
(defun check-diamond-step (diamond direction)
  (unless diamond
    (error "Scanned ~(~A~) off the end of a string.  This should not happen." direction))
  (when (eq diamond :deleted)
    (error "Scanned ~(~A~) into deleted portion of a string.  This should not happen." direction))
  )


;;Call function for all diamonds that we know about, which means:
;; Diamonds on the calendar
;; Final diamonds that are waiting to be output
;; Diamonds that have been read from dumps (If there are these, there are not others.)
;;In this order for things to work properly, it's important that this function really find all the diamonds,
;;i.e. that the diamonds we find don't link to others that we miss.
;;We do not guarantee to avoid duplicate calls.
(defun map-all-diamonds (function)
  (flet ((do-diamond (diamond) (funcall function diamond)))
    (cond (*reading-dumps*			;Reading dumps, not doing simulation
	   (map-read-diamonds #'do-diamond))	;Just do those those read in from dumps.
	  (t (map-main-calendar
	      #'(lambda (time diamond)	;Find any diamond on calendar
		  (declare (ignore time))     
		  (when (diamondp diamond) ;Ignore things that are not diamonds
		    (do-diamond diamond))))
	     (map-final-diamonds #'do-diamond) ;Now do those those awaiting transmission
	     ))))

;;Call LOOP-FUNCTION once for each loop with some diamond in the loop
;;and OPEN-FUNCTION once for each open string with arguments:
;; START-DIAMOND, START-JUNCTION, END-DIAMOND, END-JUNCTION
;;If START-DIAMOND and END-DIAMOND are the same, either START-JUNCTION is to the left of END-JUNCTION, in which case
;;this is the only diamond or START-JUNCTION is to the right, in which case there is a whole string between.
;;If LOOP-FUNCTION is not given, call the same function, but with different arguments
;;If *DUMP-TIME* is set, do only those strings that should be dumped.  See notes.
(defun map-string-paths (open-function &optional (loop-function open-function))
  (using-diamond-processed
   (unwind-protect
       (map-all-diamonds #'(lambda (diamond) (map-string-paths-1 diamond open-function loop-function)))
     (map-all-diamonds	 ;Now go back and clear processed bits
      #'(lambda (diamond) 
	  (setf (diamond-processedp diamond) nil))))))

;;Here with every diamond that we know about
;;If any strings start in the diamond, we handle those.  If the diamond is part of a loop, we handle that.
;;The processedp bit means that any string starting in the diamond or any loop including the diamond is done
;;We also set it as we look for a loop, but if we don't find one, the open string including this segment is not done
(defun map-string-paths-1 (diamond open-function loop-function)
  (unless (diamond-processedp diamond)	;If set, any string starting here or loop including this has been done.  Exit.
    (let ((starts (diamond-left-junctions diamond)))
      (cond (starts ;Any strings start here?
	     (dolist (start starts)
	       (map-string-paths-open diamond start open-function)) ;Handle them
	     ;;This diamond is all done.  Mark to protect against duplicate call with same diamond.
	     (setf (diamond-processedp diamond) t))
	    (t
	     (map-string-paths-maybe-loop diamond loop-function)) ;Might be loop: handle that
	    ))))

;;Handle one string that is not a loop.  This function is called once for each such string.  We ignore
;;processedp, which means only that later diamonds might have been considered for loops.
;;We mark the diamonds as we go, including the start diamond, unless the string is only one diamond long.
;;We don't mark the end diamond, because more strings might start there, except in the "loopback" case,
;;where we marked the same diamond when we left it earlier.
;;If a string has no dump junction, then we forget it unless it overlaps the dump surface.  See notes.
(defun map-string-paths-open (start-diamond start-junction function)
  (let ((d start-diamond)		;Get ready to scan eastward
	(junction start-junction)
	(must-output (or (null *dump-time*) ;If not dumping, we output everything
			 (if (rejoining-junction-p start-junction) ;Starts with real junction
			     ;;If dump junction, dump it
			     (or (= (rejoining-junction-left-direction start-junction) dump-destination)
				 ;;Or if clearly before dump time, then we dump it.
				 (< (diamond-position-time start-diamond :a (rejoining-junction-a start-junction)
							   :b (rejoining-junction-b start-junction))
				    (- (local-time *dump-time*) fudge-global-coordinates)))
			   (output-loop-p start-diamond)))) ;VV-junction or :deleted.  See where dump enters diamond.
	end)
    (when (loop				;Scan string for end.  Return T if function should be called
	   (setq end (string-end-right d junction)) ;See if string ends with this diamond.
	   (when (eq end :invalid) (return nil)) ;Invalidated: give up now
	   (when end			;Found end of string
	     (return (or (and (rejoining-junction-p end) ;Ends with a dump-junction?  Must be valid.
			      (= (rejoining-junction-right-direction end) dump-destination))
			 must-output)))
	   (setf (diamond-processedp d) t) ;Stepping to right, so mark this one as done
	   (setq d (diamond-e d))	;Step to right
	   (check-diamond-step d :east)	;Check for data structure problems
	   (setq junction nil))		;Forget about starting junction for now
      ;;Here unless invalid
      ;;D is the last diamond in the string.  We did not mark it above, except if it is also the starting diamond
      (funcall function start-diamond start-junction d end)
      )))

;;Consider a string starting in this diamond at junction START, or entering from the left if START as NIL.
;;If the string ends at a junction in this diamond, return that junction.  If the entering string is invalidated by
;;a starting dump-junction before any ending junction, return :INVALID.
(defun string-end-right (diamond start)
  (let ((end (find-if #'(lambda (end)	;Scan possible ending junctions
			  (or (null start) ;Came in through right edge?  Then leftmost right junction is the end
			      (junction-left-p start end))) ;Otherwise leftmost that is to right of start
		      (diamond-right-junctions diamond) ;Places where string might end, in order from right
		      :from-end t))			;Scan through junctions from left
	(restart nil))
    (when *dump-time*			;If dumping, find first dump junction at this time if any
      (dolist (junction (received-left-dump-junctions diamond))
	(when (and (fudge= (rejoining-junction-dump-time junction) *dump-time* fudge-global-coordinates)
		   (or (null restart) ;Keep first of these junctions
		       (junction-left-p junction restart)))
	  (setq restart junction))))
    ;;END is the ending junction, if any.  RESTART is the first dump-starting junction, if any.
    (if (and restart			;If there is a restart, the entry is invalid if we can reach the restart
	     (junction-left-p restart end) ;End is later (or NIL, :DELETED, etc), so we reach restart
	     (junction-left-p start restart)) ;And start is earlier, so we to start after it.
	:invalid			;and so the original start is invalid
      end)				;If there is no restart, or the restart is to the right of the end, we are OK.
    ))

;;Called with a diamond that has (?) been examined before and has no starts, so it might be part if a loop.
;;Scan eastward, marking all diamonds.  If we find a loop, call function if loop has a diamond overlapping the dump
;;surface, or if not dumping.  Otherwise leave diamonds marked
;;so that they don't need scanning again.  If we find a marked diamond, then we know we previously examined it
;;for a loop and didn't find one, so we stop.
(defun map-string-paths-maybe-loop (start function)
  (loop with diamond = start
	until (diamond-processedp diamond) ;If marked already, not loop.
	until (string-end-right diamond nil) ;String ends in this diamond?  Not a loop.  Exit.  Do not mark.
	when (diamond-left-junctions diamond) ;We should never walk into a diamond with starts but no ends
	do (error "Scanned right into a diamond with starts and no ends")
	do (setf (diamond-processedp diamond) t) ;No junctions.  Mark to prevent repetitive scanning.
	do (setq diamond (diamond-e diamond)) ;Step east across diamond
	do (check-diamond-step diamond :east) ;Check for trouble
	when (eq diamond start)		;Looped back to first diamond (which is already marked)
	return (when (output-loop-p diamond)
		 (funcall function start))) ;Call function on loop.  
  )

;;Decide whether to output a loop or a string beginning with VV-junction or :deleted
(defun output-loop-p (diamond)
  (or (null *dump-time*) ;If not dumping, always output
      (null *size*)	 ;If infinite volume, always output
      ;;Doing dump.  We dump this diamond unless valid part of the diamond definitely starts in the future
      ;;of the dump
      (let ((local-dump-time (local-time *dump-time*)))
	(and
	 ;;Overall start of diamond must be before dump.
	 (< (4vector-t (diamond-start diamond)) (+ local-dump-time fudge-global-coordinates))
	 ;;If end is before dump also, then definitely dump
	 (or (< (4vector-t (diamond-end diamond)) (+ local-dump-time fudge-global-coordinates))
	     ;;Otherwise find where dump enters diamond and make sure it is after every past boundary
	     ;;Time can be slightly before diamond start because of fudge-global-coordinates
	     ;;If it is within fudge-coordinates, then find-left-edge-position will give a=b=0,
	     ;;but if not there can be trouble
	     (let ((entry
		    (if (< local-dump-time (4vector-t (diamond-start diamond))) ;If time is earlier than start
			(diamond-start diamond) ;use start
		      (multiple-value-bind (a b)
			  (find-left-edge-position diamond local-dump-time)
			(diamond-position diamond :a a :b b)))))
	       (loop for value across (xtoi entry)
		     always (> value (- fudge-ijkl)))))))))

(mirror-images
;;Return left junctions of diamond in a new list in order from left to right
;;Ignore the beginning rejoining junction if we're doing a dump and the first junction is in the future of it
(defun diamond-left-junctions (diamond)
  (let ((first (left-rejoining-junction diamond)) ;Starting junction for whole diamond, or :deleted
	(dumps
	 (and *dump-time*			;If dumping, get dump junctions for this time
	      (sort
	       (remove-if-not
		#'(lambda (junction)
		    (fudge= (rejoining-junction-dump-time junction) *dump-time* fudge-global-coordinates))
		(received-left-dump-junctions diamond))
	       #'junction-left-p))))	;Get in order from left to right
    (cond ((null first)			;No overall starting junction?
	   dumps)			;Just dumps, if any
	  ((null dumps)			;No dumps
	   (list first))		;Just first junction
	  (t (cons first dumps))	;First, then dump junctions
	  )))
)

;;Count number of diamonds in loop.  NIL if not loop.
(defun loop-count (diamond)
  (loop with start = diamond
	do (setq diamond (diamond-e diamond))
	;unless diamond return nil
	unless (diamondp diamond) return nil  ;;modified by SJW
	count t
	until (eq diamond start)))

;;Return number of diamonds in each string in 3 lists: loops, regular open strings, and partly deleted strings
(defun string-counts ()
  (let ((loop-lengths nil)
	(open-lengths nil)
	(deleted-lengths nil))
    (map-string-paths
     #'(lambda (start-diamond start-junction end-diamond end-junction)
	 (loop for diamond = start-diamond then (diamond-e diamond)
	       count t into count
	       until (eq diamond end-diamond)
	       finally 	 (if (or (eq start-junction :deleted)
				 (eq end-junction :deleted))
			     (push count deleted-lengths)
			   (push count open-lengths))))
     #'(lambda (diamond)
	 (push (loop-count diamond) loop-lengths)))
    (values loop-lengths open-lengths deleted-lengths)))

;;Get strings without worrying about time
(defun get-all-strings ()
  (let ((data nil))
    (map-string-paths
     #'(lambda (diamond &rest rest)
	 (push (apply #'get-string-diamonds diamond rest) 
	       data)))
    data))

;;Get list of paths, not including those that are partially deleted
(defun get-paths (&optional (time (current-time)))
  (let ((data nil))
    (map-string-paths
     #'(lambda (diamond &rest rest)
	 (let ((path (apply #'get-path time diamond rest)))
	   (when path			;NIL means :deleted
	     (push path data)))))
    data))

;;Get list of closed paths
(defun get-closed-paths (&optional (time (current-time)))
  (let ((data nil))
    (map-string-paths
     #'(lambda (&rest ignore) (declare (ignore ignore)))
     #'(lambda (diamond)
	 (push (get-path time diamond) data)))
    data))

;; Get list of loops longer than current conformal time 
(defun get-superhorizon-closed-paths (&optional (time (current-time)) (min-length 1.0))
  (let ((closed-paths (get-closed-paths time)))
    (loop for path in closed-paths
	  if (> (path-length path) (* min-length time))
	  collect path)))
    

;;Get list of open paths
(defun get-open-paths (&optional (time (current-time)))
  (let ((data nil))
    (map-string-paths
     #'(lambda (start-diamond start-junction end-diamond end-junction)
	 (push (get-path time start-diamond start-junction end-diamond end-junction) data))
     #'(lambda (&rest ignore) (declare (ignore ignore))))
    data))

;;Return length of segment from get-path.  See old notes.text
(defun path-segment-length (diamond a1 b1 a2 b2)
  (declare (ignore b1 b2))
  (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond))))

;;Give the size of a path from get-path.
(defun path-length (path)
  (loop for entry in path
	sum (apply #'path-segment-length entry)))

;;Return list of path lengths
(defun path-lengths (&optional (paths (get-paths (current-time))))
  (loop for string in paths
	for length = (path-length string)
	collect length))


(defun path-segment-bh-length (diamond a1 b1 a2 b2)
  (declare (ignore b1 b2))
  (cond ((and (eq (diamond-w diamond) :BH) (eq (diamond-e diamond) :BH)) 
	 (list (diamond-bh diamond) (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	((and (eq (diamond-w diamond) :BHeatit) (eq (diamond-e diamond) :BHeatit))
         (list (diamond-bh diamond) (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	((and (eq (diamond-w diamond) :BH) (eq (diamond-e diamond) :BHeatit))
         (list (diamond-bh diamond) (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	((and (eq (diamond-w diamond) :BHeatit) (eq (diamond-e diamond) :BH))
         (list (diamond-bh diamond) (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	 ;(list (diamond-w diamond) (diamond-e diamond) 1.0))   
	;((and (eq (diamond-w diamond) :BHdeleted) (eq (diamond-e diamond) :BHdeleted))
	 ;(list (diamond-w diamond) (diamond-e diamond) (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
	 ;(list (diamond-w diamond) (diamond-e diamond) 1.0))  
	;((and (eq (diamond-w diamond) :BH) (eq (diamond-e diamond) :BHdeleted))
         ;(list (diamond-bh diamond) nil (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
	 ;(list (diamond-w diamond) (diamond-e diamond) 1.0))
	;((and (eq (diamond-w diamond) :BHdeleted) (eq (diamond-e diamond) :BH))
         ;(list nil (diamond-bh diamond) (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
	 ;(list (diamond-w diamond) (diamond-e diamond) 1.0))
	;((and (eq (diamond-w diamond) :BHpropdel) (eq (diamond-e diamond) :BHpropdel))
	 ;(list (diamond-w diamond) (diamond-e diamond) (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
	;((and (eq (diamond-w diamond) :BH) (eq (diamond-e diamond) :BHpropdel))
         ;(list (diamond-bh diamond) nil (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
        ;((and (eq (diamond-w diamond) :BHpropdel) (eq (diamond-e diamond) :BH))
         ;(list nil (diamond-bh diamond) (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
        ;((and (eq (diamond-w diamond) :BHpropdel) (eq (diamond-e diamond) :BHdeleted))
         ;(list (diamond-w diamond) (diamond-e diamond) (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
        ;((and (eq (diamond-w diamond) :BHdeleted) (eq (diamond-e diamond) :BHpropdel))
         ;(list (diamond-w diamond) (diamond-e diamond) (* 2 (abs (- a1 a2)) (4vector-t (diamond-p diamond)))))
	((eq (diamond-w diamond) :BH)
         (list (diamond-bh diamond) nil (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	 ;(list (diamond-w diamond) nil 1.0))
	((eq (diamond-e diamond) :BH)
	 (list nil (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	((eq (diamond-w diamond) :BHeatit)
         (list (diamond-bh diamond) nil (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	((eq (diamond-e diamond) :BHeatit)
         (list nil (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	 ;(list nil (diamond-e diamond) 1.0))
	((eq (diamond-w diamond) :BHdeleted)
         (list (diamond-bh diamond) nil (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	 ;(list (diamond-w diamond) nil 1.0))
	((eq (diamond-e diamond) :BHdeleted)
	 (list nil (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	 ;(list nil (diamond-e diamond) 1.0))
	((eq (diamond-w diamond) :BHpropdel)
         (list (diamond-bh diamond) nil (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	((eq (diamond-e diamond) :BHpropdel)
         (list nil (diamond-bh diamond) (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))
	(t
	 (list nil nil (* 2 (- a1 a2) (4vector-t (diamond-p diamond)))))))
         ;(list nil nil 1.0))))


(defun path-bh-length (path)
  (let* ((length 0.0)
	 bh-w
	 bh-e
	 reads)
    (loop for entry in path
	  do (setq reads (apply #'path-segment-bh-length entry))
	  do (when (car reads)
	       (setq bh-w (car reads)))
	  do (when (car (cdr reads))
	       (setq bh-e (car (cdr reads))))
	  do (when (> (car (cdr (cdr reads))) 5.0)
	       (warn "~%Possible cranky diamond: L: ~S~%" (car (cdr (cdr reads)))))
	  do (setq length (+ length (car (cdr (cdr reads))))))
    (list bh-w bh-e length)))
			     
(defun path-bh-lengths-free (&optional (paths (get-paths (current-time))))
    (loop for string in paths
	  for bh-length = (path-bh-length string)
	  unless (and (car bh-length) (car (cdr bh-length)))
	  sum (car (cdr (cdr bh-length))) into length-free
	  finally (return (list length-free))))

(defun path-bh-lengths-bh (&optional (paths (get-paths (current-time))))
    (loop for string in paths
          for bh-length = (path-bh-length string)
          when (and (car bh-length) (car (cdr bh-length)))
          sum (car (cdr (cdr bh-length))) into length-bh
          finally (return (list length-bh))))


;;Return if path has a BH in one of the ends
(defun path-bh (path)
  (let ((bh-str 0.0))
    (loop for entry in path
	  do (if (apply #'path-segment-bh entry)
		 (setf bh-str 1.0)))
    bh-str))

;; Checks if any neighbour of the diamond is a BH
(defun path-segment-bh (diamond a1 b1 a2 b2)
  (declare (ignore a1 b1 a2 b2))
  (if (or (eq (diamond-e diamond) :BH)
	  (eq (diamond-w diamond) :BH))
      (return-from path-segment-bh t)))

;;Return list of loop (or open string) sizes, and the number of kinks on them, sorted by kink number
(defun loop-sizes-kinks (&optional (paths (get-paths (current-time))))
  (let ((unsorted
	 (loop for path in paths
	       for length = (path-length path)
	       collect (list length (length path)))))
    (sort unsorted #'(lambda (x y)
		       (< (nth 1 x) (nth 1 y))))))

(defun total-length (&optional (paths (get-paths (current-time))))
  (reduce #'+ (path-lengths paths)))

;;Fraction of total length in sub-horizon loops
(defun loop-length-fraction (&optional (cutoff (current-time)) (paths (get-paths (current-time))))
  (loop for length in (path-lengths paths)
	sum length into total
	when (< length cutoff)
	sum length into short
	finally (return (/ short total))))

(defun print-loop-length-fractions (directory start end step)
  (loop with *allow-missing-predecessor-files* = t
	for time from start to end by step
	do (read-dumps directory :time time)
	do (format t "~A: ~4$~%" (format-float-reasonably time) (loop-length-fraction))))

(defun bh-length (&optional (paths (get-paths (current-time))))
  (let (ttl-free
	ttl-bh)
    (loop for length-free in (path-bh-lengths-free paths)
	  sum length-free into total-free
	  finally (setq ttl-free total-free))
    (loop for length-bh in (path-bh-lengths-bh paths)
	  sum length-bh into total-bh
	  finally (setq ttl-bh total-bh))
    (list ttl-free ttl-bh)))

(defun bh-length-new ()
  (let ((ttl-free 0.0d0)
	(ttl-bh 0.0d0)
	(ttl-all 0.0d0))
    (loop for string in *dumped-segments*
	  do (when (and (or (eq (dump-segment-start-junction string) :BH)
			    (eq (dump-segment-start-junction string) :BHeatit))
			(or (eq (dump-segment-end-junction string) :BH)
			    (eq (dump-segment-end-junction string) :BHeatit)))
	       (setf ttl-bh (+ ttl-bh (dump-segment-length string))))
	  do (when (and (eq (dump-segment-start-junction string) :loop)
			(eq (dump-segment-end-junction string) :loop))
	       (setf ttl-free (+ ttl-free (dump-segment-length string)))))
    (setf ttl-all (+ ttl-bh ttl-free))
    (list ttl-free ttl-bh ttl-all)))

(defun print-bh-length (directory start end step)
  (let* ((stream (open (concatenate 'string directory "/bhlengths.txt") :direction :output :if-does-not-exist :create :if-exists :supersede)))
	 (loop with *allow-missing-predecessor-files* = t
	       for time from start to end by step
	       do (setf *read-dump-diamonds* nil)
	       do (read-dumps directory :time time)
	       do (format stream "~A ~{~a~^ ~} ~%" (format-float-reasonably time) (bh-length-new)))
	 (close stream)))
	 
(defun bh-stats (dir bhn)
  (with-open-file (stream (concatenate 'string dir "/blackholes.dat") :if-does-not-exist nil)
		  (let* ((center (make-3vector 0.0d0 0.0d0 0.0d0))
                         (size 0.d0)
                         (bh (make-blackhole :center center :size size))
                         (bhs (make-hash-table :test 'equalp))
                         (nmb 0)
                         (zero 0)
                         (one 0)
                         (two 0)
                         (three 0)
                         (four 0)
                         (five 0)
                         (six 0)
                         (seven 0)
                         (eight 0)
                         (nine 0)
                         (ten 0)
			 (eleven 0)
                         (twelve 0)
                         (thirteen 0)
                         (fourteen 0)
			 (fifteen 0)
			 (sixteen 0)
			 (seventeen 0)
			 (eightteen 0)
			 (nineteen 0)
			 (twenty 0)
			 (twentyone 0)
			 (twentytwo 0)
			 (twentythree 0)
			 (twentyfour 0)
			 (twentyfive 0)
			 (twentysix 0)
			 (twentyseven 0)
			 (twentyeight 0)
			 (twentynine 0)
			 (therty 0)
                         bhkey)
		    (loop for i below bhn
                          do (setf bh (read stream nil))
                          do (setf center (blackhole-center bh))
                          do (setf (gethash center bhs) 0)
                          do (push center bhkey))
		    (loop for string in *dumped-segments*
			  do (when (or (eq (dump-segment-start-junction string) :BHdeleted)
				       (eq (dump-segment-end-junction string) :BHdeleted)
				       (eq (dump-segment-start-junction string) :BHpropdel)
                                       (eq (dump-segment-end-junction string) :BHdpropdel))
				       ;(eq (dump-segment-start-junction string) :deleted)
                                       ;(eq (dump-segment-end-junction string) :deleted))
			       (format t "Seg: ~S~%" string))
			  do (when (and (or (eq (dump-segment-start-junction string) :BH)
					    (eq (dump-segment-start-junction string) :BHeatit))
					(or (eq (dump-segment-end-junction string) :BH)
					    (eq (dump-segment-end-junction string) :BHeatit)))
			       (when (null (dump-segment-start-bh string))
				 (format t "Str: ~S~%" string))
			       (setf center (blackhole-center (dump-segment-start-bh string)))
                               (setf nmb (+ (gethash center bhs) 1))
                               (setf (gethash center bhs) nmb))
                          do (when (and (or (eq (dump-segment-start-junction string) :BH)
                                            (eq (dump-segment-start-junction string) :BHeatit))
                                        (or (eq (dump-segment-end-junction string) :BH)
                                            (eq (dump-segment-end-junction string) :BHeatit)))
                               (setf center (blackhole-center (dump-segment-end-bh string)))
                               (setf nmb (+ (gethash center bhs) 1))
                               (setf (gethash center bhs) nmb))
                          )
                  (loop for center in bhkey
                        do (setf nmb (gethash center bhs))
                        do (cond ((eq nmb 0)
                                  (setf zero (+ zero 1)))
                                 ((eq nmb 1)
                                  (setf one (+ one 1)))
                                 ((eq nmb 2)
                                  (setf two (+ two 1)))
                                 ((eq nmb 3)
                                  (setf three (+ three 1)))
                                 ((eq nmb 4)
                                  (setf four (+ four 1)))
                                 ((eq nmb 5)
                                  (setf five (+ five 1)))
                                 ((eq nmb 6)
                                  (setf six (+ six 1)))
                                 ((eq nmb 7)
                                  (setf seven (+ seven 1)))
                                 ((eq nmb 8)
                                  (setf eight (+ eight 1)))
                                 ((eq nmb 9)
                                  (setf nine (+ nine 1)))
                                 ((eq nmb 10)
                                  (setf ten (+ ten 1)))
                                 ((eq nmb 11)
				  (setf eleven (+ eleven 1)))
                                 ((eq nmb 12)
                                  (setf twelve (+ twelve 1)))
                                 ((eq nmb 13)
                                  (setf thirteen (+ thirteen 1)))
                                 ((eq nmb 14)
                                  (setf fourteen (+ fourteen 1)))
				 ((eq nmb 15)
                                  (setf fifteen (+ fifteen 1)))
				 ((eq nmb 16)
                                  (setf sixteen (+ sixteen 1)))
				 ((eq nmb 17)
                                  (setf seventeen (+ seventeen 1)))
				 ((eq nmb 18)
                                  (setf eightteen (+ eightteen 1)))
				 ((eq nmb 19)
                                  (setf nineteen (+ nineteen 1)))
				 ((eq nmb 20)
                                  (setf twenty (+ twenty 1)))
				 ((eq nmb 21)
                                  (setf twentyone (+ twentyone 1)))
				 ((eq nmb 22)
                                  (setf twentytwo (+ twentytwo 1)))
				 ((eq nmb 23)
                                  (setf twentythree (+ twentythree 1)))
				 ((eq nmb 24)
                                  (setf twentyfour (+ twentyfour 1)))
				 ((eq nmb 25)
                                  (setf twentyfive (+ twentyfive 1)))
				 ((eq nmb 26)
                                  (setf twentysix (+ twentysix 1)))
				 ((eq nmb 27)
                                  (setf twentyseven (+ twentyseven 1)))
				 ((eq nmb 28)
                                  (setf twentyeight (+ twentyeight 1)))
				 ((eq nmb 29)
                                  (setf twentynine (+ twentynine 1)))
				 ((eq nmb 30)
                                  (setf therty (+ therty 1)))
                                 ))
                  (list zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eightteen nineteen twenty twentyone twentytwo twentythree twentyfour twentyfive twentysix twentyseven twentyeight twentynine therty)))) 


(defun bh-stats-o (directory bhn &optional (paths (get-paths (current-time))))
  (with-open-file (stream (concatenate 'string directory "/blackholes.dat") :if-does-not-exist nil)
                  (let* ((center (make-3vector 0.0d0 0.0d0 0.0d0))
                         (size 0.d0)
                         (bh (make-blackhole :center center :size size))
			 (bhs (make-hash-table :test 'equalp))
			 (nmb 0)
			 (zero 0)
			 (one 0)
			 (two 0)
			 (three 0)
			 (four 0)
			 (five 0)
			 (six 0)
			 (seven 0)
			 (eight 0)
			 (nine 0)
			 (ten 0)
			 (eleven 0)
			 (twelve 0)
			 (thirteen 0)
			 (fourteen 0)
			 bhkey
			 bhlist
			 (count 1))
                    (loop for i below bhn
                          do (setf bh (read stream nil))
			  do (setf center (blackhole-center bh))
			  do (setf (gethash center bhs) 0)
			  do (push center bhkey))
		    (loop for string in paths
			  for bh-length = (path-bh-length string)
			  do (when (and (car bh-length) (car (cdr bh-length)))
			       (push bh-length bhlist)
			       (setf count (+ count 1)))
			  do (when (and (car bh-length) (car (cdr bh-length)))
			       (setf center (blackhole-center (car bh-length)))
			       (setf nmb (+ (gethash center bhs) 1))
			       (setf (gethash center bhs) nmb))
			  do (when (and (car bh-length) (car (cdr bh-length)))
			       (setf center (blackhole-center (car (cdr bh-length))))
			       (setf nmb (+ (gethash center bhs) 1))
			       (setf (gethash center bhs) nmb))
			  )
		  (loop for center in bhkey
			do (setf nmb (gethash center bhs))
			do (cond ((eq nmb 0)
				  (setf zero (+ zero 1)))
				 ((eq nmb 1)
				  (setf one (+ one 1)))
				 ((eq nmb 2)
				  (setf two (+ two 1)))
				 ((eq nmb 3)
				  (setf three (+ three 1)))
				 ((eq nmb 4)
				  (setf four (+ four 1)))
				 ((eq nmb 5)
				  (setf five (+ five 1)))
				 ((eq nmb 6)
				  (setf six (+ six 1)))
				 ((eq nmb 7)
                                  (setf seven (+ seven 1)))
				 ((eq nmb 8)
                                  (setf eight (+ eight 1)))
				 ((eq nmb 9)
                                  (setf nine (+ nine 1)))
				 ((eq nmb 10)
                                  (setf ten (+ ten 1)))
				 ((eq nmb 11)
                                  (setf eleven (+ eleven 1)))
				 ((eq nmb 12)
                                  (setf twelve (+ twelve 1)))
				 ((eq nmb 13)
                                  (setf thirteen (+ thirteen 1)))
				 ((eq nmb 14)
                                  (setf fourteen (+ fourteen 1)))
				 ))
		  bhlist)))
		  ;(list zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen))))

(defun bh-single-old (directory bha bhn str-out time)
  (with-open-file (stream (concatenate 'string directory "/blackholes.dat") :if-does-not-exist nil)
		  (let* ((center (make-3vector 0.0d0 0.0d0 0.0d0))
			 (centera (make-3vector 0.0d0 0.0d0 0.0d0))
			 (size 0.d0)
			 (bh (make-blackhole :center center :size size))
			 (bhs (make-hash-table :test 'equalp))
			 (segments nil)
			 (ns 1)
			 (bhnn 1)
			 (struc nil))
		    (loop for i below bhn
                          do (setf bh (read stream nil))
                          do (setf center (blackhole-center bh))
                          do (setf (gethash center bhs) bhnn)
                          do (when (eq bhnn bha)
			       (setf centera (blackhole-center bh)))
			  do (setf bhnn (+ bhnn 1)))
		    (loop for str in *dumped-segments*
			  do (when (dump-segment-start-bh str)
			       (when (null (dump-segment-end-bh str))
				 (error "This segment: ~S~%" str)))
			  do (when (dump-segment-end-bh str)
			       (when (null (dump-segment-start-bh str))
				 (error "This segment: ~S~%" str)))
			  do (when (or (and (dump-segment-start-bh str) (zerop (3vector-distance centera (blackhole-center (dump-segment-start-bh str)))))
				       (and (dump-segment-end-bh str) (zerop (3vector-distance centera (blackhole-center (dump-segment-end-bh str))))))
			       (push str segments)))
		    (format str-out "~%#### Analysing BH ~S at time ~S ####~%" (gethash centera bhs) time)
		    (loop for sg in segments
			  do (if (zerop (3vector-distance (blackhole-center (dump-segment-start-bh sg)) (blackhole-center (dump-segment-end-bh sg))))
				 (progn
				   (push (gethash centera bhs) struc)
				   (format str-out "Segment ~S is a Loop with length ~S~%" ns (dump-segment-length sg)))
			       (if (zerop (3vector-distance centera (blackhole-center (dump-segment-start-bh sg))))
				   (progn
				    (push (gethash (blackhole-center (dump-segment-end-bh sg)) bhs) struc)
				    (format str-out "Segment ~S links BH ~S with length ~S (~S)~%" ns (gethash (blackhole-center (dump-segment-end-bh sg)) bhs) (dump-segment-length sg) (3vector-distance-wrap centera (blackhole-center (dump-segment-end-bh sg)))))
				 (progn
				   (push (gethash (blackhole-center (dump-segment-start-bh sg)) bhs) struc)
				   (format str-out "Segment ~S links BH ~S with length ~S (~S)~%" ns (gethash (blackhole-center (dump-segment-start-bh sg)) bhs) (dump-segment-length sg) (3vector-distance-wrap centera (blackhole-center (dump-segment-start-bh sg)))))))
			  do (setf ns (+ ns 1)))
		    (loop for i in struc
			  do (format str-out "~S<->~S," (gethash centera bhs) i)))))


(defun bh-single (directory bha bhn str-out time)
  (with-open-file (stream (concatenate 'string directory "/blackholes.dat") :if-does-not-exist nil)
                  (let* ((center1 (make-3vector 0.0d0 0.0d0 0.0d0))
                         (center2 (make-3vector 0.0d0 0.0d0 0.0d0))
                         (size 0.d0)
                         (bh (make-blackhole :center center1 :size size))
                         (bhs (make-hash-table :test 'equalp))
                         (bh-labels (make-hash-table :test 'equalp))
                         (bhnn 1)
                         list1
                         list2
                         bhkey
                         structures)
		    (format str-out "~% ####### Structure Information at t=~S ####~%" (current-time))
                    (loop for i below bhn
                          do (setf bh (read stream nil))
                          do (setf center1 (blackhole-center bh))
                          do (setf (gethash center1 bhs) '())
                          do (setf (gethash center1 bh-labels) bhnn)
                           do (setf bhnn (+ bhnn 1))
                           do (push center1 bhkey))
                    (loop for string in *dumped-segments*
			  do (when (and (or (eq (dump-segment-start-junction string) :BH)
                                            (eq (dump-segment-start-junction string) :BHeatit))
                                        (or (eq (dump-segment-end-junction string) :BH)
                                            (eq (dump-segment-end-junction string) :BHeatit)))
                               (setf center1 (blackhole-center (dump-segment-start-bh string)))
                               (setf center2 (blackhole-center (dump-segment-end-bh string)))
                               (setf list1 (push (list (dump-segment-start-bh string) (dump-segment-end-bh string) (dump-segment-length string)) (gethash center1 bhs)))
                               (setf list2 (push (list (dump-segment-start-bh string) (dump-segment-end-bh string) (dump-segment-length string)) (gethash center2 bhs)))
                               (setf (gethash center1 bhs) list1)
                               (setf (gethash center2 bhs) list2)))
                    (loop for center1 in bhkey
                          do (when (gethash center1 bhs)
                               (let ((segs (remove-duplicates (gethash center1 bhs) :test #'eqseg))
                                     (nn 0))
                                 (loop for sg in segs
                                       do (if (zerop (3vector-distance (blackhole-center (car sg)) (blackhole-center (car (cdr sg)))))
                                              (setf nn (+ nn 2))
                                            (setf nn (+ nn 1))))
                                 (when (oddp nn)
                                   (format t "Center: ~S ~S~%" center1 (gethash center1 bh-labels))
                                   (format t "Seg: ~S~%" segs)))))
                    (loop for k below 2
                          do (loop for center1 in bhkey
                                   do (when (gethash center1 bhs)
					(loop for center2 in bhkey
                                              do (when (and (null (zerop (3vector-distance center1 center2)))
                                                            (gethash center2 bhs))
                                                   (loop for str1 in (gethash center1 bhs)
                                                         do (loop for str2 in (gethash center2 bhs)
                                                                  do (when (or (zerop (3vector-distance (blackhole-center (car str1)) (blackhole-center (car str2))))
                                                                               (zerop (3vector-distance (blackhole-center (car str1)) (blackhole-center (car (cdr str2)))))
                                                                               (zerop (3vector-distance (blackhole-center (car (cdr str1))) (blackhole-center (car str2))))
                                                                               (zerop (3vector-distance (blackhole-center (car (cdr str1))) (blackhole-center (car (cdr str2))))))
                                                                       (setf list1 (append (gethash center1 bhs) (gethash center2 bhs)))
                                                                       (setf list2 nil)
                                                                       (setf (gethash center1 bhs) list1)
                                                                       (setf (gethash center2 bhs) list2)))))))))
                    (loop for center1 in bhkey
                          do (when (gethash center1 bhs)
                               (push (remove-duplicates (gethash center1 bhs) :test #'eqseg) structures)))
                    (loop for struc in structures
                          do (let* ((strubhs (make-hash-table :test 'equalp))
                                    (nloop 0)
                                    (lloop 0.0d0)
                                    (bh-list nil)
                                    (nbh 1)
                                    (nstr 0)
                                    (lstr 0.0d0)
                                    (loops nil)
				    (print-struc nil))
                               (loop for str in struc
                                     do (if (zerop (3vector-distance (blackhole-center (car str)) (blackhole-center (car (cdr str)))))
                                            (progn
                                              (push (list time (gethash (blackhole-center (car str)) bh-labels) (car (cdr (cdr str)))) loops)
					      (setf nloop (+ nloop 1))
                                              (setf lloop (+ lloop (car (cdr (cdr str))))))
                                          (progn
                                            (setf nstr (+ nstr 1))
                                            (setf lstr (+ lstr (car (cdr (cdr str)))))
					    (push (3vector-length (blackhole-center (car str))) bh-list)
                                            (push (3vector-length (blackhole-center (car (cdr str)))) bh-list)))
                                     do (when (null (gethash (blackhole-center (car str)) strubhs))
                                          (setf (gethash (blackhole-center (car str)) strubhs) nbh)
                                          (setf nbh (+ nbh 1)))
                                     do (when (null (gethash (blackhole-center (car (cdr str))) strubhs))
                                          (setf (gethash (blackhole-center (car (cdr str))) strubhs) nbh)
                                          (setf nbh (+ nbh 1))))
                               (setf bh-list  (remove-duplicates bh-list))
			       (loop for str in struc
				     do (when (or (eq (gethash (blackhole-center (car str)) bh-labels) bha)
						  (eq (gethash (blackhole-center (car (cdr str))) bh-labels) bha))
					  (setf print-struc t))
				     until print-struc)
			       (when print-struc
				 (format str-out "~%# STRUCTURE #~%")
				 (format str-out "~S  loops with total length: ~S~%" nloop lloop)
				 (format str-out "Estructure with ~S bhs and ~S segments and total length ~S ~%" (length bh-list) nstr lstr)
				 (loop for str in struc
				       do (format str-out "~S<->~S, " (gethash (blackhole-center (car str)) bh-labels) (gethash (blackhole-center (car (cdr str))) bh-labels)))))))))
                    

(defun bh-struct (directory bhn str-out str-out2 str-loopbh str-seg str-loop inf-seg time tbh) ;&optional (paths (get-paths (current-time))))
  (with-open-file (stream (concatenate 'string directory "/blackholes.dat") :if-does-not-exist nil)
                  (let* ((center1 (make-3vector 0.0d0 0.0d0 0.0d0))
			 (center2 (make-3vector 0.0d0 0.0d0 0.0d0))
                         (size 0.d0)
                         (bh (make-blackhole :center center1 :size size))
                         (bhs (make-hash-table :test 'equalp))
			 (bh-labels (make-hash-table :test 'equalp))
			 (bhnn 1)
			 list1
			 list2
			 bhkey
			 structures
			 (Nfti 0)
			 (Nft1 0)
			 (Nft2 0)
			 (Nft3 0)
			 (Nft4 0)
			 (done nil))
		    (format str-out "~% ####### Structure Information at t=~S ####~%" (current-time)) 
                    (loop for i below bhn
                          do (setf bh (read stream nil))
                          do (setf center1 (blackhole-center bh))
                          do (setf (gethash center1 bhs) '())
			  do (setf (gethash center1 bh-labels) bhnn)
			  do (setf bhnn (+ bhnn 1))
			  do (push center1 bhkey))
		    (loop for string in *dumped-segments*;paths
                          ;for bh-length = (path-bh-length string)
			  do (when (and (eq (dump-segment-start-junction string) :loop)
					(eq (dump-segment-end-junction string) :loop))
			       (format str-loop "~S ~S~%" (dump-segment-dump-time string) (dump-segment-length string)))
			  do (when (and (or (eq (dump-segment-start-junction string) :BH)
                                            (eq (dump-segment-start-junction string) :BHeatit))
                                        (or (eq (dump-segment-end-junction string) :BH)
                                            (eq (dump-segment-end-junction string) :BHeatit)))
			       (setf center1 (blackhole-center (dump-segment-start-bh string)))
			       (setf center2 (blackhole-center (dump-segment-end-bh string)))
                               (setf list1 (push (list (dump-segment-start-bh string) (dump-segment-end-bh string) (dump-segment-length string)) (gethash center1 bhs)))
			       (setf list2 (push (list (dump-segment-start-bh string) (dump-segment-end-bh string) (dump-segment-length string)) (gethash center2 bhs)))
                               (setf (gethash center1 bhs) list1)
			       (setf (gethash center2 bhs) list2)))
		    (loop for center1 in bhkey
			  do (when (gethash center1 bhs)
			       (let ((segs (remove-duplicates (gethash center1 bhs) :test #'eqseg))
				     (nn 0))
				 (loop for sg in segs
				       do (if (zerop (3vector-distance (blackhole-center (car sg)) (blackhole-center (car (cdr sg)))))
					      (setf nn (+ nn 2))
					    (setf nn (+ nn 1))))
				 (when (oddp nn)
				   (format t "Center: ~S ~S~%" center1 (gethash center1 bh-labels))
				   (format t "Seg: ~S~%" segs)))))
		    (loop
		     do (setf done t)
		     do (loop for center1 in bhkey
			      do (when (gethash center1 bhs)
				   (loop for center2 in bhkey
					 do (when (and (null (zerop (3vector-distance center1 center2)))
						       (gethash center2 bhs))
					      (loop for str1 in (gethash center1 bhs)
						    do (loop for str2 in (gethash center2 bhs)
							     do (when (or (zerop (3vector-distance (blackhole-center (car str1)) (blackhole-center (car str2)))) 
									  (zerop (3vector-distance (blackhole-center (car str1)) (blackhole-center (car (cdr str2)))))
									  (zerop (3vector-distance (blackhole-center (car (cdr str1))) (blackhole-center (car str2))))
									  (zerop (3vector-distance (blackhole-center (car (cdr str1))) (blackhole-center (car (cdr str2))))))
								  (setf done nil)
								  (setf list1 (append (gethash center1 bhs) (gethash center2 bhs)))
								  (setf list2 nil)
								  (setf (gethash center1 bhs) list1)
								  (setf (gethash center2 bhs) list2))))))))
		     until done)
		    (loop for center1 in bhkey
			  do (when (gethash center1 bhs)
			       (push (remove-duplicates (gethash center1 bhs) :test #'eqseg) structures)))
		    (loop for struc in structures
			  do (let* ((strubhs (make-hash-table :test 'equalp))
				    (nloop 0)
				    (lloop 0.0d0)
				    (bh-list nil)
				    (nbh 1)
				    (nstr 0)
				    (lstr 0.0d0)
				    (loops nil))
			       (loop for str in struc
				     do (if (zerop (3vector-distance (blackhole-center (car str)) (blackhole-center (car (cdr str)))))
					    (progn
					      (push (list time (gethash (blackhole-center (car str)) bh-labels) (car (cdr (cdr str)))) loops)
					      ;(format str-loop "~S ~S ~S~%" time (gethash (blackhole-center (car str)) bh-labels) (car (cdr (cdr str))))
					      (setf nloop (+ nloop 1))
					      (setf lloop (+ lloop (car (cdr (cdr str))))))
					  (progn
					    (setf nstr (+ nstr 1))
					    (setf lstr (+ lstr (car (cdr (cdr str)))))
					    ;(format str-seg "~S ~S ~S ~S ~S~%" time (gethash (blackhole-center (car str)) bh-labels) (gethash (blackhole-center (car (cdr str))) bh-labels)
						    ;(car (cdr (cdr str))) (- (3vector-distance-wrap (blackhole-center (car str)) (blackhole-center (car (cdr str)))) (* 2.0 (/ (* (blackhole-size (car str)) tbh) time)))) 
					    (push (gethash (blackhole-center (car str)) bh-labels) bh-list)
					    (push (gethash (blackhole-center (car (cdr str))) bh-labels) bh-list)))
				     do (when (null (gethash (blackhole-center (car str)) strubhs))
					  (setf (gethash (blackhole-center (car str)) strubhs) nbh)
					  (setf nbh (+ nbh 1)))
				     do (when (null (gethash (blackhole-center (car (cdr str))) strubhs))
					  (setf (gethash (blackhole-center (car (cdr str))) strubhs) nbh)
					  (setf nbh (+ nbh 1))))
			       (setf bh-list  (remove-duplicates bh-list))
			       (loop for str in struc
				     do (format str-seg "~S ~S ~S ~S ~S ~S~%" time (length bh-list)  (gethash (blackhole-center (car str)) bh-labels) (gethash (blackhole-center (car (cdr str))) bh-labels)
                                                    (car (cdr (cdr str))) (- (3vector-distance-wrap (blackhole-center (car str)) (blackhole-center (car (cdr str)))) (* 2.0 (/ (* (blackhole-size (car str)) tbh) time)))))
			       (when (zerop nstr)
				 (loop for lp in loops
				       do (format str-loopbh "~S ~S ~S~%" (car lp) (car (cdr lp)) (car (cdr (cdr lp))))))
			       (when (> lstr time)
				 (loop for str in struc
				       do (format inf-seg "~S ~S ~S ~S ~S~%" time (gethash (blackhole-center (car str)) bh-labels) (gethash (blackhole-center (car (cdr str))) bh-labels)
					       (car (cdr (cdr str))) (- (3vector-distance-wrap (blackhole-center (car str)) (blackhole-center (car (cdr str)))) (* 2.0 (/ (* (blackhole-size (car str)) tbh) time))))))
			       (if (zerop nstr)
				   (if (> lloop time)
				       (setf Nfti (+ Nfti 1))
				     (setf Nft1 (+ Nft1 1)))
				 (if (> lstr time)
				     (setf Nfti (+ (length bh-list) Nfti))
				   (cond ((eq (length bh-list) 2)
					  (setf Nft2 (+ (length bh-list) Nft2)))
					 ((eq (length bh-list) 3)
					  (setf Nft3 (+ (length bh-list) Nft3)))
					 ((eq (length bh-list) 4)
					  (setf Nft4 (+ (length bh-list) Nft4))))))
			       (format str-out "~%# STRUCTURE #~%")
			       (format str-out "~S  loops with total length: ~S~%" nloop lloop)
			       (format str-out "Estructure with ~S bhs and ~S segments and total length ~S ~%" (length bh-list) nstr lstr)
			       (loop for str in struc
				     do (when (or (null (gethash (blackhole-center (car str)) bh-labels))
						  (null (gethash (blackhole-center (car (cdr str))) bh-labels)))
					  (format t "Stru: ~S~%" str))
                                     do (format str-out "~S<->~S, " (gethash (blackhole-center (car str)) bh-labels) (gethash (blackhole-center (car (cdr str))) bh-labels)))))
		    (format str-out2 "~S ~S ~S ~S ~S ~S~%" time Nfti Nft1 Nft2 Nft3 Nft4)))) 

			       
(defun eqseg (seg1 seg2)
  (if (eq (car (cdr (cdr seg1))) (car (cdr (cdr seg2))))
      (return-from eqseg t)
     (return-from eqseg nil)))
  
(defun bh-stats-c (&optional (paths (get-paths (current-time))))
  (loop for string in paths
	for bh-length = (path-bh-length string)
	do (when (or (car bh-length) (car (cdr bh-length)))                                                                                                                               
	     (format t "STR: ~S~%" bh-length))))  
    

			 
(defun bh-stats-check (directory start end step)
  (loop with *allow-missing-predecessor-files* = t
	for time from start to end by step
	do (read-dumps directory :time time)
	do (bh-stats-c)))
			  
(defun print-bh-stats (directory start end step bhn)
  (let* ((stream (open (concatenate 'string directory "/bhstats.txt") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (bhlist nil))
    (loop with *allow-missing-predecessor-files* = t
	  for time from start to end by step
	  do (setf *read-dump-diamonds* nil)
	  do (read-dumps directory :time time)
          ;do (setf bhlist (bh-stats directory bhn)))
	  do (format stream "~A ~{~a~^ ~} ~%" (format-float-reasonably time) (bh-stats directory bhn)))
    (close stream)
    bhlist))

(defun print-bh-structures (directory start end step bhn tbh)
  (let* ((stream (open (concatenate 'string directory "/bhstruct.txt") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (stream2 (open (concatenate 'string directory "/nft.txt") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (str-loopbh (open (concatenate 'string directory "/bhloops.txt") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (str-seg (open (concatenate 'string directory "/segments.txt") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (str-loop (open (concatenate 'string directory "/loops.txt") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (inf-seg (open (concatenate 'string directory "/infsegments.txt") :direction :output :if-does-not-exist :create :if-exists :supersede)))
    (loop with *allow-missing-predecessor-files* = t
	  for time from start to end by step
	  do (setf *read-dump-diamonds* nil)
	  do (read-dumps directory :time time)
	  do (bh-struct directory bhn stream stream2 str-loopbh str-seg str-loop inf-seg time tbh))
    (close stream)
    (close stream2)
    (close str-loopbh)
    (close str-seg)
    (close str-loop)
    (close inf-seg)))

(defun print-bh-single (directory start end step bhn bha)
  (let* ((stream (open (concatenate 'string directory "/bhsingle.txt") :direction :output :if-does-not-exist :create :if-exists :supersede)))
    (loop with *allow-missing-predecessor-files* = t
	  for time from start to end by step
	  do (setf *read-dump-diamonds* nil)
	  do (read-dumps directory :time time)
	  do (bh-single directory bha bhn stream time))
    (close stream)))

(defun bh-length-stats (&optional (paths (get-paths (current-time))))
  (let* ((bhloop 0)
	 (bhlong 0)
	 (free 0)
	 (bhloop-l 0.0d0)
	 (bhlong-l 0.0d0)
	 (free-l 0.0d0))
    (loop for string in paths
	  for bh-length = (path-bh-length string)
	  do (if (and (car bh-length) (car (cdr bh-length)))
		 (if (equalp (blackhole-center (car bh-length)) (blackhole-center (car (cdr bh-length))))
		     (progn
		       (setf bhloop (+ bhloop 1))
		       (when (> (car (cdr (cdr bh-length))) bhloop-l)
			 (setf bhloop-l (car (cdr (cdr bh-length))))))
		   (progn
		     (setf bhlong (+ bhlong 1))
		     (when (> (car (cdr (cdr bh-length))) bhlong-l)
                       (setf bhlong-l (car (cdr (cdr bh-length)))))))
	       (progn
		 (setf free (+ free 1))
		 (when (> (car (cdr (cdr bh-length))) free-l)
		   (setf free-l (car (cdr (cdr bh-length))))))))
    (list bhloop bhloop-l bhlong bhlong-l free free-l)))
		   
		 
  
(defun print-bh-length-stats (directory start end step)
  (let* ((stream (open (concatenate 'string directory "/bhlengthstats.txt") :direction :output :if-does-not-exist :create :if-exists :supersede)))
    (loop with *allow-missing-predecessor-files* = t
	  for time from start to end by step
	  do (read-dumps directory :time time)
	  do (format stream "~A ~{~a~^ ~} ~%" (format-float-reasonably time) (bh-length-stats)))
    (close stream)))

(defun bh-check (directory bhn)
  (with-open-file (stream (concatenate 'string directory "/blackholes.dat") :if-does-not-exist nil)
		  (let* ((center (make-3vector 0.0d0 0.0d0 0.0d0))
			 (size 0.d0)
			 (bh (make-blackhole :center center :size size))
			 (bhs nil))
		    (loop for i below bhn
			  do (setf bh (read stream nil))
			  do (push bh bhs))
		    (loop for bh1 in bhs
			  do (loop for bh2 in bhs
				   do (when (and (< (3vector-distance-wrap (blackhole-center bh1) (blackhole-center bh2)) (* 2.5 (blackhole-size bh1)))
						 (> (3vector-distance-wrap (blackhole-center bh1) (blackhole-center bh2)) 0.0))
					(format t "BH1: ~S~%" bh1)
					(format t "BH2: ~S~%" bh2)
					(error "Wrong setting of bhs")))))))

(defun bh-plot (directory bhn)
  (with-open-file (stream (concatenate 'string directory "/blackholes.dat") :if-does-not-exist nil)
                  (let* ((center (make-3vector 0.0d0 0.0d0 0.0d0))
                         (size 0.d0)
                         (bh (make-blackhole :center center :size size))
			 (str-out (open (concatenate 'string directory "/bhplot.txt") :direction :output :if-does-not-exist :create :if-exists :supersede)))
		    (loop for i below bhn
                          do (setf bh (read stream nil))
                          do (setf center (blackhole-center bh))
			  do (setf size (blackhole-size bh))
			  do (format str-out "~S ~S ~S ~S~%" (3vector-x center) (3vector-y center) (3vector-z center) size))
		    (close str-out))))

(defun point-bhs (directory i1 i2)
  (let* ((center (make-3vector 0.0d0 0.0d0 0.0d0))
	 (size 0.d0)
	 (outbhs (make-hash-table :test 'equalp))
	 (bh (make-blackhole :center center :size size))
	 (stream (open (concatenate 'string directory "/blackholes.dat") :direction :output :if-does-not-exist :create :if-exists :supersede)))
    (loop for i from i1 to i2
	  do (when (probe-file (concatenate 'string directory "/" (write-to-string i) "-blackholes.dat")) 
	       (with-open-file (stream (concatenate 'string directory "/" (write-to-string i) "-blackholes.dat") :if-does-not-exist nil)
			     (loop
			      do (setf bh (read stream nil))
			      until (null bh)
			      do (setf center (blackhole-center bh))
			      do (setf size (blackhole-size bh))
			      do (when (null (gethash (blackhole-center bh) outbhs))
				   (setf (gethash (blackhole-center bh) outbhs) bh))))))
    (format t "Number of BHS: ~S~%" (hash-table-count outbhs))
    (loop for k being each hash-key of outbhs
	  do (write (gethash k outbhs) :readably t :stream stream)
	  do (terpri stream))
    (close stream)))

;;;Dumping

(defvar *successor-streams* (make-array 5 :initial-element nil)) ;Streams to our 4 successors, or dump stream the fifth slot

(defstruct (dump-request (:include timed-event (function 'do-dump) (report-string "dump")))
  step)					;Dump goes in step-N.dat.  NIL if only lengths to be dumped

;;Request dump with given parameters.  STEP NIL if only lengths to be dumped.
(defun request-dump (step time)
  (main-calendar-add (local-time time) (make-dump-request :time time :step step)))

(define-simulate-variable *dump-average-velocities* nil) ;If set, include velocity in dump file

;;If this is set, do not include in the length any diamond that has nonzero countup, and therefore is in a loop
;;that survived one oscillation without self-intersection
(define-simulate-variable *dump-length-excludes-loops* t)

;;Bound to last explicitly given dump time for preservation code to see.  See notes.
(defvar *last-dump-time* nil)

;;Bound during dump
(defvar *dump-step*)			;Step number, or NIL for dumping only length
(defvar *dump-length*)			;Comoving invariant length of string
(defvar *dump-count*)			;Count of diamonds being dumped
(defvar *dump-total-velocity*)		;Integral of v^2 dl over strings dumped
(defvar *dump-diamonds* t)              ;Dump all the diamonds in the string



;;Do the dump according to the request
(defun do-dump (structure)
  (let ((*dump-time* (dump-request-time structure)) ;Say that we are dumping at this time
	(*dump-step* (dump-request-step structure)))
    (report-progress "[Dump ~A:" (format-float-reasonably *dump-time*))
    (if *dump-step* (dump-strings)	;Real dump?
      (dump-string-1))			;No, only wants length: just do that
    (report-progress "OK]")))

;;Dump all strings to temporary files, then rename.  This allows you to watch for the file to appear and then
;;read it without getting partial files.
(defun dump-strings ()
  (with-open-file (stream (temporary-dump-file *job-number* *dump-step*)
			  :direction :output :element-type '(unsigned-byte 8))
    (unwind-protect
	(progn
	  (setf (aref *successor-streams* dump-destination) stream) ;Install stream to enable dumping
	  (send-dump-time stream *dump-time*)
	  (dump-string-1)
	  )
      (setf (aref *successor-streams* dump-destination) nil)))
  ;;Rename to final place
  (rename-file (temporary-dump-file *job-number* *dump-step*)
	       (merge-pathnames (dump-file *job-number* *dump-step*))))


;;Process all string.  If you previously have set the dump stream in *successor-streams*, this dumps the string.
;;Otherwise the diamonds are not actually written out, but we still compute the length.
(defun dump-string-1 ()
  (let ((*dump-length* 0.0)
	(*dump-total-velocity* 0.0)
	(*dump-count* 0))
    (map-string-paths #'write-open-string #'write-string-loop)
    (report-progress "total length ~$:" *dump-length*)
    (write-string-length *dump-time* *dump-length* *dump-total-velocity*)))



;;Write out a single loop contained entirely within one processor
(defun dump-one-loop (diamond file &key overwrite)
  (with-open-file (stream file :direction :output :element-type '(unsigned-byte 8) :if-exists (if overwrite :supersede :error))
     (unwind-protect
	(progn
	  (setf (aref *successor-streams* dump-destination) stream) ;Install stream to enable dumping
	  (let ((*dump-time* (current-time))
		(*dump-length* 0.0)				    ;Make called functions happy
		(*dump-total-velocity* 0.0)
		(*dump-count* 0))
	    (send-dump-time stream *dump-time*)		     ;Start with time so reader doesn't have to know it
	    (write-string-loop-1 diamond dump-destination))) ;Dump the loop
       (setf (aref *successor-streams* dump-destination) nil))))

(defun dump-bhs ()
  (let ((transmit-bhs (if (and *all-bhs*
			     (and (> (global-job-end *total-size* *time-offset* *job-number*) *bh-start*)
				  (< (global-job-start *total-size* *time-offset* *job-number*) *bh-start*)))
			 *job-bhs*
		       *my-bhs*)))
    (format t "~%Dumping ~S bhs~%" (hash-table-count *my-bhs*))
    (loop for k being each hash-key of transmit-bhs
	  do (loop for i below 4
		do (send-myBH (aref *successor-streams* i))
		do (send-bh (aref *successor-streams* i) (gethash k transmit-bhs))))))
      

(defstruct (successor-output (:include timed-event (function 'write-final-strings) (report-string " END"))
			     (:constructor make-successor-output (time))))

;;Open streams for successor files.  Record files in *successor-files-written*.
;;We are no longer concerned with previously-existing output files.  The cost of dealing with them correctly
;;in the face of termination was high, and the value of worrying about them seemed low.
;;We use :supersede instead of :overwrite in the hope of reducing race conditions associated with
;;packets that are still on the way to the file server to output to the file.
(defun open-successor-streams ()
  (dotimes (index 4)
    (let ((file (successor-file index)))
      (setf (aref *successor-streams* index) (open file :direction :output
						   :if-does-not-exist :create
						   :if-exists :supersede
						   :element-type '(unsigned-byte 8))))))

;;Close streams.  If nothing was written, abort empty file.
(defun close-successor-streams (abort)
  (dotimes (index 4)
    (when (aref *successor-streams* index)
      (close (aref *successor-streams* index)
	     :abort (or abort (zerop (aref *send-count-diamonds* index)))))
    (setf  (aref *successor-streams* index) nil)))

(define-timer :write)

;;Write out strings at end of run.  If called from the successor-output event, we get the event
;;itself as an unneeded argument.
(defun write-final-strings (&optional event)
  (declare (ignore event))
  (find-small-segments (* monster-length-scale 0.01) 50)
  (maybe-time :output
    (account-time :write
      (write-final-strings-1))))

(defun write-final-strings-1 ()
  (let ((abort t))			;If thrown through, e.g., worker terminated during output, abort writing
    (unwind-protect
	(let (*send-previous-diamond*
	      *send-first-diamond*
	      (*dump-time* nil))	;Say not dumping
	  (fill *send-count-diamonds* 0)
	  (open-successor-streams)
	  (write-vv-output)					   ;Write out initial condition information
	  (dump-bhs)
	  (map-string-paths #'write-open-string #'write-string-loop) ;Write out strings
	  (loop for successor from 0 below 4
		for count across *send-count-diamonds*
		for stream across *successor-streams*
		do (format t "~&~D diamonds written to ~A~%" count (namestring stream))
		when *successor-file-flags*
		do (setf (logbitp successor *successor-file-flags*) (plusp count)))
	  (setq abort nil))
      (close-successor-streams abort)))
  (format t "Returning resources...") (force-output)
  (return-final-diamonds))
  
;;Write a string loop.  There are no pre-existing junctions.
(defun write-string-loop (diamond)
  (loop with start = diamond				      ;First see if single destination
	for destinations = (diamond-destinations diamond nil nil) ;See who needs to get this diamond
	when (cdr destinations)				      ;Change of destination inside this diamond
	return (write-string-1 diamond nil destinations nil nil nil) ;OK.  Start writing this to second destination
	do (setq diamond (diamond-e diamond))		      ;go eastward to next diamond
	when (eq diamond start)		;Got back to first diamond without finding any change of destination
	return (write-string-loop-1 diamond (car destinations)))) ;Loop to a single destination

;;Write a string loop to a single destination.
(defun write-string-loop-1 (diamond destination)
  (setq *send-previous-diamond* nil	;Get ready to send.  This isn't done in transmit-diamond, because
	*send-first-diamond* nil)	;we don't give it a junction
  (when (null *dump-diamonds*)
    (setf *dump-length* 0.0d0))
  (loop with start = diamond
	do (transmit-diamond destination diamond)
	do (setq diamond (diamond-e diamond))
	until (eq diamond start))
  (when (and (= destination dump-destination) (null *dump-diamonds*))
    (send-dump-length (aref *successor-streams* destination) *dump-length*)
    (setf *dump-length* 0.0))
  (transmit-loop destination)
  )


;;Write out an open string to our successors.
(defun write-open-string (start-diamond start-junction end-diamond end-junction)
  ;;Under normal circumstances, if start and end are the same, there is just this one diamond.
  ;;But it is possible to have START to the right of END.  In that case, we must loop around
  ;;and return to this diamond later.  We call that a loopback.
  (let ((loopback (and (eq start-diamond end-diamond)
		       (rejoining-junction-p start-junction) ;This doesn't happen with VV junctions
		       (rejoining-junction-p end-junction) ;because they go at the corners of the diamond.
		       (> (rejoining-junction-b start-junction) (rejoining-junction-b end-junction))))) ;start to right
    (when (null *dump-diamonds*)
      (setf *dump-length* 0.0d0))
    (write-string-1 start-diamond start-junction
		    (diamond-destinations start-diamond start-junction
					  (and (eq start-diamond end-diamond) ;If only one diamond, stop at end
					       (not loopback)
					       end-junction))
		    end-diamond end-junction loopback)))

;;Write out a string that is not a loop going to a single destination.
;;We are called with the first diamond and the whole destinations and junctions list for that diamond.
;;For open strings, START-JUNCTION is the junction at which the string starts.  Is not included in 
;;INITIAL-DESTINATIONS.
;;We start with it and then go on to the first destination.  END-DIAMOND and END-JUNCTION give the
;;end of the open string.  If DIAMOND and END-DIAMOND are the same then the string is only in this diamond, unless
;;LOOPBACK is set, in which case we exit the diamond on the right and return to it later.
;;For loops that are not all to a single destination, START-JUNCTION is NIL and END-DIAMOND is NIL.  We ignore
;;the first destination and start with the first junction in the destination's list.  Eventually we loop around
;;to the same diamond and end with the junction that we first wrote
(defun write-string-1 (start-diamond start-junction initial-destinations end-diamond end-junction loopback)
  (let ((diamond start-diamond)		;Start with first diamond
	(destinations initial-destinations) ;Start with initial list of destinations
	(first-diamond t)
	(destination nil)		;We can't write anything until we write a junction
	(done nil)			;Flag to request exit after next transmission
	this-end-junction)
    (unless start-junction		;No starting junction, so
      (pop destinations)		;Discard first destination.  We'll use it at the end.
      (unless destinations
	(error "Must either supply starting junction or have multiple destinations"))
      (setq start-junction (pop destinations))) ;Get first junction
    (setq destination (pop destinations)) ;Get first destination
    ;;State in this loop: we are about to write DIAMOND to DESTINATION.  Remaining junctions and diamonds, if any,
    ;;are on DESTINATIONS.
    (loop
     (setq this-end-junction (pop destinations)) ;Get ending junction for present transmission, or NIL
     (when (and (null this-end-junction) ;This destination is last for this diamond?
		(eq diamond end-diamond) ;And This is the last diamond?
		(not (and loopback first-diamond))) ;But if loopback we need to encounter it a second time
       (setq this-end-junction end-junction ;Yes.  Finish with last junction
	     done t))			;Done after this transmission
     (transmit-diamond destination diamond start-junction this-end-junction) ;Send this diamond piece to the successro   
     (when (and (= destination dump-destination) this-end-junction (null *dump-diamonds*))
       (setf *dump-length* 0.0))
     (cond (done (return t))		;Exit on request
	   ((and (eq diamond start-diamond) ;We're doing the first diamond again and we've already done the first piece
		 (not first-diamond))
	    (return t))			;so we're done.
	   (this-end-junction		;More work for this diamond?
	    (setq start-junction this-end-junction) ;Ending junction of last piece is start of next
	    (setq destination (pop destinations)) ;Destination for next piece
	    )				;Loop to send this diamond to the new destination
	   (t				;Go on to next diamond
	    (let ((last diamond))	;Save for error report
	      (setq diamond (diamond-e diamond))
	      (check-diamond-step diamond :east)
	      (setq destinations	;Get its destinations
		    (if (and (eq diamond start-diamond)	;Looped back to first diamond?
			     (not loopback)) ;If loopback, it's a different part of it
			initial-destinations ;Already have destinations
		      (diamond-destinations diamond nil (and (eq diamond end-diamond) end-junction)))) ;find new list
	      ;;The next diamond should start with the same destination as the old one ended with
	      (unless (eq destination (first destinations)) ;Going to new successor?
		;;There's at least one bug that can produce this error when the intersections of the successor
		;;regions lie exactly on the diamond edge.
		(error "~S started with destination ~D, whereas its predecessor, ~S ended with ~D"
		       diamond (first destinations) last destination))
	      (setq destination (pop destinations) ;New destination; remainder if any in destinations.
		    first-diamond nil
		    start-junction nil)) ;Started at beginning of diamond
	    )))))


(defun start-send-string ()
  (setq *send-first-diamond* nil
	*send-previous-diamond* nil))

;;Transmit a junction at the beginning of the string.
(defun transmit-start-junction (stream junction)
  (cond ((eq junction :deleted) 
	 (send-start-deleted stream))
        ;added by SJW
        ((eq junction :BH)
         (send-start-BH stream))
	((eq junction :BHdeleted)
	 (send-start-BHdeleted stream))
	((eq junction :BHeatit)
	 (send-start-BHeatit stream))
	((eq junction :BHpropdel)
         (send-start-BHpropdel stream))
	((vv-junction-p junction)
	 (send-start-vv-junction stream junction))
	((initial-junction-p junction)
	 (send-start-initial-junction stream junction))
	(t
	 (send-start-junction stream junction))))

(defun transmit-end-junction (stream junction)
  (cond ((eq junction :deleted)
	 (send-end-deleted stream))
        ;added by SJW
        ((eq junction :BH)
         (send-end-BH stream))
	((eq junction :BHdeleted)
	 (send-end-BHdeleted stream))
	((eq junction :BHeatit)
	 (send-end-BHeatit stream))
	((eq junction :BHpropdel)
         (send-end-BHpropdel stream))
	((vv-junction-p junction)
	 (send-end-vv-junction stream junction))
	((initial-junction-p junction)
	 (send-end-initial-junction stream junction))
	(t
	 (send-end-junction stream junction))))

(defun transmit-loop (destination)
  (let ((stream (aref *successor-streams* destination)))
    (when stream (send-loop stream))))
  

;;Transmit a diamond or piece of one.  If left-junction is set, than we're starting a new string here.
;;If right-junction is set, we end here.  (Argument names chosen to make mirror-images happy)
(defun transmit-diamond (destination diamond &optional left-junction right-junction)
  (when (and *dump-time*
	     (> (global-time (4vector-t (diamond-start diamond))) (+ *dump-time* fudge-global-coordinates)))
    (error "~S starts at later global time than dump time ~F" diamond *dump-time*))
  (let ((stream (aref *successor-streams* destination))
	(*send-offset* (destination-offset destination))) ;Offset into his coordinates
    (when stream			;Actually sending?
      (when left-junction
	(transmit-start-junction stream left-junction)
	(setq *send-previous-diamond* nil ;String begins here, so no previous
	      *send-first-diamond* nil)) ;The first diamond has not yet been sent
      (when (or (null (= destination dump-destination)) *dump-diamonds*)
	(maybe-transmit-tag destination (diamond-tag diamond)) ;Send loop tag structure if needed
	(send-diamond stream
		      (diamond-start diamond)
		      (diamond-left diamond)
		      (diamond-right diamond)
		      (diamond-end diamond)
		      (diamond-tag diamond)
		      (unless (and (diamond-sw diamond) ;If it has a SW neighbor, then the a kink has already been sent
				   (not left-junction))	;unless this diamond was the first transmitted
			(diamond-a-kink-created diamond))
		      (unless (and (diamond-nw diamond) ;If it has a NW neighbor, then the b kink has already been sent
				   (not left-junction))
			(diamond-b-kink-created diamond))
		      (diamond-countup diamond)
		      (diamond-inertp diamond)
		      )
	)
      (when (or (null (= destination dump-destination)) (and (null *dump-diamonds*) (diamond-bh diamond)) *dump-diamonds*)
	(send-bh stream (diamond-bh diamond))) ; send the bh information of the diamond
      (unless *send-first-diamond*
	(setq *send-first-diamond* diamond)) ;Now first diamond can be used by later diamonds
      (setq *send-previous-diamond* diamond)
      (if (= destination dump-destination) ;Writing to actual dump?
	  (when (> (diamond-countup diamond) 1)	;Known loop being dumped
	    (setf (tag-dumped (diamond-tag diamond)) t)) ;Tell tag about it for *loop-preservation-dump-x*
	(flet ((junction-ours-p (junction) ;See if this dump-junction is within the range we are transmitting
		 (mirror-images		;Return NIL if too far to right or left
		  (when (and right-junction
			     (rejoining-junction-p right-junction) ;Nothing can be beyond a VV junction
			     (> (rejoining-junction-b junction) (rejoining-junction-b right-junction)))	;This to R?
		    (return-from junction-ours-p nil)))
		 t))			;OK: give T
	  (mirror-images
	   (dolist (junction (created-right-dump-junctions diamond)) ;Switch from dumping to successor here
	     (when (= (rejoining-junction-right-direction junction) destination) ;Dump joined to this successor
	       (if (junction-ours-p junction)
		   ;;Tell junction to successor, so he can use in dumping.  His dump output will start here, so
		   ;;it is a left-junction for him
		   (send-note-left-junction stream junction)
		 (warn "not sending junction")
		 ))))))
      (when (and (= destination dump-destination) right-junction (null *dump-diamonds*))
	(length-diamond-transmitted destination diamond left-junction right-junction)
	(send-dump-length (aref *successor-streams* dump-destination) *dump-length*))
      (when right-junction
	(transmit-end-junction stream right-junction)
	)
      (incf (aref *send-count-diamonds* destination))))
  (when (= destination dump-destination) ;Sending (or fake sending) to dump file?
    (mirror-images
     (when left-junction		;Delete received dump-junctions when used.  Then they do not
					;lead to confusion when we try to rescale them on intersection
       (deletef left-junction (received-left-dump-junctions diamond))))) ;Remove from list if there
  (account-diamond-transmitted destination diamond left-junction right-junction) ;Count length being sent
  )


(mirror-images
;;See if one junction is to the left of another.  Real junctions that are EQ fail this test.
;;If both are rejoining junctions, we just compare the positions.
;;If JUNCTION-1 is anything but a rejoining-junction (i.e., NIL, :DELETED or a VV-JUNCTION), we take it to be at the
;;left of the diamond, whereas if JUNCTION-2 is not a rejoining-junction a VV-junction we take it to be at the right. 
;;Thus we return true in any of those cases, including both :DELETED
;;If junction-1 is a left junction and junction-2 is a right junction, these rules do what you'd expect.
;;But if junction-2 is a left junction, then if it is :DELETED, we will think it is on the extreme right,
;;so you should make sure such cases do not occur.
;;For JUNCTION-RIGHT-P, junction-1 should normally be a right junction.
(defun junction-left-p (junction-1 junction-2)
  (or (not (rejoining-junction-p junction-1)) ;If either is at end, then true
      (not (rejoining-junction-p junction-2))
      (and (not (eq junction-1 junction-2)) ;Two real junctions gives NIL
	   (< (rejoining-junction-b junction-1) (rejoining-junction-b junction-2)))))
)

(defun length-diamond-transmitted (destination diamond start-junction end-junction)
  (when (= destination dump-destination)
    (let ((start-a
	   (cond ((or (null start-junction)
		      (eq start-junction :deleted)
		      (eq start-junction :BH)
		      (eq start-junction :BHeatit)
		      (eq start-junction :BHpropdel)
		      (eq start-junction :BHdeleted))
		  (find-left-edge-position diamond (local-time *dump-time*)))
		 ((vv-junction-p start-junction)
		  1.0)
		 (t
		  (rejoining-junction-a start-junction))))
	  (end-a
	   (cond ((or (null end-junction)
		      (eq end-junction :deleted)
		      (eq end-junction :BH)
		      (eq end-junction :BHeatit)
		      (eq end-junction :BHpropdel)
		      (eq end-junction :BHdeleted))
		  (find-right-edge-position diamond (local-time *dump-time*)))
		 ((vv-junction-p end-junction)
		  0.0)
		 (t
		  (rejoining-junction-a end-junction)))))
      (let ((length 0.0d0))
	(setf length (* 2 (- start-a end-a) (4vector-t (diamond-p diamond))))
	(incf *dump-length* length)))))

;;The given range of the given diamond will be sent to the given destination
;;Included in total length, etc.
(defun account-diamond-transmitted (destination diamond start-junction end-junction)
  (when (= destination dump-destination) ;If dumping, keep track of length dumped and average velocity
    (let ((start-a			;Find valid part of diamond
	   (cond ((or (null start-junction) ;Starts at left of diamond
		      (eq start-junction :deleted) ;No previous diamond: still start at left
		      (eq start-junction :BH) ; Previous diamond BH: still start at left
		      (eq start-junction :BHeatit) ; Previous diamond BHeatit: still start at left 
		      (eq start-junction :BHpropdel) ; Previous diamond BHpropdel: still start at left   
		      (eq start-junction :BHdeleted)) ; Previous diamond BHdeleted: still start at left 
		  (find-left-edge-position diamond (local-time *dump-time*))) ;Returns values A, B
		 ((vv-junction-p start-junction) ;VV-junctions are at the left corner
		  1.0)
		 (t			;Starts with real junction: use its coordinates
		  (rejoining-junction-a start-junction))))
	  (end-a
	   (cond ((or (null end-junction) ;Ends at right of diamond
		      (eq end-junction :deleted) ;No next diamond: still ends at right
		      (eq end-junction :BH) ;Next diamond BH: still ends at right
		      (eq end-junction :BHeatit) ;Next diamond BHeatit: still ends at right 
		      (eq end-junction :BHpropdel) ;Next diamond BHpropdel: still ends at right 
		      (eq end-junction :BHdeleted)) ;Next diamond BHdeleted: still ends at right 
		  (find-right-edge-position diamond (local-time *dump-time*))) ;Returns values A, B
		 ((vv-junction-p end-junction) ;VV-junctions are at the right corner
		  0.0)
		 (t			;Ends with real junction: use its coordinates
		  (rejoining-junction-a end-junction)))))
      ;;Account length, except if this is part of an excluded loop
      (when (and (or (null *dump-length-excludes-loops*) (zerop (diamond-countup diamond))) *dump-diamonds*)
	;;Compute length in flat space approx.  See old-notes.text
	(let ((length 0.0d0)
	      ;;Compute velocity:  v = (unit p + unit q)/2.  Normalize each to 1/2, then add
	      (v (make-zero-3vector)))
	  (if (and (zerop (3vector-length (diamond-p diamond)))
		  (zerop (3vector-length (diamond-q diamond))))
	      (set v (make-zero-3vector))
	    (if (zerop (3vector-length (diamond-p diamond)))
		(setf v (3vector-normalize (diamond-q diamond) 0.5))
	      (if (zerop (3vector-length (diamond-q diamond)))
		  (setf v (3vector-normalize (diamond-p diamond) 0.5)) 
		(setf v (3vector+ (3vector-normalize (diamond-p diamond) 0.5)
				  (3vector-normalize (diamond-q diamond) 0.5))))))
	  (setf length (* 2 (- start-a end-a) (4vector-t (diamond-p diamond))))
	  (incf *dump-length* length)
	  (incf *dump-total-velocity* (* (3vector-dot v v) length)) ;We're integrating v^2 dl
	  (incf *dump-count*)
	  ))
      (when (null *dump-diamonds*)
	(let ((length 0.0d0))
	  (setf length (* 2 (- start-a end-a) (4vector-t (diamond-p diamond))))
	  (incf *dump-length* length))))))

;;This tag is about to be transmitted to DESTINATION.   Make sure it has a handle and that the destination knows 
;;about the tag.
;;Returns handle.
(defun maybe-transmit-tag (destination tag)
  (let ((stream (aref *successor-streams* destination)))
    (when stream
      (unless (logbitp destination (tag-successor-flags tag)) ;If not already sent
	(send-tag stream (get-tag-handle tag) (tag-created-position tag) ;Tell successor about it
		  (tag-last-position tag) (tag-xi tag) (tag-minimum-inert-count tag)
		  (tag-dumped tag) (tag-bh tag) )))))



  
(defun write-string-length (time length velocity)
  (when *length-output*
    (write-double time *length-output*)
    (write-double length *length-output*)
    (when *dump-average-velocities*
      (write-double velocity *length-output*))))

;;Write out information about run
(defun write-run-info-file (&key start end split-factor total-size dump-times length-times bh-size bh-number bh-times)
  (ensure-directories-exist (run-info-file))
  (with-open-file (stream (run-info-file) :direction :output )
    (prin1 (make-run-info :era *era*
			  :total-size total-size
			  :split-factor split-factor
			  :start-time start
			  :end end
			  :time-offset *time-offset*
			  :dump-times dump-times
			  :length-times length-times
			  :bh-size bh-size
			  :bh-number bh-number
			  :bh-times bh-times
			  :loop-preservation-threshold *loop-preservation-threshold*
			  :jobs *manager-jobs*)
	   stream)))

;;Find the distance between two points
(defun 3vector-distance-wrap (p1 p0)
  (3vector-distance (standardize-position p1 p0) p0))

