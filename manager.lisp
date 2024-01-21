(in-package "CL-USER")

(eval-when (:compile-toplevel :load-toplevel :execute)
  (require :sb-posix))			;For SIGPWR

;;Manager code

(defvar *manager-report* nil)		;Report everything that we do

(defvar *manager-group*)		;Group number for current run.  NIL if we are not submitting jobs
(defvar *manager-directory*)		;Directory of current run
(defvar *manager-port*)			;UDP port we communicate on

(defvar *manager-workers*)		;Number of workers.  Split factor in *split-factor*

(defvar *manager-report-time*)		;Internal time at which next report is needed
(defparameter manager-status-interval (* 5 internal-time-units-per-second))

(defvar *worker-adjustment-time*)	;Internal time at which to adjust number of workers.
(defparameter worker-adjustment-interval (* 30 internal-time-units-per-second))

(defvar *idle-workers*) 		;Count of workers in :sleep
(defvar *minimum-idle-workers*)		;Minimum number of sleeping workers in the last interval
(defvar *outstanding-acknowledgments*)	;Count of workers whose acknowledgments we are awaiting

(defparameter retransmit-interval (* 2 internal-time-units-per-second))

;;Time before exited worker can be restarted.  We enforce this only by waiting to resubmit the batch job,
;;so it does not apply in systems that requeue terminated jobs automatically.
(defparameter worker-restart-delay (* 20 internal-time-units-per-second))

(defvar *first-existing-worker*)	;First worker that still exists

(defvar *manager-failing*)		;True if any job has failed.
(defvar *manager-exiting*)		;True if the run is done

(defvar *count-jobs-run*)
(defvar *count-trivial-jobs*)

(defvar *total-worker-time*)		;Total internal-real-time used by workers

(defvar *last-jobs-time*)
(defvar *last-jobs-done*)

(defvar *last-layer-time*)
(defvar *last-layer-done*)

;;A heap of jobs, for a priority queue.   The next job to run is at the top of the heap.
;;Why use a heap here instead of a calendar queue?  A calendar queue should not have too many entries on each day.
;;Since many consecutive jobs might be runnable, the day length should thus be about 1.  But we should not be 
;;enqueueing many years in advance, so the year length should be about N^4, where N is the split factor.
;;But that is not much better than just storing everything in an array.
(defvar *ready-heap*)

(defun ready-jobs () (heap-count *ready-heap*))

;;Queue of jobs that were terminated.  We can't restart them until there is enough time for the file system changes
;;to propagate to the server.  Entries in the queue are (time job), where time is the internal real-time at which
;;this job can be run.  Jobs in the queue list: DEFERRED, and the worker (which has been killed)
;;no longer points to the job.

(defvar *deferred-jobs* nil)
(defvar *deferred-jobs-tail* nil)	;Last cons in above list, for appending.  Undefined if no deferred jobs.

;;When a worker is terminated, don't try to run the same job again for three minutes.  See above.
(defparameter terminated-job-deferral-time (* 180 internal-time-units-per-second))

;;Set if terminated jobs automatically get requeued.  This happens under condor.   Currently (4/16) Tufts does not do it
(defparameter termination-requeues nil)

;;Information about a job
(defstruct (job
	    (:print-object print-job))
  number			     ;Number of this job
  status			     ;number of unfinished predecessors (0 = ready), :running, :done, :failed, :deferred
  worker			     ;Worker who is doing or did do the job
  host				     ;host (IP address) on which job was done
  (random-seed 0)		     ;32-bit random seed used for this job
  (predecessor-flags 0)		     ;Bits set if corresponding predecessor has a file for us
  )

(defun print-job (job stream)
  (if *print-escape*
      (print-unreadable-object (job stream :type t)
	(format stream "~D" (job-number job)))
    (format stream "job ~D" (job-number job))
    ))

(defvar *jobs*)

(defvar *waiting-jobs*)			;Jobs waiting on predecessors
(defvar *done-jobs*)
(defvar *minimum-ready-jobs*)		;Minimum number of ready jobs in the last interval
(defvar *status-change*)		;Set on status change to trigger a new report

(defvar *last-job-first* nil)		;If set, attempt to do later layers before earlier ones

;;Set worker status, maintain variables.  Set the job of this worker to the argument JOB, or NIL if not given.
(defun set-worker-status (worker status &optional job)
  (let ((old (worker-status worker)))
    (when (eq old :sleep)
      (decf *idle-workers*)		;No longer sleeping
      (setq *minimum-idle-workers* (min *minimum-idle-workers* *idle-workers*))) ;Maintain minimum
    (when (eq status :sleep) (incf *idle-workers*))	;Now sleeping
    (when (member old '(:starting :exiting))
      (decf *outstanding-acknowledgments*) ;Exit from state needing acknowledgment
      (setf (worker-retransmit-time worker) nil))
    (when (member status '(:starting :exiting))
      (incf *outstanding-acknowledgments*) ;New state needs acknowledgment
      (setf (worker-retransmit-time worker) ;Say when to retransmit
	    (+ (get-internal-real-time) retransmit-interval)))
    (setf (worker-status worker) status
	  (worker-job worker) job
	  *status-change* t)))

(defvar *enqueue-ready-jobs* t)		;If NIL, newly ready jobs don't go on the heap.  Used in restore.

;;Set job status, maintain variables.
;;If keywords are non-NIL, set corresponding field.  If WORKER given, default HOST. 
;;If the job is newly started it has already been removed from the heap
(defun set-job-status (job status &key worker (host (and worker (first (worker-address worker)))) random-seed)
  (let ((old (job-status job)))
    (unless (eql status old)		;Actually changing?
      (cond ((eql status 0)		;Newly ready
	     (when (numberp old)	;Normal case: previously waiting.  Could also be running if job terminated
	       (decf *waiting-jobs*))
	     (when *enqueue-ready-jobs*
	       (heap-insert *ready-heap* (job-number job))))
	    ((eql old 0)		;Job previously ready
	     (when (numberp status)	;So now it should be starting
	       (error "Job dependency count should not increase"))
	     (setq *minimum-ready-jobs* (min *minimum-ready-jobs* (ready-jobs))) ;Maintain minimum number seen
	     )
	    ((and (numberp old)		;Straight from waiting to done.  Can happen with trivial jobs
		  (eq status :done))
	     (decf *waiting-jobs*)))
      (cond ((eq status :done)		;Newly done
	     (incf *done-jobs*)))))
  (setf (job-status job) status)
  (when worker
    (setf (job-worker job) worker
	  (job-host job) host))
  (when random-seed (setf (job-random-seed job) random-seed)))

;;Tell whether we would prefer to do job-1 before job-2
(defun job-sooner-p (job-1 job-2)
  (if *last-job-first*			;This really means to do later layers first
      (multiple-value-bind (layer-1 sub-1) (truncate job-1 (expt *split-factor* 3))
	(multiple-value-bind (layer-2 sub-2) (truncate job-2 (expt *split-factor* 3))
	  (or (> layer-1 layer-2)	;Later layer is better
	      (and (= layer-1 layer-2) (< sub-1 sub-2))))) ;Within layer, earlier job number is better
    (< job-1 job-2)))

(defvar *combination-data*)		;Array of *jobs* arrays
(defvar *combination-layer*)		;Layer at which combination will take place (i.e., layers previously done)

;;Set of initial status of jobs.  If RESTORE is T read old state.
;;If RESTORE is :PARTIAL, don't complain if there are more jobs in the file
;;If REPRODUCE set, read seeds from old run
;;This does not update dependencies
(defun setup-job-status (jobs split-factor &key restore reproduce)
  (let ((*enqueue-ready-jobs* (not restore))) ;If restoring, don't put job in heap until it is read in, in case done
    (setq *jobs* (make-array jobs)
	  *ready-heap* (create-heap #'job-sooner-p)
	  *deferred-jobs* nil)
    (dotimes (job jobs)			;Make all jobs waiting
      (setf (aref *jobs* job)
	    (make-job :number job
		      :status 4)))	;We consider ourselves to have 4 successors even if split factor = 1.
    (setq *waiting-jobs* *manager-jobs*	;All waiting
	  *done-jobs* 0			;or done
	  *minimum-ready-jobs* 0)
    (cond (*combination-layer*		;Jobs up to this layer done, next layer ready
	   (dotimes (job (* *combination-layer* (expt split-factor 3)))
	     (set-job-status (aref *jobs* job) :done))
	   (loop for job from (* *combination-layer* (expt split-factor 3))
		 repeat (expt split-factor 3)
		 do (set-job-status (aref *jobs* job) 0)
		 do (setf (job-predecessor-flags (aref *jobs* job)) 15)) ;Since we don't know, set all flags
	   )				    
	  (t				;Normal case
	   (dotimes (job (min (expt split-factor 3) jobs)) ;First N^3 are ready to go
	     (set-job-status (aref *jobs* job) 0))))
    (cond (restore
	   (restore-manager-state :partial-ok (eq restore :partial)))
	  (reproduce
	   (let ((*manager-directory* (if (eq reproduce t) *manager-directory*
					reproduce))) ;Read from different directory
	     (restore-manager-state :partial-ok t :seeds-only t))))))
	   
  
(defun setup-worker-status (workers)
  (setq *manager-workers* workers
	*idle-workers* 0
	*minimum-idle-workers* 0
	*first-existing-worker* 0
	*outstanding-acknowledgments* 0
	*workers* (make-array workers)
	)
  (dotimes (n workers)
    (setf (aref *workers* n)
	  (make-worker :number n :status nil))))

;;Set up to combine previous jobs
(defun setup-combination (n)
  (when (>= (* n multiple-predecessor-job-offset-multiplier) (expt 256 job-number-length))
    (error "~S is too large" 'multiple-predecessor-job-offset-multiplier))
  (setq *combination-data* (make-array n)) ;Array of the *jobs* arrays for predecessors
  (let ((first-info nil))
    (dotimes (element n)
      (let* ((*manager-directory* (format nil "~A/splits/~D" *manager-directory* element))
	     (info (read-run-info-file *manager-directory*)))
	(cond ((null first-info)	;First time
	       (setq first-info info
		     *combination-layer* (/ (run-info-jobs info) (expt (run-info-split-factor info) 3)))
	       (when (> (run-info-jobs info) multiple-predecessor-job-offset-multiplier)
		 (error "~S is too small" 'multiple-predecessor-job-offset-multiplier))
	       (unless (numberp *combination-layer*)
		 (error "Previously done jobs not a multiple of layer size in ~S" info)))
	      ((equalp info first-info)) ;Not first-time: should be same
	      (t (error "New info ~S from split ~D not the same as old info ~S" info element first-info)))
	(let ((*combination-layer* nil) ;No special setup for subjobs
	      (*manager-jobs* (run-info-jobs info))) ;Subjobs
	  (setup-job-status *manager-jobs* *split-factor* :restore t)) ;Restore previous run
	(setf (aref *combination-data* element) *jobs*) ;Store job array
	))))

;;Put job on deferral list and set status to :deferred
(defun defer-job (job time)
  (set-job-status job :deferred)
  (let ((cons (list (list time job))))	;List of one entry
    (if *deferred-jobs*		;Some already deferred?
	(setq *deferred-jobs-tail*
	      (setf (cdr *deferred-jobs-tail*) cons)) ;Put on end, update pointer
      (setq *deferred-jobs-tail* (setq *deferred-jobs* cons)))))

;;Time to run first deferred job, if any
(defun deferred-job-time ()
  (and *deferred-jobs* (caar *deferred-jobs*)))

;;Remove and return first deferred job
(defun undefer-job ()
  (unless *deferred-jobs* (error "Attempting to undefer jobs when there are none deferred"))
  (second (pop *deferred-jobs*)))

;;Find a job to run.
(defun find-ready-job ()
  (cond (*manager-failing* :exit)	;If jobs have failed, don't start more
	(*manager-exiting* :exit)	;Already know we're finished
	((plusp (ready-jobs))		;Some job ready to go?
	 (let ((job (aref *jobs* (heap-remove *ready-heap*)))) ;Remove and return it
	   (cond ((not (eql (job-status job) 0)) ;Some unexpected state
		  (warn "~A on heap in state ~S.  Ignoring it."  job (job-status job))
		  (find-ready-job))
		 ((and (zerop (job-predecessor-flags job)) ;Job has no input files, so it is trivial
		       (> (job-layer (job-number job)) 3)) ;except in first four layers
		  (warn "Trivial ~A on heap" job)
		  ;;This shouldn't really happen, but it can when things are wrong in the manager state file
		  (manager-report "Skipping trivial ~A~%" job)
		  (set-job-status job :done) ;Set status straight to :done.  Don't put in heap.
		  (incf *count-trivial-jobs*)
		  (update-job-dependencies job *split-factor* 0)
		  (find-ready-job))
		 (t job))))
	((plusp *waiting-jobs*)		;No jobs ready, but some are waiting for predecessors
	 :wait)
	((or *deferred-jobs*		;Jobs that can't be restarted until a little time has passed
	     (and (eq server :rsync) (not *permanent-rsync-server*) ;Workers must wait until all are done
		  (< *done-jobs* *manager-jobs*)))
	 :wait)
	(t				;All jobs are done
	 (setq *manager-exiting* t)
	 :exit)))


(defun manager-report (&rest args)
  (when *manager-report*
    (apply #'format t args)))
	
;;Start manager to process the given number of workers doing the given number of jobs
;;If we should submit jobs, group should be the group number.
(defun manager-top-level (directory group jobs split-factor workers &key restore reproduce (port 0) combine)
  (format t "~&Job manager running on host ~A in process ~D: group ~D with ~D job~:P~%"
	  (machine-instance) (sb-unix:unix-getpid) group jobs)
  (setq *manager-directory* (merge-pathnames directory batch-root-directory)
	*manager-group* group)
  (setq *manager-failing* nil *manager-exiting* nil *combination-data* nil *combination-layer* nil
	*worker-adjustment-time* (and group 0))	;If we are submitting, set this up below
  (setq *manager-jobs* jobs)
  (setup-split-factor split-factor)	;Install globally, setup data structures
  (when combine				;Before setting up for our run, read predecessors to be combined
    (setup-combination combine))
  (setup-job-status jobs split-factor :restore restore :reproduce reproduce)
  (setup-worker-status workers)
  (when (zerop (ready-jobs))		;All jobs are done?
    (format t "There's nothing to do")
    (return-from manager-top-level nil))
  (setq *last-jobs-time* (get-universal-time) *last-jobs-done* (count-jobs-done)
	*last-layer-time* *last-jobs-time* *last-layer-done* (floor (first-undone-job) (expt *split-factor* 3)))
  (unwind-protect
      (progn
	(enable-interrupt sb-unix:sigusr1 #'manager-checkpoint)
	(enable-interrupt sb-posix:sigpwr #'manager-fail)
	(setq *peer-socket* (make-instance 'inet-socket :type :datagram :protocol :udp)) ;Get UDP socket
	(socket-bind *peer-socket* *inet-address-any* port) ;Call bind to get local port assigned, or use given
	(multiple-value-bind (address port) (socket-name *peer-socket*)	;Get port
	  (declare (ignore address))
	  (manager-report "Listening on ~A:~S~%" (machine-instance) port) (force-output)
	  (setq *manager-port* port)
	  (start-n-workers (ready-jobs))
	  (manager-1)))
    (enable-interrupt sb-unix:sigusr1 :default)	;Get rid of interrupt handlers
    (enable-interrupt sb-posix:sigpwr :default)	;Get rid of interrupt handler
    (socket-close *peer-socket*)))

;;Body of manager
(defun manager-1 ()
  (let ((start-time (get-internal-real-time))
	(start-runtime (get-internal-run-time))
	(*count-jobs-run* 0)
	(*count-trivial-jobs* 0)
	(*total-worker-time* 0))
    (setq *manager-report-time* start-time	;Say report is needed
	  *status-change* t)
    (catch 'manager-done		;Throw here when all workers have exited
      (loop				;Loop receiving packets, retransmitting, and doing periodic tasks
       (let ((now (get-internal-real-time)))
	 (when (plusp *outstanding-acknowledgments*) ;Anyone possibly needing retransmission?
	   (manager-do-retransmission now)) ;Do it as needed
	 (loop while (and *deferred-jobs* (>= now (deferred-job-time))) ;Jobs whose time has come to be retried?
	       do (set-job-status (undefer-job) 0))			;Set back to ready
	 (when (and *worker-adjustment-time*
		    (>= now *worker-adjustment-time*))
	   (adjust-workers now))
	 (when (and *status-change*	;Anything to say?
		    (>= now *manager-report-time*)) ;Time to say it?
	   (manager-status-report now))
	 ;;Receive and respond to all pending messages.  Do not block.
	 ;;If the manager is somehow swamped, we would not like to run forever in this loop and never gave the user
	 ;;any feedback about what is going on.  So we exit after a number of messages equal to the
	 ;;maximum number of workers.
	 (loop repeat *manager-workers* while (peer-receive 0))
	 (let ((event-time (manager-event-time))) ;Now see when next event due
	   (setq now (get-internal-real-time))	  ;Update current time
	   (when (or (null event-time) (> event-time now)) ;Need to wait?
	   (peer-receive (and event-time (- event-time now))) ;Wait until new message or time for event
	   )))))
    (manager-status-report (get-internal-real-time))
    (write-manager-state)
    (let* ((real (- (get-internal-real-time) start-time))
	   (real-seconds (floor real internal-time-units-per-second))
	   (run (- (get-internal-run-time) start-runtime)))
      (format t "~&Run ~:[complete~;FAILED ~] in ~A.  Manager CPU ~D%~%"
      *manager-failing* (format-seconds real-seconds)
      (round (* (/ run real) 100)))
      (format t "~D jobs run, ~D trivial jobs skipped~%" *count-jobs-run* *count-trivial-jobs*))
    (format t "We used ~$ worker hours to do this run~%"
	    (/ (float *total-worker-time* 0.0) internal-time-units-per-second 3600))))
			 
;;Internal real time of next event that we need to process
(defun manager-event-time ()
  (let ((time nil))
    (flet ((add-time (this)		;Consider this time, if any
	     (when this
	       (setq time (if time (min time this) this)))))
      (when (plusp *outstanding-acknowledgments*) ;Need to retransmit requests to workers?
	(loop for n from *first-existing-worker* below *manager-workers*
		       for worker = (aref *workers* n)
		       do (add-time (worker-retransmit-time worker))))	;First retransmit time
      (when *status-change*		;Something to say in next status report?
	(add-time *manager-report-time*)) ;Report needed at given time
      (add-time *worker-adjustment-time*) ;Time to adjust workers, if any
      (add-time (deferred-job-time)))		;terminated jobs waiting to be rescheduled?
    time))

(defun manager-status-report (now)
  (sb-int:format-universal-time t (get-universal-time) :style :short :print-weekday nil :print-timezone nil)
  (format t ": ")
  (report-manager-status)
  (force-output)
  (setq *manager-report-time* (+ now manager-status-interval)) ;Schedule next
  (setq *status-change* nil)		;Marked that we reported current status
  )

;;Find something for worker to do and respond to his message
;;SEQUENCE is the number of the message to which we are responding
(defun dispatch-work (worker sequence)
  (let ((job (find-ready-job)))		;Look for work to give him
    (cond ((eq job :exit)		;Every job has been started, or we are failing
	   (manager-report "Telling ~A to exit~%" worker)
	   (set-worker-status worker nil)) ;Tell him to exit
	  ((eq job :wait)		;No runnable jobs
	   (manager-report "Telling ~A to sleep~%" worker)
	   (set-worker-status worker :sleep)) ;Tell him to sleep
	  (t				;Job to run
	   (manager-report "Giving ~A to ~A~%" job worker)
	   (set-worker-status worker :running job) ;Tell him to do the job
	   (set-job-status job :running :worker worker)
	   ))
    (manager-transmit worker sequence) ;Actually send the message to the worker
    ))

;;Tell the worker what to do based on status
(defun manager-transmit (worker sequence)
  (let ((*peer-address* (worker-address worker)))
    (ecase (worker-status worker)	;What is he supposed to be doing?
      (:sleep (send-worker-sleep (worker-buffer worker) sequence)) ;Wait for further instructions
      ((nil)
       (send-worker-exit (worker-buffer worker)) ;Exit
       (account-worker-time worker)
       (check-manager-done))
      (:running (manager-transmit-job #'send-do-job worker sequence)) ;Do this job
      (:starting (manager-transmit-job #'send-wake-job worker sequence)) ;Wake up and do this job
      (:exiting (send-wake-exit (worker-buffer worker))) ;Wake up and exit
      )))

;;Send DO-JOB or WAKE-JOB to the worker.
(defun manager-transmit-job (function worker sequence)
  (let* ((job (worker-job worker))
	 (job-number (job-number job))
	 (predecessors (predecessor-job-numbers job-number))
	 (predecessor-hosts
	  (and predecessors		;If first layer, don't send anything
	       (if (and *combination-layer* (= (job-layer (job-number job)) *combination-layer*)) ;Combining now?
		   (loop for subjob-data across *combination-data*
			 ;;Append lists from each split.  If job was not done, we use NIL.  But
			 ;;if job produced no output, we don't know that.
			 nconc (loop for predecessor in predecessors
				     collect (job-host (aref subjob-data predecessor))))
		   (loop for predecessor in predecessors ;Normal case
			 for successor from 0
			 collect (and (logbitp successor (job-predecessor-flags job)) ;Have file?
				      (job-host (aref *jobs* predecessor)))))))) ;IP address of node
    (funcall function (worker-buffer worker) sequence job-number (job-random-seed job) predecessor-hosts)))

#| Without predecessor-hosts
;;Common code for starting a job
;;Need to allow rerun with random seed here
(defun manager-transmit-job (function worker sequence)
  (let* ((job (worker-job worker)))
    (funcall function (worker-buffer worker) sequence (job-number job)
	     0
	     (job-predecessor-flags job)
	     )))
|#

;;Retransmit any requests from us that have not been acknowledged
(defun manager-do-retransmission (time)
  (loop for n from *first-existing-worker* below *manager-workers*
	for worker = (aref *workers* n)
	when (and (member (worker-status worker) '(:starting :exiting)) ;Unacknowledged request from us
		  (>= time (worker-retransmit-time worker)))
	do (peer-retransmit worker)
	(setf (worker-retransmit-time worker) (+ (get-internal-real-time) retransmit-interval))) ;reschedule
  )

;;Exit if we are done.
(defun check-manager-done ()
  (loop for n from *first-existing-worker* below *manager-workers*
	for worker = (aref *workers* n)
	when (worker-status worker)	;Worker doing something
	return nil			;Not all done
	do (when (= n *first-existing-worker*) ;Maybe update variable
	     (incf *first-existing-worker*))
	finally				;If we get here, all workers have exited
	;;Is it possible to have only deferred jobs, and no workers (e.g., is also just been terminated).
	;;In that case we should not exit but start new jobs after the deferral time.
	(unless *deferred-jobs*
	  ;;Attempt to receive more messages in case some worker-exit messages got lost.  After
	  ;;timeout, we can safely exit.  If this receives a message from a worker's whose exit
	  ;;was lost, we will reenter this function, wait again, and eventually throw from the recursive call,
	  ;;rather than returning here.  If this receives some other old message, we might return
	  ;;here, in which case we should wait again.
	  (manager-report "~&~:[All jobs done~;Failed~].  Manager exiting.~%" *manager-failing*)
	  (loop while (peer-receive (* 2 retransmit-interval)))
	  (throw 'manager-done t))))

;;If there have persistently been idle workers, kill some.  If there have persistently been ready jobs, start more workers
;;We also check here for global disk space problems
(defun adjust-workers (time)
  (manager-check-disk-space)
  (cond (*manager-failing*
	 (setq *worker-adjustment-time* nil)) ;Failing: don't do this again
	(t
	 (format t "Adjusting workers: min-idle ~D, min-ready ~D~%" *minimum-idle-workers* *minimum-ready-jobs*)
	 (force-output)
	 (cond ((plusp *minimum-idle-workers*)
		(when (plusp *minimum-ready-jobs*)
		  (error "Both ready jobs and idle workers"))
		;;Keep workers alive if needed for rsync server
		(unless (and (eq *local-data-files* :rsync) (not *permanent-rsync-server*))
		  (kill-workers *minimum-idle-workers*)))
	       ((plusp *minimum-ready-jobs*)
		(let ((start (- *minimum-ready-jobs* ;Subtract jobs that have already been submitted but not yet started
				(loop for worker from *first-existing-worker* below *manager-workers* 
				      count (eq (worker-status (aref *workers* worker)) :not-started)))))
		  (when (plusp start) (start-n-workers start)))))
	 (setq *worker-adjustment-time* (+ time worker-adjustment-interval) ;Reschedule us
	       *minimum-ready-jobs* (ready-jobs) ;Reset minima for next interval
	       *minimum-idle-workers* *idle-workers*)
	 (manager-report-speed))))

;;Maximum number of jobs to submit in one go.  After this we go back to regular running, and the next 30-second
;;check will submit additional jobs if they are still needed.  This prevents slow startup where no job can
;;run until all have been submitted.
(defvar *maximum-simultaneous-submissions* 100)

;;Start up to n worker jobs
(defun start-n-workers (n)
  (setq n (min n *maximum-simultaneous-submissions*))
  (loop with now = (get-internal-real-time)
	while (plusp n)			;Anything to do?
	for index below *manager-workers*
	for worker = (aref *workers* index)
	when (and (null (worker-status worker))	;NIL = exited or never started
		  (or (null (worker-exited worker)) ;never started
		      (> (- now (worker-exited worker)) worker-restart-delay)))	;or finished long enough ago
	do (start-worker-process worker)
	   (decf n)))

(defun start-worker-process (worker)
  (format t "Submitting ~A~%" worker)
  (initialize-worker worker)
  (let ((host (machine-instance)))
    (when (string= host "hydra.phys.uwm.edu") ;Nemo nodes have a different name for this host
      (setq host "hydra"))
    (submit-worker *manager-directory* *manager-group* (worker-number worker) host *manager-port*))
  )

;;Setup worker that has just been submitted
(defun initialize-worker (worker)
  (set-worker-status worker :not-started) ;Get set to wait for initial message
  (setf (worker-sequence worker) 0)	;Reinitialize worker structure
  (setq *first-existing-worker* (min *first-existing-worker* (worker-number worker))))

(defun manager-report-speed ()
  (let ((done (count-jobs-done)))
    (unless (= done *last-jobs-done*)
      (let* ((now (get-universal-time))
	     (jobs (- done *last-jobs-done*))
	     (seconds (- now *last-jobs-time*)))
	(if (zerop seconds)		;How can this happen?
	    (warn "~S called with zero interval" 'manager-report-speed)
	  (let ((speed (/ (float jobs 0.0) seconds))
;;	        (eta (floor (/ (- *manager-jobs* done) ;Remaining jobs.  We don't account those currently running
;;			    speed)))
		)
	    (format t "Finished ~D jobs in ~A at ~$ jobs/sec~%" ;; , estimated remaining ~A~%"
		    jobs (format-seconds seconds) speed
		    )
	    (setq *last-jobs-done* done *last-jobs-time* now)))))))


;;;Requests from worker.  Sequence number should be as expected, except when acknowledgment of wake
;;message was lost, in which case check-sequence puts us in the right state.

;;Initial message from worker
(defun receive-worker-ready (worker sequence)
  (format t "Received worker-ready from ~A.  Seq ~D, expected ~D~%" worker sequence (worker-sequence worker))
  (when (check-sequence worker sequence)
    (unless (eq (worker-status worker) :not-started)
      (warn "Worker status is ~S instead of ~S" (worker-status worker) :not-started))
    (setf (worker-address worker) *peer-address*) ;Store address
    (unless (worker-started worker)	;Say when he actually started running
      (setf (worker-started worker) (get-internal-real-time)))
    (dispatch-work worker sequence)))

;;Worker finished job
(defun receive-job-done (worker sequence job-number random-seed flags)
  (when (check-sequence worker sequence)
    (let ((status (worker-status worker))	;What we think is going on
	  (job (aref *jobs* job-number)))
      (unless (member status '(:running :starting)) ;OK if he is still starting because we missed acknowledgment
	(error "~A says he finished ~A, but we thought his state was ~S" worker job status))
      (unless (eq (job-status job) :running)
	(error "~A says he finished ~A, but we thought that job was in state ~S" worker job (job-status job)))
      (manager-report "~:(~A~) finished ~A~%" worker job)
      (incf *count-jobs-run*)
      (set-job-status job :done :random-seed random-seed)		;Say job is done
      (update-job-dependencies job *split-factor* flags)
      (dispatch-work worker sequence)		;Find more work for this worker to do
      (wake-workers)			;Wake additional workers as needed
      )))

;;Job failed
(defun receive-job-failed (worker sequence job-number random-seed)
  (when (check-sequence worker sequence)
    (let ((job (aref *jobs* job-number)))
      (setq *manager-failing* t)
      (set-job-status job :failed :random-seed random-seed)	;Say job failed
      (warn "~:@(~A failed~) in ~A~%" job worker)
      (dispatch-work worker sequence)	;Go on to tell him to exit
      (kill-workers)
      )))


;;;Send asynchronous messages to worker.  
;;Sequence number should be the number of the next message and becomes the number of the expected reply

;;Asynchronous request for idle worker to do the job
(defun wake-worker (worker job)
  (manager-report "Waking ~A for ~A~%" worker job)
  (set-worker-status worker :starting job) ;Say that he is supposed to start job
  (set-job-status job :running :worker worker)
  ;;The sequence number in the table is the one we would next expect from the worker if we had sent
  ;;DO-JOB instead of SLEEP, so it is also the one for us to use next.
  (manager-transmit worker (worker-sequence worker))
  )

;;If there are both idle workers and ready jobs, wake the workers to do the jobs
;;Or, if we are done, wake them to exit
(defun wake-workers ()
  (if (or *manager-exiting* *manager-failing*)
      (kill-workers)
    (let ((n *first-existing-worker*))
      (loop while (and (plusp (ready-jobs))
		       (plusp *idle-workers*))
	    do
	    (loop until (eq (worker-status (aref *workers* n)) :sleep) ;Find an idle worker
		  do (incf n))
	    (wake-worker (aref *workers* n) (find-ready-job))
	    ))))

;;Tell N sleeping workers to exit
(defun kill-workers (&optional (n *idle-workers*))
  (loop for index from *first-existing-worker*
	for worker = (aref *workers* index)
	when (eq (worker-status worker) :sleep)
	do
	(manager-report "Waking ~A to exit~%" worker)
	(set-worker-status worker :exiting)
	(manager-transmit worker (worker-sequence worker))
	(decf n)
	while (plusp n)))

;;When JOB is done, update those that are waiting for it.
;;FLAGS tells which successors of this job got nonempty output from us.
;;If split factor is 1, we consider ourselves to have the same successor 4 times over, and so we decrement his
;;status from 4 to 0.
(defun update-job-dependencies (job split-factor flags &key (do-trivial-jobs t))
  (loop for successor below 4
	for successor-number = (successor-job-number successor (job-number job))
	when (< successor-number *manager-jobs*) ;If not beyond limit
	do (let ((successor-job (aref *jobs* successor-number)))
	     (setf (ldb (byte 1 successor) (job-predecessor-flags successor-job))
		   (ldb (byte 1 successor) flags)) ;Tell successor whether or not we have a file for him
	     (let ((new-dependencies (1- (job-status successor-job)))) ;Decrement dependencies
	       (cond ((and do-trivial-jobs
			   (zerop new-dependencies)			 ;Now ready to go?
			   (zerop (job-predecessor-flags successor-job)) ;But no input files, so it is trivial
			   ;;Always do first four layers, because new initial data could come in
			   (> (job-layer (job-number successor-job)) 3))
		      (manager-report "Skipping trivial ~A~%" successor-job)
		      (set-job-status successor-job :done) ;Set status straight to :done.  Don't put in heap.
		      (incf *count-trivial-jobs*)
		      (update-job-dependencies successor-job split-factor 0)) ;And he has no output files either
		     (t							      ;Regular case
		      (set-job-status successor-job new-dependencies)
		      ))))))

;;;Acknowledgments of asynchronous messages.
;;Receipt increments sequence number so we are ready for synchronous reply from worker

(defun receive-acknowledge-wake-job (worker sequence)
  (when (check-sequence worker sequence)
    (unless (eq (worker-status worker) :starting)
      (error "Unexpected status ~S for ~A" (worker-status worker) worker))
    (set-worker-status worker :running (worker-job worker))))

(defun receive-acknowledge-wake-exit (worker)
  (case (worker-status worker)		;What is he doing?
    (:exiting				;Exiting, as expected
     (set-worker-status worker nil)
     (account-worker-time worker)
     (check-manager-done))
    ((nil))				;Already exited: this message is duplicate
    (t (error "Unexpected status ~S for ~A" (worker-status worker) worker))))


;;Worker is being vacated by condor.  Forget whenever he is doing.
;;There's a great potential for races here.  For example, the worker could have finished and be deleting the
;;input files when it is interrupted.
;;Another problem is that the output files have already been written, so the job crashes trying to
;;write them again.
;;Using :rsync, it is quite likely that the input files will have been copied, and if we want to use
;;that mode, we will have to restart from there.
(defun receive-worker-terminated (worker)
  (send-acknowledge-terminated (worker-buffer worker))
  (handle-worker-terminated worker)
  (check-manager-done))			;This can happen if we are aborting

(defun handle-worker-terminated (worker)
  (let ((status (worker-status worker)))
    (format t "~&~:(~A~) is being terminated while ~S ~@[~A~]~%" worker status (worker-job worker))
    (when (member status '(:running :starting)) ;Doing a job now?
      (defer-job (worker-job worker) (+ (get-internal-real-time) terminated-job-deferral-time)))
    (account-worker-time worker)		;Count time used, set status back to NIL.
    (if termination-requeues			;terminated job goes back in queue and will start again?
	(initialize-worker worker)
      (set-worker-status worker nil))))	;Otherwise set back to not started



(defun report-manager-status ()
  (report-worker-status)
  (report-job-status)
  (report-layers))

(defun report-worker-status ()
  (loop for worker across *workers*
	as status = (worker-status worker)
	as job = (worker-job worker)
	count status into total
	when (eq status :running)
	count t into running
	and minimize (job-number job) into first-running
	and maximize (job-number job) into last-running
	count (eq status :sleep) into sleeping
	count (eq status :not-started) into never-started
	count (member status '(:starting :exiting)) into waiting
	finally (format t "~D worker~:P: ~D sleeping, ~D running jobs~:*~[~2*~:; (~D..~D)~]~[~:;~:*, ~D waiting for acknowledgments~]~[~:;~:*, ~D never started~]~%"
			total sleeping running first-running last-running waiting never-started)))

(defun report-job-status ()
  (format t "~D job~:P done, ~D ready, ~D waiting~[~:;~:*, ~D deferred~]~%"
			*done-jobs* (ready-jobs) *waiting-jobs* (length *deferred-jobs*)))

(defun report-layers ()
  (let ((layer-done (floor (or (first-undone-job) *manager-jobs*) ;if none undone, use total
			   (expt *split-factor* 3)))) ;Number of layers done
    (unless (= layer-done *last-layer-done*)
      (let* ((now (get-universal-time))
	     (seconds (- now *last-layer-time*)))
	(if (zerop seconds)		;How can this happen?
	    (warn "~S called with zero interval" 'manager-report-speed)
	  (let* ((speed (/ (- layer-done *last-layer-done*) (float seconds 0.0))) ;Layers per second
		 (total-layers (/ *manager-jobs* (expt *split-factor* 3)))
		 (eta (floor (/ (- total-layers layer-done) ;Layers left
				speed))))
	    (format t "Layer ~D/~D complete at ~$ layers/hour.  Estimated remaining ~A~%"
		    layer-done total-layers (* speed 3600)
		    (format-seconds eta))
	    (setq *last-layer-done* layer-done *last-layer-time* now)))))))

(defun workers-in-state (status)
  (loop for worker across *workers*
	when (eq status (worker-status worker))
	collect worker))

(defun count-jobs-done ()
  (loop for job-number from 0 below *manager-jobs*
	as job = (aref *jobs* job-number)
	as status = (job-status job)
	count (eq status :done)))

(defun first-undone-job ()
  (loop for job-number from 0 below *manager-jobs*
	as job = (aref *jobs* job-number)
	as status = (job-status job)
	unless (eq status :done) return job-number))

;;Mark all running jobs as failed
(defun manager-fail-immediately ()
  (loop for job across *jobs*
	when (member (job-status job) '(:running :starting))
	do (set-job-status job :failed)
	and count t))

;;When worker exits, count how much time he used
(defun account-worker-time (worker)
  (let ((now (get-internal-real-time)))
    (setf (worker-exited worker) now)
    (if (worker-started worker)
	(incf *total-worker-time* (- now (worker-started worker)))
      (warn "Terminating never-started ~A" worker))
    (setf (worker-started worker) nil)))

;;;The manager-state file

;;Give the status of the job.  If it is done, it doesn't matter what nonempty predecessors it had.  If it has
;;been started, there is a seed.
;;Since we just write out the jobs sequentially, we don't bother with the job number.
(define-communicators manager-state-communicators manager-state-send
  (finished-job (seed random-seed) (host host))
  (unfinished-job (seed random-seed) (flags byte) (host host))
  (unstarted-job (flags byte)))

(defun manager-state-send (code &rest arguments)
  (let* ((communicator (aref manager-state-communicators code)))
    (send-write-byte code)		;Say what message this is
    (loop for argument in arguments
	  for argument-data in (communicator-argument-data communicator)
	  for *debug-send-argument-number* from 0
	  ;;Call the sending function for each argument, according to its type
	  do (funcall (communicator-argument-send argument-data) argument))))

(defun manager-state-file (&optional (directory *manager-directory*))
  (format nil "~A/manager-state" directory))

(defun write-manager-state ()
  (let ((file (manager-state-file)))
    (format t "~&Writing state to ~A..." file) (force-output)
    (with-open-file (stream file :direction :output :if-exists :supersede
			    :if-does-not-exist :create :element-type '(unsigned-byte 8))
      (loop for job across *jobs*
	    do (case (job-status job)
		 (:done (send-finished-job stream (job-random-seed job) (job-host job)))
		 ((:running :failed)
		  (send-unfinished-job stream (job-random-seed job) (job-predecessor-flags job) (job-host job)))
		 (t			;Not started yet
		  (send-unstarted-job stream (job-predecessor-flags job)))))))
  (format t "done.~%"))

(defvar *restore-seeds-only* nil)	;If set, restore seeds and nothing else to reproduce run.

;;Called from setup-job-status to read previous state
;;If PARTIAL-OK, don't complain if we now have fewer jobs than before.
;;If SEEDS-ONLY, just restore seeds and nothing else to reproduce run.
(defun restore-manager-state (&key partial-OK ((:seeds-only *restore-seeds-only*) *restore-seeds-only*))
  (let ((file (manager-state-file)))
    (format t "~&Reading state from ~A..." file) (force-output)
    (with-open-file (*receive-source* file :element-type '(unsigned-byte 8))
      (loop for *job-number* from 0
	    for code = (receive-read-byte nil)
	    while code			;exit at EOF
	    when (= *job-number* *manager-jobs*) ;Got to end but file is not empty
	    do (unless partial-OK	;Unless caller said that was OK
		 (cerror "Forge ahead" "You are running ~D jobs, but we had old status for more.  Information about the remaining jobs will be lost" *manager-jobs*))
	    and return nil
	    do (let ((communicator (aref manager-state-communicators code)))
		 (loop for argument-data in (communicator-argument-data communicator)
		       collect (funcall (communicator-argument-receive argument-data)) into arguments
		       finally (apply (communicator-receive-function communicator) arguments)))
	    finally (when (< *job-number* *manager-jobs*) ;Fewer jobs in file 
		      (error "You are trying to run ~D jobs, but we only had old status for ~D"
			     *manager-jobs* *job-number*)))
      (format t "done.~%")
      (report-job-status))))

;;Job was previously completed.  Store that and save seed for later reproduction
(defun receive-finished-job (seed host)
  (let ((job (aref *jobs* *job-number*)))
    (unless *restore-seeds-only*
      (set-job-status job :done)
      (setf (job-host job) host)
      ;;Update dependencies.  Don't worry about flags of future jobs: they will be read in later.
      (update-job-dependencies job *split-factor* 0 :do-trivial-jobs nil))
    (setf (job-random-seed job) seed)))

;;Started but not finished.  Status is already 0 = ready.  Install seed and flags
(defun receive-unfinished-job (seed flags host)
  (let ((job (aref *jobs* *job-number*)))
    (setf (job-random-seed job) seed)
    (unless *restore-seeds-only*
      (setf (job-host job) host
	    (job-predecessor-flags job) flags)
      (unless (eql (job-status job) 0)
	(error "Unfinished ~A wasn't ready" job))
      (heap-insert *ready-heap* *job-number*)))) ;Enqueue it now

;;Not started yet.  Status already set to number of unfinished predecessors.  Install flags for nonempty predecessors.
(defun receive-unstarted-job (flags)
  (let ((job (aref *jobs* *job-number*)))
    (setf (job-predecessor-flags job) flags)
    (when (eql (job-status job) 0)	;Ready to go?
      (heap-insert *ready-heap* *job-number*)) ;Enqueue it now
    ))


#|
;;A job was started but not finished.  If the input files have been copied, we must copy them back
(defun restore-input-files (job host *split-factor*)
  (when (predecessors-p job)
    (loop with address = (inet-address-string host)
	  for predecessor below 4
	  for predecessor-job = (predecessor-job-number predecessor job)
	  for predecessor-host = (job-host (aref *jobs* predecessor-job))
	  for predecessor-address = (inet-address-string predecessor-host)
	  for old-file = (predecessor-file predecessor :of-job-number job :copied nil)
	  for new-file = (predecessor-file predecessor :of-job-number job :copied t)
	  do (format t "~&Copying ~A:~A -> ~A:~A" address new-file predecessor-address old-file) (force-output)
	  if (zerop (process-exit-code
		     (do-run-program "ssh" :args (list
						  predecessor-address ;Run on predecessor
						  "if" "(" "!" "-f" old-file ")" ;If old file does not still exist
						  "rsync" ;Get it back from new host
						  (format nil "~A:~A" address new-file)
						  old-file))))
	  do (do-run-program "ssh" :args (list address "rm" new-file)) ;If copy succeeded, delete old file
	  
	  else do (error "Copy failed")
	  )))
|#


;;Read an old manager run for processing
;;If NJOBS given, read only that many
(defun read-manager-run (directory &optional njobs)
  (setq *manager-directory* (merge-pathnames directory batch-root-directory)
	*combination-layer* nil)
  (let ((info (read-run-info-file *manager-directory*)))
    (setq *manager-jobs* (or njobs (run-info-jobs info)))
    (setup-split-factor (run-info-split-factor info))) ;Install split-factor
  (setup-job-status *manager-jobs* *split-factor* :restore :partial))

(defun manager-get-seed (directory job-number)
  (read-manager-run directory (1+ job-number))
  (let ((seed (job-random-seed (aref *jobs* job-number))))
    (values seed (format nil "#x~X" seed))))

;;Delete old files from scratch directories.
(defun delete-old-run (directory)
  (restart-case
      (let ((hosts (hosts-used directory)))
	(format t "~&Deleting old files on:~%")
	(ecase *local-data-files*
	  ((:rsync :nfs)
	   (dolist (host hosts)
	     (run-on-node host (format nil "hostname ; rm -r ~A~A/[0-9]*" local-root-directory directory))))
	   ;;For NFS, we can do it over NFS, but ssh may be faster.
#|	  (:nfs (dolist (host hosts)
		  (let ((dir (format nil "~A~A" (remote-local-root-directory host) directory)))
		    (format t "~A~%" dir)
		    (do-run-program "rm" :args (list "-r" dir))))) |#
	  ))
    (continue () :report "Forget about deleting old scratch files")))

(defun hosts-used (directory)
  (read-manager-run directory)
  (let ((hosts (make-hash-table :test #'equalp))
	(list nil))
    (loop for job across *jobs*
	  do (setf (gethash (job-host job) hosts) t))
    (maphash #'(lambda (host value)
		 (declare (ignore value))
		 (when host
		   (push host list)))
	     hosts)
    list))

;;Write out the state file on interrupt.  This is probably only useful if the workers are not actually running.
;;Otherwise, the file will immediately be obsolete.
(defun manager-checkpoint (&rest ignore)
  (declare (ignore ignore))
  (format t "~&Checkpoint requested: ")
  (write-manager-state))

;;Set failing, cleanup, and exit
(defun manager-fail (&rest ignore)
  (declare (ignore ignore))
  (warn "Failing on request by interrupt")
  (setq *manager-failing* t))

;;Copy output files from a manager job for debugging.  Directory must be relative to batch-root-directory
;;If USE-LOADED set, does not reload the run
(defun copy-predecessor-output-files (directory *job-number* &key (destination "test")
						(user (get-current-username))
						use-loaded)
  (unless use-loaded (read-manager-run directory (1+ *job-number*)))
  (let* ((*worker-directory* directory) ;Tell predecessor-file where to find it
	 (*predecessor-hosts*
	  (loop for predecessor below 4
		collect (job-host (aref *jobs* (predecessor-job-number predecessor))))))
    (loop for predecessor below 4
	  do (let* ((old (and (nth predecessor *predecessor-hosts*) ;NIL if trivial job not done
			      (predecessor-file predecessor :local :nfs :user user)))
		    (new (let ((*output-directory* destination))
			   (predecessor-file predecessor :local nil)))
		    (old-p old))
	       (when old		;Job might have been done but not written file
		 (unless (probe-file old) (setq old-p nil)))
	       (ensure-directories-exist new)
	       (format t "~&~A -> ~A~%"
		       (if old (if old-p old ;File is there
				 (format nil "~A (no such file)" old))
			 "(empty)")
		       new)
	       (do-run-program "cp" :args (list (if old-p old "/dev/null") new))))))

(defun kill-rsync-daemons (directory)
  (dolist (host (hosts-used directory))
    (do-run-program "ssh" :args (list (inet-address-string host)
				      (format nil "hostname ; kill `cat ~A`" (rsync-pid-file))))))

(defvar manager-minimum-disk-gb 5)

(defun manager-check-disk-space ()
  (let ((free (free-disk-gb *manager-directory*)))
    (unless (>= free manager-minimum-disk-gb)
      (warn "FAILING because only ~$GB of disk space are free in ~A" free *manager-directory*)
      (setq *manager-failing* t))))

#|

;;First read-manager-run
(defun find-cheap-jobs ()
  (loop with *worker-directory* = "smooth/matter-radiation/500"
	with *allow-missing-predecessor-files* = t
	with straight = 0
	with total = 0
	for *job-number* from 4463187 repeat (expt *split-factor* 3)
        for *predecessor-hosts* = (loop for predecessor in (predecessor-job-numbers *job-number*)
					collect (job-host (aref *jobs* predecessor)))
	when (loop for host in *predecessor-hosts*
		   for successor from 0
		   thereis (and host (probe-file (predecessor-file successor))))
	do (format t "~&Job ~D: " *job-number*)
	(initialize :total-size 500.0 :split-factor 29 :ijkl-origin zero-4vector :job-number *job-number*)
	(print *job-number*)
	(read-initial-strings)
	(let ((count (reduce #'+ (reduce #'append (multiple-value-list (string-counts))))))
	  (incf total count)
	  (when (cheap-job-p)
	    (incf straight count)))
	finally (format t "~%~D out of ~D were straight~%" straight total)))

;;See if job is cheap, meaning that there's only one string and it's pretty much straight
(defun cheap-job-p ()
  (block done
    (let ((string nil))
      (map-string-paths
       #'(lambda (start-diamond start-junction end-diamond end-junction)
	   (when string
	     (format t "Multiple strings.")
	     (return-from done nil))
	   (setq string (list start-diamond start-junction end-diamond end-junction)))
       #'(lambda (diamond)
	   (declare (ignore diamond))
	   (format t "Loop.")
	   (return-from done nil)))
      (cond ((apply #'string-straight-p string)
	     (format t "Straight.")
	     t)
	    (t (format t "Not straight.")
	       nil)))))

;;Get list of p's
(defun string-p-directions (start-diamond start-junction end-diamond end-junction)
  (loop for diamond = start-diamond then (diamond-e diamond)
	for first-time = nil then t
	while diamond
	until (eq diamond :deleted)
	collect (3vector-normalize (diamond-p diamond))
	until (and (eq diamond end-diamond) ;Stop at end of string
		   (not (and first-time start-junction end-junction
			     (junction-left-p end-junction start-junction)))))) ;If end is to left, must go through whole string once

(defun string-q-directions (start-diamond start-junction end-diamond end-junction)
  (loop for diamond = start-diamond then (diamond-e diamond) ;Loop over q's
	for first-time = nil then t
	while diamond
	until (eq diamond :deleted)
	collect (3vector-normalize (diamond-q diamond))
	until (and (eq diamond end-diamond) ;Stop at end of string
		   (not (and first-time start-junction end-junction
			     (junction-left-p end-junction start-junction))))))

(defun plot-string-directions (start-diamond start-junction end-diamond end-junction)
  (let ((p-list (string-p-directions start-diamond start-junction end-diamond end-junction))
	(q-list (string-q-directions start-diamond start-junction end-diamond end-junction)))
    (gnuplot 2 (length p-list)
	     #'(lambda (plot point)
	 	(if (eq point :title) (nth plot '("p" "q"))
		   (values-list (3vector-list (nth point (if (plusp plot) q-list p-list))))))
	     :3d t)))


(defun string-straight-p (start-diamond start-junction end-diamond end-junction)
  ;;Quick and dirty way to find a region that encloses all p's
  (let* ((p-list (string-p-directions start-diamond start-junction end-diamond end-junction))
	 (average-p (loop with total = zero-3vector
			  for p in p-list
			  do (setq total (3vector+ total p))
			  finally (return (3vector-normalize total))))
	 (min-cos (loop for p in p-list
			minimize (3vector-dot p average-p)))) ;Largest angle away from central vector
    (and (plusp min-cos)		;p's in a hemisphere
	 (loop for q in (string-q-directions start-diamond start-junction end-diamond end-junction)
	       always (< (3vector-dot q average-p) min-cos))
	 )))
|#
