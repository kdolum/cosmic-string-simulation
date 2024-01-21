;;;Worker code
(in-package "CL-USER")

(define-alien-routine "fcntl" int
    (fd int)
    (cmd int)
    (arg int))

;;From /usr/include/bits/fcntl.h
(defconstant F-SETOWN 8)		;Say which process gets signal
(defconstant O-ASYNC #o20000)		;Asynchronous mode

;;Say whether socket should be asynchronous
(defun set-socket-asynchronous (socket asynchronous-p)
  (let* ((fd (socket-file-descriptor socket))
	 (flags (fcntl fd sockint::F-GETFL 0)))
    (setq flags (if asynchronous-p (logior O-ASYNC flags) ;Set bit
		  (logandc1 O-ASYNC flags))) ;clear bit
    (fcntl fd F-SETOWN (sb-unix:unix-getpid)) ;Give signal to us
    (fcntl fd sockint::F-SETFL flags)	;Give or don't give signal when input available
    ))

(defun worker-receive-retransmit ()
  (loop until (peer-receive retransmit-interval) ;Receive message: exit unless timeout
	do (peer-retransmit *manager*))
  )

(defvar *worker-command*)		;Form to evaluate to do a job
(defvar *split-factor*)		;Our split factor

(defvar *worker-sleeping* nil)		;Bound to T when sleeping
(defvar *current-job* nil)		;Job that we are doing

(defvar *backup-directory* nil)		;A directory on a different host to write backup files
(defvar *rsync-started*)		;Bound to T if this process started the server

;;(define-timer :worker)
(define-timer :worker-1)

;;List of predecessor files from completed job.  One manager has responded to JOB-DONE, we delete these.
(defvar *predecessor-files-pending-deletion* nil)

;;Top-level function in worker process for LSF
(defun worker-top-level (worker-number manager-host manager-port)
  (let ((cwd (sb-unix:posix-getcwd))
	(*rsync-started* nil))
    (format t "~%Running worker number ~D on host ~A~%~%" worker-number (machine-instance))
    (with-group-write-access
     (format t "Directory ~A~%" cwd)
     (when (eq server :uwm)
       (rename-output-files worker-number))
     (when (eq *local-data-files* :rsync)
       ;;Start an rsync process on this node unless one is running already.
       (cond ((probe-file (rsync-pid-file))
	      (format t "Not starting rsync daemon, because it is already running~%"))
	     (t
	      (format t "Starting rsync daemon~%")
	      (if *permanent-rsync-server*
		  ;;Kluge: If we start it straightforwardly, then it will be part of our process group and condor/slurm
		  ;;will kill it when we exit.
		  (do-run-program "ssh" :args (list (machine-instance) "rsync" "--daemon"
						    (format nil "--config=~A/source/rsync.conf" cwd)
						    "<" "/dev/null")) ;This prevents daemon-over-ssh protocol
		(do-run-program "rsync" :args (list "--daemon" "--config=source/rsync.conf"))
		)
	      (setq *rsync-started* t))))
     (initialize-timers)
     (time
      (worker (relative-batch-directory cwd) manager-host manager-port worker-number))
     (report-timers)
     ;;))
     )))

(defun initialize-timers ()
  (loop for (name . timer) in *timers*
	do (set timer (cons 0 0))))

(defun report-timers ()
  (format t "Waiting times: ~%")
  (loop for (name . timer) in *timers*
	for value = (symbol-value timer)
	do (format t "~& ~$ seconds in ~(~A~)~%"
		   (/ (- (car value) (cdr value)) internal-time-units-per-second)
		   name)))

;;In the vanilla universe, condor cannot append to files.  Furthermore, if the job is vacated, it will
;;restart with the same output files, deleting the old ones.  To avoid problems, we rename the open
;;output files on startup from "output" to "output.N", and similarly for "error".
(defun rename-output-files (worker-number)
  (loop with directory = (truename (worker-subdirectory "." worker-number)) ;Absolute directory where we work
	with out = (format nil "~A/output" directory) ;Absolute path of current file
	with error = (format nil "~A/error" directory)
	for serial from 1
	for out-n = (format nil "~A.~D" out serial) ;Absolute path we want
	for error-n = (format nil "~A.~D" error serial)
	while (or (probe-file out-n) (probe-file error-n))
	finally
	(rename-file out out-n)
	(rename-file error error-n)))

(defvar last-nemo-slave 772)		;Number of last existing machine

#| Not in use.  Gives mysterious unused code warnings

;;Get directory on some other host for backup
(defun find-backup-directory (subdirectory)
  (when (eq server :uwm)		;Only at Wisconsin
    (let* ((our-host (nemo-short-host (machine-instance))) ;s....
	   (host-number (parse-integer (subseq our-host 1))) ;....
	   (number host-number))
      (loop
       (when (= (incf number) last-nemo-slave) 1) ;increment and wrap through machine names
       (when (= number host-number)	;Wrapped around?
	 (error "Could not find a backup host"))
       (let ((directory (nemo-scratch-directory (format nil "s~4,'0D" number))))
	 (when (handler-case (with-timeout 30 (probe-file directory)) ;See if directory exists
		 (timeout () nil))	;Fail after 30 seconds
	   (format t "Backup directory is ~A" directory)
	   (return (format nil "~A~A" directory subdirectory))))))))
|#

;;Be a worker: open communications with manager and follow directions
;;The directory is relative to the root directories.  Log files go on the server, but
;;output files in the scratch directory
(defun worker (*worker-directory* manager-host manager-port number)
  (let* ((directory (merge-pathnames *worker-directory* batch-root-directory))
;;	 (*backup-directory* (find-backup-directory *worker-directory*))
	 (*worker-command* (with-open-file (stream (format nil "~A/command.lisp" directory))
			    (read stream)))
	 (split-factor (get-argument *worker-command* :split-factor 1)))
    (setup-split-factor split-factor)	;Install globally
    (format t "Command is ~S~%" *worker-command*)
    (setq *manager*			;Create worker structure for manager communications
	  (make-worker :number number	;Our worker number
		       :address (list (host-ent-address (get-host-by-name manager-host)) ;Address of manager
				      manager-port)
		       :status :manager
		       ))
    (setq *workers* nil)		;Say not manager
    (with-maybe-open-file (*loops-output* (get-argument *worker-command* :log) ;If requested, open log file
					  (loop-spectrum-file directory number) :direction :output :if-exists :append
					  :if-does-not-exist :create :element-type '(unsigned-byte 32))
      (with-maybe-open-file (*bh-loops-output* (get-argument *worker-command* :log)
					       (bh-loop-spectrum-file directory number) :direction :output :if-exists :append
					       :if-does-not-exist :create :element-type '(unsigned-byte 32))
        (with-maybe-open-file (*length-output* ;Check for any keyword specifying that dumps or lengths should be done
			       (get-properties (cdr *worker-command*)
					       '(:dump-interval :dump-start :dump-end :dump-times
								:length-interval :length-start :length-end :length-times))
			       (length-file directory number) :direction :output :if-exists :append
			       :if-does-not-exist :create :element-type '(unsigned-byte 64))
	(account-time :worker-1
	  (worker-1)))))))

(defun worker-1 ()
  (setq *peer-socket* (make-instance 'inet-socket :type :datagram :protocol :udp) ;Get ready to communicate
	*peer-address* (worker-address *manager*)) ;Only send to one place
  (enable-interrupt sb-unix:sigio #'worker-message-available)
  (enable-interrupt sb-unix:sigusr1 #'worker-terminate)
  (enable-interrupt sb-unix:sigint #'worker-terminate)
  (enable-interrupt sb-unix:sigterm #'worker-terminate) ;scancel, slurm preemption send this
  (unwind-protect
      (progn
	(send-worker-ready (worker-buffer *manager*) *manager* (worker-sequence *manager*)) ;Initial message
	(catch 'exit			;Throw here when finished
	  (loop (worker-receive-retransmit))))
    (when (and *rsync-started* (not *permanent-rsync-server*)) ;If we started and not leaving it around, kill it
      (format t "Killing rsync server~%")
      (do-run-program "csh" :args (list "-c" (format nil "kill `cat ~A`" (rsync-pid-file))))
      (sleep 1))			;Maybe we need to wait for it to exit cleanly before allowing slurm to kill it
    (socket-close *peer-socket*)))

;;Here on interrupt.  Receive and process message, and go back to what we were doing.
;;Interrupts are off.
(defun worker-message-available (&rest ignore)
  (declare (ignore ignore))
  (unless (main-thread-p)		;Interrupt in finalizer thread, etc.
    (return-from worker-message-available	;Pass it on
      (interrupt-thread (main-thread) #'worker-message-available)))
  (let ((conversion-data (alien-sap (make-alien (signed 8) 8)))) ;Might be in use at higher level, so make it anew
    (unless (peer-receive 0)
      ;;Somehow we get SIGIO when worker is being terminated, so don't fail here.
      (warn "No message available on interrupt that one is available"))))

;;;Synchronous communications

;;Request to do given job
;;PREDECESSOR-HOSTS is a list of length 4*n, with n=1 unless combining.
;;DO-JOB is a reply to WORKER-READY or JOB-DONE.  We are the master already.
(defun receive-do-job (sequence job random-seed predecessor-hosts)
  (when (check-sequence *manager* sequence t)
    (delete-pending-predecessor-files)	;Delete input files from previous job
    (worker-do-job job random-seed predecessor-hosts)))

(define-timer :worker-do-job)

;;Do job, tell manager what happened
;;If this called from DO-JOB, we are the master already, interrupts are on and the socket is
;;synchronous.   We do not need to receive communication from the master during the job.
;;If this is called from WAKE-JOB, we start as slave.  The socket is asynchronous because we were interrupted from
;;sleep, but interrupts are off because we're in one.  We stay in asynchronous mode while doing the job
;;so that we can acknowledge duplicate WAKE-JOB requests that happen when our first acknowledgment was lost.
;;Set socket was job is done, since we are now going to send a JOB-DONE request and become master.
(defun worker-do-job (job *random-seed* predecessor-hosts)
 (account-time :worker-do-job
  (let ((*current-job* job)
	(*successor-file-flags* 0))
    ;;Turn on interrupts in case called from WAKE-JOB.  See header comment.
    (cond ((prog1 (with-interrupts (do-job *worker-directory* job predecessor-hosts))
	     (set-socket-asynchronous *peer-socket* nil)) ;Turn off interrupts.  We are the master now.
	   (send-job-done (worker-buffer *manager*) *manager* (worker-sequence *manager*) job ;Success.  Toll manager.
			  *random-seed* *successor-file-flags*)
	   (format t "JOB-DONE sent~%"))
	  (t
	   (send-job-failed (worker-buffer *manager*) *manager* (worker-sequence *manager*) job *random-seed*))))))

(defvar *bedtime*)

;;Response to job-done is a request to sleep.  After this we are asynchronous.
;;The next message will come from the manager with a job to do.
(defun receive-worker-sleep (sequence)
  (when (check-sequence *manager* sequence t)
    (delete-pending-predecessor-files)	;Delete input files from previous job
    (let ((*worker-sleeping* t)
	  (*bedtime* (get-internal-real-time)))
      (format t "Sleeping...~%") (force-output)
      (catch 'wake-up
	(set-socket-asynchronous *peer-socket* t) ;Arrange for interrupt when we get a message
	(sleep most-positive-fixnum))	;Wait until manager wakes us
      )))
	 
;;We are all done.  Exit now.
(defun receive-worker-exit ()
  (delete-pending-predecessor-files)	;Delete input files from previous job
  (format t "~&Exiting~%")
  (throw 'exit nil))

;;;Asynchronous communications

(define-timer :sleeping)

;;Wake up from sleeping.  Interrupts are off.
;;This can also be called from top-level if the sleep was lost.
(defun receive-wake-job (sequence job random-seed predecessor-hosts)
  (when *worker-sleeping*
    (let ((slept (- (get-internal-real-time) *bedtime*)))
      (format t "Woken after ~$ seconds~%" (/ slept internal-time-units-per-second))
      (incf (car *sleeping-timer*) slept))) ;Can't use account-time, because we don't return on wake
  ;;If we are running a job (so this is duplicate), or sleeping, then we are the slave.  If we haven't got the response
  ;;to the previous job-done, then we are still the master.
  (when (check-sequence *manager* sequence (not (or *worker-sleeping* *current-job*)))
    (format t "Acknowledging job ~D~%" job)
					;Acknowledge request
    (send-acknowledge-wake-job (worker-buffer *manager*) *manager* sequence)	;Response to wake request
    (cond (*current-job*		;Already doing the job?  OK if duplicate request.
	   (unless (= job *current-job*) ;The right job? 
	     (error "Manager told us to do ~D, but we are already doing ~D" job *current-job*))) ;No, error
	  (t				;Actual wake from sleeping
	   (worker-do-job job random-seed predecessor-hosts) ;Do the job
	   (when *worker-sleeping*	;If sleeping, wake up now.  Otherwise just return.
	     (throw 'wake-up t))
	   ))))

;;Receive command to wake up and exit.  Interrupts are on.
(defun receive-wake-exit ()
  (when *worker-sleeping*
    (format t "Woken after ~$ seconds~%" (/ (- (get-internal-real-time) *bedtime*) internal-time-units-per-second)))
  (send-acknowledge-wake-exit (worker-buffer *manager*) *manager*)
  (sleep (* 2 (truncate retransmit-interval internal-time-units-per-second))) ;Wait in case acknowledgment lost
  (throw 'exit nil))
	       
;;Receive signal saying that the scheduler is killing our job.  Tell manager
(defun worker-terminate (&rest ignore)
  (declare (ignore ignore))
  (format t "~&Terminating on request by job scheduler~%")
  (unless (main-thread-p)		;Interrupt in finalizer thread, etc.
    (return-from worker-terminate	;Pass it on
      (interrupt-thread (main-thread) #'worker-terminate)))
  ;;Sometimes backtrace fails, giving recursive errors that cause premature termination
  (send-worker-terminated (worker-buffer *manager*) *manager*) ;Tell manager we are exiting
  (handler-case (print-backtrace :count 1000)
    (error (condition)
	   (format t "~&Backtrace failed~%")
	   (ignore-errors (describe condition))))
  ;;If successor files have been opened, but we have not yet finished writing them and sent JOB-DONE, delete them
  ;;because the manager will do the job again.
  (loop					;Loop here until termination acknowledged or batch system kills us
   (worker-receive-retransmit)))

(defun receive-acknowledge-terminated ()
  (throw 'exit nil))


(define-timer :probe-input)
(define-timer :retrieve)

;;Do the work for job N.  Returns true on success
;;*RANDOM-SEED* should be bound to the seed to use if we're reproducing.
;;If it is 0, we set it to a new seed.
;;On success, *predecessor-files-pending-deletion* will have input files to delete when we know success was recorded
;;and *successor-files-written* will have output files to delete if we are unable to record success.
(defun do-job (directory job-number *predecessor-hosts*)
  (when (zerop *random-seed*)		;If seed not supplied, generate one
    (setq *random-seed* (get-really-random-number)))
  (let* ((*job-number* job-number)	;Install some things globally in advance of real initialization
	 ;;Kluge: if combining, we haven't sees the information about whether predecessors generated any output
	 ;;(although we do know if predecessor job did not exist at all, so some predecessor slots will be NIL
	 ;;but others will not, even though there are in fact no files)
	 (*allow-missing-predecessor-files* (> (length *predecessor-hosts*) 4))
	 (*count-missing-predecessor-files* 0) ;But we aren't going to check because we don't know what it should be
	 ;;Split-factor already globally set
	 (*output-directory*
	  (if *local-data-files*
	      (format nil "~A~A" local-root-directory directory) ;Write to local scratch directory
	    (format nil "~A~A" batch-root-directory directory))) ;Use global directory
	 (*global-output-directory* (format nil "~A~A" batch-root-directory directory))
	 (*predecessor-files-copied* (eq *local-data-files* :rsync)) ;Tell simulation where files are
	 )
    (block success			;Return on success
      (catch 'job-failed			;Throw here on error
	(handler-bind ((error #'do-job-error))
	  (worker-check-disk-space)
	  (ensure-directories-exist (job-filename *output-directory* *job-number* ""))
	  (case *local-data-files*
	    (:rsync
	     (let ((*report-time* (get-argument *worker-command* :time))) ;Set up now so we can time retrieval
	       (account-time :retrieve
		 (maybe-time :retrieve
		   (retrieve-predecessor-files *predecessor-hosts* directory)))))
	    (:nfs
	     (account-time :probe-input
	       (unless *allow-missing-predecessor-files*
		 (wait-for-predecessor-files *predecessor-hosts*)))))
	  (eval *worker-command*)	;Do command
	  (unless (> (length *predecessor-hosts*) 4) ;If multiple predecessors, don't delete
	    (setq *predecessor-files-pending-deletion* (predecessor-files))) ;Success.  Mark input files for deletion.
	  (return-from success t)))	;Return success to caller
      ;;Here on failure
      (format t "~&********** JOB ~D FAILED **********~%" *job-number*)
      nil)))

;;Here when job gets an error.
(defun do-job-error (condition)
  (format t "~&~%Job ~D FAILED~%~A~%~%" *job-number* condition)
  (print-backtrace :count 1000)
  (throw 'job-failed nil))

;;Check that predecessor files exist
(defun wait-for-predecessor-files (predecessor-hosts)
  (when (predecessors-p)
    (loop for host in predecessor-hosts
	  for predecessor from 0
	  when host			;Should have file
	  do (loop with file = (predecessor-file predecessor)
		   for count from 0
		   when (probe-file file) return (unless (zerop count)
						   (format t "  to  ")
						   (sb-int:format-universal-time t (get-universal-time) :style :short :print-weekday nil :print-timezone nil)
						   (format t "~%"))
		   when (zerop count) 
		     do (format t "~&No file ~A.  ~%Probing remote root directory returns ~A~%Waiting for file from "
				file (probe-file (remote-local-root-directory host)))
		     and do (sb-int:format-universal-time t (get-universal-time) :style :short :print-weekday nil :print-timezone nil)
		            (force-output)
		   do (sleep 1)
		   do (format t ".") (force-output)
		   ))))

;;Copy files from predecessor's scratch file system to ours using rsync
(defun retrieve-predecessor-files (predecessor-hosts directory)
  (when (predecessors-p)
    (format t "~&Retrieving input files~%")
    (loop for predecessor below 4
	  for host in predecessor-hosts
	  when host
	  do (let ((address (inet-address-string host))
		   ;;Get old file.  We must specify job number explicitly, because geometry is not set up.
		   (old-file (let ((*input-directory* directory)) ;use directory without prefix
			       (predecessor-file predecessor :of-job-number *job-number* :copied nil
						 :local nil)))
		   (new-file (predecessor-file predecessor :copied t)))
	       (cond ((zerop (process-exit-code ;Copy old file to here
			      (do-run-program "rsync" :args (list
							     "--port" "8006"
							     "--remove-sent-files"
							     (format nil "~A::simulation/~A~A"
								     address rsync-directory-prefix old-file)
							     new-file))))
;;		      (do-run-program "ssh" :args (list address "rm" old-file)) ;If copy succeeded, delete old file
		      )
		     (t			;Failed
		      (error "Failed to copy input file ~A:~A" address old-file))
		     )))))


(defvar worker-minimum-disk-gb 0.5)

(defun worker-check-disk-space ()
  (when *local-data-files*		;Writing to local disk?
    (ensure-directories-exist local-root-directory) ;In case overall root needs to be created on each local file system
    (let ((free (free-disk-gb local-root-directory)))
      (unless (>= free worker-minimum-disk-gb)
	(error "FAILING because only ~$GB of disk space are free in ~A" free local-root-directory)))))

  


;;(define-timer :delete)

;;List of predecessor files, skipping directions where the run was empty
(defun predecessor-files ()
  (loop for host in *predecessor-hosts*
	for direction from 0
	when host
	collect (predecessor-file direction)))

;;Delete a list of files.  Do not wait for them to be deleted.
(defun delete-files (files)
;;  (account-time :delete
  (when files
    (do-run-program "rm" :wait nil :args files)))
;;  )

;;If there are pending file deletions from the previous run, request them now
(defun delete-pending-predecessor-files ()
  (delete-files *predecessor-files-pending-deletion*)
  (setq *predecessor-files-pending-deletion* nil))

#|
(define-timer :unlink)

;;Delete file quickly.  Wait if it is not there.
(defun do-delete-file (file host)
  (loop
   (multiple-value-bind (res err)
       (account-time :unlink (sb-unix:unix-unlink file)) ;Delete directly, but it doesn't seem any faster
     (when res				;success?
       (return-from do-delete-file nil))
     (if (= err sb-unix:enoent)		;No such file?
	 (wait-for-file file host)
       (sb-impl::simple-file-perror "couldn't delete ~A" file err)))))
|#

;;Diagnose disappearing file and wait for it to reappear
(defun wait-for-file (file host)
  (let ((dir (remote-local-root-directory host)))
    (cond ((probe-file dir) ;Directory exists on server, but file is absent
	   (error "File ~A disappeared even though directory ~A is present" file dir)) ;Things are unlikely to improve
	  (t
	   (format t "~A: Directory ~A vanished.  Waiting for it to reappear."
		   (sb-int:format-universal-time nil (get-universal-time)
					       :style :short :print-weekday nil :print-timezone nil)
		 dir)
	   (force-output)
	   (loop until (probe-file dir)	;Wait forever for directory to reappear
		 do (sleep 30) (format t ".") (force-output))
	   (format t "~A: Directory ~A reappeared."
		   (sb-int:format-universal-time nil (get-universal-time)
						 :style :short :print-weekday nil :print-timezone nil)
		   dir)))))
	    
