;;;Interact with batch system

(in-package "CL-USER")

;;Flags to use for batch jobs.  A list of strings such as ("-p" "mpi") or ("--qos" "dregs")
;;Set it in local-modifications.
(defvar batch-flags nil)		;for all jobs that are not simulation workers or manager
(defvar manager-batch-flags nil)	;for manager.  Supersedes batch-flags even if nil.
(defvar worker-batch-flags nil)		;for workers.  Access through function below.
;;Redefine this function for more complex behavior
(defun worker-batch-flags (worker-number)
  (declare (ignore worker-number))
  worker-batch-flags)

(defvar *batch-dynamic-space-size* (if (eq server :uwm) 1024 2048)) ;Number of megabytes (2^20) to use
(defparameter total-non-dynamic-space-size 100)			    ;Number of extra megabytes to allow beyond above

(defvar *manager-batch-dynamic-space-size* nil) ;if set, use instead of *batch-dynamic-space-size*

(defun batch-lisp-size ()
  (+ *batch-dynamic-space-size* total-non-dynamic-space-size))

;;Time limit for batch jobs.  Quantum 1 minute.  Reasonable formats are minutes,
;;hours:minutes:seconds, days-hours, days-hours:minutes.  You can't use hours:minutes, because it would be
;;interpreted as minutes:seconds.
(defvar *batch-time-limit* "7-0")	;7 days

;;Local scratch directory
(defparameter local-root-directory
  (ecase server
    (:cosmos "/strings/")
    (:tufts-lsf (format nil "/scratch2/~A/" (get-current-username)))
    (:tufts (format nil "/scratch/~A/" (get-current-username)))
    (:uwm (format nil "/localscratch/~A/" (get-current-username)))))

;;The pathname given to rsync will have simulation::, than this, then the relative pathname (e.g., "matter/1000")
;;and filename.  So this is what comes between the module pathname given in rsync.conf and local-root-directory
(defparameter rsync-directory-prefix (format nil "~A/" (get-current-username)))

;;Cache IP->host name
(defparameter host-name-table (make-hash-table :test #'equalp))

;;Get host name used for file system
(defun host-name (address)
  (or (gethash address host-name-table)
      (setf (gethash address host-name-table) (lookup-host-name address))))

(defun lookup-host-name (address)
  (let ((count 0))
    (loop
     (handler-case
	 (return (lookup-host-name-1 address)) ;Lookup address, return if found
       (try-again-error (condition)
	(when (> (incf count) 10)
	  (error "Persistent try-again-error ~S" condition))
	(format t "~&Failed to look up address ~S: trying again" address) ;Fall through to try again
	)))))

(defun lookup-host-name-1 (address)
  (let* ((name (host-ent-name (get-host-by-address address)))
	 (position (position #\. name)))
    (if position
	(subseq name 0 position)		;Only part before first dot
      name)))

;;Convert nemo-slave0662.nemo.phys.uwm.edu to s0662.
(defun nemo-short-host (hostname)
  (unless (string-equal hostname "nemo-slave" :end1 (min (length hostname) 10))
    (error "Don't know how to NFS-mount ~A" hostname))
  (format nil "s~A" (subseq hostname 10 14)))
  

;;The same directory accessed remotely through NFS
(defun remote-local-root-directory (host &key user)
  (ecase server
    (:tufts-lsf (format nil "/cluster/scratch2/~A/~A/" (host-name host) (or user (get-current-username))))
    (:uwm (nemo-scratch-directory (nemo-short-host (host-name host))))
    (t (error "Don't know path to remote local directory on server ~S" server))))

(defun nemo-scratch-directory (nemo-short-host)
  (format nil "/netdata/~A/localscratch/~A/" nemo-short-host (get-current-username)))

;;Convert real working directory of batch job into relative directory
(defun relative-batch-directory (directory)
  (when (pathnamep directory) (setq directory (format nil "~A" directory)))
  (unless (and (>= (length directory) (length batch-root-directory))
	       (string= directory batch-root-directory :end1 (length batch-root-directory)))
    (error "~A is not under ~A" directory batch-root-directory))
  (subseq directory (length batch-root-directory)))

(defun rsync-pid-file ()
  (ecase server
    (:uwm (format nil "~rsync.pid" local-root-directory))
    (:tufts "/scratch/strings-rsync.pid")))

(define-alien-routine "umask" int
    (mode int))

;;Everything happens in a directory specific to
;;this run.  On the cosmos it will be called /strings/...,
;;and on the cluster /cluster/tufts/strings/....
;;Lisp starts in this directory.
;;Under this directory there's a directory source with the source files
;;and a directory for each job called nnnn/mmm

;;Set up directories to do FORM, which must have the form (FUNCTION . KEYWORDS-AND-ARGS).
;;Various arguments are extracted from the arguments in the form.
;;Returns number of jobs used, in case we defaulted it.
;;Directory should not exist, unless overwrite is set, in which case we delete it first
(defun setup-run (directory jobs form &key overwrite combine reproduce)
  (setq directory (merge-pathnames directory batch-root-directory))
  (check-old-directory directory overwrite combine reproduce)
  (copy-files-to-run directory)
  (write-command-file directory form)
  (format t "done~%")
  jobs)

;;Sometimes format returns a base-string.  Then when you try to print it readably, you get
;;something like #A((3) BASE-CHAR . "abc").  This is readable, but the # cause trouble in things like
;;shell scripts and batch command files.  This turns base-strings into regular strings which
;;print in the usual notation with double quotes.
(defun fix-base-strings (object)
  (typecase object
    (base-string (coerce object '(simple-array character (*))))
    (cons (mapcar #'fix-base-strings object))
    (t object)))

(defun write-readably-fix-strings (object stream &rest args)
  (apply #'write (fix-base-strings object) :stream stream args))

;;Write form for lisp in top-level directory
(defun write-command-file (directory form)
  (with-open-file (stream (format nil "~A/command.lisp" directory) :direction :output :if-exists :supersede)
    (write-readably-fix-strings form stream)))

(defun read-command-file (directory)
  (with-open-file (stream (format nil "~A/command.lisp" directory))
    (read stream)))

(defun write-manager-command-file (directory form)
  (with-open-file (stream (format nil "~A/manager-commands.lisp" directory) :direction :output
			  :if-exists :append :if-does-not-exist :create)
    (write-readably-fix-strings form stream)
    (terpri stream)))

(defun check-old-directory (directory overwrite combine reproduce)
  (when (probe-file directory)		;Output directory already exists?
    (let ((delete-files (list (format nil "~A/*" directory)))) ;In simple case, give this to unix rm.
      (when (or combine reproduce)	;Only delete some files.  Must list them
	(setq delete-files
	      (remove-if #'(lambda (path)
			     (or (and combine (equal (car (last (pathname-directory path))) "splits")) ;Save splits
				 (and reproduce (equal (pathname-name path) "manager-state"))))
			 (directory (format nil "~A/*.*" directory))))) ;*.* in lisp DIRECTORY finds all files
      (when delete-files		;Not combining or combining but extra files
	(unless overwrite
	  (cerror (if combine "delete all but splits" "delete it")
		  "Directory ~A already exists" directory))
	(when *local-data-files*
	  (let ((state-file (manager-state-file directory)))
	    (when (probe-file state-file)
	      (delete-old-run (relative-batch-directory directory)))))
	(format t "~&Deleting all files ~:[~;except splits ~]in directory ~A..." combine directory) (force-output)
	(do-run-program "csh" :args (list "-c" (format nil "rm -r ~{~A ~}" delete-files)))
	(format t "done.~%")))))

;;Copy files to source directory
(defun copy-files-to-run (directory)
  (let ((source (format nil "~A/source" directory)))
    (load "load")			;Make sure files compiled
    (ensure-directories-exist (format nil "~A/" source)) ;Create source directory
    (format t "~&Copying files...") (force-output)
    (do-run-program "csh" :args (list "-c" (format nil "cp --preserve=timestamps *.lisp *.fasl *.dat rsync.conf ~A" source)))
    (do-run-program "csh" :args (list "-c" (format nil "chmod g+w ~A/*" source))) ;Make all group-writable
    (format t "done~%")))

;;When restoring job, check for new files that should perhaps be propagated to the source
;;directory of the run to be restored
(defun check-new-sources (directory)
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((files (loop for file in all-simulation-files ;Find newer sources
		     for run-date = (file-write-date (format nil "~A/source/~A.lisp" directory file))
		     for our-date = (file-write-date (format nil "~A.lisp" file))
		     when (> our-date run-date) collect file)))
    (when files				;Any?
      (format t "~&The following files have been modified since the run was started: ~{~A~^, ~}~%"
	      files)
      (when (y-or-n-p "Would you like to compile and load and copy all sources and compiled files to ~A?" directory)
	(copy-files-to-run directory)))))
	       
;;Return a list of arguments to lisp to run a job in the root directory
;;and evaluate the given form
(defun lisp-batch-arguments (form load-file)
  `("--control-stack-size" "20"
    ,@(and *batch-dynamic-space-size*
	   `("--dynamic-space-size" ,(format nil "~D" *batch-dynamic-space-size*)))
    ,@(and (eq server :uwm)		
	   `("--core" "/home/kdo/lib/sbcl/sbcl.core")) ;Core file declared in .cshrc, but no shell here, I think
    "--disable-debugger" "--load" ,load-file "--eval"
    ,(with-output-to-string (s) (write-readably-fix-strings form s :pretty nil))))

(defvar *debug-submit* nil) 			;If set, don't really submit jobs

;;Submit job
;;Arguments are:
;;DIRECTORY -- connect to this directory, relative to batch root, to do the run.  If NIL, use current directory.
;;JOB-NAME -- what to call the job
;;OUTPUT-PREFIX -- text to put before "output" for the output file.  This file will then be relative to DIRECTORY
;;FLAGS -- list of extra arguments to submit command, such as -p, --qos
;;FORM -- form to evaluate
;;Keyword arguments
;; :LOAD-FILE -- file to logo for evaluation
;; :ARRAY -- job array information.  If a number, that many jobs starting with number 0.  Otherwise a string.
(defun do-submit (&rest arguments)
  (apply (ecase server
	   (:uwm #'condor-submit)
	   (:tufts #'slurm-submit)
	   (:tufts-lsf #'lsf-submit))
	  arguments))

(defparameter all-lsf-nodes
  (append (loop for number from 1 to 56 collect (format nil "node~2,'0D" number))
	  (loop for number from 1 to 50 collect (format nil "nodeb~2,'0D" number))))

(defparameter m3-nodes
  (loop for number from 1 to 18 collect (format nil "m3n~2,'0D" number)))
(defparameter alpha-nodes
  (loop for number in '(2 3 4 7 8 9) collect (format nil "omega~3,'0D" number)))
(defparameter omega-nodes
  (loop for number in '(2 3 4 7 8 9) collect (format nil "omega~3,'0D" number)))
(defparameter all-slurm-nodes (append m3-nodes omega-nodes))

;;Not doing contributed now
;;(loop for number from 1 to 7 collect (format nil "contribb0~D" number))))

(defvar bad-nodes nil)		;If set, avoid these nodes.  Strings for tufts, integers for nemo
;;List of resource strings for bsub -R.  Examples: "scratch2 > 5000": 5G free on scratch2, i.e. nodes w/ large disks
(defvar bsub-resources nil)

;;Characters that need quoting in shell arguments
(defparameter shell-metacharaters '(#\space #\tab #\newline #\\ #\' #\" #\` #\( #\) #\!))

;;Put a backslash before every character that needs to be protected against interpretation by csh in an argument
(defun quote-shell-metacharacters (string)
  (with-output-to-string (s)
    (loop for char across string
	  when (member char shell-metacharaters) do (write-char #\\ s)
	  do (write-char char s))))

;;Submit job to LSF.
;;The output file and other things are written to the output prefix with "output", etc. appended
;;This can be either "worker-0/" or "manager-"
(defun lsf-submit (directory job-name output-prefix flags form &key (load-file "source/load"))
  (setq directory (and directory (merge-pathnames directory batch-root-directory)))
  (with-group-write-access
   (do-run-program (if *debug-submit* "echo" "bsub")
		   :args `(,(format nil "-J~A" job-name)
			   ,@flags
			   ,@(loop for string in bsub-resources
				   collect "-R"
				   collect string)
			   ,@(loop for node in bad-nodes
				   collect "-R"
				   collect (format nil "hname!=~A" node))
			   ,@(and directory `("-cwd" ,(format nil "~A" directory))) ;Run in directory if given
			   ;;Append, because we may have been running this worker number before
			   "-o" ,(format nil "~Aoutput" output-prefix) ;Output file for this worker
			   ,lisp-program
			   ,@(lisp-batch-arguments form load-file)
			   ))))

(defun slurm-submit (directory job-name output-prefix flags form &key (load-file "source/load") array)
  (setq directory (and directory (merge-pathnames directory batch-root-directory)))
  (with-group-write-access
   (multiple-value-bind (handle stream)
       (do-run-program (if *debug-submit* "echo" "sbatch")
		       :args `(,(format nil "-J~A" job-name)
			       "--mem" ,(format nil "~D" (batch-lisp-size))
			       "-t" ,*batch-time-limit*
			       ,@(and array (list "-a" (etypecase array (integer (format nil "0-~D" (1- array))) (string array))))
			       ,@(and bad-nodes (list (format nil "--exclude=~{~A,~^~}" bad-nodes)))
			       ,@(and directory (list "-D" (format nil "~A" directory))) ;Run in directory if given
			       "--open-mode=append" ;because we may have been running this worker number before
			       ,(format nil "--output=~Aoutput" output-prefix) ;Output file for this worker
			       ,@flags	;position last allows these to override things specified in other ways
			       )
		       :input :stream
		       :wait nil)
     (when *debug-submit* (setq stream t))
     ;;sbatch is now accepting a script from standard input
     (format stream "#!/bin/csh~%")
     ;;exec here causes the shell to be replaced by the program.  Then when the job is preempted, the SIGTERM
     ;;doesn't terminate the shell but instead goes to lisp, which handles it cleanly
     (format stream "exec ~A ~{~A ~}~%" lisp-program
	     ;;Since this is being processed by shell, arguments (may) need quoting
	     (mapcar #'quote-shell-metacharacters (lisp-batch-arguments form load-file)))
     (unless *debug-submit* (close stream))
     (wait-for-program handle))))

(defun condor-submit (directory job-name output-prefix flags form &key (load-file "source/load"))
  (declare (ignore job-name))
  (when flags
    (error "flags not implemented"))
  (unless directory (error "not implemented"))
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((submit-file (format nil "~A/~Acondor" directory output-prefix))
	(output-file (format nil "~Aoutput" output-prefix))
	(error-file (format nil "~Aerror" output-prefix))
	(log-file (format nil "~Alog" output-prefix)))
    (with-open-file (condor submit-file :direction :output :if-does-not-exist :create :if-exists :supersede)
      (format condor "universe = vanilla~%")
      (when bad-nodes
	(format condor "requirements = ~{regexp(\"slave~4,'0D\", Machine) == 0~^ && ~}~%" bad-nodes))
      (format condor "initialdir = ~A~%" directory) ;Run in top-level directory
      (format condor "getenv = True~%")	;copy environment
      (format condor "executable = /home/kdo/bin/sbcl~%")
      ;;The entire string goes in double quotes.  Each argument goes in single quotes, and within
      ;;each argument, each single or double quote must be duplicated.
      (format condor "arguments = \"~{'~A'~^ ~}\"~%"
	      (loop for argument in (lisp-batch-arguments form load-file)
		    collect (with-output-to-string (stream)
			      (loop for char across argument
				    do (write-char char stream)
				    when (member char '(#\' #\"))
				    do (write-char char stream)))))
      (format condor "log = ~A~%" log-file) ;Output file relative to initial directory
      (format condor "output = ~A~%" output-file)
      (format condor "error = ~A~%" error-file) ;output and error to separate files required by Condor
      (format condor "killsig = 10~%")	;SIGUSR1
      (format condor "queue~%"))
    ;;Condor isn't able to append to the error and output files, so if these files exist already
    ;;from previous runs, we should probably rename them, but we're not doing it yet.
    (do-run-program "condor_submit" :args (list submit-file))))

(defun submit-worker (directory group worker-number host port)
  (setq directory (merge-pathnames directory batch-root-directory))
  (let ((worker-output (format nil "worker-~D/" worker-number)))
    (ensure-directories-exist (format nil "~A/~A" directory worker-output))
    (do-submit directory (format nil "~D-~D" group worker-number) worker-output
	      (worker-batch-flags worker-number) `(worker-top-level ,worker-number ,host ,port))))

;;Unique name for group of jobs  
;;We store a number in a file ~/.lisp-bsub-last-name
(defun allocate-group-number ()
  (let* ((file (merge-pathnames ".lisp-bsub-last-name" (user-homedir-pathname)))
	 (next (with-open-file (stream file :if-does-not-exist nil)
		 (if stream
		     (1+ (parse-integer (read-line stream)))	;Get previous number and increment
		   0))))		;Not there, use 0
    (with-open-file (stream file :direction :output :if-exists :supersede)
      (format stream "~D~%" next))
    next))

;;Set up run and run manager here.  Duration can be T meaning we decide.
;;DIRECTORY should be relative to the cluster strings directory
;;OVERWRITE means to delete the working directory and corresponding scratch directories without asking
;;RESTORE means to continue a previous run by reading the status file.
;;REPRODUCE means to reproduce a previous run using the same seeds or those from directory given as reproduce arg
(defun manager (directory duration max-workers form
			  &rest keys &key (group (allocate-group-number)) overwrite restore reproduce submit background
			  combine)	;Combine given number of previous runs
  (when (or (> max-workers 2000) (floatp max-workers))
    (error "~D workers?  That is not reasonable" max-workers))
  (let ((*simulate-dry-run* t))		;Don't really do anything, but check arguments
    (eval form))
  ;;This is rather a kluge.  We assume that the form was given quoted and everything else is self-evaluating
  (let ((manager-form (list* 'manager directory duration max-workers `',form keys)))
    (setq directory (merge-pathnames directory batch-root-directory))
    (when (and reproduce (not (eq reproduce t)))
      (setq reproduce (merge-pathnames reproduce batch-root-directory)))
    (when (or restore reproduce (eq form t)) ;Need directory to exist already
      (unless (probe-file (or reproduce directory))
	(error "Directory ~A does not exist: Can't ~A" (or reproduce directory)
	       (cond (restore "restore") (reproduce "reproduce") (t "use previous form")))))
    (with-group-write-access
     (cond (restore
	    (check-new-sources directory) ;Make sure new code installed
	    (let ((old-form (read-command-file directory)))
	      (cond ((eq form t)	;Not reusing
		     (setq form old-form) ;Reuse
		     (format t "~&Form is ~S" form))
		    (t (unless (equalp old-form form)
			 (warn "Using new form ~S instead of old form ~S" form old-form))
		       (write-command-file directory form)))))
	   (reproduce (when (eq form t) (setq form (read-command-file (if (eq reproduce t) directory reproduce))))
		      (format t "~&Form is ~S" form))
	   ((eq form t)
	    (setq form (read-command-file directory))
	    (format t "~&Form is ~S" form)))
     (let* ((split-factor (get-argument form :split-factor 1)) ;Get some info from the form
	    (total-size (get-argument form :size (minimum-vv-simulation-size split-factor)))
	    jobs)
       (cond (restore			;Restoring old run?
	      (if (eq duration t)	;Default duration?
		  (let ((info (read-run-info-file directory)))
		    (setq jobs (run-info-jobs info))) ;Use same number of jobs as before
		(setq jobs (duration-jobs total-size split-factor duration)))
	      (write-manager-command-file directory manager-form))
	     (t				;New run
	      (when (eq duration t)	;Default duration?
		(setq duration
		      (if (get-argument form :log) ;Get default duration based on form
			  (default-duration-loops total-size (get-argument form :start 1.0)
			    (get-argument form :loop-preservation-threshold nil) (get-argument form :era :flat))
			total-size)))	;One light-crossing
	      (format t "~&Running for duration ~$~%" duration)
	      (setq jobs (duration-jobs total-size split-factor duration))
	      (setup-run directory jobs form :overwrite overwrite :combine combine ;Set up new run
			 :reproduce (eq reproduce t)) ;If T, we need to save old manager-state file
	      (write-manager-command-file directory manager-form)
	      (let ((*simulate-just-write-info* t)
		    (*manager-jobs* jobs)
		    (*output-directory* directory))
		(eval form))))		;Write run-info file
       ;if bhs has to be created create them before workers start
       (let* ((bh-size (get-argument form :bh-size nil))
	      (bh-number (get-argument form :bh-number nil))
	      (bh-start (get-argument form :bh-start nil)))
	 (when (and bh-number (null restore))
	   (sort-bhs directory total-size bh-size bh-number bh-start))) ;create blackholes.dat
       (cond (submit
	      (let ((*batch-dynamic-space-size* (or *manager-batch-dynamic-space-size* *batch-dynamic-space-size*)))
		(do-submit directory (format nil "~D manager" group) "manager-" manager-batch-flags
			   `(manager-top-level ,(format nil "~A" directory) ,group ,jobs ,split-factor ,max-workers
					       :restore ,restore :reproduce ,(and reproduce (format nil "~A" reproduce))
					       :combine ,combine))))
	     (background
	      (with-open-file (output (format nil "~A/manager-output" directory)
				      :direction :output :if-exists :append ;In case restoring
				      :if-does-not-exist :create)
		(do-run-program lisp-program :output output :input nil :wait nil
				:args (lisp-batch-arguments
				       `(manager-top-level ,(format nil "~A" directory) ,group ,jobs
							   ,split-factor ,max-workers :restore ,restore :reproduce ,reproduce :combine ,combine)
				       "load"))))
	     (t				;Run interactively, but log output
	      (with-open-file (output (format nil "~A/manager-output" directory)
				      :direction :output :if-exists :append ;In case restoring
				      :if-does-not-exist :create)
		(with-open-stream (*standard-output* (make-broadcast-stream *standard-output* output))
		  (let ((*error-output* *standard-output*)
			(*trace-output* *standard-output*))
		    (manager-top-level directory group jobs split-factor max-workers :restore restore :reproduce reproduce :combine combine)
		    )))))))))

;;Run a shell command on a node
(defun run-on-node (host command &rest do-run-program-args)
;;  (format t "~A: " host) (force-output)
  (if (eq server :tufts)
      (apply #'do-run-program "srun" :args (list "-w" host "csh" "-c" command) do-run-program-args)
    (apply #'do-run-program "ssh" :args (list "-n" (inet-address-string host) command) do-run-program-args))
  )

;;Do something on a collection of nodes by simultaneous ssh
(defun do-nodes (command &optional (nodes all-lsf-nodes))
  (let* ((count (length nodes))
	 (processes (make-array count :initial-element nil))
	 (max-fd 0)
	 (read-fds 0))
    (unwind-protect
	(progn
	  (loop for node in nodes	;Start a ssh process to each node
		for index from 0
		for process = (run-on-node node command :wait nil :input nil :output :stream :error-too t)
		do (setf (aref processes index) process)) ;Store in table
	  (loop for index below count
		for fd = (sb-sys:fd-stream-fd (process-output (aref processes index)))
		do (setf (logbitp fd read-fds) t) ;Set bits in read-fds corresponding to file descriptors
		when (> fd max-fd)
		do (setq max-fd fd))
	  (loop while (plusp read-fds)		;Loop as long as there are processes to wait for
		do
	   (multiple-value-bind (nfds readable-fds)
	       (sb-unix:unix-select (1+ max-fd) read-fds 0 0 nil) ;Wait for some stream to be readable
	     (declare (ignore nfds))
	     (loop for index below count
		   for node in nodes
		   for process = (aref processes index)
		   for stream = (process-output process)
		   for fd = (sb-sys:fd-stream-fd stream)
		   when (logbitp fd readable-fds) ;Something to say?
		   do (loop for first = t then nil
			    for char = (read-char-no-hang stream nil t)
			    while char	;Exit if nothing available
			    when (eq char t) ;EOF?
			    do (process-close process) ;Done with process
			       (setf (logbitp fd read-fds) nil)
			       (loop-finish)
			    when first do (format t "~&~A: " node)
			    do (write-char char))))))
      (loop for index below count
	    for process = (aref processes index)
	    when (and process (process-alive-p process)) ;process started, and still alive
	    do (process-kill process sb-unix:sigterm)) ;Kill it
      )))
		     


