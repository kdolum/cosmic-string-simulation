(in-package "CL-USER")

;;Expected run time for a single segment calculation is this coefficient times N_a N_b seconds
;;This should really depend on what processor we are using
(defvar backreaction-run-time-coefficient 6e-6)

(defvar backreaction-max-jobs 400)
(defvar batch-scheduling-overhead 15)

;;Always group work to take at least this time.  For quickest real-time computation, this should be small, e.g., 5
;;seconds.  But to use the least cluster resources, it should be large, e.g., 120, so not too much time is spent 
;;in dispatch.
(defvar backreaction-min-job-time 20)

;;Return number of segments to do ourselves, number for each batch job, and time in seconds that we expect each job
;;to take
(defun backreaction-bunch-numbers (Na Nb segments)
  (let* ((threads (or *threads* 1))				       ;Number of processors we can use ourselves
	 (segment-time (* Na Nb backreaction-run-time-coefficient))    ;Time to do each segment
	 (job-segments (max 1 (round backreaction-min-job-time segment-time))) ;Number of segments to do in each job
	 (job-time (* job-segments segment-time))		       ;Actual time to expect each job to take
	 (our-time (+ job-time batch-scheduling-overhead))	       ;Real time for these jobs to run
	 (our-segments (min (* threads (floor our-time segment-time))  ;How many we can do
			    segments))
	 (batch-segments (- segments our-segments)) ;How many to do in batch jobs
	 (batch-jobs (ceiling batch-segments job-segments))) ;Number of jobs to do them
    (when (plusp batch-segments)
      (setq job-segments (ceiling batch-segments batch-jobs)))  ;Distribute work evenly among jobs
    (when (> batch-jobs backreaction-max-jobs)		     ;Too many?
      ;;For simplicity, ignore overhead in this regime
      (setq batch-jobs backreaction-max-jobs
	    job-segments (ceiling segments (+ batch-jobs threads))
	    our-segments (* threads job-segments)))
    (values our-segments (and (< our-segments segments) ;Avoid returning meaningless number
			      job-segments)
	    job-time)))

;; This function submits to the cluster the calculation of one individual segment
;; of the string. It calls "compute-backreaction-segment-a or b" depending of the type
;; of segment. These functions are also introduced in backreaction.lisp and should be
;; accesible by all the nodes of the cluster.
(mirror-images
(defun submit-one-segment-backreaction-a (directory input-worldsheet-filename output-filename segment-number )
  (ensure-directories-exist directory)
  (do-submit directory
	     (format nil "br-~A~D" :a segment-number) 
	     (format nil "~A/computing-backreaction-segment-~A-~D-" directory :a segment-number) 
	     worker-batch-flags

	     `(compute-backreaction-segment-a ,directory ,input-worldsheet-filename 
					    ,output-filename ,segment-number)
	     :load-file "/cluster/home/j/b/jblanc02/strings/parallel/load"))
)


;; This function takes the information about the a-hats, b-hats etc... of the loop
;; and submits all the jobs for each individual segments to the cluster by calling
;; submit-one-segment-backreaction-a or b.
(defun submit-backreaction-loop (input-worldsheet-filename output-directory output-filename)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) 
      (read-hats-dsigmas input-worldsheet-filename)
    (declare (ignore a-dsigmas)
	     (ignore b-dsigmas))
      (mirror-images
       (loop for segment-number below (length a-hats)
	    do(submit-one-segment-backreaction-a output-directory input-worldsheet-filename 
					       (format nil "~A-~A" output-filename :a)
					       segment-number)))))
    


(mirror-images
;; These function use the information of the worldsheet to compute the correction to 
;; a particular segment.
(defun compute-backreaction-segment-a (output-directory input-worldsheet-filename output-filename segment-number)
  (multiple-value-bind (nomirror (a-hats a-dsigmas b-hats b-dsigmas))
      (read-hats-dsigmas input-worldsheet-filename)
    (let*((output-information-filename (format nil "~A/~A-~D.lisp" output-directory output-filename segment-number))
	  (segment-correction (a-segment-change a-hats a-dsigmas b-hats b-dsigmas segment-number)))		
      (with-open-file (stream output-information-filename :direction :output 
			      :if-does-not-exist :create :if-exists :supersede)
	(prin1 segment-correction stream))))))

;; This function looks for all the files that give the corrections to all the
;; segments of the worldsheet and waits until all are done to output t.
(mirror-images
(defun wait-for-backreaction-files-a (directory input-worldsheet-filename )
  (multiple-value-bind (nomirror (a-hats a-dsigmas b-hats b-dsigmas))
	(read-hats-dsigmas 
	 (format nil "~A/~A" directory input-worldsheet-filename))
   (declare 
    (ignore b-hats)
    (ignore a-dsigmas)
    (ignore b-dsigmas)
    )
   (loop for file-number below (length a-hats)
	 do (print file-number)
	 do (loop until (probe-file (format nil "~A/segment-~A-~D.lisp" directory :a file-number))
		  do (sleep 1)
		  do (format t ".") (force-output)))
   (values t)))
)

 

;; This function reads the results of the correction of one individual segment from
;; a file and puts the result in a correct format array for the rest of the computation.
(mirror-images
(defun read-a-segment-change-file (output-directory iter-a)
  (let*((a-segment-change  (make-4vector))
	(a-from-file (with-open-file (stream (format nil "~A/segment-~A-~D.lisp" output-directory :a iter-a))
		       (read stream))))
    (loop for index below 4
	  do (setf (4vector-component a-segment-change index)  (aref a-from-file index))
	  )
    (values a-segment-change)))
)



;; This is basically the identical function to "gravitational-backreaction" but using the 
;; parallel power of the cluster. It takes the information of the worldsheet, figures out the
;; jobs that need to be done to compute the correction for all the segments, submits the 
;; jobs and waits until they are all done. Then computes the new a-hats etc... as the non-parallel
;; version of this function.
(defun parallel-gravitational-backreaction (a-hats a-dsigmas b-hats b-dsigmas NGmu output-directory input-worldsheet-filename)
  (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas 
		      (format nil "~A/~A" output-directory input-worldsheet-filename))
  (submit-backreaction-loop (format nil "~A/~A" output-directory input-worldsheet-filename)
			    output-directory
			    "segment")
  (wait-for-backreaction-files-a output-directory input-worldsheet-filename)
  (wait-for-backreaction-files-b output-directory input-worldsheet-filename)
  (mirror-image-let ((new-a-hats (make-array (length a-hats)))
		     (new-a-dsigmas (make-array (length a-hats))))
    (mirror-images
     (dotimes (iter-a (length a-hats))
       (let ((d-a-prime (4vector-scale 
			 ;;(a-segment-change a-hats a-dsigmas b-hats b-dsigmas iter-a) 
			 (read-a-segment-change-file output-directory iter-a)
			 NGmu)))
	 (setf (aref new-a-hats iter-a) (3vector-normalize (3vector+ (aref a-hats iter-a) (4to3vector d-a-prime)))
	       (aref new-a-dsigmas iter-a) (* (aref a-dsigmas iter-a) (1+ (4vector-t d-a-prime))))
	 )))
    (values new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas)))


(defun submit-backreaction-loop-in-bunches (input-worldsheet-filename output-directory output-filename bunch-number)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) 
      (read-hats-dsigmas input-worldsheet-filename)
    (declare (ignore a-dsigmas)
	     (ignore b-dsigmas))
      (mirror-images
       (loop for initial-segment-number below (length a-hats) by bunch-number
	    do(submit-several-segments-backreaction-a 
	       output-directory input-worldsheet-filename 
	       (format nil "~A-~A" output-filename :a)
	       initial-segment-number (min bunch-number (- (length a-hats) initial-segment-number)))))
      ))


(mirror-images
(defun submit-several-segments-backreaction-a (output-directory input-worldsheet-filename output-filename initial-segment-number bunch-number)
  (ensure-directories-exist output-directory)
  (do-submit "/cluster/tufts/strings/JJ" 
	     (format nil "br-~A~D" :a initial-segment-number) 
	     (format nil "~A/computing-backreaction-segment-~A-~D-" output-directory :a initial-segment-number) 
	     worker-batch-flags
	     `(compute-backreaction-several-segments-a ,output-directory ,input-worldsheet-filename 
					    ,output-filename ,initial-segment-number ,bunch-number)
	     :load-file "/cluster/home/j/b/jblanc02/strings/parallel/load"))
)

(defun compute-backreaction-several-segments-a (output-directory input-worldsheet-filename output-filename initial-segment-number bunch-number)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) 
      (read-hats-dsigmas input-worldsheet-filename)
    (loop for segment-number from initial-segment-number below (+ initial-segment-number bunch-number)
	  do(let*((output-information-filename (format nil "~A/~A-~D.lisp" output-directory output-filename segment-number))
		  (segment-correction (a-segment-change a-hats a-dsigmas b-hats b-dsigmas segment-number)))		
	      (with-open-file (stream output-information-filename :direction :output 
				      :if-does-not-exist :create :if-exists :supersede)
			      (prin1 segment-correction stream))))))

(defun compute-backreaction-several-segments-b (output-directory input-worldsheet-filename output-filename initial-segment-number bunch-number)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas) 
      (read-hats-dsigmas input-worldsheet-filename)
    (loop for segment-number from initial-segment-number below (+ initial-segment-number bunch-number)
	  do(let*((output-information-filename (format nil "~A/~A-~D.lisp" output-directory output-filename segment-number))
		  (segment-correction (b-segment-change b-hats b-dsigmas a-hats a-dsigmas segment-number)))		
	      (with-open-file (stream output-information-filename :direction :output 
				      :if-does-not-exist :create :if-exists :supersede)
			      (prin1 segment-correction stream))))))

;;Obselete
(defun optimal-bunch-number (Na Nb)
  (let*((optimal-number-of-jobs (min 400 (/ (* (* 8.5 (expt 10 -5)) (+ Na Nb) (* Na Nb)) 120))))
    (ceiling (+ Na Nb) optimal-number-of-jobs)))


(mirror-images
(defun delete-backreaction-files-a (directory input-worldsheet-filename bunch-number )
  (multiple-value-bind (nomirror (a-hats a-dsigmas b-hats b-dsigmas))
	(read-hats-dsigmas 
	 (format nil "~A/~A" directory input-worldsheet-filename))
   (declare 
    (ignore b-hats)
    (ignore a-dsigmas)
    (ignore b-dsigmas)
    )
   (loop for file-number below (length a-hats)
	 do (delete-file (format nil "~A/segment-~A-~D.lisp" directory :a file-number))
	 )
   (loop for segment-number below (length a-hats) by bunch-number
	 do (delete-file (format nil "~A/computing-backreaction-segment-~A-~D-output" directory :a segment-number)) 
	 )
   (values t)))
)

 

(defun parallel-gravitational-backreaction-in-bunches (a-hats a-dsigmas b-hats b-dsigmas NGmu output-directory input-worldsheet-filename)
  (let*((Na (length a-hats))
	(Nb (length b-hats))
	(bunch-number (optimal-bunch-number Na Nb)))
    (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas 
			(format nil "~A/~A" output-directory input-worldsheet-filename))
    (submit-backreaction-loop-in-bunches (format nil "~A/~A" output-directory input-worldsheet-filename)
					 output-directory
					 "segment"
					 bunch-number)
    (wait-for-backreaction-files-a output-directory input-worldsheet-filename)
    (wait-for-backreaction-files-b output-directory input-worldsheet-filename)
    (mirror-image-let
	((a-primes (map 'vector #'(lambda (x) (3to4vector x 1.0)) a-hats))
	 (new-a-hats (make-array (length a-hats)))
	 (new-a-dsigmas (make-array (length a-hats)))
	 (d-a-prime-lengths nil))
      (mirror-images
       (dotimes (iter-a (length a-hats))
	 (let ((d-a-prime (4vector-scale 
			   ;;(a-segment-change a-hats a-dsigmas b-hats b-dsigmas iter-a) 
			   (read-a-segment-change-file output-directory iter-a)
			   NGmu)))
	   (setf (aref new-a-hats iter-a) (3vector-normalize (4to3vector (4vector+ (aref a-primes iter-a) d-a-prime)))
		 (aref new-a-dsigmas iter-a) (* (aref a-dsigmas iter-a) (1+ (4vector-t d-a-prime))))
	   (push (3vector-length (4to3vector d-a-prime)) d-a-prime-lengths)
	   )))
      (delete-backreaction-files-a output-directory input-worldsheet-filename bunch-number)
      (delete-backreaction-files-b output-directory input-worldsheet-filename bunch-number)
      (values new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas d-a-prime-lengths d-b-prime-lengths))))

;;Computation of backreaction using ODE solver

;;Structure describing process of backreaction computation.  We write one of these for each dump.
(defstruct backreaction-info
  dump-step
  time					;NGmu evolved so far
  length				;Current length of loop
  fraction-lost				;Fraction lost since the beginning
  a-kinks
  b-kinks
  agegmu				;Age of loop time Gmu.  Increases by NGmu L/2 each step
  ;;The following are since the last dump
  (simulations-with-intersections 0)	;Number of simulations that intersected
  (intersections 0)			;Count of any kind of intersection
  (rejoinings 0)
  (length-lost-backreaction 0.0)
  (length-lost-rest-frame 0.0)		;Not used, because we cannot distinguish it now that the derivatives for the ODE
					;keep us in the rest frame
  (length-lost-simulation 0.0)
  (length-lost-simulation-rest-frame 0.0)
  (threads *threads*)			;Threads used to do backreaction
  (real-time-simulation 0.0)
  (real-time-backreaction 0.0)		;Time spent doing backreaction
  (real-time-total 0.0)			;Total real time
  (batch-jobs 0)
  (derivative-calls 0)			;Number of times that ODE called for segment changes
  tolerance				;Tolerance used for ODE solution
  maximize-errors			;Flag for using maximum rather than RMS error
  backreaction-steps			;Not used
  )

(defvar *backreaction-info*)		;Current data being accumulated
(defvar *backreaction-time*)		;Where we are in backreaction process
(defvar *backreaction-directory*)	;Directory we are using
(defvar *backreaction-dump-real-time*)	;Real time of last dump
(defvar *backreaction-job-wait-time*)	;Time spent waiting for batch jobs
(defparameter fudge-backreaction-time 1e-10)

(defvar *backreaction-derivative-count*) ;Count of derivative calls.  Reset after successful integration step

(defvar *backreaction-filename-digits*)

(defun backreaction-initial-filename (directory) ;Initial data file for batch runs
  (format nil "~A/initial.dat" directory))

(defun backreaction-filename (directory step format-ctl &rest format-args) ;General scheme for filenames
  (format nil "~A/~V,'0D-~?" directory *backreaction-filename-digits* step format-ctl format-args))

(defun backreaction-data-filename (directory step) ;File giving string to be processed
  (backreaction-filename directory step "compute-loop.dat"))

(defun backreaction-segment-change-filename (directory step type index) ;File for result of a/b-segment-change
  (backreaction-filename directory step "segment-~A-~D.dat" type index))

(defun backreaction-segment-change-filename-wildcard (directory) ;File for result of a/b-segment-change
  (format nil "~A/*segment-*-*.dat" directory))

(defun backreaction-info-filename (directory step) ;Information about run
  (backreaction-filename directory step "info.lisp"))

(defun backreaction-dump-filename (directory step) ;Dump of string position
  (backreaction-filename directory step "loop.dat"))

;;Directory for simulation files.  We base it on an overall step count to avoid races, probably needlessly.
(defun backreaction-simulate-directory (directory step)
  (backreaction-filename directory step "simulation" step))

;;String before simulation.  DUMP-STEP is the step number of the the previous dump.  Since then we advanced
;;and are now simulating.
(defun backreaction-pre-simulate-filename (directory dump-step intersection-step)
  (backreaction-filename directory dump-step "simulation-~D-input.dat" intersection-step))

(defun backreaction-post-simulate-filename (directory dump-step intersection-step) ;String after simulation if no dump
  (backreaction-filename directory dump-step "simulation-~D-result.dat" intersection-step))

(defun backreaction-simulate-info-filename (directory dump-step intersection-step) ;Information after simulation
  (backreaction-filename directory dump-step "simulation-~D-result-info.lisp" intersection-step))

(defun backreaction-momentum-filename (directory step) ;Momentum 4vector
  (backreaction-filename directory step "momentum.dat"))

(defun backreaction-success-filename (directory)
  (format nil "~A/backreaction-success" directory))

;;Short form of directory for job names, etc.  Uses last ncomponents of directory name
;;only keep last nchars characters.
(defun backreaction-directory-name (directory &optional (ncomponents 2))
  (let ((list (pathname-directory directory))) ;directory structure if trailing /
    (when (pathname-name directory)	       ;another component?
      (setq list (append list  (list (pathname-name directory)))))
    (when (keywordp (car list))		; :relative, etc.?
      (pop list))
    (format nil "~{~A~^/~}" (last list ncomponents))))

;;Decide the maximum time for which we could run and the corresponding number of digits to use in filenames
;;Each oscillation the loop loses length Gamma G mu L/2.  After n oscillations, the length is L_0 exp(- Gamma G mu n/2)
;;so the number of oscillations to get to length L obeys n G mu = ln(L_0/L) 2/Gamma
(defun choose-backreaction-duration (total-loss NGmu-step)
  (let* ((min-gamma 39.0)		;Gamma must be at least this
	 (max-ngmu (* (- (log (- 1 total-loss))) (/ 2 min-gamma)))
	 (steps (ceiling max-ngmu NGmu-step)))
    (values max-ngmu (ceiling (log (1+ steps) 10)))))

(defun ode-gravitational-backreaction (a-hats a-dsigmas b-hats b-dsigmas directory
							    total-loss NGmu-dump-step ;and simulate
							    &key submit threads (tolerance 0.01) step-size
							    (maximize-errors t) no-intersections
							    minimum-loop-count resume)
  (with-group-write-access
   (setq directory (merge-pathnames directory batch-root-directory))
   (cond (resume
	  (unless (probe-file directory)
	    (error "Directory ~A does not exist" directory))
	  (check-new-sources directory)
	  (let ((files (directory (backreaction-segment-change-filename-wildcard directory))))
	    (when files
	      (format t "Deleting ~D old segment files..." (length files))) (force-output)
	    (mapc #'delete-file files)))
	 (t (when (probe-file directory)		;Output directory already exists?
	      (cerror "Delete it" "Directory ~A already exists" directory)
	      (format t "Deleting old files...") (force-output)
	      (do-run-program "csh" :args (list "-c" (format nil "rm -r ~A/*" directory)))
	      (terpri))
	    (copy-files-to-run directory)))
   (write-command-file directory -)	;Quick and dirty way to write form that invoked us
   (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
     (rest-frame-close-hats a-hats a-dsigmas b-hats b-dsigmas)) ;Make sure closed before starting
   (cond (submit
	  (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas
			      (backreaction-initial-filename directory))
	  (do-submit directory (format nil "br ~A" (backreaction-directory-name directory))
		     "backreaction-"
		     (append batch-flags (and threads `("--cpus-per-task" ,(format nil "~D" threads))))
		     `(ode-gravitational-backreaction-1
		       nil nil nil nil nil ,total-loss ,NGmu-dump-step ,tolerance ,step-size ,maximize-errors
		       ,no-intersections ,minimum-loop-count ,resume)))
	 (t (ode-gravitational-backreaction-1
	     a-hats a-dsigmas b-hats b-dsigmas directory total-loss NGmu-dump-step
	     tolerance step-size maximize-errors no-intersections minimum-loop-count resume)))))

(defun ode-gravitational-backreaction-1 (a-hats a-dsigmas b-hats b-dsigmas *backreaction-directory*
							      total-loss dump-step-size tolerance
							      step-size maximize-errors no-intersections
							      *minimum-loop-count* resume)
  (with-group-write-access
   (format t "~&Backreaction running on host ~A in process ~D with~:[out~; ~:*~D~] threads~%"
	   (machine-instance) (sb-unix:unix-getpid) *threads*)
   (unless *backreaction-directory*			     ;NIL means connected directory
     (setq *backreaction-directory* (sb-posix:getcwd)))	     ;Using "." causes some trouble with batch job submission
   (unless a-hats					     ;Batch job?
     (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas) ;Read file
       (read-hats-dsigmas (backreaction-initial-filename *backreaction-directory*))))
   (multiple-value-bind (max-ngmu *backreaction-filename-digits*)
       (choose-backreaction-duration total-loss dump-step-size)
     (let* ((*backreaction-info* (make-backreaction-info)) ;Place to accumulate energy losses
	    (ode-time 0.0)				   ;Starting time of ODE solver
	    (intersection-step 0)			   ;Filename counter for intersections
	    (initial-length (total-sigma a-dsigmas))
	    (length initial-length)
	    (length-goal (* initial-length (- 1 total-loss))) ;Stop at first dump past this
	    (agegmu 0.0)				      ;Age of loop times Gmu
	    (dump-step 0)				      
	    (simulate-step 0)		;Step counter for simulations
	    (total-derivative-count 0)
	    (*backreaction-dump-real-time* (get-internal-real-time))
	    (start-time (get-universal-time))
	    (*backreaction-job-wait-time* 0))
       (format t "Backreaction starting at ~A~%"
	   (sb-int:format-universal-time nil start-time :style :short :print-weekday nil :print-timezone nil))
       (cond (resume
	      (loop while (probe-file (backreaction-dump-filename *backreaction-directory* (1+ dump-step)))
		    do (incf dump-step))
	      (format t "Resuming from step ~D~%" dump-step)
	      (let ((info (with-open-file (stream (backreaction-info-filename *backreaction-directory* dump-step))
			    (read stream))))
		(setq agegmu (backreaction-info-agegmu info)
		      length (backreaction-info-length info)
		      ode-time (backreaction-info-time info)))
	      (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas) ;Read old output file
		(read-hats-dsigmas (backreaction-dump-filename *backreaction-directory* dump-step))))
	     (t (dump-backreaction dump-step ode-time agegmu a-hats a-dsigmas b-hats b-dsigmas length
				    tolerance maximize-errors))) ;Dump initial conditions
       (block outer
	 (loop							     ;Return here after intersection to start over
	  (setq *na* (length a-hats))				     ;Must set globally for threads
	  (let* ((y (flatten-hats a-hats a-dsigmas b-hats b-dsigmas)) ;initial point
		 (ny (length y))				      ;3(Na+Nb)
		 (atol (make-array ny :element-type 'double-float))
		 (steps-left (ceiling (- max-ngmu ode-time) dump-step-size)) ;Steps to reach ending time
		 (end-time (+ ode-time (* steps-left dump-step-size)))	     ;Round up to integer number of steps
		 ;;Count of dumps.  Will increment before using, except if there is an intersection
		 ;;in the initial conditions, which should be impossible.
		 (intersection nil)	;flag to exit dop853
		 (*backreaction-derivative-count* 0))
	    ;;Absolute tolerance for each component is tolerance times the length of the da or db of which this
	    ;;component is a part.  We don't have a way to update this as the calculation progresses, but probably
	    ;;the lengths of these factors don't change in order of magnitude
	    (loop for i below ny by 3
		  for len = (sqrt (loop for j from i repeat 3 sum (expt (aref y j) 2))) ;length of da or db
		  do (loop for j from i repeat 3
			   do (setf (aref atol j) (* tolerance len))))
	    ;;Solve as ODE.  The solver will break up the entire NGmu into sections where it can use
	    ;;its approximation.  It calls us back at the times we specify, and then we simulate.  If there
	    ;;is an intersection, we abandon the rest of the ODE solution and start again.
	    (format t "Starting ODE solution at time ~A: Length ~S lost ~1$%~%"
		    ode-time length (* 100 (- 1.0 (/ length initial-length))))
	    (setq step-size		;Set next step size when we find an intersection and abort
		  (do-dop853-dense
		   #'ode-backreaction-derivatives ode-time y end-time
		   0.0 atol		;only absolute tolerance
		   :dense-step dump-step-size
		   :output-function
		   #'(lambda (step-number this-time y) ;Here for each step
		       ;;The first call is the starting conditions of this run, which we dumped either at the beginning
		       ;;or after we found the intersection that caused us to come back here
		       (when (plusp step-number)
			 ;;Y is the position at the next dump time.  If there is no intersection, we dump it.  If
			 ;;there is an intersection, backreaction-simulate dumps it into nnnn-simulate-... where
			 ;;nnnn is the previous dump step.  Then we dump the post-simulation result.
			 ;;Set variables to position at step
			 (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas) (unflatten-hats y))
			 ;;One oscillation takes L/2.  So age increases by N L/2 and agegmu by NGmu L/2.  The
			 ;;increase in NGmu is just the dump step.
			 (incf agegmu (* dump-step-size (/ length 2)))
			 (let ((backreaction-length (average-total-sigma a-dsigmas b-dsigmas))) ;New length
			   (incf (backreaction-info-length-lost-backreaction *backreaction-info*)
				 (- length backreaction-length))
			   (setq length backreaction-length)
			   (setf (backreaction-info-fraction-lost *backreaction-info*) (- 1.0 (/ length initial-length)))
			   )
			 (multiple-value-bind (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas new-length)
			     (and (not no-intersections)
				  (backreaction-simulate simulate-step dump-step intersection-step
							 a-hats a-dsigmas b-hats b-dsigmas this-time
							 t)) ;always dumping
			   (when new-a-hats		     ;Was there an intersection?  Must start over after dumping
			     ;;Mark this intersection step used (this would only matter if we didn't dump every time)
			     (incf intersection-step)
			     (mirror-images (setq a-hats new-a-hats a-dsigmas new-a-dsigmas)) ;Install new
			     (setq length new-length)
			     (setq intersection t))
			   (incf dump-step)
			   (setq length (average-total-sigma a-dsigmas b-dsigmas))
			   (format t "Dumping step ~D at Ngmu ~S, agegmu ~S: Length ~S lost ~1$%~%"
				   dump-step this-time agegmu length (* 100 (- 1.0 (/ length initial-length))))
			   (dump-backreaction dump-step this-time agegmu
					      a-hats a-dsigmas b-hats b-dsigmas length tolerance maximize-errors)
			   (setq intersection-step 0) ;Reset counter (this would only matter if we didn't dump every time)
			   (when (<= length length-goal) ;After dump, check if finished
			     (return-from outer nil))
			   (incf simulate-step)
			   (incf total-derivative-count *backreaction-derivative-count*)
			   (setq *backreaction-derivative-count* 0) ;Start counting anew for each segment
			   (when intersection
			     (setq ode-time this-time)
			     :abort))))	;Exit dop853 to start over after intersection
		   :maximize-errors maximize-errors :start-step-size step-size)))))
       (let* ((end-time (get-universal-time))
	      (elapsed (- end-time start-time))
	      (wait-fraction (/ *backreaction-job-wait-time* elapsed)))
	 (format t "Backreaction finished with loss fraction ~S at ~A~%" (- 1.0 (/ length initial-length))
		 (sb-int:format-universal-time nil end-time :style :short :print-weekday nil
					       :print-timezone nil))
	 (format t "Elapsed time: ~$ hours, ~$% waiting for jobs ~%" (/ (- end-time start-time) 3600.0)
		 (* wait-fraction 100)))
       (format t "Derivatives were computed ~D times" total-derivative-count)
       (with-open-file (stream (backreaction-success-filename *backreaction-directory*)
			       :direction :output)) ;Create empty file to say success
       ))))	

(defun ode-backreaction-derivatives (time y dydx) ;Store derivatives at position y
  (format t "~D: Computing derivatives at time ~S~%" (incf *backreaction-derivative-count*) time) (force-output)
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
      (unflatten-hats y)
      (mirror-images
       (loop for i from 0 to (1- (length a-dsigmas))
	     for a-dsigma = (aref a-dsigmas i)
	     do (when (< a-dsigma *create-worldsheet-minimum-edge*)
		  (warn "Very small ~S-sigma: length ~S, index ~D" :a a-dsigma i))))
      (multiple-value-bind (dd-a dd-b)
	(process-backreaction-agenda a-hats a-dsigmas b-hats b-dsigmas *backreaction-derivative-count*)
      (fudge-dd-ab a-hats a-dsigmas b-hats b-dsigmas dd-a dd-b) ;fudge errors, go to rest frame
      (mirror-images
       (dotimes (i (length a-hats))
	 (dotimes (component 3)		;Install in flat array
	   (setf (aref dydx (flatten-hats-index :a i component)) (3vector-component (aref dd-a i) component)))))))
  (incf (backreaction-info-derivative-calls *backreaction-info*)))

;;All dd-a and dd-b have been computed.  But because we only compute the backreaction at the center of each
;;segment, these data will not be consistent: they will not leave the loop closed.  So we have to fudge them.
;;In addition, we have to account for the changes due to going back to the rest frame.
;;Loop velocity is v = \sum(da)/L=\sum(db)/L for closure, with L = \sum|da| = \sum|db|.  Initially v = 0.
;;After applying changes we have dv = \sum(da)/L.  While L may change, since \sum(da)=0, it doesn't affect dv.
;;These are applied differently.  The fudging is distributed proportionately to the changes, so unchanging
;;segments won't be fudged.  The boost applies to segments in proportion to their length.
;;
;;3vectors in dd-a and dd-b will be modified
(defun fudge-dd-ab (a-hats a-dsigmas b-hats b-dsigmas dd-a dd-b)
  (mirror-image-let ((total-dd-a-length 0.0)			   ;sum(|dda|)
		     (total-dd-a-a 0.0)				   ;sum(dda . A')
		     (total-length-dd-a-a (make-zero-3vector))	   ;sum(|dda| A')
		     (total-dd-a (make-zero-3vector)))		   ;sum(dda)
    (mirror-images
     (loop for hat across a-hats		    ;A'
	   for change across dd-a		    ;dda
	   for length = (3vector-length change)	    ;|dda|
	   do (incf total-dd-a-length length)
	   do (incf total-dd-a-a (3vector-dot change hat))
	   do (3vector-incf total-length-dd-a-a (3vector-scale hat length))
	   do (3vector-incf total-dd-a change)))
    (let* ((d1 (+ total-dd-a-length total-dd-b-length)) ;See Jeremy's email with Ken's correction
	   (d2 (- total-dd-a-a total-dd-b-b))
	   (d3 (3vector+ total-length-dd-a-a total-length-dd-b-b))
	   (d4 (3vector- total-dd-a total-dd-b))
	   ;;We will adjust each dda -> dda + |dda| (c A' + f), ddb -> ddb - |ddb| (c A' + f)
	   (c (/ (- (* d1 d2) (3vector-dot d3 d4)) (- (3vector-squared-length d3) (expt d1 2))))
	   (f (3vector-scale (3vector+ d4 (3vector-scale d3 c)) (/ -1.0 d1)))
	   (a-sign 1)
	   (b-sign -1))
      (format t "Fudging with c = ~S, f = ~S~%" c f)
      (mirror-images			;Do the fudging
       (dotimes (i (length a-hats))
	 (3vector-incf (aref dd-a i)
		       (3vector-scale (3vector+ (3vector-scale (aref a-hats i) c) f)
				      (* (3vector-length (aref dd-a i)) a-sign))))))
    (check-dds a-hats dd-a b-hats dd-b nil)
    ;;Fudging is done.  Now compute center of mass acceleration
    (let ((dv (3vector-scale (3vector+ (3vector-total dd-a) (3vector-total dd-b)) ;avg change divided by L
			     (/ 1 (+ (total-sigma a-dsigmas) (total-sigma b-dsigmas))))))
      (mirror-images
       (dotimes (i (length a-hats))
	 ;;Infinitesimal Lorentz transformation is just to change da by -v |da|
	 (3vector-decf (aref dd-a i) (3vector-scale dv (aref a-dsigmas i))))))
    (check-dds a-hats dd-a b-hats dd-b t)))

;;Check that these derivatives maintain loop closure in space and time by checking that the sums of dd-a and dd-b
;;are equal and that the sums of dd-a . A' and dd-b . B' are equal.  If rest-frame-p is set, check that the sums of dd-a
;;and dd-b are each zero.
(defun check-dds (a-hats dd-a b-hats dd-b rest-frame-p)
  (mirror-image-let ((total-dd-a (make-zero-3vector))
		     (total-dd-a-a 0.0))
    (mirror-images
     (loop for hat across a-hats	;A'
	   for change across dd-a	;dda
	   do (incf total-dd-a-a (3vector-dot change hat))
	   do (3vector-incf total-dd-a change)))
    (let ((fudge (* (/ (+ (abs total-dd-a-a) (abs total-dd-a-a)) 2) 1e-12))) ;We allow deviations 1e-12 of the average
      (if rest-frame-p
	  (mirror-images (assert (< (3vector-length total-dd-a) fudge))) ;In rest frame, should vanish separately
	(assert (3vector= total-dd-a total-dd-b fudge)))		 ;If not, should be equal
      (assert (fudge= total-dd-a-a total-dd-b-b fudge)))))		 ;In any case, time components should be equal

(defun backreaction-agenda (a-hats b-hats)
  (mirror-image-let ((a-agenda (loop for index below (length a-hats)
				     collect (list :a index))))
    (nconc a-agenda b-agenda)))

;;Compute all segments.  The step is the step number from the ODE solver, and the only purpose is to prevent filename
;;conflicts.
(defun process-backreaction-agenda (a-hats a-dsigmas b-hats b-dsigmas step)
  (let ((agenda (backreaction-agenda a-hats b-hats)))
    (mirror-image-let ((dd-a (make-array (length a-hats))))
      (multiple-value-bind (our-segments job-segments job-time)
	  (backreaction-bunch-numbers (length a-hats) (length b-hats) (length agenda))
	(let* ((start-real-time (get-internal-real-time))
	       (our-agenda (subseq agenda 0 our-segments))
	       (remainder (nthcdr our-segments agenda))
	       (job-agendas (and remainder
				 (loop collect (loop repeat job-segments
						     while remainder ;Stop if we run out
						     collect (pop remainder))
				       while remainder))))
	  (format t "~S to do.  ~:[No batch jobs~2*~;~D job~:P with up to ~D each~]. We do ~S segment~:P~%"
		  (length agenda) job-agendas (length job-agendas) (length (first job-agendas)) (length our-agenda))
	  (when job-agendas		;First submit batch jobs
	    (submit-backreaction-agendas a-hats a-dsigmas b-hats b-dsigmas job-agendas step job-time))
	  (process-our-backreaction-agenda a-hats a-dsigmas b-hats b-dsigmas dd-a dd-b our-agenda) ;Do our own work
	  (when job-agendas
	    (let ((start-time (get-internal-real-time)))
	      (read-dd-files job-agendas step dd-a dd-b)
	      (incf *backreaction-job-wait-time* (/ (- (get-internal-real-time) start-time)
						    internal-time-units-per-second))))
	  (let ((elapsed (/ (- (get-internal-real-time) start-real-time)
			    (double-float internal-time-units-per-second))))
	    (format t "Done in ~S minutes~%" (/ elapsed 60))
	    (incf (backreaction-info-real-time-backreaction *backreaction-info*) elapsed))
	  (values dd-a dd-b))))))

;;Do the work we do ourselves
(defun process-our-backreaction-agenda (a-hats a-dsigmas b-hats b-dsigmas dd-a dd-b agenda)
  (setq agenda (coerce agenda 'vector))
  (thread-dispatch (length agenda) *threads*
    #'(lambda (index)
	(with-local-resources
	  (format t ".") (force-output)
	  (destructuring-bind (type slot) (aref agenda index)
	    (choose-mirror-image type
		 (let ((d-a-prime (a-segment-change a-hats a-dsigmas b-hats b-dsigmas slot)))
		   ;;Change to da = (A' change) dsigma
		   (setf (aref dd-a slot) (3vector-scale d-a-prime (aref a-dsigmas slot))))
		 )))))
  (terpri))

;;Submit batch jobs to do computations
(defun submit-backreaction-agendas (a-hats a-dsigmas b-hats b-dsigmas agendas step job-time)
  (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas
		      (backreaction-data-filename *backreaction-directory* step)) ;Pass loop to batch jobs
  (loop with *batch-time-limit* = (format nil "~D" (ceiling (* job-time 3) ;allow 3 times expected runtime
							    60))		    ;Limit is in minutes
	for agenda in agendas
	for index from 0
	do (do-submit *backreaction-directory*		;Connected directory of batch job
		      (format nil "br-~A-~D" (backreaction-directory-name *backreaction-directory*) index)
		      (format nil "backreaction-~D-" index)
		      (cons "--requeue" worker-batch-flags) ;Requeue if suspended
		      `(compute-backreaction-agenda ',agenda ,step ,*backreaction-filename-digits*)))
  (incf (backreaction-info-batch-jobs *backreaction-info*) (length agendas)))

;;Tell his job has been restarted and give number of restarts
(defun job-restarted-p ()
  (sb-ext:posix-getenv "SLURM_RESTART_COUNT")) ;NIL if not there, or number

;;Run in batch job to do requested computations
(defun compute-backreaction-agenda (agenda step *backreaction-filename-digits*
					   &optional *backreaction-directory*) ;Defaults to connected directory
  (format t "Computing backreaction.  Agenda:~%~S ~%" agenda)
   (unless *backreaction-directory*			;NIL means connected directory
     (setq *backreaction-directory* (sb-posix:getcwd)))	     ;Using "." causes some trouble with batch job submission
   (when (job-restarted-p)				     ;If job restarted, skip items previously computed
     (format t "Job was restarted ~%")
     (loop while agenda
	   for (type index) = (car agenda)
	   while (probe-file (choose-mirror-image type
			       (backreaction-segment-change-filename *backreaction-directory* step :a index)))
	   do (pop agenda))
     (format t "Remaining agenda:~%~S ~%" agenda))
   (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
       (read-hats-dsigmas (backreaction-data-filename *backreaction-directory* step))
     ;;      (format t "Read file ~A~%" (backreaction-data-filename *backreaction-directory* step))
     (format t "Computing: ")
     (loop for (type index) in agenda
	   do (format t "~S ~D " type index) (force-output)
	   do (choose-mirror-image type
		(let* ((output-filename
			(backreaction-segment-change-filename *backreaction-directory* step :a index))
		       (segment-correction (a-segment-change a-hats a-dsigmas b-hats b-dsigmas index))
		       (dd-a (3vector-scale segment-correction (aref a-dsigmas index))))
		  (with-open-file (stream output-filename :direction :output
					  :element-type '(unsigned-byte 64))
		    (write-4vector dd-a stream)
		    ;;		     (format t "Wrote to ~A~%" stream)
		    ))))
     (terpri)))

;;Open file.  If this yields a FILE-DOES-NOT-EXIST error, wait 1 second and try again forever
;;When file is successfully opened, do body.  But if this signals END-OF-FILE, wait 1 second
;;and return the whole process from the beginning.
;;Caveats: The body must not have side effects that cannot be repeated.
;;         Unhandled errors of the above types will come here, even if they involve other files
;;         opened inside the body.
(defmacro with-open-file-wait-retry (binding &body body)
  `(loop
    (handler-case
	(return (with-open-file ,binding ,@body))
      ((or file-does-not-exist end-of-file) ()
       (sleep 1)))))

(defun read-dd-files (agendas step dd-a dd-b)
  (format t "Reading segment changes:~%")
  (dolist (agenda agendas)
    (destructuring-bind (type1 index1) (first agenda)
      (destructuring-bind (type2 index2) (car (last agenda))
	(if (eq type1 type2) (format t "~A ~D..~D" type1 index1 index2) ;All the same
	  (loop for last-index = nil then index				;Find last index of first type
		for (type index) in agenda
		while (eq type type1)
		finally (format t "~@{~A ~D..~D~^ ~}" type1 index1 last-index type2 0 index2)))))
    (force-output)
    (loop for (type index) in agenda
	  for file = (backreaction-segment-change-filename *backreaction-directory* step type index)
	  do (with-open-file-wait-retry (stream file :element-type '(unsigned-byte 64))
	       (choose-mirror-image type
		 (setf (aref dd-a index) (read-4vector stream))))
	  do (format t ".") (force-output))
    ;;Wait until the whole job has been read before deleting previous files, so a restarted job knows what to do
    (loop for (type index) in agenda
	  do (delete-file (backreaction-segment-change-filename *backreaction-directory* step type index)))
    (terpri))
  (delete-file (backreaction-data-filename *backreaction-directory* step)))

(defun dump-backreaction (step time agegmu a-hats a-dsigmas b-hats b-dsigmas length tolerance maximize-errors)
  (let ((now (get-internal-real-time)))
    (setf (backreaction-info-dump-step *backreaction-info*) step
	  (backreaction-info-time *backreaction-info*) time
	  (backreaction-info-agegmu *backreaction-info*) agegmu
	  (backreaction-info-length *backreaction-info*) length
	  (backreaction-info-a-kinks *backreaction-info*) (length a-hats)
	  (backreaction-info-b-kinks *backreaction-info*) (length b-hats)
	  (backreaction-info-real-time-total *backreaction-info*) (/ (- now *backreaction-dump-real-time*)
								     (double-float internal-time-units-per-second))
	  (backreaction-info-tolerance *backreaction-info*) tolerance
	  (backreaction-info-maximize-errors *backreaction-info*) maximize-errors)
    (setq *backreaction-dump-real-time* now))
  (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas
		      (backreaction-dump-filename *backreaction-directory* step))
  (with-open-file (stream (backreaction-info-filename *backreaction-directory* step) :direction :output)
    (prin1 *backreaction-info* stream)
    (terpri stream))
  (setq *backreaction-info* (make-backreaction-info))) ;New structure to accumulate new data

;;Do simulation. DUMP-STEP is the step of the last dump.  INTERSECTION-STEP is the steps since then.
(defun backreaction-simulate (simulate-step dump-step intersection-step
					    a-hats a-dsigmas b-hats b-dsigmas time will-dump)
  (format t "Simulating at time ~A~%" time)
  (let* ((start-real-time (get-internal-real-time))
	 (start-length (average-total-sigma a-dsigmas b-dsigmas))
	 ;;Simulation volume boundaries start at size 0 and retreat at the speed of light.  So if we start at time
	 ;;length/2, the loop will fit in the simulation volume.  (This could be smaller if we centered it.)
	 (start-time (/ start-length 2))
	 ;;Now we need to make the volume big enough so that we won't run into any future boundaries in any reasonable
	 ;;number of oscillations.
	 (duration (+ (* 2 start-time) (* start-length 5))) ;corner to start + end to final corner + 10 oscillations
	 (size (/ duration (job-duration 1.0)))) ;How much larger must it be to get desired duration?
    (setq *intersections-performed* 0 *rejoinings-performed* 0
	  *job-number* 0)		;Must be set outside simulate to avoid problems with evolve-until
    (simulate
     #'(lambda () (create-worldsheet-at-time a-hats a-dsigmas b-hats b-dsigmas start-time :install t :tags t))
     :count-rejoinings t :output-directory (backreaction-simulate-directory *backreaction-directory* simulate-step)
     :size size
     :start start-time :end (+ start-time (/ start-length 2)))	;One oscillation
    (when (plusp *intersections-performed*)	 ;If no intersections, done
      (loop do (evolve-until (+ (current-time) (/ start-length 10))) ;If intersections, it may take longer.
	    until (non-selfintersect-p)))		       ;Evolve in chunks until we find non-self-intersecting
    ;;Now delete directory
    (do-run-program "rm" :args (list "-r" (backreaction-simulate-directory *backreaction-directory* simulate-step)))
    (multiple-value-prog1
	(when (plusp *intersections-performed*)	;Return NIL if no intersections
	  (incf (backreaction-info-simulations-with-intersections *backreaction-info*))
	  (incf (backreaction-info-intersections *backreaction-info*) *intersections-performed*)
	  (incf (backreaction-info-rejoinings *backreaction-info*) *rejoinings-performed*)
	  (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas ;Write situation before simulation
			      (backreaction-pre-simulate-filename *backreaction-directory* dump-step intersection-step))
	  (multiple-value-bind (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas)
	      (get-ab-data-dsigmas (longest-loop))
	    (multiple-value-bind (rest-a-hats rest-a-dsigmas rest-b-hats rest-b-dsigmas)
		;;Make sure everything is closed and in the rest frame, though this seems always to be the case already
		(rest-frame-close-hats new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas)
	      (unless will-dump		;Unless about to do it as a dump
		(write-hats-dsigmas rest-a-hats rest-a-dsigmas rest-b-hats rest-b-dsigmas ;Write state after simulation
				    (backreaction-post-simulate-filename *backreaction-directory* dump-step intersection-step))
		(with-open-file (stream (backreaction-simulate-info-filename *backreaction-directory* dump-step intersection-step)
					:direction :output)
		  (prin1 *backreaction-info* stream)
		  (terpri stream)))
	      (let* ((new-length (average-total-sigma new-a-dsigmas new-b-dsigmas))
		     (rest-length (average-total-sigma rest-a-dsigmas rest-b-dsigmas))
		     (new-lost  (- start-length new-length))
		     (rest-lost (- new-length rest-length)))
		(format t "~&~D intersection~:P, ~D rejoining~:P.  New length ~S, ~S lost in sim, ~S to rest.~%"
			*intersections-performed* *rejoinings-performed* rest-length new-lost rest-lost)
		(incf (backreaction-info-length-lost-simulation *backreaction-info*) new-lost)
		(incf (backreaction-info-length-lost-simulation-rest-frame *backreaction-info*) rest-lost)
		(values rest-a-hats rest-a-dsigmas rest-b-hats rest-b-dsigmas rest-length))
	      )))
      (incf (backreaction-info-real-time-simulation *backreaction-info*)
	    (/ (- (get-internal-real-time) start-real-time) (double-float internal-time-units-per-second))))))

;; Check if the loop is non-selfintersecting by checking that it has diamonds 
;; with the countup > 1
(defun non-selfintersect-p ()
  (let ((start (longest-loop)))
    (loop for d = (diamond-e start) then (diamond-e d)
	  thereis (> (diamond-countup d) 1)
	  until (eq d start))))

(defun find-backreaction-digits (directory) ;See how many in use in filenames
  (loop for digits from 1 below 10
	when (let ((*backreaction-filename-digits* digits)) ;See if file is there
	       (probe-file (backreaction-dump-filename directory 0)))
	return digits
	finally (error "~A does not seem to contain any backreaction datafiles" directory)))
 
;;Summarize what happened with a loop, including when there was an ejection of the loop that produced
;;significant recoil velocity.  "Ejection velocity" here means the velocity of the main loop after the small
;;loop was ejected.
(defun total-backreaction-info (directory &key (major-intersection-threshold 0.05)
					  (ejection-velocity 2e-3) all-ejections)
  (setq directory (merge-pathnames directory batch-root-directory))
  (loop with *backreaction-filename-digits* = (find-backreaction-digits directory)
	with initial-length = nil
	with last-length = nil
	with last-a-kinks = nil
	with last-b-kinks = nil
	with max-v = nil
	with max-v-step
	with ejected-with-loss = nil	;Amount of loss at time of ejection
	with major-loss = nil
	with ending-time = 0.0
	initially (multiple-value-bind (a-hats a-dsigmas b-hats) (read-hats-dsigmas (backreaction-dump-filename directory 0))
		    (declare (ignore a-dsigmas))
		    (format t "N_A = ~D, N_B = ~D~%" (length a-hats) (length b-hats)))
	for step from 0
	for info =  (with-open-file (stream (backreaction-info-filename directory step) :if-does-not-exist nil)
		      (and stream (read stream)))
	while info
	if initial-length
	do (let ((lost (/ (backreaction-info-length-lost-simulation info) initial-length)))
	       (when (> lost major-intersection-threshold)
		 (format t "Major intersection in step ~D: ~S of initial length lost in simulation~%" step lost)
		 (when (or (null major-loss) (> lost major-loss))
		   (setq major-loss lost))))
	else do (setq initial-length (backreaction-info-length info)) ;First time
	;;we check the rest-frame length lost because an intersection--rejoining event can cause this value to be
	;;negative due to numerical error---the correct answer would be 0.0
	when (plusp (backreaction-info-length-lost-simulation-rest-frame info)) ;Intersected?
	;;Boost factor that it had before going to rest frame is (energy moving)/(rest energy) =
	;;1+dE/L, where dE is the energy lost going to the rest frame and L the rest length 
	do (let* ((gamma (+ 1.0 (/ (backreaction-info-length-lost-simulation-rest-frame info)
				   (backreaction-info-length info))))
		  (v (sqrt (- 1.0 (expt gamma -2))))
		  (loss (/ (- initial-length last-length) ;fractional loss before this intersection
			   initial-length)))
	     (when (or (null max-v) (> v max-v))
	       (setq max-v v max-v-step step))
	     (when  (> v ejection-velocity)
	       (when (or all-ejections (null ejected-with-loss))
		 ;;Intersection produces two kinks in each of A and B.  The number in the ejected loop is thus what
		 ;;we had before plus 2, minus what we have afterward.
		 (mirror-image-let ((loop-a-kinks (- (+ last-a-kinks 2) (backreaction-info-a-kinks info))))
		   (format t "~:[First e~;E~]jection at step ~D (~4$% evaporation), velocity ~D in ~D As and ~D Bs~%"
			   all-ejections step (* loss 100) v loop-a-kinks loop-b-kinks))
		 (unless ejected-with-loss
		   (setq ejected-with-loss loss)))))
	sum (backreaction-info-intersections info) into intersections
	sum (backreaction-info-simulations-with-intersections info) into simulations-with-intersections
	sum (backreaction-info-rejoinings info) into rejoinings
	sum (backreaction-info-length-lost-backreaction	info) into length-lost-backreaction
	;; sum (backreaction-info-length-lost-rest-frame info) into length-lost-rest-frame
	sum (backreaction-info-length-lost-simulation-rest-frame info) into length-lost-simulation-rest-frame
	sum (backreaction-info-length-lost-simulation info) into length-lost-simulation
	sum (backreaction-info-real-time-backreaction info) into real-time-backreaction
	sum (backreaction-info-real-time-simulation info) into real-time-simulation
	sum (backreaction-info-real-time-total info) into real-time-total
	sum (backreaction-info-derivative-calls info) into derivative-calls
	do (setq last-length (backreaction-info-length info)
		 ending-time (backreaction-info-time info))
	do (mirror-images (setq last-a-kinks (backreaction-info-a-kinks info)))
	finally
	(when max-v
	  (format t "Fastest speed after simulation ~S at step ~D~%" max-v max-v-step))
	(format t "~D intersection~:P (~D rejoining~:P) in ~D batches.~%"
		intersections rejoinings simulations-with-intersections)
	(format t "Evolved to NGmu ~F~%" ending-time)
	(format t "Initial length ~S.  Fractions lost:~%  backreaction: ~S~%  simulation: ~S~%  rest frame after simulation: ~S~%"
		initial-length
		(/ length-lost-backreaction initial-length)
		;; (/ length-lost-rest-frame initial-length)
		(/ length-lost-simulation initial-length)
		(/ length-lost-simulation-rest-frame initial-length))
	(format t "Total real time to compute: ~A.  Backreaction: ~A, simulation: ~A, other: ~A.~%"
		(format-seconds (floor real-time-total))
		(format-seconds (floor real-time-backreaction))
		(format-seconds (floor real-time-simulation))
		(format-seconds (floor (- real-time-total real-time-simulation real-time-backreaction)))
		)
	(format t "Derivative calls: ~D~%" derivative-calls)
	(return (values initial-length intersections major-loss max-v ejected-with-loss))))


;;List of pre- and post-simulation files and next dump
(defun backreaction-simulation-file-list (directory step)
  (append
   (loop for intersection-step from 0
	 for pre = (backreaction-pre-simulate-filename directory step intersection-step)
	 for post = (backreaction-post-simulate-filename directory step intersection-step)
	 while (probe-file pre)		;An(other) intersection this interval?
	 collect pre
	 when (probe-file post) collect post)
   (let ((dump (backreaction-dump-filename directory (1+ step))))
     (and (probe-file dump) (list dump))) ;Return NIL if ran out of files
   ))
		  
;;Call function with step, hats and sigmas, and rotation matrix for this step and for total,
;;giving transformation of shape from initial conditions to these.
(defun scan-backreaction-dumps (function directory &optional (start 0) end (step 1))
  (loop with *backreaction-filename-digits* = (find-backreaction-digits directory)
	with total-r = (make-identity-float 3)
	with old-a-hats = nil and old-b-hats
	with a-hats and a-dsigmas and b-hats and b-dsigmas
	for this-step from start by step
	until (and end (>= (1+ this-step) end)) ;When this-step = end-2 we are actually dumping end-1
	initially
	(multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas) ;Read initial state
	  (read-hats-dsigmas (backreaction-dump-filename directory start)))
	(mirror-images (setq old-a-hats a-hats))		 ;Set up for first rotation
	(funcall function 0 a-hats a-dsigmas b-hats b-dsigmas total-r total-r) ;Call with initial state
	for files = (backreaction-simulation-file-list directory this-step)
	while files			;Exit if we have processed all the files
	do
	(format t "~D " this-step) (force-output)
	;;Go through files and accumulate rotation.  When done, hats are in both a/b-hats and old-a/b-hats
	(loop with this-r = (make-identity-float 3)
	      for file in files
	      for backreaction = t then (not backreaction) ;Alternate steps are backreaction and simulation
	      do 
	      (multiple-value-setq (a-hats a-dsigmas b-hats b-dsigmas)
		(read-hats-dsigmas file)) ;Read next position
	      (when backreaction	  ;A backreaction step?
		;;Accumulate total rotation to transform position into new
		(let ((one-r (find-backreaction-rotation old-a-hats a-dsigmas a-hats
							 old-b-hats b-dsigmas b-hats)))
		  (if one-r (setq this-r (dotmm one-r this-r))
		    (warn "Failed to find rotation at step ~D" this-step))))
	      (mirror-images (setq old-a-hats a-hats))
	      finally
	      (setq total-r (dotmm this-r total-r))
	      (funcall function (+ this-step step) ;It is really the next step
		       a-hats a-dsigmas b-hats b-dsigmas this-r total-r))))

(defun describe-backreaction-rotation (directory &optional (start 0) end)
  (let ((first-axis nil)
	first-angle
	(total-rotation 0))
    (scan-backreaction-dumps
     #'(lambda (step a-hats a-dsigmas b-hats b-dsigmas this-r total-r)
	 (declare (ignore a-hats a-dsigmas b-hats b-dsigmas total-r))
	 (unless (= step start)	;First callback has no rotation
	   (multiple-value-bind (angle axis)
	       (rotation-matrix-info this-r)
	     (unless first-axis		;Second callback, first rotation?
	       (setq first-axis axis first-angle angle))
	     (incf total-rotation (abs angle))
	     (format t ": ~S (diff ~$%) around ~S (diff ~S)~%"
		     angle 
		     (* 100 (/ (- (abs angle) (abs first-angle))
			       (abs first-angle)))
		     axis
		     (let ((difference-angle (spherical-angle axis first-axis)))
		       (min difference-angle (- pi difference-angle)))))))
     directory start end)
    total-rotation))

;;Make movie of a' and b' points moving under backreaction, with overall rotation removed
(defun plot-ab-movie (directory &rest gnuplot-keys
				&key (start 0) end (step 1)
				(output-directory directory) 
				(style :points)	;Style for a' and b'
				interpolate ;Draw smooth curves with this number of points
				(frame-rate 30)
				no-rotate
				circles	;Show segments sizes with circles
				mark-kinks
				show-starts ;Circle starting points of a' and b'
				&allow-other-keys)
  (unless (probe-file directory)
    (error "~A does not exist" directory))
  (with-group-write-access
   (ensure-directories-exist (format nil "~A/" output-directory))
   (setq output-directory (namestring (truename output-directory))) ;Get real name, not relative or ~
   (let ((frames 0))
     (scan-backreaction-dumps
      #'(lambda (this-step a-hats a-dsigmas b-hats b-dsigmas this-r total-r)
	  (declare (ignore this-r))
	  (when no-rotate (setq total-r (make-identity-float 3)))
	  (mirror-image-let ((output-file (backreaction-filename output-directory this-step "plot-ab.png"))
			     (smooth-a-hats (if interpolate (interpolate-hats a-hats :split-count interpolate)
					      a-hats)))
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
		   :axes (list (dotmv33 total-r (make-3vector 0.0 0.0 1.0)) ;rotated original Z
			       (dotmv33 total-r (make-3vector 1.0 0.0 0.0))) ;rotated original X
		   :show-axes t
		   :mark-locations
		   (append (and mark-kinks
				(append (mark-large-kink-locations a-hats mark-kinks)
					(mark-large-kink-locations b-hats mark-kinks)))
			   (and show-starts
				(list (list (format nil "set object circle radius 0.1") (aref a-hats 0))
				      (list (format nil "set object circle radius 0.1") (aref b-hats 0)))))
		   :png output-file
		   :title (format nil "Step ~D" this-step)
		   gnuplot-keys)
	    (incf frames)))
      directory start end step)
     (ffmpeg-movie output-directory "plot-ab.png" (find-backreaction-digits directory)
		 start step frames  (format nil "~A/plot-ab.mp4" output-directory)
		 :frame-rate frame-rate))))

;;Make movie of a' and b' points moving under backreaction, with overall rotation removed
;;Specially mark a range of b hats which are given as inputs to study during debugging of the reduced smoothing in b
(defun plot-ab-movie-interesting-b (directory first-b-hat last-b-hat &rest gnuplot-keys
					      &key (start 0) end (step 1)
					      (output-directory directory)
					      (style :points)	;Style for a' and b'
					      interpolate ;Draw smooth curves with this number of points
					      (frame-rate 30)
					      no-rotate
					      circles	;Show segments sizes with circles
					      &allow-other-keys)
  (with-group-write-access
   (ensure-directories-exist (format nil "~A/" output-directory))
   (let ((output-files nil)
	 (frames 0))
     (scan-backreaction-dumps
      #'(lambda (step a-hats a-dsigmas b-hats b-dsigmas this-r total-r)
	  (declare (ignore this-r))
	  (when no-rotate (setq total-r (make-identity-float 3)))
					;currently coded for the specific interesting b hats we were using
					;may need to be updated to be more universal	  
	  (let* ((intrst-b-hats (coerce (coerce (subseq b-hats first-b-hat last-b-hat) 'list) 'vector))
		 (intrst-b-dsigmas (coerce (coerce (subseq b-dsigmas first-b-hat last-b-hat) 'list) 'vector)))
	    (mirror-image-let ((output-file (backreaction-filename output-directory step "plot-ab.png"))
			       (smooth-a-hats (if interpolate (interpolate-hats a-hats :split-count interpolate)
						a-hats)))
	      (let* ((smooth-interesting-b-hats (if interpolate (interpolate-hats intrst-b-hats :split-count interpolate)      
						  intrst-b-hats)))
		(apply #'plot-tangent-vectors
		       (append (list (list a-hats "A" style)
				     (list b-hats "B" style)
				     (list intrst-b-hats "B*" style))
			       (and interpolate
				    (list (list smooth-a-hats nil "lines lt 1")
					  (list smooth-b-hats nil "lines lt 2")
					  (list smooth-interesting-b-hats nil "lines lt 3")))
			       (and circles
				    (list (list (hats-circle-data a-hats a-dsigmas) nil "circles lt 1")
					  (list (hats-circle-data b-hats b-dsigmas) nil "circles lt 2")
					  (list (hats-circle-data intrst-b-hats intrst-b-dsigmas) nil "circles lt 3"))))
		       :loop t
		       :axes (list (dotmv33 total-r (make-3vector 0.0 0.0 1.0)) ;rotated original Z
				   (dotmv33 total-r (make-3vector 1.0 0.0 0.0))) ;rotated original X
		       :show-axes t
		       :png output-file
		       :title (format nil "Step ~D" step)
		       gnuplot-keys)
		(push output-file output-files))
	      (incf frames))))
      directory start end step)
     (ffmpeg-movie output-directory "plot-ab-intrst-b.png" (find-backreaction-digits directory)
		 start step frames  (format nil "~A/plot-ab-intrst-b.mp4" output-directory)
		 :frame-rate frame-rate))))

;using strongest-cusp-per-timestamp creates a list of the strongest cusps to be referenced by
;plot-ab-movie-strong-cusp
(defun strongest-cusps (directory)
  (let* ((output nil)
	 (strongest-cusps (strongest-cusp-per-timestamp directory)))
    (loop for timestamp in strongest-cusps
	  do (push (list (nth 0 (nth 1 timestamp)) (nth 2 (nth 1 timestamp)) (nth 12 (nth 1 timestamp))) output))
					;strongest cusp a-hat, strongest cusp b-hat, intersection (location)
    (setf output (reverse output))
    (let ((my-array (make-array 4 :element-type 'double-float :initial-element 0.0)))
      (push (list 0 0 my-array) output))
    (coerce output 'vector)))

;;Make movie of a' and b' points moving under backreaction, with overall rotation removed
;;Specifically marks the location of the strongest cusp at each timestamp
(defun plot-ab-movie-strong-cusps (directory &rest gnuplot-keys
					     &key (start 0) end (step 1)
					     (output-directory directory)
					     (style :linespoints)
					     interpolate
					     (frame-rate 30)
					     no-rotate
					     circles
					     &allow-other-keys)
  (let* ((strongest-cusps (strongest-cusps directory)))
    (with-group-write-access
     (ensure-directories-exist (format nil "~A/" output-directory))
     (let ((output-files nil)
	   (frames 0))
       (scan-backreaction-dumps
	#'(lambda (step a-hats a-dsigmas b-hats b-dsigmas this-r total-r)
	    (declare (ignore this-r))
	    (when no-rotate (setq total-r (make-identity-float 3)))
	    (let* ((intrst-a-low-bound (nth 0 (aref strongest-cusps step)))
		   (intrst-a-high-bound (cond ((< (+ (nth 0 (aref strongest-cusps step)) 2)
						  (length a-hats))
					       (+ (nth 0 (aref strongest-cusps step)) 2))
					      (t
					       (length a-hats))))
		   (intrst-a-hats (coerce (coerce (subseq a-hats intrst-a-low-bound intrst-a-high-bound)
						  'list)
					  'vector))
		   (intrst-a-dsigmas (coerce (coerce (subseq a-dsigmas intrst-a-low-bound
							     intrst-a-high-bound)
						     'list)
					     'vector))
		   (intrst-b-low-bound (nth 1 (aref strongest-cusps step)))
		   (intrst-b-high-bound (cond ((< (+ (nth 1 (aref strongest-cusps step)) 2)
						  (length b-hats))
					       (+ (nth 1 (aref strongest-cusps step)) 2))
					      (t
					       (length b-hats))))
		   (intrst-b-hats (coerce (coerce (subseq b-hats intrst-b-low-bound intrst-b-high-bound)
						  'list)
					  'vector))
		   (intrst-b-dsigmas (coerce (coerce (subseq b-dsigmas intrst-b-low-bound
							     intrst-b-high-bound)
						     'list)
					     'vector)))
	      (mirror-image-let ((output-file
				  (backreaction-filename output-directory step "plot-ab-strong-cusps.png"))
				 (smooth-a-hats (if interpolate (interpolate-hats a-hats :split-count interpolate)
						  a-hats))
				 (smooth-interesting-a-hats (if interpolate (interpolate-hats intrst-a-hats :split-count
											 interpolate)
							      intrst-a-hats)))
		(apply #'plot-tangent-vectors
		       (cond (interpolate
			      (append
			       (list (list smooth-a-hats nil "lines lt 1")
				     (list smooth-b-hats nil "lines lt 2")
				     (list smooth-interesting-a-hats nil "lines lt 3")
				     (list smooth-interesting-b-hats nil "lines lt 4"))
			       (and circles
				    (list (list (hats-circle-data a-hats a-dsigmas) nil "circles lt 1")
					  (list (hats-circle-data b-hats b-dsigmas) nil "circles lt 2")
					  (list (hats-circle-data intrst-a-hats intrst-a-dsigmas) nil "circles lt 3")
					  (list (hats-circle-data intrst-b-hats intrst-b-dsigmas)
						nil "circles lt 4")))))
			     (t
			      (append (list (list a-hats "A" style)
					    (list b-hats "B" style)
					    (list intrst-a-hats "A*" style)
					    (list intrst-b-hats "B*" style))
				      (and circles
					   (list (list (hats-circle-data a-hats a-dsigmas) nil "circles lt 1")
						 (list (hats-circle-data b-hats b-dsigmas) nil "circles lt 2")
						 (list (hats-circle-data intrst-a-hats intrst-a-dsigmas) nil "circles lt 3")
						 (list (hats-circle-data intrst-b-hats intrst-b-dsigmas)
						       nil "circles lt 4"))))))
		       ;:no-gnuplot t
		       :mark-locations (list (list "set object circle radius 0.05"
					     (nth 2 (aref strongest-cusps step))))
		       :loop t
		       :axes (list (dotmv33 total-r (make-3vector 0.0 0.0 1.0))
				   (dotmv33 total-r (make-3vector 1.0 0.0 0.0)))
		       :show-axes t
		       :png output-file
		       :title (format nil "Step ~D" step)
		       gnuplot-keys)
		(push output-file output-files)
		(incf frames))))
	directory start end step)
       (ffmpeg-movie output-directory "plot-ab-strong-cusps.png" (find-backreaction-digits directory)
		 start step frames  (format nil "~A/plot-ab-strong-cusps.mp4" output-directory)
		 :frame-rate frame-rate)))))


;;Convert a set of image files into a video.  Sadly I don't know how do this simply in general, so we use a kluge
(defun ffmpeg-files (files output-file &key (frame-rate 30))
  (let ((type (pathname-type (first files)))) ;Must be uniform
    (loop for file in files
	  for index from 0
	  do (do-run-program "ln" :args (list "-sf" (native-namestring (pathname file)) ;resolve ~
					      (format nil "/tmp/ffmpeg-frame-~D.~A" index type))))
    (do-run-program "ffmpeg" :args (list "-y"
					 "-r" (format nil "~D" frame-rate)
					 `"-i" (format nil "/tmp/ffmpeg-frame-%d.~A" type)
					 (native-namestring (pathname output-file))))
    (do-run-program "csh" :args (list "-c" (format nil "rm /tmp/ffmpeg-frame-*.~A" type)))))


;;Return commands for marking location of large kinks with their angle
(defun mark-large-kink-locations (hats threshold)
  (loop for index from (1- (length hats)) downto 0
	    for hat = (aref hats index)
	    for next-hat = (aref hats (next-array-index-wrapping index hats))
	    for angle = (spherical-angle hat next-hat)
	    when (> angle threshold)
	    collect (list (format nil "set label '~5$'" angle)
			  (spherical-interpolation hat next-hat 0.5))))

;;Convert a set of image files into a video.  The files are located in the given DIRECTORY and have names
;;with the frame number in DIGITS digits, a hypen, and the FILENAME.
(defun ffmpeg-movie (directory filename digits start step frames output-file &key (frame-rate 30))
  (if (= step 1)			;ffmpeg can't handle steps
      (do-run-program "ffmpeg" :args (list "-y"
					   "-r" (format nil "~D" frame-rate)
					   "-start_number" (format nil "~D" start)
					   "-i" (format nil "~A/%~Dd-~A" directory digits filename)
					   "-frames:v" (format nil "~D" frames)
					   (format nil "~A" output-file)) ;Convert to string
		      :input nil)	   ;Avoid unintended interruptions of ffmpeg
    ;;Can't do it by ffmpeg.  Use links.
    (let ((type (pathname-type filename))) ;Must be uniform
      (loop for this-step from start by step repeat frames
	    for index from 0
	    do (do-run-program "ln" :args (list "-sf" (format nil "~A/~V,'0D-~A" directory digits this-step filename)
						(format nil "/tmp/ffmpeg-frame-~D.~A" index type))))
      (do-run-program "ffmpeg" :args (list "-y"
					   "-r" (format nil "~D" frame-rate)
					   `"-i" (format nil "/tmp/ffmpeg-frame-%d.~A" type)
					   output-file)
		      :input nil)
      (do-run-program "csh" :args (list "-c" (format nil "rm /tmp/ffmpeg-frame-*.~A" type))))))
  

;;Compare different runs for the same loop
(defun plot-ab-movie-compare (directory1 directory2 &rest gnuplot-keys
				&key (start1 0) end1 (step1 1) (start2 0) end2 (step2 1)
				(output-directory directory1) ;Location for png files	
				(output-file "")	      ;mp4 file to write.
				(frame-rate 30)
				&allow-other-keys)
  (setq output-file			;Default to plot-ab-compare.mp4 in directory1
	(merge-pathnames output-file (format nil "~A/plot-ab-compare.mp4" output-directory)))
  (format t "Output file ~A~%" output-file)
  (ensure-directories-exist output-file)
  (loop with digits1 = (find-backreaction-digits directory1)
	with digits2 = (find-backreaction-digits directory2)
	for frames from 0
	for frame1 from start1 by step1
	for frame2 from start2 by step2
	until (or (and end1 (>= frame1 end1)) (and end2 (>= frame2 end2)))
	for file1 = (let ((*backreaction-filename-digits* digits1))
		      (backreaction-dump-filename directory1 frame1))
	for file2 = (let ((*backreaction-filename-digits* digits2))
		      (backreaction-dump-filename directory2 frame2))
	while (and (probe-file file1) (probe-file file2))
	do (multiple-value-bind (a-hats-1 a-dsigmas-1 b-hats-1 b-dsigmas-1)
	       (read-hats-dsigmas file1)
	     (declare (ignore a-dsigmas-1 b-dsigmas-1))
	     (multiple-value-bind (a-hats-2 a-dsigmas-2 b-hats-2 b-dsigmas-2)
		 (read-hats-dsigmas file2)
	       (declare (ignore a-dsigmas-2 b-dsigmas-2))
	       (let ((timestring
		      (format nil "~&Time ~6$"
			      (with-open-file (stream (let ((*backreaction-filename-digits* digits1))
							(backreaction-info-filename directory1 frame1)))
				(backreaction-info-time (read stream))))))
		 (format t "~&~A.~%" timestring)
		 (apply #'plot-tangent-vectors
			(list (list a-hats-1 (format nil "A ~A" directory1) "points")
			      (list b-hats-1 (format nil "B ~A" directory1) "points")
			      (list a-hats-2 (format nil "A ~A" directory2) "points pt 4") ;open squares to compare
			      (list b-hats-2 (format nil "B ~A" directory2) "points pt 17"))
			:axes nil
			:png (let ((*backreaction-filename-digits* digits1))
			       (backreaction-filename output-directory frame1 "plot-ab-compare.png"))
			:title timestring
			:key "above font ',8'"
			gnuplot-keys))))
	finally
	(ffmpeg-movie (truename output-directory) "plot-ab-compare.png" (find-backreaction-digits directory1)
		      start1 step1 frames output-file
		      :frame-rate frame-rate)))

(defun backreaction-Jacobian (a-hats a-dsigmas b-hats b-dsigmas &key (epsilon 1e-8))
  (mirror-image-let* ((n-a (length a-hats))
		      (n (nomirror (+ n-a n-b)))
		      (a-segment-changes (make-array n-a))
		      (Jacobian (make-array (list n n)))
		      ;;J_ij gives the effect of changes to hat j on the segment change of segment i
		      (j 0)
		      new-a-hats	;Used below for data with one slot modified
		      new-a-dsigmas)
    (mirror-images
     (dotimes (index n-a)	;Compute unmodified segment changes
       (setf (aref a-segment-changes index) (a-segment-change a-hats a-dsigmas b-hats b-dsigmas index))
       (format t ".") (force-output)))
    (mirror-images
     (loop for index below n-a
	   for d-a-prime = (4vector-scale (aref a-segment-changes index) ;Ignore original length.  Set to epsilon
					  (/ epsilon (4vector-Euclidian-length (aref a-segment-changes index))))
	   do
	   ;;First copy put both a and b data into new variables.  We will then update ours and leave the other alone.
	   (nomirror (mirror-images (setq new-a-hats a-hats new-a-dsigmas a-dsigmas)))
	   (setq new-a-hats (copy-seq a-hats) new-a-dsigmas (copy-seq a-dsigmas)) ;We will modify this
	   (setf (aref new-a-hats index) (3vector-normalize ;Update just the one slot
					      (3vector+ (aref a-hats index) (4to3vector d-a-prime)))
		 (aref new-a-dsigmas index) (* (aref a-dsigmas index) (1+ (4vector-t d-a-prime))))
	   (multiple-value-setq (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas)
	     (close-4vectors new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas))
	   (backreaction-fill-Jacobian
	    (nomirror new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas a-segment-changes b-segment-changes)
	    Jacobian j epsilon)
	   (incf j)))
    Jacobian))

;;Return the rms average of the 4vector euclidian change in the given segment change when the hats and sigmas
;;are perturbed by TRIES random fractional changes of order epsilon
(defun backreaction-sensitivity (a-hats a-dsigmas b-hats b-dsigmas ab index &key (tries 10) (epsilon 1e-16))
  (mirror-image-let* ((n-a (length a-hats))
		      ;;After we perturbed, we normalize.  Repeated normalization may lead to oscillations,
		      ;;so to be sure the effect is from the perturbation and not from the normalization,
		      ;;we normalize here
		      (normalized-a-hats (map 'vector #'(lambda (x) (3vector-normalize x)) a-hats))
		      normalized-a-dsigmas
		      (new-a-hats (copy-seq a-hats))
		      (new-a-dsigmas (copy-seq a-dsigmas)))
    (multiple-value-setq (normalized-a-hats normalized-a-dsigmas normalized-b-hats normalized-b-dsigmas) ;Must reclose
      (rest-frame-close-hats normalized-a-hats a-dsigmas normalized-b-hats b-dsigmas))
    (loop with original = (choose-mirror-image ab ;Result of normalizing but not otherwise perturbing
			    (a-segment-change normalized-a-hats normalized-a-dsigmas normalized-b-hats
					      normalized-b-dsigmas index))
	  repeat tries
	  with change
	  do (mirror-images
	      (dotimes (i-a n-a)
		(setf (aref new-a-hats i-a)
		      (3vector-normalize (3vector+ (aref a-hats i-a) ;Perturb the original hat
						   (make-3vector (- (random epsilon) (/ epsilon 2))
								 (- (random epsilon) (/ epsilon 2))
								 (- (random epsilon) (/ epsilon 2))))))
		(setf (aref new-a-dsigmas i-a) (+ (aref a-dsigmas i-a)
						    (* epsilon (aref a-dsigmas i-a) (- (random 1.0) 0.5))))))
	  do (multiple-value-setq (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas) ;Must reclose
	       (rest-frame-close-hats new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas))
	  do (setq change (4vector-euclidian-distance
			   (choose-mirror-image ab (a-segment-change new-a-hats new-a-dsigmas
								     new-b-hats new-b-dsigmas  index))
			   original))
	  sum (expt change 2) into total
	  finally (return (sqrt (/ total tries))))))

;;Fill in slot j of Jacobian by computing segment changes with new data and comparing with old segment changes
(defun backreaction-fill-Jacobian (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas
					      a-segment-changes b-segment-changes Jacobian j epsilon)
  (format t "~&~D" j) (force-output)
  (let ((i 0))				;Affected slot
    (mirror-images
     (dotimes (index (length a-segment-changes))
       (setf (aref Jacobian i j)
	     (4vector-scale (4vector-
			     (a-segment-change new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas index)
			     (aref a-segment-changes index))
			    (/ 1 epsilon)))
       (format t ".") (force-output)
       (incf i)))))
  
;;Make sure that Jacobian properly predicts changes to segment changes
(defun test-backreaction-Jacobian (a-hats a-dsigmas b-hats b-dsigmas NGmu)
  (let ((Jacobian (backreaction-Jacobian a-hats a-dsigmas b-hats b-dsigmas)))
    (mirror-image-let* ((n-a (length a-hats))
			(n (nomirror (+ n-a n-b)))
			(a-segment-changes (make-array n-a))
			;;J_ij gives the effect of changes to hat j on the segment change of segment i
			(change-lengths (make-array n)) ;Euclidian lengths of segment changes
			(a-segment-change-changes (make-array n-a))
			(our-segment-changes (make-array n))
			(j 0))
      (mirror-images
       ;;Jacobian gives effect of change with unit Euclidian length.  Instead, we have whatever length comes
       ;;from the calculation, which we want to scale by NGmu
       (dotimes (index n-a)
	 (setf (aref change-lengths j)
	       (* (4vector-Euclidian-length
		   (setf (aref a-segment-changes index) (a-segment-change a-hats a-dsigmas b-hats b-dsigmas index)))
		  NGmu))
	 (format t ".") (force-output)
	 (incf j)))
      (dotimes (i n)
	(setf (aref our-segment-changes i)
	      (loop with total = (make-zero-4vector)
		    for j below n
		    do (4vector-incf total (4vector-scale (aref Jacobian i j) (aref change-lengths j)))
		    finally (return total))))
      ;;Never have a prediction of what the changes should be after the first step
      (multiple-value-bind (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas) ;Do first step
	  (gravitational-backreaction a-hats a-dsigmas b-hats b-dsigmas NGmu)
	(multiple-value-setq (new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas)
	  (close-4vectors new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas))
	(mirror-images
	 (dotimes (index n-a)
	   (setf (aref a-segment-change-changes index) ;Difference in changes from new position
		 (4vector- (a-segment-change new-a-hats new-a-dsigmas new-b-hats new-b-dsigmas index)
			   (aref a-segment-changes index)))))
	(values
	 (concatenate 'vector a-segment-change-changes b-segment-change-changes) ;Actual changes
	 our-segment-changes)))))
	    
(defun describe-backreaction-Jacobian (a-hats a-dsigmas b-hats b-dsigmas)
  (let* ((Jacobian (backreaction-Jacobian a-hats a-dsigmas b-hats b-dsigmas))
	 (n (array-dimension Jacobian 0)))
    (dotimes (i n)
      (dotimes (j n)
	(setf  (aref Jacobian i j)  (4vector-Euclidian-length (aref Jacobian i j)))))
    (plot-matrix Jacobian)
    Jacobian))

;;Compute gravitational 4-momentum spectrum.  See also compute-radiated-energy
(defun backreaction-compute-radiation (directory &key (start 0) end (step 1) steps
						 (direct-n (expt 2 14)) (total-n (expt 2 40)) (bin-size 2.0)
						 (split-levels 0) submit)
   (cond (submit
	  (copy-files-to-run directory)
	  (do-submit directory (format nil "rad ~A" (backreaction-directory-name directory))
		     "radiation-"
		     batch-flags
		     `(backreaction-compute-radiation-1
		       nil ,start ,end ,step ',steps ,direct-n ,total-n ,bin-size ,split-levels)))
	 (t (backreaction-compute-radiation-1 directory start end step steps direct-n total-n bin-size split-levels))))

;;Compute from START to END in steps of STEP.  Or if STEPS is set, use those steps instead.  If END is not
;;set, stop when we find a step whose data is not there
(defun backreaction-compute-radiation-1 (directory start end step steps direct-n total-n bin-size split-levels)
  (format t "~%Running on host ~A~%~%" (machine-instance))
  (with-group-write-access
   (unless directory			     ;NIL means connected directory
     (setq directory (sb-posix:getcwd)))     ;Using "." causes some trouble with batch job submission
   (let ((start-time (get-internal-real-time))
	 (*backreaction-filename-digits* (find-backreaction-digits directory)))
     (if steps
	 (loop for this-step in steps
	       do (backreaction-compute-radiation-2 directory this-step direct-n total-n bin-size split-levels))
       (loop for this-step from start by step
	     until (and end (>= this-step end)) ;Stop if end reached
	     while (or end	;If end not set
		       (probe-file (backreaction-dump-filename directory this-step))) ;Go until end of run
	     do (backreaction-compute-radiation-2 directory this-step direct-n total-n bin-size split-levels)))
     (format t "~&Computation took ~A."
	     (format-seconds (floor (- (get-internal-real-time) start-time) internal-time-units-per-second)))
     )))

(defun backreaction-compute-radiation-2 (directory step direct-n total-n bin-size split-levels)
  (format t "Step ~D: " step)
  (with-open-file (stream (backreaction-momentum-filename directory step)
			  :direction :output :if-exists :supersede)
    (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
	(read-hats-dsigmas (backreaction-dump-filename directory step))
      (mirror-image-let ((a-sigmas (dsigma-to-sigma a-dsigmas)))
	(let ((*print-readably* t))	;Compute and write result for this step
	  (format stream 
		  ";;4-momentum computed with split-levels ~D, direct-n ~D, total-n ~D, bin-size ~S~@[, bpd-split ~D~]~%"
		  split-levels direct-n total-n bin-size (and (> *bpd-split-factor* 1) *bpd-split-factor*))
	  (format stream "~S~%"
		  (total-gravitational-4vector-spectrum a-hats a-sigmas b-hats b-sigmas
							direct-n :total-n total-n
							:bin-size bin-size :split-levels split-levels))
	  )))))


;;Plot previously computed energy, momentum, and rocket fraction over time.
;;If DIVIDE-SPECTRUM given it should be a list of in numbers and we show also the spectrum divided at
;;these positions.
(defun plot-rocket-fraction (directory &rest gnuplot-keys
				       &key (start 0) end (step 1) show
				       divide-spectrum ;Split up at given bins
				       (style :lines)
				       &allow-other-keys)
  (multiple-value-bind (momenta bins intersections lifetime)
      (get-backreaction-radiation-data
       directory :start start :end end :step step :show show :divide-spectrum divide-spectrum)
    (let ((division-count (if divide-spectrum (1+ (length divide-spectrum)) 0)))
      ;;0 = intersections.  1,2,3 = overall.  Others are divisions
      (apply #'gnuplot (1+ (* 3 (1+ division-count)))
	     (max (length momenta) (* 3 (length intersections)))
	     #'(lambda (plot point)
		 (if (zerop plot)
		     ;;Show intersections
		     (if (eq point :title) "intersections"
		       (multiple-value-bind (index part) (floor point 3)
			 (let ((data (nth index intersections)))
			   (and data
				(ecase part
				  (0 (values (nth index intersections) 0.0))
				  (1 (values (nth index intersections) 0.1))
				  (2 nil)))))) ;break line
		   (multiple-value-bind (line division) (floor (1- plot) (1+ division-count))
		     (if (eq point :title)
			 (format nil "~A~A"
				 (cond ((plusp division)
					(format nil "~D..~D "
						(if (= division 1) 0 (nth (- division 2) divide-spectrum))
						(1- (or (nth (1- division) divide-spectrum) bins)))) ;last bin
				       (divide-spectrum "all ")
				       (t "")) ;Don't label if not dividing
				 (nth line '("{/Symbol G}" "{/Symbol G}_P" "{/Symbol G}_P/{/Symbol G")))
		       (let ((data (nth point momenta)))
			 (and data
			      (destructuring-bind (age momentum) data
				(values (if age (/ age lifetime) (* point step))
					(let ((this-momentum (aref momentum division)))
					  (ecase line
					    (0 (4vector-t this-momentum))
					    (1 (3vector-length this-momentum))
					    (2 (/ (3vector-length this-momentum) ;Fraction of total
						  (4vector-t (aref momentum 0))))))))))
		       ))))
	     :prelude (format nil "set ytics nomirror~%set y2tics~%set y2label 'rocket fraction'~%")
	     :xlabel (if lifetime "fraction of lifetime" "step")
	     :ylabel "{/Symbol G}"
	     ;;Kluge: axes is not really part of style.  y2 means to use the axes for rocket fraction
	     :styles (append '("lines lt -1 axes x1y2") ;Intersections with fraction axes
			     (loop for division to division-count
				   collect (format nil "~(~A~) lt ~D" style (1+ division))) ;Gamma: solid
			     (loop for division to division-count
				   collect (format nil "~(~A~) lt ~D dt 2" style (1+ division))) ;Gamma_P Dashed
			     (loop for division to division-count			  ;fraction: dotted
				   collect (format nil "~(~A~) lt ~D dt 3 axes x1y2" style (1+ division))
				   ))
	     :title directory
	     gnuplot-keys))))

;;Plot radiation for several runs
(defun plot-backreaction-radiation (directories &rest gnuplot-keys
						&key (start 0) end (step 1) (styles :lines)
						&allow-other-keys)
  (let ((data
	 (loop for directory in directories
	       collect (multiple-value-bind (momenta bins intersections lifetime)
			   (get-backreaction-radiation-data directory :start start :end end :step step)
			 (declare (ignore intersections bins))
			 (list momenta lifetime)))))
    (apply #'gnuplot (length directories)
	   (loop for (momenta) in data maximize (length momenta))
	   #'(lambda (plot point)
	       (if (eq point :title)
		   (format nil "~A" (nth plot directories))
		 (destructuring-bind (momenta lifetime) (nth plot data)
		   (let ((entry (nth point momenta)))
		     (and entry
			  (destructuring-bind (age momentum) entry
			    (values (/ age lifetime)
				    (let ((this-momentum (aref momentum 0)))
				      (4vector-t this-momentum))))
			  )))))
	   :xlabel "fraction of lifetime"
	   :ylabel "{/Symbol G}"
	   :styles styles
	   gnuplot-keys)))

;;Read in and process spectra from backreaction-compute-radiation
;;Returns values (momenta bins intersections lifetime)
;;momenta is a list of (age data) and the data is a vector of 4vectors giving the total power at the given age
;;and then the power in each division.  bins is the number of bins (not divisions) in the datafiles.
;;intersections is a list of evaporation fractions where there were intersections and lifetime is the
;;extrapolated total lifetime of the loop.
(defun get-backreaction-radiation-data (directory &key (start 0) end (step 1) show divide-spectrum)
  (loop with *backreaction-filename-digits* = (find-backreaction-digits directory)
	with initial-info = nil
	with bins
	for this-step from start by step
	until (and end (>= this-step end))
	for last-step = nil then step
	for momentum-file = (backreaction-momentum-filename directory this-step)
	for info-file = (backreaction-info-filename directory this-step)
	while (probe-file momentum-file)
	for momentum = (with-open-file (stream momentum-file) (read stream)) 
	for info = (with-open-file (stream info-file) (read stream))
	unless initial-info do (setq initial-info info ;Save first time
				     bins (length momentum)) ;Save bin count
	collect (list (backreaction-info-agegmu info)
		      (sum-gravitational-4vectors momentum divide-spectrum))
	into momenta			;A list of (age momentum-data)
	finally
	(unless momenta (error "No momentum files found in ~A" directory))
	(unless (cdr momenta) (error "Only one momentum file found in ~A.  Is step right?" directory))
	(let* ((final-age (first (car (last momenta))))
	       (final-fraction-lost (- 1.0 (/ (backreaction-info-length info)
					      (backreaction-info-length initial-info))))
	       (lifetime (and final-age	;old runs lack this
			      (/ final-age final-fraction-lost))) ;Total age of loop before complete evaporation
	       (intersections
		(loop for step in (backreaction-intersection-steps directory)
		      for age = (backreaction-info-agegmu
				  (with-open-file (stream (backreaction-info-filename directory step))
				    (read stream)))
		      collect (if age (/ age lifetime) step))))
	  (when show			;Show table
	    (format t "~:[Step~;Age~] Gamma Gamma_P Fraction~%" (first (first momenta)))
	    (loop for (age momentum) in momenta
		  for this-step from 0 by step
		  do (format t "~S ~S ~S ~S~%"
			     (if age (/ age lifetime) this-step)
			     (4vector-t momentum)
			     (3vector-length momentum)
			     (/ (3vector-length momentum) (4vector-t momentum)))))
	  (return (values momenta bins intersections lifetime)))))


(defun plot-backreaction-radiation-spectrum (directory step &rest gnuplot-keys)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (momenta (with-open-file (stream (backreaction-momentum-filename directory step))
		    (read stream)))
	 (spectrum (map 'vector #'4vector-t momenta)))
    (apply #'plot-binned-gravitational-spectrum
	   spectrum
	   2.0
	   gnuplot-keys)))

(defun get-cusp-power (directory step)
  (total-cusp-gravitational-power-coefficient
   (let ((*backreaction-filename-digits* (find-backreaction-digits directory)))
     (multiple-value-call #'find-cusp-parameters
			  (read-hats-dsigmas (backreaction-dump-filename directory step))))))

;;Integrate the cusp spectrum over a bin
(defun cusp-power-in-bin (bin cusp-coefficient)
  (* 3 (- 1 (expt 2.0 -1/3)) (expt 2.0 (* bin -1/3)) ;\int_n^{2n} j^{-4/3} dj with n = 2^bin
     cusp-coefficient))

;;It doesn't work to start with step 0, because there are zero-angle kinks due to splitting of hats
(defun plot-cusp-coefficient (directory &key start end (step 1))
  (unless start (error "You must specify :start.  Starting at 0 doesn't work because of duplicate points"))
  (plot-data-list
   (loop with *backreaction-filename-digits* = (find-backreaction-digits directory)
	 for this-step from start by step
	 until (and end (>= this-step end))
	 while (probe-file (backreaction-dump-filename directory this-step))
	 collect (list this-step (get-cusp-power directory this-step)))))

(defun plot-backreaction-radiation-spectrum-cusps (directory step &rest gnuplot-keys)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (momenta (with-open-file (stream (backreaction-momentum-filename directory step))
		    (read stream)))
	 (spectrum (map 'vector #'4vector-t momenta))
	 (cusp-coefficient (get-cusp-power directory step)))
    (apply #'gnuplot 2 (length spectrum)
	   #'(lambda (plot point)
	       (and (numberp point)
		    (let ((n (expt 2.0 point))) ;start of bin
		      (values n
			      (ecase plot
				(0 (aref spectrum point))
				(1 (cusp-power-in-bin point cusp-coefficient)))))))
	   :logscale '(:x :y)
	   :styles '(:points :lines)
	   :prelude (format nil "set format '%h'~%")
	   gnuplot-keys)))

;;Find bin where cusp radiation matters
(defun find-backreaction-radiation-cusp-significance (directory step)
  (let* ((*backreaction-filename-digits* (find-backreaction-digits directory))
	 (file (backreaction-momentum-filename directory step))
	 (momenta (with-open-file (stream (cond ((probe-file file) file)
						((= step 1) ;I'm not sure why the regular file isn't there
						 (format nil "~A/0.0-percent-momentum.dat" directory))
						(t file))) ;Get error
			(read stream)))
	 (spectrum (map 'vector #'4vector-t momenta))
	 (cusp-coefficient (get-cusp-power directory step)))
    (loop for bin downfrom (1- (length spectrum)) to 0
	  when (< (cusp-power-in-bin bin cusp-coefficient) (aref spectrum bin)) ;cusp power less
	  return (1+ bin))))
	  

(defun backreaction-intersection-steps (directory)
  (loop with *backreaction-filename-digits* = (find-backreaction-digits directory)
	for file in (directory (format nil "~A/*simulation*input.dat" directory))
	collect (parse-integer (pathname-name file) :end *backreaction-filename-digits*)))
											  
;;Total 4vectors in array.  If divide-spectrum is given, also split them at given bins.
;;Always returns vector of sums
(defun sum-gravitational-4vectors (spectrum &optional divide-spectrum)
  (let* ((divisions (if divide-spectrum (+ 2 (length  divide-spectrum)) 1))
	 (result (coerce (loop repeat divisions collect (make-zero-4vector)) 'vector))
	 (divide-remaining divide-spectrum))
    (loop with index = 1		;first entry in division section.
	  for bin from 0 below (length spectrum)
	  do (4vector-incf (aref result 0) (aref spectrum bin)) ;Entry 0 is grand total
	  when divide-spectrum do
	  (when (and divide-remaining (>= bin (first divide-remaining)))
	    (incf index)		;On to next division
	    (pop divide-remaining))
	  (4vector-incf (aref result index) (aref spectrum bin)))
    result))

;;Test code

(defun run-jj-backreaction (set number
				&key submit threads (loss 0.5) (dump 1e-5) (tolerance 1e-4) (maximize-errors t) 
				;;The system that guesses the step-size usually gets it much too large.
				(step-size 1e-3) ;This would mean 100 dump steps for each ODE step
				(minimum-loop-count 7) ;No loops with fewer segments than this
				(root-directory "backreaction")
				;;split with c coefficient.  To turn off this feature you must give this as NIL explicitly
				(split-c 4e-3)
				min-length   ;Split up segments to have this minimum length
				big-kink-pieces ;If set, split kinks larger than threshold into this many pieces
				kink-threshold
				max-kink ;Split kinks so angles could be smoothed to more than this
				no-intersections
				dry-run		   ;just a split, don't run
				resume		   ;restart with last dumped file
				(use-0 nil)) ;Use the loop-N-0.dat file
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
      (read-jj-loop-hats set number :use-0 use-0)
    (format t "Read in N_A = ~D, N_B = ~D~%" (length a-hats) (length b-hats))
    (when big-kink-pieces
	(mirror-images (multiple-value-setq (a-hats a-dsigmas)
			 (split-large-kinky-hats a-hats a-dsigmas big-kink-pieces kink-threshold))))
    (when max-kink
	(mirror-images (multiple-value-setq (a-hats a-dsigmas)
			 (split-kinky-hats a-hats a-dsigmas max-kink))))
    (when min-length
      (mirror-images (multiple-value-setq (a-hats a-dsigmas) (split-long-hats a-hats a-dsigmas min-length))))
    (when split-c
      (mirror-images (multiple-value-setq (a-hats a-dsigmas) (split-hats-c a-hats a-dsigmas split-c))))
    (format t "N_A = ~D, N_B = ~D~%" (length a-hats) (length b-hats))
    (unless dry-run
      (ode-gravitational-backreaction a-hats a-dsigmas b-hats b-dsigmas
				      (format nil "~A/~D/~D" root-directory set number)
				      loss dump :submit submit :threads threads
				      :tolerance tolerance :maximize-errors maximize-errors
				      :step-size step-size :no-intersections no-intersections
				      :minimum-loop-count minimum-loop-count :resume resume))))

;;Read given JJ loop, return hats and sigmas in rest frame
(defun read-jj-loop-hats (set number
			       &key use-0 ;is this still needed?
			       hats)	  ;New convention
  (cond (hats
	 (multiple-value-call #'rest-frame-close-hats
			      (read-hats-dsigmas
			       (format nil "/cluster/tufts/strings/backreaction/inputs/~D/loop-hats-~D.dat" 
				       set number))))
	(t
	 (read-one-loop (format nil "/cluster/tufts/strings/backreaction/inputs/~D/loop-~D~:[~;-0~].dat" 
				set number use-0))
	 (multiple-value-call #'rest-frame-close-hats (get-ab-data-dsigmas (longest-loop))))))


(defun read-jj-loop-ratios ()
  (with-open-file (stream "/cluster/tufts/strings/backreaction/inputs/loops-sorted-by-ratio.text")
    (loop for entry =  (loop repeat 6 for field = (read stream nil) while field collect field)
	  while entry
	  collect entry)))

(defun do-jj-backreaction-run (skip count &rest arguments)
  (loop repeat count for (set number) in (nthcdr skip (read-jj-loop-ratios))
	do (format t "********** ~D/~D **********~%" set number)
	do (apply #'run-jj-backreaction set number arguments)))
						 

(defun jj-compute-radiation (root-directory skip count &rest arguments &key (start 0) &allow-other-keys)
  (loop repeat count for (set number) in (nthcdr skip (read-jj-loop-ratios))
	for directory = (format nil "~A/~D/~D" root-directory set number)
	for *backreaction-filename-digits* = (find-backreaction-digits directory)
	do (format t "~D/~D: " set number)
	do (cond ((not (probe-file (backreaction-success-filename directory)))
		  (format t "Backreaction running or failed~%"))
		 ((probe-file (backreaction-momentum-filename directory start)) ;If first file is there
		  (format t "Backreaction succeeded and momenta already computed~%"))
		 (t
		  (format t "Backreaction succeeded.  Computing~%")
		  (apply #'backreaction-compute-radiation directory arguments)
		  ))))

(defun jj-loop-effort ()
  (loop for (set number ratio n-a n-b length) in (read-jj-loop-ratios)
	sum (* n-a n-b (+ n-a n-b))))	;Number of segment changes to compute
	

(defun jj-backreaction-directories (root-directory start end)
  (loop for (set number) in (subseq (read-jj-loop-ratios) start end)
	collect (format nil "~A/~D/~D" root-directory set number)))

(defun jj-backreaction-info (root-directory count)
  (loop repeat count
	for (set number) in (read-jj-loop-ratios)
	with intersections and max-v and ejected and major-loss and initial-length
	do (format t "~%*** ~D/~D~%" set number)
	do (multiple-value-setq (initial-length intersections major-loss max-v ejected)
	     (total-backreaction-info (format nil "~A/~D/~D" root-directory set number)))
	collect intersections into all-intersections
	collect max-v into all-max-v
	collect ejected into all-ejected
	collect major-loss into all-major-loss
	finally (flet ((compare (x y)
				(or (null x) ;if X NIL, put first
				    (and y   ;if Y NIL, that first
					 (< x y)))))
		  (return (values (sort all-intersections #'compare)
				  (sort all-major-loss #'compare)
				  (sort all-max-v #'compare)
				  (sort all-ejected #'compare))))))

(defun jj-sort-backreaction-info ()
  (let ((data (read-jj-loop-ratios)))	;list of (set number ratio n-a n-b length)
    (setq data
	  (sort data #'<
		:key #'(lambda (x)
			 (let ((n-a (fourth x))
			       (n-b (fifth x)))
			   (* n-a n-b (+ n-a n-b)) ;difficulty
			   ))))
    (with-open-file (stream "/cluster/tufts/strings/backreaction/inputs/loops-sorted-by-difficulty.text"
			    :direction :output :if-exists :supersede)
      (format stream ";;Set	Number	shortest/longest seg	N_a	N_b	Length~%")
      (dolist (entry data)
	(format stream "~{~S~^	~}~%" entry)))))


(defun jj-status (root-directory start end &key (threshold 60)) ;Number of seconds to consider stalled
  (loop with summary = nil
	for index from start
	for (set number) in (subseq (read-jj-loop-ratios) start end)
	do (format t "~&~D: " index)
	do (incf (getf summary (jj-status-1 root-directory set number :threshold threshold) 0))
	finally (format t "~{~D ~(~A~)~^, ~}~%" (reverse summary))))

(defun jj-status-1 (root-directory set number &key (threshold 60))
  (let* ((directory (format nil "~A/~D/~D" root-directory set number))
	 (log-file (format nil "~A/backreaction-output" directory))
	 (success-file (backreaction-success-filename directory)))
    (format t "~D/~D: " set number)
    (cond ((probe-file success-file)
	   (format t "success file exists~%")
	   :succeeded)
	  ((probe-file log-file)
	   (let* ((written (file-write-date log-file))
		  (ago (- (get-universal-time) written))) ;Number of seconds ago that file was written
	     (cond ((< ago threshold)
		    (format t "Running.  Log file written ~D seconds ago~%" ago)
		    :running)
		   (t (format t "Failed or stalled.  Log file written ~:[~D seconds ago~*~;~*~A~], ending:~%"
			      (> ago 180)
			      ago
			      (sb-int:format-universal-time nil written))
		      (do-run-program "tail" :args (list "-3" log-file))
		      (terpri)
		      (do-run-program "squeue" :args (list "-hn" (format nil "br ~D/~D" set number))) ;show job if any
		      :stopped))))
	  (t (format t "not started~%")
	     :not-started))))

;;Compute radiation of loop from simulation as though it were first frame of backreaction
;;If hats is set, hats are already stored in files.  If not, we compute them.
(defun compute-radiation-without-backreaction (set number &key (split-levels 0) submit
						   (hats (>= set 7))) ;At this point we stored hats and dsigmas
  (multiple-value-bind (a-hats a-dsigmas b-hats b-dsigmas)
      (read-jj-loop-hats set number :hats hats)	;Read loop, call rest-frame-close-hats
    (let ((*backreaction-filename-digits* 1)		    ;Since not really doing it
	  (directory (format nil "~A/backreaction/~D/~D" batch-root-directory set number)))
      (ensure-directories-exist (format nil "~A/" directory))
      (write-hats-dsigmas a-hats a-dsigmas b-hats b-dsigmas ;Write as initial step
			(backreaction-dump-filename directory 0))
      (backreaction-compute-radiation directory :start 0 :end 1 :split-levels split-levels :submit submit))))

;;Returns a list of (percentage (run-number step-number) (run-number step-number) ...)
;;Mapping percentages to the steps of backreaction
(defun read-steps-file (&optional (file "/cluster/tufts/strings/backreaction/6-loops_step-evap_no-majors_250to400.dat"))
  (with-open-file (stream file)
    (loop with table = (make-hash-table :test #'eql)
	  with result = nil
	  for line = (read-line stream nil)
	  while line
	  for (set number step percentage) =
	      (with-input-from-string (scanner line)
		(loop repeat 4 collect (read scanner)))
	  do (push (list number step) (gethash percentage table))
	  finally
	  (maphash #'(lambda (percentage list)
		       (push (cons percentage (reverse list)) result))
		   table)
	  (return (sort result #'< :key #'car)))
	  ))

(defun list-backreaction-cusp-significance (&optional (directory "/cluster/tufts/strings/backreaction/6"))
  (loop for (percentage . items) in (read-steps-file)
	do (format t "~&~$: " percentage)
	do (loop for (number step) in items
		 for significance = 
		 (handler-case (find-backreaction-radiation-cusp-significance
				(format nil "~A/~D" directory number)  step)
		   (file-does-not-exist ()
                     (format t "No file for ~D step ~D:~%" number step)
		     :failed))
		 when (or (null significance)
			  (and (numberp significance) (< significance 15)))
		 do (format t "~D step ~D: ~A~%" number step significance)
		 )))

