(in-package "CL-USER")
(use-package "SB-THREAD")

(defvar *thread-work-counter* 0)

;;Use nthreads threads to do nwork pieces of work.  The function is called with
;;one argument going from 0 to nwork-1.
;;If nthreads is NIL, just call the function nwork times
(defun thread-dispatch (nwork nthreads function)
  (cond (nthreads
	 (setq *thread-work-counter* 0)
	 (mapc #'join-thread			;Wait for threads to finish afterwards
	       (loop for thread below (min nthreads nwork)
		     collect (make-thread #'thread-dispatch-worker :name (format nil "Thread ~D" thread)
					  :arguments (list thread nwork function)))))
	(t				;Non-thread case
	 (dotimes (work nwork)
	   (funcall function work)))))

;;The top-level function in the thread.  Gets the next piece of work in a thread-safe way
;;and calls FUNCTION with the number of the work.
(defun thread-dispatch-worker (thread nwork function)
  (declare (ignorable thread))		;For debugging
  (loop for work = *thread-work-counter*
	while (< work nwork)		;Exit if nothing to do
	;;In sbcl-1.4.0 can just say (compare-and-swap *thread-work-counter* ...)
	when (= (compare-and-swap (symbol-value '*thread-work-counter*) work (1+ work)) work) ;Try to take work
	;;Success: do work.  Otherwise loop to trying again
;;	do (format t "Thread ~D doing work ~D~%" thread work) and
	do (funcall function work)))
	
(defun kill-other-threads ()
  (loop for thread in (list-all-threads)
	unless (eq thread *current-thread*)
	do (terminate-thread thread)))

(defun setup-thread-count ()
  (let ((threads (posix-getenv "SLURM_CPUS_PER_TASK"))) ;Must be given by user with -c for it to work
    (setq *threads* (and threads (parse-integer threads)))))

(setup-thread-count)			;Set up on load
