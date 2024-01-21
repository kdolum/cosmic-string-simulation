;;;Uniform interface to run a subprocess
;;;Some of this is taken from from with-wish.lisp by Matthias Lindner.
;;;Copyright (c) 1998 Ken Olum.

;;This function handles differences in user interface between various lisps.
(defun do-run-program (cmd &key args (wait t) (input t) (output t) error-too)
  "Starts a subprogram.  Arguments:
   :WAIT -- if set (default), waits for program to finish.
   :INPUT and :OUTPUT.  If T (default), inherits stream from Lisp.
      If :STREAM, makes a stream to the process.
      If NIL, use /dev/null
      If stream, use it.
   :ERROR-TOO -- if set, error goes to output stream.  Only implemented on sbcl.
Values:
   HANDLE -- some kind of handle on the process, possibly NIL if not
     waiting.
   ISTREAM -- stream for process's input, or NIL if none
   OSTREAM -- stream for process's output, or NIL if none
    The two streams might be EQ."
  (check-type cmd string)
  (check-type args list)
  (check-type input (or (member t nil :stream) stream))
  (check-type output (or (member t nil :stream) stream))
  #+clisp(when (and wait
		    (not (eq input :stream))
		    (not (eq output :stream)))
	   (error "CLISP can't wait unless streams are used"))

  ;;CLISP has a different name for this.
  #+clisp(when (eq input t) (setq input :terminal))
  #+clisp(when (eq output t) (setq output :terminal))

  #+clisp (let ((stream
		 (run-program cmd :arguments args
			      :input input :output :output)))
	    (values stream stream stream))

  #+(or cmu sbcl)
    (let ((proc (run-program cmd args :input input :output output
				 :error (if (and (not error-too)
						 (member output '(:stream nil)))
					    t ;if pipe or discard, show errors
					  :output) ;but if file, send output also
				 :wait wait
				 #+sbcl :search #+sbcl t)))
	  (values proc (process-input proc) (process-output proc)))

  #-(or clisp cmu sbcl) (error "DO-RUN-PROGRAM not implemented in this Lisp.  Check the source"))
	  
#| This is some code from with-wish and other places.  Since I don't have these implementations I don't have a way to make them work in the above program.

  #+lucid (lcl:run-program cmd
			   :arguments arguments
			   :input     :stream
			   :output    :stream
			   :wait      wait)
  #+allegro (excl:run-shell-command (apply #'vector cmd cmd arguments)
				     :input        :stream
				     :output       :stream
				     :error-output t
				     :wait         wait)
  #+kcl (run-program cmd arguments)

  #+lispworks (if wait 
                (ffi:open-pipe "gnuplot" :direction ...) ...)

|#

(defun do-kill-program (handle)
  "Kill an inferior program giving the handle returned by DO-RUN-PROGRAM"

  #+clisp (progn
	    (close handle)		;Just close the stream
	    ;;Running some null program causes child to be waited for
	    ;;so it doesn't hang in <defunct>.
	    (run-program "echo"))

  #+(or cmu sbcl)
    (let ((in (process-input handle))
	      (out (process-output handle)))
	  (when out (close out))
	  (when in (close in))
	  (process-kill handle #+cmu unix:sigterm #+sbcl sb-unix:sigterm)
	  (process-wait handle))		;Wait for it to die
  t)

(defun wait-for-program (handle)
  "Wait for inferior program to finish"
  #+(or cmu sbcl)(process-wait handle)
  #-(or cmu sbcl)(error "Don't know how to wait in this lisp"))

(defun getpid ()
  #+sbcl(sb-unix:unix-getpid)
  #-sbcl(error "Don't know how to get process ID"))