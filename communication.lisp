;;;Low level file I/O, for communication with successor processors and snapshots
(in-package "CL-USER")

;;Bound during sending
(defvar *send-destination*)		;Stream or vector for output

(defvar *send-offset* zero-4vector)	;Vector which is subtracted from coordinates being output
(defvar *receive-handle-job-offset* 0)	;Added to handles being read from predecessors

;;Bound during reception
(defvar *receive-source* nil)		;Stream or vector for input
(defvar *receive-index* 0)		;Next byte to be read from array
(defvar *allow-missing-predecessor-files* nil) ;Don't give error on missing predecessor files (e.g., in reading dumps)
(defvar *count-missing-predecessor-files*)     ;If file missing, count here

(declaim (type fixnum *send-index* *receive-index* *receive-length*))

;;Package for symbols created by our macros.
(defconstant communicator-package (find-package "CL-USER"))

;;Add a byte to message being composed
(defun send-write-byte (byte)
  (assert (>= byte 0))
  (cond ((typep *send-destination* 'vector) ;Writing to array?
	 (unless (vector-push byte *send-destination*)
	   (error "Ran out of space in output buffer")))
	(t				;No: file
	 (write-byte byte *send-destination*))))

(defun receive-read-byte (&optional (eof-error-p t))
  (cond ((typep *receive-source* 'vector) ;Reading from array?
	 (cond ((< *receive-index* (length *receive-source*))	;something left to read?
		(prog1 (aref *receive-source* *receive-index*)
		  (incf *receive-index*)))
	       (eof-error-p
		(locally (declare (muffle-conditions compiler-note)) ;Avoid unused code warning with constant flag
		  (error "Read over end of message")))
	       (t nil)))		  ;Error not requested -- return NIL
	(t				;No: file
	 (read-byte *receive-source* eof-error-p)
	 )))

(defstruct communicator-argument
  type
  send
  receive
  )
  
(defmacro define-communicator-argument (name arglist send-form receive-form)
  (assert (and (listp arglist) (= (length arglist) 1)))
  (let ((receiver (intern (format nil "RECEIVE-ARGUMENT-~S" name) communicator-package))
	(sender (intern (format nil "SEND-ARGUMENT-~S" name) communicator-package))
	(argument (intern (format nil "COMMUNICATOR-ARGUMENT-~S" name) communicator-package)))
  `(progn
     (defun ,sender ,arglist
       ,send-form)
     (defun ,receiver ()
       ,receive-form)
     (defparameter ,argument
       (make-communicator-argument :type ',name :send ',sender :receive ',receiver)))))


;;;Here's the system which reads and writes the data from and to the files

#-(or little-endian big-endian)#.(error "endianness not known")

(defparameter conversion-data
  (alien-sap (make-alien (signed 8) 8))) ;Static pointer to 8 bytes

(declaim (inline write-unsigned-n write-convert-n read-unsigned-n read-convert-n
		 read-single read-double write-single write-double))

;;read N bytes into conversion data, LSB first.
(defun read-convert-n (n)
  (declare (optimize (speed 3))
	   (fixnum n))
  (setf (sap-ref-64 conversion-data 0) 0)	;Clear old.  May not exist on 32-bit machine.
  (loop for byte below n
	;;Put in right place according to endianness
	do (setf (sap-ref-8 conversion-data #+little-endian byte #+big-endian (- 7 byte))
		 (receive-read-byte))))

(defun read-unsigned-n (n)
  (declare (optimize (speed 3) (safety 0))) ;Don't check conversion to fixnum below
  (read-convert-n n)
  (the fixnum (sap-ref-64 conversion-data 0)))

;;functions for file I/O
(defun read-single (stream)
  (declare (muffle-conditions compiler-note)) ;Avoid warning about consing return value
  (declare (optimize (speed 3)))
  (setf (sap-ref-32 conversion-data 0) (read-byte stream))
  (sap-ref-single conversion-data 0))

(defun read-double (stream)
  (declare (muffle-conditions compiler-note)) ;Avoid warning about consing return value
  (declare (optimize (speed 3)))
  (setf (sap-ref-64 conversion-data 0) (read-byte stream))
  (sap-ref-double conversion-data 0))

;;Read elements of 4vector
(defun read-4vector (stream)
  (let ((result (make-4vector)))
    (dotimes (i 4)
      (setf (4vector-component result i) (read-double stream)))
    result))

;;Write a single-float as binary output
;;If arguement is not already a single-float, it is coerced to be
(defun write-single (single stream)
  ;(declare (optimize (speed 3)))
  (setf (sap-ref-single conversion-data 0) (single-float single))
  (write-byte (sap-ref-32 conversion-data 0) stream))

;;Write a double-float as binary output
;;If arguement is not already a double-float, it is coerced to be
(defun write-double (double stream)
  ;(declare (optimize (speed 3)))
  (setf (sap-ref-double conversion-data 0) (double-float double))
  (write-byte (sap-ref-64 conversion-data 0) stream))

;;Write elements of 4vector
(defun write-4vector (vector stream)
  (dotimes (i 4)
    (write-double (4vector-component vector i) stream)))

;;Write N bytes out of conversion data, LSB first.
(defun write-convert-n (n)
  (declare (optimize (speed 3))
	   (fixnum n))
  (loop for byte below n
	;;Get from right place according to endianness
	do (send-write-byte (sap-ref-8 conversion-data
				       #+little-endian byte #+big-endian (- 7 byte)))))

;;Write a positive integer as N bytes
(defun write-unsigned-n (n value)
  (declare (optimize (speed 3))
	   (fixnum n value))
  (assert (>= value 0))
  (when (< n 8)
    (locally (declare (type (integer 1 7) n))
      (assert (< value (the fixnum (ash 1 (* 8 n)))))))
  (setf (sap-ref-64 conversion-data 0) value)
  (write-convert-n n))

;;Write boolean as 0 = NIL, 1 = true.
(define-communicator-argument boolean (boolean)
  (send-write-byte (if boolean 1 0))
  (not (zerop (receive-read-byte))))

(define-communicator-argument byte (byte)
  (send-write-byte byte)
  (receive-read-byte))

(define-communicator-argument fixnum (value)
  (progn
    (setf (sap-ref-64 conversion-data 0) value)
    (write-convert-n 8))
  (progn
    (read-convert-n 8)
    (sap-ref-64 conversion-data 0)))

(define-communicator-argument 4vector (4vector)
  (dotimes (i 4)
    (send-argument-float (- (aref 4vector i) (aref *send-offset* i))))
  (let ((4vector (make-4vector)))
    (dotimes (i 4) (setf (aref 4vector i) (receive-argument-float)))
    4vector))

(define-communicator-argument global-4vector (4vector)
  (dotimes (i 4)
    (send-argument-float (aref 4vector i)))
  (let ((4vector (make-4vector)))
    (dotimes (i 4) (setf (aref 4vector i) (receive-argument-float)))
    4vector))

(define-communicator-argument global-4vector-or-nil (4vector)
  (cond ((null 4vector)			;None
	 (send-argument-byte 0))	;say that
	(t
	 (send-argument-byte 1)		;say we have 4vector
	 (send-argument-global-4vector 4vector))) ;send it
  (ecase (receive-argument-byte)	 ;0 = nil, 1 = number
    (0 nil)
    (1 (receive-argument-global-4vector))))

#| Not used
(define-communicator-argument 3vector (3vector)
  (dotimes (i 3) (send-argument-float (- (aref *send-offset* i) (aref 3vector i))))
  (let ((3vector (make-3vector)))
    (dotimes (i 3) (setf (aref 3vector i) (receive-argument-float)))
    3vector))
|#

;;Number of bytes to write quantities
(defconstant handle-length 4)		;handle code
(defconstant job-number-length 4)
;;Use high bits to distinguish predecessor sets
;;This must be larger than number of jobs/run being combined, and smaller than 256^job-number-length/combine
(defparameter multiple-predecessor-job-offset-multiplier 100000000)
(defconstant vv-index-length 2)		;index into VV array
(defconstant site-length 6)		;site = 3 VV indices
(defconstant string-number-length 4)	;String and diamond numbers for explicit initial conditions
(defconstant diamond-number-length 4)	;String and diamond numbers for explicit initial conditions

(define-communicator-argument float (value)
  (progn
    (setf (sap-ref-double conversion-data 0) value)
    (write-convert-n 8))
  (progn
    (read-convert-n 8)
    (sap-ref-double conversion-data 0)))

(define-communicator-argument float-or-nil (number)
  (cond ((null number)			;None
	 (send-argument-byte 0))	;say that
	(t
	 (send-argument-byte 1)		;say we have number
	 (send-argument-float number))) ;send it
  (ecase (receive-argument-byte)	 ;0 = nil, 1 = number
    (0 nil)
    (1 (receive-argument-float))))

(define-communicator-argument job-number (number)
  (write-unsigned-n job-number-length number)
  (read-unsigned-n job-number-length))

(define-communicator-argument handle (handle)
  (progn
    (debug-trace-handle handle "sending")
    (send-argument-job-number (handle-creator handle)) ;Job who created handle
    (write-unsigned-n handle-length (handle-code handle))	;handle code
    )
  (debug-trace-handle
   (make-handle (+ (receive-argument-job-number) *receive-handle-job-offset*)
		(read-unsigned-n handle-length))
   "received"))

(define-communicator-argument rejoining-junction (junction)
  (progn
    (send-argument-handle (rejoining-junction-handle junction))
    (send-argument-byte (rejoining-junction-left-direction junction))
    (send-argument-byte (rejoining-junction-right-direction junction))
    (send-argument-float (rejoining-junction-a junction))
    (send-argument-float (rejoining-junction-b junction))
    (send-argument-float (rejoining-junction-a0 junction))
    (send-argument-float (rejoining-junction-b0 junction))
    (send-argument-float-or-nil (rejoining-junction-dump-time junction))
    )
  (make-rejoining-junction
   :handle (receive-argument-handle)
   :left-direction (receive-argument-byte)
   :right-direction (receive-argument-byte)
   :a (receive-argument-float) :b (receive-argument-float)
   :a0 (receive-argument-float) :b0 (receive-argument-float)
   :dump-time (receive-argument-float-or-nil)))

(define-communicator-argument tag-count (count)
  (write-unsigned-n tag-count-bytes count)
  (read-unsigned-n tag-count-bytes))

(define-communicator-argument loop-tag (tag)
  (send-argument-handle (tag-handle tag))
  (handle-object (receive-argument-handle)))

;;VV junction, which is really just a face specification.
(define-communicator-argument vv-junction (junction)
  (progn
    (send-argument-site (vv-junction-site junction))
    (send-argument-byte (vv-junction-axis1 junction))
    (send-argument-byte (vv-junction-axis2 junction))
    )
  (make-vv-junction
   :site (receive-argument-site) :axis1 (receive-argument-byte) :axis2 (receive-argument-byte)))

;;Send information about the bh, center and the radius 
(define-communicator-argument blackhole-point (blackhole)
  (progn
    (send-argument-global-4vector (blackhole-center blackhole))
    (send-argument-float   (blackhole-size blackhole))
    )
  (make-blackhole
   :center (receive-argument-global-4vector) :size (receive-argument-float)))

;;Checks if there is BH information and send nil or the BH info
(define-communicator-argument blackhole-or-nil (blackhole)
  (cond ((null blackhole)                 ;None
         (send-argument-byte 0))        ;say that
	(t
	 (send-argument-byte 1)         ;say we have 4vector
         (send-argument-blackhole-point blackhole))) ;send it
  (ecase (receive-argument-byte)         ;0 = nil, 1 = number                                                                                                                                              
    (0 nil)
    (1 (receive-argument-blackhole-point))))


;;Junction from explicit initial conditions.  It is like a VV-junction, but instead it has string and diamond numbers.
(define-communicator-argument initial-junction (junction)
  (progn
    (write-unsigned-n string-number-length (initial-junction-string junction)) ;Number of string
    (write-unsigned-n diamond-number-length (initial-junction-diamond junction))) ;Number of diamond
  (make-initial-junction
   :string (read-unsigned-n string-number-length)
   :diamond (read-unsigned-n string-number-length)
   ))

;;Face is a list (SITE DIRECTION1 DIRECTION2).  The orientation does not have to be canonical.
(define-communicator-argument face (face)
  (progn
    (send-argument-site (first face))
    (send-argument-byte (second face))
    (send-argument-byte (third face)))
  (list (receive-argument-site) (receive-argument-byte) (receive-argument-byte)))

;;One of the corners of a diamond.  We send it  is also a zip code saying that it is is an existing corner of
;;*send-previous-diamond* or *send-first-diamond*.
(define-communicator-argument diamond-corner (point)
  (let ((code (cond ((and *send-previous-diamond* (eq point (diamond-right *send-previous-diamond*))) 1)
		    ((and *send-previous-diamond* (eq point (diamond-start *send-previous-diamond*))) 2)
		    ((and *send-previous-diamond* (eq point (diamond-end *send-previous-diamond*))) 3)
					;left of previous would not make sense
		    ((and *send-first-diamond* (eq point (diamond-left *send-first-diamond*))) 4) ;Left corner of first
		    ((and *send-first-diamond* (eq point (diamond-start *send-first-diamond*))) 5)
		    ((and *send-first-diamond* (eq point (diamond-end *send-first-diamond*))) 6)
		    (t 0))))
    (send-argument-byte code)		;Send code
    (when (zerop code)			;If we found it we're done.  Otherwise send point.
      (send-argument-4vector point)))
  (ecase (receive-argument-byte)
    (0 (if *reading-dumps* (standardize-position (receive-argument-4vector))
	 (receive-argument-4vector)))
    (1 (diamond-right *receive-previous-diamond*))
    (2 (diamond-start *receive-previous-diamond*))
    (3 (diamond-end *receive-previous-diamond*))
    (4 (diamond-left *receive-first-diamond*))
    (5 (diamond-start *receive-first-diamond*))
    (6 (diamond-end *receive-first-diamond*))))

;;A site is a triple index into Vachaspati-Vilenkin array.  Sent as global position.
(define-communicator-argument site (site)
  (write-unsigned-n site-length (globalize-vv site))
  (localize-vv (read-unsigned-n site-length)))



(defstruct communicator
  name
  receive-function			;function to call with arguments
  send-function
  arglist					;names of arguments
  argument-data)				;list of communicator-arguments structures

;;Defines functions
;;(SEND-name destination . args) which sends a message to destination in such a way that
;;destination will call (RECEIVE-name source . args).   DESTINATION and SOURCE
;;are neighbor structures.  The table of data goes in the variable TABLE.  Sender
;;functions call COMMON-SENDER with the code and arguments.
(defmacro define-communicators (table common-sender &body body)
  `(progn
     (defvar ,table)
     (setq ,table (make-array ,(length body)))
     ,@(loop for (name . arguments) in body
	     for code from 0
	     collect `(define-communicators-1 ,table ,common-sender ,code ,name ,arguments))))

(defmacro define-communicators-1 (table common-sender code-value name arguments)
  (let ((receiver (intern (format nil "RECEIVE-~S" name) communicator-package))
	(sender (intern (format nil "SEND-~S" name) communicator-package))
	(code (intern (format nil "CODE-~S" name) communicator-package))
	(arglist (mapcar #'car arguments))
	(argument-data (mapcar #'(lambda (spec)
				   (intern (format nil "COMMUNICATOR-ARGUMENT-~S" (second spec))
						   communicator-package))
			       arguments)))
    `(progn
       (defconstant ,code ,code-value)
       (defun ,sender (*send-destination* ,@arglist)
	 (,common-sender ,code ,@arglist))
       (setf (aref ,table ,code)
	     (make-communicator :name ',name
				:receive-function ',receiver
				:send-function ',sender
				:arglist ',arglist
				:argument-data (list ,@argument-data))))))


;;The list of messages and their arguments.  This is used both for snapshots and
;;communications to successor jobs.
(define-communicators communicators communicator-send
  ;;Successor communication codes 
  (diamond (start diamond-corner)
	   (left diamond-corner)
	   (right diamond-corner)
	   (end diamond-corner)
	   (tag loop-tag)
	   (a-kink-created global-4vector-or-nil)
	   (b-kink-created global-4vector-or-nil)
	   (countup tag-count)
	   (inertp boolean)
	   )
  (start-junction (junction rejoining-junction)) ;Next diamond begins with this junction
  (end-junction (junction rejoining-junction))   ;Previous diamond ends with this junction
  (note-left-junction (junction rejoining-junction)) ;Dump joined to us (moving right) in previous diamond
  (note-right-junction (junction rejoining-junction)) ;We joined to dump (moving right) in previous diamond
  (loop)			     ;Previous diamond has first diamond to its right
  (start-deleted)		     ;first diamond has :deleted to its left
  (end-deleted)			     ;previous diamond has :deleted to its right
  (tag (handle handle) (created-position global-4vector) ;Send loop tag info for later use
       (last-position global-4vector-or-nil) (xi float)
       (minimum-inert-count tag-count) (dumped boolean) (bh blackhole-or-nil))
  (dump-time (time float))	     ;Global time of this dump
  ;;Say that we took care of a cell
  (vv-cell-done (site site))
  ;;Install a phase (0, 1, or 2) at site (i, j, k) in the virtual global lattice with boundary.
  (vv-phase (phase byte) (site site))
  ;;Give location of (perturbed) point at center of face.  The indices are global, but the location is pre-offset
  (vv-face-point (location 4vector) (site site) (axis1 byte) (axis2 byte))
  ;;Information about a diamond that was designed but not really created by sender or a predecessor
  ;;Faces are outward directed from the diamond.
  (pseudo-diamond (start 4vector) (west-face face) (east-face face))
  (start-vv-junction (junction vv-junction)) ;Say that string starts in the middle of given VV face.
  (end-vv-junction (junction vv-junction)) ;Say that string ends at given face
  (dump-length (length float))		;Total invariant length of string dumped
  (start-initial-junction (junction initial-junction))
  (end-initial-junction (junction initial-junction))
  (start-BH) ;communicator for the bh diamond
  (end-BH) ;communicator for the bh diamond
  (start-BHdeleted) ;communicators for the BHdeleted case
  (end-BHdeleted)   ;communicators for the BHdeleted case
  (start-BHeatit) ;communicator for the case in which we remove the bh diamond
  (end-BHeatit) ;communicator for the case in which we remove the bh diamond
  (start-BHpropdel)
  (end-BHpropdel)
  (bh (blackhole blackhole-or-nil)) ; communicator to send and receive bh information of each diamond
  (myBH)
  )


;;SEND-... comes here with code and arguments.
;;Start outputting a message to file.  *SEND-DESTINATION* is the stream.
(defun communicator-send (code &rest arguments)
  (let* ((communicator (aref communicators code))
	 (*debug-send-communicator* communicator)) ;Store globally for debugging
    (send-write-byte code)		;Say what message this is
    (loop for argument in arguments
	  for argument-data in (communicator-argument-data communicator)
	  for *debug-send-argument-number* from 0
	  ;;Call the sending function for each argument, according to its type
	  do (funcall (communicator-argument-send argument-data) argument))))

;;Do reception: decode arguments and call receiver function
(defun communicator-receive (code)
  (let* ((communicator (aref communicators code))
	 (*debug-receive-communicator* communicator))	;Store globally
    (loop for argument-data in (communicator-argument-data communicator)
	  for *debug-receive-argument-number* from 0
	  collect (funcall (communicator-argument-receive argument-data)) into arguments
	  finally (apply (communicator-receive-function communicator) ;Call receiver
			 arguments))))

;;Instead of actually receiving message, we just tell what is in it
(defun communicator-describe (code)
  (let ((communicator (aref communicators code))
	(*print-pretty* nil))
    (format t "~&~A:~%" (communicator-name communicator))
    (loop for argument-data in (communicator-argument-data communicator)
	  for argument-name in (communicator-arglist communicator)
	  do (format t "  ~A(~A): " argument-name (communicator-argument-type argument-data))
	  do (handler-case	 ;Print argument, but if we can't parse it, don't crash
		 (format t "~S~%"
			 (funcall (communicator-argument-receive argument-data)))
	       (error (condition)
		(format t "~A~%" condition)))
	  )))

;;Terse version
(defun communicator-describe-terse (code)
  (let ((communicator (aref communicators code)))
    (apply (intern (format nil "DESCRIBE-~A" (communicator-name communicator)) communicator-package)
	   (loop for argument-data in (communicator-argument-data communicator)
		 collect (ignore-errors
			   (funcall (communicator-argument-receive argument-data)))
		 ))))

(defun describe-diamond (&rest ignore)
  (declare (ignore ignore))
  (format t "d"))
(defun describe-start-junction (junction)
  (format t "~&~A" (junction-description-terse junction)))
(defun describe-end-junction (junction)
  (format t " ~A~%" (junction-description-terse junction)))
(defun describe-note-left-junction (junction)
  (format t " L=~A" (junction-description-terse junction)))
(defun describe-note-right-junction (junction)
  (format t " R=~A" (junction-description-terse junction)))
(defun describe-tag (&rest ignore)
  (declare (ignore ignore))
  (format t "t"))

(defun describe-bh (blackhole)
  (declare (ignore blackhole))
  (format t " BH~%"))

(defun junction-description-terse (junction)
   (format nil "~D(~D)~@[@~A~]" (handle-code (rejoining-junction-handle junction))
	   (handle-creator (rejoining-junction-handle junction))
	   (and (rejoining-junction-dump-time junction)
		(format-float-reasonably (rejoining-junction-dump-time junction)))))

(defun describe-loop ()
  (format t " LOOP~%"))

(defun describe-start-deleted ()
  (format t " Start-Deleted~%"))

(defun describe-end-deleted ()
  (format t " End-Deleted~%"))

(defun describe-dump-time (time)
  (format t "Time ~A: " (format-float-reasonably time)))

(defun describe-dump-length (length)
  (format t "Total length ~$: " length))

(defun describe-vv-phase (phase site)
  (declare (ignore site))
  (format t "~D" phase))

(defun describe-vv-face-point (location site axis1 axis2)
  (declare (ignore location site axis1 axis2))
  (format t "*"))

(defun describe-pseudo-diamond (start west east)
  (declare (ignore start))
  (format t "~&P~A+~D+~D:~A~D+~D"
	  (site-list (first west)) (second west) (third west)
	  (site-list (first east)) (second east) (third east)))

(defun describe-vv-cell-done (site)
  (declare (ignore site))
  (format t "."))
  
(defun describe-start-vv-junction (junction)
  (format t "~&~A" (vv-junction-description-terse junction)))
(defun describe-end-vv-junction (junction)
  (format t " ~A~%" (vv-junction-description-terse junction)))

(defun vv-junction-description-terse (junction)
   (format nil "~D+~D+~D"
	   (site-list (vv-junction-site junction))
	   (vv-junction-axis1 junction) (vv-junction-axis2 junction)))



;;;Receiving

(define-timer :read)

;;Reads all records in file and calls receiver functions
(defun read-file (filename &key (verbose t))
  (account-time :read
    (with-open-file (*receive-source* filename 
				      :element-type '(unsigned-byte 8) 
				      :if-does-not-exist (unless *allow-missing-predecessor-files* :error))
      (cond (*receive-source*		;Got file?
	     (loop for code = (receive-read-byte nil)
		   while code		       ;exit at EOF
		   do (communicator-receive code) ;Decode arguments and call receiver
		   ))
	    (t				;No file and that is allowed
	     (incf *count-missing-predecessor-files*)
	     (when verbose (format t "No such file: ")))))))

(defun read-file-segments (filename &key (verbose t))
  (account-time :read
		(with-open-file (*receive-source* filename
						  :element-type '(unsigned-byte 8)
						  :if-does-not-exist (unless *allow-missing-predecessor-files* :error))
				(cond (*receive-source*
				       (loop for code = (receive-read-byte nil)
					     while code
					     do (communicator-receive code)
					     )
				       )
				      (t
				       (incf *count-missing-predecessor-files*)
				       (when verbose (format t "No such file: ")))))))

(defun describe-file (filename &optional terse)
  (with-open-file (*receive-source* filename :element-type '(unsigned-byte 8))
    (loop for code = (receive-read-byte nil)
	  while code			    ;exit at EOF
	  do (if terse
		 (communicator-describe-terse code) ;Describe tersely
	       (communicator-describe code) ;Regular 
	  ))))

