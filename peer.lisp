;;;Common code for manager-worker communication
(in-package "CL-USER")

(defvar *peer-socket*)			;Socket for worker-manager communications
(defvar *peer-address*)			;Bound to address that we are communicating with


;;Buffers for manager-worker communication
(defparameter peer-buffer-length 400)	;This is enough for :combine 20
(defun make-peer-buffer () (make-array peer-buffer-length :element-type '(unsigned-byte 8) :fill-pointer 0))
(defvar *peer-buffer* (make-peer-buffer)) ;Buffer to use on reception


;;Information about a worker
(defstruct (worker
	    (:print-object print-worker))
  number				;Number of this worker
  status				; :manager -- we're a worker and this is our handle on the manager
					; NIL -- exited or never started
					; :not-started --  job submitted but he hasn't contacted us
					; :running -- running a job
					; :sleep -- waiting for us to give him work
					; :starting -- we are waiting for him to acknowledge wake-job
					; :exiting -- we are waiting for him to acknowledge wake-exit
  job					;Job he is doing if status is :running or :starting
  (buffer (make-peer-buffer))		;Output buffer, containing last message sent
  address				;His IP address and port
  (sequence 0)				;Sequence number of next message to send or what we expect to receive
  (retransmit-time nil)			;When to retransmit to him
  (started nil)				;internal-real-time when worker started
  (exited nil)				;internal-real-time when worker exited
  )

(defun print-worker (worker stream)
  (if *print-escape*
      (print-unreadable-object (worker stream :type t)
	(format stream "~D" (worker-number worker)))
    (format stream "worker ~D" (worker-number worker))
    ))

(defconstant worker-number-length 2)
(defconstant sequence-number-length 4)
(defconstant random-seed-length 4)

(defvar *manager*)			;If worker, the peer structure.  If we are manager, NIL
(defvar *workers*)			;If manager, array of worker structures


;;What a worker does is to evaluate some form, which must look like (FUNCTION . KEYS-AND-VALUES)
;;Here we extract something from the keyword list.  If it is not there, we evaluate DEFAULT
(defmacro get-argument (form key &optional default)
  (let ((indicator (gensym))
	(value (gensym)))
    `(multiple-value-bind (,indicator ,value) (get-properties (cdr ,form) '(,key))
       (if ,indicator ,value		;found it
	 ,default)			;otherwise use default
       )))

(define-communicator-argument worker (worker)
  (write-unsigned-n worker-number-length (worker-number worker))
  (progn (unless *workers* (error "Message intended for manager received by worker"))
	 (aref *workers* (read-unsigned-n worker-number-length))))

(define-communicator-argument sequence-number (number)
  (write-unsigned-n sequence-number-length number)
  (read-unsigned-n sequence-number-length))

(define-communicator-argument random-seed (number)
  (write-unsigned-n random-seed-length number)
  (read-unsigned-n random-seed-length))

;;All zeroes represents no host
(defparameter no-host
  (make-array 4 :element-type '(unsigned-byte 8) :initial-element 0))

;;IP address, as defined by bsd-sockets, or NIL for no host
(define-communicator-argument host (address)
  (progn (unless address (setq address no-host))
	 (dotimes (index 4) (send-write-byte (aref address index))))
  (let ((address (make-array 4 :element-type '(unsigned-byte 8))))
    (dotimes (index 4) (setf (aref address index) (receive-read-byte (aref address index))))
    (unless (equalp address no-host)	;NIL if no-host
      address)))

;;Variable-length list of hosts
(define-communicator-argument host-list (hosts)
  (progn (send-write-byte (length hosts)) ;Send count
	 (dolist (host hosts)
	   (send-argument-host host)))
  (loop repeat (receive-read-byte)	;How many follow
	collect (receive-argument-host)))

(defun inet-address-string (address)
  (format nil "~A.~A.~A.~A" (aref address 0) (aref address 1) (aref address 2) (aref address 3)))


;;Worker-Manager protocol
;;Senders are called with the first argument being the buffer in which to write
(define-communicators peer-communicators peer-send
  ;;Requests from worker to manager
  (worker-ready (worker worker) (sequence sequence-number)) ;Initial message from worker
  ;;Job complete.  Bits set in flags tell whether there was any output to corresponding successors.
  (job-done (worker worker) (sequence sequence-number) (job job-number) (random random-seed) (flags byte))
  (job-failed (worker worker) (sequence sequence-number) (job job-number) (random random-seed))
  ;;Replies from manager to worker
  ;;Run job.
  (do-job (sequence sequence-number) (job job-number) (random random-seed) ;Seed or 0 to generate it
	  (predecessor-hosts host-list))
  (worker-exit)
  (worker-sleep (sequence sequence-number)) ;After this the next message comes from the manager
  ;;Requests from manager to worker
  (wake-exit)
  (wake-job (sequence sequence-number) (job job-number) (random random-seed) (predecessor-hosts host-list))
  ;;Replies from worker to manager
  (acknowledge-wake-exit (worker worker))
  (acknowledge-wake-job (worker worker) (sequence-number sequence-number))
  (worker-terminated (worker worker))
  (acknowledge-terminated))

;;Send one message.  *SEND-DESTINATION* has already been bound to the buffer.
(defun peer-send (code &rest arguments)
  (setf (fill-pointer *send-destination*) 0) ;Start writing at beginning of buffer
  (let ((communicator (aref peer-communicators code)))
    (send-write-byte code)		;Say what message this is
    (loop for argument in arguments
	  for argument-data in (communicator-argument-data communicator)
	  do (funcall (communicator-argument-send argument-data) argument)))
  (peer-send-it))

(define-timer :send-it)

;;Actually send message in *SEND-DESTINATION* to address in *PEER-ADDRESS*
;;This can be called again for retransmission
(defun peer-send-it ()
;;  (cond ((zerop (random 3))		;Simulate loss
;;	 (format t "Packet dropped~%"))
;;	(t
  (account-time :send-it
    (socket-send *peer-socket* *send-destination* (length *send-destination*) :address *peer-address*))
;;	 ))
	)


;;Retransmit last message to worker or manager
(defun peer-retransmit (worker)
  (let ((*send-destination* (worker-buffer worker)) ;Re-create state
	(*peer-address* (worker-address worker)))
    (peer-send-it)))

(define-timer :select)

;;Receive and process a message, or time out.  True on success, NIL on timeout.
;;TIMEOUT in internal time units
(defun peer-receive (&optional timeout)
  (when (and timeout (minusp timeout))
    (error "Timeout must be positive"))
  (multiple-value-bind (seconds remainder) (and timeout (truncate timeout internal-time-units-per-second))
    (let* ((microseconds (and remainder (* remainder (truncate 1000000 internal-time-units-per-second))))
	   (fd (socket-file-descriptor *peer-socket*)))
      ;;Using select here not only allows us to have a timeout but also protects us from problems when
      ;;you type ^C inside recvfrom
      (with-alien ((fds (struct sb-unix:fd-set)))
	(sb-unix:fd-zero fds)
	(sb-unix:fd-set fd fds)
	(let ((nset
	       (account-time :select
		 (sb-unix:unix-fast-select (1+ fd) (addr fds) nil nil seconds microseconds)))) ;Wait for packet
	  (and nset			;If interrupted, it can return NIL
	       (plusp nset)		;If none, return NIL
	       (peer-receive-1)))))))

(define-timer :receive-1)

;;Actually receive message when we know that there is one available.
(defun peer-receive-1 ()
 (account-time :receive-1
  (setf (fill-pointer *peer-buffer*)  (array-dimension *peer-buffer* 0)) ;Set to max length for reception
  (multiple-value-bind (buffer length source-ip source-port) ;Get message and bind variables to length and message source
      (without-interrupts ;Avoid reentrance to foreign functions and such.  Should be safe because select said available
       (socket-receive *peer-socket* *peer-buffer* (length *peer-buffer*)))
    (setf (fill-pointer buffer) length)	;Set to actual length received
    (let* ((*receive-source* buffer)	;Bind variables for reception & reply
	   (*peer-address* (list source-ip source-port)))
      (peer-receive-2) ;Decode arguments and call receiver
      ;;This function can be reentered here if the manager seems to be done
      t)			;Success
    )))

(define-timer :receive-2)

;;Accumulate arguments and dispatch incoming message
(defun peer-receive-2 ()
  (account-time :receive-2
    (setq *receive-index* 0)		;Get ready to read message
    (let* ((code (receive-read-byte))	;Only one message per UDP.  Read code
	   (communicator (aref peer-communicators code))) ;Look up message
      (loop for argument-data in (communicator-argument-data communicator)
	    collect (funcall (communicator-argument-receive argument-data)) into arguments
	    finally (apply (communicator-receive-function communicator) arguments)))))

;;Sequence number checking.
;;If this is a new message from the peer, increment the expected sequence number and return true
;;If MASTER-P is set, then accept also a message with one higher number.  This happens when the reply to
;;our previous request was lost and we should be the slave but don't know it yet.  We ignore the lost reply and
;;set up to process the current request.
;;If MASTER-P is not set, then we should be receiving requests.  If we receive a duplicate one, we retransmit
;;our last message.
;;The default MASTER-P is for manager calls
(defun check-sequence (worker sequence &optional (master-p (member (worker-status worker) '(:starting :exiting))))
  (let ((next (worker-sequence worker))) ;This is the message we expect
    (cond ((= sequence next)		;Got it
	   (incf (worker-sequence worker)) ;Advance to next
	   t)				;Process message
	  ((and master-p		;We are waiting for him to acknowledge our request?
		(= sequence (1+ next)))	;Acknowledgment lost, but here is his request.
	   (setf (worker-sequence worker) (1+ sequence)) ;OK.  Set sequence number for that state
	   t)				;Processing message
	  ((and (not master-p)		;We didn't send him a request, so he should be sending us one
		(= sequence (1- next))) ;Duplicate of last message?
	   (peer-retransmit worker)	;Retransmit to him
	   nil)				;But don't do anything else
	  ((< sequence next) nil)	;Earlier duplicate message: ignore
	  (t (if (member (worker-status worker) '(nil :not-started)) ;Old packet after termination?
		 (format t "~&Ignoring packet for not-started (terminated) ~A~%" worker)
	       (error "Invalid future sequence number ~D > ~D" sequence (worker-sequence worker)))))))




(define-alien-type statfs
  (struct statfs
	  (type long)
	  (bsize long)			;Block size
	  (blocks long)			;Total blocks
	  (bfree long)			;Total free blocks
	  (bavail long)			;Free blocks to non-superuser
	  (files long)			;Total inodes
	  (ffree long)			;Free inodes
	  (fsid long)
	  (namelen long)			;Max filename length
	  (frsize long)
	  (spare (array long 4))
	  ))

(declaim (inline statfs))
(define-alien-routine statfs int
  (path c-string)
  (data (* statfs)))
	  
;;Return number of bytes free on the file system with the given directory
(defun free-disk-space (file)
  (when (pathnamep file) (setq file (namestring file)))
  (with-alien ((data statfs))
    (unless (zerop (statfs file (addr data))) ;Get information about file system
      (error "Can't get disk space for ~A" file))
    (* (slot data 'bsize) (slot data 'bavail))))
    
(defun free-disk-gb (file)
  (/ (free-disk-space file) (expt 2.0 30)))
