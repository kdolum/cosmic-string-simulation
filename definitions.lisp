;;General definitions
(in-package "CL-USER")

;;Package manipulation
(require :sb-sprof)
(unuse-package "SB-PROFILE")
(use-package "SB-SPROF")		;statistical profiling

(require :sb-bsd-sockets)
(use-package "SB-BSD-SOCKETS")
(defparameter *inet-address-any* (vector 0 0 0 0)) ;In 1.4 we imported this from SB-BSD-SOCKETS, but in 2.1 it is not there

(use-package "SB-SYS")			;for SAP-REF-...

(use-package "SB-THREAD")

(eval-when (:compile-toplevel :load-toplevel :execute)
  (use-package "SB-MOP"))			;For clear-defstruct

;;For unknown reasons, this does not work when put in gc.lisp, even with eval-when
(proclaim '(disable-package-locks sb-kernel:dynamic-usage sb-kernel::get-bytes-consed))

;;Configuration
;;Figure out which site we're using
(eval-when (:compile-toplevel :load-toplevel :execute)
(defparameter server
  (let ((name (machine-instance)))
    (cond ((string-equal "cosmos.phy.tufts.edu" name) :cosmos)
	  ((string-equal "cosmos" name) :cosmos)
	  ((search ".uwm.edu" name) :uwm)
	  ((probe-file "/cluster/tufts") ;Tufts?
	   (cond ((posix-getenv "SLURM_JOB_ID") ;Running under slurm
		  :tufts)
		 ((posix-getenv "LSB_JOBID") ;Running under lsf
		  :tufts-lsf)
		 ((string-equal name "login" :end1 (min (length name) 5)) ;head node
		  :tufts)
		 (t (error "Running on unknown tufts machine ~S" name))))
	  (t (error "Running on unknown server ~S" name))))))

;;Find lisp executable.  Run whatever version we are running
(defvar lisp-program
  (ecase server
    (:tufts (format nil "/cluster/tufts/strings/sbcl-~A/run-sbcl.sh" (lisp-implementation-version)))))

;;Macro definitions
;;Once-only macro that binds the extra variables only if they are needed.
;;Gratuitous bindings would be optimized away by the compiler, but they interfere with the optimizers in vector.lisp
;;Modified from Peter Norvig, Paradigms of AI Programming.
(defmacro once-only (variables &rest body)
  "Returns the code built by BODY.  If any of VARIABLES
  might have side effects, they are evaluated once and stored
  in temporary variables that are then passed to BODY."
  (assert (every #'symbolp variables))
  (let ((temps nil))
    (dotimes (i (length variables)) (push (gensym) temps))
    `(if (every #'atom (list ,@variables)) ;Code going in your macro: see if variables are all atoms
	 (progn ,@body)			   ;Yes: just return body
       (list 'let			   ;No: generate code to bind temporaries to the forms
	     ,`(list ,@(mapcar #'(lambda (tmp var)
				   `(list ',tmp ,var))
			       temps variables))
	     (let ,(mapcar #'(lambda (var tmp) `(,var ',tmp)) ;Bind forms to temporary symbols and evaluate body
			   variables temps)
	       ,@body)))))

;;Define a constant object which requires consing.  The new one wouldn't be EQ to the old,
;;giving an error.  So if it is bound already, we use the old value.
(defmacro define-constant (name value)
  `(defconstant ,name (if (boundp ',name) (symbol-value ',name) ,value)))

;;Destructively delete object from sequence, as recommended by Peter Norvig and Kent Pitman
;;in "Tutorial on Good Lisp Programming Style"
(defmacro deletef (item sequence &rest keys &environment environment)
  (multiple-value-bind (temps vals stores store-form access-form)
      (get-setf-expansion sequence environment)
    (assert (= (length stores) 1))
    (let ((item-var (gensym "ITEM")))
      `(let* ((,item-var ,item)
	       ,@(mapcar #'list temps vals)
	       (,(first stores)
		(delete ,item-var ,access-form ,@keys)))
	  ,store-form))))

;;Data type for coordinates
(deftype coordinate () 'double-float)

(defconstant zero-coordinate 0.0)

(defconstant euler-constant 0.577215664901532)

;;Largest size of a diamond in any given dimension.
(defvar diamond-span)

;;Intersections are not permitted within this distance of an edge, because that would generate
;;a diamond that is too narrow to evolve accurately.
(defparameter minimum-diamond-width 1e-10)
;;Intersections are not permitted within a distance of an edge which is a scaled value of this
;;parameter. The scaling depends on the cosmological era and difference between the current time
;;and *overall-end-time*.
(defparameter minimum-ending-diamond-width 1e-8)
(defvar *overall-end-time* nil) ;Given a value in SIMULATE-1 based on what variables are set

;;Fudge factors for floating point comparisons

;;Fudge factor for comparing coordinates in a single simulation volume.  The volume could not be larger than about 100,
;;so the low order bit is about 1e-14.
(defparameter fudge-coordinates 1e-12)

;;Fudge factor for global coordinates.  These could perhaps be as large as 10000.
(defparameter fudge-global-coordinates 1e-10)

;;Fudge factor for ijkl coordinates involving a single cube.  The coordinates run roughly from 0 to 1.
;;numbers smaller by that factor. 
(defparameter fudge-ijkl 1e-14)

;;Fudge factor for global ijkl coordinates.  There could be perhaps 1000 layers of cubes
(defparameter fudge-ijkl 1e-11)

(declaim (double-float diamond-span minimum-diamond-width minimum-ending-diamond-width
		       fudge-coordinates fudge-global-coordinates fudge-ijkl fudge-global-ijkl	))

(declaim (inline fudge=))
;;Equal within fudge factor
(defun fudge= (x y fudge-factor)
  (<= (abs (- x y)) fudge-factor))

(defun fudge=global (x y)
  (fudge= x y fudge-global-coordinates))

;;Equal within fudge factor modulo divisor
;;We can't just use (fudge= (mod ...)) because, e.g., (mod -0.001 1.0) = 0.999
(defun fudge=mod (x y divisor fudge-factor)
  (multiple-value-bind (quotient remainder) (round (- x y) divisor)
    (declare (ignore quotient))
    (<= (abs remainder) fudge-factor)))

;;Representation specific versions of the function FLOAT
(declaim (inline double-float single-float))
(defun double-float (x)
  (coerce x 'double-float))
(defun single-float (x)
  (coerce x 'single-float))

;;Quickly take the square root of a positive real number.  We have to tell
;;the compiler that the argument won't be negative.
(defmacro real-sqrt (x)
  `(sqrt
    (the (double-float 0.0 *)		;positive double-float
      ,x)))

;;This is harder than SQRT, because the compiler doesn't know that when the argument is positive
;;the result will not be complex.
(defmacro real-log (&rest x)
  `(locally (declare (optimize (safety 0))) ;Don't check type in THE
      (the double-float (log ,@x))))
;;There was previously a without-compiler-notes here, commented "Avoid mysterious "note: doing
;;SAP to pointer coercion" having to do with alien call for log", but as of 1.4.0 it doesn't seem needed

;;A double float that is no larger than a fixnum could accommodate.  The limits are chosen
;;so that floating point errors (in type checking) or floor or ceiling
;;can never take the number outside the fixnum range
(deftype fixnum-size-float ()
  (list 'double-float
	(* (- 1.0 (* double-float-epsilon 10)) (double-float most-negative-fixnum))
	(* (- 1.0 (* double-float-epsilon 10)) (double-float most-positive-fixnum))))

;;one-argument, one-value floor without worrying about the possibility that the result is too large for a fixnum
(defmacro fixnum-floor (x)
  `(prog1 (floor (the fixnum-size-float ,x))))

(defmacro fixnum-ceiling (x)
  `(prog1 (ceiling (the fixnum-size-float ,x))))

;;Round float, giving fixnum
(defmacro fixnum-round (x)
  `(fixnum-floor (+ ,x 0.5)))

;;Avoid consing to call out for ffloor by using above
(defmacro fast-ffloor (x)
  `(float (fixnum-floor ,x) 0.0))

;;Determine endianness
(eval-when (:compile-toplevel :load-toplevel :execute)
  (with-alien ((char-pointer (array char 4)))
    (declare (optimize (speed 0)))	;Suppress optimization warnings
    (let ((int-pointer (cast char-pointer (array (unsigned 32) 1))))
      (setf (deref int-pointer 0) #x01020304)
      (let ((little (case (deref char-pointer 0)
		      (4 t)
		      (1 nil)
		      (t (error "Can't determine endianness")))))
	(if little
	    (if (member :big-endian *features*)
		(error "Wrong endianness feature already present")
	      (pushnew :little-endian *features*))
	  (if (member :little-endian *features*)
		(error "Wrong endianness feature already present")
	    (pushnew :big-endian *features*))
	    )))))

;;Index wrapping in arrays
(declaim (inline next-index-wrapping previous-index-wrapping
		 next-array-index-wrapping previous-array-index-wrapping
		 east-array-index-wrapping west-array-index-wrapping))

;;Next index, wrapping at specified limit
(defun next-index-wrapping (index limit)
  (mod (1+ index) limit))

;;Previous index
(defun previous-index-wrapping (index limit)
  (mod (1- index) limit))

;;Get length from array and wrap
(defun next-array-index-wrapping (index array)
  (mod (1+ index) (length array)))

(defun previous-array-index-wrapping (index array)
  (mod (1- index) (length array)))

;;Mirror-image compatible versions
(defun east-array-index-wrapping (index array)
  (next-array-index-wrapping index array))
(defun west-array-index-wrapping (index array)
  (previous-array-index-wrapping index array))


;;Generate symmetrical code for right and left.  The body will be expanded twice.  The second time
;;we will make replacements WEST<->EAST, RIGHT<->LEFT, NE,->NW, SE<->SW, E<->W, A<->B when these strings appear as
;;complete atoms or delimited by hyphens.
;;This does not work for things inside commas (inside backquotes), because the sbcl reader does not render these
;;as lists.
;;A piece of code (NOMIRROR ...) will be replaced by ... without mirroring.  If ... is multiple forms they will be
;;spliced into the lists this form is part of.
(defmacro mirror-images (&body body)
  `(progn
     ,@(mirror-image body nil)
     ,@(mirror-image body t)))

;;Expand let bindings in the same manner.  Bindings which have anything to mirror appear twice,
;;others are left alone
(defmacro mirror-image-let (bindings &body body)
  `(let ,(mirror-let-bindings bindings)
     ,@body))
(defmacro mirror-image-let* (bindings &body body)
  `(let* ,(mirror-let-bindings bindings)
     ,@body))

;;code needed at compile time for mirror-image-destruct definitions
(eval-when (:compile-toplevel :load-toplevel :execute)

;;Bind variables with mirror images where needed.  This needs to handle NOMIRROR specially so
;;that it can be wrapped around several bindings.  We can't just pass the binding list to MIRROR-IMAGES
;;because then they would come out in the wrong order.
(defparameter mirror-images
  (loop for (from to) on '("WEST" "EAST" "RIGHT" "LEFT" "NE" "NW" "SE" "SW" "E" "W" "A" "B" "A1" "B1" "A2" "B2") by #'cddr
	collect (cons from to)		;substitute both ways
	collect (cons to from)))

;;Turn thing into its mirror image, or if reflect-p is nil, just strip nomirror
(defun mirror-image (thing reflect-p)
  (typecase thing
    (null nil)				;Trivial case of symbol NIL
    (cons
     (when (eq (car thing) 'nomirror)
       (error "NOMIRROR should not be at top-level in something to be mirrored"))
     (loop for form in thing
	   append
	   (if (and (consp form) (eq (car form) 'nomirror)) ;form is (NOMIRROR ...)
	       (cdr form)				    ;splice in ...
	     (list (mirror-image form reflect-p)))))	    ;accumulate reflected forms
    (symbol 
     (if reflect-p
	 (let* ((old-name (symbol-name thing))
		(new-name (mirror-image-name old-name)))
	   (if (eq new-name old-name)	;no changes
	       thing
	     (intern new-name (symbol-package thing)))) ;make new symbol
       thing))
    (t thing)))

;;Scan string replacing chiral words with their mirror images
(defun mirror-image-name (name)
  (loop for position = 0 then (1+ next) ;find components between hyphens
	for next = (position-if #'(lambda (char) (find char "-*")) ;Positiong of delimiter, or NIL
				name :start position)
	for cons = (assoc (subseq name position next) mirror-images
			  :test #'string-equal)
	when cons			;Found something to change
	  do (setq name (concatenate 'string ;make change in string
				   (subseq name 0 position)
				   (cdr cons)
				   (and next (subseq name next)))
		   next (and next (+ position (length (cdr cons))))) ;adjust to keep at hyphen
		   
	while next)
  name							   ;return possibly modified name
  )

(defun mirror-let-bindings (bindings)
  (loop for binding in bindings
	append (mirror-let-binding binding)))

;;Mirror one let binding.  Or if this is (NORMIRROR ... bindings ...) just strip nomirror
(defun mirror-let-binding (binding)
  (if (and (consp binding) (eq (car binding) 'nomirror)) (cdr binding) ;not mirrored: just collect it
    (let ((original (mirror-image binding nil))
	  (mirror (mirror-image binding t)))
      (if (equal original mirror)	;no change on mirroring
	  (list original)		;don't bind twice
	(list original mirror)))))

;;Mirror keywords for defstruct construction.  A lot like mirror-let-bindings
(defun mirror-keyword-list (list)
  (loop for keyword = (pop list)
	while keyword
	append (if (and (consp keyword) (eq (car keyword) 'nomirror))
		   (cdr keyword)
		 (let* ((pair (list keyword (pop list))) ;(:key value)
			(original (mirror-image pair nil))
			(mirror (mirror-image pair t)))
		   (if (equal original mirror) ;no change on mirroring
		       original		       ;just use original pair
		     (append original mirror))))))
)					;eval-when

;;additional slots mirroring those that are given, and a constructor called MIRROR-IMAGE-(real constructor)
;;that mirror images the slot names and their initial values
(defmacro mirror-image-defstruct (name-and-options &body variables)
  (let* ((name (if (consp name-and-options) (car name-and-options) name-and-options))
	 (constructor
	  (let ((entry (and (consp name-and-options) (assoc :constructor (cdr name-and-options)))))
	    (if entry (second entry)					    ;If given use it
	      (intern (format nil "MAKE-~A" name) (symbol-package name)))))) ;Default like defstruct
    `(progn
       (defstruct ,name-and-options ,@(mirror-let-bindings variables))
       ,@(and constructor
	      (list
	       `(defmacro ,(intern (format nil "MIRROR-IMAGES-~A" constructor) (symbol-package constructor))
		    (&rest slots)
		  `(,',constructor ,@(mirror-keyword-list slots))))))))


;;Dispatch on type flag :A/:EAST or :B/:WEST.  If the first, do the body as given.  If the second, the mirror image.
(defmacro choose-mirror-image (type &body a-body)
  (once-only (type)
    `(block choose-mirror-image
       (mirror-images
	(when (member ,type '(:a :east))
	  (return-from choose-mirror-image (progn ,@a-body))))
       (error "Type was ~S not :A or :B" ,type))))

;;CMUCL has this but not SBCL
(defun required-argument ()
  (error "A required defstruct slot was not supplied."))    

(defmacro without-compiler-notes (&body body)
  `(locally (declare (muffle-conditions compiler-note)) ;Avoid notes about not optimizing
     ,@body))

(defmacro with-compiler-notes (&body body)
  `(locally (declare (unmuffle-conditions compiler-note)) ;Re-enable notes
     ,@body))

;;7/14/21:  Recent versions of sbcl have the code that conses floats to return associated with the outermost
;;part of the function instead of with the valuation of the object to be returned.  Put this macro around
;;your defun to suppress the return consing note.  Probably it suppresses notes from function entry too.
(defmacro without-compiler-notes-return (defun)
  (unless (eq (car defun) 'defun)	;Do we understand what we're doing?
    (error "~S should have been wrapped around a defun" 'without-compiler-notes-return))
  (destructuring-bind (name arglist . body) (cdr defun)
    (let ((first (list name arglist)))
      (when (stringp (first body))	;documentation string?
	(setq first (append first (list (pop body))))) ;Move it to first
      `(defun ,@first
	   (declare (muffle-conditions compiler-note)) ;Insert this declaration before body
	 (with-compiler-notes ,@body)))))				       ;Turn declarations back on locally

(defconstant vv-index-size 16)
(deftype vv-index () `(unsigned-byte ,vv-index-size))	;Index into initial VV lattice.

(deftype site () `(unsigned-byte ,(* 3 vv-index-size)))

;;Axis number.  Cannot exceed maximum rank 3.  However,
;;loop for axis of-type axis from 1 to rank increments variable one extra time,
;;so we have to allow for that.
(deftype axis () `(integer 1 4))

;;Direction number.  Like axis, but can be negative
(deftype direction () `(integer -4 4))



;;The function SIMULATE takes many arguments, most of which serve only to bind special variables,
;;so we generate the argument list from a set of declarations

(eval-when (:compile-toplevel :load-toplevel :execute)

;;Code and variables for SIMULATE system
(defvar simulate-keywords nil)
(defvar simulate-coercions nil)

;;Add code to coerce variable to given type to SIMULATE-COERCIONS
(defun add-simulate-coercion-code (variable type)
  (when type
    (pushnew `(when ,variable
		(setq ,variable (coerce ,variable ',type)))
	     simulate-coercions :test #'equal)))

)					;eval-when


;;Define an argument to simulate with, optionally, a default value and a type to coerce the argument to
;;if it is not NIL.
(defmacro define-simulate-argument (variable &optional default type)
  `(eval-when (:compile-toplevel :load-toplevel :execute)
     (pushnew '(,variable ,default) simulate-keywords :test #'equal)
     (add-simulate-coercion-code ',variable ',type)))

;;Create special VARIABLE with DEFVAR, and set its default value to DEFAULT
;;The function SIMULATE will then take this variable as an argument, with the default value being whatever 
;;the variable is already bound to.
;;If TYPE is set, the value will be coerced to that type, unless it is NIL.
(defmacro define-simulate-variable (variable default &optional type)
  (let* ((name (symbol-name variable))
	 (keyword (intern (if (char= (char name 0) (char name (1- (length name))) #\*)
			     (subseq name 1 (1- (length name))) ;Strip asterisks from variable name to make keyword
			   name)
			 "KEYWORD")))
    `(progn
       (defvar ,variable ,default)
       (eval-when (:compile-toplevel :load-toplevel :execute)
	 (pushnew '((,keyword ,variable) ,variable) ;Bind to itself by default
		  simulate-keywords :test #'equal)
	 (add-simulate-coercion-code ',variable ',type)))))

;;Code for coercing simulate arguments
(defmacro simulate-coerce ()
  `(progn ,@simulate-coercions))




;;Group modifications during which our data structures might be inconsistent.
;;This suppresses calls to (check-data-structures) in the body and calls it at the end.
(define-simulate-variable *check-data-structures* nil)	;If set, call check-data-structures periodically

;;This doesn't use use unwind-protect, because we don't want to check when there was an error.
(defmacro with-modification-group (&body body)
  `(prog1
       (let ((*check-data-structures* nil))
	 ,@body)
     (when *check-data-structures*	;Unless a suppressed or inside outer group
       (check-data-structures))	  ;Everything should be consistent now
     ))


;;;Global variables

(defvar *calendar*)			;Main calendar of events

;;Hash table giving a set of diamonds to pass to successors
(defvar *final-diamonds*)

(define-simulate-variable *job-number* nil fixnum) ;Number of our job.  Set by setup-geometry or inside worker code

(defvar *size* nil)			;The periodicity distance of this processor, if it were the only one.
					;The temporal cube length is larger by sqrt{2/3}

(defvar *total-size* nil) ;Total size (periodicity distance) of this simulation.  NIL if infinite volume.
(declaim (fixnum *split-factor*))
(defvar *split-factor* 1)		;Number of pieces that the lattice is split into along each of the 4 dimensions.
(defvar *global-location* nil)		;4vector giving the position of our starting corner
(defvar *global-ijkl* nil)		;ijkl coordinates of our starting vertex

(define-simulate-variable *random-seed* nil) ;If set, use this instead of random-seed files
(define-simulate-variable *time-offset* 0.0 double-float) ;Global time of ijkl origin


;;If set, print dots, etc.
(define-simulate-variable *print-progress* nil)

;;Time part or all of the run:  :setup, :input, :evolution, :output, :log, :overall, or a list of these
(defvar *report-time* nil)		;The corresponding simulate argument is :time.
(defmacro maybe-time (when &body body)
  `(if (or (eq *report-time* t) (eq *report-time* ,when)
	   (and (listp *report-time*) (member ,when *report-time*)))
       (time (progn ,@body))
     (progn ,@body)))

(define-simulate-variable *era* :flat)	; :FLAT, :RADIATION, or :MATTER

(defvar *initializing* nil)		;Flag to say that we are doing initial conditions
;;Local time start of simulation.  The later of the initial condition time or the beginning of the volume
(defvar *initial-time* 0.0)

;;Bound to list of predecessor hosts, or NIL just to use local files.
(defvar *predecessor-hosts* nil)
;;If set, output will set bits to tell whether successor files were nonempty.
(defvar *successor-file-flags* nil)

;;Normally we do not delete loops from the simulation if they are in self-intersecting trajectories but have not
;;actually reconnected because of low *intercommutation-probability.  If this variable is set, we treat such
;;loops (and those that have failed to rejoin for similar reasons) as non-self-intersecting and delete them
;;at the usual time.  Probably you should also set *loop-preservation-threshold* if you use this, and
;;you should not also be looking at the loop distribution.
(define-simulate-variable *delete-unlucky-loops* nil)

;;Loops larger than this fraction of the horizon size will be retained in case they might rejoin.
;;0 = always retain; never delete.
;;NIL = always delete
(define-simulate-variable *loop-preservation-threshold* nil double-float)

;;Time to start preserving loops so that they are dumped.
(define-simulate-variable *loop-preservation-dump-start* 0.0 double-float)

;;Smallest scaling size loop we keep around until the next dump.  A loop is considered dumped when any
;;diamond on it has countup >= 2 at the dump time.
;;NIL = don't do it
(define-simulate-variable *loop-preservation-dump-x* nil double-float)

(defvar *loops-found* nil)		;Array to store deleted loops, or NIL if not keeping
(defvar *loops-output* nil)		;Stream to output deleted loops

(defvar *bh-loops-found* nil) ;Array to store deleted bh loops, or NIL if not keeping
(defvar *bh-loops-output* nil) ;Stream to output deleted bh loops

(define-simulate-variable *log-loop-positions* nil) ;If set, write positions of loops as well
(define-simulate-variable *log-loop-velocities* nil) ;Write 3vector velocity

(defvar *length-output* nil)		;Stream to store length of string
(defvar *information-flow-output* nil)	;Stream to store information flow distance

;;Prevent or detect reentrant use of the diamond-processed bit
(defvar *using-diamond-processed* nil)

(defmacro using-diamond-processed (&body body)
  `(progn
     (when *using-diamond-processed*
       (error "reentrance to ~S" 'using-diamond-processed))
     (let ((*using-diamond-processed* t))
       ,@body)))

(defvar *send-first-diamond*)	      ;First diamond in current string
(defvar *send-previous-diamond*)	;Last diamond sent
(defvar *send-count-diamonds* (make-array 5)) ;Count number of diamonds sent.  If none, abort file.

(defvar *receive-start-junction*)     ;Junction at start of string
(defvar *receive-first-diamond*)	;First diamond in current string
(defvar *receive-previous-diamond*)	;Last diamond received
(defvar *receive-count-diamonds*)     ;Count number of diamonds received

;;Destinations are numbers giving the direction in which the information is to be output:
;; i = 0, j = 1, k = 2, l = 3, dump = 4
(defconstant dump-destination 4)

(defvar *dump-time* nil)		;Bound to global time of current dump, or NIL if not dumping,
					;or :unknown if we're reading a dump and haven't yet learned the time

(defvar *reading-dumps* nil)		;NIL or global time for which dumps are being or have been read
					;or :unknown if we're reading a dump and haven't yet learned the time
(defvar *read-and-store* nil) ;T if reading dumps or pre-reading from predecessor.  Just store diamond in *read-diamonds*

(defvar *longest-edge* nil)		;Longest p or q in initial data, or NIL.

;;Hash table of diamonds read from dumps
(defvar *read-diamonds*)

(defvar monster-length-scale)	 ;Conformal length scale to determine monsters using largest known diamond as reference.

;;How to handle datafiles: NIL = store globally in batch-root-directory 
;;If set, use local-root-directory
;; :RSYNC -- copy files with rsync server
;; :NFS -- NFS-mount scratch on source node

;;New tufts cluster does not support scratch-over-NFS.
;;(defparameter *local-data-files* (if (eq server :tufts) :rsync :nfs))
(defparameter *local-data-files* (if (eq server :tufts) nil :nfs))

;;On some systems is possible to run a backrub process which persists after all batch jobs have exited.
;;On others (e.g., Tufts), you cannot do that, so we have to resort to running the server only one worker is alive
;;and consequently we must keep the worker alive until the run is finished.
(defparameter *permanent-rsync-server* (not (eq server :tufts)))

;;The top-level directory (without worker-N) in which we're working.
;;If *LOCAL-DATA-FILES* is :NFS, this is where we look for input files.
;;If NIL, not inside worker
(defvar *worker-directory* nil)

;;Maintain timing information for a task.  Variable should be set to NIL, in which case no timing
;;is collected, or (real-time . run-time).
(defvar *timers* nil)

(defun timer-variable (key)			;Convert timer keyword to symbol
  (intern (format nil "*~A-TIMER*" key)))

(defmacro define-timer (key)
  `(progn
     (defvar ,(timer-variable key) nil)
     (pushnew (cons ,key ',(timer-variable key)) *timers* :key #'car)))

(defmacro account-time (key &body body)
  (let ((var (timer-variable key)))
    `(if ,var
	 ,(let ((start-real (gensym))
		(start-run (gensym)))
	    `(let ((,start-real (get-internal-real-time))
		   (,start-run (get-internal-run-time)))
	       (multiple-value-prog1 (progn ,@body)
		 (incf (car ,var) (- (get-internal-real-time) ,start-real))
		 (incf (cdr ,var) (- (get-internal-run-time) ,start-run)))))
       (progn ,@body))))


;;Information about cusps.  Not all uses use all slots.  The parameter used in x-a and in differentiation is
;;x = 0..1 for smooth strings, sigma = 0..L otherwise.
(mirror-image-defstruct cusp-info
  direction				;Direction of the cusp (a'=b')
  x-a					;Parameter value in a where cusp occurs
  i-a					;Index of a-1, the first side of the cusp.
  a-1					;hat before segment that crosses
  a-2					;hat after segment that crosses
  a-3					;next hat after a-2 for check-cusp-form-paper
  a-dsigma-1				;Amount of parameter spent at a-1
  a-dsigma-2				;Amount of parameter spent at a-2
  a-dsigma-3				;Amount of parameter spent at a-3 for check-cusp-form-paper
  sigma-a-1-to-cusp			;We split the average of a-dsigma-1 and a-dsigma-2 into
  sigma-cusp-to-a-2			;before and after the cusp by angular position
  ;;Second derivatives at the cusp with respect to the parameter, times L (but not divided by 2 pi!)
  a-pp
)

;;Structure to store information about run in file
(defstruct run-info			;Information about a run
  era
  total-size				;Total size of simulation
  split-factor				;Number of sub-cubes in each direction
  time-offset				;Global time of first cube corner
  start-time				;conformal time at start
  end					;conformal time of (desired) end
  dump-times				;Time of real dumps
  length-times				;Time length was dumped
  bh-size                               ;size of the bhs
  bh-number                             ;number of bhs
  bh-times                              ;creation time of bhs
  loop-preservation-threshold
  jobs					;Number of jobs in run
  ;;Old slots to read old runs
  dump-start				;Global time of first dump
  dump-interval
  dump-only-length			;Only lengths were dumped
  )

;;An event that happens at a certain time, rarely enough that efficiency is not important
(defstruct timed-event
  time
  (report-string "")			;What to print when we find one
  function				;Function called with event as argument
  )

;;Like with-open-file, but only if condition is met.  We used to with-open-stream, but now that tries to close NIL
(defmacro with-maybe-open-file ((stream condition &rest arguments) &body body)
  (let ((aborted (gensym)))
    `(let ((,stream (and ,condition (open ,@arguments)))
	   (,aborted t))
       (unwind-protect
	   (multiple-value-prog1 (progn ,@body)
	     (setq ,aborted nil))
	 (when ,stream (close ,stream :ABORT ,aborted))))))

(defvar *previous-umask*)

(defparameter batch-root-directory
  (ecase server
    (:cosmos "/strings/")
    ((:tufts :tufts-lsf) "/cluster/tufts/strings/")
    (:uwm "/home/kdo/strings/")))

;;Allow other members of the simulation group to modify files
(defmacro with-group-write-access (&body body)
  (case server
    ((:tufts :tufts-lsf :cosmos)
     `(let (*previous-umask*)		;Rebind this variable so macro can be used recursively
	(setq *previous-umask* (umask 2)) ;Allow group access to files written by lisp.  Remember previous mask
	(unwind-protect (progn ,@body)
	  (umask *previous-umask*))))
    (t 
     `(progn ,@body) 			;No effect elsewhere
    )))

(defvar *manager-jobs* nil)		;Number of jobs being managed.  NIL if no manager.  Used in run-info

(defvar *threads* nil)			;Number of threads to use for multithreaded computation in this process

;;Bound by communicators; used for debugging
(defvar *debug-receive-communicator* nil)	;Bound in receiving
(defvar *debug-receive-argument-number* nil)
(defvar *debug-send-communicator* nil)	;Bound in receiving
(defvar *debug-send-argument-number* nil)

;;Reset all slots in structure to their default values
;;This uses the meta-object protocol.  If something goes wrong with it, one could change the code to
;;just clear the slots by hand.
;;The structure definition needs to be available at compile time.
(defmacro clear-defstruct (type object)
  `(progn
     ,@(loop for slot in (class-slots (find-class type))
	     for name = (slot-definition-name slot)
	     for initform = (slot-definition-initform slot)
	     collect `(setf (slot-value ,object ',name) ,initform))))


;;Definitions for number of BH, creation time for BH and the size of BH, at this time all the BH                                                                                                           
;;has the same size
(defvar *bh-number* nil) ;Number of BH created                                                                                                                                                             
(defvar *bh-size* nil)   ; Radius of the BHs                                                                                                                                                               
(defvar *bh-start* nil) ; first detection of BH intersections    
(defvar *pointbh* nil) ;when set to t pointlike bhs are created 
(defvar *bh-probability* nil) ;probability for creating a BH in a given segement

;; Structures for BH analysis
;;black hole structure. It contains a center which is a 3vector and the radius of the BH
(defstruct blackhole
  center ;position of the blackhole
  size ;radius of the blackhole
)

;; BH-INTERSECTIONS                                                                                                                                                                                        
;;BH intersection structure. Information about the diamond in which a intersection takes place, the a and b and the spacetime position.
(defstruct bh-intersec
  diamond ;diamond in which the intersection is located
  bh ;blackhole which is intersecting string
  a ;a of the intersection point
  b ;b of the intersection point
  spacetime ;spacetime position of the intersection
  )



;;Variable and structure for the segment length dumping


(defvar *dumped-segments* nil)
(defvar *start-dumped-segments* (make-hash-table :test 'equalp))
(defvar *end-dumped-segments* (make-hash-table :test 'equalp))
(defvar *read-dump-diamonds* t)
(defvar *segment-start-junction* nil)
(defvar *receive-start-bh* nil)
(defvar *receive-end-bh* nil)

(defvar *receive-mybh* nil)


(defvar *all-bhs* nil)

(defvar *job-bhs* (make-hash-table :test 'equalp))
(defvar *my-bhs* (make-hash-table :test 'equalp))

(defstruct dump-segment
  dump-time
  start-junction
  start-bh
  end-junction
  end-bh
  length
  )
