;;;Handle processing.
(in-package "CL-USER")

(defstruct (handle
	    (:constructor make-handle (creator code))
	    (:print-object print-handle))
  creator				;Job number of creator
  code)					;his code

(defun print-handle (handle stream)
  (print-unreadable-object (handle stream)
   (format stream "HANDLE ~D on " (handle-code handle))
   (let ((creator (handle-creator handle)))
     (cond ((eq creator *job-number*)
	    (format stream "ME"))
	   (t
	    (format stream "job ~D" creator))))))

;;We use integer codes starting at 1
;;0 is reserved.
(defvar *last-handle-code*)

;;Allocate and returne a new handle in our space
(defun create-local-handle ()
  (incf *last-handle-code*)
  (debug-trace-handle (make-handle *job-number* *last-handle-code*) "created"))

(defvar *object-handle*)		;object -> handle structure
(defvar *handle-objects*)		;handle -> list of objects

;;Install or destroy correspondence between object and handle and return handle
(defsetf object-handle (object) (handle)
  `(progn
     (if ,handle
	 (connect-handle ,object ,handle)
       ,(let ((old (gensym)))
	  `(let ((,old (object-handle ,object nil)))
	     (when ,old (disconnect-handle ,object ,old)))))
     ,handle))

;;Get handle, if any, for given object.
;;If we created the object, the handle will be in our space.  If we received it
;;from the neighbor, the handle will be in the space of the object creator.
;;If no handle exists already, action depends on :if-does-non-exist
;;  If :error (the default), signal an error
;;  If :create, create a new local handle
;;  If NIL, return NIL
(defun object-handle (object &optional (if-does-not-exist :error))
  (or (gethash object *object-handle*)
      (ecase if-does-not-exist
	(:create	      ;Create, install, and return new local handle
	 (setf (object-handle object) (create-local-handle)))
	(:error (error "Object ~S has no handle" object))
	((nil) nil))))

;;Find the objects corresponding to a given handle.
;;If if-does-non-exist nil, returns NIL if not known,
;;otherwise gives an error.
(defun handle-objects (handle &optional (if-does-not-exist :error))
  (or (gethash handle *handle-objects*)
      (ecase if-does-not-exist
	(:error (error "Handle ~S is not known" handle))
	((nil) nil))))

;;Get object with this handle.  Error if more than 1
(defun handle-object (handle &optional (if-does-not-exist :error))
  (let ((objects (handle-objects handle if-does-not-exist)))
    (when (cdr objects)
      (error "More than one object has handle ~S" handle))
    (car objects)))

;;Establishes correspondence between an object and a handle
(defun connect-handle (object handle)
  (unless object (error "Can't give NIL a handle"))
  (unless handle (error "NIL is not a valid handle"))
  (when (gethash object *object-handle*)
      (error "~S already has a handle" object))
  (setf (gethash object *object-handle*) handle)
  (push object (gethash handle *handle-objects*)))

;;Destroys correspondence between an object and a handle
(defun disconnect-handle (object handle)
  (remhash object *object-handle*)
  (deletef object (gethash handle *handle-objects*)))

(defun handle= (handle1 handle2)
  (and (= (handle-creator handle1) (handle-creator handle2))
       (= (handle-code handle1) (handle-code handle2))))

;;Compare handles in a way which, while arbitrary, will always returns the same result on the same handles
(defun compare-handles (handle1 handle2)
  (or (< (handle-creator handle1) (handle-creator handle2))
      (and (= (handle-creator handle1) (handle-creator handle2))
	   (< (handle-code handle1) (handle-code handle2)))))
  
;;Initialize data structures
(defun initialize-handles ()
  (setq *last-handle-code* 0
	*object-handle* (make-hash-table :test #'eq)
	*handle-objects* (make-hash-table :test #'equalp)))
