(in-package "CL-USER")

(defstruct (resource
	    (:print-object print-resource))
  (name "anonymous")
  (constructor #'error :type function)	;Function of no arguments to create object
  (objects (make-array 16) :type simple-vector)
  (pointer 0 :type fixnum)	   ;Number of objects currently stored
  (max 1000000 :type fixnum)	;Maximum number of objects to remember
  (constructed 0)		;Number of objects created
  (allocated 0))		;Number of objects given out

(defun print-resource (resource stream)
  (print-unreadable-object (resource stream :identity t)
    (format stream "RESOURCE ~A storing ~D" (resource-name resource) (resource-pointer resource))))

;;Do or don't monitor.  Decide at compile time
(defmacro monitoring-resource-allocation (&body body)
;;  `(without-compiler-notes ,@body)	;Don't complain about optimization
  (declare (ignore body))
  )

;;Allocate an object from the resource.  Create one if not available.
(defun allocate (resource)
  (declare (type resource resource)
	   (optimize speed))
  (monitoring-resource-allocation (incf (resource-allocated resource)))
  (locally (declare (optimize (safety 0))) ;Check argument on entry, but otherwise don't check anything
    (let ((pointer (resource-pointer resource)))
      (cond ((plusp pointer)		;Have one available?
	     (prog1 
		 (aref (resource-objects resource) (decf pointer)) ;Give it out
	       (setf (resource-pointer resource) pointer) ;Store new pointer
	       (setf (aref (resource-objects resource) pointer) nil))) ;Clear slot so object not protected from GC
	    (t
	     (monitoring-resource-allocation (incf (resource-constructed resource)))
	     (funcall (resource-constructor resource)))))))

;;Return object to resource.  Do nothing if max number are already stored.
(defun deallocate-resource-object (resource object)
  (declare (type resource resource)
	   (optimize speed))
  ;;  (loop for index below (resource-pointer resource)
  ;;	when (eq object (aref (resource-objects resource) index))
  ;;	do (error "~S is already in ~S" object resource))
  (locally (declare (optimize (safety 0))) ;Check that argument is correct, but then don't check anything
    (let ((pointer (resource-pointer resource)))
      (unless (>= pointer (resource-max resource))		  ;Stored enough already: just drop
	(when (>= pointer (length (resource-objects resource)))	  ;Array too small
	  (expand-resource resource))
	(setf (aref (resource-objects resource) pointer) object)
	(setf (resource-pointer resource) (1+ pointer)))
      )))

;;Return object or objects to resource
(defmacro deallocate (resource &rest objects)
  `(progn
     ,@(loop for object in objects collect `(deallocate-resource-object ,resource ,object))))

;;Make a larger array for the resource
;;Not reentrant.
(defun expand-resource (resource)
  (let* ((old (resource-objects resource))
	 (new (make-array (* (length old) 2))))
    (loop for index below (resource-pointer resource) ;Copy saved objects
	  do (setf (aref new index) (aref old index)))
    (setf (resource-objects resource) new)))
  
;;Reset counts, forget objects
(defun reset-resource (resource)
  (dotimes (index (resource-pointer resource))
    (setf (aref (resource-objects resource) index) nil))
  (setf (resource-pointer resource) 0
	(resource-constructed resource) 0
	(resource-allocated resource) 0))    

;;Allocate a resource from RESOURCE and bind OBJECT to it.
;;On (normal or abnormal) exit, deallocated, unless OBJECT has been set to NIL
(defmacro using-resource ((object resource) &body body)
  `(let ((,object (allocate ,resource)))
     (unwind-protect
	 (locally ,@body)	    ;Allow user to declare object type
       (when ,object			;If set to NIL, do not deallocate
	 (deallocate ,resource ,object)))))

;;Serially binds several variables to objects allocated from resources
(defmacro using-resources (bindings &body body)
  (if bindings
      `(using-resource ,(car bindings)
	 (using-resources ,(cdr bindings)
	   ,@body))
    `(locally ,@body)))
			  
