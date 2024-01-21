;;;Loop tags
(in-package "CL-USER")

(defconstant tag-count-bytes 2)
(defconstant tag-count-bits (* tag-count-bytes 8))
(deftype tag-count () `(unsigned-byte ,tag-count-bits))

(defstruct (loop-tag
	    (:conc-name tag-)
	    (:print-object print-loop-tag))
  (created-position (required-argument) :type 4vector) ;Global location of intersection
  (last-position nil)			;Position of last engulfment
  (xi 0.0 :type double-float)		;x_i = l_i/t_i, computed at time of first engulfment.  0.0 until known.
  (minimum-inert-count 2 :type tag-count) ;diamond may not be set inert or deleted unless its count is at least this
  (dumped nil :type boolean)		;Loop with this tag has been written out to a dump.
  (successor-flags 0 :type (unsigned-byte 8)) ;Bits set if this tag has been transmitted to successors
  (handle nil)				;Handle for this tag: NIL unless transmitted to successor
  (bh nil))
  
(defun print-loop-tag (tag stream)
  (print-unreadable-object (tag stream :identity t :type t)
    (when (tag-created-position tag)
      (print-4vector (tag-created-position tag) stream))))

;;Create a new tag for the intersection occurring at a given position
;;HANDLE is NIL -- we will create one if we ever have to transmit the tag
(defun create-loop-tag (position)
  (make-loop-tag :created-position position))

;;Get a handle for tag, creating if necessary.
(defun get-tag-handle (tag)
  (or (tag-handle tag)	;existing handle
      (setf (tag-handle tag) (create-local-handle)))) ;or create one

(defun check-tag-count (count)		;Check if integer fits in tag count
  (when (>= count (expt 2 tag-count-bits))
    (error
     "COUNT ~D is too big to fit in type. Is *loop-preservation-threshold* or *intersection-probability* too small?" 
     count)))

