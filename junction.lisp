;;;Junctions for reassembling strings
(in-package "CL-USER")

;;Junction between diamonds: common slots
(defstruct junction
  (right-diamond nil)			;Diamond to right of this junction.  
  (left-diamond nil)			;Diamond to left.
  )

;;Junctions for rejoining diamonds that were split to pass to multiple successors
(defstruct (rejoining-junction (:include junction))
  handle
  left-direction				;Direction of communication to left: 0=I,...,3=K, 4=t.
  right-direction				;Direction of communication to right.
  (dump-time nil)		      ;Global time of dump with which this junction is associated, or NIL
  a					;Location of junction within its diamond
  b
  a0					;Original location of junction within uncut version of its diamond
  b0)


;;Junctions for Vachaspati-Vilenkin initial conditions
(defstruct (vv-junction (:include junction))
  (site (required-argument) :type site)			;Specify face where diamonds join.
  (axis1 (required-argument) :type axis)		;Face is oriented eastward along the string
  (axis2 (required-argument) :type axis))

;;Junctions for explicit initial conditions
(defstruct (initial-junction (:include junction))
  string				;Number of string containing diamond
  diamond)				;Number of diamond inside string on past side of junction

;;There are 3 kinds of junctions: ones passed on two successors for rejoinings, once created by us during dumping
;;and passed to our successor, and ones received from our successor and written out into dumps.

;;The meaning of right is as follows:  For rejoining-junctions, the right junction of the diamond is the one where
;;the string comes to a end as we move rightward.  It will be eventually joined to another junction in the same
;;(or cut versions of the same) diamond.

;;When we write a dump, and the dump comes to an end in a junction, we call that a right-junction
;;for the last diamond, because the dump string ends here.  However, when we pass this same junction
;;to our successor, it becomes a left junction for him, because he will start his dump there.

(mirror-images

(defun right-rejoining-junction (diamond)
  (or (diamond-right-rejoining-junction diamond)
      (and (keywordp (diamond-e diamond)) ;:deleted, :bh, :bhdeleted, ...
	   (diamond-e diamond))))

(defsetf right-rejoining-junction (diamond) (junction)
  `(progn
     (when (and ,junction ;Check there is not one already
		(diamond-right-rejoining-junction ,diamond))
       (error "~S already has a ~(~A~) rejoining junction" ,diamond :right))
     (setf (diamond-right-rejoining-junction ,diamond) ,junction)))

(defmacro created-right-dump-junctions (diamond)
  `(diamond-created-right-dump-junctions ,diamond))
(defmacro received-right-dump-junctions (diamond)
  `(diamond-received-right-dump-junctions ,diamond))

)					;mirror-images

;;Create a new junction for output to successors.  It does not need to go in our data structures unless
;;it involves dumping.  In that case, we need to remember it to dump it later.
(defun create-junction (&optional diamond left-destination right-destination a b)
  (let ((junction (make-rejoining-junction :handle (create-local-handle)
				:left-direction left-destination :right-direction right-destination
				:dump-time *dump-time* :a a :b b :a0 a :b0 b)))
    (mirror-images
     (when (eq left-destination dump-destination) ;Switch from dumping to successor
       (push junction (created-right-dump-junctions diamond)))) ;Save as that
    junction))


;;Remove junction from any diamonds it is in and from hash tables
(defun discard-rejoining-junction (junction)
  (mirror-images
   (let ((diamond (junction-right-diamond junction)))
     (when diamond
       (setf (left-rejoining-junction diamond) nil))))
  (setf (object-handle junction) nil))

(defun initialize-junctions ()
  )					;Nothing to do since we no longer use hash tables

;;Map over junctions that are in the *handle-objects* system.
;;Call function with a list of the junctions with a given handle
(defun map-junctions (function)
  (maphash
   #'(lambda (handle objects)
       (declare (ignore handle))
       (when (and objects
		  (junction-p (car objects))) ;If one is a junction, they should all be
	 (funcall function objects)))
   *handle-objects*))

