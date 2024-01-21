(in-package "CL-USER")

(define-simulate-variable *intersection-probability* 1.0)
;;Ordered list of (start-global-time probability).  Before first use *intersection-probability*.
(define-simulate-variable *intersection-probabilities* nil)
;;If set, check every intersection to see if it is a rejoining.  Only counts those where one loop is entirely in the simulation volume.
(define-simulate-variable *count-rejoinings* nil)

(Defun intersection-probability (global-time)
  (loop with probability = *intersection-probability*
	for (start p) in *intersection-probabilities*
	when (> global-time start)
	do (setq probability p)	;Time has come to use this
	finally (return probability)))

(defstruct (vector-parameters (:type (vector double-float)))
  a1 b1
  (a2 0.0)				;OK not to specify these if no second diamond
  (b2 0.0))

;;DEFSTRUCT with :TYPE doesn't automatically given a named type.
(deftype vector-parameters () `(simple-array double-float (4)))

;; Definition of the intersection structure
(declaim (disable-package-locks intersection))
(defstruct (intersection
	     (:print-object print-intersection)
	    )
   diamond-1
   diamond-2 ;This can no longer be NIL
   (vector-parameters (required-argument) :type vector-parameters) ;Array A1, B1, A2, B2
   spacetime
   (performed nil)			;intersection has been performed?
   )

(defun print-intersection (intersection stream)
  (print-unreadable-object (intersection stream :identity t)
			   (format stream "INTERSECTION ")
			   (if (intersection-spacetime intersection)
			       (print-4vector (intersection-spacetime intersection) stream)
			     (format stream "uninitialized"))))

(defvar *intersections-performed*)	;Any type of intersection
(defvar *rejoinings-performed*)		;Only used if :count-rejoinings. Only counts cases that can be definitively identified.
(defvar *suppressed-rejoinings*)	;Rejoinings suppressed to avoid monsters
(defvar *intersections-unlucky*)	;Number of intersections skipped because of *intersection-probability*


(defvar *loop-counter* 0)
(defvar *bh-loop-counter* 0)

;;Bound to list of intersections currently being performed to prevent premature erasure
(defvar *intersections-being-performed* nil)

;;Fudge factors for the intersection code
(declaim (double-float solution-accuracy))
(defparameter solution-accuracy 1e-10)

;;And this diamonds to the cells structure and check it for intersections with all other diamonds
;;Record all intersections to be performed later.
(defun check-for-intersections (diamond &optional predecessor-diamond)
   (mark-exempt-diamonds diamond predecessor-diamond t) ;Set processed bits to avoid diamonds that couldn't make a loop
   (do-diamond-cells-diamonds (d diamond :add t)
       ;;(when (diamond-abandoned-p d) (error "~D is abandoned" d))
       ;;should not compare diamonds just produced at intersection (they share start pt.)
       (let ((intersection (handle-possible-intersection diamond d)))
       ;;If intersection, make sure it is inside the other diamond.  We don't have to check that it is inside
       ;;our diamond, because that diamond was just created (by us or someone else), so isn't cut.
            (when intersection
                  (let* ((d1 (intersection-diamond-1 intersection))
                         (d2 (intersection-diamond-2 intersection))
                         (vector-a-b (intersection-vector-parameters intersection))
                         (a1 (vector-parameters-a1 vector-a-b))
                         (b1 (vector-parameters-b1 vector-a-b))
                         (a2 (vector-parameters-a2 vector-a-b))
                         (b2 (vector-parameters-b2 vector-a-b)))
                        (when (and (not (BH-intersection-wrong-side-p d1 :a a1 :b b1))
                                   (not (BH-intersection-wrong-side-p d2 :a a2 :b b2))
                                   (check-intersection-in-my-future intersection))  ;is it in my volume and in the future ?
                              (new-intersection intersection)  ;Arrange for it to get done
                        )
                  )
            )
       )
   )
   (mark-exempt-diamonds diamond predecessor-diamond nil) ;Reset processed bits
)


;defined by SJW to see if intersection point is at the ficititious part of BH-diamond
#||
(defun fake-intersection (intersection)
 (let* ((d1 (intersection-diamond-1 intersection))
        (d2 (intersection-diamond-2 intersection))
        (vector-a-b (intersection-vector-parameters intersection))
        (a1 (vector-parameters-a1 vector-a-b))
        (b1 (vector-parameters-b1 vector-a-b))
        (a2 (vector-parameters-a2 vector-a-b))
        (b2 (vector-parameters-b2 vector-a-b)))
        ;;unless at least one of d1 and d2 are BH-diamond, otherwise return nil for fake-intersection
        (unless (or (eq (diamond-e d1) :BH)
                    (eq (diamond-e d2) :BH)
                    (eq (diamond-w d1) :BH)
                    (eq (diamond-w d2) :BH))
                (return-from fake-intersection nil))
        ;;what if one of d1 and d2 is BH-diamond, return t if intersection point lies inside BH
        (when (eq (diamond-e d1) :BH) 
              (if (< b1 a1) (return-from fake-intersection nil) (return-from fake-intersection t)))
        (when (eq (diamond-e d2) :BH) 
              (if (< b2 a2) (return-from fake-intersection nil) (return-from fake-intersection t)))
        (when (eq (diamond-w d1) :BH) 
              (if (< a1 b1) (return-from fake-intersection nil) (return-from fake-intersection t)))
        (when (eq (diamond-w d1) :BH) 
              (if (< a2 b2) (return-from fake-intersection nil) (return-from fake-intersection t)))
 )
)
||#

(defun BH-intersection-wrong-side-p (diamond &key a b)
 (if (or (and (eq (diamond-e diamond) :BH) (>= b a))
         (and (eq (diamond-w diamond) :BH) (>= a b))
	 (and (eq (diamond-e diamond) :BHdeleted) (>= b a))
	 (and (eq (diamond-w diamond) :BHdeleted) (>= a b))
	 )
     t nil)
)

;;Mark or unmark diamonds that we could never intersect with, because they are too close to us to make a loop
(defun mark-exempt-diamonds (diamond predecessor mark-p)
  (setf (diamond-processedp diamond) mark-p) ;Can't intersect ourselves
  (when predecessor (setf (diamond-processedp predecessor) mark-p)) ;Can't intersect predecessor
  (mirror-images (mark-exempt-diamonds-e diamond mark-p)))

(mirror-images
(defun mark-exempt-diamonds-e (diamond mark-p)
  (declare (optimize speed))
  (let ((a 1)
	(b 1))
    (declare (type (unsigned-byte 20) a b))
    (loop 
     (let ((ne (diamond-ne diamond)))	;Try to go NE
       ;(cond (ne
        (cond ((and ne (diamondp ne))			;New diamond has same a, new b
	      (incf b)
	      (setq diamond ne))
	     (t				;No NE
	      (setq diamond (diamond-se diamond)) ;Try SE
              ;(unless diamond (return nil))
	      (unless (and diamond (diamondp diamond)) (return nil)) ;Ran out: nothing more to mark
	      (incf a)))		;Same b, new a
       ;;In order to be a candidate for intersection, we must have at least 2 a's and 3 b's, or
       ;;2 b's and 3 a's
       (when (and (>= a 2) (>= b 2) (not (= a b 2)))
	 (return nil))			;If intersection possible, we are done
       (setf (diamond-processedp diamond) mark-p))))))  

;;Decide whether to do an intersection.
;;Arrange for intersection to be performed, by putting it in the queue
(defun new-intersection (intersection)
  (let ((probability (intersection-probability (global-time (current-time)))))
    (if (or (= probability 1.0) (< (random 1.0) probability))
	(record-intersection intersection)
      (skip-unlucky-intersection intersection))))

;; Adds a intersection to the calendar(s) and the hash tables
;; for both diamonds involved in the intersection.
(defun record-intersection (intersection)
  (let*((d1 (intersection-diamond-1 intersection))
	(d2 (intersection-diamond-2 intersection)))
    (mirror-images
     (when (or (eq (diamond-e d1) d2)
	       (eq (diamond-e d2) d1))
       (error "About to record an intersection on adjacent diamonds ~S and ~S" d1 d2)))
    (push intersection (diamond-pending-intersections d1)) ;Add to list of pending work for diamonds
    (push intersection (diamond-pending-intersections d2))
    (add-to-calendar intersection))	;Put intersection in calendar to be done when its time comes
  )

;;The intersection will not be done because of *intersection-probability*
;;Update the tags so that we don't think a loop is non-self-intersecting prematurely
(defun skip-unlucky-intersection (intersection)
  (incf *intersections-unlucky*)
  (unless *delete-unlucky-loops*	;Unless disabled by user
    (flet ((fix-tag (diamond)
		    (let ((tag (diamond-tag diamond))
			  (new-count (+ (diamond-countup diamond) 3))) ;See notes
		      (check-tag-count new-count)
		      (setf (tag-minimum-inert-count tag) (max new-count (tag-minimum-inert-count tag))))))
      (fix-tag (intersection-diamond-1 intersection))
      (fix-tag (intersection-diamond-2 intersection)))))

;; Delete a particular intersection from the hash table of a diamond, if it is there
(defun erase-intersection (diamond intersection)
  (setf (diamond-pending-intersections diamond)
	(remove intersection (diamond-pending-intersections diamond)))) ;; Remove is better than delete

(defvar monster-kink-number 50)	;Number of small segments in a row to identify a monster
;;Conformal length scale used to determine monster-hood is this factor times the largest edge
(defvar monster-length-multiplier 0.1)

;; Returns T if a (large or small) loop is rejoining a segment which is so kinky that
;; a monster is probably forming.
;; Species one: a small loop wants to rejoin a segment kinky in both directions
;; Species two: a demi-kinky loop wants to rejoin another demi-kinky loop
(defun monster-maker-p (time d1 a1 b1 d2 a2 b2)
  (let ((data (list (small-left-structure-p d1 d2 monster-kink-number monster-length-scale time a1 b1)
		    (small-right-structure-p d1 d2 monster-kink-number monster-length-scale time a1 b1)
		    (small-left-structure-p d2 d1 monster-kink-number monster-length-scale time a2 b2)
		    (small-right-structure-p d2 d1 monster-kink-number monster-length-scale time a2 b2))))
    (when (and (find-if #'floatp data)			     ;Many kinks in a short region in any search?
	       (not (find :self-intersection data)))	     ;No self-intersection found?
      (let ((rigorous (rigorous-self-intersection-p d1 d2))) ;self intersection looking more carefully?
	(unless (eq rigorous t)				     ;It's definitely self intersection, don't do it
	  (let ((count1 (loop-count d1))
		(count2 (loop-count d2)))
	    (cond (rigorous		;Don't know for sure if it is a self-intersection
		   (warn "Not suppressing a possible monster with loop lengths ~D and ~D because rigorous = ~S~%"
			  count1 count2 rigorous))
		  (t			;Definitely a rejoining.  Suppress it.
		   (format t "~&Suppressing a possible monster with loop lengths ~D and ~D.~%"
			   count1 count2)
		   (format t " Intersection position ~S~%" (diamond-position d1 :a a1 :b b2))
		   (report-progress "%")	;Found a monster
		   (incf *suppressed-rejoinings*)
		   t))))))))

(mirror-images
;; decides if a segment has many kinks within a small conformal length SIZE
;; output is length encountered before NUMBER kinks, or number of kinks before length limit was reached
;; or :self-intersection if other diamond is found on same string 
(defun small-left-structure-p (d1 d2 number size time
				  parameter1 parameter2) ;(a1, b1) regardless of mirror image
  (let ((w1 (diamond-w d1))
	(measured (multiple-value-call #'path-segment-length d1 parameter1 parameter2
				       (find-left-edge-position d1 time t)))) ;a, b of edge position
    (unless (diamondp w1)
      (return-from small-left-structure-p :end-point))
    (loop with next
	  with entry
	  for n from 1 to number
	  for previous = d1 then this
	  for this = w1 then next
	  do (setq next (diamond-w this)) ;move leftward
	  unless (diamondp next)
	     do (return-from small-left-structure-p :end-point)
	  do (setq entry (cons this (nconc (multiple-value-list (find-right-edge-position this time nil))
					   (multiple-value-list (find-left-edge-position this time nil)))))
	  when (position nil entry) ;exit if the boundary of a string is found
	     do (return-from small-left-structure-p :end-point)
	  do (setq measured (+ measured (* (apply #'path-segment-length entry)))) ;add physical length of segment
	  when (> measured size)
	     do (return-from small-left-structure-p n) ;too long, not enough kinks found
	  when (eq d1 this)		;e.g. two loops colliding
	     do (return-from small-left-structure-p n) ;loop is found with not enough kinks
	  when (eq d2 this)
	     do (return-from small-left-structure-p :self-intersection)) ;self-intersection
    measured))


;; Gathers (number) distinct right-movers to the right of diamond 
(defun gather-a-to-right (d1 e1 d2 e2 number)
  (let ((set nil))
    (unless (and (diamondp d1) (diamondp e1))
      ;(format t "~%No neighbor found from coordinates! ~%")
      (return-from gather-a-to-right nil)) ;can't gather
    (loop
        with next
	for previous = d1 then this
	for this = e1 then next
	do (setq next (diamond-e this))
	while (diamondp next)
	when (eq next (diamond-se this))
	   do (push (diamond-a next) set)
	when (and (eq next e2) (eq this d2))
	   do (return-from gather-a-to-right nil) ;self-intersection, no need to gather
	when (and (eq next e1) (eq this d1))
	   do (push :loop set)
	   and do (return-from gather-a-to-right set)
	while (< (length set) number)
	)
    (push :not-loop set)
    set))
)					;mirror images

;;Is this a self-intersection.  Returns:
;; T -- we walked from one diamond to the other, so it is definitely a self-intersection
;; NIL -- we walked from one diamond around in a loop without finding the other, so it is definitely not.
;; :MAYBE -- The two strings are open and disjoint, so we don't know.
(defun rigorous-self-intersection-p (d1 d2)
  (assert (not (eq d1 d2)))
  (flet ((test (this other)		;Go eastward from one diamond looking for the other
	   (loop for d = (diamond-e this) then (diamond-e d) ;Go eastward around loop
		 while (diamondp d)	;If we ran off the end, exit this loop
		 when (eq d this)	;Looped back around?
		    do (return-from rigorous-self-intersection-p nil) ;Definitely not self: return NIL from function
		 when (eq d other)	;Found other diamond
		   do (return-from rigorous-self-intersection-p t)))) ;Definitely a self-intersection: return T
    (test d1 d2)			;Try starting at first diamond
    (test d2 d1)			;If the first test ran off end, try starting at second one
    :maybe))				;Both tests ran off end.

;;Transforms the numerical values of the parameters on the plane for an intersection
;; into spacetime position vector.
(defun solution-to-spacetime-position (solution d1 d2)
  (let* ((d1-position (diamond-position d1 :a (vector-parameters-a1 solution) :b (vector-parameters-b1 solution)))
	 (d2-position (diamond-position d2 :a (vector-parameters-a2 solution) :b (vector-parameters-b2 solution)))
	 (error (4vector-Euclidian-distance d1-position d2-position)))
    (when (> error solution-accuracy)
      (describe d1)
      (describe d2)
      (error "intersecton point is split by too much! ~D > ~D" error solution-accuracy))
    (prog1 (4vector-scale (4vector+ d1-position d2-position) 0.5d0)
      (deallocate 4vectors d1-position d2-position))))

;; Find the sign of (s1 x s2) . (v1 - v2), which is the homology number of the intersection
;; The order of d1 and d2 does not matter.  Assumes the diamonds are infinite sheets.
(defun old-flat-intersection-number (d1 d2)
  (if (eq d1 d2) (return-from old-flat-intersection-number 0))
  (let* ((p1 (4vector-scale (diamond-p d1) (/ 0.5 (4vector-t (diamond-p d1)))))
	 (q1 (4vector-scale (diamond-q d1) (/ 0.5 (4vector-t (diamond-q d1)))))
	 (p2 (4vector-scale (diamond-p d2) (/ 0.5 (4vector-t (diamond-p d2)))))
	 (q2 (4vector-scale (diamond-q d2) (/ 0.5 (4vector-t (diamond-q d2)))))
	 (v1 (4to3vector (4vector+ p1 q1)))
	 (v2 (4to3vector (4vector+ p2 q2)))
	 (s1 (4to3vector (4vector- q1 p1))) ;ok if not a unit vector, normalization is arbitrary
	 (s2 (4to3vector (4vector- q2 p2)))
	 (value (3vector-dot (3vector-cross-product s1 s2) (3vector- v1 v2))))
    (cond ((plusp value)
	   1)
	  ((minusp value)
	   -1)
	  (t
	   0))))
    
(defparameter uv-cutoff-x nil)		;forbid loop formation below this scaling size

;; decides if an intersection would form a loop too small to be allowed.
(defun uv-cutoff (intersection)
  (let* ((d1 (intersection-diamond-1 intersection))
	 (d2 (intersection-diamond-2 intersection))
	 (intersection-point (intersection-spacetime intersection))
	 (time (global-time (4vector-t intersection-point)))
	 (vector-a-b (intersection-vector-parameters intersection))
	 (a1 (vector-parameters-a1 vector-a-b))
	 (b1 (vector-parameters-b1 vector-a-b))
	 (uv-length (* time uv-cutoff-x))
	 (l-length (small-left-structure-p d1 d2 1000 uv-length time a1 b1))
	 (r-length (small-right-structure-p d1 d2 1000 uv-length time a1 b1)))
    (if (or (eq l-length :self-intersection)
	    (eq r-length :self-intersection))
	t    ;(format t "%")				;ignore this intersection if its a small loop forming
      nil     ;(format t "!")
      )))

;;Don't make loops that (are known to) have fewer than this number of segments
(defvar *minimum-loop-count* nil)

;;Return true if and this intersection would formal loop with fewer than minimum-loop-count segments
(defun check-minimum-loop-count (intersection count)
  (let ((d1 (intersection-diamond-1 intersection))
	(d2 (intersection-diamond-2 intersection)))
    (mirror-images
     (loop for d = (diamond-e d1) then (diamond-e d) ;D scans over diamond to right of D1
	   for index from 1 below count		     ;If we found COUNT of these, loop is big enough, so exit now
	   while (diamondp d)			     ;Exit if we found the end of the string
	   when (eq d d2)		;If we found the other side of the intersection before COUNT diamonds
	   do (return-from check-minimum-loop-count t)  ;then this a too-short loop: return T
	   ))))

(defun perform-intersection (intersection)
  (unless (intersection-performed intersection)
    (with-modification-group
      (let* ((d1 (intersection-diamond-1 intersection))
	     (d2 (intersection-diamond-2 intersection))
	     (intersection-point (intersection-spacetime intersection))
	     (vector-a-b (intersection-vector-parameters intersection)))
;;	(let ((count1 (loop-count d1)))
;;	  (when (and count1 (> count1 25))
;;	    (format t "~&Intersection between length ~S, count ~S, and length ~S, count ~S~% at ~S"
;;		   (loop-length d1) count1 (loop-length d2) (loop-count d2) intersection-point)))

	
	(if (or (null (diamond-start d1)) (null (diamond-start d2))) ;after removing a string we could have intersections that are not present anymore
	    (progn
	      (format t "~%Unitalized intersection~%")
	      (discard-object intersection)
	      (setf (intersection-performed intersection) t)
	      (return-from perform-intersection t)))

	(when (and *bh-start* (null *pointbh*))
	  (when (and (check-point-inbh intersection-point) (> (current-time) (local-time *bh-start*)))
	    (format t "Int inside BH~%")
	    (discard-object intersection)
	    (setf (intersection-performed intersection) t)
	    (return-from perform-intersection t)))

	(if (null (handle-possible-intersection d1 d2)) ;This is also a similar case
	    (progn
	      (format t "Null Intersection")
	      (discard-object intersection)
	      (setf (intersection-performed intersection) t)
	      (return-from perform-intersection t))
	  )

	(mirror-images
	 (when (or (eq (diamond-e d1) :BHdeleted)
		   (eq (diamond-e d2) :BHdeleted)
		   (eq (diamond-e d1) :BHpropdel)
                   (eq (diamond-e d2) :BHpropdel)
		   )
	   (discard-object intersection)
	   (setf (intersection-performed intersection) t)
	   (return-from perform-intersection t)))


	
	(mirror-images
	 (when (or (eq (diamond-e d1) d2)
		   (eq (diamond-e d2) d1))
	   (error "About to perform an intersection on adjacent diamonds ~S and ~S" d1 d2)))
	;;Remove intersection from pending-intersection lists and calendar.  It does not reuse the data structure,
	;;so it is OK that we use this object even after discarding it.
	(discard-object intersection)
	(unless (or (and uv-cutoff-x (uv-cutoff intersection))	;explicitly forbidding loop formation below some length
		    (and *minimum-loop-count* (check-minimum-loop-count intersection
									*minimum-loop-count*)) ;Too few segments?
		    (and (plusp monster-length-scale)
			 (monster-maker-p (4vector-t intersection-point)
					  d1 (vector-parameters-a1 vector-a-b) (vector-parameters-b1 vector-a-b)
					  d2 (vector-parameters-a2 vector-a-b) (vector-parameters-b2 vector-a-b))))
	  ;;We have decided to actually perform the intersection
	  (incf *intersections-performed*)
	  (when *count-rejoinings*
	    (when (eq (rigorous-self-intersection-p d1 d2) nil) ;If definitively a rejoining
	      (incf *rejoinings-performed*)))
	  (let* ((global-point (globalize-position intersection-point))
		 (a1 (vector-parameters-a1 vector-a-b))
		 (b1 (vector-parameters-b1 vector-a-b))
		 (a2 (vector-parameters-a2 vector-a-b))
		 (b2 (vector-parameters-b2 vector-a-b))
		 ;;New diamonds made from quarters of old.  Links to and from d1 are destroyed here
		 (e1 (east-quarter d1 intersection-point global-point :a a1 :b b1))
		 (w1 (west-quarter d1 intersection-point global-point :a a1 :b b1))
		 (e2 (east-quarter d2 intersection-point global-point :a a2 :b b2))
		 (w2 (west-quarter d2 intersection-point global-point :a a2 :b b2)))
	    (divide-pending-intersections d1 :a a1 :b b1 :east e1 :west w1) ;discard or move intersections to west or east quarters
	    (divide-pending-intersections d2 :a a2 :b b2 :east e2 :west w2)
	    (rescale-right-junctions d1 e1 :a a1 :b b1)
	    (rescale-right-junctions d2 e2 :a a2 :b b2)
	    (rescale-left-junctions d1 w1 :a a1 :b b1)
	    (rescale-left-junctions d2 w2 :a a2 :b b2)
	    (discard-object d1)		;this makes neighboring diamonds link to nil, so better not have any.
	    (discard-object d2)   
	    (propagate-cut-ne e1 a1)	
	    (propagate-cut-nw w1 b1)
	    (propagate-cut-ne e2 a2)
	    (propagate-cut-nw w2 b2)
	    (handle-new-diamond e1 :predecessor :meiosis) ;predecessor = :meiosis prevents rechecking for intersections
	    (handle-new-diamond w1 :predecessor :meiosis)
	    (handle-new-diamond e2 :predecessor :meiosis)
	    (handle-new-diamond w2 :predecessor :meiosis)
	    (advance-cut-point e1 w1 e2 w2 global-point) ;Create new reconnected future diamonds
	    )
	  (setf (intersection-performed intersection) t))
	))))




(mirror-images
 ;; Creates east quarter from diamond, and links in/out east neighbors
(defun east-quarter (diamond point global-start &key a b)
  (let ((new (make-diamond :start (if (or (eq (diamond-e diamond) :BH)
					  (eq (diamond-e diamond) :BHdeleted)) 
                                      (diamond-position diamond :a b :b b)
				    (diamond-position diamond :a 0.0 :b b)
				    )
			   :left point
			   :right (if (or (eq (diamond-e diamond) :BH)
					  (eq (diamond-e diamond) :BHdeleted))
                                      (diamond-position diamond :a b :b a)
				    (diamond-right diamond)
                                  )
			   :end (if (or (eq (diamond-e diamond) :BH)
					(eq (diamond-e diamond) :BHdeleted))
                                    (diamond-position diamond :a a :b a)
				  (diamond-position diamond :a a :b 1.0)
                                )
			   :se (if (or (eq (diamond-e diamond) :BH)
				       (eq (diamond-e diamond) :BHdeleted))
				   nil
				 (diamond-se diamond))
			   :ne (if (or (eq (diamond-e diamond) :BH)
				       (eq (diamond-e diamond) :BHdeleted))
				   (if (eq (diamond-e diamond) :BH)
				       :BH
				     :BHdeleted)
				 (diamond-ne diamond))
			   :b-kink-created global-start
			   :a-kink-created (diamond-a-kink-created diamond)
			   :tag (diamond-tag diamond)
			   :countup (diamond-countup diamond)
			   :bh (diamond-bh diamond)
			   )))		;Can't be inert, so there's no need to set that
    ;(if (diamond-ne diamond) (setf (diamond-sw (diamond-ne diamond)) new)) ;link in to east neighbor
    ;(if (diamond-se diamond) (setf (diamond-nw (diamond-se diamond)) new))
    (if (diamondp (diamond-ne diamond)) (setf (diamond-sw (diamond-ne diamond)) new)) ;modified by SJW
    (if (diamondp (diamond-se diamond)) (setf (diamond-nw (diamond-se diamond)) new)) ;modifed  by SJW
    (setf (diamond-ne diamond) nil)
    (setf (diamond-se diamond) nil)
    new))				;output new diamond


;; takes an eastern quarter of a diamond and propagates the cut to all ne neighbors
;; pending intersections in ne neighbors are moved or discarded, 
(defun propagate-cut-ne (diamond a)
  (let ((east (diamond-ne diamond)))
    ;(when east
    (when (diamondp east)    ;;modified by SJW
    ;(when (or (eq diamond :BH) (eq east :BH)) (print "BH-diamond detected"))  ;;tested by SJW
      (maybe-delete-diamond-cells east)	;Before reshaping, get rid of of stuff that depends on shape
      (maybe-delete-from-calendar east)
      (when (diamond-finalp east)
	(setf (diamond-finalp east) nil)
	(delete-final-diamond east))
      (report-progress "+")
      (let* ((start-1 (diamond-start east))
	     (start-2 (diamond-right diamond))
	     (start-time (4vector-t start-1))
	     (left-time (4vector-t (diamond-left east)))
	     (time (4vector-t (diamond-end diamond)))
	     (new-a (if (eq start-1 start-2) a         ;ordinary case
		      (/ (- time start-time) 
			 (- left-time start-time)))))  ;case when east is itself a west quarter from some earlier intersection 
	(setf (diamond-left east) (diamond-end diamond))
	(unless (= new-a 1.0)
	  (if (or (eq (diamond-e east) :BH)
		  (eq (diamond-e east) :BHdeleted))
	      (progn
		(format t  "~%cut propagates to ~S diamond  ~S~%" (diamond-e east) east)
		;(setf (diamond-right east) (3to4vector (4to3vector (3vector- (diamond-start east) (3vector- (diamond-left east) (diamond-start east)))) (4vector-t (diamond-left east))))
		(if (or (eq *era* :radiation)
			(eq *era* :radiation-smooth)
			(eq *era* :matter)
			(eq *era* :power))
		    (progn
		     (setf (diamond-end east) (compute-bh-diamond-end east (diamond-left east)))
		     (setf (diamond-right east) (4vector- (4vector+ (diamond-start east) (diamond-end east)) (diamond-left east))))
		  (progn
		    (setf (diamond-right east) (3to4vector (4to3vector (3vector- (diamond-start east) (3vector- (diamond-left east) (diamond-start east)))) (4vector-t (diamond-left east))))
		    (setf (diamond-end east) (compute-bh-diamond-end east (diamond-left east))))   
		  )
		)   ;;added by SJW to propagate a cut to BH-diamond with new diamond-end
	    (setf (diamond-end east) (diamond-position-wrap-dumps east :a new-a :b 1.0)) ;If reading dump, must wrap
	    )
	  (divide-pending-intersections east :a new-a :b 0.0 :east east :west nil)
	  (rescale-right-junctions east east :a new-a :b 0.0)
	  (rescale-left-junctions east east :a new-a :b 0.0))
	(when (and (not (eq start-1 start-2))
		   (4vector= start-1 start-2 fudge-coordinates))
	  (error "shared corner close but not eq: difference = ~D" (4vector- start-1 start-2)))
	(handle-new-diamond east :predecessor :meiosis)
	(propagate-cut-ne east new-a)  ;continue until no ne neighbor 
	))))
)					;mirror-images




;;Tell if two closed loops are actually the same loop
(defun closed-loops-same-p (loop-1 loop-2)
  (let ((d11 (car (first loop-1)))	;First diamond in loop-1
	(d12 (car (second loop-1))))	;Next diamond
    (loop for (d21) in loop-2		;Go through diamonds of loop-2
	  for (d22) in (cdr loop-2)	;Next diamond 
	  do (when (null d22)		;Ran out of loop
	       (setq d22 (car (first loop-2))))	;First is after last
	  thereis (and (eq d11 d21) (eq d12 d22)))))


;; Creates two new diamonds starting at intersection point, links them in to neighbors
;; and puts them in calendar and cells
(defun advance-cut-point (east-1 west-1 east-2 west-2 global-start)
  (let* ((left-1 (diamond-end west-1))
	 (right-1 (diamond-end east-1))
	 (left-2 (diamond-end west-2))
	 (right-2 (diamond-end east-2))
	 (start (diamond-left east-1))
	 (west-tag (create-loop-tag global-start))
	 (east-tag (create-loop-tag global-start))
	 (west (make-diamond :start start
			     :left left-1
			     :right right-2
			     :tag west-tag
			     :a-kink-created global-start :b-kink-created global-start
			     :sw west-1
			     :se east-2
			     :bh (diamond-bh west-1)))
	 (east (make-diamond :start start
			     :left left-2
			     :right right-1
			     :tag east-tag
			     :a-kink-created global-start :b-kink-created global-start
			     :sw west-2
			     :se east-1
			     :bh (diamond-bh east-1))))
    (mirror-images
     (setf (diamond-end west) (compute-diamond-end west)) ;Compute the end points of the sons
     (setf (diamond-ne west-1) west)
     (setf (diamond-nw east-2) west)
     )				;mirror images
    ;;Do everything that needs doing with new diamonds.  Predecessor is nil because there is no diamond ending where
    ;;this one begins, only diamonds to its sw and se.
    (handle-new-diamond west) ;only locations of :sw-cut etc. are needed.
    (handle-new-diamond east)
    ))
	  
;; Assigns pending intersections to be in west quarter, east quarter, or nowhere
;; Original is the diamond split into east and west quarters, and a,b is the splitting location
;; For slicing, one of east, west is nil, and the other is eq original
(defun divide-pending-intersections (original &key a b east west)
  (let ((pending-intersection-list (diamond-pending-intersections original)))
    (when pending-intersection-list   
      (setf (diamond-pending-intersections original) nil) ;remove all intersections from original, they will be put back as needed.
      (mapcar #'(lambda (intersection) (divide-pending-intersection intersection original a b east west))
	      pending-intersection-list))))

;; a-b parameters are modified, pending-intersection-list is modified, spacetime is not
(defun divide-pending-intersection (intersection original a b east west)
  (let* ((d1 (intersection-diamond-1 intersection))
	 (d2 (intersection-diamond-2 intersection))
	 (other (if (eq d1 original) d2 d1))
	 (original-diamond-number (if (eq d1 original) 0 2))
	 (vector-a-b (intersection-vector-parameters intersection))
	 (ai (aref vector-a-b original-diamond-number))
	 (bi (aref vector-a-b (1+ original-diamond-number))))
    (if (or (null (diamond-start d1))
	    (null (diamond-start d2)))
	(progn
	  (format t "Pending in in ~S and ~S" d1 d2)
	  (discard-object intersection)
	  (return-from divide-pending-intersection))) ;modified for the cases after creating BHs
    (when east				;i.e. not a nw slicing
      (setf (aref vector-a-b original-diamond-number) (/ ai a))	;modify solution parameters
      (setf (aref vector-a-b (1+ original-diamond-number)) (1+ (/ (- bi 1.0)
								  (- 1.0 b))))
      (if (eq d2 other) (setf d1 east)	;replace original diamond with east
	(setf d2 east))
      (when (check-solution-inside-diamond vector-a-b d1 d2);east other) ;still a good intersection, but with east
	(setf (intersection-diamond-1 intersection) d1) ;change diamonds accordingly
	(setf (intersection-diamond-2 intersection) d2) 
	(setf (intersection-vector-parameters intersection) vector-a-b) ;and vector parameters
	(push intersection (diamond-pending-intersections east)) ;record intersection
	(return-from divide-pending-intersection)))
    (when west				;not a ne slicing, and east didn't contain the intersection
      (setf (aref vector-a-b original-diamond-number) (1+ (/ (- ai 1)
							     (- 1 a))))	;destructively modify solution parameters
      (setf (aref vector-a-b (1+ original-diamond-number)) (/ bi b))
      (if (eq d2 other) (setf d1 west)
	(setf d2 west))
      (when (check-solution-inside-diamond vector-a-b d1 d2);west other) ;still a good intersection, but with west
	(setf (intersection-diamond-1 intersection) d1) ;change diamonds accordingly
	(setf (intersection-diamond-2 intersection) d2)
	(setf (intersection-vector-parameters intersection) vector-a-b)
	(push intersection (diamond-pending-intersections west))
	(return-from divide-pending-intersection)))
    (discard-object intersection)))		;neither east nor west contained this intersection


(mirror-images
;; Rescale the a and b.of existing junctions in original to be compatible with new diamond. 
;; (a,b) will eventually be the left or right corner of new.
;; Original may be eq to new in which case any junctions in new should be overwritten instead of appended.
;; (a,b) will either be to the left of all right-junctions (i.e. an intersection point)
;; or to the right of all right-junctions (from merge-diamonds).  
(defun rescale-right-junctions (original new &key a b)
  (let ((original-received-dump-junctions (received-right-dump-junctions original))
	(original-created-dump-junctions (created-right-dump-junctions original))
	(original-rejoining-junction (right-rejoining-junction original))
	(new-rejoining-junction (right-rejoining-junction new))
	(directions nil))		;list of directions to junctions.  Should all be the same.
    (unless (eq (length original-created-dump-junctions)
		(length (remove-duplicates original-created-dump-junctions)))
      (error "Found duplicate created dump junctions ~D -> ~D" 
	     (length original-created-dump-junctions) 
	     (length (remove-duplicates original-created-dump-junctions))))
    (when (eq original new) 		;If we're modifiying a diamond, 
      (if (rejoining-junction-p new-rejoining-junction) (setf (right-rejoining-junction new) nil))
      (setf (received-right-dump-junctions new) nil)
      (setf (created-right-dump-junctions new) nil)) ;don't append, overwrite.
    (cond ((rejoining-junction-p original-rejoining-junction)
	   (push (rescale-junction original-rejoining-junction :a a :b b) directions)) ;modify a,b and collect direction
	  ((vv-junction-p original-rejoining-junction)
	   (let ((structure (vv-face-ref (vv-junction-site original-rejoining-junction)
					 (vv-junction-axis1 original-rejoining-junction)
					 (vv-junction-axis2 original-rejoining-junction))))
	     (cond ((eq (face-point-left structure) original) ;Replace with old with new
		    (setf (face-point-left structure) new))
		   ((eq (face-point-right structure) original) ;Replace with old with new
		    (setf (face-point-right structure) new))
		   (t (error "Couldn't find ~S in face of ~S" original original-rejoining-junction)))))
	  (t (assert (or (null original-rejoining-junction)
			 (eq original-rejoining-junction :BHdeleted)
			 (eq original-rejoining-junction :BHeatit)
			 (eq original-rejoining-junction :BHpropdel)
                         (eq original-rejoining-junction :BH)  ;;added by SJW
			 (eq original-rejoining-junction :deleted)))))
    (when (junction-p original-rejoining-junction)
      (setf (junction-left-diamond original-rejoining-junction) new)	;junction should point to new diamond
      (setf (right-rejoining-junction new) original-rejoining-junction)) ;vv- or rejoining-junction goes into new
    (dolist (junction original-received-dump-junctions)	;only care about future dumps
      (push (rescale-junction junction :a a :b b) directions) ;modify a,b and collect direction
      (setf (junction-left-diamond junction) new))
    (dolist (junction original-created-dump-junctions)
      (push (rescale-junction junction :a a :b b) directions)) ;modify a,b and collect direction
    (when (and (position :right directions)    ;check that all were on same side or in past
	       (position :left directions))
      (error "junctions are in both directions from cut point. ~D" directions))
    (setf (received-right-dump-junctions new) (append (received-right-dump-junctions new)
						      original-received-dump-junctions))
    (setf (created-right-dump-junctions new) (append (created-right-dump-junctions new)
						     original-created-dump-junctions))
    ))
)					;mirror-images

;; Modifies the a b values of a junction in a diamond that was sliced at location (a,b).
(defun rescale-junction (junction &key a b)
  (let* ((old-a (rejoining-junction-a junction))
	 (old-b (rejoining-junction-b junction))
	 (direction (cond
		     ((and (<= old-a a)
			   (<= b old-b))
		      :right)		;junction is to the right of (a,b)
		     ((and (<= a old-a)
			   (<= old-b b))
		      :left)
		     (t
		      (error "active junction (a,b) (~D,~D) time-like from pivot location (a,b) = ( ~D, ~D )" old-a old-b a b))))
	 (new-a (if (eq direction :right) (/ old-a a)
		  (/ (- old-a a)
		     (- 1.0 a))))
	 (new-b (if (eq direction :right) (/ (- old-b b)
					     (- 1.0 b))
		  (/ old-b b))))
    (assert (<= 0.0 new-a 1.0))
    (assert (<= 0.0 new-b 1.0))
    (setf (rejoining-junction-a junction) new-a)
    (setf (rejoining-junction-b junction) new-b)
    direction				;output direction toward junction
    ))


;;Checks that the intersection is to the future of (current-time), in our volume,
;; and compatible with junctions e.g. not to the right of the right-junction.
(defun check-intersection-in-my-future (intersection)
  (let* ((now (current-time))
	 (here (intersection-spacetime intersection))
	 (time (4vector-t here))
	 (d1 (intersection-diamond-1 intersection))
	 (d2 (intersection-diamond-2 intersection))
	 (r-junction-1 (right-rejoining-junction d1))
	 (l-junction-1 (left-rejoining-junction d1))
	 (r-junction-2 (right-rejoining-junction d2))
	 (l-junction-2 (left-rejoining-junction d2))
	 (vector-a-b (intersection-vector-parameters intersection))
	 (a1 (vector-parameters-a1 vector-a-b))
	 (b1 (vector-parameters-b1 vector-a-b))
	 (a2 (vector-parameters-a2 vector-a-b))
	 (b2 (vector-parameters-b2 vector-a-b)))
    (cond ((< time now) 		;in past?
	   nil)				;skip
	  ((not (point-mine here t)) 	;not in our volume?
	   nil)				;skip
	  (t
	   (and (left-of-right-junction r-junction-1 :a a1 :b b1) ;no junction excludes this point?
		(left-of-right-junction r-junction-2 :a a2 :b b2)
		(right-of-left-junction l-junction-1 :a a1 :b b1)
		(right-of-left-junction l-junction-2 :a a2 :b b2)))))) ;then intersection is good
	   
(mirror-images
 ;; Returns t if the specified (a,b) location is within the region of the diamond
 ;; not cut off by the rejoining junction.
(defun left-of-right-junction (junction &key a b)
  (if (rejoining-junction-p junction)	
      (let ((a-j (rejoining-junction-a junction))
	    (b-j (rejoining-junction-b junction)))
	(if (and (<= a-j a)
		 (<= b b-j))
	    t				;left of the right-junction, so returnn t
	  nil))			;not left of the right-junction, return nil
    t))				;no junction, so return t
)					;mirror-images

;;Erases the intersection from the hash tables for the diamonds involved
;;and the calendars.  Intersection might not be in calenders if it is associated with an end-point
;;This is called from perform-intersection, which uses the intersection object even after discarding it, so
;;it's important that we don't put it in a resource here.
(defun discard-intersection (intersection)
  (let*((d1 (intersection-diamond-1 intersection))
	(d2 (intersection-diamond-2 intersection)))
    (erase-intersection d1 intersection) ;Remove from diamond pending-intersection lists.  Might not be necessary
    (erase-intersection d2 intersection) ;if it has been already removed in perform-intersection
    (delete-from-calendar intersection)
    ))

(defun print-intersection-details (intersection)
  (describe intersection)
  (let ((d1 (intersection-diamond-1 intersection))
	(d2 (intersection-diamond-2 intersection))
	(location (intersection-spacetime intersection)))
    (format t "Location has handle ~S~%" (object-handle location nil))
    (format t "Diamond 1 has handle ~S~%" (object-handle d1 nil))
    (describe d1)
    (when d2
      (format t "Diamond 2 has handle ~S~%" (object-handle d2 nil))
      (describe d2))))
