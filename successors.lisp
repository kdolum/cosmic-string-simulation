;;;Passing things on to successors
(in-package "CL-USER")

;;Where diamond intersects the line between successors on the future boundary
(defstruct successor-intersection
  a
  b
  successor-1
  successor-2
  ijkl)

;;List of successors whose volumes this point is in.  Conservative means only if we're sure
(defun point-successors (location conservative-p)
  ;;Possible successors are those whose ijkl coordinate of the endpoint is > 1.0.
  (let ((ijkl (xtoi location)))
    (loop for index below 4
	  when (> (aref ijkl index) (if conservative-p (+ 1.0 fudge-ijkl) (- 1.0 fudge-ijkl)))
	  collect index into result
	  finally
	  (when (and *dump-time* (> (global-time (4vector-t location))
				    (if conservative-p (+ *dump-time* fudge-global-coordinates)
				       (- *dump-time* fudge-global-coordinates))))
	    (push 4 result))			;If time > *dump-time*, Include dump surface
	  (return result))))


;;Return a list of successors and junctions for this diamondg: (successor junction successor ... junction successor)
;;If left-junction is given and isn't :deleted, start there.  Otherwise start at left of diamond.
;;Similarly for right-junction.  We don't include the left and right junctions themselves.
(defun diamond-destinations (diamond left-junction right-junction)
   (mirror-images
   (unless (rejoining-junction-p left-junction) (setq left-junction nil))) ;Other types of junction are not interesting to us
  (let* ((left (diamond-left diamond))	;Left and right unless modified by junctions
	 (right (diamond-right diamond))
	 start end
	 (amin 0.0) (amax 1.0)		;Valid range of diamond
	 (bmin 0.0) (bmax 1.0)
	 (treatBH nil)
	 )
    (when left-junction			;Junctions restrict the valid range of the diamond
      (setq amax (rejoining-junction-a left-junction)
	    bmin (rejoining-junction-b left-junction)
	    left (diamond-position diamond :a amax :b bmin))) ;New left corner position
    (when right-junction
      (setq amin (rejoining-junction-a right-junction)
	    bmax (rejoining-junction-b right-junction)
	    right (diamond-position diamond :a amin :b bmax))) ;New right corner position

       
    (when (or (eq (diamond-e diamond) :BH)
	      (eq (diamond-e diamond) :BHdeleted))
      (when (and left-junction (null right-junction))
	(setf bmax amax)
	(setf amin bmin)
	(setf right (diamond-position diamond :a (/ (+ amax amin) 2.0) :b (/ (+ amax amin) 2.0)))
	(setf treatBH t))
      (when (and (null left-junction) (null right-junction))
	(setf right (diamond-position diamond :a 0.5 :b 0.5))
	(setf treatBH t)))

    (when (or (eq (diamond-w diamond) :BH)
	      (eq (diamond-w diamond) :BHdeleted))
      (when (and (null left-junction) right-junction)
	(setf amax bmax)
	(setf bmin amin)
	(setf left (diamond-position diamond :a (/ (+ amax amin) 2.0) :b (/ (+ amax amin) 2.0)))
	(setf treatBH t))
      (when (and (null left-junction) (null right-junction))
	(setf left (diamond-position diamond :a 0.5 :b 0.5))
	(setf treatBH t)))

       
    (setq start (diamond-position diamond :a amin :b bmin)       ;Now that range has been determined, compute start and end
	  end (diamond-position diamond :a amax :b bmax))
    
    (when (point-successors start t)
     (error "The starting point of the (possibly reduced) diamond is already past some successors"))
    (let ((tries (point-successors end nil))) ;Which do we need to consider?
      (unless tries
	(error "~S has no valid successors" diamond))
      (if (cdr tries)			 ;More than one possible successor
	  (let* ((successor (find-first-successor start left end left-junction))
		 (result (list successor))
		 (intersections		;Get set of points where successor changes
		  (loop for rest on tries ;Loop over pairs of possible successors
			for first = (car rest)
			nconc (loop for second in (cdr rest)
				    ;;Look for intersection between the diamond and the place where
				    ;;we touch the 2 successors
				    nconc (loop for successor-intersection-element in ;;allow for more than one intersection  
						(if treatBH ;;when bh diamond check successors using successor-intersection-bh 
						    (successor-intersection-bh diamond amin amax bmin bmax
									       start (diamond-position diamond :a amax :b bmin) (diamond-position diamond :a amin :b bmax) end first second)
						  (successor-intersection diamond amin amax bmin bmax
									  start left right end first second)
						  )
						when successor-intersection-element
						collect it))))
		 (last-junction right-junction))
	    (setq intersections (sort intersections #'successor-left-p)) ;Sort from left to right
	    (loop for intersection in intersections ;Step through intersections making list of successors and junctions
		  for next = (cond ((= (successor-intersection-successor-1 intersection) successor)
				    (successor-intersection-successor-2 intersection))
				   ((= (successor-intersection-successor-2 intersection) successor)
				    (successor-intersection-successor-1 intersection))
				   (t
				    (error "Can't walk through successors in ~S" diamond)))
		  do (setq last-junction (create-junction diamond successor next 
							  (successor-intersection-a intersection)
							  (successor-intersection-b intersection)))
		  do (push last-junction result)
		  do (push (setq successor next) result)
		  finally (when right-junction ;Must agree with junction if any
			    (unless (= successor (rejoining-junction-right-direction right-junction))
			      (warn "Last successor does not agree with final junction.  Creating fake junction")
			      (unless last-junction (error "No previous junction"))
			      ;;Put it midway between last one and this one
			      (push (create-junction diamond successor (rejoining-junction-right-direction right-junction)
						     (/ (+ (rejoining-junction-a last-junction) (rejoining-junction-a right-junction)) 2)
						     (/ (+ (rejoining-junction-b last-junction) (rejoining-junction-b right-junction)) 2))
				    result)
			      (push (rejoining-junction-right-direction right-junction) result))))
	    (nreverse result))
	tries)	;Only one successor (usual case): just return it
      )))






(defun successor-left-p (intersection-1 intersection-2)
  (cond ((< (successor-intersection-b intersection-1) (- (successor-intersection-b intersection-2) fudge-ijkl))
	 t)				;#1 clearly to left
	((> (successor-intersection-b intersection-1) (+ (successor-intersection-b intersection-2) fudge-ijkl))
	 nil)				;clearly to right
	(t				;can't tell
	 (error "Intersections too close to sort"))))

;;Find the first successor for this diamond by looking along the left edge
(defun find-first-successor (start left end left-junction)
  (if left-junction			;Started with junction?
      (rejoining-junction-left-direction left-junction)	;Then we must pass it across the remaining surface
    (let ((tries (point-successors left t)))
      (if tries				;crossed before left?
	  (setq end left)		;OK.  End of segment is LEFT
	(setq start left		;Try northwest edge
	      tries (point-successors end nil)))
      (unless tries (error "Couldn't find first successor"))
      (if (cdr tries)			;more than one?
	  (let ((togo (make-togo start)) ;What do we need in each of the 5 directions to reach surface?
		(path (make-pq5 (4vector- end start)))) ;What direction are we going?
	    (loop with best = nil
		  with best-successor
		  for successor in tries
		  for distance = (/ (aref togo successor) (aref path successor)) ;How far along edge to this boundary?
		  when (or (null best) (< distance best))
		  do (setq best distance best-successor successor)
		  finally (return best-successor)))
	(car tries))			;Only one successor: return it
      )))

(defun find-first-successor-bh (left end)
  (let ((tries (point-successors end t)))
    (unless tries (error "Couldn't find first successor"))
    (if (cdr tries)                   ;more than one?
	(let ((togo (make-togo left)) ;What do we need in each of the 5 directions to reach surface?
	      (path (make-pq5 (4vector- end left)))) ;What direction are we going?
	  (loop with best = nil
		with best-successor
		for successor in tries
		for distance = (/ (aref togo successor) (aref path successor)) ;How far along edge to this boundary?
		when (or (null best) (< distance best))
		do (setq best distance best-successor successor)
		finally (return best-successor)))
      (car tries))                    ;Only one successor: return it 
      ))
    
(defvar successor-solution-accuracy 1e-10)

(defstruct (5vector
	    (:type (vector coordinate))
	    (:constructor nil)
	    (:copier nil)
	    (:include 4vector))
  total)

	    
(deftype 5vector () '(simple-array coordinate (5)))

(defun construct-5vector ()
  (make-array 5 :element-type 'coordinate))

(defvar 5vectors (make-resource :name "5vector" :constructor #'construct-5vector :max 100))

(defmacro make-5vector (&optional x y z time fifth)
  (when x (unless fifth (error "Must supply all slots or none")))
  `(let ((result (allocate 5vectors)))
     ,@(when x
	 `((setf (5vector-x result) ,x
		 (5vector-y result) ,y
		 (5vector-z result) ,z
		 (5vector-t result) ,time
		 (5vector-t result) ,fifth)))
     (the 5vector result)))

(defun make-pq5 (pq)
  (let ((pqi (xtoi pq)))
    (if *dump-time*
	(let ((pq5 (make-5vector)))
	  (dotimes (i 4)
	    (setf (aref pq5 i) (aref pqi i))) ;Advance rate in each ijkl coord
	  (setf (aref pq5 dump-destination) (4vector-t pq)) ;Advance rate in actual time
	  pq5)
      pqi)))

;;Vector from start point to IJKL all 1.0
(defun make-togo (start)
  (let ((togo (4vector- unit-4vector (xtoi start))))
    (if *dump-time*
	(let ((togo5 (make-5vector)))
	  (dotimes (i 4)
	    (setf (aref togo5 i) (aref togo i)))
	  ;;Actual time from start to dump time
	  (setf (aref togo5 dump-destination) (- (local-time *dump-time*) (4vector-t start)))
	  togo5)
      togo)))
	  
;;PiQj - PjQi
(defun vector-slot-determinant-parts (vector1 vector2 slot1 slot2)
  (values (* (aref vector1 slot1) (aref vector2 slot2))
	  (* (aref vector2 slot1) (aref vector1 slot2))))

(defun vector-slot-determinant (vector1 vector2 slot1 slot2)
  (multiple-value-bind (x y) (vector-slot-determinant-parts vector1 vector2 slot1 slot2)
    (- x y)))

;;Return determinat, or NIL if it is too small
(defun vector-slot-determinant-check (vector1 vector2 slot1 slot2)
  (multiple-value-bind (x y) (vector-slot-determinant-parts vector1 vector2 slot1 slot2)
    (let ((det (- x y)))
      (and (> (abs det) (* successor-solution-accuracy (max (abs x) (abs y)))) ;Not too much cancellation in subtraction
	   det))))

;;See if the diamond crosses the point where we touch both successors
;;The difference with successor-intersection possible-ab-junctions has to be computed using "solve-for-junctions" 
(defun successor-intersection-bh (diamond amin amax bmin bmax start left right end first second)
  (let*((togo (make-togo start)) ;Distance in IJKL coordinates to each surface 
        (p (make-pq5 (4vector- left start)))
        (q (make-pq5 (4vector- right start)))
        (pqdet (vector-slot-determinant-check p q first second)) ;Get determinant, unless it is too small
        (possible-ab-junctions)
        (successor-intersections))
    (if (eq *era* :flat)
	(when pqdet                     ;If determinant too small, give up                                                                                                                                 
          (setq possible-ab-junctions (list (list (/ (vector-slot-determinant togo q first second) pqdet)
                                                  (/ (vector-slot-determinant p togo first second) pqdet)))))
      (let ((delta (make-pq5 (4vector- (4vector+ left right) (4vector+ start end))))) ;end - flat-end                                                                                                      
        (setq possible-ab-junctions
              (solve-for-junctions p q delta togo first second)))) ;Solve the 2 quadratics for a and b simultaneously
     (loop for ab-list in possible-ab-junctions ;allow for more than one junction                                                                                                                          
          when ab-list
          do (let* ((a (first ab-list))
                    (b (second ab-list)))
               (when (and (<= 0.0 a 1.0) ;;Valid range of a and b?                                                                                                                                         
                          (<= 0.0 b 1.0)) ;;fudge factor?                                                                                                                                                  
                 ;;a and b above are relative to the subregion.  Convert to overall coordinates for the diamond                                                                                            
		 (setq a (+ (* amin (- 1.0 a)) (* amax a))
                       b (+ (* bmin (- 1.0 b)) (* bmax b)))
                 (let*((x-position (diamond-position diamond :a a :b b))
                       (i-position (xtoi x-position)))
                   (when (loop for index below (if *dump-time* 5 4) ;Check that this is really a solution                                                                                                  
                               for value = (if (= index dump-destination) (4vector-t x-position) ;actual time                                                                                              
                                             (aref i-position index))
                               for limit = (if (= index dump-destination) (local-time *dump-time*) 1.0)
                               always (if (or (= index first) (= index second)) ;Given index?                                                                                                              
                                          ;;Not really the right fudge factor for time                                                                                                                     
                                          (fudge= value limit successor-solution-accuracy) ;Actually solves equations?                                                                                     
                                        (< value limit) ;fudge factor?                                                                                                                                     
                                        ))
		     (when (or (and (eq (diamond-e diamond) :BH) (>= a b))
			       (and (eq (diamond-w diamond) :BH) (>= b a))
			       (and (eq (diamond-e diamond) :BHdeleted) (>= a b))
                               (and (eq (diamond-w diamond) :BHdeleted) (>= b a)))
		       (push (make-successor-intersection :a a :b b :successor-1 first :successor-2 second
							  :ijkl i-position) successor-intersections)
		       )
		   )))))
    successor-intersections))

;;See if the diamond crosses the point where we touch both successors
;; That happens when there's a place in the diamond with 2 coordinates = 1.0 and the other 2 < 1.0
(defun successor-intersection (diamond amin amax bmin bmax start left right end first second)
  (let*((togo (make-togo start)) ;Distance in IJKL coordinates to each surface
	(p (make-pq5 (4vector- left start)))
	(q (make-pq5 (4vector- right start)))
	(pqdet (vector-slot-determinant-check p q first second)) ;Get determinant, unless it is too small
	(possible-ab-junctions)
	(successor-intersections))
    (if (eq *era* :flat)
	(when pqdet			;If determinant too small, give up
	  (setq possible-ab-junctions (list (list (/ (vector-slot-determinant togo q first second) pqdet) 
						  (/ (vector-slot-determinant p togo first second) pqdet)))))
      (let ((delta (make-pq5 (4vector- (4vector+ left right) (4vector+ start end))))) ;end - flat-end
	(setq possible-ab-junctions
	      (solve-for-junctions p q delta togo first second)))) ;Solve the 2 quadratics for a and b simultaneously
    (loop for ab-list in possible-ab-junctions ;allow for more than one junction
	  when ab-list
	  do (let* ((a (first ab-list))
		    (b (second ab-list)))
	       (when (and (<= 0.0 a 1.0) ;;Valid range of a and b?
			  (<= 0.0 b 1.0)) ;;fudge factor?
		 ;;a and b above are relative to the subregion.  Convert to overall coordinates for the diamond
		 (setq a (+ (* amin (- 1.0 a)) (* amax a))
		       b (+ (* bmin (- 1.0 b)) (* bmax b)))
		 (let*((x-position (diamond-position diamond :a a :b b))
		       (i-position (xtoi x-position)))
		   (when (loop for index below (if *dump-time* 5 4) ;Check that this is really a solution
			       for value = (if (= index dump-destination) (4vector-t x-position) ;actual time
					     (aref i-position index))
			       for limit = (if (= index dump-destination) (local-time *dump-time*) 1.0)
			       always (if (or (= index first) (= index second)) ;Given index?
					  ;;Not really the right fudge factor for time
					  (fudge= value limit successor-solution-accuracy) ;Actually solves equations?
					(< value limit) ;fudge factor?
					))
		     (push (make-successor-intersection :a a :b b :successor-1 first :successor-2 second
							:ijkl i-position) successor-intersections)
		     )))))
    successor-intersections))
  
;;This function solves the 2 quadratics of the form a*p + b*q + a*b*delta = togo
;;and produces a list of the (a b) solution. Some or both could be nil
;; in which case they should be discarded later on.
(defun solve-for-junctions (p q delta togo first second)
  (let*((qd (vector-slot-determinant q delta first second))
	(pd (vector-slot-determinant p delta first second))
	(dt (vector-slot-determinant delta togo first second))
	(tq (vector-slot-determinant togo q first second))
	(tp (vector-slot-determinant togo p first second))
	(pq (vector-slot-determinant p q first second))
	(as (quadratic-solver pd (- dt pq) tq)) ;;Solve the quadratic for a
	(bs (quadratic-solver qd (+ dt pq) tp)) ;;Solve the quadratic for b
	;; we have to group the solutions for a and b in the correct order
	(a1 (first as))   
	(b1 (second bs))
	(a2 (second as)) 
	(b2 (first bs)))
    (list     
     (if  (and a1 b1) ;;Either a or b could be nil in which case we should not consider it
	 (list a1 b1) 
       nil)
     (if (and a2 b2)
	 (list a2 b2)
       nil))))

;; This function returns the possible 2 solutions of a generic quadratic
;; it tries to handle all the cases. It uses the Numerical recepies method
;; to evaluate the solutions.
(defun quadratic-solver (a b c)
  (let*((bsign (signum b))
	(discriminant (- (* b b) (* 4.0 a c)))
	(q (/ (- 0.0 (+ b (* bsign (sqrt discriminant)))) 2.0)))
    (cond 
     ((= a 0.0) ;; Takes care of the case where the quadratic is in fact linear
      (cond     ;; The order is important to pair them later on with the other solution for a
       ((= bsign 1) (list (/ (- 0.0 c) b) nil)) 
       ((= bsign -1) (list nil (/ (- 0.0 c) b)))
       (t (list nil nil))) ;; if b is zero too then there is no solution
       )
      ((= b 0.0) ;;Takes care of the case where there is no linear term
       (cond 
	((= c 0.0) (list 0.0 0.0))
	((minusp (* a c)) (list (- 0.0 (sqrt (- 0.0 (/ c a)))) (sqrt (- 0.0 (/ c a)))))
	(t (list nil nil))
	)
       )
      ((>= discriminant 0.0)  ;; Generic case
       (if (= bsign 1)
	   (list (/ c q) (/ q a))	
	 (list (/ q a) (/ c q))))
      (t (list nil nil))
      )))
