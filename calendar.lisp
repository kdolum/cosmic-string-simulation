(in-package "CL-USER")

;;Calendar implementation of priority queues.
;;Copyright (c) 2007 Ken Olum
;;See R. Brown, "Calendar queues: a fast 0(1) priority queue implementation
;;for the simulation event set problem", Communications of the ACM, 31, 1220 
 
(defstruct (calendar
	    (:constructor construct-calendar)
	    (:print-object print-calendar))
  (name nil)
  (days 0 :type fixnum)			;number of days in the calendar
  (current-day 0 :type fixnum)	    ;Day which had last-dequeued event
  (current-time 0.0d0 :type double-float) ;Time of last-dequeued event.
  (this-year 0 :type fixnum)		;The number of this year
  (day-length 0.0d0 :type double-float)
  (event-count 0 :type fixnum)		;Number of events stored
  (pages #() :type simple-vector)			;An array of length DAYS
  ;;Whether to allow queuing before current time.  T = always, NIL = never, number = by that amount at most
  (backups-allowed nil :type (or (member nil t) double-float))
  )

(defun print-calendar (calendar stream)
  (print-unreadable-object (calendar stream :identity t)
    (format stream "CALENDAR ~S ~D events"
	    (calendar-name calendar) (calendar-event-count calendar))))

;;PAGES is an array indexed by time/day-length given a list of
;;(time . thing) with the earliest time first.

(declaim (inline calendar-year calendar-day calendar-empty-p))

(defun calendar-empty-p (calendar)
  (zerop (calendar-event-count calendar)))

;;Compute the number of the day where a given time would be stored.  
(defun calendar-day (calendar time)
  (declare (optimize speed)
	   (double-float time))
  ;;First find integer day number on which to store this event.
  ;;Subsequent calls to this function will always get the same integer here,
  ;;so we never lose the event by round-off error.
  (let ((day (fixnum-floor (/ time (calendar-day-length calendar)))))
    (locally (declare (optimize (safety 0))) ;Don't type-check THE
      (the fixnum (mod day (calendar-days calendar))))))

;;Compute the year to which a given time belongs
(defun calendar-year (calendar time)
  (declare (optimize speed)
	   (double-float time))
  (let ((day (fixnum-floor (/ time (calendar-day-length calendar)))))
    (declare (fixnum day))
    (values (floor day (calendar-days calendar)))))

;;Make calendar for storing events with a period of time given by year
(defun make-calendar (year-length &key (days 4) name backups-allowed (current-time 0.0))
  (let ((calendar (construct-calendar
		   :days days
		   :day-length (/ year-length days)
		   :name name
		   :backups-allowed backups-allowed
		   :current-time current-time
		   :pages (make-array days :initial-element nil))))
    ;;Set this year and day to correspond to current time
    (setf (calendar-this-year calendar) (calendar-year calendar current-time)
	  (calendar-current-day calendar) (calendar-day calendar current-time))
    calendar))

;;Add event to calendar
(defun calendar-add (calendar time thing)
  (declare (optimize speed)
	   (double-float time))
  (when (< time (calendar-current-time calendar))
    (let ((allowed (calendar-backups-allowed calendar)))
      (typecase allowed
	(null				;No backups?
	 (error "Trying to add a past event to a calendar"))
	(double-float			;Limited backups?
	 (unless (<= (- (calendar-current-time calendar) time) allowed) ;Within allowable step
	   (without-compiler-notes	;Don't care about consing for error
	    (error "Trying to add an event ~E in the past, (allowed is ~E)" (- (calendar-current-time calendar) time) allowed))))))
    (return-from calendar-add (calendar-add-past calendar time thing)))
  (let* ((index (calendar-day calendar time))
	 ;; List of the old events in chronological order
	 (later (aref (calendar-pages calendar) index))
	 (earlier nil))
    ;;We need to find any events earlier than the one we're trying to insert
    (loop while later
	  for event = (car later)	;First possibly later event
	  for found-time double-float = (car event)	;its time
	  until (< time found-time)	;Actually is later: exit
	  when (and (eq (cdr event) thing) ;duplicate event?
		    (= time found-time))
	  do (error "~S is already in ~S at time ~S" thing calendar time)
	  ;;Found earlier event, remember it
	  do (setq earlier later)	;Mark it as earlier
	  do (pop later)		;and remove from later list
	  )
    ;;Now later is the list of events (possibly NIL) later than ours.
    ;;while earlier is the cons (if any) which stores the last earlier event
    (push (cons time thing) later)	;Now LATER has our event and all later
    (if earlier (setf (cdr earlier) later) ;link it in if there is a cons
      (setf (aref (calendar-pages calendar) index) later) ;store if not
      )
    ;;Count event, maybe resize
    (incf-calendar-count calendar)))

;;Add an event to a calendar before the current time
(defun calendar-add-past (calendar time thing)
  (let ((index (calendar-day calendar time)))
    (push (cons time thing)		;Goes at front of list, since it is earliest
	  (aref (calendar-pages calendar) index))
    ;;Back up data structures to point to new first event
    (setf (calendar-this-year calendar) (calendar-year calendar time)
	  (calendar-current-day calendar) index
	  (calendar-current-time calendar) time))
  (incf-calendar-count calendar))

;;Delete event from calendar.
(defun calendar-delete (calendar time thing)
  (let* ((index (calendar-day calendar time))
	 (pages (calendar-pages calendar)))
    ;;This is like DELETE, except it fails if the event isn't there
    (loop for last = nil then cons	;Previous cons
	  for cons on (aref pages index)
	  for (this-time . this-thing) = (car cons)
	  when (and (= this-time time) (eq this-thing thing)) ;Found it?
	  ;;Delete it and return
	  return (if last			;have previous?
		     (setf (cdr last) (cdr cons)) ;get rid of current
		   (setf (aref pages index) (cdr cons))) ;delete first
	  finally (error "~S was not found in calendar ~S at time ~S"
			 thing calendar time)))
  (decf-calendar-count calendar))

(defun object-in-calendar-p (calendar time thing)
  (loop for cons in (aref (calendar-pages calendar) (calendar-day calendar time))
	thereis (and (= (car cons) time) (eq (cdr cons) thing))
	))

;;Find next event to process.  Remove it unless PEEK is set.
;;In any case, it is too late after this to enqueue events earlier than
;;the one that we return.
;;Returns TIME and THING
(defun calendar-next (calendar &optional peek)
  (let* ((current-day (calendar-current-day calendar))
	 (pages (calendar-pages calendar))
	 (year (calendar-this-year calendar))
	 (days (calendar-days calendar))
	 (index current-day))
    (loop
     (let ((events (aref pages index)))
       (when events
	 (let ((event (car events)))
	   (when (= (calendar-year calendar (car event)) year) ;Right year?
	     ;;Found event.  Remove and return it
	     (setf (calendar-current-day calendar) index
		   (calendar-current-time calendar) (car event))
	     (unless peek
	       (setf (aref pages index) (cdr events))
	       (decf-calendar-count calendar))
	     (return (values (car event) (cdr event))) ;time, thing
	     ))))
     (incf index)
     (when (= index days)
       (setq index 0)			;Wrap around array
       (incf year)			;Go to subsequent year
       (setf (calendar-this-year calendar) year))
     (when (= index current-day) ;Did whole year without finding anything
       (return (calendar-next-direct calendar peek))) ;Use linear search
     )))

;;When we scan an entire year without finding anything, we come here to
;;do it by direct search
(defun calendar-next-direct (calendar peek)
  (loop with pages = (calendar-pages calendar)
	with first-index = nil
	with first-time
	for index from 0 below (calendar-days calendar)
	as events = (aref pages index)
	when events
	  do (let ((this-time (caar events)))
	       (when (or (null first-index) ;first seen
			 (< this-time first-time)) ;this event sooner?
		 (setq first-index index
		       first-time this-time)))
	finally
	(unless first-index
	  (error "Attempted to dequeue from an empty calendar"))
	(let* ((events (aref pages first-index))
	       (event (car events)))
	  ;;Point calendar to the place that we found
	  (setf (calendar-this-year calendar) (calendar-year calendar (car event)))
	  (setf (calendar-current-day calendar) first-index)
	  (setf (calendar-current-time calendar) (car event))
	  (unless peek
	    (setf (aref pages first-index) (cdr events)) ;remove event
	    (decf-calendar-count calendar))
	  (return (values (car event) (cdr event))))))

;;Increment count of events and maybe resize.  You should not have anything
;;cached when you do this, because we might remake everything.
(defun incf-calendar-count (calendar)
  (when (>= (incf (calendar-event-count calendar))
	    (* 2 (calendar-days calendar)))
    (resize-calendar calendar)))

;;Decrement count of events.  Resize if half the days are empty, but never
;;below size 4.
(defun decf-calendar-count (calendar)
  (when (<= 4
	    (decf (calendar-event-count calendar))
	    (truncate (calendar-days calendar) 2))
    (resize-calendar calendar)))

;;Remake the calendar with the number of days equal to the number of events
;;At the moment, we're not adjusting the year length.
;;We retain the time of the last-dequeued event, so that we can add new
;;events anytime after this.
(defun resize-calendar (calendar)
  (let ((old-pages (calendar-pages calendar))
	(old-time (calendar-current-time calendar))
	(old-days (calendar-days calendar))
	(old-day-length (calendar-day-length calendar))
	(new-days (calendar-event-count calendar))) ;New size
    (setf (calendar-days calendar) new-days)
    (setf (calendar-day-length calendar) ;New length of day
	  (/ (* old-days old-day-length) new-days))
    (setf (calendar-pages calendar)	;New data array
	  (make-array new-days :initial-element nil))
    ;;Set pointers right for old last-event time.
    (setf (calendar-this-year calendar) (calendar-year calendar old-time)
	  (calendar-current-day calendar) (calendar-day calendar old-time))
    ;;Now the calendar is empty.  Add all the old events.
    (setf (calendar-event-count calendar) 0)
    (loop for index below (length old-pages)
	  do (loop for (time . thing) in (aref old-pages index)
		   do (calendar-add calendar time thing)))
    (fill old-pages 0)			;Discourage conservative GC problems by clearing old array
    ))

;;Call function with time and thing for all events in the calendar
(defun map-calendar (function calendar)
  (loop for index below (calendar-days calendar)
	do (loop for (time . thing) in (aref (calendar-pages calendar) index)
		 do (funcall function time thing))))
