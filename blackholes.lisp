#||
============================================================
background spacetime coordinates X^mu as function
of string worldsheet coordiantes tau and sigma
   X^mu(tau,sigma) = ( a^mu(v) + b^mu(u) ) / 2
where we work with the temporal gauge
   tau = X^0 = t
and lightcone coordinates
   u =  tau + sigma
   v =  tau - sigma
therefore, X^0 could be obtained by
   X^0 = ( u + v ) / 2 = t
namely,
   a^0(v) = v
   b^0(u) = u
the spatial part X^i(t,sigma) from a^i(v) and b^i(u)
satisfies following EOM and constraints
   Xddot - Xpprm = 0
   Xdot * Xprm = ( b'(u)^2 - a'(v)^2 ) / 4 = 0
   Xdot^2 + Xprm^2 = ( a'(v)^2 + b'(u)^2 ) / 2 = 1
from
   Xdot = ( b'(u) + a'(v) ) / 2
   Xprm = ( b'(u) - a'(v) ) / 2
   Xddot = ( b''(u) + a''(v) ) / 2
   Xpprm = ( b''(u) + a''(v) ) / 2
if
   a'(v)^2 = b'(u)^2 = 1
============================================================
||#


;given a'(v) and b'(u)



(defvar *A* 0.5)
(defvar *w* 0.1)

(defvar *astart*
  (3to4vector (make-3vector 2.0d0
			    2.0d0
			    2.0d0)
	      10.0d0)
)

(defvar *bstart*
  (3to4vector (make-3vector 2.0d0
			    2.0d0
			    2.0d0)
	      10.0d0)
)


(defun aprime0 (v)
  (let ((normalizer (sqrt (+ 1.0d0 (* *A* *A* (sin (* *w* v)) (sin (* *w* v))))))) 
       (make-4vector (/ (sin v) normalizer) 
                     (/ (* -1.0d0 (cos v)) normalizer) 
                     (/ (* *A* (sin (* *w* v) )) normalizer)
                     1.0d0)
  )
)
        
(defun bprime0 (u)
  (let ((normalizer (sqrt (+ 1.0d0 (* *A* *A* (sin (* *w* u)) (sin (* *w* u)))))))
       (make-4vector (/ (cos u) normalizer)
                     (/ (sin u) normalizer)
                     (/ (* *A* (sin (* *w* u) )) normalizer)
                     1.0d0)
  )
)
 

;compute the diamond-right of a diamond
;similar to compute-diamond-end-flat
(mirror-images

(defun compute-diamond-right-flat (diamond)
       (4vector- (4vector+ (diamond-start diamond) (diamond-end diamond))
                 (diamond-left diamond))
)

)


;compute xvector from given aprime-function and bprime-function
;along positive/negative v and u directions
; a(v0+dv)=a(v0)+a'(v0)dv

(defun avector-aprime (input-avector aprime-function input-v dv)
;       (if (> dv 0)
           (4vector+ input-avector (4vector-scale (funcall aprime-function input-v) dv))
;           (4vector+ input-avector (4vector-scale (funcall aprime-function (+ input-v dv)) dv))
;       )
)

(defun bvector-bprime (input-bvector bprime-function input-u du)
;       (if (> du 0)
           (4vector+ input-bvector (4vector-scale (funcall bprime-function input-u) du))
;           (4vector+ input-bvector (4vector-scale (funcall bprime-function (+ input-u du)) du))
;       )
)


;numerical integration of a'(v) from v0 to v=v0+ndv*dv to get a(v0+ndv*dv)
(defun avector-from-aprime (aprime-function a0 v0 ndv dv)
  (let ((partial-sum (make-zero-4vector)))
       (when (> ndv 1)
             (loop for i from 1 below ndv
                do (setf partial-sum (4vector+ partial-sum (funcall aprime-function (+ v0 (* i dv)))))
             )
       )
       (4vector+ (if (= ndv 0) (make-zero-4vector)
                     (4vector-scale (4vector+ (funcall aprime-function v0)
                                              (funcall aprime-function (+ v0 (* ndv dv)))
                                    )
                                    (* dv 0.5d0)
                     )
                 )
                 (4vector-scale partial-sum dv)
                 a0
       )
  )
)

;numerical integration of b'(u) from u0 to u=u0+ndu*du to get b(u0+ndu*du)
(defun bvector-from-bprime (bprime-function b0 u0 ndu du)
  (let ((partial-sum (make-zero-4vector)))
       (when (> ndu 1)
             (loop for i from 1 below ndu
                do (setf partial-sum (4vector+ partial-sum (funcall bprime-function (+ u0 (* i du)))))
             )
       )
       (4vector+ (if (= ndu 0) (make-zero-4vector)
                     (4vector-scale (4vector+ (funcall bprime-function u0)
                                              (funcall bprime-function (+ u0 (* ndu du)))
                                    )
                                    (* du 0.5d0)
                     )
                 )
                 (4vector-scale partial-sum du)
                 b0
       )
  )
)



(defun xvector (avector bvector)
      (4vector-scale (4vector+ avector bvector) 0.5d0)
)


;compute the diamond on the southeast (dv<0) of input-diamond with (ustart, vstart) indicated
(defun compute-diamond-se-aprime (input-diamond bstart astart aprime-function vstart dv)
  (let ((output-diamond (make-diamond)))
       (setf (diamond-end   output-diamond) (diamond-right input-diamond))
       (setf (diamond-left  output-diamond) (diamond-start input-diamond))
       (setf (diamond-start output-diamond) (xvector (avector-aprime astart aprime-function vstart dv) bstart))
       (setf (diamond-right output-diamond) (compute-diamond-right-flat output-diamond))
       (setf (diamond-se input-diamond) output-diamond)
       (setf (diamond-nw output-diamond) input-diamond)
       output-diamond
  )
)

;compute the diamond on the southeast (dv<0) of input-diamond, with (ndvstart,ndustart) indicated for output-diamond
(defun compute-diamond-se-from-abprime (input-diamond aprime-function bprime-function a0 b0 v0 u0 ndvstart ndustart dv du)
  (let ((output-diamond (make-diamond)))
       (setf (diamond-end   output-diamond) (diamond-right input-diamond))
       (setf (diamond-left  output-diamond) (diamond-start input-diamond))
       (setf (diamond-start output-diamond) (xvector (avector-from-aprime aprime-function a0 v0 ndvstart dv)
                                                     (bvector-from-bprime bprime-function b0 u0 ndustart du)))
       (setf (diamond-right output-diamond) (compute-diamond-right-flat output-diamond))
       (setf (diamond-se input-diamond) output-diamond)
       (setf (diamond-nw output-diamond) input-diamond)
       output-diamond
  )
)

;compute the diamond on the northeast (du>0) of input-diamond with (uend, vend) indicated
(defun compute-diamond-ne-bprime (input-diamond aend bend bprime-function uend du)
  (let ((output-diamond (make-diamond)))
       (setf (diamond-start output-diamond) (diamond-right input-diamond))
       (setf (diamond-left  output-diamond) (diamond-end   input-diamond))
       (setf (diamond-end   output-diamond) (xvector aend (bvector-bprime bend bprime-function uend du)))
       (setf (diamond-right output-diamond) (compute-diamond-right-flat output-diamond))
       (setf (diamond-ne input-diamond) output-diamond)
       (setf (diamond-sw output-diamond) input-diamond)
       output-diamond
  )
)

;compute the diamond on the northeast (du>0) of input-diamond, with (ndvstart,ndustart) indicated for output-diamond
(defun compute-diamond-ne-from-abprime (input-diamond aprime-function bprime-function a0 b0 v0 u0 ndvstart ndustart dv du)
  (let ((output-diamond (make-diamond)))
       (setf (diamond-start output-diamond) (diamond-right input-diamond))
       (setf (diamond-left  output-diamond) (diamond-end   input-diamond))
       (setf (diamond-end   output-diamond) (xvector (avector-from-aprime aprime-function a0 v0 (- ndvstart 1) dv)
                                                     (bvector-from-bprime bprime-function b0 u0 (+ ndustart 1) du)))
       (setf (diamond-right output-diamond) (compute-diamond-right-flat output-diamond))
       (setf (diamond-ne input-diamond) output-diamond)
       (setf (diamond-sw output-diamond) input-diamond)
       output-diamond
  )
)




;create a string from given aprime and bprime
;if sigma-length equals sigma-period, then the two ends of string are identical
;if closep is true, the the two ends of string are connected with each other
(defun create-string-abprime (astart bstart aprime-function bprime-function diamond-number sigma-length closep)
 (let* ((initial-diamond (make-diamond))
        (previous-diamond (make-diamond))
        (current-diamond (make-diamond))
        (dsigma-number (if (evenp diamond-number) (/ diamond-number 2) (/ (- diamond-number 1) 2)))
        (dsigma (/ sigma-length dsigma-number))
        (unext 0.0d0) 
        (vnext 0.0d0)
        (anext astart)
        (bnext bstart)
       )
       ;set-up the initial-diamond
       (setf (diamond-start initial-diamond) (xvector anext bnext))
       (setf (diamond-left  initial-diamond) (xvector (avector-aprime anext aprime-function vnext dsigma) bnext))
       (setf (diamond-right initial-diamond) (xvector anext (bvector-bprime bnext bprime-function unext dsigma)))
       (setf (diamond-end   initial-diamond) (compute-diamond-end-flat initial-diamond))
       ;prepare to go to next diamond
       (setf previous-diamond initial-diamond)
       (loop for i below (- diamond-number 1) ;i=0,1,2,...,diamond-number - 2
          do (when (evenp i)
                   (setf current-diamond (compute-diamond-se-aprime previous-diamond bnext anext aprime-function vnext (* -1.0d0 dsigma)))
                   (setf bnext (bvector-bprime bnext bprime-function unext dsigma))
                   (setf unext (+ unext dsigma))
             )
          do (when (oddp i)
                   (setf current-diamond (compute-diamond-ne-bprime previous-diamond anext bnext bprime-function unext dsigma))
                   (setf anext (avector-aprime anext aprime-function vnext (* -1.0d0 dsigma)))
                   (setf vnext (- vnext dsigma))
             )
          do (setf previous-diamond current-diamond)
       )
       (when closep
             (setf (diamond-ne current-diamond) initial-diamond)
             (setf (diamond-sw initial-diamond) current-diamond)
       )
       current-diamond
  )
)      

(defun create-string-from-abprime (a0 b0 aprime-function bprime-function diamond-number sigma-length closep)
 (let* ((initial-diamond (make-diamond)) 
        (previous-diamond (make-diamond)) 
        (current-diamond (make-diamond))
        (dsigma-number (if (evenp diamond-number) (/ diamond-number 2) (/ (- diamond-number 1) 2)))
        (dsigma (/ sigma-length dsigma-number)) 
        (v0 0.0d0) 
        (u0 0.0d0) 
        (ndvstart 0) 
        (ndustart 0) 
        (dv (* -1.0d0 dsigma)) 
        (du dsigma)
       )
       (setf (diamond-start initial-diamond) (xvector a0 b0))
       (setf (diamond-left  initial-diamond) (xvector (avector-from-aprime aprime-function a0 v0 (+ ndvstart 1) dsigma) b0))
       (setf (diamond-right initial-diamond) (xvector a0 (bvector-from-bprime bprime-function b0 u0 (+ ndustart 1) dsigma)))
       (setf (diamond-end   initial-diamond) (compute-diamond-end-flat initial-diamond))
       (setf ndvstart (+ ndvstart 1))
       (setf previous-diamond initial-diamond)
       (loop for i from 1 below diamond-number
          do (when (oddp i)
                   (setf current-diamond (compute-diamond-se-from-abprime previous-diamond aprime-function bprime-function a0 b0 v0 u0 ndvstart ndustart dv du))
                   (setf ndustart (+ ndustart 1))
             )
          do (when (evenp i)
                   (setf current-diamond (compute-diamond-ne-from-abprime previous-diamond aprime-function bprime-function a0 b0 v0 u0 ndvstart ndustart dv du))
                   (setf ndvstart (+ ndvstart 1))
             )
          do (setf previous-diamond current-diamond)
       )
       (when (and (evenp diamond-number) closep)
             (setf (diamond-ne current-diamond) initial-diamond)
             (setf (diamond-sw initial-diamond) current-diamond)
       )
       current-diamond
  )
)


;reset the spatial components of diamond-end of a BH-diamond to its diamond-start
(defun set-BH-diamond (BH-diamond)
  (let (  (spatial-start (4to3vector (diamond-start BH-diamond)))  )
       (setf (4vector-x (diamond-end BH-diamond)) (3vector-x spatial-start))
       (setf (4vector-y (diamond-end BH-diamond)) (3vector-y spatial-start))
       (setf (4vector-z (diamond-end BH-diamond)) (3vector-z spatial-start))
       BH-diamond
  )
)


;create a BH-string from given abprime
(defun create-BH-string-abprime (astart bstart aprime-function bprime-function diamond-number sigma-length)
  (let ((return-diamond (make-diamond)) 
;        (current-diamond (create-string-abprime astart bstart aprime-function bprime-function diamond-number sigma-length nil))
        (current-diamond (create-string-from-abprime astart bstart aprime-function bprime-function diamond-number sigma-length nil))
       )
       ;set the diamond with BH on the right
       (setf (diamond-ne current-diamond) :BH)
       (setf (diamond-se current-diamond) :BH)
       (set-BH-diamond current-diamond)
       (setf return-diamond current-diamond)
       ;go to the left end of the string
       (loop while (diamond-w current-diamond)
          do (setf current-diamond (diamond-w current-diamond))
       )
       ;set the diamond with BH on the left
       (setf (diamond-nw current-diamond) :BH)
       (setf (diamond-sw current-diamond) :BH)
       (set-BH-diamond current-diamond)
       return-diamond
  )
)

(defun create-bh-tag (position bh)
  (make-loop-tag :created-position position :bh bh))


(defun update-tag-time (tag time)
  (let ((position (tag-created-position tag)))
    (setf position (3to4vector (4to3vector position) time))
    (setf (tag-created-position tag) position)
    tag))

(defun handle-bh-tags (east predecessor)
  (let* ((east-tag (diamond-tag east))
	 (tag (diamond-tag predecessor))
	 (east-count (diamond-countup east))
	 (count (diamond-countup predecessor))
	 (bh (diamond-bh predecessor))
	 (east-bh (tag-bh east-tag))
	 (east-tag-position (tag-created-position east-tag))
	 (bh-position (globalize-position (diamond-end predecessor)))
	 (new-position (3to4vector (4to3vector east-tag-position) (4vector-t bh-position)))
	 (new-bh (make-blackhole :center (3to4vector (4to3vector bh-position) (4vector-t east-tag-position)) :size 0.0))
	 (periods 1)
	 )
    (cond (east-bh
	   (cond ((< (3vector-distance (blackhole-center bh) (blackhole-center east-bh)) 0.0001)  
		  (cond ((eq (get-tag-handle tag) (get-tag-handle east-tag))
			 (cond ((and (> east-count count) (eq east-count periods))
				(values tag (+ east-count 1) t))
			       ((and (> east-count count) (> east-count periods))
				(values tag (+ east-count 1) t))
			       ((> east-count count)
				(values tag (+ east-count 1) nil))
			       (t
				(values tag count (diamond-inertp predecessor)))))
			(t
			 (setf (tag-last-position east-tag) bh-position)
			 (values east-tag 1 t))))
		 (t
		  (values (create-bh-tag new-position new-bh) 0 nil))))
	  (t
	   (values (create-bh-tag new-position new-bh) 0 nil)))))
		 


(defun record-bh-loop (predecessor east)
  (format t "~% Recording loop ")
  (when *bh-loops-found*                   ;Keeping these?
    (let* ((delete-event (globalize-position (diamond-end predecessor)))
           (delete-time (4vector-t delete-event))
           (tag (diamond-tag east))
           (inert-event (tag-last-position tag)) ;Position of last engulfment, where it was made inert
           (inert-time (4vector-t inert-event)) ;Time it was made inert.
           (ic-time (4vector-t (blackhole-center (tag-bh (diamond-tag predecessor))))) ;Time of event creating the loop
           (4-momentum (standardize-position (4vector- delete-event inert-event))) ;comoving 4-momentum averaged over 1 period
           (velocity (3vector-scale 4-momentum (/ 1 (4vector-t 4-momentum))))
           (speed (3vector-length velocity))
           (pf (/ speed                 ;this is the momentum between inert and delete times
                  (sqrt (- 1.0 (expt speed 2.0)))))
           (energy (loop-length-i ic-time inert-time delete-time))  ;energy is the comoving energy of the loop at time of formation.
           (xi (/ energy ic-time))
           (p-i (/ pf (scale-factor-ratio ic-time (/ (+ delete-time inert-time) 2.0)))))        ;correct for redshifting since formation
      (when (or (null *loop-record-start*) ;Don't record loops at early times
                (> ic-time *loop-record-start*))
        (vector-push-extend xi *bh-loops-found*)
        (vector-push-extend (realpart p-i) *bh-loops-found*)
        (vector-push-extend ic-time *bh-loops-found*)
        (when *log-loop-positions*
          (let ((position (globalize-position (diamond-left east)))) ;The point at which the first :deleted diamond began
            (vector-push-extend (4vector-x position) *bh-loops-found*)
            (vector-push-extend (4vector-y position) *bh-loops-found*)
            (vector-push-extend (4vector-z position) *bh-loops-found*)))
        (when *log-loop-velocities*
          (vector-push-extend (3vector-x velocity) *bh-loops-found*)
          (vector-push-extend (3vector-y velocity) *bh-loops-found*)
          (vector-push-extend (3vector-z velocity) *bh-loops-found*))
        ))))

(defun write-recorded-bh-loops ()
    (format t "Writing ~D bh-loops..." (/ (length *bh-loops-found*) (+ 3
                                                                 (if *log-loop-positions* 3 0))))
    (loop for x across *bh-loops-found*
	  do (write-single x *bh-loops-output*))
    (force-output *bh-loops-output*)
    (format t "done.~%"))

;advance the diamond with BH on the left

(defun bin-bh-loops-files (directory
                        &key (t-bins 24) (x-bins 120) (p-bins 1)
                        (t-min 50.0) t-max
                        (x-min 1d-8) x-max
                        (p-min 1d-8) ;Set very low since we're not using it
                        (p-max 1d3)
                        (scaling t)     ;If NIL, use l instead of x, if :mass use x = m/t instead of E/t
                        (errorbar-split 1) ;Compute error bars by splitting
                        (integrated nil) ;if NIL compute just f(x)dxdt binning, if T compute n(x,t)dx
                        )         ;Submit batch jobs to do work.  # to submit or T
  (unless scaling (error "Only scaling prebinning is currently implemented"))
  (setq directory (merge-pathnames directory batch-root-directory))
  (unless (probe-file directory)
    (error "No such directory ~A" directory))
  (format t "~&Binning loops...") (force-output)
  (let* ((run-info (read-run-info-file directory))
         (*era* (run-info-era run-info)) ;Used for cosmological data
         (total-size (double-float (run-info-total-size run-info)))
         (start (run-info-start-time run-info))
         (initial-time (run-info-start-time run-info))
         ;;The following may be wrong if you didn't actually run until the default ending time
         (t-end (default-simulation-end total-size initial-time nil)) ;FIX
         (t-min (or t-min start))
         (t-max (cond (t-max t-max)     ;last time when a loop can be produced which is relevant
                      (t (+ start (light-crossing-time total-size)))))
         (x-max (or x-max (/ (loop-length-i t-max t-max t-end) t-max 3))) ;FIX
         (prebin-info (make-prebin-info :era *era* :physical nil :scaling scaling
                                        :p-bins p-bins
                                        :p-min p-min
                                        :p-max p-max
                                        :x-bins x-bins :t-bins t-bins :errorbar-split errorbar-split
                                        :x-min x-min :x-bin-size (/ (log (/ x-max x-min)) x-bins)
                                        :t-min t-min
                                        :t-bin-size (/ (log (/ t-max t-min)) t-bins)
                                        :total-size total-size)))
    (write-prebin-info prebin-info directory integrated)
    (let ((end (loop for worker from 0                                          ;Find how many files there are
                     unless (probe-file (worker-subdirectory directory worker)) ;don't care if worker has a file or not
                     return worker)))
      (dotimes (worker end)
	(format t "~D " worker) (force-output)
	(prebin-bh-file directory worker prebin-info integrated))))
  (combine-prebin-bh-files directory :integrated integrated))

(defun prebin-bh-file (directory worker prebin-info integrated)
  (write-prebin-bh-data directory
                     worker
                     (bin-loops-file (bh-loop-spectrum-file directory worker) prebin-info integrated)
                     prebin-info
                     integrated))

(defun write-prebin-bh-data (directory worker data prebin-info integrated)
  (with-group-write-access
   (with-open-file (stream (prebin-bh-data-file directory :integrated integrated :worker worker)
                           :direction :output :element-type '(unsigned-byte 32) :if-does-not-exist :create
                           :if-exists :supersede)
     (let* ((errorbar-split (prebin-info-errorbar-split prebin-info))
            (p-bins (prebin-info-p-bins prebin-info))
            (x-bins (prebin-info-x-bins prebin-info))
            (t-bins (+ (if integrated 1 0) (prebin-info-t-bins prebin-info))))
       (dotimes (split (expt errorbar-split 3))
         (dotimes (t-bin t-bins)
           (dotimes (x-bin x-bins)
             (dotimes (p-bin p-bins)
               (write-single (aref data split t-bin x-bin p-bin) stream)))))))))

(defun combine-prebin-bh-files (directory &key integrated)
  (setq directory (merge-pathnames directory batch-root-directory))
  (combine-prebin-bh-files-1 directory
                          (read-prebin-info directory integrated)
                          integrated))

(defun combine-prebin-bh-files-1 (directory info integrated)
  (let* ((*era* (prebin-info-era info))
         (p-bins (prebin-info-p-bins info))
         (x-bins (prebin-info-x-bins info))
         (t-bins (+ (if integrated 1 0) (prebin-info-t-bins info))) ;One more "bin" if integrated                                                                                                           
         (errorbar-split (prebin-info-errorbar-split info))
         (result (make-array (list (expt errorbar-split 3) t-bins x-bins p-bins)
                             :element-type 'double-float :initial-element 0.0)))
    (declare (fixnum p-bins x-bins t-bins errorbar-split)
             (type (simple-array double-float 4) result)
;;           (optimize (speed 3))                                                                                                                                                                           
             )
    (format t "~&Combining prebin files...") (force-output)
    (loop for worker from 0
          as file = (prebin-bh-data-file directory :worker worker :integrated integrated)
          while (probe-file file)       ;Stop when file doesn't exist.  Even never-started workers should have one.                                                                                         
          do (format t "~D " worker) (force-output)
          do (with-open-file (stream file :element-type '(unsigned-byte 32))
               (dotimes (split (expt errorbar-split 3))
                 (dotimes (t-bin t-bins)
                   (dotimes (x-bin x-bins)
                     (dotimes (p-bin p-bins)
                       (incf (aref result split t-bin x-bin p-bin) (read-single-double stream))))))))
    (write-prebin-bh-data directory nil result info integrated)))


(mirror-images

 (defun advance-diamond-BH-left (BH-diamond)
   (multiple-value-bind (tag count inert) (handle-bh-tags (diamond-e BH-diamond) BH-diamond)
   (let* ((next-diamond (make-diamond))
	  (east (diamond-e BH-diamond))
	  (west (diamond-w BH-diamond))
	  (bh (diamond-bh BH-diamond)))
      
     (when (and (eq east :BH) (eq west :BH))
       (return-from advance-diamond-BH-left nil))

     (when (and (eq east :BHeatit) (eq west :BHeatit))
       (return-from advance-diamond-BH-left nil))
     
     (when (and (or (eq east :BHeatit) (eq west :BHeatit))
		(or (eq east :BH) (eq west :BH))
		(or (eq east :BHdeleted) (eq west :BHdeleted)))
       (return-from advance-diamond-BH-left nil))
     
     (if (or (eq east :BHdeleted) (eq west :BHdeleted))
	 (return-from advance-diamond-BH-left nil))
     (when (null bh)
       (error "Null bh in evolution"))

     
      
     (setf (diamond-nw next-diamond) :BH)
     (setf (diamond-start next-diamond) (diamond-end BH-diamond))
     (setf (diamond-right next-diamond) (diamond-end (diamond-ne BH-diamond)))

     (setf (diamond-bh next-diamond) bh)

     
     (when (and
	    (null *pointbh*)
	    (or (eq *era* :radiation) (eq *era* :radiation-smooth) (eq *era* :matter) (eq *era* :power))
	    (check-point-inbh (diamond-right next-diamond)))
       (setf (diamond-bh (diamond-ne BH-diamond)) bh)
       (setf (diamond-nw (diamond-ne BH-diamond)) :BHeatit)
       (return-from advance-diamond-BH-left nil))

     
     (when (null *pointbh*)
       (when  (< (3vector-distance (diamond-right next-diamond) (diamond-start next-diamond)) 1E-10)
	 (setf (diamond-bh (diamond-ne BH-diamond)) bh)
	 (setf (diamond-nw (diamond-ne BH-diamond)) :BHeatit)
	 (return-from advance-diamond-BH-left nil)))
     
      (setf (diamond-left next-diamond) (3to4vector (4to3vector (3vector- (diamond-start next-diamond) (3vector- (diamond-right next-diamond) (diamond-start next-diamond)))) (4vector-t (diamond-right next-diamond))))
     
   
      (if (or (eq *era* :radiation)
	      (eq *era* :radiation-smooth)
	      (eq *era* :matter)
	      (eq *era* :power))
	 (progn
	   (setf (diamond-end next-diamond) (compute-bh-diamond-end next-diamond (diamond-right next-diamond)))
	   (setf (diamond-left next-diamond) (4vector- (4vector+ (diamond-start next-diamond) (diamond-end next-diamond)) (diamond-right next-diamond))))
       (progn
	 (setf (diamond-end next-diamond) (compute-bh-diamond-end next-diamond (diamond-right next-diamond)))))

     (setf (diamond-se next-diamond) (diamond-ne BH-diamond))
     (setf (diamond-nw (diamond-se next-diamond)) next-diamond)

         
      ;detach
     (setf (diamond-sw (diamond-ne BH-diamond)) nil)
     (setf (diamond-ne BH-diamond) nil)


     (setf (diamond-countup next-diamond) count)
     (setf (diamond-tag next-diamond) tag)
     (setf (diamond-inertp next-diamond) inert)

     (when (and (eq count 2) inert *pointbh*)
       (format t "~%Removing bh-loop~%")
       (record-bh-loop BH-diamond east)
       (incf *bh-loop-counter*)
       (setf (diamond-nw next-diamond) :deleted))
         
     (handle-new-diamond next-diamond :predecessor BH-diamond)
     next-diamond
     )
   ))
)


(defun bh-radius (global-t)
  (let ((r *bh-size*))
    (when (or (eq *era* :radiation)
	      (eq *era* :radiation-smooth))
      (setf r (/ (* *bh-size* *bh-start*) global-t))) ;physical shrinking
      ;(setf r (/ (* *bh-size* (* (*  *bh-start* *bh-start*) *bh-start*)) (* (* global-t global-t) global-t)))) ; 1/t^3 shrinking
    r))

(defun compute-future-time-old (de teg) ; 1/t^3 shrinking  
  (let ((dr (bh-radius teg))
	(tl)
	(a)
	(b)
	(c)
	(solutions))
    (setf a 1.0)
    (if (> de dr)
	(setf b (* -1.0 (+ de teg)))
      (setf b (- de teg)))
    (if (> de dr)
	(setf c (* *bh-size* (* (*  *bh-start* *bh-start*) *bh-start*)))
      (setf c (* -1.0 (* *bh-size* (* (*  *bh-start* *bh-start*) *bh-start*)))))
    (setf tl (solve-quartic a b 0.0 0.0 c))
    (loop for sl in tl
	  do (when (and (realp sl) (> sl teg))
	       (push sl solutions)))
    (if (or (< (length solutions) 1.0) (> (length solutions) 1.0))
	(error "Not correct tg")
      (car solutions))))
        
       
    
(defun compute-future-time (de teg) ; physical shrinking
  (let ((dr (bh-radius teg))
	(tg1)
	(tg2))
    (if (> de dr)
	(progn
	  (setq tg1 (/ (+ (+ de teg) (sqrt (- (* (+ de teg) (+ de teg)) (* 4.0d0 (* *bh-size* *bh-start*))))) 2.0d0))
	  (setq tg2 (/ (- (+ de teg) (sqrt (- (* (+ de teg) (+ de teg)) (* 4.0d0 (* *bh-size* *bh-start*))))) 2.0d0)))
      (progn
	(setq tg1 (/ (+ (- teg de) (sqrt (+ (* (- teg de) (- teg de)) (* 4.0d0 (* *bh-size* *bh-start*))))) 2.0d0))
	(setq tg2 (/ (- (- teg de) (sqrt (+ (* (- teg de) (- teg de)) (* 4.0d0 (* *bh-size* *bh-start*))))) 2.0d0))))
    (if (> tg1 teg)
	tg1
      (if (> tg2 teg)
	  tg2
	(error "te1 and te2 smaller than tedge")))))
		 

(defun compute-bh-diamond-end (diamond edge)
  (if (or (eq *era* :radiation)
	  (eq *era* :radiation-smooth)
	  (eq *era* :matter)
	  (eq *era* :power))
      (if *pointbh*
	  (let* ((de (3vector-distance edge (diamond-start diamond)))
		 (te (4vector-t edge))
		 (teg (global-time te)))
	    (3to4vector (4to3vector (diamond-start diamond)) (local-time (+ teg de))))
	(let* ((bh (diamond-bh diamond))
	       (center (standardize-position (localize-position (blackhole-center bh)) edge))
	       (dir (3vector-normalize (3vector- edge center)))
	       (te (4vector-t edge))
	       (teg (global-time te))
	       (de (3vector-distance edge center))
	       (tl))
	  (setf tl (local-time (compute-future-time de teg)))
	  (3to4vector (3vector+ (4to3vector center) (3vector-scale dir (bh-radius (global-time tl)))) tl)))	
    (3to4vector (4to3vector (diamond-start diamond)) (4vector-t (compute-diamond-end-flat diamond)))))

	      

;put initial BH-string in the simulation
(defun set-BH-string-abprime (astart bstart aprime-function bprime-function diamond-number sigma-length)
  (let ((right-BH-diamond (create-BH-string-abprime astart bstart aprime-function bprime-function diamond-number sigma-length))
       )
       ;(initialize :total-size 20)
       (loop for this = right-BH-diamond then (diamond-w this)
          do (add-diamond-tags this)
          do (handle-new-diamond this)
          until (eq (diamond-w this) :BH)
       )
       right-BH-diamond
  )
)




(defun do-movie-jpg (dt nt filename)
       (plot-time-slice 
                        :xrange (list -1.5 1.5) 
                        :yrange (list -1.5 3.0) 
                        :zrange (list -1.0 1.0) 
                       ; :prelude (format nil "set view 60, 60~%") 
                       ; :epsf (concatenate 'string filename "/0.eps")
                        :prelude (progn (format nil "set view 60, 30~%")
                                        (format nil "set output '~A/0.jpg'~%set terminal jpeg~%" filename)
                                 )
                        :exit t         
       )
       (loop for i from 1 to nt by 1
          do (evolve-until (* i dt))
          do (plot-time-slice 
                              :xrange (list (+ -1.5 (/ (* i 1.5) (+ nt 30))) 
                                            (- 1.5  (/ (* i 1.5) (+ nt 30)))
                                      ) 
                              :yrange (list (+ -1.5 (/ (* i 1.5) (+ nt 30)))
                                            (- 3.0  (/ (* i 3.0) (+ nt 30)))
                                      ) 
                              :zrange (list (+ -1.0 (/ (* i 1.0) (+ nt 30)))
                                            (- 1.0  (/ (* i 1.0) (+ nt 30)))
                                      ) 
                             ; :prelude (format nil "set view 60, 60~%")
                             ; :epsf (concatenate 'string filename "/" (write-to-string i) ".eps")
                              :prelude (progn (format nil "set view 60, 30~%")
                                              (format nil "set output '~A/~D.jpg'~%set terminal jpeg~%" filename i)
					      )
                              :exit t
             )
       )
)


(defun do-movie-eps (dt nt filename)
       (plot-time-slice 
                        :xrange (list -1.5 1.0) 
                        :yrange (list -1.5 3.0) 
                        :zrange (list -1.0 1.0) 
                        :prelude (format nil "set view 60, 60~%") 
                        :epsf (concatenate 'string filename "/0.eps")
                        ;:prelude (progn (format nil "set view 60, 60~%")
                        ;                (format nil "set output '~A/0.jpg'~%set terminal jpeg~%" filename)
                        ;         )
                        ;:exit t         
       )
       (loop for i from 1 to nt by 1
          do (evolve-until (* i dt))
          do (plot-time-slice 
                              :xrange (list (+ -1.5 (/ (* i 1.5) (+ nt 1))) 
                                            (- 1.0  (/ (* i 1.0) (+ nt 1)))
                                      ) 
                              :yrange (list (+ -1.5 (/ (* i 1.5) (+ nt 1)))
                                            (- 3.0  (/ (* i 3.0) (+ nt 1)))
                                      ) 
                              :zrange (list (+ -1.0 (/ (* i 1.0) (+ nt 1)))
                                            (- 1.0  (/ (* i 1.0) (+ nt 1)))
                                      ) 
                              :prelude (format nil "set view 60, 60~%")
                              :epsf (concatenate 'string filename "/" (write-to-string i) ".eps")
                             ; :prelude (progn (format nil "set view 60, 60~%")
                             ;                 (format nil "set output '~A/~D.jpg'~%set terminal jpeg~%" filename i)
                             ;          )
                             ; :exit t
             )
       )
)


;;Functions written to create a network of strings attached to Black Holes. The idea of this functions is to create
; firstly a string nwtwork which is in scaling using the Vachaspati Vilenkin initial conditions and then create
; the blackholes in random positions. The process to create such a network is as follows:
; 1) Create a scaling network initializing the system with VV initial conditions
; 2) Create a blackholes.dat file with blackhole positions and the sizes. The positions are chosen randomly 
; 3) At the desire Blackhole creation time intersections between strings and blackhole surfaces are detected
; 4) The strings are split at the intersection points making them start/end into a blackhole diamond
; 5) Once all intersections are performed every string in the netwrok is analysed to check if it is inside a BH or not
; 6) Strings which lie inside a BH are deleted
; 7) The network evolve until the end time without checking more BH intersections 

(defun bh-schedule (start interval end)
  (let ((times nil)
	(length (+ (/ (- end start) interval) 1)))
    (loop for i from 1 to length
	  do (push (+ start (* (- i 1) interval)) times ))
    (setf times (reverse times))
    times))
    


;This function creates black-holes in random position with a given radius and store them in 
;blackholes.dat 
;(defun sort-bhs ( )
(defun sort-bhs (directory size bh-size bh-number bh-start)
  (format t "Creating ~S blackholes with size ~S in a ~S box~%" bh-number bh-size size)
  (cond (size                           ;Size given: normal case.  It is the periodicity distance.
         (initialize :total-size size :split-factor 1
                     :bh-size bh-size :bh-number bh-number :bh-start bh-start))
        (t                              ;Infinite volume
         (initialize :total-size nil)))
  (format t "File created in ~S~%" (concatenate 'string (namestring directory) "/blackholes.dat"))
  (format t "Furthest point at : ~S~%" (round (3vector-length (itox (make-4vector 1.0 1.0 0.0 0.0)))))
  (let* ((stream (open (concatenate 'string (namestring directory) "/blackholes.dat") :direction :output :if-does-not-exist :create :if-exists :supersede))
 	 (center (make-zero-3vector))
	 (vec (make-zero-3vector))
	 (rdn 1.0d0)
	 (rdnN nil)
	 (Per-vec (periodicity-vectors))
	 (per-vec-w (periodicity-vectors))
	 (size bh-size)
	 (bh (make-blackhole :center center :size size))
	 (bh-list (make-hash-table :test 'equalp))
	 (safe t)
	 (tcom (get-universal-time))
	 (prev 0)
	 (bhn 0)
	 )
    (loop
     do (setf safe t)
     do (setf per-vec-w per-vec)
	  do (loop for i below 12
		   do (setf vec (car per-vec-w))
		   do (setf per-vec-w (cdr per-vec-w))
		   do (setf rdnN (random rdn))
		   do (setf center (3vector+ (3vector-scale vec rdnN) center)))
	  do (setf center (standardize-position center))
	  do (loop for k being each hash-key of bh-list
		   do (when (< (3vector-distance-wrap center (blackhole-center (gethash k bh-list))) (* size 2.5))
			(setf safe nil))
		   until (< (3vector-distance-wrap center (blackhole-center (gethash k bh-list))) (* size 2.5)))
	  do (when safe
	       (setf bh (make-blackhole :center center :size size))
	       (setf (gethash center bh-list) bh)
	       (setf bhn (+ bhn 1)))
	  do (setf center (make-zero-3vector))
	  do (when (and (zerop (rem bhn 10000))
			(< prev bhn))
	       (format t "Number of bhs created: ~S~%" bhn)
	       (format t "Last 10000 created in ~S s~%" (- (get-universal-time) tcom))
	       (setq tcom (get-universal-time))
	       (setf prev bhn))
	  until (eq bhn bh-number)
	  )
    (loop for k being each hash-key of bh-list
	  do (when (gethash k bh-list)
		     (write (gethash k bh-list) :readably t :stream stream)
		     (terpri stream)))
    (close stream)
    ))



(defun my-bhs ()
  (format t "~%Starting time of job number ~S size ~S and toff ~S is ~S~%" *job-number* *total-size* *time-offset* (global-job-start *total-size* *time-offset* *job-number*)) 
  (format t "BH-star is: ~S~%" *bh-start*)
  (clrhash *job-bhs*)
  (when (or *all-bhs*
	    (and (> (global-job-end *total-size* *time-offset* *job-number*) *bh-start*)
		 (< (global-job-start *total-size* *time-offset* *job-number*) *bh-start*)))
   (with-open-file (str-in "blackholes.dat" :if-does-not-exist nil)
   ;(with-open-file (str-in (concatenate 'string *input-directory* "/blackholes.dat") :if-does-not-exist nil) 
		    (let* ((bh nil)
			   (center nil)
			   (dis (3vector-length (itox (make-4vector 1.0 1.0 0.0 0.0)))))
		      (loop
		       do (setf bh (read str-in nil))
		       until (null bh)
		       do (setf center (standardize-position (localize-position (blackhole-center bh))))
		       do (when (and (< (3vector-length center) (* 1.5 dis))
				     (null (gethash (blackhole-center bh) *job-bhs*)))
			    (setf (gethash (blackhole-center bh) *job-bhs*) bh)))))
    (format t "Bhs read from file N=~S~%" (hash-table-count *job-bhs*))))

(defun clean-mybhs ()
  (let* ((bh nil)
	 (center nil)
	 (dis (3vector-length (itox (make-4vector 1.0 1.0 0.0 0.0)))))
    (loop for k being each hash-key of *my-bhs*
          do (setf bh (gethash k *my-bhs*))
	  do (setf center (standardize-position (localize-position (blackhole-center bh))))
	  do (when (> (3vector-length center) (* 1.5 dis))
	       (remhash k *my-bhs*))))
  (format t "My bhs are ~S~%" (hash-table-count *my-bhs*)))
    

(defun get-mybhs ()
  (if (or (and  *all-bhs*
		(< (global-job-start *total-size* *time-offset* *job-number*) *bh-start*)
		(> (global-job-end *total-size* *time-offset* *job-number*) *bh-start*))
	  (and (null *all-bhs*)
	       (< (global-job-start *total-size* *time-offset* *job-number*) *bh-start*)
	       (> (global-job-end *total-size* *time-offset* *job-number*) *bh-start*)
	       (< (abs (- (local-time *bh-start*) (current-time))) 1E-5)))
      *job-bhs*
    *my-bhs*))

(defun create-point-bhs ()
  (let ((starts nil))
    (map-string-paths
     #'(lambda (st &rest args)
	 (push (list st args) starts)))
    (loop for it in starts
	  do (apply #'create-point it)))
  )

(defun create-point (start &rest args)
  (let ((current-diamond start)
        (start-junction nil)
        (end nil)
        (end-junction nil)
        (intersections nil)
        (int nil)
	;(prev nil)
	(create t)
	(min-dis 5.0))
    (setq args (car args))
    (setf
     start-junction (car args)
     end (car (cdr args))
     end-junction (car (cdr (cdr args))))
    (if (null end)
        (setf end start))
    (if (or (eq start-junction :deleted)
            (eq end-junction :deleted))
        (return-from create-point nil))
    (if (or (diamond-inertp start)
            (diamond-inertp end))
        (return-from create-point nil))
    (loop ;loop over the diamonds in the string
     do (setf create t)
     do (when (diamond-inertp current-diamond)
	  (setf create nil))
     do (setf int (check-bhpoint current-diamond))
     do (when (and int intersections)
	  (loop for i in intersections
		do (when (< (3vector-distance (bh-intersec-spacetime int) (bh-intersec-spacetime i)) min-dis)
		     (setf create nil))
		until (null create)))
     do (when (and int *my-bhs*)
	  (loop for k being each hash-key of *my-bhs*
		do (when (< (3vector-distance (standardize-position (localize-position (blackhole-center (gethash k *my-bhs*)))) (standardize-position (localize-position (blackhole-center (bh-intersec-bh int))))) min-dis)
		     (setf create nil))
		until (null create)))
     do (when (and int *job-bhs*)
	  (loop for k being each hash-key of *job-bhs*
		do (when (< (3vector-distance (standardize-position (localize-position (blackhole-center (gethash k *job-bhs*)))) (standardize-position (localize-position (blackhole-center (bh-intersec-bh int))))) min-dis)
		     (setf create nil))
		 until (null create)))
     do (when (and int
		   ;(null prev)
		   create
		   (null (diamond-finalp current-diamond))
                   (null (diamond-inertp current-diamond))) ; Inert diamonds are not of interest
	  (push int intersections)
	  ;(setf prev t)
	  )
     ;do (when (null int)
	  ;(setf prev nil))
     do (if (and (eq start end)
		 (null (eq (diamond-e current-diamond) :BH))
		 (diamond-e current-diamond))
	    (setf current-diamond (diamond-e current-diamond)))
     until (eq end current-diamond)
     do (if (null (eq start end))
            (setf current-diamond (diamond-e current-diamond))))
    (when intersections ;if any intersection is detected
      (loop for i in intersections
	    do (perform-bh-point i)))
    ))

(defun check-bhpoint (d)
  (when (and (< (4vector-t (diamond-start d)) (current-time))
	     (> (4vector-t (diamond-end d)) (current-time)))
    (when (eq (random (round (/ 1.0 *bh-probability*))) 0)
      (let* ((return-int (make-bh-intersec))
	     (abs nil)
	     (a nil)
	     (b nil)
	     (center (make-zero-3vector))
	     (bh (make-blackhole :center center :size 0.0))
	     (intersection-point nil)
	     (time (local-time *bh-start*)))
      
	(if (> (4vector-t (diamond-position d :a 0.5 :b 0.5)) time)                                                                                                        
	    (setf abs (multiple-value-list (interpolate-find-coordinate d 3 time :a1 0.0 :b1 0.0 :a2 0.5 :b2 0.5 :error t)))
	  (setf abs (multiple-value-list (interpolate-find-coordinate d 3 time :a1 0.5 :b1 0.5 :a2 1.0 :b2 1.0 :error t))))


	
	(format t "abs: ~S~%" abs)
	
	(setf a (car  abs))
	(setf b (car (cdr abs)))

	(setf intersection-point (diamond-position d :a a :b b))

	(when (null (point-mine intersection-point))
	  (format t "Not my point~%")
	  (return-from check-bhpoint nil))
	
	(setf (blackhole-center bh) (globalize-position intersection-point))
      
	(setf (bh-intersec-diamond return-int) d)
	(setf (bh-intersec-bh return-int) bh)
	(setf (bh-intersec-a return-int) a)
	(setf (bh-intersec-b return-int) b)
	(setf (bh-intersec-spacetime return-int) intersection-point)

	(format t "At global time: ~S~%" (global-time (4vector-t intersection-point)))
	(format t "Created point bh at: ~S~%" return-int)
	
	return-int
	)
      )
    )
  )


(defun perform-bh-point (int)
   (format t "#########################~%")
  (format t "Performing int~%")
  (let* ((d (bh-intersec-diamond int))
         (bh (bh-intersec-bh int))
	 (a (bh-intersec-a int))
         (b (bh-intersec-b int))
         (intersection-point (bh-intersec-spacetime int))
         (global-point (globalize-position intersection-point))
	 (e (east-quarter d intersection-point global-point :a a :b b))
	 (w (west-quarter d intersection-point global-point :a a :b b)))
	    
 
    ;(unless (gethash (blackhole-center bh) *my-bhs*)
    (format t "Pushing to hash~%")
    (setf (gethash (blackhole-center bh) *my-bhs*) bh)


    (divide-pending-intersections d :a a :b b :east e :west w)
    (rescale-right-junctions d e :a a :b b)
    (rescale-left-junctions d w :a a :b b)
    (discard-object d)
    (propagate-cut-ne e a)
    (propagate-cut-nw w b)
    (handle-new-diamond e :predecessor :meiosis)
    (handle-new-diamond w :predecessor :meiosis)
    (advance-cut-point-BH e  w global-point bh nil)
    ))    


(defun writepointbhs ()
  (let* ((stream (open (concatenate 'string (write-to-string *job-number*)  "-blackholes.dat") :direction :output :if-does-not-exist :create :if-exists :supersede))
	 (bh nil))
    (loop for k being each hash-key of *my-bhs*
          do (setf bh (gethash k *my-bhs*))
	  do (write bh :readably t :stream stream)
	  do (terpri stream))
    (close stream)))


	  
;Function to be called ones the string network is created and in scaling. This function is called at BH creation time 
;It maps all the strings in the network and stores the starting diamond of all of them in starts. Then, each one of
;the starting diamonds Is passed to detec-intersec, the function that creates the BH intersections. Once all of them
;are created the new network is ones again mapped and the new starting diamonds are stored once again in starts. 
;finally the starting points are passed to remove-strings, a function that decides if the string is inside a BH and
;deletes them.
(defun create-bhs ()
  (format t "~%Analysing BH-Ints~%")
  (let ((starts nil))
 ;   (format t "Mapping")
    (map-string-paths
     #'(lambda (st &rest args)
	 (push (list st args) starts)))
    (loop for it in starts
	  do (apply #'detec-intersec it)))
  (format t "Deleting Strings~%")
  (let ((starts nil))
    (map-string-paths
     #'(lambda (st &rest args)
	 (push (list st args) starts)))
    (loop for it in starts
	  do (apply #'remove-strings it)))
  (format t "Strings deleted~%"))
 


;This function receives a starting diamond and analyses every single diamond of each string 
;A dimond is checked by check-bh-intersec function to see if it intersecs with a BH and if this is 
;the case the diamond is stored in intersections. Once the string is analysed the diamonds with intersections
;are passed to perform-bh-intersec to ccreate the BH intersection.
(defun detec-intersec (start &rest args)
  (let ((current-diamond start)
	(start-junction nil)
	(end nil)
	(end-junction nil)
	(intersections nil)
	(int nil))
    (setq args (car args)) 
    (setf 
     start-junction (car args)
     end (car (cdr args))
     end-junction (car (cdr (cdr args))))

    (if (null end) ;if the end diamond is null it means that the string is a loop end=start
	(setf end start))
    
    ;(if (and (eq end start)   ; prevent the analysis of strings with only one diamond
	;     (or (null (diamond-e start))
	;	 (eq (diamond-e start) :deleted)
	;	 (eq (diamond-e start) :BH)
	;	 (eq (diamond-e start) :BHdeleted)
	;	 (eq (diamond-e start) :BHpropdel)))
	;(return-from detec-intersec nil))
    (if (or (eq start-junction :deleted) ; the string is deleting so we are not going to analyse it
	    (eq end-junction :deleted))
	(return-from detec-intersec nil))
    (if (or (diamond-inertp start) ;Inert diamonds are not analyse 
	    (diamond-inertp end))
	(return-from detec-intersec nil))
      (loop ;loop over the diamonds in the string
;     do (format t "Current-diamond ~S~%" current-diamond)
;       do (when (diamond-inertp current-diamond)
					;	    (setf inert t))
       do (setf int (check-bh-intersec current-diamond))
       do (when (and int
;		   (possible-int current-diamond)
;		   (null (diamond-finalp current-diamond))
		   (null (diamond-inertp current-diamond))) ; Inert diamonds are not of interest
	    (push int intersections))
       do (if (and (eq start end) (diamond-e current-diamond)
		   (null (eq (diamond-e current-diamond) :BH))
		   (null (eq (diamond-e current-diamond) :BHdeleted))
		   (null (eq (diamond-e current-diamond) :BHpropdel))
		   (null (eq (diamond-e current-diamond) :BHeatit)))
	    (setf current-diamond (diamond-e current-diamond)))
     until (eq end current-diamond)
     do (if (null (eq start end))
	    (setf current-diamond (diamond-e current-diamond))))
      (when intersections ;if any intersection is detected
	;(setf intersections (nreverse intersections))
	(loop for intlist in (nreverse intersections) ;loop around them and perform the intersection
	      do (loop for int in intlist
		       do (when int
			    (loop for i in int
				  do (perform-bh-intersec i))))))
    ))

(defun out-strings (start &rest args)
  (let ((start-junction nil)
        (end nil)
        (end-junction nil))
    (setq args (car args))
    (setf
     start-junction (car args)
     end (car (cdr args))
     end-junction (car (cdr (cdr args))))
    (if (null end)
	(setf end start))
    (when (or (eq start-junction :BH) (eq end-junction :BH))
      (format t "String: st ~S, st-J ~S, end ~S, end-J ~S ~%" start start-junction end end-junction))))


(defun remove-bh-loop (start end)
  (cond ((and (eq (diamond-w start) :BH) (eq (diamond-e end) :BH) (< (3vector-distance (diamond-end start) (diamond-end end)) 1E-8))
         t)
        (t
         nil)))

;Function that analyses if the string starting at start is inside a BH or not.
;If the string is inside a BH it is deleted. To decide if a string is inside a BH
;each diamond in the string is analysed to see if it lies inside the BH. Then, 
;the number of diamonds inside/outside BH are counted is the ratio between
;diamonds out against diamonds in is smaller than 0.8 the string is removed.
(defun remove-strings (start &rest args)
  (let ((current-diamond start)
	(end nil)
	(del nil)
	(d-e nil)
	(d-w nil)
	(d-nenw nil)
	(d-nesw nil)
	(d-senw nil)
	(d-sesw nil)
	(d-ne nil)
	(d-se nil)
	(d-nw nil)
	(d-sw nil))
    (setq args (car args))
    (setf
     end (car (cdr args)))
    (if (null end) ; if the end diamond is null it means that the string is a loop start=end
	(setf end start))
    (if (and (eq end start) ;prevent the analysis of strings with only one diamond
	     (or (null (diamond-e start))
		 (eq (diamond-e start) :deleted)
		 (eq (diamond-e start) :BH)
		 (eq (diamond-e start) :BHdeleted)
		 (eq (diamond-e start) :BHpropdel)
		 (eq (diamond-e start) :BHeatit)))
	(return-from remove-strings nil))
    (if (or (eq (diamond-w start) :deleted) ;this string is already deleting do not do anything
	    (eq (diamond-e end) :deleted))
	(return-from remove-strings nil))
    (if (or (diamond-inertp start) ;Inert diamonds are not analyse 
            (diamond-inertp end))
        (return-from remove-strings nil))
    (setf current-diamond start)
    (if t ;(remove-this-string start end)
    	(progn
	  (loop
	   do (if (or (diamond-finalp current-diamond) ;final diamonds are not discarded, therefore neighbors are analysed and one of them is discarded it is replaced with BHdeleted
		      (>= (4vector-t (diamond-start current-diamond)) (+ (current-time) 1E-10))
		      (< (4vector-t (diamond-end current-diamond)) (- (current-time) 1E-10))
		      (null (check-inbh current-diamond)) 
		      (eq (diamond-e current-diamond) :BH)
		      (eq (diamond-w current-diamond) :BH)
		      )
		  (let ((ne (diamond-ne current-diamond))
			(nw (diamond-nw current-diamond))
			(se (diamond-se current-diamond))
			(sw (diamond-sw current-diamond)))
		    (if (and ne se)
			(push current-diamond d-e)
		      (if (and nw sw)
			  (push current-diamond d-w) 
			(if (and ne nw)
			    (push current-diamond d-nenw)
			  (if (and ne sw)
			      (push current-diamond d-nesw)
			    (if (and se nw)
			      (push current-diamond d-senw)
			      (if (and se sw)
				  (push current-diamond d-sesw)
				(if ne
				    (push current-diamond d-ne)
				  (if se
				      (push current-diamond d-se)
				    (if nw
					(push current-diamond d-nw)
				      (if sw
					(push current-diamond d-sw))))))))))))
		  (push current-diamond del))
	     do (if (and (eq start end) (diamond-e current-diamond))
		    (setf current-diamond (diamond-e current-diamond)))
	     until (eq end current-diamond)
	     do (if (null (eq start end))
		    (setf current-diamond (diamond-e current-diamond))))))
    
    (if del ;remove the diamonds
	(loop for di in del
	      do (if (diamond-start di)
		     (progn
		       (setf (left-rejoining-junction di) nil)
		       (setf (right-rejoining-junction di) nil)
		       (discard-object di)))))
    ;neighbor analysis to replace them correctly with BHdeleted
    (if d-e
	(progn
	  (loop for di in d-e
		do (if (null (diamond-se di))
		       (setf (diamond-se di) :BHpropdel))
		do (if (null (diamond-ne di))
		       (setf (diamond-ne di) :BHpropdel)))))
	
    (if d-w
	(progn
	  (loop for di in d-w
		do (if (null (diamond-sw di))
		       (setf (diamond-sw di) :BHpropdel))
		do (if (null (diamond-nw di))
		       (setf (diamond-nw di) :BHpropdel)))))

    (if d-nenw
	(progn
	  (loop for di in d-nenw
		do (if (null (diamond-ne di))
		       (setf (diamond-ne di) :BHpropdel))
		do (if (null (diamond-nw di))
		       (setf (diamond-nw di) :BHpropdel)))))

	  

    (if d-nesw
	(progn
	  (loop for di in d-nesw
		do (if (null (diamond-ne di))
		       (setf (diamond-ne di) :BHpropdel))
		do (when (null (diamond-sw di))
		     (setf (diamond-nw di) :BHpropdel)))))

    (if d-senw
	(progn
	  (loop for di in d-senw
		do (when (null (diamond-se di))
		     (setf (diamond-ne di) :BHpropdel))
		do (if (null (diamond-nw di))
		       (setf (diamond-nw di) :BHpropdel)))))

    (if d-sesw
	(progn
	  (loop for di in d-sesw
		do (when (null (diamond-se di))
		     (setf (diamond-ne di) :BHpropdel))
		do (when (null (diamond-sw di))
		     (setf (diamond-nw di) :BHpropdel)))))


    (if d-ne
	(progn
	  (loop for di in d-ne
		do (if (null (diamond-ne di))
		       (setf (diamond-ne di) :BHpropdel)))))

    (if d-se
	(progn
	  (loop for di in d-se
		do (when (null (diamond-se di))
		     (setf (diamond-ne di) :BHpropdel)))))

    (if d-nw
	(progn
	  (loop for di in d-nw
		do (if (null (diamond-nw di))
		       (setf (diamond-nw di) :BHpropdel)))))

    (if d-sw
	(progn
	  (loop for di in d-sw
		do (when (null (diamond-sw di))
		     (setf (diamond-nw di) :BHpropdel)))))

    ))


(defun remove-this-string (start end)
  (cond ((and (eq (diamond-e end) :BHdeleted) (eq (diamond-w start) :BHdeleted)
	      (<= (4vector-t (diamond-start end)) (+ (current-time) 1E-10))
	      (<= (4vector-t (diamond-start start)) (+ (current-time) 1E-10))
	      (> (4vector-t (diamond-end end)) (- (current-time) 1E-10))
	      (> (4vector-t (diamond-end start)) (- (current-time) 1E-10)))
	 t)
	((or (and (eq (diamond-e end) :BHdeleted)
		  (<= (4vector-t (diamond-start end)) (+ (current-time) 1E-10))
		  (> (4vector-t (diamond-end end)) (- (current-time) 1E-10)))
	     (and (eq (diamond-w start) :BHdeleted)
		  (<= (4vector-t (diamond-start start)) (+ (current-time) 1E-10))
		  (> (4vector-t (diamond-end start)) (- (current-time) 1E-10))))
	 t)
	((or (and (eq (diamond-w start) :BHdeleted) (null (eq (diamond-e end) :BHdeleted)) (<= (4vector-t (diamond-start start)) (+ (current-time) 1E-10))
		  (> (4vector-t (diamond-end start)) (- (current-time) 1E-10)))
	     (and (eq (diamond-e end) :BHdeleted) (null (eq (diamond-w start) :BHdeleted)) (<= (4vector-t (diamond-start end)) (+ (current-time) 1E-10))
		  (> (4vector-t (diamond-end end)) (- (current-time) 1E-10))))
	 t)
	((and (eq (diamond-w start) :BH) (eq (diamond-e end) :BH) (< (3vector-distance (diamond-end start) (diamond-end end)) 1E-8))
	 t)
	(t
	 nil)))
; Function that checks if a specific diamond is inside a BH. This function reads the
; positions of the BH from blacholes.dat and diamond and bh info is tramsmitted to check-diamond-inbh function
(defun check-inbh (diamond)
  (if (eq diamond :BH)
      t
    (if (null diamond)
	nil
      (if (eq diamond :deleted)
	  nil
	(let* ((mybhs (get-mybhs))
	       (center (make-3vector 0.0d0 0.0d0 0.0d0))
	       (size 0.d0)
	       (bh (make-blackhole :center center :size size))
	       (inbh nil)
	       (sta-time (4vector-t (diamond-start diamond)))
	       (end-time (4vector-t (diamond-end diamond))))
	  (if (or (>= (- (current-time) 1E-10) end-time)
		  (< (+ (current-time) 1E-10) sta-time))
	      (return-from check-inbh nil))
	  (loop for k being each hash-key of mybhs
	   do (setf bh (gethash k mybhs))
	   do (if (check-diamond-inbh diamond bh)
		  (setf inbh t))
	   )
	  (if inbh
	      t
	    nil))))))



(defun check-inbh-out (diamond)
  ;(format t "Analysing: ~S~%" diamond)
;  (format t "E: ~S~%" (diamond-e diamond))
 ; (format t "W: ~S~%" (diamond-w diamond))
 ; (format t "Current time: ~S~%" (current-time))
 ; (format t "St-t: ~S~%" (4vector-t (diamond-start diamond)))
 ; (format t "En-t: ~S~%" (4vector-t (diamond-end diamond)))
 ; (format t "Int: ~S~%" (check-bh-intersec diamond))
  (if (eq diamond :BH)
      t
    (if (null diamond)
        nil
      (if (eq diamond :deleted)
          nil
	(let* ((mybhs (get-mybhs))
	       (center (make-3vector 0.0d0 0.0d0 0.0d0))
	       (size 0.d0)
	       (bh (make-blackhole :center center :size size))
	       (inbh nil)
	       (sta-time (4vector-t (diamond-start diamond)))
	       (end-time (4vector-t (diamond-end diamond))))
	  (if (or (>= (- (current-time) 1E-10) end-time)
		  (< (+ (current-time) 1E-10) sta-time))
	      (progn
		(format t "Time out of range~%")
		(return-from check-inbh-out nil)))
	  (loop for k being each hash-key of mybhs
		do (setf bh (gethash k mybhs))
		do (if (check-diamond-inbh-out diamond bh)
		       (setf inbh t))
		)
	  (if inbh
	      t
	    nil))))))

(defun check-point-inbhhalo (p)
  (let* ((mybhs (get-mybhs))
	 (center (make-zero-3vector))
	 (size 0.0d0)
	 (bh (make-blackhole :center center :size size))
	 (inbh nil)
	 (short 1000.0d0)
	 (bhn 0)
	 )
    (loop for k being each hash-key of mybhs
	  do (setf bh (gethash k mybhs))
	  do (setf center (standardize-position (localize-position (blackhole-center bh)) p))
	  do (setf size (blackhole-size bh))
	  do (setf bhn (+ bhn 1))
	  do (when (< (3vector-distance p center) short)
	       (setf short (3vector-distance p center)))
	  do (if (<= (3vector-distance p center) (* size 1.25))
		 (setf inbh t)))
    (format t "Number of bhs in the region: ~S~%" bhn)
    (format t "Shortest distance from ~S to a BH center is ~S~%" p short)
    (if inbh
	t
      nil)))

(defun check-point-inbh (p)
  (let* ((mybhs (get-mybhs))
	 (center (make-zero-3vector))
	 (size 0.0d0)
	 (bh (make-blackhole :center center :size size))
	 (inbh nil)
	 (r (bh-radius (global-time (current-time))))
	 (short 1000.0d0)
	 (bhn 0)
	 )
    (loop for k being each hash-key of mybhs
          do (setf bh (gethash k mybhs))
	  do (setf center (standardize-position (localize-position (blackhole-center bh)) p))
	  do (setf size (blackhole-size bh))
	  do (setf bhn (+ bhn 1))
	  do (when (< (3vector-distance p center) short)
	       (setf short (3vector-distance p center)))
	  do (if (<= (3vector-distance p center) r)
		 (setf inbh t)))
    (if inbh
	t
      nil)))

(defun possible-int (d bh)
  (let* ((ti (current-time))
	 (ar (car (multiple-value-list (find-right-edge-position d ti nil))))
	 (br (car (cdr (multiple-value-list (find-right-edge-position d ti nil)))))
	 (al (car (multiple-value-list (find-left-edge-position d ti nil))))
	 (bl (car (cdr (multiple-value-list (find-left-edge-position d ti nil)))))
	 (start (diamond-start d))
	 (center (standardize-position (localize-position (blackhole-center bh)) start))
	 (r (bh-radius (global-time (current-time))))
	 (dl nil)
	 (dr nil)
	 (inbh nil)
	 (res 1E-10)
         (ledge nil)
         (redge nil))
    (when (or (null ar) (null br) (null al) (null bl))
      (return-from possible-int nil))
    (setf dl (3vector-distance (diamond-position d :a al :b bl) center))
    (setf dr (3vector-distance (diamond-position d :a ar :b br) center))
    (when (and (< (abs (- al 1.0)) res)
               (< bl res))
      (setf ledge t))
    (when (and (< ar res)
               (< (abs (- br 1.0)) res))
      (setf redge t))
    (when (and (null ledge) (null redge)
               (or (and (> dl r) (< dr r))
                   (and (< dl r) (> dr r))))
      (setf inbh t))
    (if inbh
	t
      nil)))
			  
(defun possible-int-out (d bh)
  (let* ((ti (current-time))
         (ar (car (multiple-value-list (find-right-edge-position d ti nil))))
         (br (car (cdr (multiple-value-list (find-right-edge-position d ti nil)))))
         (al (car (multiple-value-list (find-left-edge-position d ti nil))))
         (bl (car (cdr (multiple-value-list (find-left-edge-position d ti nil)))))
         (start (diamond-start d))
         (center (standardize-position (localize-position (blackhole-center bh)) start))
	 (r (bh-radius (global-time (current-time))))
         (dl nil)
         (dr nil)
         (inbh nil)
	 (res 1E-10)
	 (ledge nil)
	 (redge nil))
    (when (or (null ar) (null br) (null al) (null bl))
      (return-from possible-int-out nil))
    (setf dl (3vector-distance (diamond-position d :a al :b bl) center))
    (setf dr (3vector-distance (diamond-position d :a ar :b br) center))
    (when (and (< (abs (- al 1.0)) res)
	       (< bl res))
      (setf ledge t))
    (when (and (< ar res)
	       (< (abs (- br 1.0)) res))
      (setf redge t))
    (when (and (null ledge) (null redge)
	       (or (and (> dl r) (< dr r))
		   (and (< dl r) (> dr r))))
      (format t "PR: ~S~%" (diamond-position d :a ar :b br))
      (format t "PL: ~S~%" (diamond-position d :a al :b bl))
      (format t "ar ~S br ~S~%" ar br)
      (format t "al ~S bl ~S~%" al bl)
      (format t "R: ~S~%" r)
      (format t "DR: ~S DL: ~S~%" dr dl)
      (setf inbh t))
    (if inbh
        t
      nil)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Function that loops around the BH in blackholes.dat and transmit the information about the specific BH and the diamond being analysed
;to check-bh-intersec-diamond
(defun check-bh-intersec (diamond)
  (let* ((mybhs (get-mybhs))
	 (center (make-3vector 0.0d0 0.0d0 0.0d0))
	 (size 0.d0)
	 (bh (make-blackhole :center center :size size))
	 (return-intersection nil)
	 (int nil))
    (loop for k being each hash-key of mybhs
          do (setf bh (gethash k mybhs))
	  do (setf center (standardize-position (localize-position (blackhole-center bh)) (diamond-start diamond)))
	  do (when (< (3vector-distance (diamond-start diamond) center) (* 5.0 (blackhole-size bh)))
	       (setf int (check-bh-intersec-diamond diamond bh nil))) 
	  do (if int
		 (push int return-intersection))
	  do (setf int nil))
    (if return-intersection
	return-intersection)))


(defun check-bh-intersec-diamond (d bh re)
  (let* ((max 0.1)
	 (res 1E-8)
	 (start (diamond-start d))
	 (left (diamond-left d))
	 (right (diamond-right d))
	 (end (diamond-end d))
	 (center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
	 (bhtime (current-time))
	 (reads nil)
	 (solutions nil)
	 (outbs nil)
	 (a1 nil)
	 (b1 nil)
	 (a2 nil)
	 (b2 nil)
	 (ac nil)
	 (int1 (make-bh-intersec))
	 (int2 (make-bh-intersec))
	 (bh-diamond-width 1E-100)
	 (Intp nil)
	 (intf nil)
	 (reint nil)
	 (d1 0.0d0)
	 (d2 0.0d0)
	 (r (bh-radius (global-time (current-time))))
	 )
    (setq reads (apply #'solve-intersection (parameters-quartic d bh re)))
    (loop for rds in reads
	  do (if (realp rds)
	      (push rds solutions)
	    (if (< (abs (imagpart rds)) max)
		(push (realpart rds) solutions))))
    (setq solutions (remove-duplicates solutions))
    (when (possible-int d bh)
      (setf intp t))
    (loop for bc in solutions 
	  do (when (and (realp bc) (<= bc 1.0) (>= bc 0.0))
	       (setq ac (/ (- (+ bhtime (* bc (4vector-t start))) (+ (4vector-t start) (* bc (4vector-t right)))) (- (+ (4vector-t left) (+ (* bc (4vector-t start)) (* bc (4vector-t end)))) (+ (4vector-t start) (+ (* bc (4vector-t left)) (* bc (4vector-t right)))))))
	       (when (and (<= ac 1.0) (>= ac 0.0))
		 (setf intf t)
		 (when re
		   (format t "INT ac: ~S bc: ~S~%" ac bc)))
	       (when (or (and (realp ac) (<= ac 1.0) (>= ac 0.0) (point-mine (diamond-position d :a ac :b bc)) (null (eq (diamond-w d) :BH)) (null (eq (diamond-e d) :BH)) (null (eq (diamond-w d) :BHdeleted)) (null (eq (diamond-e d) :BHdeleted)))
			 (and (eq (diamond-w d) :BH) (realp ac) (< ac 1.0) (>= ac 0.0) (<= bc 1.0) (> bc 0.0) (> (- bc ac) res) (> (3vector-distance (diamond-position d :a 0.5 :b 0.5) (diamond-position d :a ac :b bc)) bh-diamond-width) (point-mine (diamond-position d :a ac :b bc)))
			 (and (eq (diamond-e d) :BH) (realp ac) (<= ac 1.0) (> ac 0.0) (< bc 1.0) (>= bc 0.0) (>  (- ac bc) res) (> (3vector-distance (diamond-position d :a 0.5 :b 0.5) (diamond-position d :a ac :b bc)) bh-diamond-width) (point-mine (diamond-position d :a ac :b bc)))
			 (and (eq (diamond-w d) :BHdeleted) (realp ac) (< ac 1.0) (>= ac 0.0) (<= bc 1.0) (> bc 0.0) (> (- bc ac) res) (> (3vector-distance (diamond-position d :a 0.5 :b 0.5) (diamond-position d :a ac :b bc)) bh-diamond-width) (point-mine (diamond-position d :a ac :b bc)))
                         (and (eq (diamond-e d) :BHdeleted) (realp ac) (<= ac 1.0) (> ac 0.0) (< bc 1.0) (>= bc 0.0) (> (- ac bc) res) (> (3vector-distance (diamond-position d :a 0.5 :b 0.5) (diamond-position d :a ac :b bc)) bh-diamond-width) (point-mine (diamond-position d :a ac :b bc))))
		 (push bc outbs)
		 )))
    (when (and intp (null intf) (null re))
      (possible-int-out d bh)
      (format t "Params: ~S~%" (parameters-quartic d bh nil))
      (format t "Solutions: ~S~%" reads)
      (warn "Rechecking int~%")
      (setf reint (check-bh-intersec-diamond d bh t)))
      

    (when (and intp (null intf) re)
      (warn "Not int possible case~%"))


    (when (and (> (length outbs) 0) (<= (length outbs) 1))
      (setq b1 (car outbs))
      (setq a1 (/ (- (+ bhtime (* b1 (4vector-t start))) (+ (4vector-t start) (* b1 (4vector-t right)))) (- (+ (4vector-t left) (+ (* b1 (4vector-t start)) (* b1 (4vector-t end)))) (+ (4vector-t start) (+ (* b1 (4vector-t left)) (* b1 (4vector-t right))))))))

    (when (and (> (length outbs) 1) (<= (length outbs) 2))
      (if (< (car outbs) (car (cdr outbs)))
	  (progn
	    (setq b1 (car outbs))
	    (setq b2 (car (cdr outbs))))
	(progn
	  (setq b2 (car outbs))
	  (setq b1 (car (cdr outbs)))))
      (setq a1 (/ (- (+ bhtime (* b1 (4vector-t start))) (+ (4vector-t start) (* b1 (4vector-t right)))) (- (+ (4vector-t left) (+ (* b1 (4vector-t start)) (* b1 (4vector-t end)))) (+ (4vector-t start) (+ (* b1 (4vector-t left)) (* b1 (4vector-t right)))))))
      (setq a2 (/ (- (+ bhtime (* b2 (4vector-t start))) (+ (4vector-t start) (* b2 (4vector-t right)))) (- (+ (4vector-t left) (+ (* b2 (4vector-t start)) (* b2 (4vector-t end)))) (+ (4vector-t start) (+ (* b2 (4vector-t left)) (* b2 (4vector-t right))))))))

    (when (> (length outbs) 2)
      (error "More than two solutions for the intersection"))

    (when (and a2
	       (< (abs (- a1 a2)) 1E-8)
	       (< (abs (- b1 b2)) 1E-8))
      (format t "TWO INTS TOO CLOSE~%")
      (return-from check-bh-intersec-diamond nil))
    
        
    (when (and (> (length outbs) 0) (<= (length outbs) 2))
      (setf d1 (3vector-distance (diamond-position d :a a1 :b b1) center))
      (setf (bh-intersec-diamond int1) d)
      (setf (bh-intersec-bh int1) bh)
      (setf (bh-intersec-a int1) a1)
      (setf (bh-intersec-b int1) b1)
      (setf (bh-intersec-spacetime int1) (diamond-position d :a a1 :b b1))
      (when a2
	(setf d2 (3vector-distance (diamond-position d :a a2 :b b2) center))
	(setf (bh-intersec-diamond int2) d)
	(setf (bh-intersec-bh int2) bh)
	(setf (bh-intersec-a int2) a2)
	(setf (bh-intersec-b int2) b2)
	(setf (bh-intersec-spacetime int2) (diamond-position d :a a2 :b b2))
	)
      )

    (if (and intp (null intf) (null re))
	reint
      (when (and (> (length outbs) 0) (<= (length outbs) 2))
	(cond ((and (< (abs (- (/ d1 r) 1.0)) 0.1)
		    (< (abs (- (/ d2 r) 1.0)) 0.1))
	       (list int1 int2))
	      ((and (< (abs (- (/ d1 r) 1.0)) 0.1)
		    (null a2))
	       (list int1)))))
    ))
       


(defun parameters-quartic (d bh re)
  (let* ((st (4vector-t (diamond-start d)))
         (sx (4vector-x (diamond-start d)))
         (sy (4vector-y (diamond-start d)))
         (sz (4vector-z (diamond-start d)))
         (rt (4vector-t (diamond-right d)))
         (rx (4vector-x (diamond-right d)))
         (ry (4vector-y (diamond-right d)))
         (rz (4vector-z (diamond-right d)))
         (lt (4vector-t (diamond-left d)))
         (lx (4vector-x (diamond-left d)))
         (ly (4vector-y (diamond-left d)))
         (lz (4vector-z (diamond-left d)))
         (et (4vector-t (diamond-end d)))
         (ex (4vector-x (diamond-end d)))
         (ey (4vector-y (diamond-end d)))
         (ez (4vector-z (diamond-end d)))
         (center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
         (cx (3vector-x center))
         (cy (3vector-y center))
         (cz (3vector-z center))
	 (r (bh-radius (global-time (current-time))))
         (tbh (current-time))
         (t1 (- st tbh))
         (t2 (- lt st))
         (t3 (- rt st))
         (t4 (+ (- st lt) (- et rt)))
         (x1 (- sx cx))
         (x2 (- lx sx))
         (x3 (- rx sx))
         (x4 (+ (- sx lx) (- ex rx)))
         (y1 (- sy cy))
         (y2 (- ly sy))
         (y3 (- ry sy))
         (y4 (+ (- sy ly) (- ey ry)))
         (z1 (- sz cz))
         (z2 (- lz sz))
         (z3 (- rz sz))
         (z4 (+ (- sz lz) (- ez rz)))
	 (p0 (sum (list
		   (list -1.0 r r t2 t2) (list 1.0 t2 t2 x1 x1) (list 1.0 t1 t1 x2 x2) (list 1.0 t2 t2 y1 y1) (list 1.0 t1 t1 y2 y2) (list 1.0 t2 t2 z1 z1) (list 1.0 t1 t1 z2 z2) (list -2.0 t1 t2 x1 x2) (list -2.0 t1 t2 y1 y2) (list -2.0 t1 t2 z1 z2))))
	 (p1 (sum (list
		   (list -2.0 r r t2 t4) (list -2.0 t1 t4 x1 x2) (list 2.0 t1 t3 x2 x2) (list 2.0 t2 t2 x1 x3) (list 2.0 t1 t1 x2 x4) (list -2.0 t1 t4 y1 y2) (list 2.0 t1 t3 y2 y2) (list 2.0 t1 t1 y2 y4) (list -2.0 t1 t4 z1 z2) (list 2.0 t1 t3 z2 z2) (list 2.0 t2 t2 y1 y3) (list 2.0 t2 t2 z1 z3) (list 2.0 t1 t1 z2 z4) (list 2.0 t2 t4 x1 x1) (list 2.0 t2 t4 y1 y1) (list 2.0 t2 t4 z1 z1) (list -2.0 t2 t3 x1 x2) (list -2.0 t2 t3 y1 y2) (list -2.0 t2 t3 z1 z2) (list -2.0 t2 t1 x2 x3) (list -2.0 t2 t1 x1 x4) (list -2.0 t2 t1 y2 y3) (list -2.0 t2 t1 y1 y4) (list -2.0 t2 t1 z2 z3) (list -2.0 t2 t1 z1 z4))))
	 (p2 (sum (list
		   (list -1.0 r r t4 t4) (list 1.0 t4 t4 x1 x1) (list -2.0 t3 t4 x1 x2) (list 1.0 t3 t3 x2 x2) (list -2.0 t1 t4 x2 x3) (list -2.0 t1 t4 x1 x4) (list 4.0 t1 t3 x2 x4) (list 1.0 t1 t1 x4 x4) (list 1.0 t4 t4 y1 y1) (list -2.0 t3 t4 y1 y2) (list 1.0 t3 t3 y2 y2) (list -2.0 t1 t4 y2 y3) (list -2.0 t1 t4 y1 y4) (list 4.0 t1 t3 y2 y4) (list 1.0 t1 t1 y4 y4) (list 1.0 t4 t4 z1 z1) (list -2.0 t3 t4 z1 z2) (list 1.0 t3 t3 z2 z2) (list -2.0 t1 t4 z2 z3) (list 1.0 t2 t2 x3 x3) (list 1.0 t2 t2 y3 y3) (list 1.0 t2 t2 z3 z3) (list -2.0 t1 t4 z1 z4) (list 4.0 t1 t3 z2 z4) (list 1.0 t1 t1 z4 z4) (list 4.0 t2 t4 x1 x3) (list 4.0 t2 t4 y1 y3) (list 4.0 t2 t4 z1 z3) (list -2.0 t2 t3 x2 x3) (list -2.0 t2 t3 x1 x4) (list -2.0 t2 t3 y2 y3) (list -2.0 t2 t3 y1 y4) (list -2.0 t2 t3 z2 z3) (list -2.0 t2 t3 z1 z4) (list -2.0 t2 t1 x3 x4) (list -2.0 t2 t1 y3 y4) (list -2.0 t2 t1 z3 z4))))
	 (p3 (sum (list
		   (list 2.0 t4 t4 x1 x3) (list -2.0 t3 t4 x2 x3) (list -2.0 t3 t4 x1 x4) (list 2.0 t3 t3 x2 x4) (list -2.0 t1 t4 x3 x4) (list 2.0 t1 t3 x4 x4) (list 2.0 t4 t4 y1 y3) (list -2.0 t3 t4 y2 y3) (list -2.0 t3 t4 y1 y4) (list 2.0 t3 t3 y2 y4) (list -2.0 t1 t4 y3 y4) (list 2.0 t1 t3 y4 y4) (list 2.0 t4 t4 z1 z3) (list -2.0 t3 t4 z2 z3) (list -2.0 t3 t4 z1 z4) (list 2.0 t3 t3 z2 z4) (list -2.0 t1 t4 z3 z4) (list 2.0 t1 t3 z4 z4) (list 2.0 t2 t4 x3 x3) (list 2.0 t2 t4 y3 y3) (list 2.0 t2 t4 z3 z3) (list -2.0 t2 t3 x3 x4) (list -2.0 t2 t3 y3 y4) (list -2.0 t2 t3 z3 z4))))
	 (p4 (sum (list
		   (list 1.0 t4 t4 x3 x3) (list -2.0 t3 t4 x3 x4) (list 1.0 t3 t3 x4 x4) (list 1.0 t4 t4 y3 y3) (list -2.0 t3 t4 y3 y4) (list 1.0 t3 t3 y4 y4) (list 1.0 t4 t4 z3 z3) (list -2.0 t3 t4 z3 z4) (list 1.0 t3 t3 z4 z4))))
	 )
    (if re
	(list 0.0 p3 p2 p1 p0)
      (list p4 p3 p2 p1 p0))
    ))

(defun solve-intersection (a b c d e)
  (let ((max4 1E-15) ;1E-15 ; if not numerical errors create ficticius bh-int
	(max3 1E-12) ;1E-12
	(max21 1E-30))
    (if (and (< (abs a) max4)
	     (< (abs b) max3)
	     (< (abs c) max21)
	     (< (abs d) max21)
	     (< (abs e) max21))
	(return-from solve-intersection nil)
      (if (and (< (abs a) max4)
	       (< (abs b) max3)
	       (< (abs c) max21)
	       (< (abs d) max21))
	  (return-from solve-intersection nil)
	(if (and (< (abs a) max4)
		 (< (abs b) max3)
		 (< (abs c) max21))
	    (solve-linear d e)
	  (if (and (< (abs a) max4)
		   (< (abs b) max3))
	      (solve-cuadratic c d e)
	    (if (< (abs a) max4)
		(solve-cubic b c d e)
	      (solve-quartic a b c d e))))))))


(defun solve-linear (d e)
  (let ((bs1 (* -1.0 (/ e d)))) 
    (list bs1)))

(defun solve-cuadratic (c d e)
  (let* ((bs1 (/ (* -1.0 (+ d (* 1.0 (sqrt (- (* d d) (* 4.0 (* c e))))))) (* 2.0 c)))
	 (bs2 (/ (* -1.0 (+ d (* -1.0 (sqrt (- (* d d) (* 4.0 (* c e))))))) (* 2.0 c))))
    (list bs1 bs2)))

(defun solve-cubic (b c d e)
  (let* ((max 1E-3)
	 (p1 (/ c (* 3.0 b)))
	 (p2 (- (* 3.0 (* b d)) (* c c)))
	 (p31 (- (* 9.0 (* b (* c d))) (+ (* 2.0 (* c (* c c))) (* 27 (* b (* b e))))))
	 (p3 (expt (+ p31 (sqrt (+ (* 4.0 (* p2 (* p2 p2))) (* p31 p31)))) (/ 1.0 3.0)))
	 (bs1 (if (and (zerop (* 3.0 (* (expt 2 1/3) b))) (zerop (* 3.0 (* b p3))))
		  (* -1.0 p1)
		(if (zerop (* 3.0 (* b p3)))
		    (- (/ p3 (* 3.0 (* (expt 2 1/3) b))) p1)
		  (if (zerop (* 3.0 (* (expt 2 1/3) b))) 
		      (* 1.0 (+ p1 (/ (* (expt 2 1/3) p2) (* 3.0 (* b p3)))))
		    (- (/ p3 (* 3.0 (* (expt 2 1/3) b))) (+ p1 (/ (* (expt 2 1/3) p2) (* 3.0 (* b p3)))))))))
	 (bs2 (if (and (zerop (* 6.0 (* (expt 2.0 1/3) b))) (zerop (* 3.0 (* (expt 2.0 2/3) (* b p3))))) 
		  (* -1.0 p1)
		(if (zerop (* 3.0 (* (expt 2.0 2/3) (* b p3))))
		    (* -1.0 (+ p1 (/ (* (complex 1.0 (* -1.0 (sqrt 3.0))) p3) (* 6.0 (* (expt 2.0 1/3) b)))))
		  (if (zerop (* 6.0 (* (expt 2.0 1/3) b)))
		      (- (/ (* (complex 1 (sqrt 3.0)) p2) (* 3.0 (* (expt 2.0 2/3) (* b p3)))) p1)
		    (- (/ (* (complex 1 (sqrt 3.0)) p2) (* 3.0 (* (expt 2.0 2/3) (* b p3)))) (+ p1 (/ (* (complex 1.0 (* -1.0 (sqrt 3.0))) p3) (* 6.0 (* (expt 2.0 1/3) b)))))))))
	 (bs3 (if (and (zerop (* 6.0 (* (expt 2.0 1/3) b))) (zerop (* 3.0 (* (expt 2.0 2/3) (* b p3)))))
		  (* -1.0 p1)
		(if (zerop (* 3.0 (* (expt 2.0 2/3) (* b p3))))
		    (* -1.0 (+ p1 (/ (* (complex 1.0 (sqrt 3.0)) p3) (* 6.0 (* (expt 2.0 1/3) b)))))
		  (if (zerop (* 6.0 (* (expt 2.0 1/3) b)))
		      (- (/ (* (complex 1 (* -1.0 (sqrt 3.0))) p2) (* 3.0 (* (expt 2.0 2/3) (* b p3)))) p1)
		    (- (/ (* (complex 1 (* -1.0 (sqrt 3.0))) p2) (* 3.0 (* (expt 2.0 2/3) (* b p3)))) (+ p1 (/ (* (complex 1.0 (sqrt 3.0)) p3) (* 6.0 (* (expt 2.0 1/3) b))))))))))
    (when (< (abs (imagpart bs1)) max)
      (setf bs1 (realpart bs1)))
    (when (< (abs (imagpart bs2)) max)
      (setf bs2 (realpart bs2)))
    (when (< (abs (imagpart bs3)) max)
      (setf bs3 (realpart bs3)))
    (list bs1 bs2 bs3)))

(defun solve-quartic (a b c d e)
  (let* ((p1 (/ b (* 4.0 a)))
	 (p2 (- (+ (* c c) (* 12.0 (* a e))) (* 3.0 (* b d))))
	 (p3 (- (+ (* 2.0 (* c (* c c))) (+ (* 27.0 (* a (* d d))) (* 27.0 (* e (* b b))))) (+ (* 9.0 (* b (* c d))) (* 72.0 (* a (* c e))))))
	 (p4 (- (/ (* 4.0 (* b c)) (* a a)) (+ (/ (* b (* b b)) (* a (* a a))) (/ (* 8.0 d) a))))
	 (p5 (- (/ (* b b) (* 4.0 (* a a))) (/ (* 2.0 c) (* 3.0 a))))
	 (p6 (- (/ (* b b) (* 2.0 (* a a))) (/ (* 4.0 c) (* 3.0 a))))
	 (p10 (expt (+ p3 (sqrt (- (* p3 p3) (* 4.0 (* p2 (* p2 p2)))))) 1/3))
	 (p7 (/ (* (expt 2.0 1/3) p2) (* 3.0 (* a p10))))
	 (p8 (/ p10 (* 3.0 (* (expt 2.0 1/3) a))))
	 (p9 (* 4.0 (sqrt (+ p5 (+ p7 p8)))))
	 (p11 (* 0.5 (sqrt (+ p5 (+ p7 p8)))))
	 (p12 (- p6 (+ p7 p8)))
	 (bs1 nil)
	 (bs2 nil)
	 (bs3 nil)
	 (bs4 nil))
    (if (zerop p9)
	(setq
	 bs1 (* -1.0 (+ p1 (+ p11 (* 0.5 (sqrt (- p12 0.0))))))
	 bs2 (- (* 0.5 (sqrt (- p12 0.0))) (+ p1 p11))
	 bs3 (- p11 (+ p1 (* 0.5 (sqrt (+ p12 0.0)))))
	 bs4 (- (+ p11 (* 0.5 (sqrt (+ p12 0.0)))) p1))
      (setq
	 bs1 (* -1.0 (+ p1 (+ p11 (* 0.5 (sqrt (- p12 (/ p4 p9)))))))
	 bs2 (- (* 0.5 (sqrt (- p12 (/ p4 p9)))) (+ p1 p11))
	 bs3 (- p11 (+ p1 (* 0.5 (sqrt (+ p12 (/ p4 p9))))))
	 bs4 (- (+ p11 (* 0.5 (sqrt (+ p12 (/ p4 p9))))) p1)))
    (list bs1 bs2 bs3 bs4)))

(defun mult (f a b c d)
    (* f (* a (* b (* c d)))))

(defun sum (lst)
  (let ((sumv 0.0d0))
    (loop for x in lst
	  do (setq sumv (+ sumv  (apply #'mult x))))
    sumv))


(defun check-diamond-inbh (d bh)
  (let* ((tbh (current-time))
         (center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
	 (radius (bh-radius (global-time (current-time))))
         (al nil)
         (bl nil)
         (ar nil)
         (br nil)
         (disr nil)
         (disl nil))
    
    (setf ar (car (multiple-value-list (find-right-edge-position d tbh))))
    (setf br (car (cdr (multiple-value-list (find-right-edge-position d tbh)))))
    (setf al (car (multiple-value-list (find-left-edge-position d tbh))))
    (setf bl (car (cdr (multiple-value-list (find-left-edge-position d tbh)))))

    (setq disr (3vector-distance (diamond-position d :a ar :b br) center))
    (setq disl (3vector-distance (diamond-position d :a al :b bl)  center))

    
    (if (or
         (and (< disr radius) (< disl radius))
         (and (eq (diamond-e d) :BH) (< disl radius))
         (and (eq (diamond-w d) :BH) (< disr radius))
         )
        t
      nil)))
    
(defun check-diamond-inbh-out (d bh)
  (let* ((tbh (current-time))
         (center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
	 (al nil)
         (bl nil)
         (ar nil)
         (br nil)
         (disr nil)
         (disl nil))

    (setf ar (car (multiple-value-list (find-right-edge-position d tbh))))
    (setf br (car (cdr (multiple-value-list (find-right-edge-position d tbh)))))
    (setf al (car (multiple-value-list (find-left-edge-position d tbh))))
    (setf bl (car (cdr (multiple-value-list (find-left-edge-position d tbh)))))

    (setq disr (3vector-distance (diamond-position d :a ar :b br) center))
    (setq disl (3vector-distance (diamond-position d :a al :b bl)  center))

    (format t "DL: ~S DR: ~S~%" disl disr)))
   


(defun check-diamond-inbh-old (d bh)
  (let* ((precision 1E-10)
	 (tbh (current-time))
	 (center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
	 (radius (bh-radius (global-time (current-time))))
	 (start (diamond-start d))
	 (left (diamond-left d))
	 (right (diamond-right d))
	 (end (diamond-end d))
	 (alc nil)
	 (arc nil)
	 (blc nil)
	 (brc nil)
	 (al nil)
	 (bl nil)
	 (ar nil)
	 (br nil)
	 (disr nil)
	 (disl nil))

    (if (> precision (abs (- tbh (4vector-t start))))
	(setq alc 0.0
	      arc 0.0
	      blc 0.0
	      brc 0.0)
      (if (and (> precision (abs (- tbh (4vector-t left)))) (> precision (abs (- tbh (4vector-t right)))))
	  (setq alc 1.0
		blc 0.0
		arc 0.0
		brc 1.0)
	(if (> precision (abs (- tbh (4vector-t left))))
	    (setq alc 1.0
		  blc 0.0
		  arc (/ (- (+ tbh (* 1.0 (4vector-t start))) (+ (4vector-t start) (* 1.0 (4vector-t right)))) (- (+ (4vector-t left) (+ (* 1.0 (4vector-t start)) (* 1.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 1.0 (4vector-t left)) (* 1.0 (4vector-t right))))))
		  brc (/ (- (+ tbh (* 0.0 (4vector-t start))) (+ (4vector-t start) (* 0.0 (4vector-t left)))) (- (+ (4vector-t right) (+ (* 0.0 (4vector-t start)) (* 0.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 0.0 (4vector-t left)) (* 0.0 (4vector-t right)))))))
	  (if (> precision (abs (- tbh (4vector-t right))))
	      (setq alc (/ (- (+ tbh (* 0.0 (4vector-t start))) (+ (4vector-t start) (* 0.0 (4vector-t right)))) (- (+ (4vector-t left) (+ (* 0.0 (4vector-t start)) (* 0.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 0.0 (4vector-t left)) (* 0.0 (4vector-t right))))))
		    blc (/ (- (+ tbh (* 1.0 (4vector-t start))) (+ (4vector-t start) (* 1.0 (4vector-t left)))) (- (+ (4vector-t right) (+ (* 1.0 (4vector-t start)) (* 1.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 1.0 (4vector-t left)) (* 1.0 (4vector-t right))))))
		    arc 0.0
		    brc 1.0)
	    (if (> precision (abs (- tbh (4vector-t end))))
		(setq alc 1.0
		      arc 1.0
		      blc 1.0
		      brc 1.0)
	      (setq arc (/ (- (+ tbh (* 1.0 (4vector-t start))) (+ (4vector-t start) (* 1.0 (4vector-t right)))) (- (+ (4vector-t left) (+ (* 1.0 (4vector-t start)) (* 1.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 1.0 (4vector-t left)) (* 1.0 (4vector-t right))))))
		    brc (/ (- (+ tbh (* 0.0 (4vector-t start))) (+ (4vector-t start) (* 0.0 (4vector-t left)))) (- (+ (4vector-t right) (+ (* 0.0 (4vector-t start)) (* 0.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 0.0 (4vector-t left)) (* 0.0 (4vector-t right))))))
		    alc (/ (- (+ tbh (* 0.0 (4vector-t start))) (+ (4vector-t start) (* 0.0 (4vector-t right)))) (- (+ (4vector-t left) (+ (* 0.0 (4vector-t start)) (* 0.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 0.0 (4vector-t left)) (* 0.0 (4vector-t right))))))
		    blc	(/ (- (+ tbh (* 1.0 (4vector-t start))) (+ (4vector-t start) (* 1.0 (4vector-t left)))) (- (+ (4vector-t right) (+ (* 1.0 (4vector-t start)) (* 1.0 (4vector-t end)))) (+ (4vector-t start) (+ (* 1.0 (4vector-t left)) (* 1.0 (4vector-t right))))))))))))


    (if (and (zerop alc) (zerop blc))
	(setq al 0.0
	      bl 0.0)
      (if (and (eq alc 1.0) (eq blc 1.0))
	  (setq al 1.0
		bl 1.0)
	(if (and (>= alc 0.0) (<= alc 1.0))
	    (setq al alc
		  bl 0.0)
	  (if (and (>= blc 0.0) (<= blc 1.0))
	      (setq al 1.0
		    bl blc)))))

    (if	(and (zerop arc) (zerop brc))
        (setq ar 0.0
              br 0.0)
      (if (and (eq arc 1.0) (eq brc 1.0))
          (setq	ar 1.0
                br 1.0)
	(if (and (>= arc 0.0) (<= arc 1.0))
	    (setq ar arc
		  br 1.0)
	  (if (and (>= brc 0.0) (<= brc 1.0))
	      (setq ar 0.0
		    br brc)))))
    
      
    (setq disr (3vector-distance (diamond-position d :a ar :b br) center))
    (setq disl (3vector-distance (diamond-position d :a al :b bl)  center))
        
    (if (or
	 (and (< disr (+ radius precision)) (< disl (+ radius precision)))
	 (and (eq (diamond-e d) :BH) (< disl (+ radius precision)))
	 (and (eq (diamond-w d) :BH) (< disr (+ radius precision)))
	 )
	t
      nil)))
    


;This function receives a diamond and perform a intersection on it. The function only receives a diamond because
;due to an earlier intersection it could be the case that the intersection in this diamond is not more there. This is
;the reason why the first thing that the function does is to check for intersections
(defun perform-bh-intersec (int)
  (format t "#########################~%")
  (format t "Performing int~%")
  (let* ((err 1E-8)
	 (d (bh-intersec-diamond int))
	 (bh (bh-intersec-bh int))
	 (int-now nil)
	 (a (bh-intersec-a int))
	 (b (bh-intersec-b int))
	 (intersection-point (bh-intersec-spacetime int))
	 (global-point (globalize-position intersection-point))
	 (center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
	 (e nil)
	 (w nil)
	 lrv
	 rrv
	 cpv
	 ar
	 br
	 al
	 bl
	 (radius (bh-radius (global-time (current-time))))
	 distance
	 (inner nil)
	 (innerc nil)
	 (int (check-bh-intersec-diamond d bh nil)))


    
    (if int
	(setf int-now int)
      (if (or (eq (diamond-e d) :BH) (eq (diamond-w d) :BH)
	      (eq (diamond-e d) :BHdeleted) (eq (diamond-w d) :BHdeleted))
	  (if (or (eq (diamond-e d) :BH) (eq (diamond-e d) :BHdeleted))
	      (setf int-now (check-bh-intersec-diamond (diamond-w d) bh nil))
	    (setf int-now (check-bh-intersec-diamond (diamond-e d) bh nil)))
	(setf int-now (check-bh-intersec-diamond (diamond-w d) bh nil))))



    (when (null (possible-int d bh))
      (if (or (> (length int-now) 1)
	      (eq (diamond-e d) :BH) (eq (diamond-w d) :BH)
	      (eq (diamond-e d) :BHdeleted) (eq (diamond-w d) :BHdeleted))
	  (warn "Double Int")
	(progn
	  (warn "Int in strange place")
	  (return-from perform-bh-intersec nil))))  
    
    (when (> (length int-now) 1)
      (warn "Double Int"))

   
    
    (if (> (length int-now) 1)
	(if (< (3vector-distance intersection-point (bh-intersec-spacetime (car int-now))) err)
	    (setf int-now (car int-now))
	  (setf int-now (car (cdr int-now))))
      (setf int-now (car int-now)))




     (when (null int-now)
       (warn "NULL INT~%")
      (return-from perform-bh-intersec nil))


     (when (null (diamond-start (bh-intersec-diamond int-now)))
        (format t "NULL D~%")
        (return-from perform-bh-intersec nil))
     
     (when (and (or (eq (diamond-w (bh-intersec-diamond int-now)) :BH) (eq (diamond-e (bh-intersec-diamond int-now)) :BH)
		    (eq (diamond-w (bh-intersec-diamond int-now)) :BHeatit) (eq (diamond-e (bh-intersec-diamond int-now)) :BHeatit))
		(or (eq *era* :radiation)
		    (eq *era* :radiation-smooth)
		    (eq *era* :matter)
		    (eq *era* :power))
		)
      (format t "Int in BH~%")
      (return-from perform-bh-intersec nil))

     (setf a (bh-intersec-a int-now)
	   b (bh-intersec-b int-now)
	   d (bh-intersec-diamond int-now)
	   bh (bh-intersec-bh int-now)
	   center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d)))
     
     (setf distance (3vector-distance (diamond-position d :a a :b b) center))
     (format t "Dis: ~S~%" distance)
     (format t "R: ~S~%" radius)
     (format t "P: ~S~%" (globalize-position (diamond-position d :a a :b b)))

   
     (when (or (> (/ distance radius) (+ 1 1E-3))
	       (< (/ distance radius) (- 1 1E-3)))
       (warn "FICTICIUS INT~%")
       (return-from perform-bh-intersec nil))

     (format t "Possible Place: ~S~%" (possible-int d bh))

     
     (setf ar (car (multiple-value-list (find-right-edge-position d (4vector-t (diamond-position d :a a :b b))))))
     (setf br (car (cdr (multiple-value-list (find-right-edge-position d (4vector-t (diamond-position d :a a :b b)))))))
     (setf al (car (multiple-value-list (find-left-edge-position d (4vector-t (diamond-position d :a a :b b))))))
     (setf bl (car (cdr (multiple-value-list (find-left-edge-position d (4vector-t (diamond-position d :a a :b b)))))))

     (setf cpv (3vector-normalize (3vector- (standardize-position (localize-position (blackhole-center bh)) (diamond-position d :a a :b b)) (diamond-position d :a a :b b))))
     (setf lrv (3vector-normalize (3vector- (diamond-position d :a al :b bl) (diamond-position d :a a :b b))))
     (setf rrv (3vector-normalize (3vector- (diamond-position d :a ar :b br) (diamond-position d :a a :b b))))
     
     (when (< (abs (3vector-dot lrv cpv)) 0.01)
       (format t "Small Angle~%"))
     
     (when (< (abs (3vector-dot lrv cpv)) 1E-3)
       (format t "DE: ~S~%" (3vector-distance (diamond-position d :a ar :b br) center))
       (format t "DW: ~S~%" (3vector-distance (diamond-position d :a al :b bl) center))
       (format t   "Surface parallel to string ~S ~%" (3vector-dot lrv cpv))
       (return-from perform-bh-intersec nil))
     
    
    (setf a (bh-intersec-a int-now)
	  b (bh-intersec-b int-now)
	  d (bh-intersec-diamond int-now)
	  bh (bh-intersec-bh int-now)
	  center (standardize-position (localize-position (blackhole-center bh)) (diamond-start d))
	  intersection-point (bh-intersec-spacetime int-now)
	  global-point (globalize-position intersection-point)
	  e (east-quarter d intersection-point global-point :a a :b b)
	  w (west-quarter d intersection-point global-point :a a :b b))


    (when (and (null *all-bhs*)
	       (< (abs (- (local-time *bh-start*) (current-time))) 1E-5))
      (unless (gethash (blackhole-center bh) *my-bhs*)
	(format t "Pushing to hash~%")
	(setf (gethash (blackhole-center bh) *my-bhs*) bh)))
         
    (when (eq (diamond-w w) :BH)
      (setf (diamond-nw d) :BH))

    (when (eq (diamond-e e) :BH)
      (setf (diamond-ne d) :BH))
   
  
    
    (format t "lrv: ~S~%" lrv)
    (format t "rrv: ~S~%" rrv)
    (format t "cpv: ~S~%" cpv)
    
    
    (when (>= (3vector-dot lrv cpv) 0.0)
      (setf inner t))

    (when (< (3vector-distance (diamond-position d :a al :b bl) center) (3vector-distance (diamond-position d :a a :b b) center))
      (setf innerc t))

        
 ;   (format t "Inner: ~S~%" (3vector-dot lrv cpv))
    (format t "DisE: ~S~%" (3vector-distance (diamond-position d :a ar :b br) center))
    (format t "DisW: ~S~%" (3vector-distance (diamond-position d :a al :b bl) center))
    
    (format t "L-Cos: ~S~%" (3vector-dot lrv cpv))
    (format t "R-Cos: ~S~%" (3vector-dot rrv cpv))
    (divide-pending-intersections d :a a :b b :east e :west w)
    (rescale-right-junctions d e :a a :b b)
    ;(rescale-left-junctions d e :a a :b b)
    (rescale-left-junctions d w :a a :b b)
    ;(rescale-right-junctions d w :a a :b b)
    (format t "Neighs w ~S~%" (diamond-w w))
;    (format t "WW: ~S~%" (diamond-w (diamond-w w)))
    (format t "Neighs e ~S~%" (diamond-e e))
;    (format t "EE: ~S~%" (diamond-e (diamond-e e)))
    ;(format t "Neighs e ~S~%" (diamond-ne e))
    (format t "a: ~S~%" a)
    (format t "b: ~S~%" b)
 ;   (format t "GP: ~S~%" (globalize-position (diamond-position d :a a :b b)))
 ;   (format t "Mine: ~S~%" (point-mine (diamond-position d :a a :b b)))
    ;(format t "Pos: ~S~%" (diamond-position d :a a :b b))
    ;(format t "T: ~S~%" (current-time))
    ;(format t "e: ~S~%" (diamond-end e))
    ;(format t "w: ~S~%" (diamond-end w))
    ;(format t "T: ~S~%" (current-time)) 
    (discard-object d)
    (propagate-cut-ne e a)
    (propagate-cut-nw w b)
    (handle-new-diamond e :predecessor :meiosis)
    (handle-new-diamond w :predecessor :meiosis)
    (advance-cut-point-BH e  w global-point bh inner)
    (when (or (and (null inner) innerc)
              (and inner (null innerc)))
      (warn "Bh deletion going out"))
;    (format t "#####################~%")
    ))

;This function evolved the intersection point. It receives the east and the west diamonds from the intersection point    
;and creates two new diamonds. This two diamonds are north of the intersection point and divided by a imaginary line    
;which is a BH.
(defun advance-cut-point-BH (east west global-start bh inner)
  (when (null bh)
    (error "NULL BH!!!!"))
  (let* ((left (diamond-end west))
	 (right (diamond-end east))
	 (start (diamond-left east))
	 (tag-position (3to4vector (4to3vector (blackhole-center bh)) (4vector-t global-start)))
	 (bh-size (blackhole-size bh))
	 (new-bh (make-blackhole :center tag-position :size bh-size))
	 (north-east (make-diamond :start start
				   :right right
				   :left  (3to4vector (4to3vector (3vector- start (3vector- right start))) (4vector-t right))
				   :nw    (if *pointbh*
					      :BH
					    (if inner
						:BH
					      :BHdeleted))
				   :sw    nil
				   :se    east  
				   :tag   (create-bh-tag global-start new-bh)
				   :a-kink-created global-start :b-kink-created global-start
				   :bh bh
				   ))
	 (north-west (make-diamond :start start
				   :left  left
				   :right (3to4vector (4to3vector (3vector- start (3vector- left start))) (4vector-t left))
				   :ne    (if *pointbh*
					      :BH
					      (if inner
						  :BHdeleted
						:BH))
				   :se    nil
				   :sw    west
				   :tag   (create-bh-tag global-start new-bh)
				   :a-kink-created global-start :b-kink-created global-start
				   :bh bh
				   )))
    (if (or (eq *era* :radiation)
	    (eq *era* :radiation-smooth)
	    (eq *era* :matter)
	    (eq *era* :power))
        (progn
          (setf (diamond-end north-east) (compute-bh-diamond-end north-east right))
          (setf (diamond-end north-west) (compute-bh-diamond-end north-west left))
	  (setf (diamond-left north-east) (4vector- (4vector+ start (diamond-end north-east)) right))
	  (setf (diamond-right north-west) (4vector- (4vector+ start (diamond-end north-west)) left)))
      (progn
        (setf (diamond-end north-east) (compute-bh-diamond-end north-east right))
        (setf (diamond-end north-west) (compute-bh-diamond-end north-west left))))
    
    (setf (diamond-nw east) north-east)
    (setf (diamond-ne west) north-west)
    (handle-new-diamond north-east)
    (handle-new-diamond north-west)
    ;(format t "#######################~%")
    (format t "START-G: ~S~%" (globalize-position start))
    (format t "Start: ~S~%" start)
    (format t "BH: ~S~%" (blackhole-center bh))
    (format t "NE: ~D~%" north-east)
    (format t "TSNE: ~S~%" (4vector-t (diamond-start north-east)))
    (format t "TENE: ~S~%" (4vector-t (diamond-end north-east)))
    ;(format t "DSR: ~S~%" (3vector-distance start right))
    ;(format t "TSR: ~S~%" (- (4vector-t right) (4vector-t start)))
    ;(format t "DRE: ~S~%" (3vector-distance right (diamond-end north-east)))
    ;(format t "TRE: ~S~%" (- (4vector-t (diamond-end north-east)) (4vector-t right)))
    (format t "NON: ~S~%" (diamond-nw north-east))
    ;(format t "DSL: ~S~%" (3vector-distance start left))
    ;(format t "TSL: ~S~%" (- (4vector-t left) (4vector-t start)))
    ;(format t "DLE: ~S~%" (3vector-distance left (diamond-end north-west)))
    ;(format t "TLE: ~S~%" (- (4vector-t (diamond-end north-west)) (4vector-t left)))
    (format t "NW: ~D~%" north-west)
    (format t "TSNW: ~S~%" (4vector-t (diamond-start north-west)))
    (format t "TENW: ~S~%" (4vector-t (diamond-end north-west)))
    
    (format t "NON: ~S~%"(diamond-ne north-west))

    (format t "#######################~%")
    north-east
    north-west))



(defun check-new-dump (t1 t2 dt dir1 dir2 nbh)
  (let ((old nil)
	(new nil)
	(min 1.0)
	(ss1 nil)
	(ss2 nil))
    (loop for time from t1 to t2 by dt
	  do (setf new nil)
	  do (setf *read-dump-diamonds* t)
	  do (setf old (print-bh-stats dir1 time time dt nbh))
	  do (setf *read-dump-diamonds* nil)
	  do (read-dumps dir2 :time time)
	  ;do (loop for s1 in old
		   ;do (format t "O: ~S~%" s1))
	  do (loop for seg in *dumped-segments*
		do (when (and (or (eq (dump-segment-start-junction seg) :BH)
				  (eq (dump-segment-start-junction seg) :BHeatit))
			      (or (eq (dump-segment-end-junction seg) :BH)
				  (eq (dump-segment-end-junction seg) :BHeatit)))
		     (push seg new)))
	  do (format t "T: ~S~%" time)
	  do (format t "LO: ~S~%" (length old))
	  do (format t "LN: ~S~%" (length new))
	  do (loop for s1 in old
		   do (setf min 1.0)
		   do (setf ss1 s1)
		   do (loop for s2 in new
			    do (when (and (or (and (zerop (3vector-distance (blackhole-center (car s1)) (blackhole-center (dump-segment-start-bh s2))))
						   (zerop (3vector-distance (blackhole-center (car (cdr s1))) (blackhole-center (dump-segment-end-bh s2)))))
					      (and (zerop (3vector-distance (blackhole-center (car s1)) (blackhole-center (dump-segment-end-bh s2))))
						   (zerop (3vector-distance (blackhole-center (car (cdr s1))) (blackhole-center (dump-segment-start-bh s2))))))
					  (< (abs (- (/ (car (cdr (cdr s1))) (dump-segment-length s2)) 1.0)) min))
				 (setf min (abs (- (/ (car (cdr (cdr s1))) (dump-segment-length s2)) 1.0)))
				 (setf ss2 s2)))
		   ;do (format t "Min: ~S~%" min)
		   do (when (> min 1E-3)
			(warn "Not match in time ~S with L1: ~S and L2: ~S~%" time (car (cdr (cdr ss1))) (dump-segment-length ss2))))
	  do (when (null (zerop (- (length old) (length new))))
	       (error "Different number of segments in time ~S~%" time)))
    ))



(defun check-new-read (t1 t2 dt dir)
  (loop for time from t1 to t2 by dt
	do (setf *read-dump-diamonds* nil)
	do (read-dumps dir :time time)))
