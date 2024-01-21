(in-package "CL-USER")

;;This is exactly the same as in Numerical Recipes in Lisp except:
;; (optimize speed) declaration
;; (floor (/ n 2)) -> (floor n 2), etc.
;; more type declarations
;; (optimize (safety 0)) to avoid check-fixnum calls
(defun four1 (data nn &key (isign 1))
 (declare (type (simple-array double-float (*)) data)
	  (type fixnum nn)
	  (type (integer -1 1) isign)
	  (optimize speed))
 (locally (declare (optimize (safety 0)))	;Avoid checking that various numbers still are fixnums
   (prog ((wr 0d0) (wi 0d0) (wpr 0d0) (wpi 0d0) (wtemp 0d0) 
	  (theta 0d0) (tempr 0d0) (tempi 0d0) (j 0) (n 0) (m 0) 
	  (mmax 0) (istep 0))
	 (declare (type double-float wr wi wpr wpi wtemp theta tempr tempi)) 
	 (declare (type fixnum j n m mmax istep))

	 (setf n (* 2 nn)) 
	 (setf j 1) 
	 (do ((i 1 (+ i 2)))
	     ((> i n) t)
	   (declare (type fixnum i))
	   (when (> j i) 
	     (setf tempr (aref data (1- j)))
	     (setf tempi (aref data j)) 
	     (setf (aref data (1- j)) (aref data (1- i)))
	     (setf (aref data j) (aref data i)) 
	     (setf (aref data (1- i)) tempr)
	     (setf (aref data i) tempi))
	   (setf m (floor n 2))
	   label1
	   (when (and (>= m 2) (> j m))
	     (setf j (- j m)) (setf m (floor m 2))
	     (go label1))
	   (setf j (+ j m))) 

	 (setf mmax 2) 
	 label2 
	 (when (> n mmax)
	   (setf istep (* 2 mmax))
	   (setf theta (/ 6.28318530717959d0 (* isign mmax)))
	   (setf wpr (* -2.0d0 (expt (sin (* 0.5d0 theta)) 2)))
	   (setf wpi (sin theta)) (setf wr 1.0d0) (setf wi 0.0d0)
	   (do ((m 1 (+ m 2)))
	       ((> m mmax) t)
	     (declare (type fixnum m))
	     (do ((i m (+ i istep)))
		 ((> i n) t)
	       (declare (type fixnum i))
	       (setf j (+ i mmax))
	       (setf tempr (+ (* wr (aref data (1- j)))
			      (* (* -1d0 wi) (aref data j))))
	       (setf tempi (+ (* wr (aref data j))
			      (* wi (aref data (1- j)))))
	       (setf (aref data (1- j)) (+ (aref data (1- i)) (- tempr)))
	       (setf (aref data j) (+ (aref data i) (* -1d0 tempi)))
	       (setf (aref data (1- i)) (+ (aref data (1- i)) tempr))
	       (setf (aref data i) (+ (aref data i) tempi)))
	     (setf wtemp wr)
	     (setf wr (+ (+ (* wr wpr) (* (- wi) wpi)) wr))
	     (setf wi (+ (+ (* wi wpr) (* wtemp wpi)) wi)))
	   (setf mmax istep)
	   (go label2)) 
   
	 (return data))))

(defun realft (data n &key (isign 1))
 (declare (type (simple-array double-float (*)) data)
	  (type fixnum n isign)
	  (optimize speed))

 (prog ((wr 0d0) (wi 0d0) (wpr 0d0) (wpi 0d0) (wtemp 0d0) 
        (theta 0d0) (c1 0d0) (n2p3 0) (c2 0d0) (i1 0) (i2 0) (i3 0) (i4 0)
        (wrs 0d0) (wis 0d0) (h1r 0d0) (h1i 0d0) (h2r 0d0) (h2i 0d0))
  (declare (type double-float wr wi wpr wpi wtemp theta c1 c2
                              wrs wis h1r h1i h2r h2i))
  (declare (type fixnum  n2p3 i1 i2 i3 i4))

 (locally (declare (optimize (safety 0))) ;Avoid checking that various numbers still are fixnums

					;  (setq n (/ (array-dimension data 0) 2))
 
   (setf theta (/ 3.141592653589793d0 (dfloat n))) 
   (setf c1 0.5d0) 
   (cond 
    ((= isign 1) 
     (setf c2 -0.5d0)
     (setq data (four1 data n :isign 1)))
    (t 
     (setf c2 0.5d0)
     (setf theta (- theta)))) 
   (setf wpr (* (- 2.0d0) (expt (sin (* 0.5d0 theta)) 2)))
 
   (setf wpi (sin theta)) 
   (setf wr (1+ wpr)) 
   (setf wi wpi) 
   (setf n2p3 (+ (the fixnum (* 2 n)) 3)) 

   (do ((i 2 (+ i 1)))
       ((> i (floor n 2)) t)
     (declare (type fixnum i))
     (setf i1 (1- (* 2 i)))
     (setf i2 (+ i1 1))
     (setf i3 (+ n2p3 (- i2)))
     (setf i4 (1+ i3))

     (setf wrs wr)
     (setf wis wi)
     (setf h1r (* c1 (+ (fref data i1) (fref data i3))))
     (setf h1i (* c1 (- (fref data i2) (fref data i4))))
     (setf h2r (* (- c2) (+ (fref data i2) (fref data i4))))
     (setf h2i (* c2 (- (fref data i1) (fref data i3))))
     (setf (fref data i1) (- (+ h1r (* wrs h2r)) (* wis h2i)))
     (setf (fref data i2) (+ (+ h1i (* wrs h2i)) (* wis h2r)))
     (setf (fref data i3) (+ (- h1r (* wrs h2r)) (* wis h2i)))
     (setf (fref data i4) (+ (+ (- h1i) (* wrs h2i)) (* wis h2r)))
     (setf wtemp wr)
     (setf wr (+ (+ (* wr wpr) (* (- wi) wpi)) wr))
     (setf wi (+ (+ (* wi wpr) (* wtemp wpi)) wi))) 

   (cond 
    ((= isign 1) 
     (setf h1r (fref data 1))
     (setf (fref data 1) (+ h1r (fref data 2)))
     (setf (fref data 2) (+ h1r (- (fref data 2)))))
    (t
     (setf h1r (fref data 1)) 
     (setf (fref data 1) (* c1 (+ h1r (fref data 2))))
     (setf (fref data 2) (* c1 (+ h1r (- (fref data 2)))))
     (setq data (four1 data n :isign -1)))) 
   
   (return data))))

(defun cosft (y &key (isign 1))
 (declare (type (simple-array double-float (*)) y))
 (declare (type fixnum isign))
 
 (prog ((n 0) (wr 0d0) (wi 0d0) (wpr 0d0) (wpi 0d0) (wtemp 0d0)
        (theta 0d0) (sum 0d0) (m 0d0) (y1 0d0) (y2 0d0) (even 0d0) 
        (odd 0d0) (enf0 0d0) (sumo 0d0) (sume 0d0))
  (declare (type double-float wr wi wpr wpi wtemp theta 
                              y1 y2 even odd enf0 sumo sume))
  (declare (type fixnum n)) 

  (setq n (array-dimension y 0))

  (setf theta (/ 3.141592653589793d0 (dfloat n))) 
  (setf wr 1.0d0) 
  (setf wi 0.0d0) 
  (setf wpr (* -2.0d0 (expt (sin (* 0.5d0 theta)) 2))) 
  (setf wpi (sin theta)) 
  (setf sum (fref y 1)) 
  (setf m (floor (/ n 2))) 
  (do ((j 1 (+ j 1)))
      ((> j (- m 1)) t)
      (declare (type fixnum j))
    (setf wtemp wr)
    (setf wr (+ (+ (* wr wpr) (* (- wi) wpi)) wr))
    (setf wi (+ (+ (* wi wpr) (* wtemp wpi)) wi))
    (setf y1 (* 0.5d0 (+ (fref y (+ j 1)) (fref y (+ (- n j) 1)))))
    (setf y2 (+ (fref y (+ j 1)) (- (fref y (+ (- n j) 1)))))
    (setf (fref y (+ j 1)) (+ y1 (* (- wi) y2)))
    (setf (fref y (+ (- n j) 1)) (+ y1 (* wi y2)))
    (setf sum (+ sum (* wr y2)))) 

  (setq y (realft y (/ n 2) :isign 1)) 
  (setf (fref y 2) sum) 

  (do ((j 4 (+ j 2)))
      ((> j n) t)
      (declare (type fixnum j))
    (setf sum (+ sum (fref y j)))
    (setf (fref y j) sum)) 

  (when (= isign -1)
   (setf even (fref y 1)) 
   (setf odd (fref y 2))

   (do ((i 3 (+ i 2)))
       ((> i (- n 1)) t)
      (declare (type fixnum i))
     (setf even (+ even (fref y i)))
     (setf odd (+ odd (fref y (+ i 1)))))

   (setf enf0 (* 2d0 (- even odd)))
   (setf sumo (- (fref y 1) enf0))
   (setf sume (- (/ (* 2d0 odd) (dfloat n)) sumo))
   (setf (fref y 1) (* 0.5d0 enf0)) 
   (setf (fref y 2) (- (fref y 2) sume))

   (do ((i 3 (+ i 2)))
       ((> i (- n 1)) t)
       (declare (type fixnum i))
     (setf (fref y i) (- (fref y i) sumo))
     (setf (fref y (+ i 1)) (- (fref y (+ i 1)) sume)))) 

  (return y)))
