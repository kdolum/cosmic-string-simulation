(in-package "CL-USER")
;;This is the same as what's in Numerical Recipes in Lisp, except:
;; (optimize speed) declaration
;; Instead of creating INDX (output array of N bytes for index permutation
;; and VV (working array of size N) they must be supplied.
;; Byte size in index arrays is 16 instead of 8.
;; We don't maintain or return the perturbation sign D.  On success we return T.
;; Instead of signaling an error or using tiny = 1e-20 when we encounter
;; singularities, we just return NIL.

(defun ludcmp (a indx vv)
 (declare (type (simple-array double-float (* *)) a)
	  (type (simple-array (unsigned-byte 16) (*)) indx)
	  (type (simple-array double-float (*)) vv))
 (prog* (
  (n (array-dimension a 0))
  (imax 0) (aamax 0d0)
  (sum 0d0) (dum 0d0))

;;  (declare (optimize speed))		;Declare only here to suppress warning about return float
  (declare (type fixnum imax))
  (declare (type double-float aamax sum dum))
 
  (do ((i 0 (+ i 1)))
      ((> i (1- n)) t)
      (declare (type fixnum i))
    (setf aamax 0d0)
    (do ((j 0 (+ j 1)))
        ((> j (1- n)) t)
        (declare (type fixnum j))
      (if (> (abs (aref a i j)) aamax) 
          (setf aamax (abs (aref a i j)))))
    (if (= aamax 0d0) (return-from ludcmp nil))	;Zero column
    (setf (aref vv i) (/ 1d0 aamax))) 

  (do ((j 0 (+ j 1)))
      ((> j (1- n)) t)
      (declare (type fixnum j))
    (do ((i 0 (+ i 1)))
        ((> i (1- j)) t)
        (declare (type fixnum i))
      (setf sum (aref a i j))
      (do ((k 0 (+ k 1)))
          ((> k (1- i)) t)
        (setf sum (- sum (* (aref a i k) (aref a k j)))))
      (setf (aref a i j) sum))

    (setf aamax 0d0)
    (do ((i j (+ i 1)))
        ((> i (1- n)) t)
        (declare (type fixnum i))
      (setf sum (aref a i j))

      (do ((k 0 (+ k 1)))
          ((> k (1- j)) t)
          (declare (type fixnum k))
        (setf sum (- sum (* (aref a i k) (aref a k j)))))

      (setf (aref a i j) sum)
      (setf dum (* (aref vv i) (abs sum)))
      (when 
       (>= dum aamax) 
       (setf imax i)
       (setf aamax dum)))

    (when 
     (not (= j imax))
     (do ((k 0 (+ k 1)))
         ((> k (1- n)) t)
         (declare (type fixnum k))
       (setf dum (aref a imax k))
       (setf (aref a imax k) (aref a j k))
       (setf (aref a j k) dum))
     (setf (aref vv imax) (aref vv j)))

    (setf (aref indx j) imax)
    (if (= (aref a j j) 0) (return-from ludcmp nil))
    (when 
     (not (= (1- n) j))
     (setf dum (/ 1 (aref a j j)))
     (do ((i (1+ j) (+ i 1)))
         ((> i (1- n)) t)
         (declare (type fixnum i))
       (setf (aref a i j) (* (aref a i j) dum))))) 
  (return t)))

;-----------------------------------------------------------------------------
(defun lubksb (a indx b)
 (declare (type (simple-array double-float (* *)) a)) 
 (declare (type (simple-array (unsigned-byte 16) (*)) indx)) ; refers to 0 based array
 (declare (type (simple-array double-float (*)) b)) 

 (prog ((n 0) (ii 0) (sum 0d0) (ll 0))
  (declare (optimize speed))
  (declare (type fixnum n ii ll))
  (declare (type double-float sum))

  (setq n (array-dimension a 0))
  (setf ii 0) 

  (do ((i 1 (+ i 1)))
      ((> i n) t)
      (declare (type fixnum i))
    (setf ll (1+ (aref indx (1- i))))
    (setf sum (aref b (1- ll)))
    (setf (aref b (1- ll)) (aref b (1- i)))
    (cond 
     ((not (= ii 0)) 
      (do ((j ii (+ j 1)))
          ((> j (1- i)) t)
          (declare (type fixnum j))
        (setf sum (- sum (* (aref a (1- i) (1- j)) (aref b (1- j)))))))
     ((not (= sum 0d0))
      (setf ii i)))
    (setf (aref b (1- i)) sum)) 

  (do ((i (1- n) (1- i)))
      ((< i 0) t)
      (declare (type fixnum i))
    (setf sum (aref b i))
    (do ((j (1+ i) (+ j 1)))
        ((> j (1- n)) t)
        (declare (type fixnum j))
      (setf sum (- sum (* (aref a i j) (aref b j)))))
    (setf (aref b i) (/ sum (aref a i i)))) 
   
  (return b)))

