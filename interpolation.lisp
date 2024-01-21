;;Some interpolation routines from numerical recipes
(in-package "CL-USER")

(defun polint (xa ya n x)
 (declare (type (simple-array double-float (*)) xa)) 
 (declare (type (simple-array double-float (*)) ya)) 
 (declare (type double-float x))
 (declare (type fixnum n))

 (prog* ((y 0d0) (dy 0d0) 
        (c (make-array n :element-type 'double-float :initial-element 0d0))
        (d (make-array n :element-type 'double-float :initial-element 0d0))
        (dif 0d0) (dift 0d0)
        (ns 0) (ho 0d0) (hp 0d0) (w 0d0) (den 0d0))
  (declare (type (simple-array double-float (*)) c)) 
  (declare (type (simple-array double-float (*)) d))
  (declare (type double-float y dy dif dift ho hp w den))
  (declare (type fixnum ns))


  (setf ns 1) 
  (setf dif (abs (- x (aref xa 0)))) 
  (do ((i 1 (+ i 1)))
      ((> i n) t)
      (declare (type integer i))
    (setf dift (abs (- x (aref xa (1- i)))))
    (when 
      (< dift dif)
      (setf ns i) 
      (setf dif dift))
    (setf (aref c (1- i)) (aref ya (1- i)))
    (setf (aref d (1- i)) (aref ya (1- i)))) 

  (setf y (aref ya (1- ns))) 
  (setf ns (1- ns)) 
  (do ((m 1 (+ m 1)))
      ((> m (1- n)) t)
      (declare (type integer m))
    (do ((i 1 (+ i 1)))
        ((> i (- n m)) t)
        (declare (type integer i))
      (setf ho (- (aref xa (1- i)) x))
      (setf hp (- (aref xa (1- (+ i m))) x))
      (setf w (- (aref c i) (aref d (1- i))))
      (setf den (- ho hp))
      (if (= den 0d0) (error " den = 0d0 in polint "))
      (setf den (/ w den))
      (setf (aref d (1- i)) (* hp den))
      (setf (aref c (1- i)) (* ho den)))

    (cond 
     ((< (* 2 ns) (- n m)) 
      (setf dy (aref c ns)))
     (t 
      (setf dy (aref d (1- ns)))
      (setf ns (1- ns))))
    (setf y (+ y dy))) 
  
  (return (values y dy))))
