;;Integration routes from Numerical Recipes
(in-package "CL-USER")

;;TRAPZD exactly as in numerical recipes
;;This is a closure, so that it can store the level of refinement
;;previously used and the result.
(defvar trapzd 
  (let ((it 0) (s 0d0))
    (declare (type fixnum it))
    (declare (type double-float s))

    #'(lambda (func a b n)
	(declare (type double-float a b))
	(declare (type fixnum n))

	(prog ((sum 0d0) (tnm 0) (del 0d0) (x 0d0))
	      (declare (type double-float sum del x))
	      (declare (type fixnum tnm))

	      (cond 
	       ((= n 1) 
		(setf s (* (* 0.5d0 (- b a)) 
			   (+ (dfloat-check (funcall func a)) 
			      (dfloat-check (funcall func b)))))
		(setf it 1))
	       (t 
		(setf tnm it)
		(setf del (/ (- b a) tnm)) 
		(setf x (+ a (* 0.5d0 del))) 
		(setf sum 0d0)
		(do ((j 1 (+ j 1)))
		    ((> j it) t)
		  (declare (type fixnum j))
		  (setf sum (+ sum (dfloat-check (funcall func x))))
		  (setf x (+ x del)))
		(setf s (* 0.5d0 (+ s (/ (* (- b a) sum) tnm)))) 
		(setf it (* 2 it)))) 
   
	      (return (the double-float s))))))


(defun qtrap (func a b &key (eps 1.0d-6) (jmax 20)) 
  (declare (type double-float a b eps))
  (declare (type fixnum jmax))

  (prog ((olds 0d0) (s 0d0))
	(declare (type double-float s olds))
	(declare (special trapzd))

	(setf olds -1.0d30) 
	(do ((j 1 (+ j 1)))
	    ((> j jmax) t)
	  (declare (type fixnum j))
	  (setq s (dfloat-check (funcall trapzd func a b j)))
    
	  (if (< (abs (- s olds)) (* eps (abs olds))) (go end))
	  (setf olds s))

	(error " too many steps in qtrap. ")
	end
	(return (the double-float s))
	))

;;Romberg integration exactly as an numerical recipes, except that the default value of
;;jmaxp is jmax+1.
(defun qromb (func a b
		   &key (eps 1.0d-6) (jmax 20) (jmaxp (1+ jmax)) (k 5) (km (1- k)))
  (declare (type double-float a b eps))
  (declare (type fixnum jmax jmaxp k km))

  (prog* (
	  (s (make-array jmaxp :element-type 'double-float :initial-element 0d0))
	  (h (make-array jmaxp :element-type 'double-float :initial-element 0d0))
	  (hj 
	   (make-array  k :element-type 'double-float :initial-element 0d0))
	  (sj 
	   (make-array k :element-type 'double-float :initial-element 0d0))
	  (ss 0d0) (dss 0d0))


	 (declare (type (simple-array double-float (*)) s)) 
	 (declare (type (simple-array double-float (*)) h)) 
	 (declare (type (simple-array double-float (*)) sj)) 
	 (declare (type (simple-array double-float (*)) hj)) 
	 (declare (type double-float ss dss))


	 (setf (aref h 0) 1d0)

	 (do ((j 1 (+ j 1)))
	     ((> j jmax) t)
	   (declare (type fixnum j))
	   (setf (aref s (1- j)) (dfloat-check (funcall trapzd func a b j)))
	   (when 
	       (>= j k)
	     (do ((i 1 (1+ i)))
		 ((> i k) t)
	       (declare (type fixnum i))
	       (setf (aref sj (1- i)) (aref s (- (+ i (- j km)) 2)))
	       (setf (aref hj (1- i)) (aref h (- (+ i (- j km)) 2))))


	     (multiple-value-setq (ss dss) (polint hj sj k 0d0)) 
	     (if (< (abs dss) (* eps (abs ss)))
		 (go end)))
	   (setf (aref s j) ss)		; (aref s (1- j)))
	   (setf (aref h j) (* 0.25d0 (aref h (1- j))))) 

	 (error " too many steps in qromb. ")
	 end
	 (return (the double-float ss))))

