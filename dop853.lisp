;;; Hairer-Nørsett-Wanner Dormand-Price Runge-Kutta solution of
;;; differential equations Translated from the FORTRAN code of Hairer,
;;; Nørsett, and Wanner, downloaded 6 April 2021 from
;;; http://www.unige.ch/~hairer/software.html.  See Ernst Hairer, Syvert
;;; Paul Nørsett, Gerhard Wanner, Solving Ordinary Differential
;;; Equations I. Nonstiff Problems.  Springer Series in
;;; Computational Mathematics, Vol. 8, Springer-Verlag 1987, Second revised
;;; edition 1993.  Translation by Ken Olum with the assistance of f2cl.

;;I have rewritten this code to match usual Lisp calling conventions:
;;  Work arrays eliminated in favor of storage allocated by code.
;;  Parameters in keyword arguments.  Defaults given by not specifying them.
;;  "Scalar" versus "vector" tolerances determined by datatype of
;;    tolerance arguments instead of special code.
;;
;;Features eliminated:
;; Parameter arrays to be passed back to user functions.  This can be done
;;   with closures or dynamically-bound variables.
;; Data passed to SOLOUT only for the purpose of passing back to CONTD8.
;; Special return from SOLOUT to abort the solution.  This can be handled
;;   by nonlocal exit
;; Special return from SOLOUT to indicate that the solution has altered,
;;   and so DP86CO should recompute the derivatives at the current point.
;;   I don't know what feature is provided here.  It seems that one could
;;   simply restart DOP853.  But it could be put back if needed.
;;
;;Some FORTRAN-style conventions are still present, though.

(in-package "CL-USER")

;; ----------------------------------------------------------
;;     Numerical solution of a system of first order
;;     ordinary differential equations  y'=f(x,y).
;;     this is an explicit runge-kutta method of order 8(5,3)
;;     due to dormand & prince (with stepsize control and
;;     dense output)
;;
;;     authors: e. hairer and g. wanner
;;              universite de geneve, dept. de mathematiques
;;              ch-1211 geneve 24, switzerland
;;              e-mail:  ernst.hairer@unige.ch
;;                       gerhard.wanner@unige.ch
;;
;;     this code is described in:
;;         e. hairer, s.p. norsett and g. wanner, solving ordinary
;;         differential equations i. nonstiff problems. 2nd edition.
;;         springer series in computational mathematics,
;;         springer-verlag (1993)
;;
;;     Translated from FORTRAN version of october 11, 2009
;;      (new option iout=3 for sparse dense output)
;;
;;     input parameters
;;     ----------------
;;     n           Dimension of the system
;;
;;     fcn         Function to compute computing f(x,y), called
;;                    with arguments (n x y f)
;;
;;     x           Initial x-value
;;
;;     y           Initial array of y values.  It will be modified as
;;                 the solution evolves.
;;
;;     xend        Final x-value (xend-x may be positive or negative)
;;
;;     solout      Function given the numerical solution during integration.
;;                 Unless iout = 0, it is called once with the initial
;;                 conditions and once for each successful step
;;                 with arguments (nr xold x y ifun):
;;                 y is the solution at the nr-th
;;                    grid-point x (thereby the initial value is
;;                    at nr=0)
;;                 xold is the preceeding grid-point.
;;                 ifun is a function to call with arguments (i s)
;;                 to get an approximation to the i-th
;;                 component of the solution at the point s.  The value
;;                 s should lie in the interval [xold,x].
;;                 If solout returns :abort, dop853 will exit.  This allows
;;                 you to get the predicted next step size.
;;                 Otherwise, if iout=3, the return value from solout is a
;;                 value xout giving the next x value for which we want to
;;                 perform dense output.
;;
;;     iout        switch for calling function solout:
;;                    iout=0: function is never called
;;                    iout=1: function is called after every successful step
;;                    iout=2: as for iout=1, but structures are set up
;;                            to allow ifun to be called for dense output
;;                    iout=3: dense output is performed in steps defined
;;                            by the user (see return from solout above).
;;
;;-----------------------------------------------------------------------
;;
;;     sophisticated setting of parameters
;;     -----------------------------------
;;              several keyword arguments to dop853 allow
;;              to adapt the code to the problem and to the needs of
;;              the user.
;;
;;    uround    the rounding unit, default double-float-epsilon
;;
;;    safe      the safety factor in step size prediction,
;;              default 0.9d0.
;;
;;    fac1, fac2  parameters for step size selection
;;              the new step size is chosen subject to the restriction
;;                 fac1 <= hnew/hold <= fac2
;;              default values: fac1=0.333d0, fac2=6.d0
;;
;;    beta      the beta for stabilized step size control
;;              (see section iv.2). positive values of beta ( <= 0.04 )
;;              make the step size control more stable.
;;              default 0.0d0.
;;
;;    hmax      maximal step size, default xend-x.
;;
;;    h         initial step size.  If not given an initial guess
;;              is computed with help of the function hinit
;;
;;    nmax      Maximal number of allowed steps.  Default 100000.
;;
;;    nstiff    test for stiffness is activated after step number
;;              j*nstiff (j integer).  If negative, the stiffness test is
;;              never activated; default value is 1000
;;
;;    nrdens    number of components, for which dense output
;;              is required.
;;
;;    icomp     array of component numbers for dense output if
;;              for   0 < nrdens < n.  For nrdens=n this is supplied
;;              aby the code.
;;
;;

;;----------------------------------------------------------------------
;;
;;     output parameters
;;     -----------------
;;
;;     y           numerical solution at final point
;;
;;     h           predicted step size of the last accepted step
;;
;;   nfcn          number of function evaluations
;;   nstep         number of computed steps
;;   naccpt        number of accepted steps
;;   nrejct        number of rejected steps (due to error test),
;;                 (rejections in the first step are not counted)
;;
;; You should expect 1 call for the initial conditions,
;; plus 1 if you didn't supply an initial stepsize,
;; 12 for each successful step (+ 3 more if you use dense output)
;; and 11 for each unsuccessful step
(defun dop853 (n fcn x y xend rtol atol solout iout
		 &key (uround double-float-epsilon)
		 (safe 0.9)
		 (fac1 0.333)
		 (fac2 6.0)
		 (beta 0.0)
		 (hmax (- xend x))
		 h			;NIL or 0.0 means guess using hinit
		 (nmax 100000)
		 (nstiff 1000)
		 nrdens
		 icomp maximize-errors)
  (declare (type (array double-float (*)) y)
           (type double-float xend x)
	   (optimize debug))
  ;;Check parameters
  (assert (> nmax 0))
  (unless h
    (setq h 0.0))			;Default to guess using hinit
  (when (< nstiff 0)			;Negative means never activate
    (setf nstiff (+ nmax 10)))
  (unless nrdens			;Default nrdens
    (setq nrdens (if (>= iout 2) n 0)))	;All components if dense output
  (when (or (< nrdens 0) (> nrdens n))
    (error "Invalid nrdens ~S" nrdens))
  (unless icomp
    (cond ((= nrdens n)			;If array not supplied, create
	   (setf icomp (make-array n :element-type '(signed-byte 32)))
	   (dotimes (i nrdens)
	     (setf (aref icomp i) i)
	     ))
	  ((zerop nrdens))		;Nothing to do
	  (t (error "Must supply icomp if nrdens given < n"))))
  (assert (typep rtol '(or double-float vector)))
  (assert (typep atol '(or double-float vector)))
  (assert (> 1.0 uround 1e-35))
  (assert (> 1.0 safe 1e-4))
  (assert (<= beta 0.2))
  ;; -------- call to core integrator ------------
  (dp86co n fcn x y xend hmax h rtol atol solout iout
	  nmax uround nstiff safe beta fac1 fac2 icomp nrdens
	  :maximize-errors maximize-errors))

(defparameter dop853-c2 0.05260015195876773)
(defparameter dop853-c3 0.0789002279381516)
(defparameter dop853-c4 0.1183503419072274)
(defparameter dop853-c5 0.2816496580927726)
(defparameter dop853-c6 0.3333333333333333)
(defparameter dop853-c7 0.25)
(defparameter dop853-c8 0.3076923076923077)
(defparameter dop853-c9 0.6512820512820513)
(defparameter dop853-c10 0.6)
(defparameter dop853-c11 0.8571428571428571)
(defparameter dop853-c14 0.1)
(defparameter dop853-c15 0.2)
(defparameter dop853-c16 0.7777777777777778)
(defparameter dop853-b1 0.054293734116568765)
(defparameter dop853-b6 4.450312892752409)
(defparameter dop853-b7 1.8915178993145003)
(defparameter dop853-b8 (- 5.801203960010585))
(defparameter dop853-b9 0.3111643669578199)
(defparameter dop853-b10 (- 0.1521609496625161))
(defparameter dop853-b11 0.20136540080403034)
(defparameter dop853-b12 0.04471061572777259)
(defparameter dop853-bhh1 0.2440944881889764)
(defparameter dop853-bhh2 0.7338466882816118)
(defparameter dop853-bhh3 0.022058823529411766)
(defparameter dop853-er1 0.01312004499419488)
(defparameter dop853-er6 (- 1.2251564463762044))
(defparameter dop853-er7 (- 0.4957589496572502))
(defparameter dop853-er8 1.6643771824549864)
(defparameter dop853-er9 (- 0.35032884874997366))
(defparameter dop853-er10 0.3341791187130175)
(defparameter dop853-er11 0.08192320648511571)
(defparameter dop853-er12 (- 0.022355307863886294))
(defparameter dop853-a21 0.05260015195876773)
(defparameter dop853-a31 0.0197250569845379)
(defparameter dop853-a32 0.0591751709536137)
(defparameter dop853-a41 0.02958758547680685)
(defparameter dop853-a43 0.08876275643042054)
(defparameter dop853-a51 0.2413651341592667)
(defparameter dop853-a53 (- 0.8845494793282861))
(defparameter dop853-a54 0.924834003261792)
(defparameter dop853-a61 0.037037037037037035)
(defparameter dop853-a64 0.17082860872947386)
(defparameter dop853-a65 0.12546768756682242)
(defparameter dop853-a71 0.037109375)
(defparameter dop853-a74 0.17025221101954405)
(defparameter dop853-a75 0.06021653898045596)
(defparameter dop853-a76 (- 0.017578125))
(defparameter dop853-a81 0.03709200011850479)
(defparameter dop853-a84 0.17038392571223998)
(defparameter dop853-a85 0.10726203044637328)
(defparameter dop853-a86 (- 0.015319437748624402))
(defparameter dop853-a87 0.008273789163814023)
(defparameter dop853-a91 0.6241109587160757)
(defparameter dop853-a94 (- 3.3608926294469414))
(defparameter dop853-a95 (- 0.868219346841726))
(defparameter dop853-a96 27.59209969944671)
(defparameter dop853-a97 20.154067550477894)
(defparameter dop853-a98 (- 43.48988418106996))
(defparameter dop853-a101 0.47766253643826434)
(defparameter dop853-a104 (- 2.4881146199716677))
(defparameter dop853-a105 (- 0.590290826836843))
(defparameter dop853-a106 21.230051448181193)
(defparameter dop853-a107 15.279233632882423)
(defparameter dop853-a108 (- 33.28821096898486))
(defparameter dop853-a109 (- 0.020331201708508627))
(defparameter dop853-a111 (- 0.9371424300859873))
(defparameter dop853-a114 5.186372428844064)
(defparameter dop853-a115 1.0914373489967295)
(defparameter dop853-a116 (- 8.149787010746927))
(defparameter dop853-a117 (- 18.52006565999696))
(defparameter dop853-a118 22.739487099350505)
(defparameter dop853-a119 2.4936055526796523)
(defparameter dop853-a1110 (- 3.0467644718982196))
(defparameter dop853-a121 2.273310147516538)
(defparameter dop853-a124 (- 10.53449546673725))
(defparameter dop853-a125 (- 2.0008720582248625))
(defparameter dop853-a126 (- 17.9589318631188))
(defparameter dop853-a127 27.94888452941996)
(defparameter dop853-a128 (- 2.8589982771350235))
(defparameter dop853-a129 (- 8.87285693353063))
(defparameter dop853-a1210 12.360567175794303)
(defparameter dop853-a1211 0.6433927460157636)
(defparameter dop853-a141 0.056167502283047954)
(defparameter dop853-a147 0.25350021021662483)
(defparameter dop853-a148 (- 0.2462390374708025))
(defparameter dop853-a149 (- 0.12419142326381637))
(defparameter dop853-a1410 0.15329179827876568)
(defparameter dop853-a1411 0.00820105229563469)
(defparameter dop853-a1412 0.007567897660545699)
(defparameter dop853-a1413 (- 0.008298))
(defparameter dop853-a151 0.03183464816350214)
(defparameter dop853-a156 0.028300909672366776)
(defparameter dop853-a157 0.053541988307438566)
(defparameter dop853-a158 (- 0.05492374857139099))
(defparameter dop853-a1511 (- 1.0834732869724932e-4))
(defparameter dop853-a1512 3.825710908356584e-4)
(defparameter dop853-a1513 (- 3.4046500868740456e-4))
(defparameter dop853-a1514 0.1413124436746325)
(defparameter dop853-a161 (- 0.42889630158379194))
(defparameter dop853-a166 (- 4.697621415361164))
(defparameter dop853-a167 7.683421196062599)
(defparameter dop853-a168 4.06898981839711)
(defparameter dop853-a169 0.3567271874552811)
(defparameter dop853-a1613 (- 0.0013990241651590145))
(defparameter dop853-a1614 2.9475147891527724)
(defparameter dop853-a1615 (- 9.15095847217987))
(defparameter dop853-d41 (- 8.428938276109013))
(defparameter dop853-d46 0.5667149535193777)
(defparameter dop853-d47 (- 3.0689499459498917))
(defparameter dop853-d48 2.38466765651207)
(defparameter dop853-d49 2.117034582445028)
(defparameter dop853-d410 (- 0.871391583777973))
(defparameter dop853-d411 2.2404374302607883)
(defparameter dop853-d412 0.6315787787694688)
(defparameter dop853-d413 (- 0.08899033645133331))
(defparameter dop853-d414 18.148505520854727)
(defparameter dop853-d415 (- 9.194632392478356))
(defparameter dop853-d416 (- 4.436036387594894))
(defparameter dop853-d51 10.427508642579134)
(defparameter dop853-d56 242.28349177525817)
(defparameter dop853-d57 165.20045171727028)
(defparameter dop853-d58 (- 374.5467547226902))
(defparameter dop853-d59 (- 22.113666853125306))
(defparameter dop853-d510 7.733432668472264)
(defparameter dop853-d511 (- 30.674084731089398))
(defparameter dop853-d512 (- 9.332130526430229))
(defparameter dop853-d513 15.697238121770845)
(defparameter dop853-d514 (- 31.139403219565178))
(defparameter dop853-d515 (- 9.35292435884448))
(defparameter dop853-d516 35.81684148639408)
(defparameter dop853-d61 19.985053242002433)
(defparameter dop853-d66 (- 387.0373087493518))
(defparameter dop853-d67 (- 189.17813819516758))
(defparameter dop853-d68 527.8081592054236)
(defparameter dop853-d69 (- 11.57390253995963))
(defparameter dop853-d610 6.8812326946963)
(defparameter dop853-d611 (- 1.0006050966910838))
(defparameter dop853-d612 0.7777137798053443)
(defparameter dop853-d613 (- 2.778205752353508))
(defparameter dop853-d614 (- 60.19669523126412))
(defparameter dop853-d615 84.32040550667716)
(defparameter dop853-d616 11.99229113618279)
(defparameter dop853-d71 (- 25.69393346270375))
(defparameter dop853-d76 (- 154.18974869023643))
(defparameter dop853-d77 (- 231.5293791760455))
(defparameter dop853-d78 357.6391179106141)
(defparameter dop853-d79 93.40532418362432)
(defparameter dop853-d710 (- 37.45832313645163))
(defparameter dop853-d711 104.0996495089623)
(defparameter dop853-d712 29.8402934266605)
(defparameter dop853-d713 (- 43.53345659001114))
(defparameter dop853-d714 96.32455395918828)
(defparameter dop853-d715 (- 39.17726167561544))
(defparameter dop853-d716 (- 149.72683625798564))
(declaim
 (double-float 
  dop853-c2 dop853-c3 dop853-c4 dop853-c5 dop853-c6 dop853-c7
  dop853-c8 dop853-c9 dop853-c10 dop853-c11 dop853-c14 dop853-c15
  dop853-c16 dop853-b1 dop853-b6 dop853-b7 dop853-b8 dop853-b9
  dop853-b10 dop853-b11 dop853-b12 dop853-bhh1 dop853-bhh2 dop853-bhh3
  dop853-er1 dop853-er6 dop853-er7 dop853-er8 dop853-er9 dop853-er10
  dop853-er11 dop853-er12 dop853-a21 dop853-a31 dop853-a32 dop853-a41
  dop853-a43 dop853-a51 dop853-a53 dop853-a54 dop853-a61 dop853-a64
  dop853-a65 dop853-a71 dop853-a74 dop853-a75 dop853-a76 dop853-a81
  dop853-a84 dop853-a85 dop853-a86 dop853-a87 dop853-a91 dop853-a94
  dop853-a95 dop853-a96 dop853-a97 dop853-a98 dop853-a101 dop853-a104
  dop853-a105 dop853-a106 dop853-a107 dop853-a108 dop853-a109
  dop853-a111 dop853-a114 dop853-a115 dop853-a116 dop853-a117
  dop853-a118 dop853-a119 dop853-a1110 dop853-a121 dop853-a124
  dop853-a125 dop853-a126 dop853-a127 dop853-a128 dop853-a129
  dop853-a1210 dop853-a1211 dop853-a141 dop853-a147 dop853-a148
  dop853-a149 dop853-a1410 dop853-a1411 dop853-a1412 dop853-a1413
  dop853-a151 dop853-a156 dop853-a157 dop853-a158 dop853-a1511
  dop853-a1512 dop853-a1513 dop853-a1514 dop853-a161 dop853-a166
  dop853-a167 dop853-a168 dop853-a169 dop853-a1613 dop853-a1614
  dop853-a1615 dop853-d41 dop853-d46 dop853-d47 dop853-d48 dop853-d49
  dop853-d410 dop853-d411 dop853-d412 dop853-d413 dop853-d414
  dop853-d415 dop853-d416 dop853-d51 dop853-d56 dop853-d57 dop853-d58
  dop853-d59 dop853-d510 dop853-d511 dop853-d512 dop853-d513 dop853-d514
  dop853-d515 dop853-d516 dop853-d61 dop853-d66 dop853-d67 dop853-d68
  dop853-d69 dop853-d610 dop853-d611 dop853-d612 dop853-d613 dop853-d614
  dop853-d615 dop853-d616 dop853-d71 dop853-d76 dop853-d77 dop853-d78
  dop853-d79 dop853-d710 dop853-d711 dop853-d712 dop853-d713 dop853-d714
  dop853-d715 dop853-d716))

;;atol or rtol is either an array of tolerances
;;or a single tolerance to be used for all
(defun get-tolerance (tol i)
  (if (typep tol 'double-float) tol
    (aref tol i)))

(defun dp86co (n fcn x y xend hmax h rtol atol solout iout nmax
		 uround nstiff safe beta fac1 fac2 icomp nrd 
		 &key maximize-errors) ;Maximize error instead of taking RMS average.  Added by Ken Olum
  (let ((nfcn 0)
	(nstep 0)
	(naccpt 0)
	(nrejct 0)
	(k1 (make-array n :element-type 'double-float)) ;Make our own working storage
	(k2 (make-array n :element-type 'double-float))
	(k3 (make-array n :element-type 'double-float))
	(k4 (make-array n :element-type 'double-float))
	(k5 (make-array n :element-type 'double-float))
	(k6 (make-array n :element-type 'double-float))
	(k7 (make-array n :element-type 'double-float))
	(k8 (make-array n :element-type 'double-float))
	(k9 (make-array n :element-type 'double-float))
	(k10 (make-array n :element-type 'double-float))
	(y1 (make-array n :element-type 'double-float))
	(cont (make-array (* 8 n) :element-type 'double-float)))
    (declare (type double-float fac2 fac1 beta safe uround h hmax xend x)
	     (type (signed-byte 32) nrejct naccpt nstep nfcn nrd nstiff
		   nmax iout n))
    (let ((xold 0.0)
	  (bspl 0.0) (ydiff 0.0) (nonsti 0) (stden 0.0)
	  (stnum 0.0) (hnew 0.0) (fac 0.0) (fac11 0.0) (deno 0.0)
	  (erri 0.0) (sk 0.0) (err2 0.0) (err 0.0) (xph 0.0)
	  (xout 0.0) (iord 0) (iasti 0) (hlamb 0.0)
	  (posneg 0.0) (facc2 0.0) (facc1 0.0)
	  (expo1 0.0) (facold 0.0) (event nil) (last nil) (reject nil))
      (declare (type boolean event last reject)
	       (type (signed-byte 32) iasti iord nonsti)
	       (type double-float xold facold expo1 facc1 facc2 posneg
		     hlamb xph err err2 sk erri deno fac11 fac hnew
		     stnum stden ydiff bspl))
      ;; ----------------------------------------------------------
      ;;     core integrator for dop853
      ;; *** *** *** *** *** *** ***
      ;;  initialisations
      ;; *** *** *** *** *** *** ***
      (setf facold 1.0e-4)
      (setf expo1 (- (/ 1.0 8.0) (* beta 0.2)))
      (setf facc1 (/ 1.0 fac1))
      (setf facc2 (/ 1.0 fac2))
      (setf posneg (signum (- xend x)))
      ;; --- initial preparations
      (setf last nil)
      (setf hlamb 0.0)
      (setf iasti 0)
      (funcall fcn n x y k1)		;Initial derivative evaluation
      (incf nfcn)
      (setf hmax (abs hmax))
      (setf iord 8)
      (when (= h 0.0)
	(setf h (hinit n fcn x y posneg k1 iord hmax
		       atol rtol))
	;;Count one call if we called hinit.  This is wrong
	;;in the original code
	(incf nfcn))
      (setf reject nil)
      (setf xold x)
      (unless (= iout 0)
	;;Initial call to SOLOUT with no solution yet.  The only purpose
	;;of this is to allow setting the initial XOUT
	(setf xout (funcall solout 0 xold x y n nil)))
      (loop				;Loop until we reach xend
       ;; --- basic integration step
       (when (> nstep nmax)
	 (error "More than nmax = ~S steps are needed" nmax))
       (when (<= (* 0.1 (abs h)) (* (abs x) uround))
	 (error "Step size too small, h=~S at x=~S" h x))
       (when (> (* (+ x (* 1.01 h) (- xend)) posneg) 0.0)
	 (setf h (- xend x))
	 (setf last t))
       (incf nstep)
       (format t "Trying DOP853 step size ~S~%" h)
       ;;Eleven derivative evaluations.
       ;;I removed here a piece of code which recomputes K1 in case SOLOUT
       ;;returned a code indicating that the solution had been altered
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h dop853-a21 (aref k1 i)))))
       (funcall fcn n (+ x (* dop853-c2 h)) y1 k2)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a31 (aref k1 i))
			(* dop853-a32 (aref k2 i)))))))
       (funcall fcn n (+ x (* dop853-c3 h)) y1 k3)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a41 (aref k1 i))
			(* dop853-a43 (aref k3 i)))))))
       (funcall fcn n (+ x (* dop853-c4 h)) y1 k4)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a51 (aref k1 i))
			(* dop853-a53 (aref k3 i))
			(* dop853-a54 (aref k4 i)))))))
       (funcall fcn n (+ x (* dop853-c5 h)) y1 k5)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a61 (aref k1 i))
			(* dop853-a64 (aref k4 i))
			(* dop853-a65 (aref k5 i)))))))
       (funcall fcn n (+ x (* dop853-c6 h)) y1 k6)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a71 (aref k1 i))
			(* dop853-a74 (aref k4 i))
			(* dop853-a75 (aref k5 i))
			(* dop853-a76 (aref k6 i)))))))
       (funcall fcn n (+ x (* dop853-c7 h)) y1 k7)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a81 (aref k1 i))
			(* dop853-a84 (aref k4 i))
			(* dop853-a85 (aref k5 i))
			(* dop853-a86 (aref k6 i))
			(* dop853-a87 (aref k7 i)))))))
       (funcall fcn n (+ x (* dop853-c8 h)) y1 k8)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a91 (aref k1 i))
			(* dop853-a94 (aref k4 i))
			(* dop853-a95 (aref k5 i))
			(* dop853-a96 (aref k6 i))
			(* dop853-a97 (aref k7 i))
			(* dop853-a98 (aref k8 i)))))))
       (funcall fcn n (+ x (* dop853-c9 h)) y1 k9)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a101 (aref k1 i))
			(* dop853-a104 (aref k4 i))
			(* dop853-a105 (aref k5 i))
			(* dop853-a106 (aref k6 i))
			(* dop853-a107 (aref k7 i))
			(* dop853-a108 (aref k8 i))
			(* dop853-a109 (aref k9 i)))))))
       (funcall fcn n (+ x (* dop853-c10 h)) y1 k10)
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a111 (aref k1 i))
			(* dop853-a114 (aref k4 i))
			(* dop853-a115 (aref k5 i))
			(* dop853-a116 (aref k6 i))
			(* dop853-a117 (aref k7 i))
			(* dop853-a118 (aref k8 i))
			(* dop853-a119 (aref k9 i))
			(* dop853-a1110 (aref k10 i)))))))
       (funcall fcn n (+ x (* dop853-c11 h)) y1 k2)
       (setf xph (+ x h))
       (dotimes (i n) 
	 (setf (aref y1 i)
	       (+ (aref y i)
		  (* h
		     (+ (* dop853-a121 (aref k1 i))
			(* dop853-a124 (aref k4 i))
			(* dop853-a125 (aref k5 i))
			(* dop853-a126 (aref k6 i))
			(* dop853-a127 (aref k7 i))
			(* dop853-a128 (aref k8 i))
			(* dop853-a129 (aref k9 i))
			(* dop853-a1210 (aref k10 i))
			(* dop853-a1211 (aref k2 i)))))))
       (funcall fcn n xph y1 k3)
       (incf nfcn 11)
       (dotimes (i n) 
	 (setf (aref k4 i)
	       (+ (* dop853-b1 (aref k1 i))
		  (* dop853-b6 (aref k6 i))
		  (* dop853-b7 (aref k7 i))
		  (* dop853-b8 (aref k8 i))
		  (* dop853-b9 (aref k9 i))
		  (* dop853-b10 (aref k10 i))
		  (* dop853-b11 (aref k2 i))
		  (* dop853-b12 (aref k3 i))))
	 (setf (aref k5 i)
	       (+ (aref y i)
		  (* h (aref k4 i)))))
       ;; --- error estimation
       (setf err 0.0)
       (setf err2 0.0)
       (cond
	(maximize-errors	  ;Maximum error rather than RMS average
	 (let (err1)
	   (dotimes (i n)
	     (setf sk
		   (+ (get-tolerance atol i)
		      (* (get-tolerance rtol i)
			 (max
			  (abs (aref y i))
			  (abs (aref k5 i))))))
	     (setf erri
		   (- (aref k4 i)
		      (* dop853-bhh1 (aref k1 i))
		      (* dop853-bhh2 (aref k9 i))
		      (* dop853-bhh3 (aref k3 i))))
;;	     (format t "err2 = ~S, " (/ erri sk))
	     (setf err2 (expt (/ erri sk) 2)) ;err2 from this component only
	     (setf erri
		   (+ (* dop853-er1 (aref k1 i))
		      (* dop853-er6 (aref k6 i))
		      (* dop853-er7 (aref k7 i))
		      (* dop853-er8 (aref k8 i))
		      (* dop853-er9 (aref k9 i))
		      (* dop853-er10 (aref k10 i))
		      (* dop853-er11 (aref k2 i))
		      (* dop853-er12 (aref k3 i))))
;;	     (format t "err1 = ~S, " (/ erri sk))
	     (setf err1 (expt (/ erri sk) 2)) ;err1 from this component only
	     (setf deno (+ err1 (* 0.01 err2)))	;Regular calculation for just this component
	     (when (zerop deno)	       ;Was <=, but that caused warnings
	       (setf deno 1.0))
	     (let ((this (* (abs h) err1 (sqrt (/ 1.0 deno)))))
;;	       (format t "this = ~S~%" this)
	       (setf err (max err this)))))) ;Maximize over all components
	(t				;Use RMS error
	 (dotimes (i n) 
	   (setf sk
		 (+ (get-tolerance atol i)
		    (* (get-tolerance rtol i)
		       (max
			(abs (aref y i))
			(abs (aref k5 i))))))
	   (setf erri
		 (- (aref k4 i)
		    (* dop853-bhh1 (aref k1 i))
		    (* dop853-bhh2 (aref k9 i))
		    (* dop853-bhh3 (aref k3 i))))
;;	   (format t "err2 = ~S, " (/ erri sk))
	   (setf err2 (+ err2 (expt (/ erri sk) 2)))
	   (setf erri
		 (+ (* dop853-er1 (aref k1 i))
		    (* dop853-er6 (aref k6 i))
		    (* dop853-er7 (aref k7 i))
		    (* dop853-er8 (aref k8 i))
		    (* dop853-er9 (aref k9 i))
		    (* dop853-er10 (aref k10 i))
		    (* dop853-er11 (aref k2 i))
		    (* dop853-er12 (aref k3 i))))
;;	   (format t "err1 = ~S~%" (/ erri sk))
	   (setf err (+ err (expt (/ erri sk) 2))))
;;	 (format t "sum-squared err2 = ~S, err = ~S; " err2 err)
	 (setf deno (+ err (* 0.01 err2)))
	 (when (<= deno 0.0)
	   (setf deno 1.0))
	 (setf err (* (abs h) err (sqrt (/ 1.0 (* n deno)))))
;;	 (format t "RMS error ~S~%" err))
	 ))
       (format t "Error estimate ~S~%" err)
       ;; --- computation of hnew
       (setf fac11 (expt err expo1))
       ;; --- lund-stabilization
       (setf fac (/ fac11 (expt facold beta)))
       ;; --- we require  fac1 <= hnew/h <= fac2
       (setf fac (max facc2 (min facc1 (/ fac safe))))
       (setf hnew (/ h fac))
       (cond
	((<= err 1.0) ; --- step is accepted.
	 (setf facold (max err 1.0e-4))
	 (incf naccpt)
	 (funcall fcn n xph k5 k4)	;Derivatives at final point
	 (incf nfcn)
	 ;; ------- stiffness detection
	 (when (or (= (mod naccpt nstiff) 0) (> iasti 0)) (setf stnum 0.0)
	   (setf stden 0.0)
	   (dotimes (i n) 
	     (setf stnum
		   (+ stnum
		      (expt
		       (- (aref k4 i) (aref k3 i))
		       2)))
	     (setf stden
		   (+ stden
		      (expt
		       (- (aref k5 i) (aref y1 i))
		       2))))
	   (when (> stden 0.0)
	     (setf hlamb (* (abs h) (sqrt (/ stnum stden)))))
	   (cond
	    ((> hlamb 6.1) (setf nonsti 0)
	     (incf iasti)
	     (when (= iasti 15)
	       (error "The problem seems to become stiff at X = ~S" x)))
	    (t (incf nonsti)
	       (when (= nonsti 6)
		 (setf iasti 0)))))
	 ;; ------- final preparation for dense output
	 (setf event (and (= iout 3) (<= xout xph)))
	 (when (or (= iout 2) event)
	   ;; ----    save the first function evaluations
	   (dotimes (j nrd)
	     (let ((i (aref icomp j)))
	       (setf (aref cont j) (aref y i))
	       (setf ydiff
		     (- (aref k5 i) (aref y i)))
	       (setf (aref cont (+ j nrd))
		     ydiff)
	       (setf bspl
		     (- (* h (aref k1 i)) ydiff))
	       (setf (aref cont (+ j (* nrd 2)))
		     bspl)
	       (setf (aref cont (+ j (* nrd 3)))
		     (- ydiff
			(* h (aref k4 i))
			bspl))
	       (setf (aref cont (+ j (* nrd 4)))
		     (+ (* dop853-d41 (aref k1 i))
			(* dop853-d46 (aref k6 i))
			(* dop853-d47 (aref k7 i))
			(* dop853-d48 (aref k8 i))
			(* dop853-d49 (aref k9 i))
			(* dop853-d410 (aref k10 i))
			(* dop853-d411 (aref k2 i))
			(* dop853-d412 (aref k3 i))))
	       (setf (aref cont (+ j (* nrd 5)))
		     (+ (* dop853-d51 (aref k1 i))
			(* dop853-d56 (aref k6 i))
			(* dop853-d57 (aref k7 i))
			(* dop853-d58 (aref k8 i))
			(* dop853-d59 (aref k9 i))
			(* dop853-d510 (aref k10 i))
			(* dop853-d511 (aref k2 i))
			(* dop853-d512 (aref k3 i))))
	       (setf (aref cont (+ j (* nrd 6)))
		     (+ (* dop853-d61 (aref k1 i))
			(* dop853-d66 (aref k6 i))
			(* dop853-d67 (aref k7 i))
			(* dop853-d68 (aref k8 i))
			(* dop853-d69 (aref k9 i))
			(* dop853-d610 (aref k10 i))
			(* dop853-d611 (aref k2 i))
			(* dop853-d612 (aref k3 i))))
	       (setf (aref cont (+ j (* nrd 7)))
		     (+ (* dop853-d71 (aref k1 i))
			(* dop853-d76 (aref k6 i))
			(* dop853-d77 (aref k7 i))
			(* dop853-d78 (aref k8 i))
			(* dop853-d79 (aref k9 i))
			(* dop853-d710 (aref k10 i))
			(* dop853-d711 (aref k2 i))
			(* dop853-d712 (aref k3 i))))))
	   ;;Three function evaluations to set up dense output
	   (dotimes (i n) 
	     (setf (aref y1 i)
		   (+ (aref y i)
		      (* h
			 (+ (* dop853-a141 (aref k1 i))
			    (* dop853-a147 (aref k7 i))
			    (* dop853-a148 (aref k8 i))
			    (* dop853-a149 (aref k9 i))
			    (* dop853-a1410 (aref k10 i))
			    (* dop853-a1411 (aref k2 i))
			    (* dop853-a1412 (aref k3 i))
			    (* dop853-a1413 (aref k4 i)))))))
	   (funcall fcn n (+ x (* dop853-c14 h)) y1 k10)
	   (dotimes (i n) 
	     (setf (aref y1 i)
		   (+ (aref y i)
		      (* h
			 (+ (* dop853-a151 (aref k1 i))
			    (* dop853-a156 (aref k6 i))
			    (* dop853-a157 (aref k7 i))
			    (* dop853-a158 (aref k8 i))
			    (* dop853-a1511 (aref k2 i))
			    (* dop853-a1512 (aref k3 i))
			    (* dop853-a1513 (aref k4 i))
			    (* dop853-a1514 (aref k10 i)))))))
	   (funcall fcn n (+ x (* dop853-c15 h)) y1 k2)
	   (dotimes (i n) 
	     (setf (aref y1 i)
		   (+ (aref y i)
		      (* h
			 (+ (* dop853-a161 (aref k1 i))
			    (* dop853-a166 (aref k6 i))
			    (* dop853-a167 (aref k7 i))
			    (* dop853-a168 (aref k8 i))
			    (* dop853-a169 (aref k9 i))
			    (* dop853-a1613 (aref k4 i))
			    (* dop853-a1614 (aref k10 i))
			    (* dop853-a1615 (aref k2 i)))))))
	   (funcall fcn n (+ x (* dop853-c16 h)) y1 k3)
	   (incf nfcn 3)
	   ;; ---     final preparation
	   (dotimes (j nrd)
	     (let ((i (aref icomp j)))
	       (setf (aref cont (+ j (* nrd 4)))
		     (* h
			(+ (aref cont (+ j (* nrd 4)))
			   (* dop853-d413 (aref k4 i))
			   (* dop853-d414 (aref k10 i))
			   (* dop853-d415 (aref k2 i))
			   (* dop853-d416 (aref k3 i)))))
	       (setf (aref cont (+ j (* nrd 5)))
		     (* h
			(+ (aref cont (+ j (* nrd 5)))
			   (* dop853-d513 (aref k4 i))
			   (* dop853-d514 (aref k10 i))
			   (* dop853-d515 (aref k2 i))
			   (* dop853-d516 (aref k3 i)))))
	       (setf (aref cont (+ j (* nrd 6)))
		     (* h
			(+ (aref cont (+ j (* nrd 6)))
			   (* dop853-d613 (aref k4 i))
			   (* dop853-d614 (aref k10 i))
			   (* dop853-d615 (aref k2 i))
			   (* dop853-d616 (aref k3 i)))))
	       (setf (aref cont (+ j (* nrd 7)))
		     (* h
			(+ (aref cont (+ j (* nrd 7)))
			   (* dop853-d713 (aref k4 i))
			   (* dop853-d714 (aref k10 i))
			   (* dop853-d715 (aref k2 i))
			   (* dop853-d716 (aref k3 i))))))))
	 (dotimes (i n) 
	   (setf (aref k1 i) (aref k4 i))
	   (setf (aref y i) (aref k5 i)))
	 (setf xold x)
	 (setf x xph)
	 (when (or (= iout 1) (= iout 2) event)
	   (setf xout
		 (funcall solout naccpt xold x y n
			  ;;Function to find component ii of Y at point s
			  #'(lambda (ii s)
			      (contd8 ii xold s h cont icomp nrd)))))
	 (when (or last			;------- normal exit
		   (eq xout :abort))	;exit requested
	   (return (values hnew nfcn nstep naccpt nrejct)))
	 (when (> (abs hnew) hmax)
	   (setf hnew (* posneg hmax)))
	 (when reject
	   (setf hnew (* posneg (min (abs hnew) (abs h)))))
	 (setf reject nil))
	(t				;--- step is rejected
	 (setf hnew (/ h (min facc1 (/ fac11 safe))))
	 (setf reject t)
	 (when (>= naccpt 1)		
	   ;;Only count rejected steps that are not in the initial tries
	   (incf nrejct))
	 (setf last nil)))
       (setf h hnew))			;Loop for next segment
      )))

;;Guess the initial stepsize.  F0 gives the derivatives at the
;;at the initial position.
;;This code does not know about maximize-errors
(defun hinit (n fcn x y posneg f0 iord hmax atol rtol)
  (declare (type (array double-float (*)) f0 y)
           (type double-float hmax posneg x)
           (type (signed-byte 32) iord n))
  (let ((h1 0.0) (der12 0.0) (der2 0.0) (h 0.0) (sk 0.0)
	(dny 0.0) (dnf 0.0)
	(f1 (make-array n :element-type 'double-float))
	(y1 (make-array n :element-type 'double-float)))
    (declare (type double-float dnf dny sk h der2 der12 h1))
    ;; ----------------------------------------------------------
    ;; ----  computation of an initial step size guess
    ;; ----------------------------------------------------------
    ;; ---- compute a first guess for explicit euler as
    ;; ----   h = 0.01 * norm (y0) / norm (f0)
    ;; ---- the increment for explicit euler is small
    ;; ---- compared to the solution
    (setf dnf 0.0)
    (setf dny 0.0)
    (dotimes (i n) 
      (setf sk (+ (get-tolerance atol i)
		  (* (get-tolerance rtol i)
		     (abs (aref y i)))))
      (setf dnf (+ dnf (expt (/ (aref f0 i) sk) 2)))
      (setf dny (+ dny (expt (/ (aref y i) sk) 2))))
    (cond ((or (<= dnf 1.0e-10) (<= dny 1.0e-10))
	   (setf h 1.0e-6))
	  (t (setf h (* (sqrt (/ dny dnf)) 0.01))))
    (setf h (min h hmax))
    (setf h (* h posneg))  ;Now h is a step in the direction we're going
    ;; ---- perform an explicit euler step
    (dotimes (i n) 
      (setf (aref y1 i)
	    (+ (aref y i)
	       (* h (aref f0 i)))))
    (funcall fcn n (+ x h) y1 f1)
    ;; ---- estimate the second derivative of the solution
    (setf der2 0.0)
    (dotimes (i n) 
      (setf sk
	    (+ (get-tolerance atol i)
	       (* (get-tolerance rtol i)
		  (abs (aref y i)))))
      (setf der2
	    (+ der2
	       (expt
		(/
		 (- (aref f1 i) (aref f0 i))
		 sk)
		2))))
    (setf der2 (/ (sqrt der2) h))
    ;; ---- step size is computed such that
    ;; ----  h**iord * max ( norm (f0), norm (der2)) = 0.01
    (setf der12 (max (abs der2) (sqrt dnf)))
    (cond ((<= der12 1.0e-15) (setf h1 (max 1.0e-6 (* (abs h) 0.001))))
	  (t (setf h1 (expt (/ 0.01 der12) (/ 1.0 iord)))))
    (setf h (min (* 100 (abs h)) h1 hmax)) ;This is the step we will use
    (* h posneg)))			   ;Return with right sign

;;SOLOUT can call this with call this function via the closure passed to it
;;to compute component ii of the solution at point s.
(defun contd8 (ii xold x h con icomp nd)
  (declare (type (array (signed-byte 32) (*)) icomp)
           (type (array double-float (*)) con)
           (type double-float x)
           (type (signed-byte 32) nd ii))
  (let ((conpar 0.0) (s1 0.0) (s 0.0) (i nil))
    (declare (type double-float s s1 conpar))
    ;; ----------------------------------------------------------
    ;;     this function can be used for continuous output in connection
    ;;     with the output-subroutine for dop853. it provides an
    ;;     approximation to the ii-th component of the solution at x.
    ;; ----------------------------------------------------------
    ;; ----- compute place of ii-th component
    (dotimes (j nd) 
      (when (= (aref icomp j) ii)
	(setf i j)))
    (unless i
      (error "No dense output available for comp. ~D" ii))
    (setf s (/ (- x xold) h))
    (setf s1 (- 1.0 s))
    (setf conpar (+ (aref con (+ i (* nd 4)))
		    (* s (+ (aref con (+ i (* nd 5)))
			    (* s1 (+ (aref con (+ i (* nd 6)))
				     (* s (aref con (+ i (* nd 7))))))))))
    (+ (aref con ii)			;Compute and return value
       (* s (+ (aref con (+ i nd))
	       (* s1 (+ (aref con (+ i (* nd 2)))
			(* s (+ (aref con (+ i (* nd 3)))
				(* s1 conpar))))))))))

;;;Simplified interface by Ken Olum

;;No parameters.  Usually in lisp you do this with special variables or closures.

;;FCN is the function to evaluate the derivatives
;;X is the starting point for the independent variable and XEND the end
;;Y is the starting (vector) position.  It will be modified.
;;RTOL and ATOL are relative and absolute tolerances, either as numbers or vectors
;;If DENSE-OUTPUT-FUNCTION is set, we request dense output and this function will be called after each
;;successful solution with arguments (XOLD X FN), giving are the range of the solution and a function
;;FN to call with II and X to get the IIth component of the solution at that point.
;;
;;FUNCTION will be called with (x y f).  X is the independent variable,
;;Y the current (vector) position, and F the array into which you should put the derivatives.
;;F may have extra elements, in which case they should be ignored
(defun do-dop853 (function xstart y xend rtol atol 
			   &optional dense-output-function maximize-errors
			   start-step-size)
  (let* ((n (length y)))
    (multiple-value-bind (hnew nfcn nstep naccpt nrejct)
	(dop853 (length y)
		#'(lambda (n x y f)	;Function called by DOP853
		    ;;Y may be a displaced array into the workspace with the wrong length
		    (setq y (make-array n :element-type 'double-float :displaced-to y))
		    (funcall function x y f)) ;Call provided function to fill in derivatives
		xstart y xend rtol atol
		#'(lambda (nr xold x y n interpolate-function) ;Called to record solution
		    (declare (ignore y n))
		    (unless (zerop nr)	;First call is to set up sparse output.  Don't compute anything.  Return NIL
		      (funcall dense-output-function xold x interpolate-function)))
		(if dense-output-function 2 0) ;dense output if requested
		:nrdens (if dense-output-function n 0)
		:h start-step-size
		:maximize-errors maximize-errors)
      (format t "DOP853 succeeded.  ~D steps (~D accepted, ~D rejected) next step ~S. ~D derivative calls~%"
	      nstep naccpt nrejct hnew nfcn)
      (values hnew nfcn nstep naccpt nrejct))))

;;Call DO-DOP853 for dense output.  Positions will be logged in steps of DENSE-STEP starting at time DENSE-START
;;by calling OUTPUT-FUNCTION with (step x y).  If the function returns :ABORT, DOP853 will exit.
;;If XEND is credibly close to a step, we use it.
;;"sparse dense output" is not implemented.  We could easily not do the computations needed for dense output
;;if dense-start is significantly after xstart but we have not.
;;Returns the last stepsize used that didn't end at XEND
(defun do-dop853-dense (function xstart y xend rtol atol
				 &key (dense-start xstart) dense-step
				 output-function maximize-errors start-step-size)
  (declare (optimize debug))
  (assert (>= dense-start xstart))
  (let* ((n (length y))
	 ;;Find number of steps.  It is just the length of the interval divided by the step size, except that
	 ;;if this is infinitesimally below an integer we round up, to make sure the final step is included.
	 (steps (floor (* (/ (- xend dense-start) dense-step) (+ 1.0 (* double-float-epsilon 100)))))
	 (next-step 0))
    (do-dop853 function xstart y xend rtol atol
	       #'(lambda (xold x interpolate-function)
		   (let ((end-step (if (>= x xend) steps ;Integration finished: be sure not to miss final point
				     (floor (- x dense-start) dense-step))))
		     (format t "Solution from ~S to ~S so we do steps ~D through ~D.~%"
			     xold x next-step end-step)
		     (loop for index from next-step to end-step
			   for s = (+ dense-start (* index dense-step))
			   for position = (make-array n :element-type 'double-float)
			   do (dotimes (component n)
				(setf (aref position component)
				      (funcall interpolate-function component s)))
			   do (when (eq (funcall output-function index s position) :abort)
				(return :abort)) ;pass on :abort.  Otherwise NIL.
			   finally (setq next-step (1+ end-step)))))
	       maximize-errors start-step-size)))
    
;;Return path
(defun do-dop853-dense-path (function xstart y xend rtol atol
				      &rest args &key (dense-start xstart) dense-step &allow-other-keys)
  (let* ((steps (floor (* (/ (- xend dense-start) dense-step) (+ 1.0 (* double-float-epsilon 100))))) ;Fudge as above
	(path (make-array (1+ steps))))
    (apply #'do-dop853-dense function xstart y xend rtol atol
	   :output-function #'(lambda (index x y)
				(declare (ignore x))
				(setf (aref path index) y)
				nil)
	   args)
    path))

;;;Test code

;;A simple test of y''=-y, implemented as y1' = y2, y2' = -y1
(defun test-dop853 (&key (xmax (* 2 pi))
			 (tolerance 0.01)
			 (step (/ (* 2 pi) 30))
			 (dense-start 0.0)
			 (plot t)
			 maximize-errors
			 start-step-size)
  (let ((y (make-array 2 :element-type 'double-float :initial-contents '(1.0 0.0)))
	;; (count 0)
	)
    (flet ((derivatives (x y f)
	     (declare (ignore x))
	     ;;(format t "~&~S computing derivatives at ~S~%" (incf count) x)
	     (setf (aref f 0) (aref y 1)
		   (aref f 1) (- (aref y 0)))))
      (let ((path (do-dop853-dense-path #'derivatives 0.0 y xmax
					0.0 ;no relative tolerance
					tolerance ;absolute tolerance
					:dense-step step :dense-start dense-start
					:start-step-size start-step-size
					:maximize-errors maximize-errors)))
	(when plot (plot-dop853-output path dense-start step))
	(describe-dop853-output path dense-start step)
	))))

(defun plot-dop853-output (path start step)
  (gnuplot 4 (array-dimension path 0)
	   #'(lambda (plot point)
	       (if (eq point :title) (nth plot '("y" "y'" "cos" "sin"))
		 (let ((x (+ start (* point step))))
		   (values x
			   (ecase plot
			     ((0 1) (aref (aref path point) plot))
			     (2 (cos x))
			     (3 (- (sin x))))))))
	   :styles '(:points :points :lines :lines)))

;;Return max error
(defun describe-dop853-output (path start step)
  (loop for index below (array-dimension path 0)
	for x = (+ start (* index step))
	maximize (max (abs (- (aref (aref path index) 0) (cos x)))
		      (abs (- (aref (aref path index) 1) (- (sin x)))))))
