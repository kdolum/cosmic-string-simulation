;;;Modify GC policy for simulation
(in-package "CL-USER")
(use-package "SB-ALIEN")
(eval-when (:compile-toplevel :load-toplevel :execute)
  (import (intern "OS-VM-SIZE-T"  "SB-KERNEL")))
(define-alien-variable "gencgc_verbose" int)

;;For verbose output, do this:
;;(setq gencgc-verbose 2)

(define-alien-routine "print_generation_stats" void
    (verbose int))

(define-alien-type page-index-t long)

(define-alien-type generation
  (struct generation
    (alloc-start-page page-index-t)
    (alloc-unboxed-start-page page-index-t)
    (alloc-large-start-page page-index-t)
    (alloc-large-unboxed-start-page page-index-t)
    (bytes-allocated os-vm-size-t)		;Size of space
    (gc-trigger os-vm-size-t)			;Size at which to collect (if age is high enough)
    (bytes-consed-between-gc os-vm-size-t)	;Trigger is set to this plus size after collection
    (number-of-gcs int)			;Number of times this generation collected since last raise
    (number-of-gcs-before-promotion int)		;Number of GCs after which we raise
    (cum-sum-bytes-allocated os-vm-size-t)	;Total of sizes after garbage collection since last raise
    (min-av-mem-age double)		;Minimum average age at which to collect
    ))

(defconstant num-generations 8)	;6 normal generations, pseudo-static, and scratch

(define-alien-variable generations (array generation #.num-generations))

(defmacro generation-slot (generation slot)
  `(slot (deref generations ,generation) ',slot))

;;Automatic GCs during simulation do not raise objects higher than this.
(defparameter highest-simulation-generation 2)

(defvar highest-generation-gc-count)	;Last value of number-of-gcs in this generation

;;Set up GC policy during simulation.  Objects are raised through our highest generation
;;by the normal process.  But then, we schedule garbage collection whenevifer the new
;;number of bytes raised into our generation exceeds half the amount that remain after the
;;previous GC.
(defun setup-simulation-gc-policy ()
  (declare (optimize (debug 0)))	;Avoid alien optimization problems
  (setf (bytes-consed-between-gcs) (* 200 (expt 2 20)))	;200MB to reduce scavenging
  ;;Never automatically raise beyond this
  (setf (generation-slot highest-simulation-generation number-of-gcs-before-promotion) 1000000000)
  ;;Don't do a GC immediately upon filling this region, but after that do it if trigger was reached
  (setf (generation-slot highest-simulation-generation min-av-mem-age) 0.0001d0)
  (pushnew 'maybe-schedule-next-gc *after-gc-hooks*) ;Maintain trigger for our policy
  (without-gcing
   (setq highest-generation-gc-count 0)
   (setf (generation-slot highest-simulation-generation number-of-gcs) 0)
   (schedule-next-gc))
  )

;;If our region was collected, scheduled next
(defun maybe-schedule-next-gc ()
   ;;If generation collected, NUMBER-OF-GCS increases by 1.  If raised (only happens explicitly), set to 0
   (unless (= (generation-slot highest-simulation-generation number-of-gcs) highest-generation-gc-count)
     (schedule-next-gc)))

;;We garbage collect next when our generation has grown by half.
(defun schedule-next-gc ()
  (declare (optimize (debug 0)))	;Avoid alien optimization problems
   (let (next)
     (without-gcing
      (setq next (max (expt 2 21)		;2MB at least
		      (ceiling (* (generation-slot highest-simulation-generation bytes-allocated) 3) 2)))
      (setf (generation-slot highest-simulation-generation gc-trigger) next)
      (setq highest-generation-gc-count (generation-slot highest-simulation-generation number-of-gcs)))
     (format t "~&Next GC in region ~D at ~DMB~%" highest-simulation-generation (round (/ next (expt 2 20))))))

(setup-simulation-gc-policy)

(define-alien-type page
  (struct page
    (region-start-offset unsigned-long)
    (bytes-used unsigned-long)
    (flags char)
    (gen char)))

(define-alien-variable "page_table" (* page))

(define-alien-type alloc-region
  (struct alloc-region
    (free-pointer (* t))
    (end-addr (* t))
#-sbcl2    (first-page page-index-t)
    (last-page page-index-t)
    (start-addr (* t))))

#-sbcl2 (define-alien-variable "boxed_region" alloc-region)
#-sbcl2 (define-alien-variable "unboxed_region" alloc-region)

;;In later versions, these things are in a single C array
#+sbcl2 (define-alien-variable "gc_alloc_region" (array alloc-region 3))
#+sbcl2 (defparameter boxed-region (deref gc-alloc-region 0))
#+sbcl2 (defparameter unboxed-region (deref gc-alloc-region 1))

(define-alien-variable "bytes_allocated" unsigned-long)

(locally (declare (muffle-conditions compiler-note)) ;Avoid redefinition warnings

;;Redefine this that it includes the current allocation regions
;;Unfortunately this doesn't work perfectly, because the previous definition was declared inline.
(defun sb-kernel:dynamic-usage ()
  (+ bytes-allocated
     (sap- (alien-sap (slot unboxed-region 'free-pointer)) (alien-sap (slot unboxed-region 'start-addr)))
     (sap- (alien-sap (slot boxed-region 'free-pointer)) (alien-sap (slot boxed-region 'start-addr)))))


(defun get-bytes-consed ()
  (+ (sb-kernel:dynamic-usage)
     sb-kernel::*n-bytes-freed-or-purified*))

)
