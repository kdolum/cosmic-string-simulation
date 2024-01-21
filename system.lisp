;;;Compile and load files in order, without confusion by multiprocessing
(in-package "CL-USER")

(eval-when (:compile-toplevel)
  (error "Don't compile this file.  Load the source."))

;;Load a system.  We do all the work in the directory of the pathname DIRECTORY
(defun load-system (directory files &optional compile)
  (with-compilation-unit ()
   (loop with timestamp-file = (make-pathname :name "compilation-timestamp" :type nil
					      :defaults directory)
	 with compiled-anything = nil
	 for file in files
	 for source = (merge-pathnames (merge-pathnames file (make-pathname :type "lisp")) directory)
	 for binary = (compile-file-pathname source)
	 for source-time = (file-write-date source)
	 for binary-time = (and (probe-file binary) (file-write-date binary))
	 maximize source-time into latest-source-time
	 do (cond ((and (not compile)
			binary-time
			(>= binary-time latest-source-time)) ;Compiled recently enough?
		   (load binary))	;Good, just load
		  (t			;Need to compile
		   (setq compiled-anything t)
		   (compile-file source :output-file binary)
		   (load binary)))
	 finally
	 (when compiled-anything      ;Was anything compiled?
	   (with-open-file (stream timestamp-file ;write new timestamp
				   :direction :output :if-exists :supersede
				   :if-does-not-exist :create)))
	 )))
