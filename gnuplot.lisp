(in-package "CL-USER")
;;;Lisp interface to GNUPLOT
;;;This version copyright (c) 2011 Ken Olum.

;;Plot a bunch of functions or data sets using gnuplot.  All computation
;;of values is done in Lisp.
;;
;;NPLOTS gives the number of plots.  They are numbered from 0 to NPLOTS-1
;;NPTS gives the number of points, from 0 to NPTS-1.
;;
;;FUNC is called in two ways:
;;  With arguments PLOT and :TITLE, it should return a title for the given
;;   plot number.
;;  With arguments PLOT and POINT it should return values X and Y (or X, Y,
;;   and Z for a 3D plot) giving the coordinates of the point, or NIL to
;;   say that the point doesn't exist.

;; :PRINT is a file to receive postscript output instead of displaying
;;  on the screen.
;;For PostScript output, :SOLID and :COLOR set those switches, and
;; :FONT-SIZE give the font size for the key and labels
;; :OUTPUT-SIZE gives option for size in set terminal postscript

;; :TITLE is a title for the plot as a whole
;; :XRANGE and :YRANGE are (MIN MAX) or (MIN . MAX)
;;If :ERRORBARS is given, 3rd value is ydelta
;;LOGSCALE is an axis (:X or :Y) or a list of axes to log-scale
;;If :EXACT-RANGES is set then we specify the range, which prevents gnuplot
;;from expanding it to "good" values.
;; :XTICS  and :YTICS must be lists of 2 or 3 elements
;;giving start, incr, [end] for where axes should be labelled, or else
;;lists of lists ("label" pos) for all the tics.  TICS positions are printed
;;with ~A so you can format the numbers yourself.
;;If not given, tics are computed automatically by gnuplot.
;; If :3D is given, then each call to FUNC should return (X, Y, Z).
;; Other possibilities also exist, according to gnuplot options that you can set in the prelude.
;;:KEY should be a two-element list giving the key location in plot coordinates,
;; a string for "set key", or NIL for no key.
;;STYLES is either a keyword or string for the style to use or a list of such
;;You can include such things as "lines linetype 3"
;;REUSE means to reuse the last gnuplot process (and window).
;;If ANIMATE is set, we make a series of calls for the entire set of plots,
;; put them in files and then run through them in gnuplot.  The value of
;; ANIMATE is the number of frames.  The frame number is passdd as the first
;; argument to FUNC.
;;If :REUSE-DATA is set then we reuse old data files and don't call
;;FUNC except for titles.
;;If NO-GNUPLOT is set, we don't actually call gnuplot, but rather type out
;; what would be done.
;;RETURN-ARGS -- Return as a list the arguments with which gnuplot was called, instead of doing anything.
;;DATA-FILE-PREFIX specifies the directory and filename prefix for datafiles,
;; superseding usual serial-number system.
;;PRELUDE gives gnuplot code to go before plotting commands, but after
;; prelude material we include.
;;POSTLUDE gives code to go after.
;;If ADDITIONAL-PLOTS is given, it is a string to be written as part of the plot command after the last
;; plot and a comma.
;;If EXIT is set, stream will be closed after commands

;;See MULTIPLOT for grouping plots together

;;Data about a gnuplot invocation
(defstruct gnuplot
  serial				;serial number for filenames
  handle				;handle on process (e.g., process object)
  stream				;stream to process
  files					;files used in this process
  )

(defvar *gnuplots* nil)			;Extant gnuplots.  Most recent first.

(defvar *gnuplot-serial* 0)		;Next serial number to use

(defvar *multiplot-stream* nil)		;If set, stream to gnuplot inside MULTIPLOT.

(defvar *gnuplot-remote-host* nil)	;If set, ship files to this host and run gnuplot via ssh
(defvar *gnuplot-remote-user* nil)
;;Where to put files on remote host.  If this is different from where we put the locally, the filenames
;;would need remapping, which is not implemented.
(defvar *gnuplot-remote-directory* "/tmp")

(defun gnuplot (nplots npts func &rest keys
		       &key animate reuse print epsf png no-gnuplot styles no-data
		       postlude exit return-args
		       &allow-other-keys)
  (check-type nplots (integer 1 *))
  (check-type npts (integer 0 *))
  (when return-args
    (return-from gnuplot (list* nplots npts func keys)))
  (when (and reuse (or print epsf png))
    (error "Reusing gnuplot to print or output does not make sense"))
  (when (consp styles)
    (unless (and (listp styles)		;Reasonable list?
		 (= (length styles) nplots))
      (error ":STYLES should be a list of ~D style types" nplots)))
  (unless (and *gnuplots* (gnuplot-handle (car *gnuplots*)))
    (setq reuse nil))			;Nothing to reuse
  (let (gnuplot)
    (cond (reuse			;Reuse previous window
	   (setq gnuplot (car *gnuplots*))
	   (unless no-data (delete-gnuplot-files gnuplot))) ;Delete old datafiles
	  (t
	   (setq gnuplot (make-gnuplot :serial *gnuplot-serial*)) ;Allocate new
	   (incf *gnuplot-serial*)
	   (push gnuplot *gnuplots*)))
    (multiple-value-bind (files flags xrange yrange zrange)
	(if no-data (gnuplot-files gnuplot) ;Don't write data.  Reuse files.  May not work.
	  (apply #'gnuplot-write-data gnuplot animate nplots npts func keys))
      (setf (gnuplot-files gnuplot) files) ;Remember for later
      (when *gnuplot-remote-host*	   ;Need to ship them?
	(gnuplot-ship files))
      (cond (*multiplot-stream*		;multiplotting
	     (setf (gnuplot-stream gnuplot) *multiplot-stream*)) ;Send commands there
	    (no-gnuplot			   ;Just type out what we would do
	     (setf (gnuplot-stream gnuplot) *standard-output*))
	    (reuse)			   ;Previous process?  Nothing to do
	    (t (start-gnuplot gnuplot))) ;Start up gnuplot process
      (when no-gnuplot		;Header if not using process
	(format t "~&Commands to gnuplot:~%"))
      (apply #'gnuplot-header (gnuplot-stream gnuplot) xrange yrange zrange keys)
      (if animate
	  (apply #'gnuplot-animate (gnuplot-stream gnuplot) animate nplots func files keys)
	(apply #'gnuplot-plot-command (gnuplot-stream gnuplot) nil nplots func (car files) flags keys))
      (let ((undone (count 0 flags)))
	(when (plusp undone) (format t "  ~D plots were empty.~%" undone)))
      (when postlude
	(format (gnuplot-stream gnuplot) "~A" postlude))
      (cond ((or print epsf png exit)	;If printing or output requested, flush process
	     (unless (or *multiplot-stream* no-gnuplot)
	       (close (gnuplot-stream gnuplot)) 
	       (wait-for-program (gnuplot-handle gnuplot))
	       (delete-gnuplot-files gnuplot)
	       (setq *gnuplots* (remove gnuplot *gnuplots*))))
	    (t
	     (force-output (gnuplot-stream gnuplot))	;Otherwise keep it around
	     ))))
  t)

;;Do all plots in body in a single gnuplot
;;Options are
;; :ROWS, :COLUMNS -- layout plots in this way
;; :COMMAND-TEXT -- anything to put in the "set multiplot" command
;; :REUSE -- reuse gnuplot
;; :PRELUDE -- commands to go before "set multiplot"
(defmacro multiplot (multiplot-options &body body)
  `(do-multiplot #'(lambda () ,@body) ,@multiplot-options))

(defun do-multiplot (body &key rows columns command-text prelude reuse no-gnuplot)
  (when *multiplot-stream*
    (error "Cannot do recursive multiplotting"))
  (unless (and *gnuplots* (gnuplot-handle (car *gnuplots*)))
    (setq reuse nil))			;Nothing to reuse
  (let ((gnuplot nil))
    (cond (reuse (setq gnuplot (car *gnuplots*)))		;Reuse last gnuplot for this multiplot
	  (t (setq gnuplot (make-gnuplot :serial *gnuplot-serial*))
	     (incf *gnuplot-serial*)))
    (cond (no-gnuplot (setf (gnuplot-stream gnuplot) *standard-output*))
	  (reuse)			;Use old stream
	  (t (start-gnuplot gnuplot)))	;Start new process
    (let ((*multiplot-stream* (gnuplot-stream gnuplot)))
      (when prelude
	(format *multiplot-stream* "~A" prelude))
      (when (or rows columns)
	(setq command-text (format nil "layout ~D,~D~@[ ~A~]" (or rows 1) (or columns 1) command-text)))
      (format *multiplot-stream* "set multiplot~@[ ~A~]~%" command-text) ;Enter multiplot mode
      (funcall body)			;Do plots
      (format *multiplot-stream* "unset multiplot~%")
      (force-output *multiplot-stream*)
      (if reuse
	  (setq *gnuplots* (cons gnuplot (delete gnuplot *gnuplots*))) ;Move to front
	(push gnuplot *gnuplots*))	;Save now, so that it is before subsidiaries
      t)))

;;Send something to the multiplot stream
(defun multiplot-format (ctl-string &rest args)
  (format *multiplot-stream* ctl-string args))
    
;;Execute body with each variable in VARS expanded to the same
;;variable with X, Y or Z on the beginning if variable 3d is set,
;;or X or Y if not.  Also binds variable
;;XYZ to "x", "y" or "z".
(defmacro xyz (vars &body body)
  `(progn
     ,@(loop for xyz in '("X" "Y" "Z")
	    as form = `(let ((xyz ,(string-downcase xyz))
			   ,@(loop for var in vars
				   collect `(,var ,(intern
						    (format nil "~A~A"
							    xyz var)))))
		       ,@body)
	    if (string-equal xyz "Z")
	    collect `(when 3d ,form)
	    else collect form)))
    
;;Give all the commands that come before the plot command
(defun gnuplot-header (stream xlimits ylimits zlimits
		       &key 3d exact-ranges ((:title overall-title)) styles
		            print epsf png portrait color solid font-size output-size
			    (key :default)
			    logscale size-scale
			    xlabel ylabel zlabel xrange yrange zrange
			    xtics ytics ztics ticslevel
			    uniform-scale prelude
		       &allow-other-keys)
  (check-type key (or list symbol string))
  (when (and styles (atom styles))		;All same style?
    (format stream "set style data ~(~A~)~%" styles)) ;Just do it
  ;;Multiple style handled in GNUPLOT-PLOT-COMMAND
  (when (and print epsf)
    (error "Can't set both :PRINT and :EPSF"))
  (when (eq print t)			;Default printing:
    (setq print "|lp"))		;Send straight to printer 
  ;;Set parametric for 3d plots.  It has various beneficial effects.
  (when 3d (format stream "set parametric~%"))
  (when overall-title
    (format stream "set title '~A'~%" overall-title))
  (cond ((eq key :default))		;Let it pick
	((null key)			;suppress
	 (format stream "set nokey~%"))
	((listp key)			;specify coordinates
	 (apply #'format stream "set key at ~F, ~F~%" key))
	(t				;string
	 (format stream "set key ~A~%" key)))
  (when logscale
    (unless (listp logscale) (setq logscale (list logscale)))
    (format stream "set logscale ~(~{~A~}~)~%" logscale))
  (xyz (label) (when label
		 (format stream "set ~Alabel \"~A\"~%" xyz label)))
  (xyz (tics) (when tics
		(format stream "set ~Atics ~A" xyz (tics-format tics))))
  (when ticslevel
    (format stream "set ticslevel ~A~%" ticslevel))
  (cond (uniform-scale
	 (let ((range (max (- (cdr xlimits) (car xlimits))
			   (- (cdr ylimits) (car ylimits))
			   (- (cdr zlimits) (car zlimits)))))
	   (setq xrange (cons (- (/ (+ (cdr xlimits) (car xlimits)) 2) range)
			      (+ (/ (+ (cdr xlimits) (car xlimits)) 2) range))
		 yrange (cons (- (/ (+ (cdr ylimits) (car ylimits)) 2) range)
			      (+ (/ (+ (cdr ylimits) (car ylimits)) 2) range))
		 zrange (cons (- (/ (+ (cdr zlimits) (car zlimits)) 2) range)
			      (+ (/ (+ (cdr zlimits) (car zlimits)) 2) range)))))

	(exact-ranges
	 (setq xrange xlimits yrange ylimits zrange zlimits)))
  (xyz (range) (when range
		 (format stream "set ~Arange [~A:~A]~%" xyz
			 (gnuplot-format-data (car range))
			 (gnuplot-format-data (if (consp (cdr range)) (second range) (cdr range))))))
  (when (or print epsf)
    (format stream "set terminal postscript ~A ~:[monochrome~;color~]~:[~; solid~]~@[ font ~D~]~@[ size ~A~]~%"
	    (cond (epsf "epsf")
		  (portrait "portrait")
		  (t "landscape"))
	    color solid font-size output-size))
  (when png
    (format stream "set terminal png~%"))
  (when (or print epsf png)
    (format stream "set output '~A'~%" (or print epsf png)))
  (when size-scale
    (if (atom size-scale)
	(setq size-scale (list size-scale size-scale)))
    (format stream "set size ~A,~A~%" (gnuplot-format-data (first size-scale))
	    (gnuplot-format-data (second size-scale))))
    (when prelude
      (format stream "~A" prelude)))

;;Emit plot command. 
(defun gnuplot-plot-command (stream frame nplots func files flags
				    &key 3d styles additional-plots &allow-other-keys)
  (loop for plot below nplots
	for file in files
	with first = t
	as title = (if frame (funcall func frame plot :title)
		     (funcall func plot :title))
	for style = (and (consp styles) (pop styles)) ;If styles specified, use next
	when (plusp (bit flags plot))	;If there are any points
	do
	(if first (format stream "~:[plot~;splot~] " 3d)
	  (format stream ", "))
	(format stream "'~A' ~:[notitle~;title ~:*\"~A\"~]~@[ with ~(~A~)~]"
		file title style)
	(setq first nil))
  (when additional-plots
    (format stream ", ~A" additional-plots))
  (terpri stream))


(defun gnuplot-animate (stream nframes nplots func files
			       &rest keys &key pause (slower 1)
			       &allow-other-keys)
  (loop for frame below nframes
	for these-files in files
	do (loop repeat slower do 
		 (apply #'gnuplot-plot-command
			stream frame nplots func these-files keys))
     (force-output stream)
     (when pause
       (when (plusp (length (read-line)))
	 (loop-finish)))))

;;Return string for tics placement.  The TICS argument can be either
;;a list of lists (LABEL COORD) or else a single list (START INCR)
;;or (START INCR END)
(defun tics-format (tics)
  (format nil (if (listp (car tics)) ;List of (label pos)
					;Print as (l p, l p, ...)
		  "(~{~{\"~A\" ~S~}~^, ~})~%"
		"~S, ~S~^, ~S~%") ;No, min, incr, [max]
	  tics))

(defun get-current-username ()
  #+cmu(cdr (assoc :user *environment-list*))
  #+sbcl(posix-getenv "USER"))

;;Calls user function and writes the data to files for gnuplot.
;; DATA-FILE-PREFIX gives location and filename for datafiles
;; The function can return any number of values, and we will write them.  We don't check what this number should be
;;Returns values:
;; FILES -- list of filenames used for plots
;; FLAGS -- bit vector telling whether there are any data points in the corresponding plot (at any time if animate)
;; followed by (minimum . maximum) of all values returned by the user's function
(defun gnuplot-write-data (gnuplot animate nplots npts func
				   &key styles
				   (data-file-prefix
				    (format nil "/tmp/~A-gnuplot-~D-~D.tmp."
					    (get-current-username)
					    (getpid)
					    (gnuplot-serial gnuplot)))
				   &allow-other-keys)
  (loop with flags = (make-array nplots :element-type 'bit)
	with limits = (make-array 10 :fill-pointer 0 :adjustable t) ;Surely we won't really need to adjust it.
	for frame below (or animate 1)
	collect
    (loop for plot below nplots
	  as style = (if (atom styles) styles (pop styles))
	  as filename = (format nil "~A~:[~*~;~D-~]~D"
				data-file-prefix animate frame plot)
	  collect filename
	  do
       (with-open-file (stream filename :direction :output :if-exists :supersede)	;Output file for this plot
	 (loop for pt below npts
	       as data = (multiple-value-list  ;List of data: (x y), (x y z), (x y color), etc.
			  (if animate (funcall func frame plot pt)
			    (funcall func plot pt)))
	       when (first data) ;All this if point is supplied.  If not, just output empty line to cause break in plot.
	       do
	       (unless (numberp (second data))
		 (error "Invalid y value ~S trying to plot ~S, ~S"
			(second data) plot pt))
	       (setf (bit flags plot) 1) ;Say that this plot is used
	       (format stream "~{~A~^ ~}" (mapcar #'gnuplot-format-data data))
	       ;;Maintain minimum and maximum values of numeric data
	       (loop for index from 0
		     for value in data
		     unless (> (fill-pointer limits) index) do (vector-push-extend nil limits) ;Make sure big enough
		     when (numberp value)
		     do (if (aref limits index)						  ;previous data?
			    (setf (car (aref limits index)) (min (car (aref limits index)) value) ;elt is (min . max)
				  (cdr (aref limits index)) (min (cdr (aref limits index)) value))
			  (setf (aref limits index) (cons value value)))) ;First time
	       do (format stream "~%"))))
    into files
    finally (return (values-list (list* files flags (coerce limits 'list))))))

;;host or user@host
(defun gnuplot-remote-spec ()
  (format nil "~@[~A@~]~A" *gnuplot-remote-user* *gnuplot-remote-host*))

;;Start up the gnuplot process
(defun start-gnuplot (gnuplot)
  (multiple-value-bind (handle stream)
      (if *gnuplot-remote-host*
	  (do-run-program "ssh" :args (list (gnuplot-remote-spec) "gnuplot" "-display" ":0")  :input :stream :wait nil)
	(do-run-program "gnuplot" :input :stream :wait nil))
    (setf (gnuplot-handle gnuplot) handle
	  (gnuplot-stream gnuplot) stream)))

(defun kill-gnuplots ()
  (loop for gnuplot in (reverse *gnuplots*) ;Close oldest first
	for handle = (gnuplot-handle gnuplot)
	when handle			;Unless didn't actually run process
	do (do-kill-program handle)
	do (delete-gnuplot-files gnuplot))
  (setq *gnuplots* nil
	*gnuplot-serial* 0))

;;Delete files on exit
(defun delete-gnuplot-files (gnuplot)
  (dolist (frame (gnuplot-files gnuplot)) ;A list of files for each animation frame, even if not animating
    (dolist (file frame)
      (handler-case (delete-file file)
	 (file-error (condition) (warn "Failed to delete ~A: ~A" file condition))))))
    
(defun last-gnuplot-files (&optional (nth 0))
  (gnuplot-files (nth nth *gnuplots*)))

;;Format a data item for gnuplot.  If it is a float, we print it in exponential notation with an E
(defun gnuplot-format-data (x)
  (typecase x
    (float (format nil "~,,,,,,VE" #\E x))
    (integer (format nil "~D" x))
    (string x)
    (t (error "Don't know how to format data item ~S for gnuplot" x))))


;;Rotate an existing plot by continuously replotting with different views
(defun gnuplot-rotate (fromxrot toxrot fromzrot tozrot steps &optional pause)
  (let ((s (gnuplot-stream (car *gnuplots*))))
    ;;Tics are distracting when rotating
    (format s "set noxtics~%set noytics~%set noztics~%")
    (loop repeat steps
	  for xrot from fromxrot by (/ (- toxrot fromxrot) steps)
	  for zrot from fromzrot by (/ (- tozrot fromzrot) steps)
	  do (format s "set view ~A, ~A~%"
		     (gnuplot-format-data xrot) (gnuplot-format-data zrot))
	     (format s "replot~%")
	     (force-output s)
	     (when pause (sleep pause)))
    (format s "set xtics~%set ytics~%set ztics~%")
    (format s "replot~%")
    (force-output s)
    ))
  
  




;;Plot a sequence of data, which should consist of pairs (X Y).
(defun plot-data-list (data &optional (title "Data list")
			    &rest gnuplot-keys)
  (apply #'gnuplot 1 (length data)
	 #'(lambda (plot point)
	     (declare (ignore plot))
	     (if (eq point :title)
		 title
	       (values-list (elt data point))))
	 gnuplot-keys))
	       
;;Plot a sequence of Y values for integer X.
(defun plot-list (list &rest title-and-gnuplot-keys)
  (apply #'plot-data-list
	 (loop for x from 1
	       for y in (coerce list 'list)
	       collect (and y (list x y)))
	 title-and-gnuplot-keys))

(defun test-plot (&rest args)
  (apply #'gnuplot 2 20
	 #'(lambda (plot pt)
	     (if (eq pt :title)
		 (ecase plot
		   (0 "sin(x)")
		   (1 "cos(x)"))
	       (let ((x (* 2 pi (/ pt 20.0))))
		 (values x
			 (ecase plot
			   (0 (sin x))
			   (1 (cos x)))))))
	 args))



;;Plot a surface given by a two-dimensional array of data.  The indices to the array are X and Y
(defun plot-data-surface (data &rest gnuplot-keys
			       &key (styles :pm3d)
			       (view-x 0)
			       (view-z 0)
			       &allow-other-keys)
  (let ((xpoints (array-dimension data 0))
	(ypoints (array-dimension data 1)))
    (apply #'gnuplot
	   1
	   (* xpoints (1+ ypoints))
	   #'(lambda (plot point)
	       (declare (ignore plot))
	       (unless (eq point :title)
		 (multiple-value-bind (x y)
		     (floor point (1+ ypoints))
		     (unless (= y ypoints) ;NIL to break at end of each line
		       (values-list (list x y (aref data x y)))))))
	   :styles styles
	   :prelude (format nil "set contour~%set view ~D, ~D~%" view-x view-z)
	   :3d t
	   gnuplot-keys)))

;;Plot the values of a matrix.  The first index is the row and the second index the column, and (0,0) is
;;in the upper left.  This is a 90 degree rotation from the default in plot-data-surface
(defun plot-matrix (m &rest gnuplot-keys)
  (apply #'plot-data-surface m :view-x 0 :view-z 90 gnuplot-keys))

;;Make sure files are deleted on exit
(pushnew 'kill-gnuplots *exit-hooks*
	 :test #'equal)

;;Remote plotting

;;Copy local files to same names (usually on /tmp) on remote host
(defun gnuplot-ship (files)
  (format t "~&Shipping ~D file~:P to ~A..." (length files) *gnuplot-remote-host*) (terpri)
  (do-run-program "scp"
		  :args (append (reduce #'append files)
				(list (format nil "~A:~A/" (gnuplot-remote-spec) *gnuplot-remote-directory*))))
  (format t "done.~%"))


