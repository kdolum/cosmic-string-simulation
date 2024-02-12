;;Recompile everything
(let ((*recompile* t))
  (declare (special *recompile*))
  (load (merge-pathnames "load" *load-pathname*)))
