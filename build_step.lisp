;(restrict-compiler-policy 'speed  0 0)
;(restrict-compiler-policy 'debug  3 3)
;(restrict-compiler-policy 'safety 3 3)
(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
(ql:quickload :cl-mpm-worker)
;; (ql:quickload :cl-mpm)
;; (ql:quickload :cl-mpm/damage)
;; (ql:quickload :cl-mpm/mpi)
(ql:quickload :cl-mpm/examples/chalk)
(ql:quickload :cl-mpm/mpi)
(in-package :cl-mpm-worker)
;; (in-package :cl-mpm/examples/chalk)
;; ;(asdf:compile-system :cl-mpm :force t)
;; ;(asdf:compile-system :cl-mpm/damage :force t)
;; ;(asdf:compile-system :cl-mpm/mpi :force t)
;; (asdf:compile-system :cl-mpm/examples/chalk :force t)
;; (in-package :cl-mpm-worker)

(sb-ext:save-lisp-and-die
   "mpi-worker"
    :executable t
    :toplevel #'main
    :save-runtime-options t)
(uiop:quit)
