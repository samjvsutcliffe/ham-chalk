;(sb-ext:restrict-compiler-policy 'speed 3 3)
;(sb-ext:restrict-compiler-policy 'debug 0 0)
;(sb-ext:restrict-compiler-policy 'safety 0 0)
;(setf *block-compile-default* t)
;(sb-ext:restrict-compiler-policy 'speed  0 0)
;(sb-ext:restrict-compiler-policy 'debug  3 3)
;(sb-ext:restrict-compiler-policy 'safety 3 3)
;(setf *load-verbose* nil)
;(setf *load-print* nil)
;(setf *compile-verbose* nil)
;(setf *compile-print* nil)


(ql:quickload "cl-mpm/examples/chalk")
(in-package :cl-mpm/examples/chalk)
;(ql:quickload "magicl")
;(ql:quickload "cl-mpm")
;;(asdf:compile-system :cl-mpm :force T)
;(ql:quickload "cl-mpm/examples/chalk")
;(ql:quickload "cl-mpm/mpi")
;(in-package :cl-mpm/examples/chalk)

(defun plot (sim)
  (cl-mpm/plotter:simple-plot
   *sim*
   :plot :deformed
   ;; :colour-func (lambda (mp) (cl-mpm/utils:get-stress (cl-mpm/particle::mp-stress mp) :zz))
   :colour-func (lambda (mp) (cl-mpm/particle::mp-damage mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-damage-ybar mp))
   ;; :colour-func (lambda (mp) (cl-mpm/particle::mp-strain-plastic-vm mp))
   ))

(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               #'cl-mpm/shape-function:make-shape-function-bspline
               ;'cl-mpm::mpm-sim-usf
               ;'cl-mpm/damage::mpm-sim-damage-nd-2
               'cl-mpm/mpi::mpm-sim-mpi-nodes
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let ()
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                ;; (mapcar (lambda (x) 0) size)
                offset
                ;; '(0 0)
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-chalk-brittle
                ;;'cl-mpm/particle::particle-chalk-delayed
                ;;:E 1d9
                ;;:nu 0.35d0
                ;;:rho 1d6
                ;;:enable-plasticity nil
                ;;:coheasion 1d5
                ;;:friction-angle 40d0
                ;;:fracture-energy 5000d0
                ;;:initiation-stress 500d3
                ;;:delay-time 1d3
                ;;:critical-damage 1d0;0.999d0
                ;;:local-length 10d0
                ;;:local-length-damaged 0.01d0

                'cl-mpm/particle::particle-mc
                :E 1d9
                :nu 0.35d0
                :psi (* 30d0 (/ pi 180))
                :phi (* 30d0 (/ pi 180))
                :c 0.10d6


                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-15)
      (let ((ms 1d5))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-2 density ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) 10.0d0)
        (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 ms))
        )

      (dotimes (i 0)
        (dolist (dir (list :x :y))
          (cl-mpm::split-mps-criteria
           sim
           (lambda (mp h)
             (when
                 (and
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     100)
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     200
                     )
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                     100
                     )
                  )
               dir
               )))))

      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             ;; (lambda (i) nil)
             ;; (lambda (i) nil)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))

      ;(defparameter *floor-bc*
      ;  (cl-mpm/penalty::make-bc-penalty-point-normal
      ;   sim
      ;   (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
      ;   (cl-mpm/utils:vector-from-list (list 00d0 (+ h-y) 0d0))
      ;   (* density 1d5)
      ;   0.9d0
      ;   ;; 1d1
      ;   ))

      ;; (setf (cl-mpm::sim-bcs-force-list sim)
      ;;       (list
      ;;        (cl-mpm/bc:make-bcs-from-list
      ;;         (list
      ;;          *floor-bc*
      ;;          ))))

      ;; (let* ((terminus-size (second block-size))
      ;;        (ocean-y (* terminus-size 1.0d0)))
      ;;   (setf (cl-mpm::sim-bcs-force-list sim)
      ;;         (list
      ;;          (cl-mpm/bc:make-bcs-from-list
      ;;           (list
      ;;            (cl-mpm/buoyancy::make-bc-buoyancy-clip
      ;;             sim
      ;;             ocean-y
      ;;             1000d0
      ;;             (lambda (pos datum)
      ;;               t)
      ;;             ))))))

      sim)))


(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (c cl-mpm/particle::mp-c))
      mp
    (let ((rho_0 0.2d6)
          (rho_1 0.0020d6)
          (soft 5d-1))
      (setf c (max rho_1
                     (* rho_0 (exp (- (* soft ps)))))))
    )
  )

(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))  
    ))

(defun setup (&key (undercut 0d0))
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))
  (format t "Setup ~%")
    (let* ((mesh-size 25)
           (mps-per-cell 2)
           (shelf-height 100)
           (soil-boundary 50)
           (shelf-aspect 1.5)
           (runout-aspect 3)
           (shelf-length (* shelf-height shelf-aspect))
           (domain-length (+ shelf-length (* runout-aspect shelf-height)))
           (shelf-height (+ shelf-height soil-boundary))
           (offset (list 0 (* 0 mesh-size))))
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height 20))
                         (list domain-length shelf-height)
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))

    (let* ((undercut-angle undercut)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1))))
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (> (magicl:tref p 1 0) soil-boundary)
                                     (cl-mpm/setup::plane-point-sdf
                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                              (magicl:tref p 1 0)) '(2 1))
                                      normal
                                      (magicl:from-list (list shelf-length soil-boundary)
                                                        '(2 1) :type 'double-float))
                                     1d0)
                                 )))
    ;; (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
    ;;                                              (magicl:from-list (list (magicl:tref p 0 0)
    ;;                                                                      (magicl:tref p 1 0)) '(2 1))
    ;;                                              (list (- shelf-length (* 0.5 shelf-height)) shelf-height)
    ;;                                              (list shelf-length soil-boundary)
    ;;                                              30d0
    ;;                                              )) 1.00d0)
    )

    ;(loop for mp across (cl-mpm:sim-mps *sim*)
  ;                 do (setf (cl-mpm/particle::mp-damage mp) (random 0.01d0) ))
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))


(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  &allow-other-keys;; (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
  ;;                                           *sim*
  ;;                                           (list :dt target-time))

  (let* ((target-time 1d2)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (when (= steps 5)
                       ;(setf (cl-mpm::sim-enable-damage *sim*) t)
                       (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                        (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-8 ms)))
                       )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (time
                      (dotimes (i substeps);)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                     (incf *sim-step*)
                     ;(plot *sim*)
                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(defun mpi-loop ()
  (format t "Starting mpi~%")
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (setup :undercut 0d0)
    ;(setup)

    ;;Square
    ;(let ((dsize (floor (sqrt (cl-mpi:mpi-comm-size)))))
    ;  (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize dsize 1)))
    (let ((dsize (floor (cl-mpi:mpi-comm-size))))
      (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize 1 1)))

    ;; (let ((dhalo-size (* 4 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
    ;;  (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))

    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Decompose~%"))
    (cl-mpm/mpi::domain-decompose *sim*)
    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*))))
    (when (= rank 0)
      (format t "Run mpi~%"))
    (run-mpi)
    (when (= rank 0)
      (format t "Done mpi~%"))
    )
  )

(defun run-mpi ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (target-time 1d3)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
      (when (= rank 0)
        (format t "CFL dt estimate: ~f~%" dt-e)
        (format t "CFL step count estimate: ~D~%" substeps-e))
      (setf substeps substeps-e))
    (when (= rank 0)
      (format t "Substeps ~D~%" substeps)
 	    (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                             *sim*
                                             (list :dt target-time)) )
    (time (loop for steps from 0 to 500
                while *run-sim*
                do
                   (progn
                     (when (= steps 5)

                       (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                        (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-8 ms)))
                       ;(setf (cl-mpm::sim-enable-damage *sim*) t)
                       ;; (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                       ;;   (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-3 ms)))
                       )
                     (when (= rank 0)
                       (format t "Step ~d ~%" steps))
                     ;(cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                     (time
                      (dotimes (i substeps);)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                     (incf *sim-step*)
                     ;(plot *sim*)
                     (swank.live:update-swank)
                     )))
    )
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

;(defun run-undercut ()
;  (let* ((target-time 1d3)
;         (dt (cl-mpm:sim-dt *sim*))
;         (substeps (floor target-time dt))
;         (dt-scale 1d0))
;
;    (cl-mpm::update-sim *sim*)
;    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
;                    (format t "CFL dt estimate: ~f~%" dt-e)
;                    (format t "CFL step count estimate: ~D~%" substeps-e)
;                    (setf substeps substeps-e))
;    (format t "Substeps ~D~%" substeps)
;    (time (loop for steps from 0 to 10
;                while *run-sim*
;                do
;                   (progn
;                     (format t "Step ~d ~%" steps)
;                     (time
;                      (dotimes (i substeps)
;                        (cl-mpm::update-sim *sim*)
;                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
;                     (incf *sim-step*)
;                     )))))
;
;(defun test-undercut ()
;
;  (setup)
;  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
;                               *sim*)
;  (loop for c in (list 0d0 10d0 20d0 30d0 40d0 50d0)
;        while *run-sim*
;        do
;           (progn
;             (setup :undercut (- c))
;             (loop for mp across (cl-mpm:sim-mps *sim*)
;                   do
;                      (setf
;                       (cl-mpm/particle::mp-initiation-stress mp) 1d20
;                       (cl-mpm/particle::mp-gf mp) 1d20))
;             (run-undercut)
;             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/chalk_~5,'0d.vtk" c)) *sim*)
;             )))

;(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
;(test-undercut)

(setf lparallel:*kernel* (lparallel:make-kernel 16 :name "custom-kernel"))
;(defparameter *run-sim* nil)
;(setup)
;(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;(run)

(format t "Running~%")
;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
(mpi-loop)


