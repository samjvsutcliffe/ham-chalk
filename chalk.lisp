(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)

(ql:quickload "magicl")
(ql:quickload "cl-mpm")
(asdf:compile-system :cl-mpm :force T)
(ql:quickload "cl-mpm/examples/chalk")
(in-package :cl-mpm/examples/chalk)

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
               ;; 'cl-mpm::mpm-sim-usf
               'cl-mpm/damage::mpm-sim-damage
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
                'cl-mpm/particle::particle-chalk-brittle
                ;; 'cl-mpm/particle::particle-elastic
                ;; 'cl-mpm/particle::particle-vm
                :E 1d9
                :nu 0.30d0
                :rho 1d7

                :coheasion 1d4
                :friction-angle 40d0

                :fracture-energy 3d3
                :initiation-stress 1d6
                ;; :initiation-stress 0d0
                ;; :damage-rate 1d-15

                :damage-rate 1d-6
                :critical-damage 0.50d0
                :local-length 10d0
                ;; :local-length-damaged 10d0
                ;; :local-length-damaged 50d0
                :local-length-damaged 1d-5
                ;; :local-length-damaged 1d-5
                ;; :damage 0.0d0
                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                ))))
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
      (let ((ms 1d4))
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
                     120)
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     200 
                     )
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                     10
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
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
             ;; (lambda (i) nil)
             ;; (lambda (i) nil)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))

      (defparameter *floor-bc*
        (cl-mpm/penalty::make-bc-penalty-point-normal
         sim
         (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
         (cl-mpm/utils:vector-from-list (list 00d0 (+ h-y) 0d0))
         (* density 1d5)
         0.9d0
         ;; 1d1
         ))

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

  (let* ((mesh-size 5)
         (mps-per-cell 4)
         (shelf-height 100)
         (soil-boundary 50)
         (shelf-aspect 3)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* 1 shelf-height)))
         (shelf-height (+ shelf-height soil-boundary))
         (offset (list 0 (* 0 mesh-size)))
         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height 100))
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
    ;; (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;       do
    ;;          (setf (cl-mpm/particle:mp-damage mp) (random 0.1d0)))
    ;; (cl-mpm/setup::damage-sdf
    ;;  *sim*
    ;;  (lambda (p)
    ;;    (cl-mpm/setup::line-sdf (magicl:from-list (list (magicl:tref p 0 0)
    ;;                                                    (magicl:tref p 1 0))
    ;;                                              '(2 1))
    ;;                            (list (- shelf-length shelf-height) shelf-height)
    ;;                            (list shelf-length soil-boundary)
    ;;                            10d0
    ;;                            )) 0.8d0)
    ;(let ((sdf
    ;        (lambda (p)
    ;          (cl-mpm/setup::line-sdf (magicl:from-list (list (magicl:tref p 0 0)
    ;                                                          (magicl:tref p 1 0))
    ;                                                    '(2 1))
    ;                                  (list (- shelf-length shelf-height) shelf-height)
    ;                                  (list shelf-length 0d0)
    ;                                  20d0
    ;                                  ))
    ;        ))
    ;  (loop for mp across (cl-mpm:sim-mps *sim*)
    ;        do (with-accessors ((pos cl-mpm/particle:mp-position)
    ;                            (damage cl-mpm/particle:mp-damage)) mp
    ;             (when (>= 0 (funcall sdf pos))
    ;               (setf damage (min 1d0 (max 0d0 (coerce (* (funcall sdf pos) -0.1d0) 'double-float)))))
    ;             )))


    )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  ;(loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))


(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  (let* ((target-time 1d3)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 500
                while *run-sim*
                do
                   (progn
                     (when (= steps 5)
                       (setf (cl-mpm::sim-enable-damage *sim*) t)
                       (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                         (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-3 ms))))
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (time
                      (dotimes (i substeps);)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                     (incf *sim-step*)
                     (plot *sim*)
                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))

(defun run-undercut ()
  (let* ((target-time 1d3)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0))

    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 10
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (time
                      (dotimes (i substeps)
                        (cl-mpm::update-sim *sim*)
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))
                     (incf *sim-step*)
                     )))))

(defun test-undercut ()

  (setup)
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                               *sim*)
  (loop for c in (list 0d0 10d0 20d0 30d0 40d0 50d0)
        while *run-sim*
        do
           (progn
             (setup :undercut (- c))
             (loop for mp across (cl-mpm:sim-mps *sim*)
                   do
                      (setf
                       (cl-mpm/particle::mp-initiation-stress mp) 1d20
                       (cl-mpm/particle::mp-gf mp) 1d20))
             (run-undercut)
             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/chalk_~5,'0d.vtk" c)) *sim*)
             )))

(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
(test-undercut)

;; (setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
;; (require :sb-sprof)
;; (defparameter *run-sim* nil)
;; (setup)
;; (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;; (run)

