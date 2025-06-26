(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/joss)
(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(declaim (optimize (debug 0) (safety 0) (speed 3)))

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  (cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  ;(cl-mpm::update-stress-kirchoff-noscale mesh mp dt fbar)
  ;(cl-mpm::scale-domain-size mesh mp)
  )
(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  (cl-mpm::update-domain-det mesh mp dt)
  ;; (cl-mpm::co-domain-corner-2d mesh mp dt)
  ;; (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;; (cl-mpm::scale-domain-size mesh mp)
  )

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  ;local-length
  )
(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     (ft cl-mpm/particle::mp-ft)
                     (fc cl-mpm/particle::mp-fc)
                     (E cl-mpm/particle::mp-e)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     ) mp
      (declare (double-float pressure damage))
      (progn
        ;(setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
        (setf damage-increment (max 0d0 (cl-mpm/damage::criterion-mohr-coloumb-stress-tensile stress (* angle (/ pi 180d0)))))
		;(setf damage-increment
        ;      (max 0d0
        ;           (cl-mpm/damage::criterion-dp-coheasion
        ;            stress
        ;            ;(magicl:scale stress (/ 1d0 (magicl:det def))) 
        ;            (* angle (/ pi 180d0)))))
        ;(setf damage-increment (cl-mpm/damage::criterion-max-principal-stress stress))

        ;(incf damage-increment
        ;      (* E (cl-mpm/particle::mp-strain-plastic-vm mp)))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))


(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; :sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes
               ;:sim-type 'cl-mpm/mpi::mpm-sim-usl-mpi-nodes-damage
               :sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes-damage
               ;:sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* (;(length-scale (* 2 h))
             ;(length-scale 0.5d0)
             (E 1d9)
             ;(length-scale (* 2 h))
             (length-scale 0.5d0)
             (init-c 26d3)
             (angle 50d0)
             ;(init-c 131d3)
             ;(angle 42d0)
             (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile init-c (* angle (/ pi 180))))
             ;(gf 48d0)
             (gf 48d0)
             ;(gf 4.8d0)
             ;(gf 10d0)
             (ductility (estimate-ductility-jirsek2004 gf length-scale init-stress E))
             (oversize (cl-mpm/damage::compute-oversize-factor 0.99d0 ductility))
             )
		(when (= (cl-mpi:mpi-comm-rank) 0) 
			(format t "Estimated ductility ~E~%" ductility)
			(format t "Estimated oversize ~E~%" oversize)
			(format t "Estimated gf ~E~%" (cl-mpm/damage::gf-from-ductility ductility length-scale init-stress E))
			(when (< ductility 1d0)
			  (error "Ductility too low ~A" ductility)))

        (cl-mpm:add-mps
         sim
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density

                'cl-mpm/particle::particle-chalk-delayed
                :E 1d9
                :nu 0.24d0
                :enable-plasticity t
                :enable-damage t

                :ft 1d0
                :fc 10d0 
                :friction-angle angle

                :kt-res-ratio 1d0
                :kc-res-ratio 0d0
                :g-res-ratio 0.51d0
                ;:g-res-ratio 0.6d0
                :peerlings-damage t

                :fracture-energy 3000d0

                :initiation-stress init-stress
                :delay-time 1d0
                :delay-exponent 2d0
                :ductility ductility

                ;:critical-damage 1d0;(- 1.0d0 1d-3)
                :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage
                :local-length length-scale
                :local-length-damaged 10d-10

                :psi (* 0d0 (/ pi 180))
                :phi (* angle (/ pi 180))
                :c (* init-c oversize)
                :softening 0d0
                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))
                ))))

      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm/damage::sim-enable-length-localisation sim) t)
      ;(setf (cl-mpm::sim-velocity-algorithm sim) :BLEND)
      (setf (cl-mpm::sim-velocity-algorithm sim) :BLEND)
      ;(setf (cl-mpm::sim-velocity-algorithm sim) :FLIP)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-10)
      ;(cl-mpm/setup::set-mass-filter sim density :proportion 1d-2)
            
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim) 
		      ;(* 0.70d0)
              (* 0.1d0 
                 (sqrt (cl-mpm::sim-mass-scale sim))
                 (cl-mpm/setup::estimate-critical-damping sim))))

      ;; (dotimes (i 2)
      ;;   (dolist (dir (list :x :y))
      ;;     (cl-mpm::split-mps-criteria
      ;;      sim
      ;;      (lambda (mp h)
      ;;        (when
      ;;            (and
      ;;             (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
      ;;                80)
      ;;             (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
      ;;                200
      ;;                )
      ;;             (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
      ;;                50
      ;;                )
      ;;             )
      ;;          dir
      ;;          )))))

      (setf (cl-mpm:sim-dt sim) (cl-mpm/setup::estimate-elastic-dt sim :dt-scale 0.1d0))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
	 (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-varfix
             (cl-mpm:sim-mesh sim)
             '(0 nil 0)
             '(0 nil 0)
             '(0 0 0)
             '(0 0 0)
             '(nil nil 0)
             '(nil nil 0)))
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

(defun setup ()
  (let* ((mesh-size 0.5)
         (mps-per-cell 3)
         (mp-refine 0)
         (shelf-height 15.5)
         (soil-boundary 1)
         (shelf-aspect 1)
         (runout-aspect 1)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
         (shelf-height (+ shelf-height soil-boundary))
         (depth 400)
         (offset (list 0 (* 0 mesh-size)
                       ;; 0
                       ))
         )
    (format t "Mesh size: ~F~%" mesh-size)
    (defparameter *sim*
      (setup-test-column (list domain-length
                               (+ shelf-height (* 5 mesh-size))
                               ;; depth
                               )
                         (list domain-length shelf-height
                               ;; depth
                               )
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))

    ;;Refine around tip
	(dotimes (i 0)
      (dolist (dir (list :x
                         :y
                         ))
        (cl-mpm::split-mps-criteria
         *sim*
         (lambda (mp h)
           (when
               (and
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length
                      15d0
                      ))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (- shelf-length 6d0))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   (+ shelf-length 2d0))
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   (- soil-boundary 2d0)
                   ))
             dir)))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 6.8d0))
           (measured-angle 78d0)
           (undercut-angle ;(- 82.5d0 90d0)
             (- measured-angle 90d0)
                           )
           ;; (undercut-angle 0d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list (- shelf-length (* (tan (- (* pi (/ undercut-angle 180d0)))) sloped-height))
                                     soil-boundary)
                               '(2 1) :type 'double-float)

             )
           (edge-refine 1)
           )
    (flet ((cutout (p)
               (if (and
                    (> (magicl:tref p 1 0) soil-boundary))
                   (if (< (magicl:tref p 1 0) (+ soil-boundary sloped-height))
                       (cl-mpm/setup::plane-point-sdf
                        (magicl:from-list (list (magicl:tref p 0 0)
                                                (magicl:tref p 1 0)) '(2 1))
                        normal
                        (magicl:from-list (list shelf-length soil-boundary)
                                          '(2 1) :type 'double-float))

                       (cl-mpm/setup::plane-point-sdf
                        (magicl:from-list (list (magicl:tref p 0 0)
                                                (magicl:tref p 1 0)) '(2 1))
                        (magicl:from-list (list 1d0 0d0) '(2 1)  :type 'double-float)
                        sloped-inflex-point)
                       )
                   1d0)))
        (dotimes (i edge-refine)
          (dolist (dir (list :x :y))
            (cl-mpm::split-mps-criteria
             *sim*
             (lambda (mp h)
               (when (and 
                      (or
                       (<= (cutout
                            (cl-mpm/fastmaths:fast-.+
                             (cl-mpm/particle:mp-position mp)
                             (cl-mpm/fastmaths:fast-.*
                              (cl-mpm/particle::mp-domain-size mp)
                              (cl-mpm/utils:vector-from-list (list 0.5d0 0.5d0 0d0))))
                            ) (* mesh-size 0d0))
                       (<= (cutout
                              (cl-mpm/fastmaths:fast-.+
                               (cl-mpm/particle:mp-position mp)
                               (cl-mpm/fastmaths:fast-.*
                                (cl-mpm/particle::mp-domain-size mp)
                                (cl-mpm/utils:vector-from-list (list 0.5d0 0d0 0d0))))
                              ) (* mesh-size 0d0))
                       (<= (cutout
                            (cl-mpm/fastmaths:fast-.+
                             (cl-mpm/particle:mp-position mp)
                             (cl-mpm/fastmaths:fast-.*
                              (cl-mpm/particle::mp-domain-size mp)
                              (cl-mpm/utils:vector-from-list (list 0.5d0 -0.5d0 0d0))))
                            ) (* mesh-size 0d0))
                       )
                          (> (cutout (cl-mpm/particle:mp-position mp)) (- (* mesh-size 1d0))))
                 dir)))))
        (cl-mpm/setup::remove-sdf *sim*
                                  #'cutout
                                  )
        )
      (let* ( 
             (notched-depth 0.0d0)
             ;(notched-depth 0.5d0) 
           ;; (undercut-angle 45d0)
           (undercut-angle 45d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list
                                (- shelf-length notched-depth)
                                soil-boundary)
                               '(2 1) :type 'double-float)))

      (flet ((cutout (p)
               (if (and
                    (> (magicl:tref p 1 0) soil-boundary))
                   (cl-mpm/setup::plane-point-sdf
                    (magicl:from-list (list (magicl:tref p 0 0)
                                            (magicl:tref p 1 0)) '(2 1))
                    normal
                    sloped-inflex-point)
                   1d0)))
        (when (> notched-depth 0d0)
          (dotimes (i edge-refine)
            (dolist (dir (list :x :y))
              (cl-mpm::split-mps-criteria
               *sim*
               (lambda (mp h)
                 (when (and
                        (or
                         (<= (cutout
                              (cl-mpm/fastmaths:fast-.+
                               (cl-mpm/particle:mp-position mp)
                               (cl-mpm/fastmaths:fast-.*
                                (cl-mpm/particle::mp-domain-size mp)
                                (cl-mpm/utils:vector-from-list (list 0.5d0 0.5d0 0d0))))
                              ) (* mesh-size 0d0))
                         (<= (cutout
                              (cl-mpm/fastmaths:fast-.+
                               (cl-mpm/particle:mp-position mp)
                               (cl-mpm/fastmaths:fast-.*
                                (cl-mpm/particle::mp-domain-size mp)
                                (cl-mpm/utils:vector-from-list (list 0.5d0 0d0 0d0))))
                              ) (* mesh-size 0d0))
                         (<= (cutout
                              (cl-mpm/fastmaths:fast-.+
                               (cl-mpm/particle:mp-position mp)
                               (cl-mpm/fastmaths:fast-.*
                                (cl-mpm/particle::mp-domain-size mp)
                                (cl-mpm/utils:vector-from-list (list 0.5d0 -0.5d0 0d0))))
                              ) (* mesh-size 0d0))
                         )
                        (> (cutout (cl-mpm/particle:mp-position mp)) (- (* mesh-size 1d0))))
                   dir)))))
          (cl-mpm/setup:remove-sdf *sim*
                                   #'cutout
                                   ))))

      (when nil
		(let ((cut-height (* 0.5d0 shelf-height-true))
              (cut-back-distance 0.15d0)
              (width (* 0.5d0 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
          (cl-mpm/setup::apply-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
                                                      (cl-mpm/utils:vector-from-list (list (magicl:tref p 0 0)
                                                                              (magicl:tref p 1 0)
                                                                              0d0))
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            (float shelf-height 0d0)
                                                            0d0)
                                                      (list (- (magicl:tref sloped-inflex-point 0 0)
                                                               (* cut-back-distance shelf-height-true))
                                                            (float (- shelf-height cut-height) 0d0)
                                                            0d0)
                                                      width))
                                   (lambda (mp v)
                                     (let ((d 
                                             ;(* 0.99d0 (exp (- (expt (/ (+ width v) width) 1))))
										      (* 0.99d0 (cl-mpm/damage::weight-func (expt (+ width v) 2) width))
                                             ))
                                       (setf (cl-mpm/particle:mp-damage mp)
                                             d)
                                       (let ((k (cl-mpm/damage::find-k-damage-mp mp d)))
                                         (setf (cl-mpm/particle::mp-history-stress mp)
                                               k)))
                                     (cl-mpm/damage::update-damage mp 1d-3)
                                     ))
          )
      ))
  
     (setf cl-mpm::*max-split-depth* 4)

     ;; (let ((ratio 1.0d0))
     ;;   (cl-mpm/setup::damage-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
     ;;                                                (magicl:from-list (list (magicl:tref p 0 0)
     ;;                                                                        (magicl:tref p 1 0)) '(2 1))
     ;;                                                (list (- shelf-length (* shelf-height ratio)) shelf-height)
     ;;                                                (list shelf-length soil-boundary)
     ;;                                                10d0
     ;;                                                )) 1.0d0))
     )
    ;; (let ((upper-random-bound 0.5d0))
    ;;   (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;         do (setf (cl-mpm/particle::mp-damage mp)
    ;;                  (reduce #'*
    ;;                          (loop for i from 0 to 2
    ;;                                collect (random upper-random-bound))))))
    (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
    ;(loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* *output-directory* )) do (uiop:delete-file-if-exists f))
    (defparameter *run-sim* t)
    (defparameter *t* 0)
    (defparameter *oobf* 0)
    (defparameter *energy* 0)
    (defparameter *sim-step* 0))


(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)

  ;; (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
  ;;                                           *sim*
  ;;                                           (list :dt target-time))

  (let* ((target-time 1d2)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1.0d0))

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
                       (setf (cl-mpm::sim-enable-damage *sim*) t)
                       (let ((ms (cl-mpm::sim-mass-scale *sim*)))
                        (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-4 ms))))
                     (format t "Step ~d ~%" steps)
                     (format t "MPs ~d ~%" (length (cl-mpm:sim-mps *sim*)))
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

(defparameter *balance-point* 1.5d0)
(defun mpi-loop ()
  (format t "Starting mpi~%")
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (setup :undercut 0d0)
    ;(setup)


    (when (typep *sim* 'cl-mpm/mpi::mpm-sim-mpi)
      ;;Square
      (let ((height 1))
        ;(when (> (cl-mpi:mpi-comm-size) 8)
        ;  (setf height 2))
        (let* ( (dsize (ceiling (cl-mpi:mpi-comm-size) height)))
          (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize height 1))))
      ;;Setup domain deomposition
      (setf cl-mpm/mpi::*prune-nodes* nil)
      ;(setf cl-mpm/mpi::*prune-nodes* t)
      (when (= rank 0)
        (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
        (format t "Decompose~%"))

     (let ((mp (aref (cl-mpm:sim-mps *sim*) 0)))
        (when (slot-exists-p mp 'cl-mpm/particle::local-length)
          (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
            ;(when (= rank 0)
            ;  (format t "Min size ~A length scale ~F~%" (mapcar (lambda (x) (abs (reduce #'- x)))  (cl-mpm/mpi::mpm-sim-mpi-domain-bounds *sim*)) dhalo-size) )
            (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))))

      (cl-mpm/mpi::setup-domain-bounds *sim*)


	  (let ((balance 1.5d0)
            (stag nil))
        (loop repeat 1000
              while (and (if balance (> balance 1.10d0) t)
                         (not stag)
                         )
              do (multiple-value-bind (bal stagnent)
                     (cl-mpm/mpi::load-balance *sim*
                                               :exchange-mps nil
                                               :step-size 1d-2
                                               ;:dims (list :x :y)
                                               )
                   (setf balance bal
                         stag stagnent)))
        (when stag
          (when (= rank 0)
            (format t "Stagnated ~%")))
        (defparameter *balance-point* balance))
      

      (cl-mpm/mpi::domain-decompose *sim*) 

      (let ((mp (aref (cl-mpm:sim-mps *sim*) 0)))
        (when (slot-exists-p mp 'cl-mpm/particle::local-length)
          (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
            (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size)))))

    (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps *sim*)))
    (when (= rank 0)
      (format t "Run mpi~%"))
    (run-mpi)
    (when (= rank 0)
      (format t "Done mpi~%"))
    )
  )

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
      (time
        (progn
          ,@body))
      (progn
        ,@body)))

(defun run-mpi ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames *output-directory* "mesh.vtk") *sim*)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (target-time 1d1)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (accelerate-target-time 1d0)
         (accelerate-mass-scale 1d4)
         (collapse-target-time 0.1d0)
         (collapse-mass-scale 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (settle-steps 0)
         (damp-steps 0)
         (sim-state :settle)
         (dt-scale 0.5d0)
         (dt-0 0d0)
         (plasticity-enabled t)
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
         (damping-0
            (* 1d-4
               (cl-mpm/setup::estimate-critical-damping *sim*)))
		 (damage-0
		    (cl-mpm/mpi:mpi-sum
		 	   (lparallel:pmap-reduce (lambda (mp)
		 		   (*
		 			 1d0
		 			 (cl-mpm/particle::mp-mass mp)
		 			 (cl-mpm/particle::mp-damage mp)))
		 		   #'+ (cl-mpm:sim-mps *sim*)
		 		   :initial-value 0d0)))
         (max-step 5000)
         (criteria-energy 1d-2)
         (criteria-oobf 1d-1)
         )

    (defparameter *data-damage* 0d0)
    (defparameter *data-energy* 0d0)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))
    (setf (cl-mpm:sim-damping-factor *sim*) 
            (* 0.1d0 
               (cl-mpm/setup::estimate-critical-damping *sim*)))

    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :dt-scale dt-scale
     :energy-crit 1d-2
     :oobf-crit 1d-1
     :substeps 50
     :conv-steps 5000
     :dt-scale dt-scale
     :post-iter-step
     (lambda (i e o)))
    (cl-mpm:iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
       (cl-mpm/fastmaths::fast-zero (cl-mpm/particle::mp-acceleration mp))))


    (let ((ms accelerate-mass-scale))
      (setf (cl-mpm::sim-mass-scale *sim*) ms) 
      (setf target-time accelerate-target-time)
      (setf (cl-mpm:sim-damping-factor *sim*) 
            (* 0.1d0 
               (cl-mpm/setup::estimate-critical-damping *sim*))))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))

    (when t
      (setf sim-state :settle)
      (setf (cl-mpm:sim-damping-factor *sim*)
            damping-0))
    (when t
      (setf (cl-mpm::sim-enable-damage *sim*) t)
      (cl-mpm::iterate-over-mps
       (cl-mpm:sim-mps *sim*)
       (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled))))

    (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
        (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) 
              (round substeps 1)))

    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))

    (when (= rank 0)
      (format t "Substeps ~D~%" substeps)
 	    (cl-mpm/output::save-simulation-parameters (merge-pathnames *output-directory* "settings.json")
                                             *sim*
                                             (list :dt target-time
                                                   :criteria-energy criteria-energy
                                                   :criteria-oobf criteria-oobf
                                                   ))
		(with-open-file (stream (merge-pathnames *output-directory* "timesteps.csv") :direction :output :if-exists :supersede)
			  (format stream "steps,time,damage,energy,oobf,step-type~%")))
    (let ((work 0d0))
        (time (loop for steps from 0 to max-step
                while *run-sim*
                do
                   (progn

                    (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
                        (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) (ceiling (* substeps 0.5))))
                     (when (= rank 0)
                       (format t "Step ~d/~D ~%" steps max-step))
                     (when (typep *sim* 'cl-mpm/mpi::mpm-sim-mpi)
                         (cl-mpm/mpi::load-balance-algo *sim*
                                                        :step-size 1d-2
                                                        ;:min-bounds 1.1d0
                                                        ;:max-bounds 1.5d0
                                                        :min-bounds *balance-point* 
                                                        :max-bounds (* 2d0 *balance-point* )
                                                        ))
                     (when (= (mod steps 1) 0)
                       (cl-mpm/output:save-vtk (merge-pathnames *output-directory* (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                       ;(cl-mpm/output::save-vtk-nodes (merge-pathnames *output-directory* (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                       )
					 (when (= rank 0)
						(with-open-file (stream (merge-pathnames *output-directory* "timesteps.csv") :direction :output :if-exists :append)
                                               (format stream "~D,~f,~f,~f,~f,~A~%"
													   steps
													   *t*
                                                       *data-damage*
                                                       *data-energy*
                                                       *oobf*
                                                       sim-state)))
                     (let ((energy-estimate 0d0)
                           ;(work 0d0)
                           (oobf 0d0)
                           )
                       (rank-0-time
                        rank
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                          (incf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                          (incf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                       ;;This is the proper way
                       ;(setf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                       ;(setf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                       (setf energy-estimate (/ energy-estimate substeps) 
                             oobf (/ oobf substeps))
                       (if (= work 0d0)
                           (setf energy-estimate 0d0)  
                           (setf energy-estimate (abs (/ energy-estimate work))))

                       (setf *oobf* oobf)
;
                       (setf *data-energy* energy-estimate)
                       (let ((damage-est
                               (- (cl-mpm/mpi:mpi-sum
                                    (lparallel:pmap-reduce (lambda (mp)
                                                             (* 
                                                               1d0
					                                           (cl-mpm/particle::mp-mass mp)
                                                               (cl-mpm/particle::mp-damage mp)))
                                                           #'+ (cl-mpm:sim-mps *sim*)
                                                           :initial-value 0d0))
                                  damage-0)))
                         (setf *data-damage* damage-est)
                         (when (= rank 0)
                           (format t "Total damage: ~E~%" damage-est)))

                       (when (= rank 0)
                         (format t "Energy estimate: ~E~%" energy-estimate)
                         (format t "OOBF estimate: ~E~%" *oobf*)
                         (format t "Work estimate ~E~%" work)
                         )

                       (when (= steps damp-steps)
                         (setf sim-state :settle)
                         (setf (cl-mpm:sim-damping-factor *sim*)
                               damping-0))
                       (when (= steps settle-steps)
                         (setf (cl-mpm::sim-enable-damage *sim*) t)
                         (cl-mpm::iterate-over-mps
                          (cl-mpm:sim-mps *sim*)
                          (lambda (mp) (setf (cl-mpm/particle::mp-enable-plasticity mp) plasticity-enabled))))
                       (when (>= steps settle-steps)
                         (if (or
                              ;; t
                              (> energy-estimate criteria-energy)
                              (> *oobf* criteria-oobf)
                              )
                             (when (not (eq sim-state :collapse))
                               (setf sim-state :collapse)
                               (setf work 0d0)
                               (when (= rank 0)
                                 (format t "Changed to collapse~%")))
                             (progn
                               (when (not (eq sim-state :accelerate))
                                 (when (= rank 0)
                                   (format t "Changed to accelerate~%"))
                                 (setf work 0d0)
                                 (setf sim-state :accelerate)
                                 (cl-mpm:iterate-over-mps
                                  (cl-mpm:sim-mps *sim*)
                                  (lambda (mp)
                                    (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp)))))))
                         (case sim-state
                           (:accelerate
                            (when (= rank 0)
                              (format t "Accelerate timestep~%"))
                            (setf
                             target-time accelerate-target-time
                             (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale))
                           (:collapse
                            (when (= rank 0)
                              (format t "Collapse timestep~%"))
                            (setf
                             target-time collapse-target-time
                             (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale)))))

                     ;(let* ((dt-e (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
                     ;       (substep-e (floor target-time dt-e)))
                     ;     (when (= rank 0)
                     ;       (format t "CFL dt estimate: ~f~%" dt-e)
                     ;       (format t "CFL step count estimate: ~D~%" substeps-e))
                     ;     (setf substeps substeps-e
                     ;           (cl-mpm:sim-dt *sim*) dt-e
                     ;           ));)
                     ;(when (not (= )))
                     ;(multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                     ; (when (= rank 0)
                     ;   (format t "CFL dt estimate: ~f~%" dt-e)
                     ;   (format t "CFL step count estimate: ~D~%" substeps-e))
                     ; (setf substeps substeps-e))

					 (let* (;(dt-est (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
					    	(dt-est (* dt-0 (sqrt (cl-mpm::sim-mass-scale *sim*))))
					 	    (substeps-est (floor target-time dt-est)))
					    (when t;(< substeps-est substeps)
					 	   (setf (cl-mpm:sim-dt *sim*) dt-est)
					 	   (setf substeps substeps-est)))

                     (when (= rank 0)
                        (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
                        (format t "CFL step count estimate: ~D~%" substeps))
                     (incf *sim-step*)))))))

(defparameter *output-directory* (merge-pathnames "/nobackup/rmvn14/ham-chalk-conv/output-0.5/"))
(format t "Outputting to ~A~%" *output-directory*)
;(defparameter *output-directory* (merge-pathnames "./output/"))
(ensure-directories-exist *output-directory*)
(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "1"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
;(defparameter *run-sim* nil)
;(setup)
;(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;(run)

(format t "Running~%")
;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
(mpi-loop)
