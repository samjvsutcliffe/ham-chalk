(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/joss)

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;(cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  )


(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
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
                     ) mp
      (declare (double-float pressure damage))
        (progn
          (when (< damage 1d0)
            (let ((cauchy-undamaged (magicl:scale stress (/ 1d0 (magicl:det def)))))
              (multiple-value-bind (s_1 s_2 s_3) (cl-mpm/damage::principal-stresses-3d stress)
                (multiple-value-bind (e_1 e_2 e_3) (cl-mpm/damage::principal-stresses-3d strain)
                  (let* ((pressure-effective (* 1d0 damage pressure))
                         ;(j2 (cl-mpm/constitutive::voigt-j2
                         ;      (cl-mpm/utils::deviatoric-voigt cauchy-undamaged)))
                         ;(p (+ s_1 s_2 s_3))
                         ;(angle (* angle (/ pi 180d0)))

                         ;;Good drucker-prager
                         ;(max 
                         ;  (s_1 (* (/ 3d0 (+ 3 (tan angle)))
                         ;          (+ (sqrt (* 3 j2)) (* 1/3 (tan angle) p))))
                         ;  (max p 0d0))
                       (et (sqrt (max 0d0
                                      (* E
                                         (+
                                          (* (max s_1 0d0) e_1)
                                          (* (max s_2 0d0) e_2)
                                          (* (max s_3 0d0) e_3)
                                          )))))
                       (ec (sqrt (max 0d0
                                      (* E
                                         (+
                                          (* (min s_1 0d0) e_1)
                                          (* (min s_2 0d0) e_2)
                                          (* (min s_3 0d0) e_3)
                                          )))))
                       (k (/ fc ft))
                       (s_1
                         et
                         ;(/
                         ; (+ (* k et) ec)
                         ; (+ k 1d0)
                         ; )
                         )

                         ;(k (/ fc ft))
                         ;(i1 (+ s_1 s_2 s_3))
                         ;(k-factor (/ (- k 1d0)
                         ;            (- 1d0 (* 2d0 nu))))
                         ;(s_1 (+ (* i1 (/ k-factor (* 2d0 k)))
                         ;       (* (/ 1d0 (* 2d0 k))
                         ;          (sqrt (+ (expt (* k-factor i1) 2)
                         ;                   (* (/ (* 12 k) (expt (- 1d0 nu) 2))j2)
                         ;                   )))))
                         )
                    (when (> s_1 0d0)
                      (setf damage-increment s_1)))))))
          (when (>= damage 1d0)
            (setf damage-increment 0d0))
          ;;Delocalisation switch
          (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
          (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)))))


(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes-damage
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
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                'cl-mpm/particle::particle-chalk-delayed
                :E 1d9
                :nu 0.24d0
                :enable-plasticity t

                ;:ft 0.70d6
                ;:fc 4.0d6
                :ft 10d3
                :fc 400d3 
                :friction-angle 60d0

                :kt-res-ratio 0d-15
                :kc-res-ratio 1d-2
                :g-res-ratio 6.5d-3
                ;:g-res-ratio 5.0d-3

                :fracture-energy 3000d0

                :initiation-stress 20d3
                :delay-time 1d0
                :ductility 4d0

                :critical-damage 1d0;(- 1.0d0 1d-3)
                :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage
                :local-length 0.5d0
                :local-length-damaged 10d-10

                :psi (* 0d0 (/ pi 180))
                :phi (* 42d0 (/ pi 180))
                :c 121d3

                :gravity -9.8d0
                :gravity-axis (cl-mpm/utils:vector-from-list (list 0d0 1d0 0d0))
                ))))
      (let* ((mp-0 (aref (cl-mpm:sim-mps sim) 0))
             (fc (cl-mpm/particle::mp-fc mp-0))
             (ft (cl-mpm/particle::mp-ft mp-0))
             (angle-d (* (/ 180 pi) (atan (* 3 (/ (- fc ft) (+ fc ft))))))
             (rc (cl-mpm/particle::mp-k-compressive-residual-ratio mp-0))
             (rs (cl-mpm/particle::mp-shear-residual-ratio mp-0))
             (angle-plastic (cl-mpm/particle::mp-phi mp-0))
             (angle-plastic-damaged (atan (* (/ rs rc) (tan angle-plastic))))
             )
        (format t "Chalk damage growth angle: ~F~%"
                angle-d)
        (format t "Chalk plastic virgin angle: ~F~%"
                (* (/ 180 pi) angle-plastic))
        (format t "Chalk plastic residual angle: ~F~%"
                (* (/ 180 pi) angle-plastic-damaged)))
      ;; (let* ((mp-0 (aref (cl-mpm:sim-mps *sim*) 0))
      ;;        (fc (cl-mpm/particle::mp-fc mp-0))
      ;;        )
      ;;   (format t "Chalk damage angle: ~F~%"
      ;;           (atan (* 3 (/ (- fc ft) (+ fc ft))))))
      ;; (cl-mpm/examples/tpb::calculate-ductility-param 1d9 200d0 1d0 200d3)
      (setf (cl-mpm:sim-allow-mp-split sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-15)
      (let ((ms 1d6))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 ms))
        ;; (setf (cl-mpm:sim-damping-factor sim) 10.0d0)
        ;; (setf (cl-mpm:sim-damping-factor sim) (* 1d-1 ms))
        )

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

      (let ((dt-scale 0.8d0))
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
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
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
      sim)))




(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (c cl-mpm/particle::mp-c))
      mp
    ;; (let ((rho_0 200d3)
    ;;       (rho_1 0.002d3)
    ;;       (soft 1d0))
    ;;   (setf c (max rho_1
    ;;                  (* rho_0 (exp (- (* soft ps)))))))
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

(defun setup ()
  (let* ((mesh-size 0.5)
         (mps-per-cell 2)
         (shelf-height 15.5)
         (soil-boundary 2)
         (shelf-aspect 1.25)
         (runout-aspect 1.25)
         (shelf-length (* shelf-height shelf-aspect))
         (domain-length (+ shelf-length (* runout-aspect shelf-height)))
         (shelf-height-true shelf-height)
         (shelf-height (+ shelf-height soil-boundary))
         (depth 400)
         (offset (list 0 (* 0 mesh-size)
                       ;; 0
                       ))
         )
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
                      ;; 3d0
                      (* 0.5d0 shelf-height)
                      ))
                (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                   shelf-length)
                (> (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                   soil-boundary
                   )
                ;; (< (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)
                ;;    (- shelf-height 13d0)
                ;;    )
                )
             dir
             )))))
    (loop for mp across (cl-mpm::sim-mps *sim*)
          do (setf (cl-mpm/particle::mp-split-depth mp) 0))

    (let* ((sloped-height (- (- shelf-height soil-boundary) 7d0))
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
           )
      (cl-mpm/setup::remove-sdf *sim*
                               (lambda (p)
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
                                     1d0)
                                 ))
      (let ((cut-height (* 0.5d0 shelf-height-true))
            (width mesh-size))
        (cl-mpm/setup::apply-sdf *sim* (lambda (p) (cl-mpm/setup::line-sdf
                                                    (magicl:from-list (list (magicl:tref p 0 0)
                                                                            (magicl:tref p 1 0)) '(2 1))
                                                    (list (- (magicl:tref sloped-inflex-point 0 0)
                                                             (* 0.15d0 shelf-height-true))
                                                          shelf-height)
                                                    (list (- (magicl:tref sloped-inflex-point 0 0)
                                                             (* 0.15d0 shelf-height-true))
                                                          (- shelf-height cut-height))
                                                    width
                                                    ))
                                 (lambda (mp v)
                                   (setf (cl-mpm/particle:mp-damage mp)
                                         (exp (- (expt (+ width v) 2)))
                                         )
                                   ))
        ;; (cl-mpm/setup::damage-sdf
        ;;  *sim*
        ;;  (lambda (p)
        ;;    (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
        ;;        (funcall
        ;;         (cl-mpm/setup::rectangle-sdf
        ;;          (list (- (magicl:tref sloped-inflex-point 0 0)
        ;;                   (* 0.15d0 shelf-height-true))
        ;;                shelf-height)
        ;;          (list (* 1.0d0 mesh-size) cut-height)
        ;;          ) p)
        ;;        0.99d0)
        ;;    ))
        ))
    (let* ((notched-depth 1.0d0)
           (undercut-angle 30d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))) '(2 1)))
           (sloped-inflex-point
             (magicl:from-list (list
                                (- shelf-length notched-depth)
                                soil-boundary)
                               '(2 1) :type 'double-float)

             )
           )
      (cl-mpm/setup:remove-sdf *sim*
                               (lambda (p)
                                 (if (and
                                      (> (magicl:tref p 1 0) soil-boundary))
                                     (cl-mpm/setup::plane-point-sdf
                                      (magicl:from-list (list (magicl:tref p 0 0)
                                                              (magicl:tref p 1 0)) '(2 1))
                                      normal
                                      sloped-inflex-point)
                                     1d0)
                                 ))
      )

   
     (let ((notch-height 0.0d0))
       (cl-mpm/setup:remove-sdf *sim*
                                (lambda (p)
                                  (if t ;(< (abs (- (magicl:tref p 2 0) (* 0.5d0 depth))) 20d0)
                                      (funcall
                                       (cl-mpm/setup::rectangle-sdf
                                        (list shelf-length (+ notch-height soil-boundary))
                                        (list (* 2d0 notch-height) notch-height)
                                        ) p)
                                      1d0)
                                  )))
    
     (setf cl-mpm::*max-split-depth* 6)

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
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./output/")) do (uiop:delete-file-if-exists f))
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
                        (setf (cl-mpm:sim-damping-factor *sim*) (* 1d-4 ms)))
                       )
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

(defun mpi-loop ()
  (format t "Starting mpi~%")
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (setup :undercut 0d0)
    ;(setup)

    ;;Square
    (let ((dsize (floor (cl-mpi:mpi-comm-size))))
      (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize 1 1)))
    ;;Setup domain deomposition
    ;(let ((dsize (floor (cl-mpi:mpi-comm-size) 2)))
    ;  (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize 2 1)))

    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Decompose~%"))
    ;(cl-mpm/mpi::domain-decompose *sim*)
	(let* ((density-ratio 0.75d0)
			  (inflection 0.60d0)
			  (di (* density-ratio inflection))
			  )
		  (cl-mpm/mpi::domain-decompose *sim* :domain-scaler (lambda (domain)
															   (format t "~A~%" domain)
															   (destructuring-bind (x y z) domain
																 (list (mapcar
																		(lambda (p)
                                                                          (expt p 1.3)
																		;  (if (> p inflection)
																		;	  (+
																		;	   (* (/ (- 1d0 di)
																		;			 (- 1d0 inflection)) p)
																		;	   (-
																		;		di
																		;		(*
																		;		 (/ (- 1d0 di)
																		;			  (- 1d0 inflection))
																		;		 inflection)
																		;		))
																		;	  (* density-ratio p)
																		;			 )
																				 ) x)
																	   y z)))))
    ;(when (= rank 0))
    (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps *sim*)))
    (let ((mp (aref (cl-mpm:sim-mps *sim*) 0)))
      (when (slot-exists-p mp 'cl-mpm/particle::local-length)
        (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
          (when (= rank 0)
            (format t "Min size ~A length scale ~F~%" (mapcar (lambda (x) (abs (reduce #'- x)))  (cl-mpm/mpi::mpm-sim-mpi-domain-bounds *sim*)) dhalo-size) )
          (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))))
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
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (target-time 1d1)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (collapse-target-time 0.1d0)
         (collapse-mass-scale 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (settle-steps 6)
         (damp-steps 5)
         (dt-scale 1.0d0)
         (dt-0 0d0)
         )

    (cl-mpm::update-sim *sim*)
    ;(multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
    ;  (when (= rank 0)
    ;    (format t "CFL dt estimate: ~f~%" dt-e)
    ;    (format t "CFL step count estimate: ~D~%" substeps-e))
    ;  (setf substeps substeps-e))

    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))

    (when (= rank 0)
      (format t "Substeps ~D~%" substeps)
 	    (cl-mpm/output::save-simulation-parameters #p"output/settings.json"
                                             *sim*
                                             (list :dt target-time)) 
		(with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :supersede)
			  (format stream "steps,time~%")) 
        )
    (time (loop for steps from 0 to 200
                while *run-sim*
                do
                   (progn
                     (when (= rank 0)
                       (format t "Step ~d ~%" steps))
                     ;(cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
					 (when (= rank 0)
						(with-open-file (stream (merge-pathnames "output/timesteps.csv") :direction :output :if-exists :append)
											   (format stream "~D,~f~%"
													   steps
													   *t*
													   )))
                     ;(let ((damage-mps (cl-mpm/mpi::mpi-sync-damage-mps *sim* (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*))))
                     ;  (cl-mpm/mpi::save-damage-vtk
                     ;   (merge-pathnames (format nil "output/sim_damage_~2,'0d_~5,'0d.vtk" rank *sim-step*))
                     ;   damage-mps)
                     ;  )
                     (let ((energy-estimate 0d0))
                       (rank-0-time
                        rank
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          (incf energy-estimate (estimate-energy-crit *sim*))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                       (setf energy-estimate (/ energy-estimate substeps))

                       (when (= rank 0)
                         (format t "Energy estimate: ~E~%" energy-estimate))
                       (when (>= steps damp-steps)
                         (setf (cl-mpm:sim-damping-factor *sim*) 0d-2))
                       (when (>= steps settle-steps)
                         (setf (cl-mpm::sim-enable-damage *sim*) t)
                         (if (> energy-estimate 1d1)
                           (progn
                             (when (= rank 0)
                               (format t "Collapse timestep~%"))
                             ;(when (not (= target-time collapse-target-time))
                                 (setf
                                  target-time collapse-target-time
                                  (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale)
                               ;)
                             
                             )
                           (progn
                             (when (= rank 0)
                               (format t "Accelerate timestep~%"))
                             (setf
                              ;dt-scale 0.8d0
                              target-time 1d0
                              (cl-mpm::sim-mass-scale *sim*) 1d2
                              )
                             ;(when (not (= target-time collapse-target-time)) 
                             ;    (setf
                             ;     target-time 1d0
                             ;     (cl-mpm::sim-mass-scale *sim*) 1d4
                             ;     ))
                             
                             )
                           
                           )))
                     ;(let* ((dt-e (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
                     ;       (substep-e (floor target-time dt-e)))
                     ;     (when (= rank 0)
                     ;       (format t "CFL dt estimate: ~f~%" dt-e)
                     ;       (format t "CFL step count estimate: ~D~%" substeps-e))
                     ;     (setf substeps substeps-e
                     ;           (cl-mpm:sim-dt *sim*) dt-e
                     ;           ));)
                     ;(when (not (= )))
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                      (when (= rank 0)
                        (format t "CFL dt estimate: ~f~%" dt-e)
                        (format t "CFL step count estimate: ~D~%" substeps-e))
                      (setf substeps substeps-e))

                     (setf (cl-mpm:sim-dt *sim*) (* dt-0 (sqrt (cl-mpm::sim-mass-scale *sim*))))

                     (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))

                     (when (= rank 0)
                        (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
                        (format t "CFL step count estimate: ~D~%" substeps))
                     (incf *sim-step*)
                     ;(plot *sim*)
                     (swank.live:update-swank)
                     )))
    )
  )

(setf lparallel:*kernel* (lparallel:make-kernel 16 :name "custom-kernel"))
;(defparameter *run-sim* nil)
;(setup)
;(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;(run)

(format t "Running~%")
;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
(mpi-loop)
