(module bl
        #:export (
                  ensemble-create
                  ensemble-push
                  ensemble-create-space
                  ensemble-get-num-ptcls
                  ensemble-set-num-ptcls
                  ensemble-get-component
                  ensemble-set-component
                  ensemble-get-phase-space
                  ensemble-get-internal-state
                  ))

(load-extension "libbl_guile" "init_ensemble")


(define (ensemble-get-component component ensemble)
  (let ((a
         (make-typed-array 'f64 0
                           (ensemble-get-num-ptcls ensemble))))
    ensemble-get-comp ensemble a))

