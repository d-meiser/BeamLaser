(load-extension "./libbl_guile" "init_ensemble")


(define (ensemble-get-component component ensemble)
  (let ((a
         (make-typed-array 'f64 0
                           (ensemble-get-num-ptcls ensemble))))
    ensemble-get-component component ensemble a))

