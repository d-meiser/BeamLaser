(define-module (test-ensemble)
  #:use-module (srfi srfi-64)
  #:use-module (beam-laser)
  )

(test-begin "ensemble")

(test-assert "Can be created"
  (ensemble-create "e1" 20 4))

(test-assert "Set number of particles"
  (let ((e (ensemble-create "e1" 20 4)))
    (begin
      (ensemble-set-num-ptcls 5 e)
      (eq? 5 (ensemble-get-num-ptcls e)))))

(test-assert "Can create space"
  (let ((e (ensemble-create "e1" 20 4)))
    (ensemble-create-space 5 e)))

(test-assert "Can create more space than original capacity"
  (let ((e (ensemble-create "e1" 20 4)))
    (ensemble-create-space 25 e)))

(test-assert "Get component"
  (let ((e (ensemble-create "e1" 20 4)))
    (begin
      (ensemble-set-num-ptcls 4 e)
      (let
          ((comp (ensemble-get-component 0 e)))
        (and
         (eq? 4 (array-length comp))
         (eq? 1 (array-rank comp)))
      ))))

(test-end "ensemble")

(exit (= (test-runner-fail-count (test-runner-current)) 0))
