(define-module (test-ensemble)
  #:use-module (srfi srfi-64)
  #:use-module (beam-laser)
  )

(test-begin "ensemble")

(test-assert "dummy test" #t)

(test-end "ensemble")

(exit (= (test-runner-fail-count (test-runner-current)) 0))
