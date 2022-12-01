FUNCTION polynomial(x, c, n)

  IMPLICIT NONE

!f2py intent(hide) n
  REAL(8), INTENT(in) :: x, c(n)
  REAL(8) :: polynomial

  INTEGER(4) :: n, i

  polynomial = 0.

  DO i = 1, n
    polynomial = polynomial + c(i) * x**(i - 1)
  END DO

END FUNCTION
