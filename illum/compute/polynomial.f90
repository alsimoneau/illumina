FUNCTION POLYNOMIAL(x, c, n) RESULT(y)

  IMPLICIT NONE

!f2py intent(hide) n
  REAL(8), INTENT(in) :: x, c(n)
  REAL(8) :: y

  INTEGER(4) :: n, i

  y = 0.

  DO i = 1, n
    y = y + c(i) * x**(i - 1)
  END DO

END FUNCTION
