MODULE MATH_M

  IMPLICIT NONE

  REAL(8), PARAMETER :: PI = 4 * ATAN(1.0D0)

CONTAINS

  FUNCTION DEG2RAD(deg) RESULT(rad)

    REAL(8), INTENT(IN) :: deg
    REAL(8) :: rad

    rad = deg * PI / 180

  END FUNCTION DEG2RAD

  FUNCTION RAD2DEG(rad) RESULT(deg)

    REAL(8), INTENT(IN) :: rad
    REAL(8) :: deg

    deg = rad * 180 / PI

  END FUNCTION RAD2DEG

  FUNCTION POLYNOMIAL(x, coeffs) RESULT(y)

    REAL(8), INTENT(IN) :: x, coeffs(:)
    REAL(8) :: y
    INTEGER(4) :: i

    y = 0

    DO i = 1, SIZE(coeffs)
      y = y + coeffs(i) * x**(i - 1)
    END DO

  END FUNCTION

END MODULE MATH_M
