MODULE MATH_M

  IMPLICIT NONE

  REAL(8), PARAMETER :: PI = 4 * ATAN(1.0D0)

CONTAINS

  REAL(8) FUNCTION DEG2RAD(deg) RESULT(rad)

    REAL(8), INTENT(IN) :: deg

    rad = deg * PI / 180

  END FUNCTION DEG2RAD

  REAL(8) FUNCTION RAD2DEG(rad) RESULT(deg)

    REAL(8), INTENT(IN) :: rad

    deg = rad * 180 / PI

  END FUNCTION RAD2DEG

  REAL(8) FUNCTION POLYNOMIAL(val, coeffs) RESULT(p)

    REAL(8), INTENT(IN) :: val, coeffs(:)
    INTEGER(4) :: i

    p = 0

    DO i = 1, SIZE(coeffs)
      p = p + coeffs(i) * val**(i - 1)
    END DO

  END FUNCTION POLYNOMIAL

END MODULE MATH_M
