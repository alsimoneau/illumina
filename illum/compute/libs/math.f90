MODULE MATH_M

  IMPLICIT NONE

  REAL(8), PARAMETER :: PI = 4 * ATAN(1.0D0)

CONTAINS

  REAL(8) FUNCTION POLYNOMIAL(val, coeffs) RESULT(p)

    REAL(8), INTENT(IN) :: val, coeffs(:)
    INTEGER(4) :: i

    p = 0

    DO i = 1, SIZE(coeffs)
      p = p + coeffs(i) * val**(i - 1)
    END DO

  END FUNCTION POLYNOMIAL

  INTEGER FUNCTION RAD2RANK(rad, N_angle) RESULT(rank)

    REAL(8), INTENT(IN) :: rad
    INTEGER(4), INTENT(IN) :: N_angle

    rank = INT(rad / PI * (N_angle - 1)) + 1

  END FUNCTION RAD2RANK

END MODULE MATH_M
