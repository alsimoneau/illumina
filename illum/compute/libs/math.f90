MODULE MATH_M

  IMPLICIT NONE
  PRIVATE
  PUBLIC PI, POINT, DEG2RAD, RAD2DEG, POLYNOMIAL

  REAL(8), PARAMETER :: PI = 4 * ATAN(1.0D0)

  TYPE POINT
    REAL(8), DIMENSION(3) :: arr = [0, 0, 0]

  CONTAINS
    PROCEDURE, PASS(self) :: X, Y, Z
    PROCEDURE, PRIVATE, PASS(self) :: PT_ADD_PT, PT_SUB_PT
    PROCEDURE, PRIVATE, PASS(self) :: PT_MUL_PT, PT_MUL_REAL, PT_DIV_REAL
    GENERIC :: OPERATOR(+) => PT_ADD_PT
    GENERIC :: OPERATOR(-) => PT_SUB_PT
    GENERIC :: OPERATOR(*) => PT_MUL_PT, PT_MUL_REAL
    GENERIC :: OPERATOR(/) => PT_DIV_REAL
  END TYPE

CONTAINS

  FUNCTION X(self)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8) :: x
    x = self % arr(1)
  END FUNCTION X

  FUNCTION Y(self)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8) :: y
    y = self % arr(2)
  END FUNCTION Y

  FUNCTION Z(self)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8) :: z
    z = self % arr(3)
  END FUNCTION Z

  FUNCTION PT_ADD_PT(self, pt) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    TYPE(POINT), INTENT(IN) :: pt
    TYPE(POINT) :: res

    res % arr = self % arr + pt % arr

  END FUNCTION PT_ADD_PT

  FUNCTION PT_SUB_PT(self, pt) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    TYPE(POINT), INTENT(IN) :: pt
    TYPE(POINT) :: res

    res % arr = self % arr - pt % arr

  END FUNCTION PT_SUB_PT

  FUNCTION PT_MUL_PT(self, pt) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    TYPE(POINT), INTENT(IN) :: pt
    REAL(8) :: res

    res = SUM(self % arr * pt % arr)

  END FUNCTION PT_MUL_PT

  FUNCTION PT_MUL_REAL(self, val) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8), INTENT(IN) :: val
    TYPE(POINT) :: res

    res % arr = self % arr * val

  END FUNCTION PT_MUL_REAL

  FUNCTION PT_DIV_REAL(self, val) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8), INTENT(IN) :: val
    TYPE(POINT) :: res

    res % arr = self % arr / val

  END FUNCTION PT_DIV_REAL

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

  FUNCTION POLYNOMIAL(val, coeffs) RESULT(p)

    REAL(8), INTENT(IN) :: val, coeffs(:)
    REAL(8) :: p
    INTEGER(4) :: i

    p = 0

    DO i = 1, SIZE(coeffs)
      p = p + coeffs(i) * val**(i - 1)
    END DO

  END FUNCTION POLYNOMIAL

END MODULE MATH_M
