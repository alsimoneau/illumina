MODULE MATH_M

  IMPLICIT NONE
  PRIVATE
  PUBLIC PI, POINT, DEG2RAD, RAD2DEG, POLYNOMIAL

  REAL(8), PARAMETER :: PI = 4 * ATAN(1.0D0)

  TYPE POINT
    REAL(8), DIMENSION(3) :: arr = [0, 0, 0]

  CONTAINS
    PROCEDURE, PASS(self) :: x, y, z
    PROCEDURE, PRIVATE, PASS(SELF) :: pt_add_pt, pt_sub_pt
    PROCEDURE, PRIVATE, PASS(self) :: pt_mul_pt, pt_mul_real, pt_div_real
    GENERIC :: OPERATOR(+) => pt_add_pt
    GENERIC :: OPERATOR(-) => pt_sub_pt
    GENERIC :: OPERATOR(*) => pt_mul_pt, pt_mul_real
    GENERIC :: OPERATOR(/) => pt_div_real
  END TYPE

CONTAINS

  FUNCTION x(self)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8) :: x
    x = self % arr(1)
  END FUNCTION x

  FUNCTION y(self)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8) :: y
    y = self % arr(2)
  END FUNCTION y

  FUNCTION z(self)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8) :: z
    z = self % arr(3)
  END FUNCTION z

  FUNCTION pt_add_pt(self, pt) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    TYPE(POINT), INTENT(IN) :: pt
    TYPE(POINT) :: res

    res % arr = self % arr + pt % arr

  END FUNCTION pt_add_pt

  FUNCTION pt_sub_pt(self, pt) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    TYPE(POINT), INTENT(IN) :: pt
    TYPE(POINT) :: res

    res % arr = self % arr - pt % arr

  END FUNCTION pt_sub_pt

  FUNCTION pt_mul_pt(self, pt) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    TYPE(POINT), INTENT(IN) :: pt
    REAL(8) :: res

    res = SUM(self % arr * pt % arr)

  END FUNCTION pt_mul_pt

  FUNCTION pt_mul_real(self, val) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8), INTENT(IN) :: val
    TYPE(POINT) :: res

    res % arr = self % arr * val

  END FUNCTION pt_mul_real

  FUNCTION pt_div_real(self, val) RESULT(res)
    CLASS(POINT), INTENT(IN) :: self
    REAL(8), INTENT(IN) :: val
    TYPE(POINT) :: res

    res % arr = self % arr / val

  END FUNCTION pt_div_real

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
