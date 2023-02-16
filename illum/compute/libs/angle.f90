MODULE ANGLE_M

  USE MATH_M, ONLY: PI
  IMPLICIT NONE

CONTAINS

  REAL(8) FUNCTION ZENITH_ANGLE(point1, point2) RESULT(angle)

    REAL(8), DIMENSION(3), INTENT(IN) :: point1, point2

    angle = ATAN2(NORM2(point2(:2) - point1(:2)), &
                  point2(3) - point1(3))

  END FUNCTION ZENITH_ANGLE

  REAL(8) FUNCTION AZIMUTH_ANGLE(point1, point2) RESULT(angle)

    REAL(8), DIMENSION(3), INTENT(IN) :: point1, point2

    ! ATAN2 is in the (-PI;PI) range, MODULO takes it into (0;2PI)
    angle = MODULO(ATAN2(point2(2) - point1(2), &
                         point2(1) - point1(1)), 2 * PI)

  END FUNCTION AZIMUTH_ANGLE

  REAL(8) FUNCTION ANGLE3PT(point1, point2, point3) RESULT(angle)
    REAL(8), DIMENSION(3), INTENT(IN) :: point1, point2, point3
    REAL(8) :: arg

    REAL(8), DIMENSION(3) :: u, v

    u = point2 - point1
    v = point3 - point2

    arg = SUM(u * v) / (NORM2(u) * NORM2(v))

    IF (arg >= 1) THEN
      angle = 0
    ELSE IF (arg <= -1) THEN
      angle = PI
    ELSE
      angle = ACOS(arg)
    END IF

  END FUNCTION ANGLE3PT

END MODULE ANGLE_M
