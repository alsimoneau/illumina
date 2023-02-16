MODULE ANGLE_M

  USE MATH_M
  IMPLICIT NONE

CONTAINS

  FUNCTION ZENITH_ANGLE(point1, point2) RESULT(angle)

    TYPE(POINT), INTENT(IN) :: point1, point2
    REAL(8) :: angle

    angle = ATAN2(NORM2(point2 % arr(:2) - point1 % arr(:2)), &
                  point2 % z() - point1 % z())

  END FUNCTION ZENITH_ANGLE

  FUNCTION AZIMUTH_ANGLE(point1, point2) RESULT(angle)

    TYPE(POINT), INTENT(IN) :: point1, point2
    REAL(8) :: angle

    ! ATAN2 is in the (-PI;PI) range, MODULO takes it into (0;2PI)
    angle = MODULO(ATAN2(point2 % y() - point1 % y(), &
                         point2 % x() - point1 % x()), 2 * PI)

  END FUNCTION AZIMUTH_ANGLE

  FUNCTION ANGLE3PT(point1, point2, point3) RESULT(angle)
    TYPE(POINT), INTENT(IN) :: point1, point2, point3
    REAL(8) :: angle, arg

    TYPE(POINT) :: u, v

    u = point2 - point1
    v = point3 - point2

    arg = u * v / (NORM2(u % arr) * NORM2(v % arr))

    IF (arg >= 1) THEN
      angle = 0
    ELSE IF (arg <= -1) THEN
      angle = PI
    ELSE
      angle = ACOS(arg)
    END IF

  END FUNCTION ANGLE3PT

END MODULE ANGLE_M
