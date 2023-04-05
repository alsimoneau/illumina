MODULE GEOM_M

  USE MATH_M, ONLY: PI
  IMPLICIT NONE

CONTAINS

  FUNCTION LINE(point1, point2, coords) RESULT(pts)

    INTEGER(4), DIMENSION(2), INTENT(IN) :: point1, point2
    REAL(8), INTENT(IN) :: coords(:, :, :) ! x/y, r, theta

    INTEGER(4), ALLOCATABLE :: pts(:, :) ! i/j, N

    REAL(8) :: angle_offset, angle_new, angle_prev, coord1(2), coord2(2)
    INTEGER(4) :: i, j, r, v, w, dir = 1, idir
    INTEGER(4), ALLOCATABLE :: found_pts(:)
    LOGICAL :: found

    IF (point1(2) == point2(2)) THEN
      ALLOCATE (pts(2, ABS(point1(1) - point2(1)) + 1))
      idir = 1
      IF (point2(1) < point1(1)) idir = -1
      pts(1, :) = [(i, i=point1(1), point2(1), idir)]
      pts(2, :) = point1(2)
      RETURN
    END IF

    ALLOCATE (found_pts(0))

    IF (point1(2) > point2(2)) THEN
      dir = -1 * dir
    END IF

    IF (ABS(point1(2) - point2(2)) > SIZE(coords, 2) / 2) THEN
      dir = -1 * dir
    END IF

    coord1 = coords(:, point1(1), point1(2))
    coord2 = coords(:, point2(1), point2(2))
    angle_offset = ATAN2(coord2(1) - coord1(1), &
                         coord2(2) - coord1(2))

    r = point1(1)
    w = point1(2)
    DO WHILE (w /= point2(2))
      v = w
      w = MODULO(v + dir - 1, SIZE(coords, 3)) + 1
      found = .FALSE.
      angle_prev = PI
      DO i = 1, SIZE(coords, 2)
        angle_new = MODULO(ATAN2(coords(1, i, w) - coord1(1), &
                                 coords(2, i, w) - coord1(2)) - angle_offset, 2 * PI)
        IF (ABS(angle_new - angle_prev) > PI) THEN
          found = .TRUE.
          EXIT
        END IF
        angle_prev = angle_new
      END DO

      IF (.NOT. found) i = 1
      PRINT *, w, i

      idir = 1
      IF (i < r) idir = -1
      found_pts = [found_pts, &
                   [([j, v], j=r, (r + i) / 2, idir)], &
                   [([j, w], j=(r + i) / 2 + idir, i + idir, idir)]]

      r = i
    END DO

    pts = RESHAPE(found_pts, [2, SIZE(found_pts) / 2])

  END FUNCTION LINE

END MODULE GEOM_M
