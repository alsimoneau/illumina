SUBROUTINE average_index(nx, ny, index_max, indices, data, average)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, index_max
  INTEGER(KIND=4), INTENT(IN) :: indices(nx, ny)
  REAL(KIND=8), INTENT(IN) :: data(nx, ny)
  REAL(KIND=8), INTENT(OUT) :: average(nx, ny)

  INTEGER :: i, j, z
  INTEGER :: N(index_max)
  REAL :: acc(index_max)

  DO i = 1, index_max
    acc(i) = 0
    N(i) = 0
  END DO

  DO j = 1, ny
    DO i = 1, nx
      z = indices(i, j)
      N(z) = N(z) + 1
      acc(z) = acc(z) + data(i, j)
    END DO
  END DO

  DO i = 1, index_max
    IF (N(i) .EQ. 0) THEN
      acc(i) = 0
    ELSE
      acc(i) = acc(i) / N(i)
    END IF
  END DO

  DO j = 1, ny
    DO i = 1, nx
      average(i, j) = acc(indices(i, j))
    END DO
  END DO

END SUBROUTINE
