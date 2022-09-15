SUBROUTINE average_index(nx, ny, index_max, indices, arr, average)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, index_max
  INTEGER(4), INTENT(IN) :: indices(nx, ny)
  REAL(8), INTENT(IN) :: arr(nx, ny)
  REAL(8), INTENT(OUT) :: average(nx, ny)

  INTEGER :: i, j
  INTEGER :: N(0:index_max)
  REAL(8) :: acc(0:index_max)

  acc = 0
  N = 0

  DO j = 1, ny
    DO i = 1, nx
      ASSOCIATE (z => indices(i, j))
        N(z) = N(z) + 1
        acc(z) = acc(z) + arr(i, j)
      END ASSOCIATE
    END DO
  END DO

  acc = MERGE(0._8, acc / N, N == 0)

  DO j = 1, ny
    DO i = 1, nx
      average(i, j) = acc(indices(i, j))
    END DO
  END DO

END SUBROUTINE
