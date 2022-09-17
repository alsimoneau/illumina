SUBROUTINE average_index(arr, indices, index_max, average, n, m)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n, m, index_max
  INTEGER(4), INTENT(IN) :: indices(n)
  REAL(8), INTENT(IN) :: arr(m, n)
  REAL(8), INTENT(OUT) :: average(m, n)

  INTEGER :: i
  INTEGER :: weight(0:index_max)
  REAL(8) :: acc(m, 0:index_max)

  acc = 0
  weight = 0

  DO i = 1, n
    ASSOCIATE (z => indices(i))
      weight(z) = weight(z) + 1
      acc(:, z) = acc(:, z) + arr(:, i)
    END ASSOCIATE
  END DO

  acc(:, 0) = 0
  DO i = 1, index_max
    acc(:, i) = acc(:, i) / weight(i)
  END DO

  DO i = 1, n
    average(:, i) = acc(:, indices(i))
  END DO

END SUBROUTINE
