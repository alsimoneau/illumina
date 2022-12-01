SUBROUTINE average_index(arr, indices, average, n, m)

  IMPLICIT NONE

!f2py intent(hide) n, m
  INTEGER(4), INTENT(IN) :: n, m
  INTEGER(4), INTENT(IN) :: indices(n)
  REAL(8), INTENT(IN) :: arr(m, n)
  REAL(8), INTENT(OUT) :: average(m, n)

  INTEGER(4) :: i, j, index_max
  INTEGER(4), ALLOCATABLE :: weight(:)
  REAL(8), ALLOCATABLE :: acc(:, :)

  index_max = MAXVAL(indices)
  ALLOCATE (weight(0:index_max))
  ALLOCATE (acc(m, 0:index_max))

  acc = 0
  weight = 0

  DO i = 1, n
    ASSOCIATE (z => indices(i))
      weight(z) = weight(z) + 1
      acc(:, z) = acc(:, z) + arr(:, i)
    END ASSOCIATE
  END DO

  DO i = 0, index_max
    acc(:, i) = acc(:, i) / weight(i)
  END DO

  DO i = 1, n
    average(:, i) = acc(:, indices(i))
  END DO

END SUBROUTINE
