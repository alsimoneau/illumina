SUBROUTINE average_index(arr, indices, average, n, m)

!$ USE OMP_LIB

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

  !$OMP PARALLEL DO REDUCTION(+:weight,acc)
  DO i = 1, n
    ASSOCIATE (z => indices(i))
      weight(z) = weight(z) + 1
      DO j = 1, m
        acc(j, z) = acc(j, z) + arr(j, i)
      END DO
    END ASSOCIATE
  END DO
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO
  DO i = 0, index_max
    DO j = 1, m
      acc(j, i) = acc(j, i) / weight(i)
    END DO
  END DO
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO
  DO i = 1, n
    DO j = 1, m
      average(j, i) = acc(j, indices(i))
    END DO
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE
