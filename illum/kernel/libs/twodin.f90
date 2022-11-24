SUBROUTINE twodin(nbx, nby, filename, bindata)
  INTEGER :: width                                                      ! matrix dimension in length/width and height
  PARAMETER(width=512)
! read double precision array in binary
  INTEGER :: nbx, nby, i, j
  REAL :: bindata(width, width)
  CHARACTER*72 :: filename
  OPEN (unit=1, form='unformatted', file=filename, action='read')
  READ (1) nbx, nby
  IF ((nbx > width) .OR. (nby > width)) THEN
    PRINT *, 'You try to use a domain larger than the maximum'
    PRINT *, 'allowed. Please restrict it to no more that 512 x 512'
    PRINT *, 'Your domain size is: ', nbx, 'x', nby
    PRINT *, 'Computation aborted'
    STOP
  END IF
  DO j = nby, 1, -1
    DO i = 1, nbx
      READ (1) bindata(i, j)
    END DO
  END DO
  CLOSE (unit=1)
  RETURN
END SUBROUTINE twodin
