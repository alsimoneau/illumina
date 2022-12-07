SUBROUTINE twodout(nbx, nby, filename, bindata)
! write double precision array in binary
  INTEGER :: nbx, nby, i, j
  REAL :: bindata(512, 512)
  CHARACTER*72 :: filename
  OPEN (unit=1, form='unformatted', file=filename, action='write')
  WRITE (1) nbx, nby
  DO j = nby, 1, -1
    DO i = 1, nbx
      WRITE (1) bindata(i, j)
    END DO
  END DO
  CLOSE (unit=1)
  RETURN
END SUBROUTINE twodout
