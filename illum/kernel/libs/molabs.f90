!  adding molecular absorption

SUBROUTINE molabs(lambda, bandw, tabs)
  REAL :: lambda, bandw, tabs, wmax, wmin, wla, ta, bwa, bwabs
  wmin = lambda - bandw / 2.0
  wmax = lambda + bandw / 2.0
  tabs = 0.0
  wla = 0.0
  bwabs = 0.0
  OPEN (unit=1, file='MolecularAbs.txt', status='old')
  READ (1, *)
  DO WHILE (wla < wmax)
    READ (1, *) wla, ta, bwa
    IF (wla > wmin) THEN
      tabs = tabs + ta * bwa
      bwabs = bwabs + bwa
    END IF
  END DO
  IF (bwabs > 0.0) tabs = tabs / bwabs
  CLOSE (unit=1)
  RETURN
END SUBROUTINE molabs
