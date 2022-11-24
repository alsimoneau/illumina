!  adding earth curvature to the topography

SUBROUTINE curvature(distc, hcur)
  REAL :: rearth, distc, hcur
  Rearth = 6371000.0
  hcur = Rearth - SQRT(Rearth**2.0 - distc**2.0)
  hcur = -1.0 * hcur
  RETURN
END SUBROUTINE curvature
