FUNCTION trans(zenith_angle, z_i, z_f, distd, scale_height, trana)

  IMPLICIT NONE

  REAL, INTENT(in) :: trana
  REAL, INTENT(in) :: angz, haer
  REAL, INTENT(in) :: z_i, z_f, z1, z2
  REAL, INTENT(in) :: distd

  REAL :: trans

  IF (distd == 0.0) THEN
    trans = 1.0
    RETURN
  END IF

  IF (z_i == z_f) THEN
    trans = EXP((LOG(trana)) * EXP(-1.0 * z_i / scale_height) * distd / scale_height)
    RETURN
  END IF

  IF (z_i > z_f) THEN
    z2 = z_i
    z1 = z_f
  ELSE
    z1 = z_i
    z2 = z_f
  END IF

  trans = EXP((LOG(trana) / ABS(COS(zenith_angle))) &
              * (EXP(-1.0 * z1 / scale_height) - EXP(-1.0 * z2 / scale_height)))

  IF (trans == 0.0) THEN
    PRINT *, 'ERREUR transa - no transmission', z_i, z_f, zenith_angle, trana, distd, scale_height
  END IF

  IF (trans > 1.0) THEN
    PRINT *, 'ERREUR avec transa', transa, z_i, z_f, zenith_angle
    STOP
  END IF

END FUNCTION

FUNCTION trans_toa(pression, wavelength) RESULT(trans)

  IMPLICIT NONE

  REAL :: m, lambm, taua, pressi, tranam, tranaa, tranal, layaod, tabs, bandw

  !  From Bodhaine et al. (1999)
  trans = EXP(-1.0 * pression / 101.325 * 0.0021520 &
              * (1.0455996 - 341.29061 * (lambm / 1000.0)**-2.0 &
                 - 0.90230850 * (lambm / 1000.0)**2.0) &
              / (1.0 + 0.0027059889 * (lambm / 1000.0)**-2.0 &
                 - 85.968563 * (lambm / 1000.0)**2.0))

END FUNCTION
