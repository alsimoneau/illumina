FUNCTION trans(zenith_angle, z_i, z_f, distd, scale_height, trana)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: trana
  REAL(8), INTENT(IN) :: zenith_angle, scale_height
  REAL(8), INTENT(IN) :: z_i, z_f
  REAL(8), INTENT(IN) :: distd

  REAL(8) :: trans, z1, z2

  IF (distd == 0) THEN
    trans = 1
    RETURN
  END IF

  IF (z_i == z_f) THEN
    trans = EXP((LOG(trana)) * EXP(-1 * z_i / scale_height) * distd / scale_height)
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
              * (EXP(-1 * z1 / scale_height) - EXP(-1 * z2 / scale_height)))

  IF (trans == 0) THEN
    PRINT *, 'ERREUR trans - no transmission', z_i, z_f, zenith_angle, trana, distd, scale_height
  END IF

  IF (trans > 1) THEN
    PRINT *, 'ERREUR avec trans', trans, z_i, z_f, zenith_angle
    STOP
  END IF

END FUNCTION

FUNCTION trans_toa(pression, wavelength) RESULT(trans)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: pression, wavelength
  REAL(8) :: trans

  !  From Bodhaine et al. (1999)
  trans = EXP(-1 * pression / 101.325D0 * 0.0021520D0 &
              * (1.0455996D0 - 341.29061D0 * (wavelength / 1000.0D0)**-2 &
                 - 0.90230850D0 * (wavelength / 1000.0D0)**2) &
              / (1 + 0.0027059889D0 * (wavelength / 1000.0D0)**-2 &
                 - 85.968563D0 * (wavelength / 1000.0D0)**2))

END FUNCTION
