MODULE ATMOSPHERE_M

  USE MATH_M, ONLY: PI, RAD2RANK
  IMPLICIT NONE

CONTAINS

  FUNCTION TRANSMITTANCE(zenith_angle, z_i, z_f, dist, scale_height, trans_toa)

    ! Computes the transmittance between two points of height 'z_i' and 'z_f'
    ! with an horizontal distance 'dist' between the two for a given
    ! scale height and vertical transmittance

    REAL(8), INTENT(IN) :: trans_toa
    REAL(8), INTENT(IN) :: zenith_angle, scale_height
    REAL(8), INTENT(IN) :: z_i, z_f
    REAL(8), INTENT(IN) :: dist
    REAL(8) :: transmittance

    IF (dist == 0) THEN
      transmittance = 1
      RETURN
    END IF

    IF (z_i == z_f) THEN
      transmittance = trans_toa**(EXP(-1 * z_i / scale_height) * dist / scale_height)
      RETURN
    END IF

    transmittance = trans_toa**ABS((EXP(-1 * z_i / scale_height) &
                                    - EXP(-1 * z_f / scale_height)) &
                                   / COS(zenith_angle))

    IF (transmittance > 1) THEN
      PRINT *, 'ERROR in trans', transmittance, z_i, z_f, zenith_angle
      STOP
    END IF

  END FUNCTION

  FUNCTION TRANS_VERTICAL(pression, wavelength)

    ! Computes the vertical atmospheric transmittance for a given
    ! pression [kPa] and wavelength [nm] based on Bodhaine et al. (1999)

    REAL(8), INTENT(IN) :: pression, wavelength
    REAL(8) :: trans_vertical, wav2

    wav2 = (wavelength / 1000)**2

    trans_vertical = EXP(-1 * pression / 101.325D0 * 0.0021520D0 &
                         * (1.0455996D0 - 341.29061D0 / wav2 - 0.90230850D0 * wav2) &
                         / (1 + 0.0027059889D0 / wav2 - 85.968563D0 * wav2))

  END FUNCTION

  FUNCTION DIFFUSION(scattering_angle, altitude, &
                     transmittance_molecular, scale_height_molecular, scattering_ratio, &
                     transmittance_aerosol, scale_height_aerosol, phase_function_aerosol, &
                     N_aerosol_layer, N_angle)

    ! Calculate scattering probability per unit of steradian for a voxel of 1x1x1m
    ! Refers to equation 1 in Aubé et al. (2020)

    ! Array sizes
    !F2PY INTENT(HIDE) N_aerosol_layer, N_angle
    INTEGER, INTENT(IN) :: N_aerosol_layer, N_angle

    ! Inputs
    REAL(8), INTENT(IN) :: scattering_angle, altitude
    REAL(8), INTENT(IN) :: transmittance_molecular
    REAL(8), INTENT(IN) :: scale_height_molecular
    REAL(8), INTENT(IN) :: phase_function_aerosol(N_aerosol_layer, N_angle)
    REAL(8), INTENT(IN) :: scattering_ratio(N_aerosol_layer)
    REAL(8), INTENT(IN) :: transmittance_aerosol(N_aerosol_layer)
    REAL(8), INTENT(IN) :: scale_height_aerosol(N_aerosol_layer)

    REAL(8) :: diffusion

    ! Internal variables
    INTEGER :: rang, layer
    REAL(8) :: angle
    REAL(8) :: phase_function_molecular

    angle = scattering_angle
    IF (angle < 0) angle = -angle
    IF (angle - PI > 0.0001) angle = PI
    rang = RAD2RANK(angle, N_angle)

    phase_function_molecular = 0.75 * (1 + COS(angle)**2) / (4 * PI)
    diffusion = phase_function_molecular &
                * (1 - TRANSMITTANCE(0.0D0, altitude, altitude, 1.0D0, &
                                     scale_height_molecular, &
                                     transmittance_molecular))

    DO layer = 1, N_aerosol_layer
      IF ((transmittance_aerosol(layer) <= 1) .AND. &
          (transmittance_aerosol(layer) > 0)) THEN
        diffusion = diffusion + scattering_ratio(layer) &
                    * phase_function_aerosol(layer, rang) &
                    * (1 - TRANSMITTANCE(0.0D0, altitude, altitude, 1.0D0, &
                                         scale_height_aerosol(layer), &
                                         transmittance_aerosol(layer)))
      END IF
    END DO

    IF (diffusion > 1) THEN
      diffusion = 1
    END IF
    IF (diffusion < 0) THEN
      diffusion = 0
    END IF

  END FUNCTION DIFFUSION

END MODULE ATMOSPHERE_M
