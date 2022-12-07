SUBROUTINE KERNEL( &
  view_elevation, view_azimuth, obs_elevation, air_pressure, &
  cloud_type, cloud_base, cloud_fraction, scat_single, scat_double, &
  fov, aod, alpha, scale_height, &
  Nview, Nobsh, Natm, Ncloud, Nsingle, Ndouble, Nfov, Naerosol &
  )

  ! Looping the kernel over various parameters:
  !   Viewing angle (az, el)
  !   Observer elevation
  !   Air pressure
  !   cloud (model, base height, fraction)
  !   single scattering
  !   double scattering
  !   direct fov
  !   Aerosols (aod, alpha, scale height)

!$ USE OMP_LIB

  IMPLICIT NONE

!F2PY INTENT(HIDE) Nview, Nobsh, Natm, Ncloud, Nsingle, Ndouble, Nfov, Naerosol
  INTEGER, INTENT(IN) :: Nview, Nobsh, Natm, Ncloud, Nsingle, Ndouble, Nfov
  INTEGER :: i_view, i_obsh, i_atm, i_cloud, i_single, i_double, i_fov

  REAL(8), INTENT(IN) :: view_elevation(Nview), view_azimuth(Nview)
  REAL(8), INTENT(IN) :: obs_elevation(Nobsh)
  REAL(8), INTENT(IN) :: air_pressure(Natm)
  INTEGER, INTENT(IN) :: cloud_type(Ncloud)
  REAL(8), INTENT(IN) :: cloud_base(Ncloud), cloud_fraction(Ncloud)
  LOGICAL, INTENT(IN) :: scat_single(Nsingle), scat_double(Ndouble)
  REAL(8), INTENT(IN) :: fov(Nfov)
  INTEGER, INTENT(IN) :: Naerosol
  REAL(8), INTENT(IN) :: aod(Naerosol), alpha(Naerosol), scale_height(Naerosol)

!$ CALL OMP_SET_NESTED(.TRUE.)

  !$OMP PARALLEL DO COLLAPSE(7)
  DO i_view = 1, Nview
    DO i_obsh = 1, Nobsh
      DO i_atm = 1, Natm
        DO i_cloud = 1, Ncloud
          DO i_single = 1, Nsingle
            DO i_double = 1, Ndouble
              DO i_fov = 1, Nfov
                CALL KERNEL_SINGLE( &
                  view_elevation(i_view), &
                  view_azimuth(i_view), &
                  obs_elevation(i_obsh), &
                  air_pressure(i_atm), &
                  cloud_type(i_cloud), &
                  cloud_base(i_cloud), &
                  cloud_fraction(i_cloud), &
                  scat_single(i_single), &
                  scat_double(i_double), &
                  fov(i_fov), &
                  aod, alpha, scale_height, Naerosol &
                  )
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE KERNEL_SINGLE( &
  view_elevation, view_azimuth, obs_elevation, air_pressure, &
  cloud_type, cloud_base, cloud_fraction, scat_single, scat_double, &
  fov, aod, alpha, scale_height, Naerosol &
  )

!$ USE OMP_LIB

  IMPLICIT NONE

  !F2PY INTENT(HIDE) Naerosol

  REAL(8), INTENT(IN) :: view_elevation, view_azimuth
  REAL(8), INTENT(IN) :: obs_elevation
  REAL(8), INTENT(IN) :: air_pressure
  INTEGER, INTENT(IN) :: cloud_type
  REAL(8), INTENT(IN) :: cloud_base, cloud_fraction
  LOGICAL, INTENT(IN) :: scat_single, scat_double
  REAL(8), INTENT(IN) :: fov
  INTEGER, INTENT(IN) :: Naerosol
  REAL(8), INTENT(IN) :: aod(Naerosol), alpha(Naerosol), scale_height(Naerosol)

  INTEGER :: layer

!$ CALL OMP_SET_NESTED(.TRUE.)

  PRINT *, "Elevation:", view_elevation
  PRINT *, "Azimuth:", view_azimuth
  PRINT *, "Observer elevation:", obs_elevation
  PRINT *, "Air pressure:", air_pressure
  PRINT *, "Cloud type:", cloud_type
  PRINT *, "Cloud base height:", cloud_base
  PRINT *, "Cloud fraction:", cloud_fraction
  PRINT *, "Single scattering:", scat_single
  PRINT *, "Double scattering:", scat_double
  PRINT *, "Field of view:", fov
  DO layer = 1, Naerosol
    PRINT *, "----------------"
    PRINT *, "Aerosol layer:", layer
    PRINT *, "Optical depth:", aod(layer)
    PRINT *, "Angstrom exponent:", alpha(layer)
    PRINT *, "Scale height:", scale_height(layer)
  END DO
  PRINT *, ""
  PRINT *, ""

END SUBROUTINE
