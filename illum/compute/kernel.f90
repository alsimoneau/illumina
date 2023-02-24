SUBROUTINE KERNEL( &
  grid_x, grid_y, ground_elevation, pixel_size, aerosols_layers, &
  wavelengths, bandwidth, reflectance, reflection_radius, stop_limit, &
  viewing_elevation, viewing_zenith, viewing_azimuth, &
  cloud_type, cloud_base_height, cloud_fraction, &
  compute_first_scattering, compute_second_scattering, &
  N_theta, N_radii, N_aerosols, N_wavelength, &
  radiance)

  IMPLICIT NONE

  ! Array lenghts
!F2PY INTENT(HIDE) N_theta, N_radii, N_aerosols, N_wavelength
  INTEGER(4), INTENT(IN) :: N_theta, N_radii, N_aerosols, N_wavelength

  ! Input parameters
  REAL(8), INTENT(IN), DIMENSION(N_theta, N_radii) :: grid_x, grid_y, ground_elevation
  REAL(8), INTENT(IN), DIMENSION(N_radii) :: pixel_size
  REAL(8), INTENT(IN), DIMENSION(3, N_aerosols) :: aerosols_layers ! optical depth, angstrom exponent, scale height
  REAL(8), INTENT(IN), DIMENSION(N_wavelength) :: wavelengths
  REAL(8), INTENT(IN) :: bandwidth, reflection_radius, stop_limit
  REAL(8), INTENT(IN), DIMENSION(N_theta, N_radii, N_wavelength) :: reflectance
  REAL(8), INTENT(IN) :: viewing_elevation, viewing_zenith, viewing_azimuth
  INTEGER(4), INTENT(IN) :: cloud_type
  REAL(8), INTENT(IN) :: cloud_base_height, cloud_fraction
  LOGICAL, INTENT(IN) :: compute_first_scattering, compute_second_scattering

  ! Return values
  REAL(8), INTENT(OUT), DIMENSION(N_wavelength) :: radiance

END SUBROUTINE KERNEL
