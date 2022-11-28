! Based on Shapiro 1982 Table 11

FUNCTION cloud_reflectance(zenith_angle, cloud_type) RESULT(refl)

  IMPLICIT NONE

  INTEGER(4), INTENT(IN) :: cloud_type
  REAL(8), INTENT(IN) :: zenith_angle
  REAL(8) :: refl

  REAL(8) :: c(4), z

  SELECT CASE (cloud_type)
  CASE (1)  ! thin cirrus & cirrostratus
    c = [0.25674, -0.18077, -0.21961, 0.25272]
  CASE (2)  ! thick cirrus & cirrostratus
    c = [0.60540, -0.55142, -0.23389, 0.43648]
  CASE (3)  ! altostratus & altocumulus
    c = [0.66152, -0.14863, -0.08193, 0.13442]
  CASE (4)  ! stratocumulus & stratus
    c = [0.71214, -0.15033, 0.00696, 0.03904]
  CASE (5)  ! cumulus & cumulonimbus
    c = [0.67072, -0.13805, -0.10895, 0.09460]
  CASE default
    c = [1.0, 0.0, 0.0, 0.0]
  END SELECT

  z = COS(zenith_angle)
  refl = c(1) + c(2) * z + c(3) * z**2 + c(4) * z**3

END FUNCTION

FUNCTION cloud_transmitance(zenith_angle, cloud_type) RESULT(trans)

  IMPLICIT NONE

  INTEGER(4), INTENT(IN) :: cloud_type
  REAL(8), INTENT(IN) :: zenith_angle
  REAL(8) :: trans

  REAL(8) :: c(4), z

  SELECT CASE (cloud_type)
  CASE (1)        ! thin cirrus & cirrostratus
    c = [0.63547, 0.35229, 0.08709, -0.22902]
  CASE (2)   ! thick cirrus & cirrostratus
    c = [0.26458, 0.66829, 0.24228, -0.49357]
  CASE (3)   ! altostratus & altocumulus
    c = [0.19085, 0.32817, -0.08613, -0.08197]
  CASE (4)   ! stratocumulus & stratus
    c = [0.13610, 0.29964, -0.14041, 0.00952]
  CASE (5)   ! cumulus & cumulonimbus
    c = [0.17960, 0.34855, -0.14875, 0.01962]
  CASE default
    c = [1.0, 0.0, 0.0, 0.0]
  END SELECT

  z = COS(zenith_angle)
  trans = c(1) + c(2) * z + c(3) * z**2 + c(4) * z**3

END FUNCTION
