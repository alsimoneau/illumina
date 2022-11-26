! Based on Shapiro 1982 Table 11

FUNCTION cloud_reflectance(zenith_angle, cloud_type) RESULT(refl)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: cloud_type
  REAL, INTENT(IN) :: zenith_angle
  REAL :: refl

  REAL :: c(4), z

  IF (cloud_type == 1) THEN       ! thin cirrus & cirrostratus
    c = [0.25674, -0.18077, -0.21961, 0.25272]
  ELSE IF (cloud_type == 2) THEN  ! thick cirrus & cirrostratus
    c = [0.60540, -0.55142, -0.23389, 0.43648]
  ELSE IF (cloud_type == 3) THEN  ! altostratus & altocumulus
    c = [0.66152, -0.14863, -0.08193, 0.13442]
  ELSE IF (cloud_type == 4) THEN  ! stratocumulus & stratus
    c = [0.71214, -0.15033, 0.00696, 0.03904]
  ELSE IF (cloud_type == 5) THEN  ! cumulus & cumulonimbus
    c = [0.67072, -0.13805, -0.10895, 0.09460]
  END IF

  z = COS(zenith_angle)
  refl = c(1) + c(2) * z + c(3) * z**2 + c(4) * z**3

END FUNCTION

FUNCTION cloud_transmitance(zenith_angle, cloud_type) RESULT(trans)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: cloud_type
  REAL, INTENT(IN) :: zenith_angle
  REAL :: trans

  REAL :: c(4), z

  IF (cloud_type == 1) THEN       ! thin cirrus & cirrostratus
    c = [0.63547, 0.35229, 0.08709, -0.22902]
  ELSE IF (cloud_type == 2) THEN  ! thick cirrus & cirrostratus
    c = [0.26458, 0.66829, 0.24228, -0.49357]
  ELSE IF (cloud_type == 3) THEN  ! altostratus & altocumulus
    c = [0.19085, 0.32817, -0.08613, -0.08197]
  ELSE IF (cloud_type == 4) THEN  ! stratocumulus & stratus
    c = [0.13610, 0.29964, -0.14041, 0.00952]
  ELSE IF (cloud_type == 5) THEN  ! cumulus & cumulonimbus
    c = [0.17960, 0.34855, -0.14875, 0.01962]
  END IF

  z = COS(zenith_angle)
  trans = c(1) + c(2) * z + c(3) * z**2 + c(4) * z**3

END FUNCTION
