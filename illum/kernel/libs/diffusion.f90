!------------------------------------------------------------------------

!=======================================================================
!  Routine diffusion

! Determine the probability of light scattering per unit of solid angle
! in the angdif direction. The scattering parameters are given
! by secdif (ratio of the effective section of scattering on the section
! of extinction), fonc_anorm (scattering phase function

! Returns the scattering probability pdif

! an Illumina routine
!-----------------------------------------------------------------------

!    Copyright (C) 2021  Martin Aube

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

!    Contact: martin.aube@cegepsherbrooke.qc.ca

SUBROUTINE diffusion(angdif, tranam, tranaa, tranal, un, secdif, secdil, &
                     fonc_a, fonc_l, haer, hlay, pdif, altit)
  REAL :: angdif, pdif, prob_a, prob_m, prob_l, secdif, secdil
  REAL :: fctmol, pi, fonc_a(181), fonc_l(181), fonc_ae, fonc_le
  REAL :: angdeg, tranam, tranaa, tranal
  REAL :: altit, un, hlay, haer
  INTEGER :: rang
  PARAMETER(pi=3.1415926)
!--------------------------------------------------------
  IF (angdif < 0.0) angdif = -angdif
  IF (angdif - pi > 0.00001) angdif = pi
  angdeg = ((angdif * 180.0) / pi)
  rang = INT(angdeg) + 1
!----------------------------------------
!  Calculate scattering probability per unit of steradian
!  The probability is for a voxel of 1x1x1m refer to equation 1 in Aubé et al.

!  Aubé, M., Simoneau, A., Muñoz-Tuñón, C., Díaz-Castro, J., &
!  Serra-Ricart, M. (2020). Restoring the night sky darkness at
!  Observatorio del Teide: First application of the model Illumina
!  version 2.0 Monthly Notices of the Royal Astronomical Society,
!  497(3), 2501-2516.0
!----------------------------------------
  IF ((tranaa <= 1.0) .AND. (tranaa > 0.0)) THEN
    ! value of the aerosol phase function
    fonc_ae = fonc_a(rang)
    ! Functions are normalized in the main code. See their division by 4pi
    prob_a = (1.0 - EXP(LOG(tranaa) * EXP(-1.0 * altit / haer) * un / haer)) * secdif * fonc_ae
  ELSE
    prob_a = 0.0
  END IF
  IF ((tranal <= 1.0) .AND. (tranal > 0.0)) THEN
    ! value of the layer phase function
    fonc_le = fonc_l(rang)
    prob_l = (1.0 - EXP(LOG(tranal) * EXP(-1.0 * altit / hlay) * un / hlay)) * secdil * fonc_le
  ELSE
    prob_l = 0.0
  END IF
  ! value of the molecule phase function
  fctmol = 0.75 * (1.0 + ((COS(angdif))**2.0)) / (4.0 * pi)
  prob_m = (1.0 - EXP(LOG(tranam) * EXP(-1.0 * altit / 8000.0) * un / 8000.0)) * fctmol

  ! This is an approximation valide if 1-transa,
  ! 1-transm et 1-transl are small
  pdif = prob_a + prob_m + prob_l
  IF (prob_a > 1.0) THEN
    PRINT *, 'prob_a>1.0'
    STOP
  END IF
  IF (prob_a < 0.0) THEN
    PRINT *, 'prob_a<0.0.'
    STOP
  END IF
  IF (prob_m > 1.0) THEN
    PRINT *, 'prob_m>1.0'
    STOP
  END IF
  IF (prob_m < 0.0) THEN
    PRINT *, 'prob_m<0.0.'
    STOP
  END IF
  IF (prob_l > 1.0) THEN
    PRINT *, 'prob_l>1.0'
    STOP
  END IF
  IF (prob_l < 0.0) THEN
    PRINT *, 'prob_l<0.0.'
    STOP
  END IF
  IF (pdif > 1.0) THEN
    pdif = 1.0
  END IF
  IF (pdif < 0.0) THEN
    pdif = 0.0
  END IF
  RETURN
END SUBROUTINE diffusion
