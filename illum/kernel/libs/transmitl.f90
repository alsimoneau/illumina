!-----------------------------------------------------------------------

!=======================================================================
! Routine transmitl

! Determine la transmittance des aerosols sur un parcours entre les
! cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
! Est fonction de l'epaisseur optique des aerosols
! Recoit des longueurs d'ondes en nanometre et les transforme en microns.
! La pression doit etre en KPa
! Retourne la transmittance transl

! pour utilisation avec Illumina
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

SUBROUTINE transmitl(angz, z_i, z_f, distd, hlay, transl, tranal)
  ! declaration des variables.
  REAL :: transl
  ! vertical transmittance of the complete atmosphere (aerosols)
  REAL :: tranal
  REAL :: angz, hlay
  REAL :: z_i, z_f, z1, z2
  REAL :: distd

  IF (z_i > z_f) THEN
    z2 = z_i
    z1 = z_f
  ELSE
    z1 = z_i
    z2 = z_f
  END IF
  IF (z1 /= z2) THEN
    transl = EXP((LOG(tranal) / ABS(COS(angz))) &
                 * (EXP(-1.0 * z1 / hlay) - EXP(-1.0 * z2 / hlay)))
  ELSE
    transl = EXP((LOG(tranal)) * EXP(-1.0 * z1 / hlay) * distd / hlay)
  END IF
  IF (transl == 0.0) THEN
    PRINT *, 'ERREUR transl - no transmission', z_i, z_f, angz
  END IF
  IF (transl > 1.0) THEN
    PRINT *, 'ERREUR avec transa', transl, z_i, z_f, angz
    STOP
  END IF
  RETURN
END SUBROUTINE transmitl
