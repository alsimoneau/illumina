!-----------------------------------------------------------------------

!=======================================================================
! Routine transmitm (Andre Morin, Alex Neron, Etienne Rousseau 2004)
! debuggee par Martin Aube 2004
! Determine la transmittance des molecules atmospheriques sur un parcours
! entre les cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
! Recoit des longueurs d'ondes en nanometre et les transforme en microns.
! La pressi doit etre en KPa
! Retourne la transmittance transm

!  *** J'ai valide le calcul zenith tout atm avec modtran et
!      cela concorde M. Aubé mars 2010

! pour utilisation avec Illumina
!-----------------------------------------------------------------------

!    Copyright (C) 2009  Martin Aube

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

SUBROUTINE transmitm(angz, z_i, z_f, distd, transm, tranam, tabs)
  ! declaration des variables.
  REAL :: transm
  ! vertical transmittance of the complete atmosphere (molecules)
  REAL :: tranam
  REAL :: angz
  REAL :: distd
  REAL :: z_i, z_f, z1, z2
  IF (z_i > z_f) THEN
    z2 = z_i
    z1 = z_f
  ELSE
    z1 = z_i
    z2 = z_f
  END IF
  IF (z1 /= z2) THEN
    transm = EXP((LOG(tranam * tabs) / ABS(COS(angz))) * &
                 (EXP(-1.0 * z1 / 8000.0) - EXP(-1.0 * z2 / 8000.0)))
  ELSE
    transm = EXP((LOG(tranam * tabs)) * EXP(-1.0 * z1 / 8000.0) * distd / 8000.0)
  END IF
  IF (distd == 0.0) transm = 1.0
  IF ((transm < 0.0) .OR. (transm > 1.0)) THEN
    PRINT *, 'ERREUR avec transm', transm, tranam, &
      z_f, z_i, distd, angz
    STOP
  END IF
  RETURN
END SUBROUTINE transmitm
