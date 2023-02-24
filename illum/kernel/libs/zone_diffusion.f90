!=======================================================================
!  Routine zone diffusion 2010

!  Determine les cellules se trouvant dans la zone de diffusion (effet)
!  des cellules (x1,y1,z1) et (x2,y2,z2)
!  Retourne la matrice des cellules diffusantes (diffusion) ainsi que le nombre
!  de cellules diffusantes (ncellule)

!------------------------------------------------------------------------

!    Copyright (C) 2010  Martin Aube

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

SUBROUTINE zone_diffusion(effet, zondif, ncell, stepdi, siz)

  IMPLICIT NONE

  REAL, INTENT(IN) :: effet, siz
  INTEGER, INTENT(IN) :: stepdi
  REAL, INTENT(OUT) :: zondif(3000000, 3)
  INTEGER, INTENT(OUT) :: ncell

  INTEGER :: i, j, k
  INTEGER :: neffet, keep
  REAL :: x0, y0, z0
  REAL :: d2, siz2, dmin2

  neffet = NINT(effet / siz)

  ncell = 0
  siz2 = siz**2
  dmin2 = effet**2
  keep = 0

  DO i = -neffet, neffet
    x0 = REAL(i) * siz
    DO j = -neffet, neffet
      y0 = REAL(j) * siz
      DO k = -neffet, neffet
        z0 = REAL(k) * siz
        d2 = x0**2.0 + y0**2.0 + z0**2.0
        IF (d2 <= dmin2) THEN
          keep = keep + 1
          IF (keep == stepdi) THEN
            keep = 0
            ncell = ncell + 1
            zondif(ncell, :) = [x0, y0, z0]
          END IF
        END IF
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE zone_diffusion
