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

SUBROUTINE zone_diffusion( &
  effet, zondif, ncell, stepdi, siz)
  IMPLICIT NONE
  INTEGER :: i, j, k
  INTEGER :: ncell, neffet, imin, imax, jmin, jmax, kmin, kmax
  INTEGER :: keep, stepdi
  REAL :: x0, y0, z0
  REAL :: effet, dmin, d
  REAL :: zondif(3000000, 3), siz
  REAL :: pi
  pi = 3.141592654

  keep = 0
  neffet = NINT(effet / siz)
  dmin = effet
  stepdi = 1

  ! limits of the calculations loops
  imin = -neffet
  imax = +neffet
  jmin = -neffet
  jmax = +neffet
  kmin = -neffet
  kmax = neffet
  ncell = 0

  DO i = imin, imax
    x0 = REAL(i) * siz
    DO j = jmin, jmax
      y0 = REAL(j) * siz
      DO k = kmin, kmax
        z0 = REAL(k) * siz
        d = SQRT(x0**2.0 + y0**2.0 + z0**2.0)
        IF (d <= dmin) THEN
          keep = keep + 1
          IF (keep == stepdi) THEN
            keep = 0
            ncell = ncell + 1
            zondif(ncell, 1) = x0
            zondif(ncell, 2) = y0
            zondif(ncell, 3) = z0
          END IF
        END IF
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE zone_diffusion
