!=======================================================================
!  Routine horizon (Martin Aube 2017)

!  Determine si la lumiere est bloquee par l'horizon

!  pour utilisation avec Illumina
!-----------------------------------------------------------------------

!    Copyright (C) 2019  Martin Aube

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

SUBROUTINE horizon(x, y, z, dx, dy, altsol, anga, zhoriz, d)
  ! matrix dimension in length/width and height
  INTEGER :: width
  PARAMETER(width=512)
  INTEGER :: x, y, nx, ny
  REAL :: dx, dy, altsol(width, width), anga, zout, pi, angaz1, ix, iy
  ! earth curvature terrain
  REAL :: hcur, distc
  REAL :: posx, posy, scalef, zhoriz, z, d, dout
  pi = 3.141592654

  angaz1 = anga
  ! viewing vector components
  ix = (COS(angaz1))
  iy = (SIN(angaz1))
  scalef = dx / 3.0
  posx = REAL(x) * dx
  posy = REAL(y) * dy
  zhoriz = pi
  d = (REAL(width)) * SQRT(1.0 + TAN(angaz1)**2.0)

  DO WHILE (((posx <= REAL(width) * dx) .AND. (posx >= 1.0 * dx)) .AND. &
            ((posy <= REAL(width) * dy) .AND. (posy >= 1.0 * dy)))
    posx = posx + ix * scalef
    posy = posy + iy * scalef
    nx = NINT(posx / dx)
    ny = NINT(posy / dy)

    ! earth curvature (first order correction)
    distc = SQRT((dx * REAL(nx - x))**2.0 + (dy * REAL(ny - y))**2.0)
    CALL curvature(distc, hcur)
    ! to forbid division by zero
    IF ((nx == x) .AND. (ny == y)) THEN
      ! reverse curvature to limit horizontal distance. curv is negative. this is a hack
      IF (z > altsol(nx, ny) - hcur) THEN
        zout = pi
        d = 0.0
      ELSE
        zout = 0.0
        d = 0.0
      END IF
    ELSE
      dout = distc
      zout = pi / 2.0 - ATAN((altsol(nx, ny) - hcur - z) / dout)
      IF (altsol(nx, ny) - hcur == z) THEN
        ! bug for zhoriz = pi, anyway in the real world pi is almost impossible
        zout = pi / 2.0 - 0.0001 * pi / 180.0
      END IF
      IF (zout < zhoriz) THEN
        zhoriz = zout
        d = dout
      END IF
    END IF
  END DO
END SUBROUTINE horizon
