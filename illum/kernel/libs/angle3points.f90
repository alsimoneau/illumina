!-----------------------------------------------------------------------

!=======================================================================
!  Routine angle3points (Andre Morin 2004)
!  debuggee et modifiee par Martin Aube 2004

!  Determine l'angle entre 3 points (x1,y1,z1), (x2,y2,z2) et (x3,y3,z3)
!  dont le sommet est au point 2
!  Retourne l'angle angle en radians

!  pour utilisation avec Illumina
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

SUBROUTINE angle3points(x1, y1, z1, x2, y2, z2, x3, y3, z3, an3pts)
  REAL :: x1, y1, x2, y2, x3, y3
  REAL :: z1, z2, z3, an3pts, argume
  REAL :: xu, yu, zu, xv, yv, zv, pi
  PARAMETER(pi=3.141592654)
  xu = x2 - x1
! Voici les composantes du vecteur u.
  yu = y2 - y1
  zu = z2 - z1
  xv = x3 - x2
! Voici les composantes du vecteur v.
  yv = y3 - y2
  zv = z3 - z2
  IF ((xv == 0.0) .AND. (yv == 0.0) .AND. (zv == 0.0)) THEN
    PRINT *, 'ERREUR vecteur sortie nul'
    PRINT *, x1, y1, z1, x2, y2, z2, x3, y3, z3
    STOP
  END IF
  IF ((xu == 0.0) .AND. (yu == 0.0) .AND. (zu == 0.0)) THEN
    PRINT *, 'ERREUR vecteur d entree nul'
    PRINT *, x1, y1, z1, x2, y2, z2, x3, y3, z3
    STOP
  END IF
  argume = (xu * xv + yu * yv + zu * zv) / (SQRT(xu**2.0 + yu**2.0 + zu**2.0) * &
                                            SQRT(xv**2.0 + yv**2.0 + zv**2.0))
  IF (argume >= 1.0) THEN
    an3pts = 0.0
  ELSE IF (argume <= -1.0) THEN
    an3pts = pi
  ELSE
    an3pts = ACOS(argume)
  END IF
  IF (an3pts < 0.0) THEN
    PRINT *, 'ERREUR an3pts < 0'
    STOP
  END IF
  RETURN
END SUBROUTINE angle3points
