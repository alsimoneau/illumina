!=======================================================================
!  Routine angleazimutal (Martin Aube 2010)

!  Determine l'angle azimutal entre les points (x1,y1,z1) et (x2,y2,z2)
!  Retourne l'angle anglezen en radians

!  pour utilisation avec Illumina
!-----------------------------------------------------------------------

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

SUBROUTINE angleazimutal(x1, y1, x2, y2, angazi)
  REAL :: x1, y1, x2, y2
  REAL :: pi, angazi
  PARAMETER(pi=3.141592654)
  IF (x2 - x1 /= 0.0) angazi = ABS(ATAN((y2 - y1) / (x2 - x1)))
  IF (((x2 - x1) == 0.0) .AND. ((y2 - y1) == 0.0)) THEN
    angazi = 0.0
  ELSE
    IF (x2 - x1 > 0.0) THEN
      IF (y2 - y1 < 0.0) THEN
        angazi = 2.0 * pi - angazi
      END IF
    ELSE IF (x2 - x1 < 0.0) THEN
      IF (y2 - y1 < 0.0) THEN
        angazi = angazi + pi
      ELSE
        angazi = pi - angazi
      END IF
    ELSE     ! i.e. x2-x1 = 0.0
      IF (y2 > y1) angazi = pi / 2.0
      IF (y2 < y1) angazi = 3.0 * pi / 2.0
    END IF
    IF ((angazi < 0.0) .OR. (angazi > 2.0 * pi)) THEN
      PRINT *, 'ERREUR angazi = ', angazi, x1, y2, x2, y2
      STOP
    END IF
!        if ((x1.0eq.x2).and.(y1.0eq.y2)) then
!          print*,'ERREUR cant compute angle between identical points!'
!          stop
!        endif
  END IF
  RETURN
END SUBROUTINE angleazimutal
