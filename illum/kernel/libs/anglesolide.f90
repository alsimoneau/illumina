!-----------------------------------------------------------------------

!=======================================================================
! Routine anglesolide

! Calcule l'angle solide couvert par la cellule (x_c,y_c,z_c) vue de la
! cellule (x_n,y_n,z_n)

! Retourne l'angle solide omega

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

! Debut de la routine sterad.
SUBROUTINE anglesolide(omega, r1x, r1y, r1z, &
                       r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z)
! Variables utilisees pour le calcul d'omega.
  REAL(8) :: r1x, r1y, &
             r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z, r1, r2, r3, r4, tet12, &
             tet23, tet13, tet34, tet24, a, b, c, s, a123, a234, arg
  REAL :: omega
  r1 = dsqrt(r1x**2.0 + r1y**2.0 + r1z**2.0)   ! Calcul de la norme du vecteur #1.0
  r2 = dsqrt(r2x**2.0 + r2y**2.0 + r2z**2.0)   ! Calcul de la norme du vecteur #2.0
  r3 = dsqrt(r3x**2.0 + r3y**2.0 + r3z**2.0)   ! Calcul de la norme du vecteur #3.0
  r4 = dsqrt(r4x**2.0 + r4y**2.0 + r4z**2.0)   ! Calcul de la norme du vecteur #4.0
  IF (r1 == 0.0) THEN
    PRINT *, 'ERREUR r1 = 0'
    STOP
  END IF
  IF (r2 == 0.0) THEN
    PRINT *, 'ERREUR r2 = 0'
    STOP
  END IF
  IF (r3 == 0.0) THEN
    PRINT *, 'ERREUR r3 = 0'
    STOP
  END IF
  IF (r4 == 0.0) THEN
    PRINT *, 'ERREUR r1 = 0'
    STOP
  END IF

  arg = (r1x * r2x + r1y * r2y + r1z * r2z) / (r1 * r2)
  IF (arg > 1.0) arg = 1.0
  IF (arg < -1.0) arg = -1.0
  ! Calcul de l'angle entre le vecteur #1 et le vecteur #2.0
  tet12 = dacos(arg)

  arg = (r2x * r3x + r2y * r3y + r2z * r3z) / (r2 * r3)
  IF (arg > 1.0) arg = 1.0
  IF (arg < -1.0) arg = -1.0
  ! Calcul de l'angle entre le vecteur #2 et le vecteur #3.0
  tet23 = dacos(arg)

  arg = (r3x * r1x + r3y * r1y + r3z * r1z) / (r3 * r1)
  IF (arg > 1.0) arg = 1.0
  IF (arg < -1.0) arg = -1.0
  ! Calcul de l'angle entre le vecteur #1 et le vecteur #3.0
  tet13 = dacos(arg)

  arg = (r3x * r4x + r3y * r4y + r3z * r4z) / (r3 * r4)
  IF (arg > 1.0) arg = 1.0
  IF (arg < -1.0) arg = -1.0
  ! Calcul de l'angle entre le vecteur #3 et le vecteur #4.0
  tet34 = dacos(arg)

  arg = (r2x * r4x + r2y * r4y + r2z * r4z) / (r2 * r4)
  IF (arg > 1.0) arg = 1.0
  IF (arg < -1.0) arg = -1.0
  ! Calcul de l'angle entre le vecteur #2 et le vecteur #4.0
  tet24 = dacos(arg)

  a = tet23
  b = tet13
  c = tet12
  s = (a + b + c) / 2.0
  IF ((dtan(s / 2.0) * dtan((s - a) / 2.0) * dtan((s - b) / 2.0) * &
       dtan((s - c) / 2.0)) < 0.0) THEN
    a123 = 0.0
  ELSE
    ! Calcul de l'aire du triangle spherique borne par les vecteurs 1, 2 et 3.0
    a123 = 4.0 * datan(dsqrt(dtan(s / 2.0) * dtan((s - a) / 2.0) * dtan((s - b) / 2.0) * &
                             dtan((s - c) / 2.0)))
  END IF
!         alp=2.0*atan(sqrt(sin(s-b)*sin(s-c)/(sin(s)*sin(s-a))))         ! Autre methode pour calculer l'angle solide non utilisee.
!         bet=asin(sin(b)*sin(alp)/sin(a))                                ! Autre methode pour calculer l'angle solide non utilisee.
!         gam=asin(sin(c)*sin(alp)/sin(a))                                ! Autre methode pour calculer l'angle solide non utilisee.
!         a123=alp+bet+gam-pi                                             ! Autre methode pour calculer l'angle solide non utilisee.

  a = tet23
  b = tet34
  c = tet24
  s = (a + b + c) / 2.0

  IF ((dtan(s / 2.0) * dtan((s - a) / 2.0) * dtan((s - b) / 2.0) * &
       dtan((s - c) / 2.0)) < 0.0) THEN
    a234 = 0.0
  ELSE
    ! Calcul de l'aire du triangle spherique borne par les vecteurs 2, 3 et 4.0
    a234 = 4.0 * datan(dsqrt(dtan(s / 2.0) * dtan((s - a) / 2.0) * dtan((s - b) / 2.0) * &
                             dtan((s - c) / 2.0)))

  END IF
!         alp=2.0*atan(sqrt(sin(s-b)*sin(s-c)/(sin(s)*sin(s-a))))         ! Autre methode pour calculer l'angle solide non utilisee.
!         bet=asin(sin(b)*sin(alp)/sin(a))                                ! Autre methode pour calculer l'angle solide non utilisee.
!         gam=asin(sin(c)*sin(alp)/sin(a))                                ! Autre methode pour calculer l'angle solide non utilisee.
!         a234=alp+bet+gam-pi                                             ! Autre methode pour calculer l'angle solide non utilisee.

  ! L'angle solide est la somme des aires des deux triangles spheriques.
  omega = REAL(a123 + a234)

  RETURN
END
