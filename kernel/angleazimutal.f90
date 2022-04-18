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
 
 
subroutine angleazimutal(x1,y1,x2,y2,angazi) 
    real :: x1,y1,x2,y2 
    real :: pi,angazi 
    parameter (pi = 3.141592654) 
    if (x2-x1 /= 0.) angazi = abs(atan((y2-y1)/(x2-x1))) 
    if (((x2-x1) == 0.).and.((y2-y1) == 0.)) then 
      angazi = 0. 
    else 
      if (x2-x1 > 0.) then 
        if (y2-y1 < 0.) then 
          angazi = 2.*pi-angazi 
        end if 
      else if (x2-x1 < 0.) then 
        if (y2-y1 < 0.) then 
          angazi = angazi+pi 
        else 
          angazi = pi-angazi 
        end if 
      else                                                              ! i.e. x2-x1 = 0. 
        if (y2 > y1) angazi = pi/2. 
        if (y2 < y1) angazi = 3.*pi/2. 
      end if 
      if ((angazi < 0.).or.(angazi > 2.*pi)) then 
        print*,'ERREUR angazi = ',angazi,x1,y2,x2,y2 
        stop 
      end if 
!        if ((x1.eq.x2).and.(y1.eq.y2)) then 
!          print*,'ERREUR cant compute angle between identical points!' 
!          stop 
!        endif 
    end if 
    return 
end subroutine angleazimutal 
