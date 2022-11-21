!======================================================================= 
!  Routine anglezenithal (Andre Morin 2004) 
 
!  debugge par Martin Aube 2004 (cette routine ne calculait pas du tout 
!  l'angle zenithal 
 
!  Determine l'angle zenithal entre les points (x1,y1,z1) et (x2,y2,z2) 
!  Retourne l'angle angzen en radians 
 
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
 
 
subroutine anglezenithal(x1,y1,z1,x2,y2,z2,angzen) 
    real :: x1,y1,x2,y2 
    real :: z1,z2,pi,angzen 
    parameter (pi = 3.141592654) 
    hdist = sqrt((x2-x1)**2. &
              +(y2-y1)**2.) 
    if (z2-z1 /= 0.) then 
      angzen = atan(hdist/abs(z2-z1)) 
      if (z2-z1 < 0.) angzen = pi-angzen 
    else 
      angzen = pi/2. 
    end if 
    if ((angzen < 0.).or.(angzen > pi)) then 
      print*,'ERREUR angzen2 = ',angzen 
      stop 
    end if 
    return 
end subroutine anglezenithal 
