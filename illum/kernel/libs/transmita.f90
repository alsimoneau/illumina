!----------------------------------------------------------------------- 
 
!======================================================================= 
! Routine transmita (Andre Morin, Alex Neron, Etienne Rousseau 2004) 
 
! Determine la transmittance des aerosols sur un parcours entre les 
! cellules (x_i,y_i,z_i) et (x_f,y_f,z_f) 
! Est fonction de l'epaisseur optique des aerosols 
! Recoit des longueurs d'ondes en nanometre et les transforme en microns. 
! La pression doit etre en KPa 
! Retourne la transmittance transa 
 
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
 
 
subroutine transmita(angz,z_i,z_f,distd,haer,transa,tranaa) 
    real :: transa                                                         ! declaration des variables. 
    real :: tranaa                                                         ! vertical transmittance of the complete atmosphere (aerosols) 
    real :: angz,haer 
    real :: z_i,z_f,z1,z2 
    real :: distd 
 
 
    if (z_i > z_f) then 
      z2 = z_i 
      z1 = z_f 
    else 
      z1 = z_i 
      z2 = z_f 
    end if 
    if (z1 /= z2) then 
      transa = exp((log(tranaa)/abs(cos(angz)))*(exp(-1.*z1/haer)- &
                exp(-1.*z2/haer))) 
    else 
      transa = exp((log(tranaa))*exp(-1.*z1/haer)*distd/haer) 
    end if 
    if (distd == 0.) transa = 1. 
    if (transa == 0.) then 
      print*,'ERREUR transa - no transmission',z_i,z_f,angz,tranaa,distd &
                ,haer 
    end if 
    if (transa > 1.) then 
      print*,'ERREUR avec transa',transa,z_i,z_f,angz 
      stop 
    end if 
    return 
end subroutine transmita 
