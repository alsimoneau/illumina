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
 
 
subroutine zone_diffusion( &
              effet,zondif,ncell,stepdi,siz) 
    implicit none 
    integer :: i,j,k 
    integer :: ncell,neffet,imin,imax,jmin,jmax,kmin,kmax,nvol 
    integer :: keep,stepdi 
    real :: x0,y0,z0 
    real :: effet,dmin,d 
    real :: zondif(3000000,3),siz 
    real :: pi 
    pi = 3.141592654 
    keep = 0 
    neffet = nint(effet/siz) 
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
    do i = imin,imax 
      x0 = real(i)*siz 
      do j = jmin,jmax 
        y0 = real(j)*siz 
        do k = kmin,kmax 
          z0 = real(k)*siz 
          d = sqrt(x0**2.+y0**2.+z0**2.) 
          if (d <= dmin) then 
            keep = keep+1 
            if (keep == stepdi) then 
              keep = 0 
              ncell = ncell+1 
              zondif(ncell,1) = x0 
              zondif(ncell,2) = y0 
              zondif(ncell,3) = z0 
            end if 
          end if 
        end do 
      end do 
    end do 
    return 
end subroutine zone_diffusion 
