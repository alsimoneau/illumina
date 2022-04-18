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
 
 
subroutine horizon(x,y,z,dx,dy,altsol,anga,zhoriz,d) 
    integer :: width                                                       ! matrix dimension in length/width and height 
    parameter (width = 512) 
    integer :: x,y,nx,ny 
    real :: dx,dy,altsol(width,width),anga,zout,pi,angaz1,ix,iy 
    real :: hcur,distc                                                           ! earth curvature terrain 
    real :: posx,posy,scalef,zhoriz,z,d,dout 
    pi = 3.141592654 
    angaz1 = anga 
    ix = (cos(angaz1))                                                  ! viewing vector components 
    iy = (sin(angaz1)) 
    scalef = dx/3. 
    posx = real(x)*dx 
    posy = real(y)*dy 
    zhoriz = pi 
    d = (real(width))*sqrt(1.+tan(angaz1)**2.) 
    do while (((posx <= real(width)*dx).and.(posx >= 1.*dx)).and. &
                ((posy <= real(width)*dy).and.(posy >= 1.*dy))) 
      posx = posx+ix*scalef 
      posy = posy+iy*scalef 
      nx = nint(posx/dx) 
      ny = nint(posy/dy) 
! earth curvature (first order correction) 
      distc = sqrt((dx*real(nx-x))**2.+(dy*real(ny-y))**2.) 
      call curvature(distc,hcur) 
      if ((nx == x).and.(ny == y)) then                                  ! to forbid division by zero 
      if (z > altsol(nx,ny)-hcur) then                               ! reverse curvature to limit horizontal distance. curv is negative. this is a hack 
      zout = pi 
      d = 0. 
    else 
      zout = 0. 
      d = 0. 
    end if 
    else 
      dout = distc 
      zout = pi/2.-atan((altsol(nx,ny)-hcur-z)/dout) 
      if (altsol(nx,ny)-hcur == z) then 
        zout=pi/2.-0.0001*pi/180.                                      ! bug for zhoriz = pi, anyway in the real world pi is almost impossible 
      end if 
      if (zout < zhoriz) then 
        zhoriz = zout 
        d = dout 
      end if 
    end if 
    end do 
end subroutine horizon 
