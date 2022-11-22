! **      http://cegepsherbrooke.qc.ca/~aubema/index.php/Prof/IllumEn?action=download&upname=intensite_lumineuse.pdf  ** 
! **                                                                                                                  ** 
! ********************************************************************************************************************** 
 
!    Copyright (C) 2015 Martin Aube 
 
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
 
 
!2345678901234567890123456789012345678901234567890123456789012345678901234 
 
subroutine cloudreflectance(angzen,cloudt,rcloud) 
 
!======================================================================= 
!     Variables declaration 
!======================================================================= 
    integer :: cloudt 
!  fitted parameters for the cloud reflectance as a function of the incident zenith angle 
!  rho(z)=a0+a1*cos z + a2 * cos^2 z + a3 * cos^3 z according to Shapiro 1982 Table 10 
! cloud no 1 are thin cirrus & cirrostratus 
! cloud no 2 are thick cirrus & cirrostratus 
! cloud no 3 are altostratus & altocumulus 
! cloud no 4 are stratocumulus & stratus 
! cloud no 5 are cumulus & cumulonimbus 
    real :: rhocld(5,4) 
    real :: angzen,rcloud 
 
    rhocld(1,1) = 0.25674 
    rhocld(1,2) = -0.18077 
    rhocld(1,3) = -0.21961 
    rhocld(1,4) = 0.25272 
    rhocld(2,1) = 0.60540 
    rhocld(2,2) = -0.55142 
    rhocld(2,3) = -0.23389 
    rhocld(2,4) = 0.43648 
    rhocld(3,1) = 0.66152 
    rhocld(3,2) = -0.14863 
    rhocld(3,3) = -0.08193 
    rhocld(3,4) = 0.13442 
    rhocld(4,1) = 0.71214 
    rhocld(4,2) = -0.15033 
    rhocld(4,3) = 0.00696 
    rhocld(4,4) = 0.03904 
    rhocld(5,1) = 0.67072 
    rhocld(5,2) = -0.13805 
    rhocld(5,3) = -0.10895 
    rhocld(5,4) = 0.09460 
    rcloud = rhocld(cloudt,1)+rhocld(cloudt,2)*cos(angzen) &
              +rhocld(cloudt,3)*(cos(angzen))**2.+rhocld(cloudt,4)* &
              (cos(angzen))**3. 
!        print*,rcloud,angzen,cos(angzen) 
    return 
end subroutine cloudreflectance 
 
