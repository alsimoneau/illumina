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
 
subroutine cloudtransmitance(angzen,cloudt,tcloud) 
 
!======================================================================= 
!     Variables declaration 
!======================================================================= 
    integer :: cloudt 
!  fitted parameters for the cloud reflectance as a function of the incident zenith angle 
!  rho(z)=b0+b1*cos z + b2 * cos^2 z + b3 * cos^3 z according to Shapiro 1982 Table 11 
    real :: thocld(5,4) 
    real :: angzen,tcloud 
 
    thocld(1,1) = 0.63547 
    thocld(1,2) = 0.35229 
    thocld(1,3) = 0.08709 
    thocld(1,4) = -0.22902 
    thocld(2,1) = 0.26458 
    thocld(2,2) = 0.66829 
    thocld(2,3) = 0.24228 
    thocld(2,4) = -0.49357 
    thocld(3,1) = 0.19085 
    thocld(3,2) = 0.32817 
    thocld(3,3) = -0.08613 
    thocld(3,4) = -0.08197 
    thocld(4,1) = 0.13610 
    thocld(4,2) = 0.29964 
    thocld(4,3) = -0.14041 
    thocld(4,4) = 0.00952 
    thocld(5,1) = 0.17960 
    thocld(5,2) = 0.34855 
    thocld(5,3) = -0.14875 
    thocld(5,4) = 0.01962 
    tcloud = thocld(cloudt,1)+thocld(cloudt,2)*cos(angzen) &
              +thocld(cloudt,3)*(cos(angzen))**2.+thocld(cloudt,4)* &
              (cos(angzen))**3. 
!        print*,rcloud,angzen,cos(angzen) 
    return 
end subroutine cloudtransmitance 
 
