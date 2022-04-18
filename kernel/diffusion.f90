!------------------------------------------------------------------------ 
 
!======================================================================= 
!  Routine diffusion 
 
! Determine the probability of light scattering per unit of solid angle 
! in the angdif direction. The scattering parameters are given 
! by secdif (ratio of the effective section of scattering on the section 
! of extinction), fonc_anorm (scattering phase function 
 
! Returns the scattering probability pdif 
 
! an Illumina routine 
!----------------------------------------------------------------------- 
 
!    Copyright (C) 2021  Martin Aube 
 
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
 
 
subroutine diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil, &
              fonc_a,fonc_l,haer,hlay,pdif,altit) 
    real :: angdif,pdif,prob_a,prob_m,prob_l,secdif,secdil 
    real :: fctmol,pi,fonc_a(181),fonc_l(181),fonc_ae,fonc_le 
    real :: angdeg,tranam,tranaa,tranal 
    real :: altit,un,hlay,haer 
    integer :: rang,na,naz 
    parameter (pi = 3.1415926) 
!-------------------------------------------------------- 
    if (angdif < 0.) angdif = -angdif 
    if (angdif-pi > 0.00001) angdif = pi 
    angdeg = ((angdif*180.)/pi) 
    rang = int(angdeg)+1 
!---------------------------------------- 
!  Calculate scattering probability per unit of steradian                 ! The probability is for a voxel of 1x1x1m refer to equation 1 in Aubé et al. 
 
!  Aubé, M., Simoneau, A., Muñoz-Tuñón, C., Díaz-Castro, J., & 
!  Serra-Ricart, M. (2020). Restoring the night sky darkness at 
!  Observatorio del Teide: First application of the model Illumina 
!  version 2. Monthly Notices of the Royal Astronomical Society, 
!  497(3), 2501-2516. 
!---------------------------------------- 
    if ((tranaa <= 1.).and.(tranaa > 0.)) then 
      fonc_ae = fonc_a(rang)                                             ! value of the aerosol phase function 
! Functions are normalized in the main code. See their division by 4pi 
      prob_a = (1.-exp(log(tranaa)*exp(-1.*altit/haer)*un/haer))* &
                secdif*fonc_ae 
    else 
      prob_a = 0. 
    end if 
    if ((tranal <= 1.).and.(tranal > 0.)) then 
      fonc_le = fonc_l(rang)                                             ! value of the layer phase function 
      prob_l = (1.-exp(log(tranal)*exp(-1.*altit/hlay)*un/hlay))* &
                secdil*fonc_le 
    else 
      prob_l = 0. 
    end if 
    fctmol = 0.75*(1.+((cos(angdif))**2.))/(4.*pi)                        ! value of the molecule phase function 
    prob_m = (1.-exp(log(tranam)*exp(-1.*altit/8000.)*un/8000.))* &
              fctmol 
 
    pdif = prob_a+prob_m+prob_l                                         ! This is an approximation valide if 1-transa, 
    ! 1-transm et 1-transl are small 
    if (prob_a > 1.) then 
      print*,'prob_a>1.' 
      stop 
    end if 
    if (prob_a < 0.) then 
      print*,'prob_a<0..' 
      stop 
    end if 
    if (prob_m > 1.) then 
      print*,'prob_m>1.' 
      stop 
    end if 
    if (prob_m < 0.) then 
      print*,'prob_m`¸^<0..' 
      stop 
    end if 
    if (prob_l > 1.) then 
      print*,'prob_l>1.' 
      stop 
    end if 
    if (prob_l < 0.) then 
      print*,'prob_l`¸^<0..' 
      stop 
    end if 
    if (pdif > 1.) then 
      pdif = 1. 
    end if 
    if (pdif < 0.) then 
      pdif = 0. 
    end if 
    return 
end subroutine diffusion 
