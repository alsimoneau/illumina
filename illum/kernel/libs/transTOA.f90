! calculation of the total vertical transmittance of the atmosphere
! (aerosols and molecules separately

! pour utilisation avec Illumina
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

SUBROUTINE transtoa(lambm, bandw, taua, layaod, pressi, tranam, tranaa, tranal, tabs)
  REAL :: m, lambm, taua, pressi, tranam, tranaa, tranal, layaod, tabs, bandw
  m = (pressi / 101.325)
!  transmittance taken from Bodhaine et al. (1999)
  tranam = EXP(-1.0 * m * 0.0021520 &
               * (1.0455996 - 341.29061 * (lambm / 1000.0)**-2.0 &
                  - 0.90230850 * (lambm / 1000.0)**2.0) &
               / (1.0 + 0.0027059889 * (lambm / 1000.0)**-2.0 &
                  - 85.968563 * (lambm / 1000.0)**2.0))
  CALL molabs(lambm, bandw, tabs)
  tranaa = EXP(-1.0 * taua)
  tranal = EXP(-1.0 * layaod)
  RETURN
END SUBROUTINE transtoa
