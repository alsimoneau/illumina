! programme pour faire la somme des flux dans un fichier .out de illumina 
 
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
 
 
 
    character*60 :: name &
              	 intger i &
              	 rea value, flux 
    open (unit=1,file='end do_illumina.in',status = 'old') 
    read(1,*) name 
    close (unit = 1) 
    name = name//'.out' 
    open(unit=2,file=name,status = 'old') &
              	   fux = 0. 
    do i = 1,1500 
      read(2,*,end = 10) bidon,bidon,bidon,bidon,value 
      flux = flux+value 
    end do 
 10   close (unit = 2) 
    open (unit=1,file='end do_illumina.out',status = 'unknown') 
    write(1,*) flux, 'Total flux for experiment ',name 
    close (unit = 1) &
              	 sto 
    	 end 
