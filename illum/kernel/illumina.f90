!                         *                          *                       iii                   *     * 
!                                                                           iiiii 
!  IIIIII    lLLLL    *    lLLLL         UUU    UUU      MMMMM      MMMMM    iii        NNNN     NN          AAAA 
!   IIII     LLLL          LLLL   *     UUU      UUU     MMMMMMM  MMMMMMM          *    NNNNN    NN        AAAaaAAA 
!   IIII     LLLL          LLLL        UUU *      UUU    MMM MMMMMMMM MMM    iii        NNNNNN   NN       AAA    AAA 
!   IIII     LLLL   *      LLLL        UUU        UUU    MMM *        MMM  iii          NNN  NNN NN     AAAAAAAAAAAAAA 
!   IIII     LLLl          LLLl        UUUu      uUUU    MMM          MMM  iiii    ii   NNN   NNNNN    AAAa        aAAA 
!   IIII    LLLLLLLLLL    LLLLLLLLLL    UUUUUuuUUUUU     MMM          MMM   iiiiiiiii   NNN    NNNN   aAAA    *     AAAa 
!  IIIIII   LLLLLLLLLLL   LLLLLLLLLLL     UUUUUUUU      mMMMm        mMMMm   iiiiiii   nNNNn    NNNn  aAAA          AAAa 
 
! ********************************************************************************************************************** 
! ** Illumina VERSION 2 - in Fortran 77                                                                               ** 
! ** Programmers in decreasing order of contribution  :                                                               ** 
! **                            Martin Aube                                                                           ** 
! **              Still having very few traces of their contributions :                                               ** 
! **                            Loic Franchomme-Fosse,  Mathieu Provencher, Andre Morin                               ** 
! **                            Alex Neron, Etienne Rousseau                                                          ** 
! **                            William Desroches, Maxime Girardin, Tom Neron                                         ** 
! **                                                                                                                  ** 
! ** Illumina can be downloaded via:   git clone https://github.com/aubema/illumina.git                               ** 
! ** To compile:                                                                                                      ** 
! **    cd hg/illumina                                                                                                ** 
! **    mkdir bin                                                                                                     ** 
! **    bash makeILLUMINA                                                                                             ** 
! **                                                                                                                  ** 
! **  Current version features/limitations :                                                                          ** 
! **                                                                                                                  ** 
! **    - Calculation of artificial sky radiance in a given line of sight                                             ** 
! **    - Calculation of the atmospheric transmittance and 1st and 2nd order of scattering                            ** 
! **    - Lambertian reflexion on the ground                                                                          ** 
! **    - Terrain slope considered (apparent surface and shadows)                                                     ** 
! **    - Angular photometry of a lamp is considered uniform along the azimuth                                        ** 
! **    - Sub-grid obstacles considered (with the mean free path of light toward ground, mean obstacle height, and    ** 
! **      obstacles transparency (filling factor)                                                                     ** 
! **    - Molecules and aerosol optics (phase function, scattering probability, aerosol absorption)                   ** 
! **    - Exponential concentrations vertical profile                                                                 ** 
! **    - Accounting for heterogeneity luminaires number, luminaires heights, luminaire spectrum,                     ** 
! **      angular photometry, obstacle properties                                                                     ** 
! **    - Wavelength dependant                                                                                        ** 
! **    - Cloud models (type and cloud base height) only the overhead clouds are considered with cloud fraction       ** 
! **    - Support direct observation of a source                                                                      ** 
! **    - Direct observation of the ground is implemented                                                             ** 
! **                                                                                                                  ** 
! ********************************************************************************************************************** 
 
!  Copyright (C) 2021 Martin Aube PhD 
 
!  This program is free software: you can redistribute it and/or modify 
!  it under the terms of the GNU General Public License as published by 
!  the Free Software Foundation, either version 3 of the License, or 
!  (at your option) any later version. 
 
!  This program is distributed in the hope that it will be useful, 
!  but WITHOUT ANY WARRANTY; without even the implied warranty of 
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
!  GNU General Public License for more details. 
 
!  You should have received a copy of the GNU General Public License 
!  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 
!  Contact: martin.aube@cegepsherbrooke.qc.ca 
 
 
 
    program illumina                                                    ! Beginning 
    implicit none 
 
!======================================================================= 
!     Variables declaration 
!======================================================================= 
 
    integer :: width,nzon                                                  ! matrix dimension in length/width and height 
    parameter (width=512,nzon = 256) 
    integer :: iun,ideux 
    real :: pi,pix4 
    real :: zero,un                                                        ! value of 0. and 1. 
    integer :: verbose                                                     ! verbose = 1 to have more print out, 0 for silent 
    parameter (pi = 3.141592654) 
    parameter (pix4 = 4.*pi) 
    character*72 :: mnaf                                                   ! terrain elevation file 
    character*72 :: diffil                                                 ! aerosol file 
    character*72 :: outfile                                                ! results file 
    character*72 :: pclf,pclgp                                             ! files containing contribution and sensitivity maps 
    character*72 :: pclimg,pcwimg 
    character*72 :: basenm                                                 ! base name of files 
    integer :: lenbase                                                     ! length of the base name of the experiment 
    real :: lambda,pressi,drefle(width,width)                              ! wavelength (nanometer), atmospheric pressure (kpa), mean free path to the ground (meter). 
    real :: reflsiz                                                        ! size of the reflecting surface 
    integer :: ntype                                                       ! number of light source types or zones considered 
    real :: largx                                                          ! width (x axis) of the modeling domain (meter) 
    real :: largy                                                          ! length (y axis) of the modeling domain (meter) 
    integer :: nbx,nby                                                     ! number of pixels in the modeling domain 
    real :: val2d(width,width)                                             ! temporary input array 2d 
    real :: altsol(width,width)                                            ! ground elevation (meter) 
    real :: srefl                                                          ! ground reflectance 
    integer :: stype                                                       ! source type or zone index 
    character*72 :: pafile,lufile,alfile,ohfile,odfile,offile              ! files related to light sources and obstacles (photometric function of the sources (sr-1), flux (w), height (m), obstacles c                                                               ! height (m), obstacle distance (m), obstacle filling factor (0-1). 
    real :: lamplu(width,width,nzon)                                       ! source fluxes 
    real :: lampal(width,width)                                            ! height of the light sources relative to the ground (meter) 
    real :: pval(181,nzon),pvalto,pvalno(181,nzon)                         ! values of the angular photometry functions (unnormalized, integral, normalized) 
    real :: dtheta                                                         ! angle increment of the photometric function of the sources 
    real :: dx,dy,dxp,dyp                                                  ! width of the voxel (meter) 
    integer :: boxx,boxy                                                   ! reflection window size (pixels) 
    real :: fdifa(181),fdifan(181)                                         ! aerosol scattering functions (unnormalized and normalized) 
    real :: extinc,scatte,anglea(181)                                      ! aerosol cross sections (extinction and scattering), scattering angle (degree) 
    real :: secdif                                                         ! contribution of the scattering to the extinction 
    real :: inclix(width,width)                                            ! tilt of the ground pixel along x (radian) 
    real :: incliy(width,width)                                            ! tilt of the ground pixel along y (radian) 
    integer :: x_obs,y_obs                                                 ! position of the observer (integer) 
    real :: rx_obs,ry_obs 
    real :: z_o                                                            ! observer height relative to the ground (meter) 
    real :: z_obs                                                          ! height of the observer (meter) to the vertical grid scale 
    integer :: ncible,icible                                               ! number of line of sight voxels, number loops over the voxels 
    integer :: x_c,y_c                                                     ! position of the line of sight voxel (integer) 
    real :: rx_c,ry_c 
    real :: z_c                                                            ! height of the line of sight voxel (meter) 
    integer :: dirck                                                       ! test for the position of the source (case source = line of sight voxel) 
    integer :: x_s,y_s,x_sr,y_sr,x_dif,y_dif,zceldi                        ! positions of the source, the reflecting surface, and the scattering voxels 
    real :: z_s,z_sr,z_dif                                                 ! heights of the source, the reflecting surface, and the scattering voxel (metre). 
    real :: rx_s,ry_s,rx_sr,ry_sr,rx_dif,ry_dif 
    real :: angzen,ouvang                                                  ! zenithal angle between two voxels (radians) and opening angle of the solid angle in degrees. 
    integer :: anglez                                                      ! emitting zenithal angle from the luminaire. 
    real :: p_dir,p_indir,p_dif1                                           ! photometric function of the light sources (direct,indirect,scattered) 
    real :: transa,transm,transl                                           ! transmittance between two voxels (aerosols,molecules,particle layer). 
    real :: tran1a,tran1m                                                  ! transmittance of the voxel (aerosols,molecules). 
    real :: taua                                                           ! aerosol optical depth @ 500nm. 
    real :: alpha                                                          ! angstrom coefficient of aerosol aod 
    real(8) xc,yc,zc,xn,yn,zn                                            ! position (meter) of the elements (starting point, final point) :: for the calculation of the solid angle. 
    real(8) :: r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! components of the vectors used in the solid angle calculation routine. 
    real :: omega,omega1                                                   ! solid angles 
    real :: fldir                                                          ! flux coming from a source (watt). 
    real :: flindi                                                         ! flux coming from a reflecting ground element (watt). 
    real :: fldiff                                                         ! flux coming from a scattering voxel (watt). 
    real :: zidif,zfdif                                                    ! initial and final limits of a scattering path. 
    real :: angdif                                                         ! scattering angle. 
    real :: pdifdi,pdifin,pdifd1,pdifd2                                    ! scattering probability (direct,indirect,1st and 2nd order of scattering 
    real :: intdir                                                         ! direct intensity toward the sensor from a scattering voxel. 
    real :: intind                                                         ! contribution of the reflecting cell to the reflected intensity toward the sensor. 
    real :: itotind                                                        ! total contribution of the source to the reflected intensity toward the sensor. 
    real :: idiff2                                                         ! contribution of the scattering voxel to the scattered intensity toward the sensor. 
    real :: itodif                                                         ! total contribution of the source to the scattered intensity toward the sensor. 
    real :: isourc                                                         ! total contribution of the source to the intensity from a line of sight voxel toward the sensor. 
    real :: itotty                                                         ! total contribution of a source type to the intensity coming from a line of sight voxel toward the sensor. 
    real :: itotci                                                         ! total intensity from a line of sight voxel toward the sensor. 
    real :: itotrd                                                         ! total intensity a voxel toward the sensor after reflexion and double scattering. 
    real :: flcib                                                          ! flux reaching the observer voxel from a line of sight voxel. 
    real :: fcapt                                                          ! flux reaching the observer voxel from all fov voxels in a given model level 
    real :: ftocap                                                         ! total flux reaching the observer voxel 
    real :: haut                                                           ! haut (negative indicate that the surface is lighted from inside the ground. i.e. not considered in the calculation 
    real :: epsilx,epsily                                                  ! tilt of the ground pixel 
    real :: flrefl                                                         ! flux reaching a reflecting surface (watts). 
    real :: irefl,irefl1                                                   ! intensity leaving a reflecting surface toward the line of sight voxel. 
    real :: effdif                                                         ! distance around the source voxel and line of sight voxel considered to compute the 2nd order of scattering. 
    real :: zondif(3000000,3)                                              ! array for the scattering voxels positions 
    integer :: ndiff,idi                                                   ! number of scattering voxels, counter of the loop over the scattering voxels 
    integer :: stepdi                                                      ! scattering step to speedup the calculation e.g. if =2 one computation over two will be done 
    integer :: ssswit                                                      ! activate double scattering (1=yes, 0 = no) 
    integer :: fsswit                                                      ! activate first scattering (1=yes, 0 = no) 
    integer :: nvis0                                                       ! starting value for the calculation along of the viewing line. 
!                                                                         ! by default the value is 1 but it can be larger 
!                                                                         ! when we resume a previous interrupted calculation. 
    real :: fldif1,fldif2                                                  ! flux reaching a scattering voxel. 
    real :: fdif2                                                          ! flux reaching the line of sight voxel after reflexion > scattering 
    real :: idif1,idif2,idif2p                                             ! intensity toward a line of sight voxel from a scattering voxel (without and with reflexion). 
    real :: portio                                                         ! ratio of voxel surface to the solid angle of the sensor field of view. 
    real :: dis_obs                                                        ! distance between the line of sight and the observer. 
    real :: ometif                                                         ! solid angle of the telescope objective as seen from the line of sight voxel 
    real :: omefov                                                         ! solid angle of the spectrometer slit. 
    real :: angvis,azim                                                    ! viewing angles of the sensor. 
!                                                                         ! Useful for the calculation of the lambertian reflectance. 
    real :: nbang                                                          ! for the averaging of the photometric function 
    real :: obsh(width,width),angmin                                       ! averaged height of the sub-grid obstacles, minimum angle under wich 
!                                                                         ! a light ray cannot propagate because it is blocked by a sub-grid obstable 
    real :: ofill(width,width)                                             ! fill factor giving the probability to hit an obstacle when pointing in its direction real 0-1 
    integer :: naz,na 
    real :: itt(width,width,nzon)                                          ! total intensity per type of lamp 
    real :: itc(width,width)                                               ! total intensity per line of sight voxel 
    real :: ftc(width,width)                                               ! fraction of the total flux at the sensor level 
    real :: fca(width,width)                                               ! sensor flux array 
    real :: lpluto(width,width)                                            ! total luminosity of the ground cell for all lamps 
    character*3 :: lampno                                                  ! lamp number string 
    integer :: imin(nzon),imax(nzon),jmin(nzon),jmax(nzon)                 ! x and y limits containing a type of lamp 
    real :: angazi                                                         ! azimuth angle between two points in rad, max dist for the horizon determination 
    real :: latitu                                                         ! approximate latitude of the domain center 
    integer :: prmaps                                                      ! flag to enable the tracking of contribution and sensitivity maps 
    integer :: cloudt                                                      ! cloud type 0=clear, 1=thin cirrus/cirrostratus, 2=thick cirrus/cirrostratus, 3 = altostratus/altocumulus, 
    ! 4=Stratocumulus/stratus, 5 = Cumulus/Cumulonimbus 
    real :: cloudslope                                                     ! slope of the radiance dependency on the cloud fraction (in percentage) according to 
    ! Sciezoras 2020 the slope vary depending on the level of LP and how it is distributed. 
    ! We decided instead to simplify this by using an average slope of -0.013. 
    ! Rad = Rad_100 * 10**(0.4*(100-cloudfrac)*cloudslope) this equation is derived from 
    ! Tomasz Sciezor, The impact of clouds on the brightness of the night sky, Journal of 
    ! Quantitative Spectroscopy & Radiative Transfer (2020), 
    ! doi: https://doi.org/10.1016/j.jqsrt.2020.106962 
    real :: cloudfrac                                                      ! cloud fraction in percentage 
    integer :: xsrmi,xsrma,ysrmi,ysrma                                     ! limits of the loop valeur for the reflecting surfaces 
    real :: rcloud                                                         ! cloud relfectance 
    real :: azencl                                                         ! zenith angle from cloud to observer 
    real :: icloud                                                         ! cloud reflected intensity 
    real :: fcloud                                                         ! flux reaching the intrument from the cloud voxel 
    real :: fccld                                                          ! correction for the fov to the flux reaching the intrument from the cloud voxel 
    real :: fctcld                                                         ! total flux from cloud at the sensor level 
    real :: totlu(nzon)                                                    ! total flux of a source type 
    real :: stoplim                                                        ! stop computation when the new voxel contribution is less than 1/stoplim of the cumulated flux 
    real :: ff,ff2,hh                                                          ! temporary obstacle filling factor and horizon blocking factor 
    real :: cloudbase,cloudtop,cloudhei                                    ! cloud base and top altitude (m), cloud layer avg height (m) 
    real :: distd                                                          ! distance to compute the scattering probability 
    real :: volu                                                           ! volume of a voxel 
    real :: scal                                                           ! stepping along the line of sight 
    real :: scalo                                                          ! previous value of scal 
    real :: siz                                                            ! resolution of the 2nd scat grid in meter 
    real :: angvi1,angaz1,angze1                                           ! viewing angles in radian 
    real :: ix,iy,iz                                                       ! base vector of the viewing (length = 1) 
    real :: dsc2,doc2                                                      ! square of the path lengths for the cloud contribution 
    real :: azcl1,azcl2                                                    ! zenith angle from the (source, refl surface, or scattering voxel) to line of path and observer to line p. 
    real :: dh,dho                                                         ! distance of the horizon limit 
    integer :: n2nd                                                        ! desired number of voxel in the calculation of the 2nd scattering 
    integer :: step                                                        ! skiping 2nd scat on 1 dim 
    real :: omemax                                                         ! max solid angle allowed 
    real :: tcloud                                                         ! low cloud transmission 
    real :: rx_sp,ry_sp                                                    ! position of a low cloud pixel 
    real :: flcld(width,width)                                             ! flux crossing a low cloud 
    real :: ds1,ds2,ds3,dss                                                ! double scattering distances 
    integer :: nss                                                         ! number of skipped 2nd scat elements 
    integer :: ndi                                                         ! number of cell under ground 
    integer :: nvol                                                        ! number of cell for second scat calc un full resolution 
    real :: diamobj                                                        ! instrument objective diameter 
    integer :: i,j,k,id,jd 
    real :: tranam,tranaa                                                  ! atmospheric transmittancess of a path (molecular, aerosol) 
    real :: zhoriz                                                         ! zenith angle of the horizon 
    real :: direct                                                         ! direct radiance from sources on a surface normal to the line of sight (no scattering) 
    real :: rdirect                                                        ! direct radiance from a reflecting surface on a surface normal to the line of sight (no scattering) 
    real :: irdirect                                                       ! direct irradiance from sources on a surface normal to the line of sight (no scattering) 
    real :: irrdirect                                                      ! direct irradiance from a reflecting surface on a surface normal to the line of sight (no scattering) 
    real :: dang                                                           ! angle between the line of sight and the direction of a source 
    real :: dzen                                                           ! zenith angle of the source-observer line 
    real :: ddir_obs                                                       ! distance between the source and the observer 
    real :: rx,ry,rz                                                       ! driving vector for the calculation of the projection angle for direct radiance. it is 20km long 
    real :: dfov                                                           ! field of view in degrees for the calculation of the direct radiance this number will be a kind of smoothing effect. the angular grid resolution to create a direct radiance panorama should be finer than that number 
    real :: fo                                                             ! flux correction factor for obstacles 
    real :: thetali                                                        ! limit angle for the obstacles blocking of viirs 
    integer :: viirs(width,width)                                          ! viirs flag 1=yes 0 = no 
    character*72 :: vifile                                                 ! name of the viirs flag file 
    real :: dh0,dhmax                                                      ! horizontal distance along the line of sight and maximum distance before beeing blocked by topography 
    character*72 :: layfile                                                ! filename of the optical properties of the particle layer 
    real :: layaod                                                         ! 500 nm aod of the particle layer 
    real :: layalp                                                         ! spectral exponent of the aod for the particle layer 
    real :: hlay                                                           ! exponential vertical scale height of the particle layer 
    real :: secdil                                                         ! scattering/extinction ratio for the particle layer 
    real :: fdifl(181)                                                     ! scattering phase function of the particle layer 
    real :: tranal                                                         ! top of atmos transmission of the particle layer 
    real :: haer                                                           ! exponential vertical scale height of the background aerosol layer 
    real :: distc,hcur                                                     ! distance to any cell and curvature  correction for the earth curvature 
    real :: bandw                                                          ! bandwidth of the spectral bin 
    real :: tabs                                                           ! toa transmittance related to molecule absorption 
    integer :: obsobs                                                      ! flag to activate the direct light obstacle blocking aroud the observer. 
    verbose=1                                                           ! Very little printout=0, Many printout = 1, even more = 2 
    diamobj = 1.                                                          ! A dummy value for the diameter of the objective of the instrument used by the observer. 
    volu = 0. 
    zero = 0. 
    un = 1. 
    ff = 0. 
    ff2 = 0. 
    step = 1 
    ncible = 1024 
    stepdi = 1 
    cloudslope = -0.013 
    cloudfrac = 100. 
    if (verbose >= 1) then 
      print*,'Starting ILLUMINA computations...' 
    end if 
! reading of the fichier d'entree (illumina.in) 
    print*,'Reading illumina.in input file' 
    open(unit=1,file='illumina.in',status = 'old') 
    read(1,*) 
    read(1,*) basenm 
    read(1,*) dx,dy 
    read(1,*) diffil 
    read(1,*) layfile, layaod, layalp, hlay 
    read(1,*) ssswit 
    read(1,*) fsswit 
    read(1,*) lambda,bandw 
    read(1,*) srefl 
    read(1,*) pressi 
    read(1,*) taua,alpha,haer 
    read(1,*) ntype 
    read(1,*) stoplim 
    read(1,*) 
    read(1,*) x_obs,y_obs,z_o 
    read(1,*) obsobs 
    read(1,*) angvis,azim 
    read(1,*) dfov 
    read(1,*) 
    read(1,*) 
    read(1,*) 
    read(1,*) reflsiz 
    read(1,*) cloudt, cloudbase, cloudfrac 
    read(1,*) 
    close(1) 
    if (angvis > 90.) then 
      print*,'Error: elevation angle larger than 90 deg' 
      stop 
    end if 
    if (angvis < -90.) then 
      print*,'Error: elevation angle smaller than -90 deg' 
      stop 
    end if 
! conversion of the geographical viewing angles toward the cartesian 
! angle we assume that the angle in the file illumina.in 
! is consistent with the geographical definition 
! geographical, azim=0 toward north, 90 toward east, 180 toward south 
! etc 
! cartesian, azim=0 toward east, 90 toward north, 180 toward west etc 
    azim = 90.-azim 
    if (azim < 0.) azim = azim+360. 
    if (azim >= 360.) azim = azim-360. 
    angvi1 = (pi*angvis)/180. 
    angze1 = pi/2.-angvi1 
    angaz1 = (pi*azim)/180. 
    ix = ( sin((pi/2.)-angvi1) ) * (cos(angaz1))                        ! viewing vector components 
    iy = ( sin((pi/2.)-angvi1) ) * (sin(angaz1)) 
    iz = (sin(angvi1)) 
    dfov = (dfov*pi/180.)/2. 
    siz = 2500. 
    if (ssswit == 0) then 
      effdif = 0. 
    else 
      effdif = 40000. 
    end if 
    scal = 19. 
    scalo = scal 
    boxx = nint(reflsiz/dx)                                               ! Number of column to consider left/right of the source for the reflection. 
    boxy = nint(reflsiz/dy)                                               ! Number of column to consider up/down of the source for the reflection. 
! omemax: exclude calculations too close (<10m) this is a sustended angle of 1 deg. 
! the calculated flux is highly sensitive to that number for a very high 
! pixel resolution (a few 10th of meters). We assume anyway that somebody 
! observing the sky will never lies closer than that distance to a 
! light fixture. This number is however somehow subjective and that means 
! that the value of sky brightness near sources will be affected by this 
! choice 
    omemax = 1./((10.)**2.) 
    if (verbose > 0) then 
      print*,'2nd order scattering grid = ',siz,'m' 
      print*,'2nd order scattering radius = ',effdif,'m' 
      print*,'Pixel size = ',dx,' x ',dy 
      print*,'Maximum radius for reflection = ',reflsiz 
    end if 
! computing the actual AOD at the wavelength lambda 
    if (verbose >= 1) print*,'500nm aod = ',taua,'500nm angstrom coeff.= &
              ',alpha 
    taua = taua*(lambda/500.)**(-1.*alpha) 
    layaod = layaod*(lambda/500.)**(-1.*layalp) 
!  determine the Length of basenm 
    lenbase = index(basenm,' ')-1 
    mnaf = basenm(1:lenbase)//'_topogra.bin'                              ! determine the names of input and output files 
    outfile = basenm(1:lenbase)//'.out' 
    pclf = basenm(1:lenbase)//'_pcl.txt' 
    pclimg = basenm(1:lenbase)//'_pcl.bin' 
    pcwimg = basenm(1:lenbase)//'_pcw.bin' 
    pclgp = basenm(1:lenbase)//'_pcl.gplot' 
! opening output file 
    open(unit=2,file=outfile,status = 'unknown') 
    write(2,*) "ILLUMINA version __version__" 
    write(2,*) 'FILE USED:' 
    write(2,*) mnaf,diffil 
    print*,'Wavelength (nm):',lambda, &
              ' Aerosol optical depth:',taua 
    write(2,*) 'Wavelength (nm):',lambda, &
              ' Aerosol optical depth:',taua 
    write(2,*) '2nd order scattering radius:',effdif,' m' 
    print*,'2nd order scattering radius:',effdif,' m' 
    write(2,*) 'Observer position (x,y,z)',x_obs,y_obs,z_o 
    print*,'Observer position (x,y,z)',x_obs,y_obs,z_o 
    write(2,*) 'Elevation angle:',angvis,' azim angle (counterclockwise &
              from east)',azim 
    print*,'Elevation angle:',angvis,' azim angle (counterclockwise &
              from east)',azim 
! Initialisation of the arrays and variables 
    if (verbose >= 1) print*,'initializing variables...' 
    if (cloudt == 0) then 
      cloudbase = 1000000000. 
    end if 
    prmaps = 1 
    iun = 0 
    ideux = 1 
    icloud = 0. 
    do i = 1,width 
      do j = 1,width 
        val2d(i,j) = 0. 
        altsol(i,j) = 0. 
        obsH(i,j) = 0. 
        viirs(i,j) = 0 
        ofill(i,j) = 0. 
        inclix(i,j) = 0. 
        incliy(i,j) = 0. 
        lpluto(i,j) = 0. 
        ITC(i,j) = 0. 
        FTC(i,j) = 0. 
        FCA(i,j) = 0. 
        flcld(i,j) = 0. 
        do k = 1,nzon 
          lamplu(i,j,k) = 0. 
          lampal(i,j) = 0. 
          ITT(i,j,k) = 0. 
        end do 
      end do 
    end do 
    do i = 1,181 
      fdifa(i) = 0. 
      fdifan(i) = 0. 
      fdifl(i) = 0. 
      anglea(i) = 0. 
      do j = 1,nzon 
        pval(i,j) = 0. 
        pvalno(i,j) = 0. 
      end do 
    end do 
    do i = 1,3000000 
      do j = 1,3 
        zondif(i,j) = 1. 
      end do 
    end do 
    idif1 = 0. 
    idif2 = 0. 
    fdif2 = 0 
    idif2p = 0. 
    fldir = 0. 
    flindi = 0. 
    fldiff = 0. 
    pdifdi = 0. 
    pdifin = 0. 
    pdifd1 = 0. 
    pdifd2 = 0. 
    intdir = 0. 
    intind = 0. 
    idiff2 = 0. 
    angmin = 0. 
    isourc = 0. 
    itotty = 0. 
    itotci = 0. 
    itotrd = 0. 
    flcib = 0. 
    flrefl = 0. 
    irefl = 0. 
    irefl1 = 0. 
    fldif1 = 0. 
    fldif2 = 0. 
    portio = 0. 
    fccld = 0. 
    fctcld = 0. 
    ometif = 0. 
    omefov = 0. 
    hh = 1. 
! determine the 2nd scattering zone 
    if (ssswit /= 0) then 
      call zone_diffusion(effdif, &
                zondif,ndiff,stepdi,siz) 
      dss = 1.*siz 
      if (verbose > 0) then 
        print*,'2nd order scattering grid points =',ndiff 
        print*,'2nd order scattering smoothing radius =',dss,'m' 
      end if 
    end if 
! determination of the vertical atmospheric transmittance 
! tranam and tranaa are the top of atmosphere transmittance (molecules and aerosols) 
    call transtoa(lambda,bandw,taua,layaod,pressi,tranam,tranaa, &
              tranal,tabs) 
 
! reading of the environment variables 
! reading of the elevation file 
    call twodin(nbx,nby,mnaf,altsol) 
! computation of the tilt of the pixels along x and along y 
    do i = 1,nbx                                                        ! beginning of the loop over the column (longitude) of the domain. 
      do j = 1,nby                                                      ! beginning of the loop over the rows (latitu) of the domain. 
        if (i == 1) then                                              ! specific case close to the border of the domain (vertical side left). 
        inclix(i,j) = atan((altsol(i+1,j)-altsol(i,j))/real(dx))      ! computation of the tilt along x of the surface. 
      else if (i == nbx) then                                        ! specific case close to the border of the domain (vertical side right). 
        inclix(i,j) = atan((altsol(i-1,j)-altsol(i,j))/(real(dx)))    ! computation of the tilt along x of the surface. 
      else 
! computation of the tilt along x of the surface. 
        inclix(i,j) = atan((altsol(i+1,j)-altsol(i-1,j))/(2. &
                  *real(dx))) 
      end if 
      if (j == 1) then                                              ! specific case close to the border of the domain (horizontal side down). 
      incliy(i,j) = atan((altsol(i,j+1)-altsol(i,j))/(real(dy)))    ! computation of the tilt along y of the surface. 
    else if (j == nby) then                                        ! specific case close to the border of the domain (horizontal side up). 
      incliy(i,j) = atan((altsol(i,j-1)-altsol(i,j))/(real(dy)))    ! computation of the tilt along y of the surface. 
    else 
! computation of the tilt along y of the surface 
      incliy(i,j) = atan((altsol(i,j+1)-altsol(i,j-1))/(2. &
                *real(dy))) 
    end if 
    end do                                                           ! end of the loop over the rows (latitu) of the domain 
    end do                                                             ! end of the loop over the column (longitude) of the domain 
! reading of the values of P(theta), height, luminosities and positions 
! of the sources, obstacle height and distance 
    ohfile = basenm(1:lenbase)//'_obsth.bin' 
    odfile = basenm(1:lenbase)//'_obstd.bin' 
    alfile = basenm(1:lenbase)//'_altlp.bin'                            ! setting the file name of height of the sources lumineuse. 
    offile = basenm(1:lenbase)//'_obstf.bin' 
    vifile = 'origin.bin' 
    dtheta = .017453293                                                 ! one degree 
! reading lamp heights 
    call twodin(nbx,nby,alfile,val2d) 
    do i = 1,nbx                                                        ! beginning of the loop over all cells along x. 
      do j = 1,nby                                                      ! beginning of the loop over all cells along y. 
        lampal(i,j) = val2d(i,j)                                        ! filling of the array for the lamp stype 
      end do                                                           ! end of the loop over all cells along y. 
    end do                                                             ! end of the loop over all cells along x. 
! reading subgrid obstacles average height 
    call twodin(nbx,nby,ohfile,val2d) 
    do i = 1,nbx                                                        ! beginning of the loop over all cells along x. 
      do j = 1,nby                                                      ! beginning of the loop over all cells along y. 
        obsH(i,j) = val2d(i,j)                                          ! filling of the array 
      end do                                                           ! end of the loop over all cells along y. 
    end do 
! reading subgrid obstacles average distance 
    call twodin(nbx,nby,odfile,val2d) 
    do i = 1,nbx                                                        ! beginning of the loop over all cells along x. 
      do j = 1,nby                                                      ! beginning of the loop over all cells along y. 
        drefle(i,j) = val2d(i,j)/2. 
        if (drefle(i,j) == 0.) drefle(i,j) = dx                         ! when outside a zone, block to the size of the cell (typically 1km) 
      end do                                                           ! end of the loop over all cells along y. 
    end do 
! reading subgrid obstacles filling factor 
    call twodin(nbx,nby,offile,val2d) 
    do i = 1,nbx                                                        ! beginning of the loop over all cells along x. 
      do j = 1,nby                                                      ! beginning of the loop over all cells along y. 
        ofill(i,j) = val2d(i,j)                                         ! Filling of the array 0-1 
      end do                                                           ! end of the loop over all cells along y. 
    end do 
! reading viirs flag 
    call twodin(nbx,nby,vifile,val2d) 
    do i = 1,nbx                                                        ! beginning of the loop over all cells along x. 
      do j = 1,nby                                                      ! beginning of the loop over all cells along y. 
        viirs(i,j) = nint(val2d(i,j))                                   ! viirs flag array 0 or 1 
      end do                                                           ! end of the loop over all cells along y. 
    end do 
! reading of the scattering parameters for background aerosols 
    open(unit = 1, file = diffil,status= 'old')                       ! opening file containing the scattering parameters 
    read(1,*)  secdif                                               ! the scattering / extinction ratio 
    read(1,*) 
    do i = 1,181 
      read(1,*) anglea(i), fdifa(i)                                 ! reading of the scattering functions 
      fdifan(i) = fdifa(i)/pix4                                       ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle 
    end do 
    close(1) 
! reading scattering parameters of particle layer 
    open(unit = 1, file = layfile,status= 'old')                     ! opening file containing the scattering parameters 
    read(1,*)  secdil                                               ! the scattering / extinction ratio of particle layer 
    read(1,*) 
    do i = 1,181 
      read(1,*) anglea(i), fdifl(i)                                 ! reading of the scattering functions of the particle layer 
      fdifl(i) = fdifl(i)/pix4                                        ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle 
    end do 
    close(1) 
! Some preliminary tasks 
    do stype = 1,ntype                                                  ! beginning of the loop 1 for the nzon types of sources. 
      imin(stype) = nbx 
      jmin(stype) = nby 
      imax(stype) = 1 
      jmax(stype) = 1 
      pvalto = 0. 
      write(lampno, '(I3.3)' ) stype                                  ! support of nzon different sources (3 digits) 
      pafile = basenm(1:lenbase)//'_fctem_'//lampno//'.dat'             ! setting the file name of angular photometry. 
      lufile = basenm(1:lenbase)//'_lumlp_'//lampno//'.bin'             ! setting the file name of the luminosite of the cases. 
! reading photometry files 
      open(UNIT=1, FILE=pafile,status = 'OLD')                          ! opening file pa#.dat, angular photometry. 
      do i = 1,181                                                    ! beginning of the loop for the 181 data points 
        read(1,*) pval(i,stype)                                     ! reading of the data in the array pval. 
! Sum of the values of the  photometric function 
        pvalto = pvalto+pval(i,stype)*2.*pi* &
                  sin(real(i-1)*dtheta)*dtheta                                ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale (i-1) x 1 degrees). 
      end do                                                         ! end of the loop over the 181 donnees of the fichier pa#.dat. 
      close(1)                                                        ! closing file pa#.dat, angular photometry. 
      do i = 1,181 
        if (pvalto /= 0.) pvalno(i,stype) = pval(i,stype)/pvalto        ! normalisation of the photometric function. 
      end do 
! reading luminosity files 
      call twodin(nbx,nby,lufile,val2d) 
      do i = 1,nbx                                                      ! beginning of the loop over all cells along x. 
        do j = 1,nby                                                    ! beginning of the loop over all cells along y. 
          if (val2d(i,j) < 0.) then                                  ! searching of negative fluxes 
          print*,'***Negative lamp flux!, stopping execution' 
          stop 
        end if 
      end do                                                         ! end of the loop over all cells along y. 
    end do 
    do i = 1,nbx                                                      ! searching of the smallest rectangle containing the zone 
      do j = 1,nby                                                    ! of non-null luminosity to speedup the calculation 
        if (val2d(i,j) /= 0.) then 
          if (i-1 < imin(stype)) imin(stype) = i-2 
          if (imin(stype) < 1) imin(stype) = 1 
          goto 333 
        end if 
      end do 
    end do 
    imin(stype) = 1 
    333    do i = nbx,1,-1 
      do j = 1,nby 
        if (val2d(i,j) /= 0.) then 
          if (i+1 > imax(stype)) imax(stype) = i+2 
          if (imax(stype) > nbx) imax(stype) = nbx 
          goto 334 
        end if 
      end do 
    end do 
    imax(stype) = 1 
    334    do j = 1,nby 
      do i = 1,nbx 
        if (val2d(i,j) /= 0.) then 
          if (j-1 < jmin(stype)) jmin(stype) = j-2 
          if (jmin(stype) < 1) jmin(stype) = 1 
          goto 335 
        end if 
      end do 
    end do 
    jmin(stype) = 1 
    335    do j = nby,1,-1 
      do i = 1,nbx 
        if (val2d(i,j) /= 0.) then 
          if (j+1 > jmax(stype)) jmax(stype) = j+2 
          if (jmax(stype) > nby) jmax(stype) = nby 
          goto 336 
        end if 
      end do 
    end do 
    jmax(stype) = 1 
    336    do i = 1,nbx 
! beginning of the loop over all cells along x. 
      do j = 1,nby 
! beginning of the loop over all cells along y. 
        lamplu(i,j,stype) = val2d(i,j) 
! remplir the array of the lamp type: stype 
! Atmospheric correction and obstacles masking corrections to the lamp 
! flux arrays (lumlp) 
        if (viirs(i,j) == 1) then 
          lamplu(i,j,stype) = lamplu(i,j,stype)/(tranam*tranaa* &
                    tranal) 
          thetali = atan2(drefle(i,j),obsH(i,j)) 
          if (thetali  <  70.*pi/180.) then 
            Fo = (1.-cos(70.*pi/180.))/(1.-ofill(i,j)*cos(thetali)+ &
                      (ofill(i,j)-1.)*cos(70.*pi/180.)) 
            lamplu(i,j,stype) = lamplu(i,j,stype)*Fo 
          else 
            Fo = 1. 
          end if 
        end if 
        totlu(stype) = totlu(stype)+lamplu(i,j,stype) 
! the total lamp flux should be non-null to proceed to the calculations 
      end do 
! end of the loop over all cells along y. 
    end do 
! end of the loop over all cells along x. 
    end do 
! end of the loop 1 over the nzon types of sources. 
    dy = dx 
    omefov = 0.00000001 
! solid angle of the spectrometer slit on the sky. Here we only need a small value 
    z_obs = z_o+altsol(x_obs,y_obs) 
! _obs = the local observer elevation plus the height of observation above ground (z_o) 
    rx_obs = real(x_obs)*dx 
    ry_obs = real(y_obs)*dy 
    if (z_obs == 0.) z_obs = 0.001 
    largx = dx*real(nbx) 
! computation of the Width along x of the case. 
    largy = dy*real(nby) 
! computation of the Width along y of the case. 
    write(2,*) 'Width of the domain [NS](m):',largx,'#cases:',nbx 
    write(2,*) 'Width of the domain [EO](m):',largy,'#cases:',nby 
    write(2,*) 'Size of a cell (m):',dx,' X ',dy 
    write(2,*) 'latitu center:',latitu 
 
 
 
 
 
 
    direct = 0. 
! initialize the total direct radiance from sources to observer 
    rdirect = 0. 
! initialize the total reflected radiance from surface to observer 
    irdirect = 0. 
! initialize the total direct irradiance from sources to observer 
    irrdirect = 0. 
! initialize the total reflected irradiance from surface to observer 
! ================================= 
! Calculation of the direct radiances 
 
    if (verbose >= 1) print*,' calculating obtrusive light...' 
    do stype = 1,ntype 
! beginning of the loop over the source types. 
      if (totlu(stype) /= 0.) then 
! check if there are any flux in that source type otherwise skip this lamp 
        if (verbose >= 1) print*,' turning on lamps',stype 
        if (verbose >= 1) write(2,*) ' turning on lamps', &
                  stype 
        do x_s = imin(stype),imax(stype) 
! beginning of the loop over the column (longitude the) of the domain. 
          do y_s = jmin(stype),jmax(stype) 
! beginning of the loop over the rows (latitud) of the domain. 
            intdir = 0. 
            itotind = 0. 
            itodif = 0. 
            itotrd = 0. 
            isourc = 0. 
            rx_s = real(x_s)*dx 
            ry_s = real(y_s)*dy 
            if (lamplu(x_s,y_s,stype)  /=  0.) then 
! if the luminosite of the case is null, the program ignore this case. 
              z_s = (altsol(x_s,y_s)+lampal(x_s,y_s)) 
! Definition of the position (metre) vertical of the source. 
 
! ********************************************************************************************************* 
! calculation of the direct radiance of sources falling on a surface perpendicular 
! to the viewing angle Units of W/nm/m2/sr 
! ********************************************************************************************************* 
              rx = rx_obs+20000.*ix 
              ry = ry_obs+20000.*iy 
              rz = z_obs+20000.*iz 
              dho = sqrt((rx_obs-rx_s)**2. &
                        +(ry_obs-ry_s)**2.) 
              if ((dho > 0.).and.(z_s /= z_obs)) then 
                call anglezenithal(rx_obs,ry_obs,z_obs &
                          ,rx_s,ry_s,z_s,dzen) 
! zenithal angle source-observer 
                call angleazimutal(rx_obs,ry_obs,rx_s, &
                          ry_s,angazi) 
! computation of the angle azimutal direct line of sight-source 
                if (dzen > pi/4.) then 
! 45deg. it is unlikely to have a 1km high mountain less than 1 
                  call horizon(x_obs,y_obs,z_obs,dx,dy, &
                            altsol,angazi,zhoriz,dh) 
                  if (dh <= dho) then 
                    if (dzen-zhoriz < 0.00001) then 
! shadow the path line of sight-source is not below the horizon => we compute 
                      hh = 1. 
                    else 
                      hh = 0. 
                    end if 
                  else 
                    hh = 1. 
                  end if 
                else 
                  hh = 1. 
                end if 
                ff = 0. 
                if (obsobs == 1) then 
! sub-grid obstacles 
                  if (dho > drefle(x_obs,y_obs)+drefle(x_s,y_s)) then 
! light path to observer larger than the mean free path -> subgrid obstacles 
                    angmin = pi/2.-atan2((altsol(x_obs,y_obs)+ &
                              obsH(x_obs,y_obs)-z_obs),drefle(x_obs, &
                              y_obs)) 
                    if (dzen < angmin) then 
! condition sub-grid obstacles direct. 
                      ff = 0. 
                    else 
                      ff = ofill(x_obs,y_obs) 
                    end if 
                  end if 
! end light path to the observer larger than mean free path 
                end if 
 
 
 
 
 
                call anglezenithal(rx_s,ry_s,z_s &
                          ,rx_obs,ry_obs,z_obs,dzen) 
! zenithal angle source-observer 
                ff2 = 0. 
                if (dho > drefle(x_s,y_s)) then 
! light path from source larger than the mean free path -> subgrid obstacles 
                  angmin = pi/2.-atan2((altsol(x_s,y_s)+ &
                            obsH(x_s,y_s)-z_s),drefle(x_s, &
                            y_s)) 
                  if (dzen < angmin) then 
! condition sub-grid obstacles direct. 
                    ff2 = 0. 
                  else 
                    ff2 = ofill(x_s,y_s) 
                  end if 
                end if                                                   ! end light path to the observer larger than mean free path 
                call anglezenithal(rx_obs,ry_obs,z_obs &
                          ,rx_s,ry_s,z_s,dzen) 
! zenithal angle source-observer 
 
 
 
 
! projection angle of line to the lamp and the viewing angle 
                call angle3points (rx_s,ry_s,z_s,rx_obs, &
                          ry_obs,z_obs,rx,ry,rz,dang) 
! scattering angle. 
                dang = pi-dang 
! computation of the solid angle of the line of sight voxel seen from the source 
                anglez = nint(180.*(pi-dzen)/pi)+1 
                P_dir = pvalno(anglez,stype) 
! computation of the flux direct reaching the line of sight voxel 
                if ((cos(dang) > 0.).and.(dang < pi/2.)) &
                          then 
                ddir_obs = sqrt((rx_obs-rx_s)**2.+ &
                          (ry_obs-ry_s)**2.+(z_obs-z_s)**2.) 
! distance direct sight between source and observer 
! computation of the solid angle 1m^2 at the observer as seen from the source 
                omega = 1.*abs(cos(dang))/ddir_obs**2. 
                call transmitm(dzen,z_obs,z_s,ddir_obs, &
                          transm,tranam,tabs) 
                call transmita(dzen,z_obs,z_s,ddir_obs, &
                          haer,transa,tranaa) 
                call transmitl(dzen,z_obs,z_s,ddir_obs, &
                          hlay,transl,tranal) 
                if (dang < dfov) then 
!check if the reflecting surface enter the field of view of the observer 
                  direct = direct+lamplu(x_s,y_s,stype)* &
                            transa*transm*transl*P_dir*omega*(1.-ff)*(1.-ff2) &
                            *hh/(pi*dfov**2.)                                      ! correction for obstacle filling factor 
                end if 
                irdirect = irdirect+lamplu(x_s,y_s,stype)* &
                          transa*transm*transl*P_dir*omega*(1.-ff)*(1.-ff2)*hh    ! correction for obstacle filling factor 
              end if 
            end if 
 
! ********************************************************************************** 
! * computation of the direct light toward the observer by the ground reflection   * 
! ********************************************************************************** 
 
            xsrmi = x_s-boxx 
            if (xsrmi < 1) xsrmi = 1 
            xsrma = x_s+boxx 
            if (xsrma > nbx) xsrma = nbx 
            ysrmi = y_s-boxy 
            if (ysrmi < 1) ysrmi = 1 
            ysrma = y_s+boxy 
            if (ysrma > nby) ysrma = nby 
            do x_sr = xsrmi,xsrma                                       ! beginning of the loop over the column (longitude) reflecting. 
              rx_sr = real(x_sr)*dx 
              do y_sr = ysrmi,ysrma                                     ! beginning of the loop over the rows (latitu) reflecting. 
                ry_sr = real(y_sr)*dy 
                irefl = 0. 
                z_sr = altsol(x_sr,y_sr) 
                if((x_sr > nbx).or.(x_sr < 1).or. &
                          (y_sr > nby).or.(y_sr < 1)) then 
                if (verbose == 2) then 
                  print*,'Ground cell out of borders' 
                end if 
              else 
                if((x_s == x_sr).and.(y_s == y_sr) &
                          .and.(z_s == z_sr)) then 
                if (verbose == 2) then 
                  print*,'Source pos = Ground cell' 
                end if 
              else 
! if haut is negative, the ground cell is lighted from below 
                haut = -(rx_s-rx_sr)*tan( &
                          inclix(x_sr,y_sr))-(ry_s- &
                          ry_sr)*tan(incliy(x_sr, &
                          y_sr))+z_s-z_sr 
                if (haut  >  0.) then                          ! condition: the ground cell is lighted from above 
! computation of the zenithal angle between the source and the surface reflectance 
! computation of the zenithal angle between the source and the line of sight voxel. 
! end of the case "observer at the same latitu/longitude than the source". 
 
                call anglezenithal(rx_s,ry_s, &
                          z_s,rx_sr,ry_sr,z_sr, &
                          angzen) 
! computation of the transmittance between the source and the ground surface 
                distd = sqrt((rx_s-rx_sr)**2. &
                          +(ry_s-ry_sr)**2.+ &
                          (z_s-z_sr)**2.) 
                call transmitm(angzen,z_s, &
                          z_sr,distd,transm,tranam,tabs) 
                call transmita(angzen,z_s, &
                          z_sr,distd,haer,transa,tranaa) 
                call transmitl(angzen,z_s,z_sr,distd, &
                          hlay,transl,tranal) 
! computation of the solid angle of the reflecting cell seen from the source 
                xc = dble(x_sr)*dble(dx)                        ! Position in meters of the observer voxel (longitude). 
                yc = dble(y_sr)*dble(dy)                        ! Position in meters of the observer voxel (latitu). 
                zc = dble(z_sr)                                 ! Position in meters of the observer voxel (altitude). 
                xn = dble(x_s)*dble(dx)                         ! Position in meters of the source (longitude). 
                yn = dble(y_s)*dble(dy)                         ! Position in meters of the source (latitu). 
                zn = dble(z_s)                                  ! Position in meters of the source (altitude). 
                epsilx = inclix(x_sr,y_sr)                      ! tilt along x of the ground reflectance 
                epsily = incliy(x_sr,y_sr)                      ! tilt along x of the ground reflectance 
                if (dx > reflsiz) then                       ! use a sub-grid surface when the reflectance radius is smaller than the cell size 
                if ((x_sr == x_s).and.(y_sr &
                          == y_s)) then 
                dxp = reflsiz 
              else 
                dxp = dx 
              end if 
            else 
              dxp = dx 
            end if 
            if (dy > reflsiz) then 
              if ((x_sr == x_s).and.(y_sr &
                        == y_s)) then 
              dyp = reflsiz 
            else 
              dyp = dy 
            end if 
          else 
            dyp = dy 
          end if 
          r1x = xc-dble(dxp)/2.-xn                        ! computation of the composante along x of the first vector. 
          r1y = yc+dble(dyp)/2.-yn                        ! computation of the composante along y of the first vector. 
! computation of the composante en z of the first vector. 
          r1z = zc-tan(dble(epsilx))* &
                    dble(dxp)/2.+tan(dble(epsily)) &
                    *dble(dyp)/2.-zn 
          r2x = xc+dble(dxp)/2.-xn                        ! computation of the composante along x of the second vector. 
          r2y = yc+dble(dyp)/2.-yn                        ! computation of the composante along y of the second vector. 
! computation of the composante en z of the second vector. 
          r2z = zc+tan(dble(epsilx))* &
                    dble(dxp)/2.+tan(dble(epsily)) &
                    *dble(dyp)/2.-zn 
          r3x = xc-dble(dxp)/2.-xn                        ! computation of the composante along x of the third vector. 
          r3y = yc-dble(dyp)/2.-yn                        ! computation of the composante along y of the third vector. 
! computation of the composante en z of the third vector. 
          r3z = zc-tan(dble(epsilx))* &
                    dble(dxp)/2.-tan(dble(epsily)) &
                    *dble(dyp)/2.-zn 
          r4x = xc+dble(dxp)/2.-xn                        ! computation of the composante along x of the fourth vector. 
          r4y = yc-dble(dyp)/2.-yn                        ! computation of the composante along y of the fourth vector. 
! computation of the composante en z of the fourth vector. 
          r4z = zc+tan(dble(epsilx))* &
                    dble(dxp)/2.-tan(dble(epsily)) &
                    *dble(dyp)/2.-zn 
! Call of the routine anglesolide to compute the angle solide. 
          call anglesolide(omega,r1x, &
                    r1y,r1z,r2x,r2y,r2z,r3x,r3y, &
                    r3z,r4x,r4y,r4z) 
          if (omega < 0.) then 
            print*,'ERROR: Solid angle of the reflecting surface < 0.' 
            stop 
          end if 
! estimation of the half of the underlying angle of the solid angle       ! this angle servira a obtenir un meilleur isime (moyenne) of 
!                                                                         ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang. 
          ouvang = sqrt(omega/pi)                         ! Angle in radian. 
          ouvang = ouvang*180./pi                         ! Angle in degrees. 
! computation of the photometric function of the light fixture toward the reflection surface 
!======================================================================= 
 
          anglez = nint(180.*angzen/pi) 
          if (anglez < 0) &
                    anglez = -anglez 
          if (anglez > 180) anglez = 360 &
                    -anglez 
          anglez = anglez+1                               ! Transform the angle in integer degree into the position in the array. 
! average +- ouvang 
          naz = 0 
          nbang = 0. 
          P_indir = 0. 
          do na = -nint(ouvang),nint(ouvang) 
            naz = anglez+na 
            if (naz < 0) naz = -naz 
            if (naz > 181) naz = 362-naz                 ! symetric function 
            if (naz == 0) naz = 1 
            P_indir = P_indir+pvalno(naz, &
                      stype)*abs(sin(pi*real(naz) &
                      /180.))/2. 
            nbang = nbang+1.*abs(sin(pi* &
                      real(naz)/180.))/2. 
          end do 
          P_indir = P_indir/nbang 
! computation of the flux reaching the reflecting surface 
          flrefl = lamplu(x_s,y_s,stype)* &
                    P_indir*omega*transm*transa*transl 
! computation of the reflected intensity leaving the ground surface 
          irefl1 = flrefl*srefl/pi                        ! The factor 1/pi comes from the normalisation of the fonction 
 
! ********************************************************************************************************* 
! calculation of the direct radiance from reflection falling on a surface perpendicular 
! to the viewing angle Units of W/nm/m2/sr 
! ********************************************************************************************************* 
          dho = sqrt((rx_obs-rx_sr)**2. &
                    +(ry_obs-ry_sr)**2.) 
          if ((dho > 0.).and.(z_s /= z_obs)) then 
! zenithal angle source-observer 
            call anglezenithal(rx_obs,ry_obs,z_obs &
                      ,rx_sr,ry_sr,z_sr,dzen) 
! computation of the angle azimutal direct line of sight-source 
            call angleazimutal(rx_obs,ry_obs,rx_sr, &
                      ry_sr,angazi) 
            if (dzen > pi/4.) then                     ! 45deg. it is unlikely to have a 1km high mountain less than 1 
            call horizon(x_obs,y_obs,z_obs,dx,dy, &
                      altsol,angazi,zhoriz,dh) 
            if (dh <= dho) then 
              if (dzen-zhoriz < 0.00001) then        ! shadow the path line of sight-source is not below the horizon => we compute 
              hh = 1. 
            else 
              hh = 0. 
            end if 
          else 
            hh = 1. 
          end if 
        else 
          hh = 1. 
        end if 
! sub-grid obstacles 
        ff = 0. 
        if (obsobs == 1) then 
! light path to observer larger than the mean free path -> subgrid obstacles 
          if (dho > drefle(x_obs,y_obs)+ &
                    drefle(x_sr,y_sr)) then 
          angmin = pi/2.-atan2((altsol(x_obs,y_obs) &
                    +obsH(x_obs,y_obs)-z_obs),drefle(x_obs, &
                    y_obs)) 
          if (dzen < angmin) then                  ! condition sub-grid obstacles direct. 
          ff = 0. 
        else 
          ff = ofill(x_obs,y_obs) 
        end if 
      end if                                       ! end light path to the observer larger than mean free path 
    end if 
 
 
! zenithal angle surface-observer 
    call anglezenithal(rx_sr,ry_sr,z_sr &
              ,rx_obs,ry_obs,z_obs,dzen) 
    ff2 = 0. 
    if (dho > drefle(x_sr,y_sr)) then                      ! light path from reflecting surface larger than the mean free path -> subgrid obstacles 
    angmin = pi/2.-atan2((altsol(x_sr,y_sr)+ &
              obsH(x_sr,y_sr)-z_sr),drefle(x_sr, &
              y_sr)) 
    if (dzen < angmin) then                              ! condition sub-grid obstacles direct. 
    ff2 = 0. 
    else 
      ff2 = ofill(x_sr,y_sr) 
    end if 
    end if 
! end light path to the observer larger than mean free path 
! zenithal angle source-observer 
    call anglezenithal(rx_obs,ry_obs,z_obs &
              ,rx_sr,ry_sr,z_sr,dzen) 
 
 
 
! projection angle of line to the lamp and the viewing angle 
! scattering angle. 
    call angle3points (rx_sr,ry_sr,z_sr, &
              rx_obs,ry_obs,z_obs,rx,ry,rz,dang) 
    dang = pi-dang 
 
! computation of the flux direct reaching the line of sight voxel 
    if ((cos(dang) > 0.).and.(dang < pi/2.)) &
              then 
! distance direct sight between source and observer 
    ddir_obs = sqrt((rx_obs-rx_sr)**2.+ &
              (ry_obs-ry_sr)**2.+(z_obs-z_sr)**2.) 
! computation of the solid angle of the line of sight voxel seen from the source 
    omega = 1.*abs(cos(dang))/ddir_obs**2. 
    call transmitm(dzen,z_obs,z_sr,ddir_obs, &
              transm,tranam,tabs) 
    call transmita(dzen,z_obs,z_sr,ddir_obs, &
              haer,transa,tranaa) 
    call transmitl(dzen,z_obs,z_sr,ddir_obs, &
              hlay,transl,tranal) 
    if (dang < dfov) then                    ! check if the reflecting surface enter the field of view of the observer 
    rdirect = rdirect+irefl1*omega*transa* &
              transm*transl*hh*(1.-ff)*(1.-ff2) &
              /(pi*dfov**2.) 
    end if 
    irrdirect = irrdirect+irefl1*omega*transa* &
              transm*transl*hh*(1.-ff)*(1.-ff2) 
    end if 
 
    end if 
    end if 
    end if 
    end if 
    end do 
    end do 
 
    end if 
    end do 
    end do 
    end if 
    end do 
 
! End of calculation of the direct radiances 
! ================================= 
 
 
 
    if (fsswit /= 0) then 
! ================================= 
! Calculation of the scattered radiances 
 
 
! temporaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      cloudtop = 100000. 
 
      if ((z_obs >= cloudbase).and.(z_obs <= cloudtop)) then 
        print*,'The observer is inside the cloud! Abort computing.', &
                  z_obs,cloudbase 
        stop 
      end if 
 1110  format(I4,1x,I4,1x,I4) 
      fctcld = 0. 
      ftocap = 0.                                                         ! Initialisation of the value of flux received by the sensor 
! calculating the distance before the line of sight beeing blocked by topography 
      call horizon(x_obs,y_obs,z_obs,dx,dy,altsol,angaz1,zhoriz, &
                dhmax) 
      rx_c = real(x_obs)*dx-ix*scal/2. 
      ry_c = real(y_obs)*dx-iy*scal/2. 
      z_c = z_obs-iz*scal/2. 
      do icible = 1,ncible                                                ! beginning of the loop over the line of sight voxels 
        rx_c = rx_c+ix*(scalo/2.+scal/2.) 
        ry_c = ry_c+iy*(scalo/2.+scal/2.) 
        dh0 = sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2) 
        if ((dh0 <= dhmax).or.((dh0 > dhmax).and.(angze1-zhoriz < &
                  0.00001))) then 
! the line of sight is not yet blocked by the topography 
        x_c = nint(rx_c/dx) 
        if (x_c < 1) x_c = 1 
        if (x_c > width) x_c = width 
        y_c = nint(ry_c/dy) 
        if (y_c < 1) y_c = 1 
        if (y_c > width) y_c = width 
        z_c = z_c+iz*(scalo/2.+scal/2.) 
        if (z_c > altsol(x_c,y_c)) then 
! stop the calculation of the viewing line when the increment is lower than 1/stoplim 
          if ((fcapt >= ftocap/stoplim).and.(z_c < cloudbase).and. &
                    (z_c < 35000.)) then                                           ! or when hitting a cloud or when z>40km (scattering probability =0 (given precision) 
          fcapt = 0. 
          do i = 1,nbx 
            do j = 1,nby 
              FCA(i,j) = 0. 
            end do 
          end do 
! Calculate the solid angle of the line of sight voxel unit voxel 
! (1 m^3) given the fixed FOV of the observer. 
! For line of sight voxel near the observer 
! we need to calculate the scattering on a part of the voxel. For far 
! voxels we may be needed to increase the solid angle since the FOV can 
! encompass more than the voxel size. This correction is done with the 
! portio parameter calculated as the ratio of the solid angle of the 
! observer FOV over the line of sight voxel solid angle as seen from the 
! observer. 
          distd = sqrt((rx_c-rx_obs)**2.+ &
                    (ry_c-ry_obs)**2.+(z_c-z_obs)**2.) 
! computation of the Solid angle of the line of sight voxel seen from the observer 
          omega = 1./distd**2. 
          if (omega > omemax) then 
            omega = 0. 
            portio = 0. 
          else 
            portio = (omefov/omega) 
          end if 
          itotci = 0.                                                     ! Initialisation of the contribution of the line of sight at the sensor level 
          do i = 1,nbx 
            do j = 1,nby 
              ITC(i,j) = 0. 
            end do 
          end do 
! Condition line of sight inside the modelling domain 
          if( (rx_c > real(nbx*dx)).or.(rx_c < dx).or. &
                    (ry_c > (nby*dy)).or.(ry_c < dy)) then 
        else 
          if (verbose >= 1) print*,'================================ &
                    ================' 
          if (verbose >= 1) print*,' progression along the line of sight :' &
                    ,icible 
          if (verbose >= 1) print*,' horizontal dist. line of sight =', &
                    sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.),' m' 
          if (verbose >= 1) print*,' vertical dist. line of sight =', &
                    abs(z_c-z_obs),' m' 
          if (verbose >= 1) write(2,*) '======================== &
                    =====================' 
          if (verbose >= 1) write(2,*) ' progression along the line of sight &
                    :',icible 
          if (verbose >= 1) write(2,*) ' horizontal dist. line of sight =', &
                    sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.),' m' 
          if (verbose >= 1) write(2,*) ' vertical dist. line of sight =', &
                    abs(z_c-z_obs),' m' 
          dis_obs = sqrt((z_c-z_obs)**2.+(ry_c-ry_obs)**2. &
                    +(rx_c-rx_obs)**2.) 
          if (dis_obs == 0.) then 
            print*,'ERROR problem with dis_obs',dis_obs 
            print*,rx_c,x_obs,y_c,y_obs,z_c,z_obs 
            stop 
          end if 
          ometif = pi*(diamobj/2.)**2./dis_obs**2. 
! beginning of the loop over the types of light sources 
          do stype = 1,ntype                                            ! beginning of the loop over the source types. 
            if (totlu(stype) /= 0.) then                              ! check if there are any flux in that source type otherwise skip this lamp 
            if (verbose >= 1) print*,' turning on lamps',stype 
            if (verbose >= 1) write(2,*) ' turning on lamps', &
                      stype 
            itotty = 0.                                               ! Initialisation of the contribution of a source types to 
            do x_s = 1,nbx                                            ! the intensity toward the sensor by a line of sight voxel. 
              do y_s = 1,nby 
                ITT(x_s,y_s,stype) = 0. 
              end do 
            end do 
            do x_s = imin(stype),imax(stype)                          ! beginning of the loop over the column (longitude the) of the domain. 
              do y_s = jmin(stype),jmax(stype)                        ! beginning of the loop over the rows (latitud) of the domain. 
                intdir = 0. 
                itotind = 0. 
                itodif = 0. 
                itotrd = 0. 
                isourc = 0. 
                rx_s = real(x_s)*dx 
                ry_s = real(y_s)*dy 
                if (lamplu(x_s,y_s,stype)  /=  0.) then             ! if the luminosite of the case is null, the program ignore this case. 
                z_s = (altsol(x_s,y_s)+lampal(x_s,y_s))             ! Definition of the position (metre) vertical of the source. 
 
! ********************************************************************************************************* 
! * computation of the scattered intensity toward the observer by a line of sight voxel from the source   * 
! ********************************************************************************************************* 
 
                dirck = 0                                           ! Initialisation of the verification of the position of the source 
! if the position of the source and the line of sight voxel are the 
                if ((rx_s == rx_c).and.(ry_s == ry_c).and. &
                          (z_s == z_c)) &
                          then 
                dirck = 1 
                if (verbose >= 1) then 
                  print*,'Source = line of sight' 
                end if 
              end if                                             ! end of the case positions x and y source and line of sight voxel identical. 
              if (dirck /= 1) then                              ! the source is not at the line of sight voxel position 
! computation of the zenithal angle between the source and the line of sight 
! computation of the horizon for the resolved shadows direct              ! horizon resolution is 1 degree 
              distd = sqrt((rx_c-rx_s)**2. &
                        +(ry_c-ry_s)**2. &
                        +(z_c-z_s)**2.) 
              dho = sqrt((rx_c-rx_s)**2. &
                        +(ry_c-ry_s)**2.) 
              call anglezenithal(rx_s,ry_s,z_s &
                        ,rx_c,ry_c,z_c,angzen)                      ! computation of the zenithal angle between the source and the line of sight voxel. 
! computation of the angle azimutal direct line of sight-source 
              call angleazimutal(rx_s,ry_s,rx_c, &
                        ry_c,angazi) 
              if (angzen > pi/4.) then                   ! 45deg. it is unlikely to have a 1km high mountain less than 1 
              call horizon(x_s,y_s,z_s,dx,dy,altsol, &
                        angazi,zhoriz,dh) 
              if (dh <= dho) then 
                if (angzen-zhoriz < 0.00001) then      ! shadow the path line of sight-source is not below the horizon => we compute 
                hh = 1. 
              else 
                hh = 0. 
              end if 
            else 
              hh = 1. 
            end if 
          else 
            hh = 1. 
          end if 
! sub-grid obstacles 
          ff = 0. 
          if (dho > drefle(x_s,y_s)) then            ! light path to observer larger than the mean free path -> subgrid obstacles 
          angmin = pi/2.-atan2((altsol(x_s,y_s)+ &
                    obsH(x_s,y_s)-z_s),drefle(x_s,y_s)) 
          if (angzen < angmin) then                ! condition sub-grid obstacles direct. 
          ff = 0. 
        else 
          ff = ofill(x_s,y_s) 
        end if 
      end if 
! computation of the transmittance between the source and the line of sight 
      call transmitm(angzen,z_s,z_c,distd, &
                transm,tranam,tabs) 
      call transmita(angzen,z_s,z_c,distd, &
                haer,transa,tranaa) 
      call transmitl(angzen,z_s,z_c,distd, &
                hlay,transl,tranal) 
! computation of the solid angle of the line of sight voxel seen from the source 
      omega = 1./distd**2. 
      if (omega > omemax) omega = 0. 
      anglez = nint(180.*angzen/pi)+1 
      P_dir = pvalno(anglez,stype) 
! computation of the flux reaching the line of sight voxel 
      fldir = lamplu(x_s,y_s,stype)*P_dir* &
                omega*transm*transa*transl*(1.-ff)*hh       ! correction for obstacle filling factor 
! computation of the scattering probability of the direct light 
! distance pour traverser la cellule unitaire parfaitement orientée 
      if (omega /= 0.) then 
! scattering angle. 
        call angle3points (rx_s,ry_s,z_s,rx_c, &
                  ry_c,z_c,rx_obs,ry_obs,z_obs, &
                  angdif) 
! scattering probability of the direct light. ############################################ secdif et un etaient inverses 
        call diffusion(angdif, &
                  tranam,tranaa,tranal,un,secdif,secdil, &
                  fdifan,fdifl,haer,hlay,pdifdi,z_c) 
      else 
        pdifdi = 0. 
      end if 
! computation of the source contribution to the scattered intensity toward the sensor by a line of sight voxel 
      intdir = fldir*pdifdi 
! contribution of the cloud reflection of the light coming directly from the source 
      if (cloudt /= 0) then                       ! line of sight voxel = cloud 
      if (cloudbase-z_c <= iz*scal) then 
        call anglezenithal(rx_c,ry_c,z_c, &
                  rx_obs,ry_obs,z_obs,azcl1)              ! zenith angle from cloud to observer 
        call anglezenithal(rx_c,ry_c,z_c, &
                  rx_s,ry_s,z_s,azcl2)                    ! zenith angle from source to cloud 
        doc2 = (rx_c-rx_obs)**2.+ &
                  (ry_c-ry_obs)**2.+(z_c-z_obs)**2. 
        dsc2 = (rx_s-rx_c)**2.+ &
                  (ry_s-ry_c)**2.+(z_s-z_c)**2. 
! cloud intensity from direct illum 
        call cloudreflectance(angzen, &
                  cloudt,rcloud) 
        icloud = icloud+ &
                  fldir/omega*rcloud*doc2*omefov* &
                  abs(cos(azcl2)/cos(azcl1))/dsc2/pi 
      end if 
    end if 
    else 
      intdir = 0. 
    end if                                             ! end of the case position source is not equal to the line of sight voxel position 
! end of the computation of the scattered intensity 
 
 
 
 
! ********************************************************************************************************************** 
! * computation of the scattered light toward the observer by a line of sight voxel lighted by the ground reflection    * 
! ********************************************************************************************************************** 
! etablissement of the conditions ands boucles 
    itotind = 0.                                    ! Initialisation of the reflected intensity of the source 
    itotrd = 0. 
    xsrmi = x_s-boxx 
    if (xsrmi < 1) xsrmi = 1 
    xsrma = x_s+boxx 
    if (xsrma > nbx) xsrma = nbx 
    ysrmi = y_s-boxy 
    if (ysrmi < 1) ysrmi = 1 
    ysrma = y_s+boxy 
    if (ysrma > nby) ysrma = nby 
    do x_sr = xsrmi,xsrma                           ! beginning of the loop over the column (longitude) reflecting. 
      rx_sr = real(x_sr)*dx 
      do y_sr = ysrmi,ysrma                         ! beginning of the loop over the rows (latitu) reflecting. 
        ry_sr = real(y_sr)*dy 
        irefl = 0. 
        z_sr = altsol(x_sr,y_sr) 
        if((x_sr > nbx).or.(x_sr < 1).or. &
                  (y_sr > nby).or.(y_sr < 1)) then 
        if (verbose == 2) then 
          print*,'Ground cell out of borders' 
        end if 
      else 
        if((x_s == x_sr).and.(y_s == y_sr) &
                  .and.(z_s == z_sr)) then 
        if (verbose == 2) then 
          print*,'Source pos = Ground cell' 
        end if 
      else 
! if haut is negative, the ground cell is lighted from below 
        haut = -(rx_s-rx_sr)*tan( &
                  inclix(x_sr,y_sr))-(ry_s- &
                  ry_sr)*tan(incliy(x_sr, &
                  y_sr))+z_s-z_sr 
        if (haut  >  0.) then              ! condition: the ground cell is lighted from above 
! computation of the zenithal angle between the source and the surface reflectance 
! computation of the zenithal angle between the source and the line of sight voxel. 
! end of the case "observer at the same latitu/longitude than the source". 
        call anglezenithal(rx_s,ry_s, &
                  z_s,rx_sr,ry_sr,z_sr, &
                  angzen) 
! computation of the transmittance between the source and the ground surface 
        distd = sqrt((rx_s-rx_sr)**2. &
                  +(ry_s-ry_sr)**2.+ &
                  (z_s-z_sr)**2.) 
        call transmitm(angzen,z_s, &
                  z_sr,distd,transm,tranam,tabs) 
        call transmita(angzen,z_s, &
                  z_sr,distd,haer,transa,tranaa) 
        call transmitl(angzen,z_s,z_sr, &
                  distd,hlay,transl,tranal) 
! computation of the solid angle of the reflecting cell seen from the source 
        xc = dble(x_sr)*dble(dx)            ! Position in meters of the observer voxel (longitude). 
        yc = dble(y_sr)*dble(dy)            ! Position in meters of the observer voxel (latitu). 
        zc = dble(z_sr)                     ! Position in meters of the observer voxel (altitude). 
        xn = dble(x_s)*dble(dx)             ! Position in meters of the source (longitude). 
        yn = dble(y_s)*dble(dy)             ! Position in meters of the source (latitu). 
        zn = dble(z_s)                      ! Position in meters of the source (altitude). 
        epsilx = inclix(x_sr,y_sr)          ! tilt along x of the ground reflectance 
        epsily = incliy(x_sr,y_sr)          ! tilt along x of the ground reflectance 
        if (dx > reflsiz) then           ! use a sub-grid surface when the reflectance radius is smaller than the cell size 
        if ((x_sr == x_s).and.(y_sr &
                  == y_s)) then 
        dxp = reflsiz 
      else 
        dxp = dx 
      end if 
    else 
      dxp = dx 
    end if 
    if (dy > reflsiz) then 
      if ((x_sr == x_s).and.(y_sr &
                == y_s)) then 
      dyp = reflsiz 
    else 
      dyp = dy 
    end if 
    else 
      dyp = dy 
    end if 
    r1x = xc-dble(dxp)/2.-xn            ! computation of the composante along x of the first vector. 
    r1y = yc+dble(dyp)/2.-yn            ! computation of the composante along y of the first vector. 
! computation of the composante en z of the first vector. 
    r1z = zc-tan(dble(epsilx))* &
              dble(dxp)/2.+tan(dble(epsily)) &
              *dble(dyp)/2.-zn 
    r2x = xc+dble(dxp)/2.-xn            ! computation of the composante along x of the second vector. 
    r2y = yc+dble(dyp)/2.-yn            ! computation of the composante along y of the second vector. 
! computation of the composante en z of the second vector. 
    r2z = zc+tan(dble(epsilx))* &
              dble(dxp)/2.+tan(dble(epsily)) &
              *dble(dyp)/2.-zn 
    r3x = xc-dble(dxp)/2.-xn            ! computation of the composante along x of the third vector. 
    r3y = yc-dble(dyp)/2.-yn            ! computation of the composante along y of the third vector. 
! computation of the composante en z of the third vector. 
    r3z = zc-tan(dble(epsilx))* &
              dble(dxp)/2.-tan(dble(epsily)) &
              *dble(dyp)/2.-zn 
    r4x = xc+dble(dxp)/2.-xn            ! computation of the composante along x of the fourth vector. 
    r4y = yc-dble(dyp)/2.-yn            ! computation of the composante along y of the fourth vector. 
! computation of the composante en z of the fourth vector. 
    r4z = zc+tan(dble(epsilx))* &
              dble(dxp)/2.-tan(dble(epsily)) &
              *dble(dyp)/2.-zn 
! Call of the routine anglesolide to compute the angle solide. 
    call anglesolide(omega,r1x, &
              r1y,r1z,r2x,r2y,r2z,r3x,r3y, &
              r3z,r4x,r4y,r4z) 
    if (omega < 0.) then 
      print*,'ERROR: Solid angle of the reflecting surface < 0.' 
      stop 
    end if 
! estimation of the half of the underlying angle of the solid angle       ! this angle servira a obtenir un meilleur isime (moyenne) of 
!                                                                         ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang. 
    ouvang = sqrt(omega/pi)             ! Angle in radian. 
    ouvang = ouvang*180./pi             ! Angle in degrees. 
! computation of the photometric function of the light fixture toward the reflection surface 
!======================================================================= 
 
    anglez = nint(180.*angzen/pi) 
    if (anglez < 0) &
              anglez = -anglez 
    if (anglez > 180) anglez = 360 &
              -anglez 
    anglez = anglez+1                   ! Transform the angle in integer degree into the position in the array. 
! average +- ouvang 
    naz = 0 
    nbang = 0. 
    P_indir = 0. 
    do na = -nint(ouvang),nint(ouvang) 
      naz = anglez+na 
      if (naz < 0) naz = -naz 
      if (naz > 181) naz = 362-naz     ! symetric function 
      if (naz == 0) naz = 1 
      P_indir = P_indir+pvalno(naz, &
                stype)*abs(sin(pi*real(naz) &
                /180.))/2. 
      nbang = nbang+1.*abs(sin(pi* &
                real(naz)/180.))/2. 
    end do 
    P_indir = P_indir/nbang 
! computation of the flux reaching the reflecting surface 
    flrefl = lamplu(x_s,y_s,stype)* &
              P_indir*omega*transm*transa* &
              transl 
! computation of the reflected intensity leaving the ground surface 
    irefl1 = flrefl*srefl/pi            ! The factor 1/pi comes from the normalisation of the fonction 
 
 
 
 
 
! ************************************************************************************** 
! * computation of the 2nd scattering contributions (2 order scattering and after reflection) 
! ************************************************************************************** 
    if (effdif > 0.) then 
      nss = 0 
      ndi = 0 
      do idi = 1,ndiff                                                      ! beginning of the loop over the scattering voxels. 
        rx_dif = zondif(idi,1)+(rx_s+rx_c)/2. 
        x_dif = nint(rx_dif/dx) 
        ry_dif = zondif(idi,2)+(ry_s+ry_c)/2. 
        y_dif = nint(ry_dif/dy) 
        z_dif = zondif(idi,3)+(z_s+z_c)/2. 
        id = nint(rx_dif/dx) 
        if (id > width) id = width 
        if (id < 1) id = 1 
        jd = nint(ry_dif/dy) 
        if (jd > width) jd = width 
        if (jd < 1) jd = 1 
        if (z_dif-siz/2. <= altsol(id,jd).or.(z_dif > 35000.).or. &
                  (z_dif > cloudbase)) then                                        ! beginning diffusing cell underground 
        ndi = ndi+1 
      else 
        ds1 = sqrt((rx_sr-rx_dif)**2.+(ry_sr-ry_dif)**2.+ &
                  (z_sr-z_dif)**2.) 
        ds2 = sqrt((rx_c-rx_dif)**2.+(ry_c-ry_dif)**2.+ &
                  (z_c-z_dif)**2.) 
        ds3 = sqrt((rx_s-rx_dif)**2.+(ry_s-ry_dif)**2.+ &
                  (z_s-z_dif)**2.) 
        if ((ds1 < dss).or.(ds2 < dss).or.(ds3 < dss)) then 
          nss = nss+1 
!       print*,ds1,ds2,ds3,'c',rx_c,ry_c,z_c,'d',rx_dif,ry_dif,z_dif 
        else 
          dho = sqrt((rx_dif-rx_sr)**2.+(ry_dif-ry_sr)**2.) 
! computation of the zenithal angle between the reflection surface and the scattering voxel 
! shadow reflection surface-scattering voxel 
          call anglezenithal(rx_sr,ry_sr, &
                    z_sr,rx_dif,ry_dif,z_dif,angzen)                              ! computation of the zenithal angle reflection surface - scattering voxel. 
          call angleazimutal(rx_sr,ry_sr,rx_dif,ry_dif,angazi)          ! computation of the angle azimutal line of sight-scattering voxel 
! horizon blocking not a matte because dif are closeby and some downward 
          hh = 1. 
! sub-grid obstacles 
          ff = 0. 
          if (dho > drefle(x_sr,y_sr)) then                            ! light path to observer larger than the mean free path -> subgrid obstacles 
          angmin = pi/2.-atan2(obsH(x_sr,y_sr),drefle(x_sr,y_sr)) 
          if (angzen < angmin) then                                  ! condition obstacle reflechi->scattered 
          ff = 0. 
        else 
          ff = ofill(x_sr,y_sr) 
        end if 
      end if 
! computation of the transmittance between the reflection surface and the scattering voxel 
      distd = sqrt((rx_dif-rx_sr)**2.+(ry_dif-ry_sr)**2.+ &
                (z_dif-z_sr)**2.) 
      call transmitm(angzen,z_sr,z_dif,distd,transm,tranam,tabs) 
      call transmita(angzen,z_sr,z_dif,distd,haer,transa,tranaa) 
      call transmitl(angzen,z_sr,z_dif, &
                distd,hlay,transl,tranal) 
! computation of the solid angle of the scattering voxel seen from the reflecting surface 
      omega = 1./distd**2. 
      if (omega > omemax) omega = 0. 
! computing flux reaching the scattering voxel 
      fldif2 = irefl1*omega*transm*transa*transl*(1.-ff)*hh 
! computing the scattering probability toward the line of sight voxel 
! cell unitaire 
      if (omega /= 0.) then 
! scattering angle. 
        call angle3points (rx_sr,ry_sr,z_sr,rx_dif,ry_dif,z_dif, &
                  rx_c,ry_c,z_c,angdif) 
! scattering probability of the direct light. 
        call diffusion(angdif,tranam,tranaa,tranal,un,secdif, &
                  secdil,fdifan,fdifl,haer,hlay,pdifd1,z_dif) 
      else 
        pdifd1 = 0. 
      end if 
      volu = siz**3. 
      if (volu < 0.) then 
        print*,'ERROR, volume 2 is negative!' 
        stop 
      end if 
! computing scattered intensity toward the line of sight voxel from the scattering voxel 
      idif2 = fldif2*pdifd1*volu 
! computing zenith angle between the scattering voxel and the line of sight voxel 
      call anglezenithal(rx_dif,ry_dif,z_dif,rx_c,ry_c,z_c,angzen)  ! computation of the zenithal angle between the scattering voxel and the line of sight voxel. 
      call angleazimutal(rx_dif,ry_dif,rx_c,ry_c,angazi)            ! computation of the azimutal angle surf refl-scattering voxel 
! subgrid obstacles 
      if ((x_dif < 1).or.(x_dif > nbx).or.(y_dif < 1).or. &
                (y_dif > nbx)) then 
      ff = 0. 
    else 
      dho = sqrt((rx_dif-rx_c)**2.+(ry_dif-ry_c)**2.) 
      ff = 0. 
      if (dho > drefle(x_dif,y_dif)) then                        ! light path to observer larger than the mean free path -> subgrid obstacles 
      angmin = pi/2.-atan2((obsH(x_dif,y_dif)+ &
                altsol(x_dif,y_dif)-z_dif),drefle(x_dif,y_dif)) 
      if (angzen < angmin) then                                ! condition subgrid obstacle scattering -> line of sight 
      ff = 0. 
    else 
      ff = ofill(x_dif,y_dif) 
    end if 
    end if 
    end if 
    hh = 1. 
! computing transmittance between the scattering voxel and the line of sight voxel 
    distd = sqrt((rx_dif-rx_c)**2.+(ry_dif-ry_c)**2.+ &
              (z_dif-z_c)**2.) 
    call transmitm(angzen,z_dif,z_c,distd,transm,tranam,tabs) 
    call transmita(angzen,z_dif,z_c,distd,haer,transa,tranaa) 
    call transmitl(angzen,z_dif,z_c, &
              distd,hlay,transl,tranal) 
! computing the solid angle of the line of sight voxel as seen from the scattering voxel 
    omega = 1./distd**2. 
    if (omega > omemax) omega = 0. 
! computation of the scattered flux reaching the line of sight voxel 
    fdif2 = idif2*omega*transm*transa*transl*(1.-ff)*hh 
! cloud contribution for double scat from a reflecting pixel 
    if (cloudt /= 0) then                                         ! line of sight voxel = cloud 
    if (cloudbase-z_c <= iz*scal) then 
      call anglezenithal(rx_c,ry_c,z_c, &
                rx_obs,ry_obs,z_obs,azcl1)                                ! zenith angle from cloud to observer 
      call anglezenithal(rx_c,ry_c,z_c, &
                rx_dif,ry_dif,z_dif,azcl2)                                ! zenith angle from source to cloud 
      doc2 = (rx_c-rx_obs)**2.+ &
                (ry_c-ry_obs)**2.+(z_c-z_obs)**2. 
      dsc2 = (rx_dif-rx_c)**2.+ &
                (ry_dif-ry_c)**2.+(z_dif-z_c)**2. 
! cloud intensity from direct illum 
      call cloudreflectance(angzen, &
                cloudt,rcloud) 
      icloud = icloud+ &
                fldif2/omega*rcloud*doc2*omefov* &
                abs(cos(azcl2)/cos(azcl1))/dsc2/pi 
    end if 
    end if 
! computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c) 
    if (omega /= 0.) then 
! scattering angle. 
      call angle3points(rx_dif,ry_dif,z_dif,rx_c,ry_c,z_c, &
                rx_obs,ry_obs,z_obs,angdif) 
! scattering probability of the direct light. 
      call diffusion(angdif,tranam,tranaa,tranal,un,secdif, &
                secdil,fdifan,fdifl,haer,hlay,pdifd2,z_c) 
    else 
      pdifd2 = 0. 
    end if 
! computing scattered intensity toward the observer from the line of sight voxel 
    idif2p = fdif2*pdifd2 
! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation 
    idif2p = idif2p*real(stepdi)*real(ndiff-ndi)/ &
              real(ndiff-ndi-nss) 
    itotrd = itotrd+idif2p 
! ******************************************************************************** 
! *  section for the calculation of the 2nd scat from the source without reflexion 
! ******************************************************************************** 
    if ((x_sr == x_s).and.(y_sr == y_s)) then                     ! beginning condition source = reflection for the computation of the source scat line of sight 
! computation of the zenithal angle between the source and the scattering voxel 
! shadow source-scattering voxel 
    call anglezenithal(rx_s,ry_s,z_s, &
              rx_dif,ry_dif,z_dif,angzen)                                 ! computation of the zenithal angle source-scattering voxel. 
! computation of the angle azimutal line of sight-scattering voxel 
    call angleazimutal(rx_s,ry_s,rx_dif, &
              ry_dif,angazi) 
! horizon blocking not a matter because some path are downward and most of them closeby 
    hh = 1. 
    angmin = pi/2.-atan2((obsH(x_s,y_s)+ &
              altsol(x_s,y_s)-z_s),drefle(x_s, &
              y_s)) 
    if (angzen < angmin) then                                  ! condition obstacle source->scattering. 
    ff = 0. 
    else 
      ff = ofill(x_s,y_s) 
    end if 
! computation of the transmittance between the source and the scattering voxel 
    distd = sqrt((rx_s-rx_dif)**2. &
              +(ry_s-ry_dif)**2. &
              +(z_s-z_dif)**2.) 
    call transmitm(angzen,z_s,z_dif, &
              distd,transm,tranam,tabs) 
    call transmita(angzen,z_s,z_dif, &
              distd,haer,transa,tranaa) 
    call transmitl(angzen,z_s,z_dif, &
              distd,hlay,transl,tranal) 
! computation of the Solid angle of the scattering unit voxel seen from the source 
    omega = 1./distd**2. 
    if (omega > omemax) omega = 0. 
    anglez = nint(180.*angzen/pi)+1 
    P_dif1 = pvalno(anglez,stype) 
! computing flux reaching the scattering voxel 
    fldif1 = lamplu(x_s,y_s,stype)*P_dif1* &
              omega*transm*transa*transl*(1.-ff)*hh 
! computing the scattering probability toward the line of sight voxel 
    if (omega /= 0.) then 
! scattering angle. 
      call angle3points (rx_s,ry_s,z_s, &
                rx_dif,ry_dif,z_dif,rx_c,ry_c,z_c, &
                angdif) 
! scattering probability of the direct light. 
      call diffusion(angdif, &
                tranam,tranaa,tranal,un,secdif,secdil, &
                fdifan,fdifl,haer,hlay,pdifd1,z_dif) 
    else 
      pdifd1 = 0. 
    end if 
    volu = siz**3. 
    if (volu < 0.) then 
      print*,'ERROR, volume 1 is negative!' 
      stop 
    end if 
! computing scattered intensity toward the line of sight voxel from the scattering voxel 
    idif1 = fldif1*pdifd1*volu 
! computing zenith angle between the scattering voxel and the line of sight voxel 
    call anglezenithal(rx_dif,ry_dif, &
              z_dif,rx_c,ry_c,z_c,angzen)                                 ! computation of the zenithal angle between the scattering voxel and the line of sight voxel. 
! computation of the azimutal angle surf refl-scattering voxel 
    call angleazimutal(rx_dif,ry_dif, &
              rx_c,ry_c,angazi) 
! subgrid obstacles 
    if ((x_dif < 1).or.(x_dif > nbx).or.(y_dif < 1).or. &
              (y_dif > nbx)) then 
    dho = sqrt((rx_dif-rx_c)**2.+(ry_dif-ry_c)**2.) 
    ff = 0. 
    else 
      ff = 0. 
      if (dho > drefle(x_dif,y_dif)) then 
        angmin = pi/2.-atan2((obsH(x_dif,y_dif) &
                  +altsol(x_dif,y_dif)-z_dif),drefle( &
                  x_dif,y_dif)) 
        if (angzen < angmin) then                             ! condition obstacles scattering->line of sight 
        ff = 0. 
      else 
        ff = ofill(x_dif,y_dif) 
      end if 
    end if 
    end if 
    hh = 1. 
! Computing transmittance between the scattering voxel and the line of sight voxel 
    distd = sqrt((rx_c-rx_dif)**2. &
              +(ry_c-ry_dif)**2. &
              +(z_c-z_dif)**2.) 
    call transmitm(angzen,z_dif,z_c, &
              distd,transm,tranam,tabs) 
    call transmita(angzen,z_dif,z_c, &
              distd,haer,transa,tranaa) 
    call transmitl(angzen,z_dif,z_c, &
              distd,hlay,transl,tranal) 
! computing the solid angle of the line of sight voxel as seen from the scattering voxel 
    omega = 1./distd**2. 
    if (omega > omemax) omega = 0. 
! computation of the scattered flux reaching the line of sight voxel 
    fldiff = idif1*omega*transm*transa*transl*(1.-ff)*hh 
! cloud contribution to the double scattering from a source 
    if (cloudt /= 0) then                                       ! line of sight voxel = cloud 
    if (cloudbase-z_c <= iz*scal) then 
      call anglezenithal(rx_c,ry_c,z_c, &
                rx_obs,ry_obs,z_obs,azcl1)                              ! zenith angle from cloud to observer 
      call anglezenithal(rx_c,ry_c,z_c, &
                rx_dif,ry_dif,z_dif,azcl2)                              ! zenith angle from source to cloud 
      doc2 = (rx_c-rx_obs)**2.+ &
                (ry_c-ry_obs)**2.+(z_c-z_obs)**2. 
      dsc2 = (rx_dif-rx_c)**2.+ &
                (ry_dif-ry_c)**2.+(z_dif-z_c)**2. 
! cloud intensity from direct illum 
      call cloudreflectance(angzen, &
                cloudt,rcloud) 
      icloud = icloud+ &
                fldiff/omega*rcloud*doc2*omefov* &
                abs(cos(azcl2)/cos(azcl1))/dsc2/pi 
    end if 
    end if 
! computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c) 
    if (omega /= 0.) then 
! scattering angle. 
      call angle3points(rx_dif,ry_dif, &
                z_dif,rx_c,ry_c,z_c,rx_obs,ry_obs, &
                z_obs,angdif) 
! scattering probability of the direct light. 
      call diffusion(angdif, &
                tranam,tranaa,tranal,un,secdif,secdil, &
                fdifan,fdifl,haer,hlay,pdifd2,z_c) 
    else 
      pdifd2 = 0. 
    end if 
! computing scattered intensity toward the observer from the line of sight voxel 
    idiff2 = fldiff*pdifd2 
! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation 
    idiff2 = idiff2*real(stepdi)*real(ndiff-ndi)/ &
              real(ndiff-ndi-nss) 
    itodif = itodif+idiff2                                        ! sum over the scattering voxels 
    end if                                                         ! end condition source = reflection for the computation of the source scat line of sight 
    end if                                                           ! end of the case scattering pos = source pos or line of sight pos 
    end if                                                             ! end diffusing celle underground 
    end do                                                               ! end of the loop over the scattering voxels. 
    end if                             ! end of the condition ou effdif > 0 
! End of 2nd scattered intensity calculations 
!=================================================================== 
 
 
 
! ********************************************************************** 
! * section refected light with single scattering 
! ********************************************************************** 
! verify if there is shadow between sr and line of sight voxel 
! zenithal angle between the reflecting surface and the line of sight voxel. 
    call anglezenithal(rx_sr,ry_sr, &
              z_sr,rx_c,ry_c,z_c,angzen) 
! computation of the azimutal angle reflect-line of sight 
    call angleazimutal(rx_sr,ry_sr, &
              rx_c,ry_c,angazi) 
    distd = sqrt((rx_sr-rx_c)**2. &
              +(ry_sr-ry_c)**2. &
              +(z_sr-z_c)**2.) 
    dho = sqrt((rx_sr-rx_c)**2. &
              +(ry_sr-ry_c)**2.) 
    if (angzen > pi/4.) then         ! 45deg. it is unlikely to have a 1km high mountain less than 1 
    call horizon(x_sr,y_sr,z_sr,dx,dy,altsol,angazi,zhoriz,dh) 
    if (dh <= dho) then 
      if (angzen-zhoriz < &
                0.00001) then                 ! the path line of sight-reflec is not below the horizon => we compute 
      hh = 1. 
    else 
      hh = 0. 
    end if                         ! end condition reflecting surf. above horizon 
    else 
      hh = 1. 
    end if 
    else 
      hh = 1. 
    end if 
    irefl = irefl1 
! case: line of sight position = Position of reflecting cell 
    if((rx_c == rx_sr).and.(ry_c == &
              ry_sr).and.(z_c == z_sr)) then 
    intind = 0. 
    else 
! obstacle 
      dho = sqrt((rx_sr-rx_c)**2. &
                +(ry_sr-ry_c)**2.) 
      ff = 0. 
      if (dho > drefle(x_sr,y_sr)) &
                then 
      angmin = pi/2.-atan2(obsH(x_sr &
                ,y_sr),drefle(x_sr,y_sr)) 
      if (angzen < angmin) then     ! condition obstacle reflected. 
      ff = 0. 
    else 
      ff = ofill(x_sr,y_sr) 
    end if 
    end if 
! computation of the transmittance between the ground surface and the line of sight voxel 
    call transmitm(angzen,z_sr, &
              z_c,distd,transm,tranam,tabs) 
    call transmita(angzen,z_sr, &
              z_c,distd,haer,transa,tranaa) 
    call transmitl(angzen,z_sr, &
              z_c,distd,hlay,transl,tranal) 
! computation of the solid angle of the line of sight voxel seen from the reflecting cell 
    omega = 1./distd**2. 
    if (omega > omemax) omega = 0. 
! computation of the flux reflected reaching the line of sight voxel 
    flindi = irefl*omega*transm* &
              transa*transl*(1.-ff)*hh        ! obstacles correction 
! cloud contribution to the reflected light from a ground pixel 
    if (cloudt /= 0) then                       ! line of sight voxel = cloud 
    if (cloudbase-z_c <= iz*scal) then 
      call anglezenithal(rx_c,ry_c,z_c, &
                rx_obs,ry_obs,z_obs,azcl1)              ! zenith angle from cloud to observer 
      call anglezenithal(rx_c,ry_c,z_c, &
                rx_sr,ry_sr,z_sr,azcl2)                 ! zenith angle from source to cloud 
! cloud intensity from direct illum 
      doc2 = (rx_c-rx_obs)**2.+ &
                (ry_c-ry_obs)**2.+(z_c-z_obs)**2. 
      dsc2 = (rx_sr-rx_c)**2.+ &
                (ry_sr-ry_c)**2.+(z_sr-z_c)**2. 
      call cloudreflectance(angzen, &
                cloudt,rcloud) 
      icloud = icloud+ &
                flindi/omega*rcloud*doc2*omefov* &
                abs(cos(azcl2)/cos(azcl1))/dsc2/pi 
    end if 
    end if 
! computation of the scattering probability of the reflected light 
    if (omega /= 0.) then 
! scattering angle. 
      call angle3points(rx_sr, &
                ry_sr,z_sr,rx_c,ry_c,z_c, &
                rx_obs,ry_obs,z_obs,angdif) 
! scattering probability of the reflected light. 
      call diffusion(angdif, &
                tranam,tranaa,tranal,un, &
                secdif,secdil,fdifan,fdifl, &
                haer,hlay,pdifin,z_c) 
    else 
      pdifin = 0. 
    end if 
! computation of the reflected intensity toward the sensor by a reflecting cell 
    intind = flindi*pdifin*(1.-ff) &
              *hh 
    end if                             ! end of the case posi reflecting cell =  line of sight voxel position 
    itotind = itotind+intind            ! Sum of the intensities of each reflecting cell. 
    end if                               ! end of the condition surface not lighted from the top. 
    end if                                   ! end of the condition reflecting cell is not on the source. 
    end if                                     ! end of the condition surface of the domain. 
    end do                                       ! end of the loop over the rows (latitu) reflecting. 
    end do                                         ! end of the loop over the column (longitude) reflecting. 
!   end of the computation of the reflected intensity 
 
!********************************************************************** 
! computation of the total intensity coming from a source to the line of sight voxel toward the sensor 
!********************************************************************** 
! In the order 1st scat; refl->1st scat; 1st scat->2nd scat, 
! refl->1st scat->2nd scat 
    isourc = intdir+itotind+itodif+itotrd           ! Sum of the intensities of a given type of source reaching the line of sight voxel. 
    isourc = isourc*scal                            ! scaling the values according to the path length in the l. of sight voxel of 1m3 
    isourc = isourc*portio                          ! correct for the field of view of the observer 
! include clouds in the total intensity 
    isourc = isourc+icloud 
 
 
    if ((itodif < 0.).or.(itotrd < 0.)) then 
      print*,intdir,itotind,itodif,itotrd 
      stop 
    end if 
 
 
 
    if (verbose == 2) then 
      print*,' Total intensity per component for type ',ntype,':' 
      print*,' source->scattering = ',intdir 
      print*,' source->reflexion->scattering = ',itotind 
      print*,' source->scattering->scattering = ',itodif 
      print*,' source->reflexion->scattering->scattering = ',itotrd 
      if (intdir*itotind*itodif*itotrd < 0.) then 
        print*,'PROBLEM! Negative intensity.' 
        stop 
      end if 
    end if 
!********************************************************************** 
! computation of the total intensity coming from all the sources of a given type 
!********************************************************************** 
    itotty = itotty+isourc                          ! Sum of the intensities all sources of the same typeand a given line of sight element 
    ITT(x_s,y_s,stype) = ITT(x_s,y_s,stype)+isourc  ! ITT stores itotty in a matrix 
    end if                                             ! end of the condition "the luminosity of the ground pixel x_s,y_s in not null". 
    end do                                               ! end the loop over the lines (latitude) of the domain (y_s). 
    end do                                                 ! end the loop over the column (longitude) of the domain (x_s). 
! end of the computation of the intensity of one source type 
    itotci = itotci+itotty                                  ! Sum of the intensities all source all type to a line of sight element 
    do x_s = imin(stype),imax(stype) 
      do y_s = jmin(stype),jmax(stype) 
        ITC(x_s,y_s) = ITC(x_s,y_s)+ITT(x_s,y_s,stype) 
      end do 
    end do 
! calculate total lamp flux matrix for all lamp types 
    do x_s = 1,nbx 
      do y_s = 1,nby 
        lpluto(x_s,y_s) = lpluto(x_s,y_s)+ &
                  lamplu(x_s,y_s,stype) 
      end do 
    end do 
    end if                                                   ! end of condition if there are any flux in that source type 
    end do                                                       ! end of the loop over the types of sources (stype). 
! end of the computation of the intensity coming from a line of sight voxel toward the sensor 
 
 
!*********************************************************************** 
! computation of the luminous flux reaching the observer 
!*********************************************************************** 
! computation of the zenithal angle between the observer and the line of sight voxel 
!======================================================================= 
    call anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs, &
              angzen)                                                   ! computation of the zenithal angle between the line of sight voxel and the observer. 
!                                                                         ! end of the case "observer at the same latitu/longitude than the source". 
! computation of the transmittance between the line of sight voxel and the observer 
    distd = sqrt((rx_c-rx_obs)**2. &
              +(ry_c-ry_obs)**2. &
              +(z_c-z_obs)**2.) 
    call transmitm(angzen,z_c,z_obs,distd,transm,tranam, &
              tabs) 
    call transmita(angzen,z_c,z_obs,distd,haer,transa, &
              tranaa) 
    call transmitl(angzen,z_c,z_obs, &
              distd,hlay,transl,tranal) 
! computation of the flux reaching the objective of the telescope from the line of sight voxel 
    fcapt = itotci*ometif*transa*transm*transl                         ! computation of the flux reaching the intrument from the line of sight voxel 
    do x_s = 1,nbx 
      do y_s = 1,nby 
        FCA(x_s,y_s) = ITC(x_s,y_s)*ometif*transa*transm* &
                  transl 
      end do 
    end do 
    if (cos(pi-angzen) == 0.) then 
      print*,'ERROR perfectly horizontal sight is forbidden' 
      stop 
    end if 
! end of the computation of the flux reaching the observer voxel from the line of sight voxel 
    ftocap = ftocap+fcapt                                       ! flux for all source all type all line of sight element 
    do x_s = 1,nbx 
      do y_s = 1,nby 
        FTC(x_s,y_s) = FTC(x_s,y_s)+FCA(x_s,y_s)                ! FTC is the array of the flux total at the sensor to identify 
        ! the contribution of each ground pixel to the total flux at the observer level 
        ! The % is simply given by the ratio FTC/ftocap 
      end do 
    end do 
! correction for the FOV to the flux reaching the intrument from the cloud voxel 
    if (cloudt /= 0) then 
! computation of the flux reaching the intrument from the cloud voxel 
      fccld = icloud*ometif*transa*transm*transl 
      fctcld = fctcld+fccld                                       ! cloud flux for all source all type all line of sight element 
    end if 
    if (verbose >= 1) print*,'added radiance =', &
              fcapt/omefov/(pi*(diamobj/2.)**2.) 
    if (verbose >= 1) print*,'radiance accumulated =', &
              ftocap/omefov/(pi*(diamobj/2.)**2.) 
    if (verbose >= 1) write(2,*) 'added radiance =', &
              fcapt/omefov/(pi*(diamobj/2.)**2.) 
    if (verbose >= 1) write(2,*) 'radiance accumulated =', &
              ftocap/omefov/(pi*(diamobj/2.)**2.) 
    end if                                                       ! end of the condition line of sight voxel inside the modelling domain 
    end if                                                           ! end condition line of sight voxel 1/stoplim 
! accelerate the computation as we get away from the sources 
    scalo = scal 
    if (scal <= 3000.)  scal = scal*1.12 
    end if 
    else 
!           print*,'End of line of sight - touching the ground' 
    end if                                                             ! line of sight not blocked by topography 
    end do                                                             ! end of the loop over the line of sight voxels. 
    fctcld = fctcld*10**(0.4*(100.-cloudfrac)*cloudslope)               ! correction for the cloud fraction (defined from 0 to 100) 
    if (prmaps == 1) then 
!          open(unit=9,file=pclf,status='unknown') 
      do x_s = 1,nbx 
        do y_s = 1,nby 
          FTC(x_s,y_s) = FTC(x_s,y_s)/ftocap                          ! Here FTC becomes the flux fraction of each pixel. The sum of FTC values over all pixels give the total flux 
        end do 
      end do 
      if (verbose == 2) then 
        print*,'Writing normalized contribution array' 
        print*,'Warning Cloud contrib. excluded from that array.' 
      end if 
!            do x_s=1,nbx 
!              do y_s=1,nby 
!                write(9,*) x_s,y_s,FTC(x_s,y_s)                           ! emettrice au sol, c'est un % par unite of watt installes 
!              enddo 
!            enddo 
      call twodout(nbx,nby,pclimg,FTC) 
!          close(unit=9) 
! creation of gnuplot file. To visualize, type gnuplot and then 
! load 'BASENAME_pcl.gplot' 
!          open(unit=9,file=pclgp,status='unknown') 
!            write(9,*) 'sand dgrid3d',nbx,',',nby 
!            write(9,*) 'sand hidden3d' 
!            write(9,*) 'sand pm3d' 
!            write(9,*) 'splot "'//basenm(1:lenbase)//'_pcl.txt" 
!     +      with dots' 
!          close(unit=9) 
    end if                                                             ! end of condition for creating contrib and sensit maps 
    end if                                                             ! end of scattered light 
 
! End of calculation of the scattered radiances 
! ================================= 
 
    if (verbose >= 1) print*,'====================================== &
              ===============' 
    print*,'         Direct irradiance from sources (W/m**2/nm)' 
    write(*,2001)  irdirect 
    print*,'       Direct irradiance from reflexion (W/m**2/nm)' 
    write(*,2001)  irrdirect 
    print*,'         Direct radiance from sources (W/str/m**2/nm)' 
    write(*,2001)  direct 
    print*,'         Direct radiance from reflexion (W/str/m**2/nm)' 
    write(*,2001)  rdirect 
    print*,'             Cloud radiance (W/str/m**2/nm)' 
    write(*,2001) fctcld/omefov/(pi*(diamobj/2.)**2.) 
    print*,'            Diffuse radiance (W/str/m**2/nm)' 
    write(*,2001) (ftocap+fctcld)/omefov/(pi*(diamobj/2.)**2.) 
    if (verbose >= 1) write(2,*) '================================== &
              =================' 
    write(2,*) '     Direct irradiance from sources (W/m**2/nm)' 
    write(2,2001)  irdirect 
    write(2,*) '     Direct irradiance from reflexion (W/m**2/nm)' 
    write(2,2001)  irrdirect 
    write(2,*) '     Direct radiance from sources (W/str/m**2/nm)' 
    write(2,2001)  direct 
    write(2,*) '     Direct radiance from reflexion (W/str/m**2/nm)' 
    write(2,2001)  rdirect 
    write(2,*) '           Cloud radiance (W/str/m**2/nm)         ' 
    write(2,2001) fctcld/omefov/(pi*(diamobj/2.)**2.) 
    write(2,*) '         Diffuse radiance (W/str/m**2/nm)          ' 
    write(2,2001) (ftocap+fctcld)/omefov/(pi*(diamobj/2.)**2.) 
    close(2) 
 2001 format('                   ',E10.3E2) 
    stop 
end program illumina                                                    ! beginning 
!*********************************************************************************************************************** 
!*                                                                                                                     * 
!*                                         end of the programme                                                        * 
!*                                                                                                                     * 
!*********************************************************************************************************************** 
