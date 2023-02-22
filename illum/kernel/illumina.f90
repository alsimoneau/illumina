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

PROGRAM illumina
  IMPLICIT NONE

!=======================================================================
!     Variables declaration
!=======================================================================

  INTEGER :: width, nzon                                                 ! matrix dimension in length/width and height
  PARAMETER(width=512, nzon=256)
  INTEGER :: iun, ideux
  REAL :: pi, pix4
  REAL :: zero, un                                                       ! value of 0.0 and 1.0
  INTEGER :: verbose                                                     ! verbose = 1 to have more print out, 0 for silent
  PARAMETER(pi=3.141592654)
  PARAMETER(pix4=4.0 * pi)
  CHARACTER*72 :: mnaf                                                   ! terrain elevation file
  CHARACTER*72 :: diffil                                                 ! aerosol file
  CHARACTER*72 :: outfile                                                ! results file
  CHARACTER*72 :: pclf, pclgp                                            ! files containing contribution and sensitivity maps
  CHARACTER*72 :: pclimg, pcwimg
  CHARACTER*72 :: basenm                                                 ! base name of files
  INTEGER :: lenbase                                                     ! length of the base name of the experiment
  REAL :: lambda, pressi, drefle(width, width)                           ! wavelength (nanometer), atmospheric pressure (kpa), mean free path to the ground (meter).
  REAL :: reflsiz                                                        ! size of the reflecting surface
  INTEGER :: ntype                                                       ! number of light source types or zones considered
  REAL :: largx                                                          ! width (x axis) of the modeling domain (meter)
  REAL :: largy                                                          ! length (y axis) of the modeling domain (meter)
  INTEGER :: nbx, nby                                                    ! number of pixels in the modeling domain
  REAL :: val2d(width, width)                                            ! temporary input array 2d
  REAL :: altsol(width, width)                                           ! ground elevation (meter)
  REAL :: srefl                                                          ! ground reflectance
  INTEGER :: stype                                                       ! source type or zone index
  CHARACTER*72 :: pafile, lufile, alfile, ohfile, odfile, offile         ! files related to light sources and obstacles (photometric function of the sources (sr-1), flux (w), height (m), obstacles c                                                               ! height (m), obstacle distance (m), obstacle filling factor (0-1).
  REAL :: lamplu(width, width, nzon)                                     ! source fluxes
  REAL :: lampal(width, width)                                           ! height of the light sources relative to the ground (meter)
  REAL :: pval(181, nzon), pvalto, pvalno(181, nzon)                     ! values of the angular photometry functions (unnormalized, integral, normalized)
  REAL :: dtheta                                                         ! angle increment of the photometric function of the sources
  REAL :: dx, dy, dxp, dyp                                               ! width of the voxel (meter)
  INTEGER :: boxx, boxy                                                  ! reflection window size (pixels)
  REAL :: fdifa(181), fdifan(181)                                        ! aerosol scattering functions (unnormalized and normalized)
  REAL :: anglea(181)                                                    ! scattering angle (degree)
  REAL :: secdif                                                         ! contribution of the scattering to the extinction
  REAL :: inclix(width, width)                                           ! tilt of the ground pixel along x (radian)
  REAL :: incliy(width, width)                                           ! tilt of the ground pixel along y (radian)
  INTEGER :: x_obs, y_obs                                                ! position of the observer (integer)
  REAL :: rx_obs, ry_obs
  REAL :: z_o                                                            ! observer height relative to the ground (meter)
  REAL :: z_obs                                                          ! height of the observer (meter) to the vertical grid scale
  INTEGER :: ncible, icible                                              ! number of line of sight voxels, number loops over the voxels
  INTEGER :: x_c, y_c                                                    ! position of the line of sight voxel (integer)
  REAL :: rx_c, ry_c
  REAL :: z_c                                                            ! height of the line of sight voxel (meter)
  INTEGER :: dirck                                                       ! test for the position of the source (case source = line of sight voxel)
  INTEGER :: x_s, y_s, x_sr, y_sr, x_dif, y_dif                          ! positions of the source, the reflecting surface, and the scattering voxels
  REAL :: z_s, z_sr, z_dif                                               ! heights of the source, the reflecting surface, and the scattering voxel (metre).
  REAL :: rx_s, ry_s, rx_sr, ry_sr, rx_dif, ry_dif
  REAL :: angzen, ouvang                                                 ! zenithal angle between two voxels (radians) and opening angle of the solid angle in degrees.
  INTEGER :: anglez                                                      ! emitting zenithal angle from the luminaire.
  REAL :: p_dir, p_indir, p_dif1                                         ! photometric function of the light sources (direct,indirect,scattered)
  REAL :: transa, transm, transl                                         ! transmittance between two voxels (aerosols,molecules,particle layer).
  REAL :: taua                                                           ! aerosol optical depth @ 500nm.
  REAL :: alpha                                                          ! angstrom coefficient of aerosol aod
  REAL(8) xc, yc, zc, xn, yn, zn                                         ! position (meter) of the elements (starting point, final point) :: for the calculation of the solid angle.
  REAL(8) :: r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z  ! components of the vectors used in the solid angle calculation routine.
  REAL :: omega                                                          ! solid angles
  REAL :: fldir                                                          ! flux coming from a source (watt).
  REAL :: flindi                                                         ! flux coming from a reflecting ground element (watt).
  REAL :: fldiff                                                         ! flux coming from a scattering voxel (watt).
  REAL :: angdif                                                         ! scattering angle.
  REAL :: pdifdi, pdifin, pdifd1, pdifd2                                 ! scattering probability (direct,indirect,1st and 2nd order of scattering
  REAL :: intdir                                                         ! direct intensity toward the sensor from a scattering voxel.
  REAL :: intind                                                         ! contribution of the reflecting cell to the reflected intensity toward the sensor.
  REAL :: itotind                                                        ! total contribution of the source to the reflected intensity toward the sensor.
  REAL :: idiff2                                                         ! contribution of the scattering voxel to the scattered intensity toward the sensor.
  REAL :: itodif                                                         ! total contribution of the source to the scattered intensity toward the sensor.
  REAL :: isourc                                                         ! total contribution of the source to the intensity from a line of sight voxel toward the sensor.
  REAL :: itotty                                                         ! total contribution of a source type to the intensity coming from a line of sight voxel toward the sensor.
  REAL :: itotci                                                         ! total intensity from a line of sight voxel toward the sensor.
  REAL :: itotrd                                                         ! total intensity a voxel toward the sensor after reflexion and double scattering.
  REAL :: flcib                                                          ! flux reaching the observer voxel from a line of sight voxel.
  REAL :: fcapt                                                          ! flux reaching the observer voxel from all fov voxels in a given model level
  REAL :: ftocap                                                         ! total flux reaching the observer voxel
  REAL :: haut                                                           ! haut (negative indicate that the surface is lighted from inside the ground. i.e. not considered in the calculation
  REAL :: epsilx, epsily                                                 ! tilt of the ground pixel
  REAL :: flrefl                                                         ! flux reaching a reflecting surface (watts).
  REAL :: irefl, irefl1                                                  ! intensity leaving a reflecting surface toward the line of sight voxel.
  REAL :: effdif                                                         ! distance around the source voxel and line of sight voxel considered to compute the 2nd order of scattering.
  REAL :: zondif(3000000, 3)                                             ! array for the scattering voxels positions
  INTEGER :: ndiff, idi                                                  ! number of scattering voxels, counter of the loop over the scattering voxels
  INTEGER :: stepdi                                                      ! scattering step to speedup the calculation e.g. if =2 one computation over two will be done
  INTEGER :: ssswit                                                      ! activate double scattering (1=yes, 0 = no)
  INTEGER :: fsswit                                                      ! activate first scattering (1=yes, 0 = no)
!                                                                        ! by default the value is 1 but it can be larger
!                                                                        ! when we resume a previous interrupted calculation.
  REAL :: fldif1, fldif2                                                 ! flux reaching a scattering voxel.
  REAL :: fdif2                                                          ! flux reaching the line of sight voxel after reflexion > scattering
  REAL :: idif1, idif2, idif2p                                           ! intensity toward a line of sight voxel from a scattering voxel (without and with reflexion).
  REAL :: portio                                                         ! ratio of voxel surface to the solid angle of the sensor field of view.
  REAL :: dis_obs                                                        ! distance between the line of sight and the observer.
  REAL :: ometif                                                         ! solid angle of the telescope objective as seen from the line of sight voxel
  REAL :: omefov                                                         ! solid angle of the spectrometer slit.
  REAL :: angvis, azim                                                   ! viewing angles of the sensor.
!                                                                        ! Useful for the calculation of the lambertian reflectance.
  REAL :: nbang                                                          ! for the averaging of the photometric function
  REAL :: obsh(width, width), angmin                                     ! averaged height of the sub-grid obstacles, minimum angle under wich
!                                                                        ! a light ray cannot propagate because it is blocked by a sub-grid obstable
  REAL :: ofill(width, width)                                            ! fill factor giving the probability to hit an obstacle when pointing in its direction real 0-1
  INTEGER :: naz, na
  REAL :: itt(width, width, nzon)                                        ! total intensity per type of lamp
  REAL :: itc(width, width)                                              ! total intensity per line of sight voxel
  REAL :: ftc(width, width)                                              ! fraction of the total flux at the sensor level
  REAL :: fca(width, width)                                              ! sensor flux array
  REAL :: lpluto(width, width)                                           ! total luminosity of the ground cell for all lamps
  CHARACTER*3 :: lampno                                                  ! lamp number string
  INTEGER :: imin(nzon), imax(nzon), jmin(nzon), jmax(nzon)              ! x and y limits containing a type of lamp
  REAL :: angazi                                                         ! azimuth angle between two points in rad, max dist for the horizon determination
  REAL :: latitu                                                         ! approximate latitude of the domain center
  INTEGER :: prmaps                                                      ! flag to enable the tracking of contribution and sensitivity maps
  INTEGER :: cloudt                                                      ! cloud type 0=clear, 1=thin cirrus/cirrostratus, 2=thick cirrus/cirrostratus, 3 = altostratus/altocumulus,
  ! 4=Stratocumulus/stratus, 5 = Cumulus/Cumulonimbus
  REAL :: cloudslope                                                     ! slope of the radiance dependency on the cloud fraction (in percentage) according to
  ! Sciezoras 2020 the slope vary depending on the level of LP and how it is distributed.
  ! We decided instead to simplify this by using an average slope of -0.013.0
  ! Rad = Rad_100 * 10**(0.4*(100-cloudfrac)*cloudslope) this equation is derived from
  ! Tomasz Sciezor, The impact of clouds on the brightness of the night sky, Journal of
  ! Quantitative Spectroscopy & Radiative Transfer (2020),
  ! doi: https://doi.org/10.1016/j.jqsrt.2020.106962
  REAL :: cloudfrac                                                      ! cloud fraction in percentage
  INTEGER :: xsrmi, xsrma, ysrmi, ysrma                                  ! limits of the loop valeur for the reflecting surfaces
  REAL :: rcloud                                                         ! cloud relfectance
  REAL :: icloud                                                         ! cloud reflected intensity
  REAL :: fccld                                                          ! correction for the fov to the flux reaching the intrument from the cloud voxel
  REAL :: fctcld                                                         ! total flux from cloud at the sensor level
  REAL :: totlu(nzon)                                                    ! total flux of a source type
  REAL :: stoplim                                                        ! stop computation when the new voxel contribution is less than 1/stoplim of the cumulated flux
  REAL :: ff, ff2, hh                                                    ! temporary obstacle filling factor and horizon blocking factor
  REAL :: cloudbase, cloudtop                                            ! cloud base and top altitude (m)
  REAL :: distd                                                          ! dizhorizstance to compute the scattering probability
  REAL :: volu                                                           ! volume of a voxel
  REAL :: scal                                                           ! stepping along the line of sight
  REAL :: scalo                                                          ! previous value of scal
  REAL :: siz                                                            ! resolution of the 2nd scat grid in meter
  REAL :: angvi1, angaz1, angze1                                         ! viewing angles in radian
  REAL :: ix, iy, iz                                                     ! base vector of the viewing (length = 1)
  REAL :: dsc2, doc2                                                     ! square of the path lengths for the cloud contribution
  REAL :: azcl1, azcl2                                                   ! zenith angle from the (source, refl surface, or scattering voxel) to line of path and observer to line p.
  REAL :: dh, dho                                                        ! distance of the horizon limit
  INTEGER :: step                                                        ! skiping 2nd scat on 1 dim
  REAL :: omemax                                                         ! max solid angle allowed
  REAL :: flcld(width, width)                                            ! flux crossing a low cloud
  REAL :: ds1, ds2, ds3, dss                                             ! double scattering distances
  INTEGER :: nss                                                         ! number of skipped 2nd scat elements
  INTEGER :: ndi                                                         ! number of cell under ground
  REAL :: diamobj                                                        ! instrument objective diameter
  INTEGER :: i, j, id, jd
  REAL :: tranam, tranaa                                                 ! atmospheric transmittancess of a path (molecular, aerosol)
  REAL :: zhoriz                                                         ! zenith angle of the horizon
  REAL :: direct                                                         ! direct radiance from sources on a surface normal to the line of sight (no scattering)
  REAL :: rdirect                                                        ! direct radiance from a reflecting surface on a surface normal to the line of sight (no scattering)
  REAL :: irdirect                                                       ! direct irradiance from sources on a surface normal to the line of sight (no scattering)
  REAL :: irrdirect                                                      ! direct irradiance from a reflecting surface on a surface normal to the line of sight (no scattering)
  REAL :: dang                                                           ! angle between the line of sight and the direction of a source
  REAL :: dzen                                                           ! zenith angle of the source-observer line
  REAL :: ddir_obs                                                       ! distance between the source and the observer
  REAL :: rx, ry, rz                                                     ! driving vector for the calculation of the projection angle for direct radiance. it is 20km long
  REAL :: dfov                                                           ! field of view in degrees for the calculation of the direct radiance this number will be a kind of smoothing effect. the angular grid resolution to create a direct radiance panorama should be finer than that number
  REAL :: fo(width, width)                                               ! flux correction factor for obstacles
  REAL :: thetali(width, width)                                          ! limit angle for the obstacles blocking of viirs
  INTEGER :: viirs(width, width)                                         ! viirs flag 1=yes 0 = no
  CHARACTER*72 :: vifile                                                 ! name of the viirs flag file
  REAL :: dh0, dhmax                                                     ! horizontal distance along the line of sight and maximum distance before beeing blocked by topography
  CHARACTER*72 :: layfile                                                ! filename of the optical properties of the particle layer
  REAL :: layaod                                                         ! 500 nm aod of the particle layer
  REAL :: layalp                                                         ! spectral exponent of the aod for the particle layer
  REAL :: hlay                                                           ! exponential vertical scale height of the particle layer
  REAL :: secdil                                                         ! scattering/extinction ratio for the particle layer
  REAL :: fdifl(181)                                                     ! scattering phase function of the particle layer
  REAL :: tranal                                                         ! top of atmos transmission of the particle layer
  REAL :: haer                                                           ! exponential vertical scale height of the background aerosol layer
  REAL :: bandw                                                          ! bandwidth of the spectral bin
  REAL :: tabs                                                           ! toa transmittance related to molecule absorption
  INTEGER :: obsobs                                                      ! flag to activate the direct light obstacle blocking aroud the observer.
  verbose = 1                                                            ! Very little printout=0, Many printout = 1, even more = 2
  diamobj = 1.0                                                           ! A dummy value for the diameter of the objective of the instrument used by the observer.
  volu = 0.0
  zero = 0.0
  un = 1.0
  ff = 0.0
  ff2 = 0.0
  step = 1
  ncible = 1024
  stepdi = 1
  cloudslope = -0.013
  cloudfrac = 100.0
  IF (verbose >= 1) THEN
    PRINT *, 'Starting ILLUMINA computations...'
  END IF
  ! reading of the fichier d'entree (illumina.in)
  PRINT *, 'Reading illumina.in input file'
  OPEN (unit=1, file='illumina.in', status='old')
  READ (1, *)
  READ (1, *) basenm
  READ (1, *) dx, dy
  READ (1, *) diffil
  READ (1, *) layfile, layaod, layalp, hlay
  READ (1, *) ssswit
  READ (1, *) fsswit
  READ (1, *) lambda, bandw
  READ (1, *) srefl
  READ (1, *) pressi
  READ (1, *) taua, alpha, haer
  READ (1, *) ntype
  READ (1, *) stoplim
  READ (1, *)
  READ (1, *) x_obs, y_obs, z_o
  READ (1, *) obsobs
  READ (1, *) angvis, azim
  READ (1, *) dfov
  READ (1, *)
  READ (1, *)
  READ (1, *)
  READ (1, *) reflsiz
  READ (1, *) cloudt, cloudbase, cloudfrac
  READ (1, *)
  CLOSE (1)
  IF (angvis > 90.0) THEN
    PRINT *, 'Error: elevation angle larger than 90 deg'
    STOP
  END IF
  IF (angvis < -90.0) THEN
    PRINT *, 'Error: elevation angle smaller than -90 deg'
    STOP
  END IF
  ! conversion of the geographical viewing angles toward the cartesian
  ! angle we assume that the angle in the file illumina.in
  ! is consistent with the geographical definition
  ! geographical, azim=0 toward north, 90 toward east, 180 toward south, etc
  ! cartesian, azim=0 toward east, 90 toward north, 180 toward west, etc
  azim = 90.0 - azim
  IF (azim < 0.0) azim = azim + 360.0
  IF (azim >= 360.0) azim = azim - 360.0
  angvi1 = (pi * angvis) / 180.0
  angze1 = pi / 2.0 - angvi1
  angaz1 = (pi * azim) / 180.0
  ! viewing vector components
  ix = (SIN((pi / 2.0) - angvi1)) * (COS(angaz1))
  iy = (SIN((pi / 2.0) - angvi1)) * (SIN(angaz1))
  iz = (SIN(angvi1))
  dfov = (dfov * pi / 180.0) / 2.0
  siz = 2500.0
  IF (ssswit == 0) THEN
    effdif = 0.0
  ELSE
    effdif = 40000.0
  END IF
  scal = 19.0
  scalo = scal
  boxx = NINT(reflsiz / dx)  ! Number of column to consider left/right of the source for the reflection.
  boxy = NINT(reflsiz / dy)  ! Number of column to consider up/down of the source for the reflection.

  ! omemax: exclude calculations too close (<10m) this is a sustended angle of 1 deg.
  ! the calculated flux is highly sensitive to that number for a very high
  ! pixel resolution (a few 10th of meters). We assume anyway that somebody
  ! observing the sky will never lies closer than that distance to a
  ! light fixture. This number is however somehow subjective and that means
  ! that the value of sky brightness near sources will be affected by this
  ! choice
  omemax = 1.0 / ((10.0)**2.0)
  IF (verbose > 0) THEN
    PRINT *, '2nd order scattering grid = ', siz, 'm'
    PRINT *, '2nd order scattering radius = ', effdif, 'm'
    PRINT *, 'Pixel size = ', dx, ' x ', dy
    PRINT *, 'Maximum radius for reflection = ', reflsiz
  END IF

  ! computing the actual AOD at the wavelength lambda
  IF (verbose >= 1) PRINT *, '500nm aod = ', taua, '500nm angstrom coeff.= ', alpha
  taua = taua * (lambda / 500.0)**(-1.0 * alpha)
  layaod = layaod * (lambda / 500.0)**(-1.0 * layalp)
  !  determine the Length of basenm
  lenbase = INDEX(basenm, ' ') - 1
  mnaf = basenm(1:lenbase)//'_topogra.bin'   ! determine the names of input and output files
  outfile = basenm(1:lenbase)//'.out'
  pclf = basenm(1:lenbase)//'_pcl.txt'
  pclimg = basenm(1:lenbase)//'_pcl.bin'
  pcwimg = basenm(1:lenbase)//'_pcw.bin'
  pclgp = basenm(1:lenbase)//'_pcl.gplot'

  ! opening output file
  OPEN (unit=2, file=outfile, status='unknown')
  WRITE (2, *) "ILLUMINA version __version__"
  WRITE (2, *) 'FILE USED:'
  WRITE (2, *) mnaf, diffil
  PRINT *, 'Wavelength (nm):', lambda, &
    ' Aerosol optical depth:', taua
  WRITE (2, *) 'Wavelength (nm):', lambda, &
    ' Aerosol optical depth:', taua
  WRITE (2, *) '2nd order scattering radius:', effdif, ' m'
  PRINT *, '2nd order scattering radius:', effdif, ' m'
  WRITE (2, *) 'Observer position (x,y,z)', x_obs, y_obs, z_o
  PRINT *, 'Observer position (x,y,z)', x_obs, y_obs, z_o
  WRITE (2, *) 'Elevation angle:', angvis, ' azim angle (counterclockwise from east)', azim
  PRINT *, 'Elevation angle:', angvis, ' azim angle (counterclockwise from east)', azim

  ! Initialisation of the arrays and variables
  IF (verbose >= 1) PRINT *, 'initializing variables...'
  IF (cloudt == 0) THEN
    cloudbase = 1000000000.0
  END IF
  prmaps = 1
  iun = 0
  ideux = 1
  icloud = 0.0
  val2d = 0.0
  altsol = 0.0
  obsH = 0.0
  viirs = 0
  ofill = 0.0
  inclix = 0.0
  incliy = 0.0
  lpluto = 0.0
  ITC = 0.0
  FTC = 0.0
  FCA = 0.0
  flcld = 0.0
  lamplu = 0.0
  lampal = 0.0
  ITT = 0.0
  fdifa = 0.0
  fdifan = 0.0
  fdifl = 0.0
  anglea = 0.0
  pval = 0.0
  pvalno = 0.0
  zondif = 1.0
  idif1 = 0.0
  idif2 = 0.0
  fdif2 = 0
  idif2p = 0.0
  fldir = 0.0
  flindi = 0.0
  fldiff = 0.0
  pdifdi = 0.0
  pdifin = 0.0
  pdifd1 = 0.0
  pdifd2 = 0.0
  intdir = 0.0
  intind = 0.0
  idiff2 = 0.0
  angmin = 0.0
  isourc = 0.0
  itotty = 0.0
  itotci = 0.0
  itotrd = 0.0
  flcib = 0.0
  flrefl = 0.0
  irefl = 0.0
  irefl1 = 0.0
  fldif1 = 0.0
  fldif2 = 0.0
  portio = 0.0
  fccld = 0.0
  fctcld = 0.0
  ometif = 0.0
  omefov = 0.0
  hh = 1.0

  ! determine the 2nd scattering zone
  IF (ssswit /= 0) THEN
    CALL zone_diffusion(effdif, &
                        zondif, ndiff, stepdi, siz)
    dss = 1.0 * siz
    IF (verbose > 0) THEN
      PRINT *, '2nd order scattering grid points =', ndiff
      PRINT *, '2nd order scattering smoothing radius =', dss, 'm'
    END IF
  END IF

  ! determination of the vertical atmospheric transmittance
  ! tranam and tranaa are the top of atmosphere transmittance (molecules and aerosols)
  CALL transtoa(lambda, bandw, taua, layaod, pressi, tranam, tranaa, tranal, tabs)

  ! reading of the environment variables

  ! reading of the elevation file
  CALL twodin(nbx, nby, mnaf, altsol)

  ! computation of the tilt of the pixels along x and along y
  ! beginning of the loop over the column (longitude) of the domain.
! Computation of tilt along x of the surface
  inclix(:, 1) = ATAN((altsol(:, 2) - altsol(:, 1)) / dx)
  inclix(:, 2:nby - 1) = ATAN((altsol(:, 3:nby) - altsol(:, 1:nby - 2)) / (2.0 * dx))
  inclix(:, nby) = ATAN((altsol(:, nby - 1) - altsol(:, nby)) / dx)

  ! Computation of tilt along y of the surface
  incliy(1, :) = ATAN((altsol(2, :) - altsol(1, :)) / dy)
  incliy(2:nbx - 1, :) = ATAN((altsol(3:nbx, :) - altsol(1:nbx - 2, :)) / (2.0 * dy))
  incliy(nbx, :) = ATAN((altsol(nbx - 1, :) - altsol(nbx, :)) / dy)

  ! reading of the values of P(theta), height, luminosities and positions
  ! of the sources, obstacle height and distance
  ohfile = basenm(1:lenbase)//'_obsth.bin'
  odfile = basenm(1:lenbase)//'_obstd.bin'
  alfile = basenm(1:lenbase)//'_altlp.bin'
  offile = basenm(1:lenbase)//'_obstf.bin'
  vifile = 'origin.bin'
  ! one degree
  dtheta = .017453293
  ! reading lamp heights
  CALL twodin(nbx, nby, alfile, lampal)
  ! reading subgrid obstacles average height
  CALL twodin(nbx, nby, ohfile, obsH)
  ! reading subgrid obstacles average distance
  CALL twodin(nbx, nby, odfile, drefle)
  drefle = drefle / 2.0
  drefle = MERGE(drefle, dx, drefle == 0.0)
  ! reading subgrid obstacles filling factor
  CALL twodin(nbx, nby, offile, ofill)
  ! reading viirs flag
  CALL twodin(nbx, nby, vifile, val2d)
  viirs = INT(val2d)

  ! reading of the scattering parameters for background aerosols
  ! opening file containing the scattering parameters
  OPEN (unit=1, file=diffil, status='old')
  ! the scattering / extinction ratio
  READ (1, *) secdif
  READ (1, *)
  DO i = 1, 181
    ! reading of the scattering functions
    READ (1, *) anglea(i), fdifa(i)
  END DO
  ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle
  fdifan = fdifa / pix4
  CLOSE (1)

  ! reading scattering parameters of particle layer
  ! opening file containing the scattering parameters
  OPEN (unit=1, file=layfile, status='old')
  ! the scattering / extinction ratio of particle layer
  READ (1, *) secdil
  READ (1, *)
  DO i = 1, 181
    ! reading of the scattering functions of the particle layer
    READ (1, *) anglea(i), fdifl(i)
  END DO
  ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle
  fdifl = fdifl / pix4
  CLOSE (1)

  ! Some preliminary tasks

  ! beginning of the loop 1 for the nzon types of sources.
  DO stype = 1, ntype
    imin(stype) = nbx
    jmin(stype) = nby
    imax(stype) = 1
    jmax(stype) = 1
    pvalto = 0.0
    ! support of nzon different sources (3 digits)
    WRITE (lampno, '(I3.3)') stype
    ! setting the file name of angular photometry.
    pafile = basenm(1:lenbase)//'_fctem_'//lampno//'.dat'
    ! setting the file name of the luminosite of the cases.
    lufile = basenm(1:lenbase)//'_lumlp_'//lampno//'.bin'

    ! reading photometry files
    ! opening file pa#.dat, angular photometry.
    OPEN (UNIT=1, FILE=pafile, status='OLD')
    DO i = 1, 181
      READ (1, *) pval(i, stype)
      ! Sum of the values of the  photometric function
      ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale (i-1) x 1 degrees).
      pvalto = pvalto + pval(i, stype) * 2.0 * pi * &
               SIN(REAL(i - 1) * dtheta) * dtheta
    END DO
    CLOSE (1)

    ! normalisation of the photometric function.
    IF (pvalto /= 0.0) pvalno(:, stype) = pval(:, stype) / pvalto

    ! reading luminosity files
    CALL twodin(nbx, nby, lufile, val2d)
    DO i = 1, nbx
      DO j = 1, nby
        IF (val2d(i, j) < 0.0) THEN
          PRINT *, '***Negative lamp flux!, stopping execution'
          STOP
        END IF
      END DO
    END DO

    ! searching of the smallest rectangle containing the zone
    ! of non-null luminosity to speedup the calculation
    imin(stype) = 1
    imax(stype) = nbx
    jmin(stype) = 1
    jmax(stype) = nby
    DO i = 1, nbx
      DO j = 1, nby
        IF (val2d(i, j) /= 0.0) THEN
          imin(stype) = MAX(imin(stype), i - 2)
          imax(stype) = MIN(imax(stype), i + 2)
          jmin(stype) = MAX(jmin(stype), j - 2)
          jmax(stype) = MIN(jmax(stype), j + 2)
        END IF
      END DO
    END DO

    lamplu(:, :, stype) = val2d(:, :)
    WHERE (viirs == 1)
      lamplu(:, :, stype) = lamplu(:, :, stype) / (tranam * tranaa * tranal)

      thetali = ATAN2(drefle, obsH)
      WHERE (thetali < 70.0 * pi / 180.0)
        Fo = (1.0 - COS(70.0 * pi / 180.0)) &
             / (1.0 - ofill * COS(thetali) &
                + (ofill - 1.0) * COS(70.0 * pi / 180.0))
        lamplu(:, :, stype) = lamplu(:, :, stype) * Fo
      END WHERE
    END WHERE
    ! end of the loop 1 over the nzon types of sources.
  END DO

  dy = dx
  ! solid angle of the spectrometer slit on the sky. Here we only need a small value
  omefov = 0.00000001
  ! z_obs = the local observer elevation plus the height of observation above ground (z_o)
  z_obs = z_o + altsol(x_obs, y_obs)
  rx_obs = REAL(x_obs) * dx
  ry_obs = REAL(y_obs) * dy
  IF (z_obs == 0.0) z_obs = 0.001

  ! computation of the Width along x of the case.
  largx = dx * REAL(nbx)
  ! computation of the Width along y of the case.
  largy = dy * REAL(nby)

  WRITE (2, *) 'Width of the domain [NS](m):', largx, '#cases:', nbx
  WRITE (2, *) 'Width of the domain [EO](m):', largy, '#cases:', nby
  WRITE (2, *) 'Size of a cell (m):', dx, ' X ', dy
  WRITE (2, *) 'latitu center:', latitu

  ! initialize the total direct radiance from sources to observer
  direct = 0.0
  ! initialize the total reflected radiance from surface to observer
  rdirect = 0.0
  ! initialize the total direct irradiance from sources to observer
  irdirect = 0.0
  ! initialize the total reflected irradiance from surface to observer
  irrdirect = 0.0

! =================================
! Calculation of the direct radiances

  IF (verbose >= 1) PRINT *, ' calculating obtrusive light...'
  ! beginning of the loop over the source types.
  DO stype = 1, ntype
    ! check if there are any flux in that source type otherwise skip this lamp
    IF (totlu(stype) /= 0.0) THEN
      IF (verbose >= 1) PRINT *, ' turning on lamps', stype
      IF (verbose >= 1) WRITE (2, *) ' turning on lamps', &
        stype
      ! beginning of the loop over the column (longitude the) of the domain.
      DO x_s = imin(stype), imax(stype)
        ! beginning of the loop over the rows (latitud) of the domain.
        DO y_s = jmin(stype), jmax(stype)
          intdir = 0.0
          itotind = 0.0
          itodif = 0.0
          itotrd = 0.0
          isourc = 0.0
          rx_s = REAL(x_s) * dx
          ry_s = REAL(y_s) * dy
          ! if the luminosite of the case is null, the program ignore this case.
          IF (lamplu(x_s, y_s, stype) /= 0.0) THEN
            ! Definition of the position (metre) vertical of the source.
            z_s = (altsol(x_s, y_s) + lampal(x_s, y_s))

! *********************************************************************************************************
! calculation of the direct radiance of sources falling on a surface perpendicular
! to the viewing angle Units of W/nm/m2/sr
! *********************************************************************************************************

            rx = rx_obs + 20000.0 * ix
            ry = ry_obs + 20000.0 * iy
            rz = z_obs + 20000.0 * iz
            dho = SQRT((rx_obs - rx_s)**2.0 &
                       + (ry_obs - ry_s)**2.0)
            IF ((dho > 0.0) .AND. (z_s /= z_obs)) THEN
              ! zenithal angle source-observer
              CALL anglezenithal(rx_obs, ry_obs, z_obs, rx_s, ry_s, z_s, dzen)
              ! computation of the angle azimutal direct line of sight-source
              CALL angleazimutal(rx_obs, ry_obs, rx_s, ry_s, angazi)

              ! 45deg. it is unlikely to have a 1km high mountain less than 1
              hh = 1.0
              IF (angzen > (pi / 4.0)) THEN
                CALL horizon(x_sr, y_sr, z_sr, dx, dy, altsol, angazi, zhoriz, dh)
                hh = MERGE(0.0, 1.0, dh <= dho .AND. angzen - zhoriz < 0.00001)
              END IF

              ff = 0.0
              IF (x_dif >= 1 .AND. x_dif <= nbx .AND. y_dif >= 1 .AND. y_dif <= nbx) THEN
                dho = SQRT((rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0)
                IF (dho <= drefle(x_dif, y_dif)) THEN
                  angmin = pi / 2.0 - ATAN2((obsH(x_dif, y_dif) + altsol(x_dif, y_dif) - z_dif), drefle(x_dif, y_dif))
                  IF (angzen >= angmin) THEN
                    ff = ofill(x_dif, y_dif)
                  END IF
                END IF
              END IF

              ! zenithal angle source-observer
              CALL anglezenithal(rx_s, ry_s, z_s, rx_obs, ry_obs, z_obs, dzen)
              ff2 = 0.0
              ! light path from source larger than the mean free path -> subgrid obstacles
              IF (dho > drefle(x_s, y_s)) THEN
                angmin = pi / 2.0 - ATAN2((altsol(x_s, y_s) + obsH(x_s, y_s) - z_s), &
                                          drefle(x_s, y_s))
                ! condition sub-grid obstacles direct.
                IF (dzen >= angmin) THEN
                  ff2 = ofill(x_s, y_s)
                END IF
                ! end light path to the observer larger than mean free path
              END IF

              ! zenithal angle source-observer
              CALL anglezenithal(rx_obs, ry_obs, z_obs, rx_s, ry_s, z_s, dzen)

              ! projection angle of line to the lamp and the viewing angle
              CALL angle3points(rx_s, ry_s, z_s, rx_obs, ry_obs, z_obs, rx, ry, rz, dang)

              ! scattering angle.
              dang = pi - dang
              ! computation of the solid angle of the line of sight voxel seen from the source
              anglez = NINT(180.0 * (pi - dzen) / pi) + 1
              P_dir = pvalno(anglez, stype)
              ! computation of the flux direct reaching the line of sight voxel
              IF ((COS(dang) > 0.0) .AND. (dang < pi / 2.0)) THEN
                ! distance direct sight between source and observer
                ddir_obs = SQRT((rx_obs - rx_s)**2.0 + (ry_obs - ry_s)**2.0 + (z_obs - z_s)**2.0)
                ! computation of the solid angle 1m^2 at the observer as seen from the source
                omega = 1.0 * ABS(COS(dang)) / ddir_obs**2.0
                CALL transmitm(dzen, z_obs, z_s, ddir_obs, transm, tranam, tabs)
                CALL transmita(dzen, z_obs, z_s, ddir_obs, haer, transa, tranaa)
                CALL transmitl(dzen, z_obs, z_s, ddir_obs, hlay, transl, tranal)
                !check if the reflecting surface enter the field of view of the observer
                IF (dang < dfov) THEN
                  ! correction for obstacle filling factor
                  direct = direct + lamplu(x_s, y_s, stype) * transa * transm * transl * P_dir * omega &
                           * (1.0 - ff) * (1.0 - ff2) * hh / (pi * dfov**2.0)
                END IF
                ! correction for obstacle filling factor
                irdirect = irdirect + lamplu(x_s, y_s, stype) * transa * transm * transl * P_dir * omega &
                           * (1.0 - ff) * (1.0 - ff2) * hh
              END IF
            END IF

! **********************************************************************************
! * computation of the direct light toward the observer by the ground reflection   *
! **********************************************************************************

            xsrmi = x_s - boxx
            IF (xsrmi < 1) xsrmi = 1
            xsrma = x_s + boxx
            IF (xsrma > nbx) xsrma = nbx
            ysrmi = y_s - boxy
            IF (ysrmi < 1) ysrmi = 1
            ysrma = y_s + boxy
            IF (ysrma > nby) ysrma = nby
            ! beginning of the loop over the column (longitude) reflecting.
            DO x_sr = xsrmi, xsrma
              rx_sr = REAL(x_sr) * dx
              ! beginning of the loop over the rows (latitu) reflecting.
              DO y_sr = ysrmi, ysrma
                ry_sr = REAL(y_sr) * dy
                irefl = 0.0
                z_sr = altsol(x_sr, y_sr)
                IF ((x_sr > nbx) .OR. (x_sr < 1) .OR. &
                    (y_sr > nby) .OR. (y_sr < 1)) THEN
                  IF (verbose == 2) THEN
                    PRINT *, 'Ground cell out of borders'
                  END IF
                ELSE
                  IF ((x_s == x_sr) .AND. (y_s == y_sr) &
                      .AND. (z_s == z_sr)) THEN
                    IF (verbose == 2) THEN
                      PRINT *, 'Source pos = Ground cell'
                    END IF
                  ELSE
                    ! if haut is negative, the ground cell is lighted from below
                    haut = -(rx_s - rx_sr) * TAN(inclix(x_sr, y_sr)) - (ry_s - ry_sr) &
                           * TAN(incliy(x_sr, y_sr)) + z_s - z_sr
                    ! condition: the ground cell is lighted from above
                    IF (haut > 0.0) THEN
                      ! computation of the zenithal angle between the source and the surface reflectance
                      ! computation of the zenithal angle between the source and the line of sight voxel.
                      ! end of the case "observer at the same latitu/longitude than the source".

                      CALL anglezenithal(rx_s, ry_s, z_s, rx_sr, ry_sr, z_sr, angzen)
                      ! computation of the transmittance between the source and the ground surface
                      distd = SQRT((rx_s - rx_sr)**2.0 + (ry_s - ry_sr)**2.0 + (z_s - z_sr)**2.0)
                      CALL transmitm(angzen, z_s, z_sr, distd, transm, tranam, tabs)
                      CALL transmita(angzen, z_s, z_sr, distd, haer, transa, tranaa)
                      CALL transmitl(angzen, z_s, z_sr, distd, hlay, transl, tranal)

                      ! computation of the solid angle of the reflecting cell seen from the source
                      xc = DBLE(x_sr) * DBLE(dx)   ! Position in meters of the observer voxel (longitude).
                      yc = DBLE(y_sr) * DBLE(dy)   ! Position in meters of the observer voxel (latitu).
                      zc = DBLE(z_sr)              ! Position in meters of the observer voxel (altitude).
                      xn = DBLE(x_s) * DBLE(dx)    ! Position in meters of the source (longitude).
                      yn = DBLE(y_s) * DBLE(dy)    ! Position in meters of the source (latitu).
                      zn = DBLE(z_s)               ! Position in meters of the source (altitude).
                      epsilx = inclix(x_sr, y_sr)  ! tilt along x of the ground reflectance
                      epsily = incliy(x_sr, y_sr)  ! tilt along x of the ground reflectance

                      ! use a sub-grid surface when the reflectance radius is smaller than the cell size
                      dxp = MIN(dxp, reflsiz)
                      IF ((x_sr /= x_s) .OR. (y_sr /= y_s)) dxp = dx
                      dyp = MIN(dyp, reflsiz)
                      IF ((x_sr /= x_s) .OR. (y_sr /= y_s)) dyp = dy

                      ! computation of the composante along x of the first vector.
                      r1x = xc - DBLE(dxp) / 2.0 - xn
                      ! computation of the composante along y of the first vector.
                      r1y = yc + DBLE(dyp) / 2.0 - yn
                      ! computation of the composante en z of the first vector.
                      r1z = zc - TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                            + TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn

                      ! computation of the composante along x of the second vector.
                      r2x = xc + DBLE(dxp) / 2.0 - xn
                      ! computation of the composante along y of the second vector.
                      r2y = yc + DBLE(dyp) / 2.0 - yn
                      ! computation of the composante en z of the second vector.
                      r2z = zc + TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                            + TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn

                      ! computation of the composante along x of the third vector.
                      r3x = xc - DBLE(dxp) / 2.0 - xn
                      ! computation of the composante along y of the third vector.
                      r3y = yc - DBLE(dyp) / 2.0 - yn
                      ! computation of the composante en z of the third vector.
                      r3z = zc - TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                            - TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn

                      ! computation of the composante along x of the fourth vector.
                      r4x = xc + DBLE(dxp) / 2.0 - xn
                      ! computation of the composante along y of the fourth vector.
                      r4y = yc - DBLE(dyp) / 2.0 - yn
                      ! computation of the composante en z of the fourth vector.
                      r4z = zc + TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                            - TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn

                      ! Call of the routine anglesolide to compute the angle solide.
                      CALL anglesolide(omega, r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z)

                      IF (omega < 0.0) THEN
                        PRINT *, 'ERROR: Solid angle of the reflecting surface < 0.0'
                        STOP
                      END IF

! estimation of the half of the underlying angle of the solid angle
! this angle servira a obtenir un meilleur isime (moyenne) of
! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                      ouvang = SQRT(omega / pi)
                      ouvang = ouvang * 180.0 / pi

! computation of the photometric function of the light fixture toward the reflection surface
!=======================================================================

                      anglez = NINT(180.0 * angzen / pi)
                      IF (anglez < 0) anglez = -anglez
                      IF (anglez > 180) anglez = 360 - anglez
                      ! Transform the angle in integer degree into the position in the array.
                      anglez = anglez + 1

                      ! average +- ouvang
                      naz = 0
                      nbang = 0.0
                      P_indir = 0.0
                      DO na = -NINT(ouvang), NINT(ouvang)
                        naz = anglez + na
                        IF (naz < 0) naz = -naz
                        IF (naz > 181) naz = 362 - naz   ! symetric function
                        IF (naz == 0) naz = 1
                        P_indir = P_indir + pvalno(naz, stype) &
                                  * ABS(SIN(pi * REAL(naz) / 180.0)) / 2.0
                        nbang = nbang + 1.0 * ABS(SIN(pi * REAL(naz) / 180.0)) / 2.0
                      END DO
                      P_indir = P_indir / nbang
                      ! computation of the flux reaching the reflecting surface
                      flrefl = lamplu(x_s, y_s, stype) * &
                               P_indir * omega * transm * transa * transl
                      ! computation of the reflected intensity leaving the ground surface
                      ! The factor 1/pi comes from the normalisation of the fonction
                      irefl1 = flrefl * srefl / pi

! *********************************************************************************************************
! calculation of the direct radiance from reflection falling on a surface perpendicular
! to the viewing angle Units of W/nm/m2/sr
! *********************************************************************************************************

                      dho = SQRT((rx_obs - rx_sr)**2.0 + (ry_obs - ry_sr)**2.0)
                      IF ((dho > 0.0) .AND. (z_s /= z_obs)) THEN
                        ! zenithal angle source-observer
                        CALL anglezenithal(rx_obs, ry_obs, z_obs, rx_sr, ry_sr, z_sr, dzen)
                        ! computation of the angle azimutal direct line of sight-source
                        CALL angleazimutal(rx_obs, ry_obs, rx_sr, ry_sr, angazi)

                        ! 45deg. it is unlikely to have a 1km high mountain less than 1
                        hh = 1.0
                        IF (angzen > (pi / 4.0)) THEN
                          CALL horizon(x_sr, y_sr, z_sr, dx, dy, altsol, angazi, zhoriz, dh)
                          hh = MERGE(0.0, 1.0, dh <= dho .AND. angzen - zhoriz < 0.00001)
                        END IF

                        ! sub-grid obstacles
                        ff = 0.0
                        IF (x_dif >= 1 .AND. x_dif <= nbx .AND. y_dif >= 1 .AND. y_dif <= nbx) THEN
                          dho = SQRT((rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0)
                          IF (dho <= drefle(x_dif, y_dif)) THEN
                            angmin = pi / 2.0 - ATAN2((obsH(x_dif, y_dif) + altsol(x_dif, y_dif) - z_dif), drefle(x_dif, y_dif))
                            IF (angzen >= angmin) THEN
                              ff = ofill(x_dif, y_dif)
                            END IF
                          END IF
                        END IF

                        ! zenithal angle surface-observer
                        CALL anglezenithal(rx_sr, ry_sr, z_sr, rx_obs, ry_obs, z_obs, dzen)
                        ff2 = 0.0
                        ! light path from reflecting surface larger than the mean free path -> subgrid obstacles
                        IF (dho > drefle(x_sr, y_sr)) THEN
                          angmin = pi / 2.0 - ATAN2((altsol(x_sr, y_sr) + obsH(x_sr, y_sr) - z_sr), &
                                                    drefle(x_sr, y_sr))
                          ! condition sub-grid obstacles direct.
                          IF (dzen >= angmin) THEN
                            ff2 = ofill(x_sr, y_sr)
                          END IF
                        END IF
                        ! end light path to the observer larger than mean free path

                        ! zenithal angle source-observer
                        CALL anglezenithal(rx_obs, ry_obs, z_obs, rx_sr, ry_sr, z_sr, dzen)

                        ! projection angle of line to the lamp and the viewing angle scattering angle.
                        CALL angle3points(rx_sr, ry_sr, z_sr, rx_obs, ry_obs, z_obs, rx, ry, rz, dang)
                        dang = pi - dang

                        ! computation of the flux direct reaching the line of sight voxel
                        IF ((COS(dang) > 0.0) .AND. (dang < pi / 2.0)) THEN
                          ! distance direct sight between source and observer
                          ddir_obs = SQRT((rx_obs - rx_sr)**2.0 + (ry_obs - ry_sr)**2.0 + (z_obs - z_sr)**2.0)
                          ! computation of the solid angle of the line of sight voxel seen from the source
                          omega = 1.0 * ABS(COS(dang)) / ddir_obs**2.0
                          CALL transmitm(dzen, z_obs, z_sr, ddir_obs, transm, tranam, tabs)
                          CALL transmita(dzen, z_obs, z_sr, ddir_obs, haer, transa, tranaa)
                          CALL transmitl(dzen, z_obs, z_sr, ddir_obs, hlay, transl, tranal)

                          ! check if the reflecting surface enter the field of view of the observer
                          IF (dang < dfov) THEN
                            rdirect = rdirect + irefl1 * omega * transa * transm * transl &
                                      * hh * (1.0 - ff) * (1.0 - ff2) / (pi * dfov**2.0)
                          END IF
                          irrdirect = irrdirect + irefl1 * omega * transa * transm * transl &
                                      * hh * (1.0 - ff) * (1.0 - ff2)
                        END IF

                      END IF
                    END IF
                  END IF
                END IF
              END DO
            END DO

          END IF
        END DO
      END DO
    END IF
  END DO

  ! End of calculation of the direct radiances
  ! =================================

  IF (fsswit /= 0) THEN
    ! =================================
    ! Calculation of the scattered radiances

    ! temporaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cloudtop = 100000.0

    IF ((z_obs >= cloudbase) .AND. (z_obs <= cloudtop)) THEN
      PRINT *, 'The observer is inside the cloud! Abort computing.', z_obs, cloudbase
      STOP
    END IF
1110 FORMAT(I4, 1X, I4, 1X, I4)

    fctcld = 0.0
    ! Initialisation of the value of flux received by the sensor
    ftocap = 0.0

    ! calculating the distance before the line of sight beeing blocked by topography
    CALL horizon(x_obs, y_obs, z_obs, dx, dy, altsol, angaz1, zhoriz, dhmax)
    rx_c = REAL(x_obs) * dx - ix * scal / 2.0
    ry_c = REAL(y_obs) * dx - iy * scal / 2.0
    z_c = z_obs - iz * scal / 2.0

    ! beginning of the loop over the line of sight voxels
    DO icible = 1, ncible
      rx_c = rx_c + ix * (scalo / 2.0 + scal / 2.0)
      ry_c = ry_c + iy * (scalo / 2.0 + scal / 2.0)
      dh0 = SQRT((rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2)

      ! the line of sight is not yet blocked by the topography
      IF ((dh0 <= dhmax) .OR. ((dh0 > dhmax) .AND. (angze1 - zhoriz < 0.00001))) THEN
        x_c = NINT(rx_c / dx)
        IF (x_c < 1) x_c = 1
        IF (x_c > width) x_c = width
        y_c = NINT(ry_c / dy)
        IF (y_c < 1) y_c = 1
        IF (y_c > width) y_c = width
        z_c = z_c + iz * (scalo / 2.0 + scal / 2.0)
        ! stop the calculation of the viewing line when the increment is lower than 1/stoplim
        IF (z_c > altsol(x_c, y_c)) THEN
          ! or when hitting a cloud or when z>40km (scattering probability =0 (given precision)
          IF ((fcapt >= ftocap / stoplim) .AND. (z_c < cloudbase) .AND. (z_c < 35000.0)) THEN
            fcapt = 0.0
            DO i = 1, nbx
              DO j = 1, nby
                FCA(i, j) = 0.0
              END DO
            END DO

! Calculate the solid angle of the line of sight voxel unit voxel
! (1 m^3) given the fixed FOV of the observer.
! For line of sight voxel near the observer
! we need to calculate the scattering on a part of the voxel. For far
! voxels we may be needed to increase the solid angle since the FOV can
! encompass more than the voxel size. This correction is done with the
! portio parameter calculated as the ratio of the solid angle of the
! observer FOV over the line of sight voxel solid angle as seen from the
! observer.
            distd = SQRT((rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0 + (z_c - z_obs)**2.0)
            ! computation of the Solid angle of the line of sight voxel seen from the observer
            omega = 1.0 / distd**2.0
            IF (omega > omemax) THEN
              omega = 0.0
              portio = 0.0
            ELSE
              portio = (omefov / omega)
            END IF

            ! Initialisation of the contribution of the line of sight at the sensor level
            itotci = 0.0
            ITC = 0.0

            ! Condition line of sight inside the modelling domain
            IF ((rx_c > REAL(nbx * dx)) .OR. (rx_c < dx) .OR. &
                (ry_c > (nby * dy)) .OR. (ry_c < dy)) THEN
            ELSE
              IF (verbose >= 1) THEN
                PRINT *, '================================================'
                PRINT *, ' progression along the line of sight :', icible
                PRINT *, ' horizontal dist. line of sight =', SQRT((rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0), ' m'
                PRINT *, ' vertical dist. line of sight =', ABS(z_c - z_obs), ' m'
                WRITE (2, *) '============================================='
                WRITE (2, *) ' progression along the line of sight :', icible
                WRITE (2, *) ' horizontal dist. line of sight =', SQRT((rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0), ' m'
                WRITE (2, *) ' vertical dist. line of sight =', ABS(z_c - z_obs), ' m'
              END IF

              dis_obs = SQRT((z_c - z_obs)**2.0 + (ry_c - ry_obs)**2.0 + (rx_c - rx_obs)**2.0)
              IF (dis_obs == 0.0) THEN
                PRINT *, 'ERROR problem with dis_obs', dis_obs
                PRINT *, rx_c, x_obs, y_c, y_obs, z_c, z_obs
                STOP
              END IF

              ometif = pi * (diamobj / 2.0)**2.0 / dis_obs**2.0
              ! beginning of the loop over the types of light sources
              DO stype = 1, ntype
                ! check if there are any flux in that source type otherwise skip this lamp
                IF (totlu(stype) /= 0.0) THEN
                  IF (verbose >= 1) THEN
                    PRINT *, ' turning on lamps', stype
                    WRITE (2, *) ' turning on lamps', stype
                  END IF
                  ! Initialisation of the contribution of a source types to
                  ! the intensity toward the sensor by a line of sight voxel.
                  itotty = 0.0
                  ITT(:, :, stype) = 0.0
                  ! beginning of the loop over the column (longitude the) of the domain.
                  DO x_s = imin(stype), imax(stype)
                    ! beginning of the loop over the rows (latitud) of the domain.
                    DO y_s = jmin(stype), jmax(stype)
                      intdir = 0.0
                      itotind = 0.0
                      itodif = 0.0
                      itotrd = 0.0
                      isourc = 0.0
                      rx_s = REAL(x_s) * dx
                      ry_s = REAL(y_s) * dy
                      ! if the luminosite of the case is null, the program ignore this case.
                      IF (lamplu(x_s, y_s, stype) /= 0.0) THEN
                        ! Definition of the position (metre) vertical of the source.
                        z_s = (altsol(x_s, y_s) + lampal(x_s, y_s))

! *********************************************************************************************************
! * computation of the scattered intensity toward the observer by a line of sight voxel from the source   *
! *********************************************************************************************************

                        ! Initialisation of the verification of the position of the source
                        dirck = 0

                        ! if the position of the source and the line of sight voxel are the same
                        IF ((rx_s == rx_c) .AND. (ry_s == ry_c) .AND. (z_s == z_c)) THEN
                          dirck = 1
                          IF (verbose >= 1) PRINT *, 'Source = line of sight'
                        END IF

                        ! the source is not at the line of sight voxel position
                        IF (dirck /= 1) THEN
                          ! computation of the zenithal angle between the source and the line of sight
                          ! computation of the horizon for the resolved shadows direct
                          ! horizon resolution is 1 degree
                          distd = SQRT((rx_c - rx_s)**2.0 + (ry_c - ry_s)**2.0 + (z_c - z_s)**2.0)
                          dho = SQRT((rx_c - rx_s)**2.0 + (ry_c - ry_s)**2.0)
                          ! computation of the zenithal angle between the source and the line of sight voxel.
                          CALL anglezenithal(rx_s, ry_s, z_s, rx_c, ry_c, z_c, angzen)
                          ! computation of the angle azimutal direct line of sight-source
                          CALL angleazimutal(rx_s, ry_s, rx_c, ry_c, angazi)

                          ! 45deg. it is unlikely to have a 1km high mountain less than 1
                          hh = 1.0
                          IF (angzen > (pi / 4.0)) THEN
                            CALL horizon(x_sr, y_sr, z_sr, dx, dy, altsol, angazi, zhoriz, dh)
                            hh = MERGE(0.0, 1.0, dh <= dho .AND. angzen - zhoriz < 0.00001)
                          END IF
                          ! sub-grid obstacles
                          ff = 0.0
                          ! light path to observer larger than the mean free path -> subgrid obstacles
                          IF (dho > drefle(x_s, y_s)) THEN
                            angmin = pi / 2.0 - ATAN2((altsol(x_s, y_s) + obsH(x_s, y_s) - z_s), &
                                                      drefle(x_s, y_s))
                            ! condition sub-grid obstacles direct.
                            IF (angzen >= angmin) THEN
                              ff = ofill(x_s, y_s)
                            END IF
                          END IF
                          ! computation of the transmittance between the source and the line of sight
                          CALL transmitm(angzen, z_s, z_c, distd, transm, tranam, tabs)
                          CALL transmita(angzen, z_s, z_c, distd, haer, transa, tranaa)
                          CALL transmitl(angzen, z_s, z_c, distd, hlay, transl, tranal)
                          ! computation of the solid angle of the line of sight voxel seen from the source
                          omega = 1.0 / distd**2.0
                          IF (omega > omemax) omega = 0.0
                          anglez = NINT(180.0 * angzen / pi) + 1
                          P_dir = pvalno(anglez, stype)
                          ! computation of the flux reaching the line of sight voxel
                          ! correction for obstacle filling factor
                          fldir = lamplu(x_s, y_s, stype) * P_dir * omega * transm * transa * transl * (1.0 - ff) * hh
                          ! computation of the scattering probability of the direct light
                          ! distance pour traverser la cellule unitaire parfaitement orientée
                          IF (omega /= 0.0) THEN
                            ! scattering angle.
                            CALL angle3points(rx_s, ry_s, z_s, rx_c, ry_c, z_c, rx_obs, ry_obs, z_obs, angdif)
                            ! scattering probability of the direct light.
                            ! ############################################ secdif et un etaient inverses
                            CALL diffusion(angdif, tranam, tranaa, tranal, un, secdif, secdil, &
                                           fdifan, fdifl, haer, hlay, pdifdi, z_c)
                          ELSE
                            pdifdi = 0.0
                          END IF

                          ! computation of the source contribution to the scattered intensity toward the sensor by a line of sight voxel
                          intdir = fldir * pdifdi
                          ! contribution of the cloud reflection of the light coming directly from the source
                          ! line of sight voxel = cloud
                          IF (cloudt /= 0) THEN
                            IF (cloudbase - z_c <= iz * scal) THEN
                              ! zenith angle from cloud to observer
                              CALL anglezenithal(rx_c, ry_c, z_c, rx_obs, ry_obs, z_obs, azcl1)
                              ! zenith angle from source to cloud
                              CALL anglezenithal(rx_c, ry_c, z_c, rx_s, ry_s, z_s, azcl2)
                              doc2 = (rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0 + (z_c - z_obs)**2.0
                              dsc2 = (rx_s - rx_c)**2.0 + (ry_s - ry_c)**2.0 + (z_s - z_c)**2.0
                              ! cloud intensity from direct illum
                              CALL cloudreflectance(angzen, cloudt, rcloud)
                              icloud = icloud + fldir / omega * rcloud * doc2 * omefov &
                                       * ABS(COS(azcl2) / COS(azcl1)) / dsc2 / pi
                            END IF
                          END IF
                        ELSE
                          intdir = 0.0
                          ! end of the case position source is not equal to the line of sight voxel position
                        END IF
! end of the computation of the scattered intensity

! **********************************************************************************************************************
! * computation of the scattered light toward the observer by a line of sight voxel lighted by the ground reflection    *
! **********************************************************************************************************************

                        ! etablissement of the conditions ands boucles
                        ! Initialisation of the reflected intensity of the source
                        itotind = 0.0
                        itotrd = 0.0
                        xsrmi = x_s - boxx
                        IF (xsrmi < 1) xsrmi = 1
                        xsrma = x_s + boxx
                        IF (xsrma > nbx) xsrma = nbx
                        ysrmi = y_s - boxy
                        IF (ysrmi < 1) ysrmi = 1
                        ysrma = y_s + boxy
                        IF (ysrma > nby) ysrma = nby
                        ! beginning of the loop over the column (longitude) reflecting.
                        DO x_sr = xsrmi, xsrma
                          rx_sr = REAL(x_sr) * dx
                          ! beginning of the loop over the rows (latitu) reflecting.
                          DO y_sr = ysrmi, ysrma
                            ry_sr = REAL(y_sr) * dy
                            irefl = 0.0
                            z_sr = altsol(x_sr, y_sr)
                            IF ((x_sr > nbx) .OR. (x_sr < 1) .OR. &
                                (y_sr > nby) .OR. (y_sr < 1)) THEN
                              IF (verbose == 2) THEN
                                PRINT *, 'Ground cell out of borders'
                              END IF
                            ELSE
                              IF ((x_s == x_sr) .AND. (y_s == y_sr) .AND. (z_s == z_sr)) THEN
                                IF (verbose == 2) PRINT *, 'Source pos = Ground cell'
                              ELSE
                                ! if haut is negative, the ground cell is lighted from below
                                haut = -(rx_s - rx_sr) * TAN(inclix(x_sr, y_sr)) - (ry_s - ry_sr) &
                                       * TAN(incliy(x_sr, y_sr)) + z_s - z_sr

                                ! condition: the ground cell is lighted from above
                                IF (haut > 0.0) THEN
                                  ! computation of the zenithal angle between the source and the surface reflectance
                                  ! computation of the zenithal angle between the source and the line of sight voxel.
                                  ! end of the case "observer at the same latitu/longitude than the source".
                                  CALL anglezenithal(rx_s, ry_s, z_s, rx_sr, ry_sr, z_sr, angzen)
                                  ! computation of the transmittance between the source and the ground surface
                                  distd = SQRT((rx_s - rx_sr)**2.0 + (ry_s - ry_sr)**2.0 + (z_s - z_sr)**2.0)
                                  CALL transmitm(angzen, z_s, z_sr, distd, transm, tranam, tabs)
                                  CALL transmita(angzen, z_s, z_sr, distd, haer, transa, tranaa)
                                  CALL transmitl(angzen, z_s, z_sr, distd, hlay, transl, tranal)

                                  ! computation of the solid angle of the reflecting cell seen from the source
                                  xc = DBLE(x_sr) * DBLE(dx)   ! Position in meters of the observer voxel (longitude).
                                  yc = DBLE(y_sr) * DBLE(dy)   ! Position in meters of the observer voxel (latitu).
                                  zc = DBLE(z_sr)              ! Position in meters of the observer voxel (altitude).
                                  xn = DBLE(x_s) * DBLE(dx)    ! Position in meters of the source (longitude).
                                  yn = DBLE(y_s) * DBLE(dy)    ! Position in meters of the source (latitu).
                                  zn = DBLE(z_s)               ! Position in meters of the source (altitude).
                                  epsilx = inclix(x_sr, y_sr)  ! tilt along x of the ground reflectance
                                  epsily = incliy(x_sr, y_sr)  ! tilt along x of the ground reflectance

                                  ! use a sub-grid surface when the reflectance radius is smaller than the cell size
                                  dxp = MIN(dxp, reflsiz)
                                  IF ((x_sr /= x_s) .OR. (y_sr /= y_s)) dxp = dx
                                  dyp = MIN(dyp, reflsiz)
                                  IF ((x_sr /= x_s) .OR. (y_sr /= y_s)) dyp = dy

                                  ! computation of the composante along x of the first vector.
                                  r1x = xc - DBLE(dxp) / 2.0 - xn
                                  ! computation of the composante along y of the first vector.
                                  r1y = yc + DBLE(dyp) / 2.0 - yn
                                  ! computation of the composante en z of the first vector.
                                  r1z = zc - TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                                        + TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn
                                  ! computation of the composante along x of the second vector.
                                  r2x = xc + DBLE(dxp) / 2.0 - xn
                                  ! computation of the composante along y of the second vector.
                                  r2y = yc + DBLE(dyp) / 2.0 - yn
                                  ! computation of the composante en z of the second vector.
                                  r2z = zc + TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                                        + TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn
                                  ! computation of the composante along x of the third vector.
                                  r3x = xc - DBLE(dxp) / 2.0 - xn
                                  ! computation of the composante along y of the third vector.
                                  r3y = yc - DBLE(dyp) / 2.0 - yn
                                  ! computation of the composante en z of the third vector.
                                  r3z = zc - TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                                        - TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn
                                  ! computation of the composante along x of the fourth vector.
                                  r4x = xc + DBLE(dxp) / 2.0 - xn
                                  ! computation of the composante along y of the fourth vector.
                                  r4y = yc - DBLE(dyp) / 2.0 - yn
                                  ! computation of the composante en z of the fourth vector.
                                  r4z = zc + TAN(DBLE(epsilx)) * DBLE(dxp) / 2.0 &
                                        - TAN(DBLE(epsily)) * DBLE(dyp) / 2.0 - zn

                                  ! Call of the routine anglesolide to compute the angle solide.
                                  CALL anglesolide(omega, r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z)
                                  IF (omega < 0.0) THEN
                                    PRINT *, 'ERROR: Solid angle of the reflecting surface < 0.0'
                                    STOP
                                  END IF

                                  ! estimation of the half of the underlying angle of the solid angle
                                  ! this angle servira a obtenir un meilleur isime (moyenne) of
                                  ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                                  ouvang = SQRT(omega / pi)
                                  ouvang = ouvang * 180.0 / pi

! computation of the photometric function of the light fixture toward the reflection surface
!=======================================================================

                                  anglez = NINT(180.0 * angzen / pi)
                                  IF (anglez < 0) &
                                    anglez = -anglez
                                  IF (anglez > 180) anglez = 360 - anglez
                                  ! Transform the angle in integer degree into the position in the array.
                                  anglez = anglez + 1

                                  ! average +- ouvang
                                  naz = 0
                                  nbang = 0.0
                                  P_indir = 0.0
                                  DO na = -NINT(ouvang), NINT(ouvang)
                                    naz = anglez + na
                                    IF (naz < 0) naz = -naz
                                    IF (naz > 181) naz = 362 - naz   ! symetric function
                                    IF (naz == 0) naz = 1
                                    P_indir = P_indir + pvalno(naz, stype) &
                                              * ABS(SIN(pi * REAL(naz) / 180.0)) / 2.0
                                    nbang = nbang + 1.0 * ABS(SIN(pi * REAL(naz) / 180.0)) / 2.0
                                  END DO

                                  P_indir = P_indir / nbang
                                  ! computation of the flux reaching the reflecting surface
                                  flrefl = lamplu(x_s, y_s, stype) * &
                                           P_indir * omega * transm * transa * &
                                           transl
                                  ! computation of the reflected intensity leaving the ground surface
                                  ! The factor 1/pi comes from the normalisation of the fonction
                                  irefl1 = flrefl * srefl / pi

! **************************************************************************************
! * computation of the 2nd scattering contributions (2 order scattering and after reflection)
! **************************************************************************************

                                  IF (effdif > 0.0) THEN
                                    nss = 0
                                    ndi = 0
                                    ! beginning of the loop over the scattering voxels.
                                    DO idi = 1, ndiff
                                      rx_dif = zondif(idi, 1) + (rx_s + rx_c) / 2.0
                                      x_dif = NINT(rx_dif / dx)
                                      ry_dif = zondif(idi, 2) + (ry_s + ry_c) / 2.0
                                      y_dif = NINT(ry_dif / dy)
                                      z_dif = zondif(idi, 3) + (z_s + z_c) / 2.0
                                      id = NINT(rx_dif / dx)
                                      IF (id > width) id = width
                                      IF (id < 1) id = 1
                                      jd = NINT(ry_dif / dy)
                                      IF (jd > width) jd = width
                                      IF (jd < 1) jd = 1
                                      ! beginning diffusing cell underground
                                      IF (z_dif - siz / 2.0 <= altsol(id, jd) .OR. &
                                          (z_dif > 35000.0) .OR. (z_dif > cloudbase)) THEN
                                        ndi = ndi + 1
                                      ELSE
                                        ds1 = SQRT((rx_sr - rx_dif)**2.0 + (ry_sr - ry_dif)**2.0 + (z_sr - z_dif)**2.0)
                                        ds2 = SQRT((rx_c - rx_dif)**2.0 + (ry_c - ry_dif)**2.0 + (z_c - z_dif)**2.0)
                                        ds3 = SQRT((rx_s - rx_dif)**2.0 + (ry_s - ry_dif)**2.0 + (z_s - z_dif)**2.0)

                                        IF ((ds1 < dss) .OR. (ds2 < dss) .OR. (ds3 < dss)) THEN
                                          nss = nss + 1
                                        ELSE
                                          dho = SQRT((rx_dif - rx_sr)**2.0 + (ry_dif - ry_sr)**2.0)
                                          ! computation of the zenithal angle between the reflection surface and the scattering voxel
                                          ! shadow reflection surface-scattering voxel
                                          ! computation of the zenithal angle reflection surface - scattering voxel.
                                          CALL anglezenithal(rx_sr, ry_sr, z_sr, rx_dif, ry_dif, z_dif, angzen)
                                          ! computation of the angle azimutal line of sight-scattering voxel
                                          CALL angleazimutal(rx_sr, ry_sr, rx_dif, ry_dif, angazi)

                                          ! horizon blocking not a matte because dif are closeby and some downward
                                          hh = 1.0
                                          ! sub-grid obstacles
                                          ff = 0.0
                                          ! light path to observer larger than the mean free path -> subgrid obstacles
                                          IF (dho > drefle(x_sr, y_sr)) THEN
                                            angmin = pi / 2.0 - ATAN2(obsH(x_sr, y_sr), drefle(x_sr, y_sr))
                                            ! condition obstacle reflechi->scattered
                                            IF (angzen >= angmin) THEN
                                              ff = ofill(x_sr, y_sr)
                                            END IF
                                          END IF
                                          ! computation of the transmittance between the reflection surface and the scattering voxel
                                          distd = SQRT((rx_dif - rx_sr)**2.0 + (ry_dif - ry_sr)**2.0 + (z_dif - z_sr)**2.0)
                                          CALL transmitm(angzen, z_sr, z_dif, distd, transm, tranam, tabs)
                                          CALL transmita(angzen, z_sr, z_dif, distd, haer, transa, tranaa)
                                          CALL transmitl(angzen, z_sr, z_dif, distd, hlay, transl, tranal)
                                          ! computation of the solid angle of the scattering voxel seen from the reflecting surface
                                          omega = 1.0 / distd**2.0
                                          IF (omega > omemax) omega = 0.0
                                          ! computing flux reaching the scattering voxel
                                          fldif2 = irefl1 * omega * transm * transa * transl * (1.0 - ff) * hh
                                          ! computing the scattering probability toward the line of sight voxel
                                          ! cell unitaire
                                          IF (omega /= 0.0) THEN
                                            ! scattering angle.
                                            CALL angle3points(rx_sr, ry_sr, z_sr, rx_dif, ry_dif, z_dif, &
                                                              rx_c, ry_c, z_c, angdif)
                                            ! scattering probability of the direct light.
                                            CALL diffusion(angdif, tranam, tranaa, tranal, un, secdif, secdil, &
                                                           fdifan, fdifl, haer, hlay, pdifd1, z_dif)
                                          ELSE
                                            pdifd1 = 0.0
                                          END IF

                                          volu = siz**3.0
                                          IF (volu < 0.0) THEN
                                            PRINT *, 'ERROR, volume 2 is negative!'
                                            STOP
                                          END IF

                                          ! computing scattered intensity toward the line of sight voxel from the scattering voxel
                                          idif2 = fldif2 * pdifd1 * volu
                                          ! computing zenith angle between the scattering voxel and the line of sight voxel
                                          ! computation of the zenithal angle between the scattering voxel and the line of sight voxel.
                                          CALL anglezenithal(rx_dif, ry_dif, z_dif, rx_c, ry_c, z_c, angzen)
                                          ! computation of the azimutal angle surf refl-scattering voxel
                                          CALL angleazimutal(rx_dif, ry_dif, rx_c, ry_c, angazi)
                                          ! subgrid obstacles
                                          ff = 0.0
                                          IF (x_dif >= 1 .AND. x_dif <= nbx .AND. y_dif >= 1 .AND. y_dif <= nbx) THEN
                                            dho = SQRT((rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0)
                                            IF (dho <= drefle(x_dif, y_dif)) THEN
                                angmin = pi / 2.0 - ATAN2((obsH(x_dif, y_dif) + altsol(x_dif, y_dif) - z_dif), drefle(x_dif, y_dif))
                                              IF (angzen >= angmin) THEN
                                                ff = ofill(x_dif, y_dif)
                                              END IF
                                            END IF
                                          END IF

                                          hh = 1.0
                                          ! computing transmittance between the scattering voxel and the line of sight voxel
                                          distd = SQRT((rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0 + (z_dif - z_c)**2.0)
                                          CALL transmitm(angzen, z_dif, z_c, distd, transm, tranam, tabs)
                                          CALL transmita(angzen, z_dif, z_c, distd, haer, transa, tranaa)
                                          CALL transmitl(angzen, z_dif, z_c, distd, hlay, transl, tranal)
                                          ! computing the solid angle of the line of sight voxel as seen from the scattering voxel
                                          omega = 1.0 / distd**2.0
                                          IF (omega > omemax) omega = 0.0
                                          ! computation of the scattered flux reaching the line of sight voxel
                                          fdif2 = idif2 * omega * transm * transa * transl * (1.0 - ff) * hh
                                          ! cloud contribution for double scat from a reflecting pixel
                                          ! line of sight voxel = cloud
                                          IF (cloudt /= 0) THEN
                                            IF (cloudbase - z_c <= iz * scal) THEN
                                              ! zenith angle from cloud to observer
                                              CALL anglezenithal(rx_c, ry_c, z_c, rx_obs, ry_obs, z_obs, azcl1)
                                              ! zenith angle from source to cloud
                                              CALL anglezenithal(rx_c, ry_c, z_c, rx_dif, ry_dif, z_dif, azcl2)
                                              doc2 = (rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0 + (z_c - z_obs)**2.0
                                              dsc2 = (rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0 + (z_dif - z_c)**2.0
                                              ! cloud intensity from direct illum
                                              CALL cloudreflectance(angzen, cloudt, rcloud)
                                              icloud = icloud + &
                                                       fldif2 / omega * rcloud * doc2 * omefov * &
                                                       ABS(COS(azcl2) / COS(azcl1)) / dsc2 / pi
                                            END IF
                                          END IF
                                          ! computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c)
                                          IF (omega /= 0.0) THEN
                                            ! scattering angle.
                                            CALL angle3points(rx_dif, ry_dif, z_dif, rx_c, ry_c, z_c, &
                                                              rx_obs, ry_obs, z_obs, angdif)
                                            ! scattering probability of the direct light.
                                            CALL diffusion(angdif, tranam, tranaa, tranal, un, secdif, &
                                                           secdil, fdifan, fdifl, haer, hlay, pdifd2, z_c)
                                          ELSE
                                            pdifd2 = 0.0
                                          END IF

                                          ! computing scattered intensity toward the observer from the line of sight voxel
                                          idif2p = fdif2 * pdifd2
                                          ! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation
                                          idif2p = idif2p * REAL(stepdi) * REAL(ndiff - ndi) &
                                                   / REAL(ndiff - ndi - nss)
                                          itotrd = itotrd + idif2p

! ********************************************************************************
! *  section for the calculation of the 2nd scat from the source without reflexion
! ********************************************************************************

                                          ! beginning condition source = reflection for the computation of the source scat line of sight
                                          IF ((x_sr == x_s) .AND. (y_sr == y_s)) THEN
                                            ! computation of the zenithal angle between the source and the scattering voxel
                                            ! shadow source-scattering voxel
                                            ! computation of the zenithal angle source-scattering voxel.
                                            CALL anglezenithal(rx_s, ry_s, z_s, rx_dif, ry_dif, z_dif, angzen)
                                            ! computation of the angle azimutal line of sight-scattering voxel
                                            CALL angleazimutal(rx_s, ry_s, rx_dif, ry_dif, angazi)
                                            ! horizon blocking not a matter because some path are downward and most of them closeby
                                            hh = 1.0
                                            angmin = pi / 2.0 - ATAN2((obsH(x_s, y_s) + altsol(x_s, y_s) - z_s), &
                                                                      drefle(x_s, y_s))
                                            ! condition obstacle source->scattering.
                                            IF (angzen < angmin) THEN
                                              ff = 0.0
                                            ELSE
                                              ff = ofill(x_s, y_s)
                                            END IF
                                            ! computation of the transmittance between the source and the scattering voxel
                                            distd = SQRT((rx_s - rx_dif)**2.0 + (ry_s - ry_dif)**2.0 + (z_s - z_dif)**2.0)
                                            CALL transmitm(angzen, z_s, z_dif, distd, transm, tranam, tabs)
                                            CALL transmita(angzen, z_s, z_dif, distd, haer, transa, tranaa)
                                            CALL transmitl(angzen, z_s, z_dif, distd, hlay, transl, tranal)
                                            ! computation of the Solid angle of the scattering unit voxel seen from the source
                                            omega = 1.0 / distd**2.0
                                            IF (omega > omemax) omega = 0.0
                                            anglez = NINT(180.0 * angzen / pi) + 1
                                            P_dif1 = pvalno(anglez, stype)
                                            ! computing flux reaching the scattering voxel
                                            fldif1 = lamplu(x_s, y_s, stype) * P_dif1 * &
                                                     omega * transm * transa * transl * (1.0 - ff) * hh
                                            ! computing the scattering probability toward the line of sight voxel
                                            IF (omega /= 0.0) THEN
                                              ! scattering angle.
                                              CALL angle3points(rx_s, ry_s, z_s, rx_dif, ry_dif, z_dif, &
                                                                rx_c, ry_c, z_c, angdif)
                                              ! scattering probability of the direct light.
                                              CALL diffusion(angdif, &
                                                             tranam, tranaa, tranal, un, secdif, secdil, &
                                                             fdifan, fdifl, haer, hlay, pdifd1, z_dif)
                                            ELSE
                                              pdifd1 = 0.0
                                            END IF

                                            volu = siz**3.0
                                            IF (volu < 0.0) THEN
                                              PRINT *, 'ERROR, volume 1 is negative!'
                                              STOP
                                            END IF

                                            ! computing scattered intensity toward the line of sight voxel from the scattering voxel
                                            idif1 = fldif1 * pdifd1 * volu
                                            ! computing zenith angle between the scattering voxel and the line of sight voxel
                                            CALL anglezenithal(rx_dif, ry_dif, z_dif, rx_c, ry_c, z_c, angzen)
                                            ! computation of the azimutal angle surf refl-scattering voxel
                                            CALL angleazimutal(rx_dif, ry_dif, rx_c, ry_c, angazi)
                                            ! subgrid obstacles
                                            IF ((x_dif < 1) .OR. (x_dif > nbx) .OR. (y_dif < 1) .OR. (y_dif > nbx)) THEN
                                              dho = SQRT((rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0)
                                              ff = 0.0
                                            ELSE
                                              ff = 0.0
                                              IF (dho > drefle(x_dif, y_dif)) THEN
                                                angmin = pi / 2.0 - ATAN2((obsH(x_dif, y_dif) + altsol(x_dif, y_dif) - z_dif), &
                                                                          drefle(x_dif, y_dif))
                                                ! condition obstacles scattering->line of sight
                                                IF (angzen >= angmin) THEN
                                                  ff = ofill(x_dif, y_dif)
                                                END IF
                                              END IF
                                            END IF

                                            hh = 1.0
                                            ! Computing transmittance between the scattering voxel and the line of sight voxel
                                            distd = SQRT((rx_c - rx_dif)**2.0 + (ry_c - ry_dif)**2.0 + (z_c - z_dif)**2.0)
                                            CALL transmitm(angzen, z_dif, z_c, distd, transm, tranam, tabs)
                                            CALL transmita(angzen, z_dif, z_c, distd, haer, transa, tranaa)
                                            CALL transmitl(angzen, z_dif, z_c, distd, hlay, transl, tranal)

                                            ! computing the solid angle of the line of sight voxel as seen from the scattering voxel
                                            omega = 1.0 / distd**2.0
                                            IF (omega > omemax) omega = 0.0
                                            ! computation of the scattered flux reaching the line of sight voxel
                                            fldiff = idif1 * omega * transm * transa * transl * (1.0 - ff) * hh
                                            ! cloud contribution to the double scattering from a source
                                            ! line of sight voxel = cloud
                                            IF (cloudt /= 0) THEN
                                              IF (cloudbase - z_c <= iz * scal) THEN
                                                ! zenith angle from cloud to observer
                                                CALL anglezenithal(rx_c, ry_c, z_c, rx_obs, ry_obs, z_obs, azcl1)
                                                ! zenith angle from source to cloud
                                                CALL anglezenithal(rx_c, ry_c, z_c, rx_dif, ry_dif, z_dif, azcl2)
                                                doc2 = (rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0 + (z_c - z_obs)**2.0
                                                dsc2 = (rx_dif - rx_c)**2.0 + (ry_dif - ry_c)**2.0 + (z_dif - z_c)**2.0
                                                ! cloud intensity from direct illum
                                                CALL cloudreflectance(angzen, cloudt, rcloud)
                                                icloud = icloud + &
                                                         fldiff / omega * rcloud * doc2 * omefov * &
                                                         ABS(COS(azcl2) / COS(azcl1)) / dsc2 / pi
                                              END IF
                                            END IF

                                            ! computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c)
                                            IF (omega /= 0.0) THEN
                                              ! scattering angle.
                                              CALL angle3points(rx_dif, ry_dif, z_dif, rx_c, ry_c, z_c, &
                                                                rx_obs, ry_obs, z_obs, angdif)
                                              ! scattering probability of the direct light.
                                              CALL diffusion(angdif, &
                                                             tranam, tranaa, tranal, un, secdif, secdil, &
                                                             fdifan, fdifl, haer, hlay, pdifd2, z_c)
                                            ELSE
                                              pdifd2 = 0.0
                                            END IF

                                            ! computing scattered intensity toward the observer from the line of sight voxel
                                            idiff2 = fldiff * pdifd2
                                            ! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation
                                            idiff2 = idiff2 * REAL(stepdi) * REAL(ndiff - ndi) / &
                                                     REAL(ndiff - ndi - nss)
                                            ! sum over the scattering voxels
                                            itodif = itodif + idiff2
                                            ! end condition source = reflection for the computation of the source scat line of sight
                                          END IF
                                          ! end of the case scattering pos = source pos or line of sight pos
                                        END IF
                                        ! end diffusing celle underground
                                      END IF
                                      ! end of the loop over the scattering voxels.
                                    END DO
                                    ! end of the condition ou effdif > 0
                                  END IF
! End of 2nd scattered intensity calculations
!===================================================================

! **********************************************************************
! * section refected light with single scattering
! **********************************************************************

                                  ! verify if there is shadow between sr and line of sight voxel
                                  ! zenithal angle between the reflecting surface and the line of sight voxel.
                                  CALL anglezenithal(rx_sr, ry_sr, z_sr, rx_c, ry_c, z_c, angzen)
                                  ! computation of the azimutal angle reflect-line of sight
                                  CALL angleazimutal(rx_sr, ry_sr, rx_c, ry_c, angazi)
                                  distd = SQRT((rx_sr - rx_c)**2.0 + (ry_sr - ry_c)**2.0 + (z_sr - z_c)**2.0)
                                  dho = SQRT((rx_sr - rx_c)**2.0 + (ry_sr - ry_c)**2.0)

                                  ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                  hh = 1
                                  IF (angzen > pi / 4.0) THEN
                                    CALL horizon(x_sr, y_sr, z_sr, dx, dy, altsol, angazi, zhoriz, dh)
                                    ! the path line of sight-reflec is not below the horizon => we compute
                                    hh = MERGE(1, 0, (dh <= dho) .AND. (angzen - zhoriz < 0.00001))
                                  END IF

                                  irefl = irefl1
                                  ! case: line of sight position = Position of reflecting cell
                                  IF ((rx_c == rx_sr) .AND. (ry_c == ry_sr) .AND. (z_c == z_sr)) THEN
                                    intind = 0.0
                                  ELSE
                                    ! obstacle
                                    dho = SQRT((rx_sr - rx_c)**2.0 + (ry_sr - ry_c)**2.0)
                                    ff = 0.0
                                    IF (dho > drefle(x_sr, y_sr)) THEN
                                      angmin = pi / 2.0 - ATAN2(obsH(x_sr, y_sr), drefle(x_sr, y_sr))
                                      ! condition obstacle reflected.
                                      IF (angzen >= angmin) THEN
                                        ff = ofill(x_sr, y_sr)
                                      END IF
                                    END IF

                                    ! computation of the transmittance between the ground surface and the line of sight voxel
                                    CALL transmitm(angzen, z_sr, z_c, distd, transm, tranam, tabs)
                                    CALL transmita(angzen, z_sr, z_c, distd, haer, transa, tranaa)
                                    CALL transmitl(angzen, z_sr, z_c, distd, hlay, transl, tranal)
                                    ! computation of the solid angle of the line of sight voxel seen from the reflecting cell
                                    omega = 1.0 / distd**2.0
                                    IF (omega > omemax) omega = 0.0
                                    ! computation of the flux reflected reaching the line of sight voxel
                                    ! obstacles correction
                                    flindi = irefl * omega * transm * transa * transl * (1.0 - ff) * hh
                                    ! cloud contribution to the reflected light from a ground pixel
                                    ! line of sight voxel = cloud
                                    IF (cloudt /= 0) THEN
                                      IF (cloudbase - z_c <= iz * scal) THEN
                                        ! zenith angle from cloud to observer
                                        CALL anglezenithal(rx_c, ry_c, z_c, rx_obs, ry_obs, z_obs, azcl1)
                                        ! zenith angle from source to cloud
                                        CALL anglezenithal(rx_c, ry_c, z_c, rx_sr, ry_sr, z_sr, azcl2)
                                        ! cloud intensity from direct illum
                                        doc2 = (rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0 + (z_c - z_obs)**2.0
                                        dsc2 = (rx_sr - rx_c)**2.0 + (ry_sr - ry_c)**2.0 + (z_sr - z_c)**2.0
                                        CALL cloudreflectance(angzen, cloudt, rcloud)
                                        icloud = icloud + &
                                                 flindi / omega * rcloud * doc2 * omefov * &
                                                 ABS(COS(azcl2) / COS(azcl1)) / dsc2 / pi
                                      END IF
                                    END IF
                                    ! computation of the scattering probability of the reflected light
                                    IF (omega /= 0.0) THEN
                                      ! scattering angle.
                                      CALL angle3points(rx_sr, &
                                                        ry_sr, z_sr, rx_c, ry_c, z_c, &
                                                        rx_obs, ry_obs, z_obs, angdif)
                                      ! scattering probability of the reflected light.
                                      CALL diffusion(angdif, &
                                                     tranam, tranaa, tranal, un, &
                                                     secdif, secdil, fdifan, fdifl, &
                                                     haer, hlay, pdifin, z_c)
                                    ELSE
                                      pdifin = 0.0
                                    END IF
                                    ! computation of the reflected intensity toward the sensor by a reflecting cell
                                    intind = flindi * pdifin * (1.0 - ff) * hh
                                    ! end of the case posi reflecting cell =  line of sight voxel position
                                  END IF
                                  ! Sum of the intensities of each reflecting cell.
                                  itotind = itotind + intind
                                  ! end of the condition surface not lighted from the top.
                                END IF
                                ! end of the condition reflecting cell is not on the source.
                              END IF
                              ! end of the condition surface of the domain.
                            END IF
                            ! end of the loop over the rows (latitu) reflecting.
                          END DO
                          ! end of the loop over the column (longitude) reflecting.
                        END DO
                        !   end of the computation of the reflected intensity

!**********************************************************************
! computation of the total intensity coming from a source to the line of sight voxel toward the sensor
!**********************************************************************

                        ! In the order 1st scat; refl->1st scat; 1st scat->2nd scat,
                        ! refl->1st scat->2nd scat
                        ! Sum of the intensities of a given type of source reaching the line of sight voxel.
                        isourc = intdir + itotind + itodif + itotrd
                        ! scaling the values according to the path length in the l. of sight voxel of 1m3
                        isourc = isourc * scal
                        ! correct for the field of view of the observer
                        isourc = isourc * portio
                        ! include clouds in the total intensity
                        ! TODO: PROBABLY AN ERROR, VALIDATE
                        ! isourc = isourc + icloud

                        IF ((itodif < 0.0) .OR. (itotrd < 0.0)) THEN
                          PRINT *, intdir, itotind, itodif, itotrd
                          STOP
                        END IF

                        IF (verbose == 2) THEN
                          PRINT *, ' Total intensity per component for type ', ntype, ':'
                          PRINT *, ' source->scattering = ', intdir
                          PRINT *, ' source->reflexion->scattering = ', itotind
                          PRINT *, ' source->scattering->scattering = ', itodif
                          PRINT *, ' source->reflexion->scattering->scattering = ', itotrd
                          IF (intdir * itotind * itodif * itotrd < 0.0) THEN
                            PRINT *, 'PROBLEM! Negative intensity.'
                            STOP
                          END IF
                        END IF

!**********************************************************************
! computation of the total intensity coming from all the sources of a given type
!**********************************************************************

                        ! Sum of the intensities all sources of the same typeand a given line of sight element
                        itotty = itotty + isourc
                        ! ITT stores itotty in a matrix
                        ITT(x_s, y_s, stype) = ITT(x_s, y_s, stype) + isourc
                        ! end of the condition "the luminosity of the ground pixel x_s,y_s in not null".
                      END IF
                      ! end the loop over the lines (latitude) of the domain (y_s).
                    END DO
                    ! end the loop over the column (longitude) of the domain (x_s).
                  END DO
                  ! end of the computation of the intensity of one source type
                  ! Sum of the intensities all source all type to a line of sight element
                  itotci = itotci + itotty
                  DO x_s = imin(stype), imax(stype)
                    DO y_s = jmin(stype), jmax(stype)
                      ITC(x_s, y_s) = ITC(x_s, y_s) + ITT(x_s, y_s, stype)
                    END DO
                  END DO

                  ! calculate total lamp flux matrix for all lamp types
                  DO x_s = 1, nbx
                    DO y_s = 1, nby
                      lpluto(x_s, y_s) = lpluto(x_s, y_s) + lamplu(x_s, y_s, stype)
                    END DO
                  END DO
                  ! end of condition if there are any flux in that source type
                END IF
                ! end of the loop over the types of sources (stype).
              END DO
              ! end of the computation of the intensity coming from a line of sight voxel toward the sensor

!***********************************************************************
! computation of the luminous flux reaching the observer
!***********************************************************************

              ! computation of the zenithal angle between the line of sight voxel and the observer.
              CALL anglezenithal(rx_c, ry_c, z_c, rx_obs, ry_obs, z_obs, angzen)
              ! end of the case "observer at the same latitu/longitude than the source".

              ! computation of the transmittance between the line of sight voxel and the observer
              distd = SQRT((rx_c - rx_obs)**2.0 + (ry_c - ry_obs)**2.0 + (z_c - z_obs)**2.0)
              CALL transmitm(angzen, z_c, z_obs, distd, transm, tranam, tabs)
              CALL transmita(angzen, z_c, z_obs, distd, haer, transa, tranaa)
              CALL transmitl(angzen, z_c, z_obs, distd, hlay, transl, tranal)

              ! computation of the flux reaching the intrument from the line of sight voxel
              fcapt = itotci * ometif * transa * transm * transl
              DO x_s = 1, nbx
                DO y_s = 1, nby
                  FCA(x_s, y_s) = ITC(x_s, y_s) * ometif * transa * transm * transl
                END DO
              END DO

              IF (COS(pi - angzen) == 0.0) THEN
                PRINT *, 'ERROR perfectly horizontal sight is forbidden'
                STOP
              END IF
              ! end of the computation of the flux reaching the observer voxel from the line of sight voxel

              ! flux for all source all type all line of sight element
              ftocap = ftocap + fcapt
              DO x_s = 1, nbx
                DO y_s = 1, nby
                  ! FTC is the array of the flux total at the sensor to identify
                  FTC(x_s, y_s) = FTC(x_s, y_s) + FCA(x_s, y_s)
                  ! the contribution of each ground pixel to the total flux at the observer level
                  ! The % is simply given by the ratio FTC/ftocap
                END DO
              END DO

              ! correction for the FOV to the flux reaching the intrument from the cloud voxel
              IF (cloudt /= 0) THEN
                ! computation of the flux reaching the intrument from the cloud voxel
                fccld = icloud * ometif * transa * transm * transl
                ! cloud flux for all source all type all line of sight element
                fctcld = fctcld + fccld
              END IF
              IF (verbose >= 1) THEN
                PRINT *, 'added radiance =', fcapt / omefov / (pi * (diamobj / 2.0)**2.0)
                PRINT *, 'radiance accumulated =', ftocap / omefov / (pi * (diamobj / 2.0)**2.0)
                WRITE (2, *) 'added radiance =', fcapt / omefov / (pi * (diamobj / 2.0)**2.0)
                WRITE (2, *) 'radiance accumulated =', ftocap / omefov / (pi * (diamobj / 2.0)**2.0)
              END IF
              ! end of the condition line of sight voxel inside the modelling domain
            END IF
            ! end condition line of sight voxel 1/stoplim
          END IF

          ! accelerate the computation as we get away from the sources
          scalo = scal
          IF (scal <= 3000.0) scal = scal * 1.12
        END IF
      ELSE
        ! print*,'End of line of sight - touching the ground'
        ! line of sight not blocked by topography
      END IF
      ! end of the loop over the line of sight voxels.
    END DO

    ! correction for the cloud fraction (defined from 0 to 100)
    fctcld = fctcld * 10**(0.4 * (100.0 - cloudfrac) * cloudslope)
    IF (prmaps == 1) THEN
      ! open(unit=9,file=pclf,status='unknown')
      DO x_s = 1, nbx
        DO y_s = 1, nby
          ! Here FTC becomes the flux fraction of each pixel. The sum of FTC values over all pixels give the total flux
          FTC(x_s, y_s) = FTC(x_s, y_s) / ftocap
        END DO
      END DO

      IF (verbose == 2) THEN
        PRINT *, 'Writing normalized contribution array'
        PRINT *, 'Warning Cloud contrib. excluded from that array.'
      END IF
!            do x_s=1,nbx
!              do y_s=1,nby
!                write(9,*) x_s,y_s,FTC(x_s,y_s)                           ! emettrice au sol, c'est un % par unite of watt installes
!              enddo
!            enddo
      CALL twodout(nbx, nby, pclimg, FTC)
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
      ! end of condition for creating contrib and sensit maps
    END IF
    ! end of scattered light
  END IF

! End of calculation of the scattered radiances
! =================================

  IF (verbose >= 1) PRINT *, '====================================================='
  PRINT *, '         Direct irradiance from sources (W/m**2/nm)'
  WRITE (*, 2001) irdirect
  PRINT *, '       Direct irradiance from reflexion (W/m**2/nm)'
  WRITE (*, 2001) irrdirect
  PRINT *, '         Direct radiance from sources (W/str/m**2/nm)'
  WRITE (*, 2001) direct
  PRINT *, '         Direct radiance from reflexion (W/str/m**2/nm)'
  WRITE (*, 2001) rdirect
  PRINT *, '             Cloud radiance (W/str/m**2/nm)'
  WRITE (*, 2001) fctcld / omefov / (pi * (diamobj / 2.0)**2.0)
  PRINT *, '            Diffuse radiance (W/str/m**2/nm)'
  WRITE (*, 2001) (ftocap + fctcld) / omefov / (pi * (diamobj / 2.0)**2.0)
  IF (verbose >= 1) WRITE (2, *) '====================================================='
  WRITE (2, *) '     Direct irradiance from sources (W/m**2/nm)'
  WRITE (2, 2001) irdirect
  WRITE (2, *) '     Direct irradiance from reflexion (W/m**2/nm)'
  WRITE (2, 2001) irrdirect
  WRITE (2, *) '     Direct radiance from sources (W/str/m**2/nm)'
  WRITE (2, 2001) direct
  WRITE (2, *) '     Direct radiance from reflexion (W/str/m**2/nm)'
  WRITE (2, 2001) rdirect
  WRITE (2, *) '           Cloud radiance (W/str/m**2/nm)         '
  WRITE (2, 2001) fctcld / omefov / (pi * (diamobj / 2.0)**2.0)
  WRITE (2, *) '         Diffuse radiance (W/str/m**2/nm)          '
  WRITE (2, 2001) (ftocap + fctcld) / omefov / (pi * (diamobj / 2.0)**2.0)
  CLOSE (2)
2001 FORMAT('                   ', E10.3E2)
  STOP
END PROGRAM illumina
!***********************************************************************************************************************
!*                                                                                                                     *
!*                                         end of the programme                                                        *
!*                                                                                                                     *
!***********************************************************************************************************************
