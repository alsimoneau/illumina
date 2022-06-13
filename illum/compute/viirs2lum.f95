subroutine viirs2lum(N, nzones, nangles, nwav, nbands, nsources, &
  nlamps, nlops, nspcts, viirs, zones, angles, wav, bands, &
  sens, lops, spcts, sources, ivtr, pixsize, reflect, lumlp)
  ! =====================================================
  ! Computes light sources power from the viirs image
  ! =====================================================

  implicit none

  ! Array lenghts
  integer, intent(in) :: N
  integer, intent(in) :: nzones
  integer, intent(in) :: nangles
  integer, intent(in) :: nwav
  integer, intent(in) :: nbands
  integer, intent(in) :: nsources
  integer, intent(in) :: nlamps
  integer, intent(in) :: nlops
  integer, intent(in) :: nspcts

  ! Inputs
  real,    intent(in) :: viirs(N,N)
  integer, intent(in) :: zones(N,N)
  real,    intent(in) :: angles(nangles)
  real,    intent(in) :: wav(nwav)
  logical, intent(in) :: bands(nwav,nbands)
  real,    intent(in) :: sens(nwav)
  real,    intent(in) :: lops(nlops,nangles)
  real,    intent(in) :: spcts(nspcts,nwav)
  integer, intent(in) :: sources(nlops)
  real,    intent(in) :: ivtr(nlamps,4)
  real,    intent(in) :: pixsize, reflect

  ! Output
  real, intent(out) :: lumlp(nbands,nsources,N,N)

  ! Internal variables
  integer :: i,j,z,a,b,wl,s
  real    :: norm
  real    :: lamps(nzones,nsources,nangles,nwav)
  real    :: mids(nangles + 1)
  real    :: sinx(nangles)
  real    :: Gdown(nzones,nsources,nwav), Gup(nzones,nsources,nwav)
  real    :: integral(nzones)
  real    :: phie(N,N)
  real    :: ratio(nbands,nzones,nsources)

  ! Defining usefull values
  real, parameter :: PI=4.D0*DATAN(1.D0)

  mids = 0
  do a = 2, nangles
    mids(a) = (angles(a-1) + angles(a))/2
  end do
  do a = 1, nangles
    sinx(a) = 2 * PI * (COS(mids(a+1)) - COS(mids(a)))
  end do

  ! Building lamp inventory
  do i = 1, nlamps
    do a = 1, nangles
      do wl = 1, nwav
        j = sources(ivtr(i,4))
        lamps(ivtr(i,1),j,a,wl) = lamps(ivtr(i,1),j,a,wl) &
          + ivtr(i,2) * spcts(ivtr(i,3),wl) * lops(ivtr(i,4),a)
      end do
    end do
  end do

  ! Ensuring normalization
  do z = 1, nzones
    do s = 1, nsources
      norm = 0
      do a = 1, nangles
        do wl = 1, nwav
          norm = norm + lamps(z,s,a,wl) * sinx(a)
        end do
      end do
      norm = norm * ( wav(2) - wav(1) )
      print *, "Zone norm:", norm
      do a = 1, nangles
        do wl = 1, nwav
          lamps(z,s,a,wl) = lamps(z,s,a,wl) / norm
        end do
      end do
    end do
  end do

  ! Inversion
  ! phie = DNB * S / int( R ( rho/pi Gdown + Gup ) ) dlambda
  norm = 0
  do a = 1, nangles
    if ( angles(a) .lt. 70 ) then
      norm = norm + sinx(a)
    end if
  end do
  do z = 1, nzones
    do s = 1, nsources
      do a = 1, nangles
        do wl = 1, nwav
          if ( angles(a) .gt. 90 ) then
            Gdown(z,s,wl) = Gdown(z,s,wl) + lamps(z,s,a,wl) * sinx(a)
          end if
          if ( angles(a) .lt. 70 ) then
            Gup(z,s,wl) = Gup(z,s,wl) + lamps(z,s,a,wl) * sinx(a) / norm
          end if
        end do
      end do
    end do
  end do

  do z = 1, nzones
    do s = 1, nsources
      do wl = 1, nwav
        integral(z) = integral(z) + sens(wl) * (wav(2) - wav(1)) &
          * ( Gdown(z,s,wl) * reflect / PI + Gup(z,s,wl) )
      end do
    end do
  end do

  do i = 1, N
    do j = 1, N
      norm = 0
      do z = 1, nzones
        if ( zones(i,j) .eq. (z - 1) ) then
          norm = norm + integral(z)
        end if
        if ( norm .eq. 0 ) then
          phie(i,j) = 0
        else
          phie(i,j) = phie(i,j) + viirs(i,j) * pixsize ** 2 / norm
        end if
      end do
    end do
  end do

  do b = 1, nbands
    norm = 0
    do wl = 1, nwav
      if ( bands(wl,b) ) then
        do z = 1, nzones
          do s = 1, nsources
            do a = 1, nangles
              ratio(b,z,s) = ratio(b,z,s) + lamps(z,s,a,wl) * sinx(a)
            end do
          end do
        end do
        norm = norm + 1
      end if
    end do
    do z = 1, nzones
      do s = 1, nsources
        ratio(b,z,s) = ratio(b,z,s) / norm
      end do
    end do
  end do

  do i = 1, N
    do j = 1, N
      do b = 1, nbands
        do s = 1, nsources
          do z = 1, nzones
            if ( zones(i,j) .eq. (z - 1) ) then
              lumlp(b,s,i,j) = lumlp(b,s,i,j) + ratio(b,z,s)
            end if
          end do
          lumlp(b,s,i,j) = lumlp(b,s,i,j) * phie(i,j)
        end do
      end do
    end do
  end do

end subroutine
