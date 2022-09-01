SUBROUTINE viirs2lum(N, nzones, nangles, nwav, nbands, nsources, &
                     nlamps, nlops, nspcts, viirs, zones, angles, wav, bands, &
                     sens, lops, spcts, sources, ivtr, pixsize, reflect, lumlp)
   ! =====================================================
   ! Computes light sources power from the viirs image
   ! =====================================================

   IMPLICIT NONE

   ! Array lenghts
   INTEGER, INTENT(in) :: N
   INTEGER, INTENT(in) :: nzones
   INTEGER, INTENT(in) :: nangles
   INTEGER, INTENT(in) :: nwav
   INTEGER, INTENT(in) :: nbands
   INTEGER, INTENT(in) :: nsources
   INTEGER, INTENT(in) :: nlamps
   INTEGER, INTENT(in) :: nlops
   INTEGER, INTENT(in) :: nspcts

   ! Inputs
   REAL, INTENT(in) :: viirs(N, N)
   INTEGER, INTENT(in) :: zones(N, N)
   REAL, INTENT(in) :: angles(nangles)
   REAL, INTENT(in) :: wav(nwav)
   LOGICAL, INTENT(in) :: bands(nbands, nwav)
   REAL, INTENT(in) :: sens(nwav)
   REAL, INTENT(in) :: lops(nlops, nangles)
   REAL, INTENT(in) :: spcts(nspcts, nwav)
   INTEGER, INTENT(in) :: sources(nlops)
   REAL, INTENT(in) :: ivtr(nlamps, 4)
   REAL, INTENT(in) :: pixsize, reflect

   ! Output
   REAL, INTENT(out) :: lumlp(nbands, nsources, N, N)

   ! Internal variables
   INTEGER :: i, j, z, a, b, wl, s
   REAL :: norm, acc
   REAL :: lamps(nzones, nsources, nangles, nwav)
   REAL :: mids(nangles + 1)
   REAL :: sinx(nangles)
   REAL :: Gdown(nzones, nsources, nwav), Gup(nzones, nsources, nwav)
   REAL :: integral(nzones)
   REAL :: phie(N, N)
   REAL :: ratio(nbands, nzones, nsources)

   ! Defining usefull values
   REAL, PARAMETER :: PI = 4.D0 * DATAN(1.D0)

   mids = 0
   DO a = 2, nangles
      mids(a) = (angles(a - 1) + angles(a)) / 2
   END DO
   DO a = 1, nangles
      sinx(a) = 2 * PI * (COS(mids(a + 1)) - COS(mids(a)))
   END DO

   ! Building lamp inventory
   DO wl = 1, nwav
      DO a = 1, nangles
         DO i = 1, nlamps
            j = sources(INT(ivtr(i, 4)))
            lamps(INT(ivtr(i, 1)), j, a, wl) = lamps(INT(ivtr(i, 1)), j, a, wl) &
                                               + ivtr(i, 2) * spcts(INT(ivtr(i, 3)), wl) * lops(INT(ivtr(i, 4)), a)
         END DO
      END DO
   END DO

   ! Ensuring normalization
   DO z = 1, nzones
      DO s = 1, nsources
         norm = 0
         DO a = 1, nangles
            DO wl = 1, nwav
               norm = norm + lamps(z, s, a, wl) * sinx(a)
            END DO
         END DO
         norm = norm * (wav(2) - wav(1))
         DO a = 1, nangles
            DO wl = 1, nwav
               lamps(z, s, a, wl) = lamps(z, s, a, wl) / norm
            END DO
         END DO
      END DO
   END DO

   ! Inversion
   ! phie = DNB * S / int( R ( rho/pi Gdown + Gup ) ) dlambda
   norm = 0
   DO a = 1, nangles
      IF (angles(a) .LT. 70) THEN
         norm = norm + sinx(a)
      END IF
   END DO
   DO wl = 1, nwav
      DO a = 1, nangles
         DO s = 1, nsources
            DO z = 1, nzones
               IF (angles(a) .GT. 90) THEN
                  Gdown(z, s, wl) = Gdown(z, s, wl) + lamps(z, s, a, wl) * sinx(a)
               END IF
               IF (angles(a) .LT. 70) THEN
                  Gup(z, s, wl) = Gup(z, s, wl) + lamps(z, s, a, wl) * sinx(a) / norm
               END IF
            END DO
         END DO
      END DO
   END DO

   DO wl = 1, nwav
      DO s = 1, nsources
         DO z = 1, nzones
            integral(z) = integral(z) + sens(wl) * (wav(2) - wav(1)) &
                          * (Gdown(z, s, wl) * reflect / PI + Gup(z, s, wl))
         END DO
      END DO
   END DO

   DO j = 1, N
      DO i = 1, N
         norm = 0
         DO z = 1, nzones
            IF (zones(i, j) .EQ. (z - 1)) THEN
               acc = integral(z)
               IF (acc .EQ. acc) THEN ! if not nan
                  norm = norm + acc
               END IF
            END IF
            IF (norm .EQ. 0) THEN
               phie(i, j) = 0
            ELSE
               phie(i, j) = phie(i, j) + viirs(i, j) * pixsize**2 / norm
            END IF
         END DO
      END DO
   END DO

   DO b = 1, nbands
      norm = 0
      DO wl = 1, nwav
         IF (bands(b, wl)) THEN
            DO a = 1, nangles
               DO s = 1, nsources
                  DO z = 1, nzones
                     ratio(b, z, s) = ratio(b, z, s) + lamps(z, s, a, wl) * sinx(a)
                  END DO
               END DO
            END DO
            norm = norm + 1
         END IF
      END DO
      DO s = 1, nsources
         DO z = 1, nzones
            ratio(b, z, s) = ratio(b, z, s) / norm
         END DO
      END DO
   END DO

   DO j = 1, N
      DO i = 1, N
         DO s = 1, nsources
            DO z = 1, nzones
               DO b = 1, nbands
                  IF (zones(i, j) .EQ. (z - 1)) THEN
                     lumlp(b, s, i, j) = lumlp(b, s, i, j) + ratio(b, z, s)
                  END IF
               END DO
               lumlp(b, s, i, j) = lumlp(b, s, i, j) * phie(i, j)
            END DO
         END DO
      END DO
   END DO

END SUBROUTINE
