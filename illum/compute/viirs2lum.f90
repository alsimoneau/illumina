SUBROUTINE viirs2lum(Nt, Nr, nzones, nangles, nwav, nbands, nsources, &
                     nlamps, nlops, nspcts, viirs, zones, angles, wav, bands, &
                     sens, lops, spcts, sources, ivtr, pixsize, reflect, lumlp)
  ! =====================================================
  ! Computes light sources power from the viirs image
  ! =====================================================

  USE constants

  IMPLICIT NONE

  ! Array lenghts
  INTEGER(4), INTENT(IN) :: Nt, Nr
  INTEGER(4), INTENT(IN) :: nzones
  INTEGER(4), INTENT(IN) :: nangles
  INTEGER(4), INTENT(IN) :: nwav
  INTEGER(4), INTENT(IN) :: nbands
  INTEGER(4), INTENT(IN) :: nsources
  INTEGER(4), INTENT(IN) :: nlamps
  INTEGER(4), INTENT(IN) :: nlops
  INTEGER(4), INTENT(IN) :: nspcts

  ! Inputs
  REAL(8), INTENT(IN) :: viirs(Nt, Nr)
  INTEGER(4), INTENT(IN) :: zones(Nt, Nr)
  REAL(8), INTENT(IN) :: angles(nangles)
  REAL(8), INTENT(IN) :: wav(nwav)
  LOGICAL, INTENT(IN) :: bands(nbands, nwav)
  REAL(8), INTENT(IN) :: sens(nwav)
  REAL(8), INTENT(IN) :: lops(nlops, nangles)
  REAL(8), INTENT(IN) :: spcts(nspcts, nwav)
  INTEGER(4), INTENT(IN) :: sources(nlops)
  REAL(8), INTENT(IN) :: ivtr(nlamps, 4)
  REAL(8), INTENT(IN) :: pixsize(Nr), reflect(nbands)

  ! Output
  REAL(8), INTENT(OUT) :: lumlp(nbands, nsources, Nt, Nr)

  ! Internal variables
  INTEGER(4) :: i, j, z, a, b, wl, s
  REAL(8) :: norm, acc
  REAL(8) :: lamps(nzones, nsources, nangles, nwav)
  REAL(8) :: mids(nangles + 1)
  REAL(8) :: sinx(nangles)
  REAL(8) :: Gdown(nzones, nsources, nwav), Gup(nzones, nsources, nwav)
  REAL(8) :: integral(nzones)
  REAL(8) :: phie(Nt, Nr)
  REAL(8) :: ratio(nbands, nzones, nsources)

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
        ASSOCIATE (l => lamps(INT(ivtr(i, 1)), j, a, wl))
          l = l + ivtr(i, 2) * spcts(INT(ivtr(i, 3)), wl) * lops(INT(ivtr(i, 4)), a)
        END ASSOCIATE
      END DO
    END DO
  END DO

  ! Ensuring normalization
  DO s = 1, nsources
    DO z = 1, nzones
      norm = 0
      DO wl = 1, nwav
        DO a = 1, nangles
          norm = norm + lamps(z, s, a, wl) * sinx(a)
        END DO
      END DO
      norm = norm * (wav(2) - wav(1))
      lamps(z, s, :, :) = lamps(z, s, :, :) / norm
    END DO
  END DO

  ! Inversion
  ! phie = DNB * S / int( R ( rho/pi Gdown + Gup ) ) dlambda
  norm = 0
  DO a = 1, nangles
    IF (angles(a) < 70) THEN
      norm = norm + sinx(a)
    END IF
  END DO
  DO wl = 1, nwav
    DO a = 1, nangles
      DO s = 1, nsources
        DO z = 1, nzones
          IF (angles(a) > 90) THEN
            Gdown(z, s, wl) = Gdown(z, s, wl) + lamps(z, s, a, wl) * sinx(a)
          END IF
          IF (angles(a) < 70) THEN
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
                      * (Gdown(z, s, wl) * reflect(wl) / PI + Gup(z, s, wl))
      END DO
    END DO
  END DO

  DO j = 1, Nr
    DO i = 1, Nt
      norm = 0
      DO z = 1, nzones
        IF (zones(i, j) == (z - 1)) THEN
          acc = integral(z)
          IF (acc == acc) THEN ! if not nan
            norm = norm + acc
          END IF
        END IF
        IF (norm == 0) THEN
          phie(i, j) = 0
        ELSE
          phie(i, j) = phie(i, j) + viirs(i, j) * pixsize(j) / norm
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
    ratio(b, :, :) = ratio(b, :, :) / norm
  END DO

  DO j = 1, Nr
    DO i = 1, Nt
      DO z = 1, nzones
        IF (zones(i, j) == (z - 1)) THEN
          lumlp(:, :, i, j) = lumlp(:, :, i, j) + ratio(:, z, :)
        END IF
      END DO
      lumlp(:, :, i, j) = lumlp(:, :, i, j) * phie(i, j)
    END DO
  END DO

END SUBROUTINE
