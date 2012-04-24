      SUBROUTINE SFCFLX(
     O                 FT,   TAUX,   TAUY,
     I                  T)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the surface forcing fluxes from the data given in the
c boundary condition files
c
c  HISTORY
c     '99.08.18  H.Hasumi
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '02.05.29  H.Nakano: tracer dimension
c     '02.06.02  H.Hasumi
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8      FT(NXYDIM, NTDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8       T(NXYDIM, NZDIM, NTDIM)

      INTEGER     IJ,      L,      N

      REAL*8  TSFC(NXYDIM, NTDIM), TDMP(NXYDIM, NTDIM)
      COMMON /WORK/ TSFC, TDMP

      CALL TMINTP(
     O              TAUX,
     I                 1)
      CALL TMINTP(
     O              TAUY,
     I                 2)
      DO L = 1, NTDIM
         N = L + 2
         CALL TMINTP(
     O                 TSFC(1, L),
     I                    N)
      END DO
      DO L = 1, NTDIM
         N = L + 2 + NTDIM
         CALL TMINTP(
     O                 TDMP(1, L),
     I                    N)
      END DO

      DO 100 IJ = IJSTR, IJEND
         FT(IJ, 1) = TDMP(IJ, 1) * (TSFC(IJ, 1) - T(IJ, KSTR, 1)) *
     &               DZ(IJ, KSTR) * AMSKT(IJ, KSTR)
         FT(IJ, 2) = TDMP(IJ, 2) * (TSFC(IJ, 2) - T(IJ, KSTR, 2)) *
     &               DZ(IJ, KSTR) / T(IJ, KSTR, 2) * AMSKT(IJ, KSTR)
  100 CONTINUE
      DO L = 3, NTDIM
         DO IJ = IJSTR, IJEND
            FT(IJ, L) = TDMP(IJ, L) * (TSFC(IJ, L) - T(IJ, KSTR, L)) *
     &                  DZ(IJ, KSTR) * AMSKT(IJ, KSTR)
         END DO
      END DO

      RETURN
      END
#ifdef OPT_BODY
c *********************************************************************

      SUBROUTINE BDYFLX(
     O                 TQ,
     I                  T)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the surface forcing fluxes from the data given in the
c boundary condition files
c
c  HISTORY
c     '99.08.25  H.Hasumi
c     '02.06.02  H.Hasumi: tracer dimension
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8      TQ(NXYDIM, NZDIM, NTDIM)
      REAL*8       T(NXYDIM, NZDIM, NTDIM)

      INTEGER     IJ,      K,      L,      N

      REAL*8  TBDY(NXYDIM, NZDIM, NTDIM), TDMB(NXYDIM, NZDIM, NTDIM)
      COMMON /WORK/ TBDY, TDMB

c      DO L = 1, NTDIM
      DO L = 1, 2
         N = L
         CALL TMINTB(
     O                 TBDY(1, 1, L),
     I                    N)
      END DO
c      DO L = 1, NTDIM
c         N = L + NTDIM
      DO L = 1, 2
         N = L + 2
         CALL TMINTB(
     O                 TDMB(1, 1, L),
     I                    N)
      END DO

c      DO 120 L = 1, NTDIM
      DO 120 L = 1, 2
         DO 110 K = KSTR, KEND
            DO 100 IJ = IJTSTR, IJTEND
               TQ(IJ, K, L) = TDMB(IJ, K, L) *
     &                        (TBDY(IJ, K, L) - T(IJ, K, L)) *
     &                        AMSKT(IJ, K)
 100        CONTINUE
 110     CONTINUE
 120  CONTINUE

      RETURN
      END
#endif
