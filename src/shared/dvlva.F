      SUBROUTINE VELVAD(
     O               UADV,   VADV,   WADV,
     I                  U,      V,      W)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Velocity for momentum advection
c
c  HISTORY
c     '99.08.17  H.Hasumi
c     '01.01.25  H.Hasumi: estimate of WADV is corrected
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '01.02.13  H.Hasumi: for incorporating BBL model
c     '01.05.02  H.Hasumi
c     '07.04.23  H.Hasumi
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.02.23  Y.Komuro: bug fix by Dr. Kurogi
c                          (loop length extended)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM)
      REAL*8       U(NXYDIM, NZDIM),      V(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM)

      INTEGER     IJ,      K

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            UADV(IJ, K) = 0.D0
            VADV(IJ, K) = 0.D0
            WADV(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO 210 K = KSTR, KSTR+KZ-1
         DO 200 IJ = IJTSTR, IJTEND
            WADV(IJ, K) = (  (W(IJ+LNE, K) + W(IJ+LN, K)) * DY(IJ+LN)
     &                     + (W(IJ+LE , K) + W(IJ   , K)) * DY(IJ)) *
c     &                    0.25D0 * RYM(IJ) * AMFVZ(IJ, K)
     &                    0.25D0 * RYM(IJ) * AMSKV(IJ, K)
            VADV(IJ, K) = (  V(IJ+LE, K) + V(IJ+LW, K)
     &                     + 2.D0 * V(IJ, K)
     &                     + V(IJ+LSE, K) + V(IJ+LSW, K)
     &                     + 2.D0 * V(IJ+LS, K)) *
     &                    0.125D0 * AMSKT(IJ, K) * AMSKT(IJ+LE, K)
            UADV(IJ, K) = (  (U(IJ+LNW, K) + U(IJ+LN, K)) * DY(IJ+LN)
     &                     + (U(IJ+LSW, K) + U(IJ+LS, K)) * DY(IJ)
     &                     + (U(IJ+LW , K) + U(IJ   , K)) * DYM(IJ) *
     &                       2.D0
     &                    ) * 0.125 * RYM(IJ) *
     &                    AMSKT(IJ, K) * AMSKT(IJ+LN, K)
  200    CONTINUE
  210 CONTINUE

      DO 310 K = KSTR+KZ, KEND
         DO 300 IJ = IJTSTR, IJTEND
            WADV(IJ, K) = (  (W(IJ+LNE, K) + W(IJ+LN, K)) * DY(IJ+LN)
     &                     + (W(IJ+LE , K) + W(IJ   , K)) * DY(IJ)) *
c     &                    0.25D0 * RYM(IJ) * AMFVZ(IJ, K)
     &                    0.25D0 * RYM(IJ) * AMSKV(IJ, K)
            VADV(IJ, K) = (  V(IJ+LE , K) * DZV(IJ+LE , K)
     &                     + V(IJ+LW , K) * DZV(IJ+LW , K)
     &                     + V(IJ    , K) * DZV(IJ    , K) * 2.D0
     &                     + V(IJ+LSE, K) * DZV(IJ+LSE, K)
     &                     + V(IJ+LSW, K) * DZV(IJ+LSW, K)
     &                     + V(IJ+LS , K) * DZV(IJ+LS , K) * 2.D0) *
     &                    0.125D0 * AMSKT(IJ, K) * AMSKT(IJ+LE, K)
            UADV(IJ, K) = (  (  U(IJ+LNW, K) * DZV(IJ+LNW, K)
     &                        + U(IJ+LN , K) * DZV(IJ+LN , K)
     &                        + U(IJ+LW , K) * DZV(IJ+LW , K)
     &                        + U(IJ    , K) * DZV(IJ    , K)) *
     &                       DY(IJ+LN)
     &                     + (  U(IJ    , K) * DZV(IJ    , K)
     &                        + U(IJ+LW , K) * DZV(IJ+LW , K)
     &                        + U(IJ+LS , K) * DZV(IJ+LS , K)
     &                        + U(IJ+LSW, K) * DZV(IJ+LSW, K)) *
     &                       DY(IJ)
     &                     ) * RYM(IJ) * 0.125D0 *
     &                     AMSKT(IJ, K) * AMSKT(IJ+LN, K)
  300    CONTINUE
  310 CONTINUE

      RETURN
      END
