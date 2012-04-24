      SUBROUTINE VELTAD(
     O               UADV,   VADV,
     I                UBT,    VBT,      U,      V)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Velocity for tracer advection
c
c  HISTORY
c     '99.08.16  H.Hasumi
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '07.04.23  H.Hasumi
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.02.23  Y.Komuro: bug fix by Dr. Kurogi
c                          (initialize UADV/VADV)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8       U(NXYDIM, NZDIM),      V(NXYDIM, NZDIM)
      REAL*8     UBT(NXYDIM),    VBT(NXYDIM)

      REAL*8    UAVR(NXYDIM),   VAVR(NXYDIM)
      REAL*8    UBAR(NXYDIM),   VBAR(NXYDIM)
      COMMON /WORK/ UAVR, VAVR, UBAR, VBAR

      INTEGER     IJ,      K

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            UADV(IJ, K) = 0.0D0
            VADV(IJ, K) = 0.0D0
         END DO
      END DO

      DO 100 IJ = IJVSTR, IJVEND
         UAVR(IJ) = U(IJ, KSTR) * DZV(IJ, KSTR) * AMSKV(IJ, KSTR)
         VAVR(IJ) = V(IJ, KSTR) * DZV(IJ, KSTR) * AMSKV(IJ, KSTR)
  100 CONTINUE
      DO 110 K = KSTR+1, KEND
         DO 120 IJ = IJVSTR, IJVEND
            UAVR(IJ) = UAVR(IJ) + U(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            VAVR(IJ) = VAVR(IJ) + V(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
  120    CONTINUE
  110 CONTINUE
      DO 130 IJ = IJVSTR, IJVEND
         UAVR(IJ) = UAVR(IJ) * RDEPV(IJ)
         VAVR(IJ) = VAVR(IJ) * RDEPV(IJ)
         UBAR(IJ) = UBT(IJ) * RDEPV(IJ)
         VBAR(IJ) = VBT(IJ) * RDEPV(IJ)
  130 CONTINUE

      DO 210 K = KSTR, KEND
         DO 200 IJ = IJVSTR, IJVEND
            UADV(IJ, K) = (UBAR(IJ) + (U(IJ, K) - UAVR(IJ))) *
     &                    AMSKV(IJ, K)
            VADV(IJ, K) = (VBAR(IJ) + (V(IJ, K) - VAVR(IJ))) *
     &                    AMSKV(IJ, K)
  200    CONTINUE
  210 CONTINUE

      RETURN
      END
