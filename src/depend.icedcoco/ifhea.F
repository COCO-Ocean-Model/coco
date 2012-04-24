      SUBROUTINE FIHEAT(
     M                FT,
     O               WAO,    WIO,    WAS,    WIL,
     I                 A,      T,     SH,
     I               QAO,    QAI,    QIO,    QII,  SWABS)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '02.10.10  H.Hasumi: from MIROC3.1
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: multicategory sea ice
c     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8     WAO(NXYDIM)
      REAL*8     WIO(NXYDIM, NIC),    WAS(NXYDIM, NIC)
      REAL*8     WIL(NXYDIM, NIC)
      REAL*8       A(NXYDIM, 0:NIC)
      REAL*8     QIO(NXYDIM, NIC),    QAI(NXYDIM, NIC)
      REAL*8     QII(NXYDIM, NIC)
      REAL*8     QAO(NXYDIM),  SWABS(NXYDIM)
      REAL*8      FT(NXYDIM, NTDIM)
      REAL*8       T(NXYDIM, NZDIM, NTDIM)
      REAL*8      SH(NXYDIM)
      REAL*8  SWCONV(NXYDIM)

      REAL*8  SWCNV1(NXYDIM)
      LOGICAL OFIRST
      SAVE SWCNV1, OFIRST

      REAL*8    WFRZ(NXYDIM),    WIB(NXYDIM),   WILM(NXYDIM)
      COMMON /WORK/ WFRZ, WIB, WILM

      REAL*8    TDEV
      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      REAL*8   FBTAB
      SAVE FBTAB
      NAMELIST /NMBTAB/ FBTAB

      DATA OFIRST / .TRUE. /
      DATA FBTAB / 0.5D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBTAB, END=99)
 99      WRITE(JFPAR, NMBTAB)
      END IF

      DO IJ = IJTSTR, IJTEND
         TDEV    = T(IJ, KSTR, 1) - dtds * T(IJ, KSTR, 2)
         WFRZ(IJ) = - rhoo * cpo * TDEV *
     &               (SH(IJ) + ZBOT) * DS(KSTR) / TS
         FT (IJ, 1) = TDEV * (SH(IJ) + ZBOT) * DS(KSTR) / TS
         WAO(IJ) = A(IJ, 0) * WFRZ(IJ)
     &           + QAO(IJ) - SWCNV1(IJ) * SWABS(IJ)
      END DO

      DO IJ = IJTSTR, IJTEND
         IF ((A(IJ, 0) .LT. 1.D0) .AND. (WAO(IJ) .LT. 0.D0)) THEN
            WIB(IJ) = FBTAB * WAO(IJ) / (1.D0 - A(IJ, 0))
            WILM(IJ) = (1.D0 - FBTAB) * WAO(IJ) / (1.D0 - A(IJ, 0))
            WAO(IJ) = 0.D0
         ELSE
            WIB(IJ) = 0.D0
            WILM(IJ) = 0.D0
         END IF
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            WAS(IJ, K) = QAI(IJ, K) - QII(IJ, K)
            WIO(IJ, K) = A(IJ, K) * (WFRZ(IJ) + WIB(IJ))
     &                 + QIO(IJ, K)
            WIL(IJ, K) = A(IJ, K) * WILM(IJ)
         END DO
      END DO

      RETURN

c =====================================================================

      ENTRY PUTSWC(
     I             SWCONV)

      DO 10 IJ = 1, NXYDIM
         SWCNV1(IJ) = SWCONV(IJ)
   10 CONTINUE

      RETURN
      END
