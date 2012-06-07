      SUBROUTINE ICETMP(
     O               EIX,
     M               TIX,
     M               QAO,    QAI,    QIO,    QII,
     I                AX,    HIX)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '07.09.30  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8     EIX(NXYDIM, 0:NIC),    TIX(NXYDIM, 0:NIC)
      REAL*8     QAO(NXYDIM)
      REAL*8     QAI(NXYDIM, NIC),    QIO(NXYDIM, NIC)
      REAL*8     QII(NXYDIM, NIC)
      REAL*8      AX(NXYDIM, 0:NIC),    HIX(NXYDIM, 0:NIC)

      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8      SI
      NAMELIST /NMISLT/ SI
      SAVE SI

!===== Define statement function 
#include "zocite.F"
!===== 

      DATA OFIRST / .TRUE. /
      DATA SI / 5.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMISLT, END=97)
 97      WRITE(JFPAR, NMISLT)
      END IF

      DO K = 0, NIC ! zero-set for K = 0 (HI = 0 for K = 0)
         DO IJ = 1, NXYDIM
            EIX(IJ, K) = EI(TIX(IJ, K), SI) * HIX(IJ, K)
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            EIX(IJ, K) = EIX(IJ, K)
     &                 + TS * (QII(IJ, K) - QIO(IJ, K)) / rhoi
            TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         QAO(IJ) = QAO(IJ) * AX(IJ, 0)
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            QAI(IJ, K) = QAI(IJ, K) * AX(IJ, K)
            QIO(IJ, K) = QIO(IJ, K) * AX(IJ, K)
            QII(IJ, K) = QII(IJ, K) * AX(IJ, K)
         END DO
      END DO

      RETURN
      END
