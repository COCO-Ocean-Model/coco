      SUBROUTINE ICTRNS(
     M                 AX,    HIX,    HSX,    EIX,    TIX,
     M                 FT,     FS)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
c      use zocite
      IMPLICIT NONE

c ---- information ----------------------------------------------------
c
c  Sea ice thickness category transfer
c
c  HISTORY
c     '03.07.30  H.Hasumi
c     '07.09.25  H.Hasumi: for COCO4
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

      REAL*8      AX(NXYDIM, 0:NIC)
      REAL*8     HIX(NXYDIM, 0:NIC)
      REAL*8     HSX(NXYDIM, 0:NIC)
      REAL*8     EIX(NXYDIM, 0:NIC),    TIX(NXYDIM, 0:NIC)
      REAL*8      FT(NXYDIM, NTDIM),     FS(NXYDIM)

      REAL*8   AXHIX(NXYDIM, 0:NIC),  AXHSX(NXYDIM, 0:NIC)
      REAL*8   AXEIX(NXYDIM, 0:NIC)
      REAL*8      CI(NXYDIM)
c      COMMON /WORK/ AXHIX, AXHSX, AXEIX, CI

      REAL*8     RRI,    RRS
      REAL*8     TMI
      LOGICAL OFIRST
      SAVE RRI, RRS, TMI, OFIRST

      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      REAL*8    AMIN,   AMAX,     SI
      INTEGER    MIC
      NAMELIST /NMAMIN/ AMIN, AMAX, MIC
      NAMELIST /NMISLT/ SI
      SAVE AMIN, SI

#include "zocite.F"

      DATA OFIRST / .TRUE. /
      DATA AMIN, AMAX / 1.0D-6, 1.D0 /
      DATA MIC / NIC /
      DATA SI / 5.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMAMIN, END=98)
 98      WRITE(JFPAR, NMAMIN)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMISLT, END=97)
 97      WRITE(JFPAR, NMISLT)

         TMI = dtds * SI
         RRI    = rhoo / rhoi
         RRS    = rhoo / rhos
      END IF

      DO IJ = IJTSTR, IJTEND
         AXHIX(IJ, 0) = AX(IJ, 0) * HIX(IJ, 0)
         AXHSX(IJ, 0) = AX(IJ, 0) * HSX(IJ, 0)
         AXEIX(IJ, 0) = AX(IJ, 0) * EIX(IJ, 0)
      END DO

      DO K = 1, NIC-1
         DO IJ = IJTSTR, IJTEND
            IF (HIX(IJ, K) .GE. HIC(K+1)) THEN
c            IF (     (AX(IJ, K) .GT. 0.D0)
c     &         .AND. (HIX(IJ, K) .GE. HIC(K+1))) THEN
               AXHIX(IJ, K+1) = AX(IJ, K+1) * HIX(IJ, K+1)
     &                        + AX(IJ, K) * HIX(IJ, K)
               AXHSX(IJ, K+1) = AX(IJ, K+1) * HSX(IJ, K+1)
     &                        + AX(IJ, K) * HSX(IJ, K)
               AXEIX(IJ, K+1) = AX(IJ, K+1) * EIX(IJ, K+1)
     &                        + AX(IJ, K) * EIX(IJ, K)
               AX(IJ, K+1) = AX(IJ, K+1) + AX(IJ, K)
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               EIX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
               HIX(IJ, K+1) = AXHIX(IJ, K+1) / AX(IJ, K+1)
               HSX(IJ, K+1) = AXHSX(IJ, K+1) / AX(IJ, K+1)
               EIX(IJ, K+1) = AXEIX(IJ, K+1) / AX(IJ, K+1)
               TIX(IJ, K+1) = TI(EIX(IJ, K+1)/HIX(IJ, K+1), SI)
            END IF
         END DO
      END DO
      DO K = NIC-1, 1, -1
         DO IJ = IJTSTR, IJTEND
            IF (HIX(IJ, K+1) .LT. HIC(K+1)) THEN
c            IF (      (HIX(IJ, K+1) .LT. HIC(K+1))
c     &          .AND. (AX(IJ, K+1) .GT. 0.D0)) THEN
               AXHIX(IJ, K) = AX(IJ, K+1) * HIX(IJ, K+1)
     &                      + AX(IJ, K) * HIX(IJ, K)
               AXHSX(IJ, K) = AX(IJ, K+1) * HSX(IJ, K+1)
     &                      + AX(IJ, K) * HSX(IJ, K)
               AXEIX(IJ, K) = AX(IJ, K+1) * EIX(IJ, K+1)
     &                      + AX(IJ, K) * EIX(IJ, K)
               AX(IJ, K) = AX(IJ, K+1) + AX(IJ, K)
               AX(IJ, K+1) = 0.D0
               HIX(IJ, K) = AXHIX(IJ, K) / AX(IJ, K)
               HSX(IJ, K) = AXHSX(IJ, K) / AX(IJ, K)
               EIX(IJ, K) = AXEIX(IJ, K) / AX(IJ, K)
               TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
               HIX(IJ, K+1) = HIC(K+1)
               HSX(IJ, K+1) = 0.D0
               EIX(IJ, K+1) = 0.D0
               TIX(IJ, K+1) = TMI
            END IF
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         IF (HIX(IJ, 1) .LT. HIC(1)) THEN
c         IF (      (AX(IJ, 1) .GT. 0.D0)
c     &       .AND. (HIX(IJ, 1) .LT. HIC(1))) THEN
            AXHIX(IJ, 1) = AX(IJ, 1) * HIX(IJ, 1)
            AXHSX(IJ, 1) = AX(IJ, 1) * HSX(IJ, 1)
            AXEIX(IJ, 1) = AX(IJ, 1) * EIX(IJ ,1)
            HIX(IJ, 1) = HIC(1)
            AX(IJ, 1) = AXHIX(IJ, 1) / HIX(IJ, 1)
            HSX(IJ, 1) = AXHSX(IJ, 1) / AX(IJ, 1)
            EIX(IJ, 1) = AXEIX(IJ ,1) / AX(IJ, 1)
*           ice temperature does not change
         END IF
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF ((AX(IJ, K) .GT. 0.D0) .AND. (AX(IJ, K) .LT. AMIN)) THEN
               AXHIX(IJ, 0) = AXHIX(IJ, 0)
     &                      + AX(IJ, K) * HIX(IJ, K)
               AXHSX(IJ, 0) = AXHSX(IJ, 0)
     &                      + AX(IJ, K) * HSX(IJ, K)
               AXEIX(IJ, 0) = AXEIX(IJ, 0)
     &                      + AX(IJ, K) * EIX(IJ, K)
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               EIX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
            END IF
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND
c         IF (HIX(IJ, 0) .GT. 0.D0) THEN
            FT(IJ, 2) = FT(IJ, 2)
     &                - (AXHIX(IJ, 0) / RRI + AXHSX(IJ, 0) / RRS)
     &                  / TS * AMSKT(IJ, KSTR)
            FS(IJ) = FS(IJ)
     &             - AXHIX(IJ, 0) / RRI / TS * SI * AMSKT(IJ, KSTR)
            FT(IJ, 1) = FT(IJ, 1)
     &                - (  AXHSX(IJ, 0) / RRS * hfus
     &                   + AXEIX(IJ, 0) / RRI) / TS / cpo *
     &                  AMSKT(IJ, KSTR)
c         END IF
      END DO

c      ENTRY IC0SET(
c     M                 AX,    HIX,    HSX)

*** initialize the category 0 ***
      DO IJ = IJTSTR, IJTEND
         AX(IJ, 0) = 1.D0
         HIX(IJ, 0) = 0.D0
         HSX(IJ, 0) = 0.D0
         EIX(IJ, 0) = 0.D0
         TIX(IJ, 0) = TMI
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AX(IJ, 0) = AX(IJ, 0) - AX(IJ, K)
         END DO
      END DO

*** because of the finit precision, AX(IJ, 0) could be a small negative
*** value, which causes some problems
      DO IJ = IJTSTR, IJTEND
         AX(IJ, 0) = MAX(0.D0, AX(IJ ,0))
      END DO

      RETURN

c =====================================================================

      ENTRY ICADJS(
     M                AX,    HIX,    HSX,    EIX)

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (AMAX .GE. 1.D0) RETURN

      DO K = 1, MIC
         DO IJ = IJTSTR, IJTEND
            AXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            AXHSX(IJ, K) = AX(IJ, K) * HSX(IJ, K)
            AXEIX(IJ, K) = AX(IJ, K) * EIX(IJ, K)
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         CI(IJ) = 0.D0
         AX(IJ, 0) = 1.D0
      END DO
      DO K = 1, MIC
         DO IJ = IJTSTR, IJTEND
            CI(IJ) = CI(IJ) + AX(IJ, K)
         END DO
      END DO

      DO K = 1, MIC
         DO IJ = IJTSTR, IJTEND
            IF ((CI(IJ) .GT. AMAX) .AND. (AX(IJ, K) .GT. 0.D0)) THEN
               AX(IJ, K) = AX(IJ, K) * AMAX / CI(IJ)
               HIX(IJ, K) = AXHIX(IJ, K) / AX(IJ, K)
               HSX(IJ, K) = AXHSX(IJ, K) / AX(IJ, K)
               EIX(IJ, K) = AXEIX(IJ, K) / AX(IJ, K)
            END IF
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AX(IJ, 0) = AX(IJ, 0) - AX(IJ, K)
         END DO
      END DO

      RETURN

      END
