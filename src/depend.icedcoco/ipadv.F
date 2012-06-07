      SUBROUTINE PADVCT(
     M                 AX,    HIX,    EIX,    HSX,    TIX,
     O                 AZ,    HIZ,    EIZ,    HSZ,
     O                FIX,    FIY,    FSX,    FSY,
     O                FEX,    FEY,
     I                UIY,    VIY)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c ---- information ----------------------------------------------------
c
c  Horizontal advection of ice category fraction, ice volume, and
c snow volume. The current choice of the advection algorithm is the
c simplest first-order upstream. No explict diffusion is taken into
c account. Advection algorithm must be positivity preserving.
c
c  HISTORY
c     '03.08.06  H.Hasumi
c     '07.09.25  H.Hasumi: for COCO4
c     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.??.??  Y.Komuro
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.05.25  Y.Komuro: CMIP5 output code included
c     '09.09.04  Y.Komuro: extra output code (FEX/FEY)
c
c ---------------------------------------------------------------------

      REAL*8      AX(NXYDIM, 0:NIC),     AZ(NXYDIM, 0:NIC)
      REAL*8     HIX(NXYDIM, 0:NIC),    HIZ(NXYDIM, 0:NIC)
      REAL*8     EIX(NXYDIM, 0:NIC),    EIZ(NXYDIM, 0:NIC)
      REAL*8     HSX(NXYDIM, 0:NIC),    HSZ(NXYDIM, 0:NIC)
      REAL*8     TIX(NXYDIM, 0:NIC)
      REAL*8     UIY(NXYDIM),    VIY(NXYDIM)
      REAL*8     FIX(NXYDIM, 0:NIC),    FIY(NXYDIM, 0:NIC)
      REAL*8     FSX(NXYDIM, 0:NIC),    FSY(NXYDIM, 0:NIC)
      REAL*8     FEX(NXYDIM, 0:NIC),    FEY(NXYDIM, 0:NIC)

      REAL*8     FAX(NXYDIM,   NIC),    FAY(NXYDIM,   NIC)
      REAL*8   AXHIX(NXYDIM, 0:NIC),  AXHSX(NXYDIM, 0:NIC)
      REAL*8   AXEIX(NXYDIM, 0:NIC)
c      COMMON /WORK/ FAX, FAY, AXHIX, AXHSX, AXEIX

      REAL*8       U,     UP,     UM
      REAL*8       V,     VP,     VM
      REAL*8     TMI
      INTEGER     IJ,      K
      INTEGER   IJLS,   IJLW,   IJLN,   IJLE
      SAVE TMI

      INTEGER  IFPAR,  JFPAR
      LOGICAL OFIRST
      REAL*8      SI
      SAVE SI, OFIRST
      NAMELIST /NMISLT/ SI

!===== Define statement function 
#include "zocite.F"
!===== 

      DATA SI / 5.D0 /
      DATA OFIRST / .TRUE. /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMISLT, END=99)
 99      WRITE(JFPAR, NMISLT)
         TMI = dtds * SI
      END IF

      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            AZ (IJ, K) = AX (IJ, K)
            HIZ(IJ, K) = HIX(IJ, K)
            EIZ(IJ, K) = EIX(IJ, K)
            HSZ(IJ, K) = HSX(IJ, K)
            AXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            AXEIX(IJ, K) = AX(IJ, K) * EIX(IJ, K)
            AXHSX(IJ, K) = AX(IJ, K) * HSX(IJ, K)
            FIX(IJ, K) = 0.D0
            FIY(IJ, K) = 0.D0
            FSX(IJ, K) = 0.D0
            FSY(IJ, K) = 0.D0
            FEX(IJ, K) = 0.D0
            FEY(IJ, K) = 0.D0
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = 1, NXYDIM
            FAX(IJ, K) = 0.D0
            FAY(IJ, K) = 0.D0
         END DO
      END DO

      DO K = 1, NIC

         DO IJ = IJTSTR, IJTEND+1
            IJLW = IJ + LW
            U  = (  UIY(IJLW  ) * HYU(IJLW  )
     &            + UIY(IJ+LSW) * HYU(IJ+LSW)) * 0.25D0
            UP = U + ABS(U)
            UM = U - ABS(U)
            FAX(IJ, K) = - (  UP * AX(IJLW, K)
     &                   + UM * AX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLW, KSTR)
            FIX(IJ, K) = - (  UP * AXHIX(IJLW, K)
     &                   + UM * AXHIX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLW, KSTR)
            FEX(IJ, K) = - (  UP * AXEIX(IJLW, K)
     &                   + UM * AXEIX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLW, KSTR)
            FSX(IJ, K) = - (  UP * AXHSX(IJLW, K)
     &                   + UM * AXHSX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLW, KSTR)
         END DO

         DO IJ = IJTSTR, IJTEND+NXDIM
            IJLS = IJ + LS
            V  = (  VIY(IJLS  ) * HXU(IJLS  )
     &            + VIY(IJ+LSW) * HXU(IJ+LSW)) * 0.25D0
            VP = V + ABS(V)
            VM = V - ABS(V)
            FAY(IJ, K) = - (  VP * AX(IJLS, K)
     &                   + VM * AX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLS, KSTR)
            FIY(IJ, K) = - (  VP * AXHIX(IJLS, K)
     &                   + VM * AXHIX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLS, KSTR)
            FEY(IJ, K) = - (  VP * AXEIX(IJLS, K)
     &                   + VM * AXEIX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLS, KSTR)
            FSY(IJ, K) = - (  VP * AXHSX(IJLS, K)
     &                   + VM * AXHSX(IJ, K)) *
     &                  AMSKT(IJ, KSTR) * AMSKT(IJLS, KSTR)
         END DO

         DO IJ = IJTSTR, IJTEND
            IJLE = IJ + LE
            IJLN = IJ + LN
            AX(IJ, K) = AX(IJ, K)
     &                + TS * (  (FAX(IJLE, K) - FAX(IJ, K)) * RX
     &                        + (FAY(IJLN, K) - FAY(IJ, K)) * RY(IJ)) *
     &                  RXT(IJ) * RYT(IJ) * AMSKT(IJ, KSTR)
            AXHIX(IJ, K) = AXHIX(IJ, K)
     &                   + TS * (  (FIX(IJLE, K) - FIX(IJ, K)) * RX
     &                         + (FIY(IJLN, K) - FIY(IJ, K)) * RY(IJ)) *
     &                     RXT(IJ) * RYT(IJ) * AMSKT(IJ, KSTR)
            AXEIX(IJ, K) = AXEIX(IJ, K)
     &                   + TS * (  (FEX(IJLE, K) - FEX(IJ, K)) * RX
     &                         + (FEY(IJLN, K) - FEY(IJ, K)) * RY(IJ)) *
     &                     RXT(IJ) * RYT(IJ) * AMSKT(IJ, KSTR)
            AXHSX(IJ, K) = AXHSX(IJ, K)
     &                   + TS * (  (FSX(IJLE, K) - FSX(IJ, K)) * RX
     &                         + (FSY(IJLN, K) - FSY(IJ, K)) * RY(IJ)) *
     &                     RXT(IJ) * RYT(IJ) * AMSKT(IJ, KSTR)
         END DO

      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AX(IJ, K) .GT. 0.D0) THEN
               HIX(IJ, K) = AXHIX(IJ, K) / AX(IJ, K)
               EIX(IJ, K) = AXEIX(IJ, K) / AX(IJ, K)
               HSX(IJ, K) = AXHSX(IJ, K) / AX(IJ, K)
               TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
            ELSE
               HIX(IJ, K) = HIC(K)
               EIX(IJ, K) = 0.D0
               HSX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
            END IF
         END DO
      END DO

c     AX(IJ, 0) can be negative
      DO IJ = IJTSTR, IJTEND
         AX(IJ, 0) = 1.D0
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AX(IJ, 0) = AX(IJ, 0) - AX(IJ, K)
         END DO
      END DO

c      CALL REWNML(IFPAR, JFPAR)
c      DO K = 1, NIC
c         DO IJ = IJTSTR, IJTEND
c            IF (AX(IJ, K) .LT. 0.D0) THEN
c               WRITE(JFPAR, *) '### NEGATIVE AREA (ipadv) ###',
c     &         IJ, K, AX(IJ, K)
c               stop
c            END IF
c            IF (     (AXHIX(IJ, K) .LT. 0.D0)
c     &          .OR. (AXHSX(IJ, K) .LT. 0.D0)) THEN
c               WRITE(JFPAR, *) '### NEGATIVE ICE/SNOW (ipadv) ###',
c     &         IJ, K, AX(IJ, K), AXHIX(IJ, K), AXHSX(IJ, K)
c               stop
c            END IF
c            IF (HIX(IJ, K) .LE. 0.D0) THEN
c               WRITE(JFPAR, *) '### INVALID ICE THICKNESS (ipadv) ###',
c     &         IJ, K, AX(IJ, K), HIX(IJ, K)
c               stop
c            END IF
c         END DO
c      END DO

      RETURN
      END
