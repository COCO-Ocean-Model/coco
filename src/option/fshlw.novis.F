      SUBROUTINE MODGXY(
     M                GXX,    GYY,
     I               UBTX,   VBTX)
      IMPLICIT NONE

c ---- information ----------------------------------------------------
c
c  HISTORY
c     '03.04.23  H.Hasumi: from COCO3.4
c     '06.01.10  H.Hasumi: bug fix (calculation of SYY)
c     '07.04.23  H.Hasumi: Formulation changed
c     '07.04.26  H.Hasumi: No explicit viscosity version
c     '07.07.30  H.Hasumi: Loop boundary for S?? calculation reduced
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8      HX(NXYDIM),     HY(NXYDIM)
      REAL*8    UBTX(NXYDIM),   UBTY(NXYDIM)
      REAL*8    VBTX(NXYDIM),   VBTY(NXYDIM)
      REAL*8     GXX(NXYDIM),    GYY(NXYDIM)
      REAL*8    PTOP(NXYDIM)
      REAL*8      FW(NXYDIM)

      REAL*8      GH(NXYDIM)
      LOGICAL OFIRST
      SAVE GH, OFIRST

      REAL*8     FUX(NXYDIM),    FUY(NXYDIM)
      REAL*8     FVX(NXYDIM),    FVY(NXYDIM)
      REAL*8     FHX(NXYDIM),    FHY(NXYDIM)
      REAL*8      GU(NXYDIM),     GV(NXYDIM)
      REAL*8     SXX(NXYDIM),    SYY(NXYDIM)
      REAL*8     SXY(NXYDIM),    SYX(NXYDIM)
      COMMON /WORK/ FUX, FUY, FVX, FVY, FHX, FHY, GU, GV,
     &              SXX, SYY, SXY, SYX

      REAL*8      CF
      INTEGER     IJ
      INTEGER   IJLE,   IJLN,   IJLW,   IJLS
      INTEGER   IJNW,   IJSE,   IJSW
      INTEGER  IFPAR,  JFPAR

      REAL*8    ACCB,    ACC
      SAVE ACCB
      NAMELIST /NMACCB/ ACCB
      NAMELIST /NMACCV/ ACC

      DATA OFIRST / .TRUE. /
      DATA ACCB / -1.D0 /
      DATA ACC / 1.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMACCB, END=99)
 99      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMACCV, END=98)
 98      CONTINUE

         IF (ACCB .LE. 0.D0) THEN
            ACCB = ACC
         END IF

         DO 10 IJ = 1, NXYDIM
            GH   (IJ) = gravit / RDEPV(IJ)
 10      CONTINUE
      END IF

      RETURN

c =====================================================================

      ENTRY SHALOW(
     M                 HX,   UBTX,   VBTX,
     I                 HY,   UBTY,   VBTY,
     I                GXX,    GYY,   PTOP,  FW)

      DO IJ = 1, NXYDIM
         FUX(IJ) = 0.D0
         FUY(IJ) = 0.D0
         FVX(IJ) = 0.D0
         FVY(IJ) = 0.D0
         FHX(IJ) = 0.D0
         FHY(IJ) = 0.D0
      END DO

      DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+2
         FHX(IJ) =
     &           - (  UBTY(IJ+LW) * HYU(IJ+LW)
     &              + UBTY(IJ+LSW) * HYU(IJ+LSW)) * 0.5D0
      END DO
      DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+NXDIM+1
         FHY(IJ) =
     &           - (  VBTY(IJ+LS) * HXU(IJ+LS)
     &              + VBTY(IJ+LSW) * HXU(IJ+LSW)) * 0.5D0
      END DO
      DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
         HX(IJ) = HX(IJ)
     &          + TSS * (  (FHX(IJ+LE) - FHX(IJ)) * RX
     &                   + (FHY(IJ+LN) - FHY(IJ)) * RY(IJ)) *
     &            RXT(IJ) * RYT(IJ) * AMSKT(IJ, KSTR)
     &          - TSS * FW(IJ) * AMSKT(IJ, KSTR)
      END DO

C      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1 ! 07.07.30
      DO IJ = IJSTR-NXDIM-1, IJEND
         GU(IJ) = GXX(IJ) + COR(IJ) * VBTX(IJ)
     &          - (  (  HY(IJ+LNE) + HY(IJ+LE)
     &                - HY(IJ+LN ) - HY(IJ   )) * GH(IJ)
     &             + (  PTOP(IJ+LNE) + PTOP(IJ+LE)
     &                - PTOP(IJ+LN ) - PTOP(IJ   ))
     &               / RDEPV(IJ) / rhoo
     &            ) * 0.5D0 * RX * RXU(IJ)
         GV(IJ) = GYY(IJ) - COR(IJ) * UBTX(IJ)
     &          - (  (  HY(IJ+LNE) + HY(IJ+LN)
     &                - HY(IJ+LE ) - HY(IJ   )) * GH(IJ)
     &             + (  PTOP(IJ+LNE) + PTOP(IJ+LN)
     &                - PTOP(IJ+LE ) - PTOP(IJ   ))
     &               / RDEPV(IJ) / rhoo
     &            ) * 0.5D0 * RYM(IJ) * RYU(IJ)
      END DO

C      DO 1800 IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1 ! 07.07.30
      DO 1800 IJ = IJSTR-NXDIM-1, IJEND
         CF = COR(IJ) * TSS / ACCB * 0.5D0
         UBTX(IJ) = UBTX(IJ)
     &            + TSS / ACCB / (1.D0 + CF * CF) *
     &              (GU(IJ) + CF * GV(IJ)) * AMSKV(IJ, KSTR)
         VBTX(IJ) = VBTX(IJ)
     &            + TSS / ACCB / (1.D0 + CF * CF) *
     &              (GV(IJ) - CF * GU(IJ)) * AMSKV(IJ, KSTR)
 1800 CONTINUE

      RETURN
      END
