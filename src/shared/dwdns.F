      SUBROUTINE WDENST(
     O                  W,
     I                 UX,     VX,
     I                 HX,     HZ)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the vertical velocity at T-points
c
c  HISTORY
c     '03.04.23  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HX(NXYDIM),            HZ(NXYDIM)
      REAL*8       W(NXYDIM, NZDIM)

      REAL*8     FTX(NXYDIM, NZDIM),    FTY(NXYDIM, NZDIM)
      REAL*8    DHDT(NXYDIM), RHZBOT(NXYDIM)
c      COMMON /WORK/ FTX , FTY, DHDT, RHZBOT

      INTEGER     IJ,      K
      INTEGER   IJLW,   IJLS,  IJLSW

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            FTX(IJ, K) = 0.D0
            FTY(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO IJ = 1, NXYDIM
         DHDT(IJ) = (HX(IJ) - HZ(IJ)) / TS
         RHZBOT(IJ) = 1.D0 / (HZ(IJ) + ZBOT)
      END DO

      DO K = KSTR, KEND
         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+NXDIM+1
            IJLW = IJ + LW
            IJLS = IJ + LS
            IJLSW = IJ + LSW
            FTX(IJ, K) = (  UX(IJLW, K) * DZV(IJLW, K) *
     &                      HYU(IJLW)
     &                    + UX(IJLSW, K) * DZV(IJLSW, K) *
     &                      HYU(IJLSW)) *
     &                   AMSKT(IJ, K) * AMSKT(IJLW, K)
            FTY(IJ, K) = (  VX(IJLS, K) * DZV(IJLS, K) *
     &                      HXU(IJLS)
     &                    + VX(IJLSW, K) * DZV(IJLSW, K) *
     &                      HXU(IJLSW)) *
     &                   AMSKT(IJ, K) * AMSKT(IJLS, K)
         END DO
      END DO

      DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
         W(IJ, KEND) = - (  (FTX(IJ+LE, KEND) - FTX(IJ, KEND)) *
     &                      RX
     &                    + (FTY(IJ+LN, KEND) - FTY(IJ, KEND)) *
     &                      RY(IJ)) * RXT(IJ) * RYT(IJ) * 0.5D0
      END DO
      DO K = KEND-1, KSTR+KZ, -1
         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
            W(IJ, K) = W(IJ, K+1)
     &               - (  (FTX(IJ+LE, K) - FTX(IJ, K)) * RX
     &                  + (FTY(IJ+LN, K) - FTY(IJ, K)) * RY(IJ)) *
     &                 RXT(IJ) * RYT(IJ) * 0.5D0
         END DO
      END DO
      K = KSTR+KZ-1
      DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
         W(IJ, K) = W(IJ, K+1) * RHZBOT(IJ)
     &            - (  (  (FTX(IJ+LE, K) - FTX(IJ, K)) * RX
     &                  + (FTY(IJ+LN, K) - FTY(IJ, K)) * RY(IJ)) *
     &                 RXT(IJ) * RYT(IJ) * 0.5D0
     &               + DHDT(IJ) * DS(K)) * RHZBOT(IJ)
      END DO
      DO K = KSTR+KZ-2, KSTR+1, -1
         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
            W(IJ, K) = W(IJ, K+1)
     &               - (  (  (FTX(IJ+LE, K) - FTX(IJ, K)) * RX
     &                     + (FTY(IJ+LN, K) - FTY(IJ, K)) * RY(IJ)) *
     &                    RXT(IJ) * RYT(IJ) * 0.5D0
     &                  + DHDT(IJ) * DS(K)) * RHZBOT(IJ)
         END DO
      END DO
      DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
         W(IJ, KSTR) = 0.D0
      END DO

      RETURN
      END
