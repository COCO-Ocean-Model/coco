      SUBROUTINE SHDIFF(
     M                 TX,     HX)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '03.04.24  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      HX(NXYDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)

      REAL*8     TSH(NXYDIM, NZDIM, NTDIM)
      REAL*8     FHX(NXYDIM),    FHY(NXYDIM)
      REAL*8     FTX(NXYDIM),    FTY(NXYDIM)
      REAL*8  RHXBOT(NXYDIM)
c      COMMON /WORK/ TSH, FHX, FHY, FTX, FTY, RHXBOT

      INTEGER     IJ,      K,      N 
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8     ASH
      SAVE ASH
      NAMELIST /NMDFSH/ ASH

      DATA OFIRST / .TRUE. /
      DATA ASH / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDFSH, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMDFSH)
      END IF

      DO N = 1, NTDIM
         DO K = KSTR, KSTR+KZ-1
            DO IJ = 1, NXYDIM
               TSH(IJ, K, N) = TX(IJ, K, N) * (HX(IJ) + ZBOT)
            END DO
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND+NXDIM
         FHX(IJ) = ASH * (HX(IJ) - HX(IJ+LW)) * RX *
     &             (HYU(IJ+LW) + HYU(IJ+LSW)) / (HXT(IJ) + HXT(IJ+LW)) *
     &             AMSKT(IJ, KSTR) * AMSKT(IJ+LW, KSTR)
         FHY(IJ) = ASH * (HX(IJ) - HX(IJ+LS)) * RYM(IJ+LS) *
     &             (HXU(IJ+LS) + HXU(IJ+LSW)) / (HYT(IJ) + HYT(IJ+LS)) *
     &             AMSKT(IJ, KSTR) * AMSKT(IJ+LS, KSTR)
      END DO
      DO IJ = IJTSTR, IJTEND
         HX(IJ) = HX(IJ)
     &          + TS * (  (FHX(IJ+LE) - FHX(IJ)) * RX
     &                  + (FHY(IJ+LN) - FHY(IJ)) * RY(IJ)) *
     &            RXT(IJ) * RYT(IJ) * AMSKT(IJ, KSTR)
      END DO
      DO IJ = IJTSTR, IJTEND
         RHXBOT(IJ) = 1.D0 / (HX(IJ) + ZBOT)
      END DO

      DO N = 1, NTDIM
         DO K = KSTR, KSTR+KZ-1
            DO IJ = IJTSTR, IJTEND+NXDIM
               FTX(IJ) =
     &         (  (FHX(IJ) + ABS(FHX(IJ))) * TX(IJ, K, N)
     &          + (FHX(IJ) - ABS(FHX(IJ))) * TX(IJ+LW, K, N)) *
     &         0.5D0 * AMSKT(IJ, K) * AMSKT(IJ+LW, K)
               FTY(IJ) =
     &         (  (FHY(IJ) + ABS(FHY(IJ))) * TX(IJ, K, N)
     &          + (FHY(IJ) - ABS(FHY(IJ))) * TX(IJ+LS, K, N)) *
     &         0.5D0 * AMSKT(IJ, K) * AMSKT(IJ+LS, K)
            END DO
            DO IJ = IJTSTR, IJTEND
               TSH(IJ, K, N) =
     &         TSH(IJ, K, N)
     &         + TS * (  (FTX(IJ+LE) - FTX(IJ)) * RX
     &                 + (FTY(IJ+LN) - FTY(IJ)) * RY(IJ)) *
     &           RXT(IJ) * RYT(IJ) * AMSKT(IJ, K)
            END DO
            DO IJ = IJTSTR, IJTEND
               TX(IJ, K, N) = TSH(IJ, K, N) * RHXBOT(IJ)
            END DO
         END DO
      END DO

      RETURN
      END
