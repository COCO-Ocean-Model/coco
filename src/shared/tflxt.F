      SUBROUTINE FLXTRC(
     O                ADT,  DIFFZ,
     I                 TX,     HX,  
     I                 TY,     HZ,
     I                 UY,     VY,
     I                  W,    AHV)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the advection and diffusion terms of the tracer equations.
c
c  HISTORY
c     '03.04.23  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      HX(NXYDIM),     HZ(NXYDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      LOGICAL OFIRST
      REAL*8     ALP,   ALPV
      SAVE ALP, ALPV, OFIRST

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      COMMON /WORK/ FTX, FTY, FTZ

      REAL*8   HZBOT(NXYDIM)
      REAL*8       U,      V,     WT
      INTEGER     IJ,      K,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8     AHH,  ALPHA, ALPHAV
      SAVE AHH
      NAMELIST /NMDIFH/ AHH
      NAMELIST /NMWUPC/ ALPHA, ALPHAV

      DATA OFIRST / .TRUE. /
      DATA AHH / 0.D0 /
      DATA ALPHA, ALPHAV / 0.5D0, 0.5D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFH, END=98)
 98      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMWUPC, END=97)
 97      CONTINUE
         WRITE(JFPAR, NMDIFH)
         WRITE(JFPAR, NMWUPC)

         ALP  = (ALPHA - 0.5D+0) * 0.5D+0
         ALPV = ALPHAV - 0.5D+0
      END IF

      DO 120 N = 1, NTDIM
         DO 110 K = 1, NZDIM
            DO 100 IJ = 1, NXYDIM
               FTX(IJ, K, N) = 0.D0
               FTY(IJ, K, N) = 0.D0
               FTZ(IJ, K, N) = 0.D0
               ADT(IJ, K, N) = 0.D0
  100        CONTINUE
  110      CONTINUE
  120 CONTINUE
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DIFFZ(IJ, K) = 0.D0
         END DO
      END DO
      DO IJ = 1, NXYDIM
         HZBOT(IJ) = HZ(IJ) + ZBOT
      END DO


      DO K = KSTR+1, KSTR+KZ-1
         DO IJ = IJTSTR, IJTEND
            DIFFZ(IJ, K) = AHV(IJ, K) / DSM(K) / HZBOT(IJ) *
     &                     AMFTZ(IJ, K)
         END DO
      END DO
      DO 220 N = 1, NTDIM
         DO 210 K = KSTR+1, KSTR+KZ-1
            DO 200 IJ = IJTSTR, IJTEND
               WT = SIGN(ALPV, W(IJ, K))
               FTZ(IJ, K, N) = 
     &         (  DIFFZ(IJ, K) * (TX(IJ, K-1, N) - TX(IJ, K, N))
     &          - W(IJ, K) * HZBOT(IJ) *
     &            (  (0.5D0 - WT) * TX(IJ, K-1, N)
     &             + (0.5D0 + WT) * TX(IJ, K  , N))) *
     &         AMFTZ(IJ, K)
  200       CONTINUE
  210    CONTINUE
  220 CONTINUE

      DO K = KSTR+KZ, KEND
         DO IJ = IJTSTR, IJTEND
            DIFFZ(IJ, K) = AHV(IJ, K) / DZM(IJ, K) * AMFTZ(IJ, K)
         ENDDO
      ENDDO
      DO 250 N = 1, NTDIM
         DO 240 K = KSTR+KZ, KEND
            DO 230 IJ = IJTSTR, IJTEND
               WT = SIGN(ALPV, W(IJ, K))
               FTZ(IJ, K, N) = 
     &         (  DIFFZ(IJ, K) * (TX(IJ, K-1, N) - TX(IJ, K, N))
     &          - W(IJ, K) * (  (0.5D0 - WT) * TX(IJ, K-1, N)
     &                        + (0.5D0 + WT) * TX(IJ, K,   N))) *
     &         AMFTZ(IJ, K)
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE    

      DO 610 N = 1, NTDIM
         DO 600 K = KSTR, KEND

            DO 300 IJ = IJTSTR, IJTEND+NXDIM
               V = VY(IJ+LS , K) * DZV(IJ+LS , K) * HXU(IJ+LS)
     &           + VY(IJ+LSW, K) * DZV(IJ+LSW, K) * HXU(IJ+LSW)
               WT = SIGN(ALP, V)
               FTY(IJ, K, N) =
     &            AHH * (TX(IJ, K, N) - TX(IJ+LS, K, N)) * RYM(IJ+LS) *
     &            (HXU(IJ+LS) + HXU(IJ+LSW)) / (HYT(IJ) + HYT(IJ+LS)) *
     &            AMFTY(IJ, K)
     &          - V * (  (0.25D0 - WT) * TX(IJ, K,    N)
     &                 + (0.25D0 + WT) * TX(IJ+LS, K, N))
 300        CONTINUE

            DO 400 IJ = IJTSTR, IJTEND+1
               U = UY(IJ+LW , K) * DZV(IJ+LW , K) * HYU(IJ+LW)
     &           + UY(IJ+LSW, K) * DZV(IJ+LSW, K) * HYU(IJ+LSW)
               WT = SIGN(ALP, U)
               FTX(IJ, K, N) =
     &           AHH * (TX(IJ, K, N) - TX(IJ+LW, K, N)) * RX *
     &           (HYU(IJ+LW) + HYU(IJ+LSW)) / (HXT(IJ) + HXT(IJ+LW)) *
     &           AMFTX(IJ, K)
     &         - U * (  (0.25D0 - WT) * TX(IJ, K, N)
     &                + (0.25D0 + WT) * TX(IJ+LW, K, N))
 400        CONTINUE

            DO 500 IJ = IJTSTR, IJTEND
               ADT(IJ, K, N) =
     &         (  (  (FTX(IJ+LE, K, N) - FTX(IJ, K, N)) * RX
     &             + (FTY(IJ+LN, K, N) - FTY(IJ, K, N)) * RY(IJ)) *
     &            RXT(IJ) * RYT(IJ)
     &          + FTZ(IJ, K, N) - FTZ(IJ, K+1, N)) / DZ(IJ, K)
 500        CONTINUE

 600     CONTINUE
 610  CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE FLXTRB(
     O                ADT,  DIFFZ,
     I                 TX, 
     I                 TY, 
     I                 UY,     VY,      W,
     I                AHV)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the advection and diffusion terms of the tracer equations
c for BBL.
c
c  HISTORY
c     '01.02.13  H.Hasumi
c     '01.05.02  H.Hasumi
c     '02.05.29  H.Nakano: tracer dimension
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      REAL*8     ALP,   ALPV
      LOGICAL OFIRST
      SAVE ALP, ALPV, OFIRST

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      COMMON /WORK/ FTX, FTY, FTZ

      REAL*8       U,      V,      WT
      INTEGER     IJ,      K,     KUP,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8  AHHBBL,  ALPHA, ALPHAV
      SAVE AHHBBL
      NAMELIST /NMBBDH/ AHHBBL
      NAMELIST /NMWUPC/ ALPHA, ALPHAV

      DATA OFIRST / .TRUE. /
      DATA AHHBBL / 0.D0 /
      DATA ALPHA, ALPHAV / 0.5D0, 0.5D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBBDH, END=98)
 98      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMWUPC, END=97)
 97      CONTINUE
         WRITE(JFPAR, NMBBDH)

         ALP  = (ALPHA - 0.5D0) * 0.5D0
         ALPV = ALPHAV - 0.5D0
      END IF

      DO N = 1, NTDIM

         DO IJ = IJTSTR, IJTEND
            K = NBOT(IJ)
            KUP = MAX(K-1, 1)
            WT = SIGN(ALPV, W(IJ, KEND))
            FTZ(IJ, KEND, N) =
     &      (  DIFFZ(IJ, K) * (TX(IJ, KUP, N) - TX(IJ, KEND, N))
     &       - W(IJ, KEND) * (  (0.5D0 - WT) * TX(IJ, KUP, N)
     &                        + (0.5D0 + WT) * TX(IJ, KEND, N))) *
     &      AMSKTB(IJ)
         END DO

         DO IJ = IJTSTR, IJTEND+NXDIM
            V = VY(IJ+LS , KEND) * DZV(IJ+LS , KEND) * HXU(IJ+LS)
     &        + VY(IJ+LSW, KEND) * DZV(IJ+LSW, KEND) * HXU(IJ+LSW)
            WT = SIGN(ALP, V)
            FTY(IJ, KEND, N) =
     &         AHHBBL * (TX(IJ, KEND, N) - TX(IJ+LS, KEND, N)) *
     &         RYM(IJ) * 
     &         (HXU(IJ+LS) + HXU(IJ+LSW)) / (HYT(IJ) + HYT(IJ+LS)) *
     &         AMFTY(IJ, KEND)
     &       - V * (  (0.25D0 - WT) * TX(IJ, KEND, N)
     &              + (0.25D0 + WT) * TX(IJ+LS, KEND, N))
         END DO

         DO IJ = IJTSTR, IJTEND+1
            U = UY(IJ+LW , KEND) * DZV(IJ+LW , KEND) * HYU(IJ+LW)
     &        + UY(IJ+LSW, KEND) * DZV(IJ+LSW, KEND) * HYU(IJ+LSW)
            WT = SIGN(ALP, U)
            FTX(IJ, KEND, N) = 
     &        AHHBBL * (TX(IJ, KEND, N) - TX(IJ+LW, KEND, N)) * RX *
     &        (HYU(IJ+LW) + HYU(IJ+LSW)) / (HXT(IJ) + HXT(IJ+LW)) *
     &        AMFTX(IJ, KEND)
     &      - U * (  (0.25D0 - WT) * TX(IJ, KEND, N)
     &             + (0.25D0 + WT) * TX(IJ+LW, KEND, N))
         END DO

         DO IJ = IJTSTR, IJTEND
            ADT(IJ, KEND, N) =
     &      (  (  (  FTX(IJ+LE, KEND, N)
     &             - FTX(IJ, KEND, N)) * RX
     &          + (  FTY(IJ+LN, KEND, N)
     &             - FTY(IJ, KEND, N)) * RY(IJ)) * RXT(IJ) * RYT(IJ)
     &       + FTZ(IJ, KEND, N)) / DZ(IJ, KEND)
         END DO

      END DO

      RETURN
      END
#endif
