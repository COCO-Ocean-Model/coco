      SUBROUTINE ADVVEL(
     O                FUX,    FUY,   FUNE,   FUSE,
     O                FVX,    FVY,   FVNE,   FVSE,
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I                 HX,     HY,
     I               UADV,   VADV,   WADV)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Advection and metric terms of the equation of motion.
c
c  HISTORY
c     '03.05.13  H.Hasumi: from COCO3.4
c     '06.04.19  T.Suzuki:
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '10.04.14  M.kurogi
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HY(NXYDIM),            HX(NXYDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM, 9)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8    FUNE(NXYDIM, NZDIM),   FVNE(NXYDIM, NZDIM)
      REAL*8    FUSE(NXYDIM, NZDIM),   FVSE(NXYDIM, NZDIM)

      REAL*8     CXN(NXYDIM, NZDIM),    CXS(NXYDIM, NZDIM)
      REAL*8     CYE(NXYDIM, NZDIM),    CYW(NXYDIM, NZDIM)
      REAL*8     CNE(NXYDIM, NZDIM),    CSE(NXYDIM, NZDIM)
      SAVE CXN, CXS, CYE, CYW, CNE, CSE

      REAL*8     FUZ(NXYDIM, NZDIM),    FVZ(NXYDIM, NZDIM)
      REAL*8    FUZU(NXYDIM, NZDIM),   FVZU(NXYDIM, NZDIM)
      REAL*8    FUZD(NXYDIM, NZDIM),   FVZD(NXYDIM, NZDIM)
      REAL*8      RZ(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
      COMMON /WORK/ FUZ, FUZU, FUZD, FVZ, FVZU, FVZD, RZ, RZM

      REAL*8     DIV(NXYDIM, KSTR:KSTR+KZ-1)
      REAL*8   HVBOT(NXYDIM)
      INTEGER     IJ,      K

      LOGICAL OFIRST
      SAVE OFIRST

      DATA OFIRST / .TRUE. /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
#ifdef OPT_BBL
         CALL RMMSKV
#endif
         DO 31 K = 1, NZDIM
            DO 30 IJ = 1, NXYDIM
               CXN(IJ, K) = 0.D0
               CXS(IJ, K) = 0.D0
               CYE(IJ, K) = 0.D0
               CYW(IJ, K) = 0.D0
               CNE(IJ, K) = 0.D0
               CSE(IJ, K) = 0.D0
 30         CONTINUE
 31      CONTINUE
         DO 41 K = KSTR, KEND
            DO 40 IJ = IJVSTR-NXDIM-1, IJVEND+NXDIM+1
            CXN(IJ, K) = AMSKV(IJ, K) * AMSKV(IJ+LW, K) *
     &                   (3.D0 - AMSKV(IJ+LS, K) - AMSKV(IJ+LSW, K)
     &                         + AMSKV(IJ+LS, K) * AMSKV(IJ+LSW, K))
            CXS(IJ, K) = AMSKV(IJ+LS, K) * AMSKV(IJ+LSW, K) *
     &                   (3.D0 - AMSKV(IJ, K) - AMSKV(IJ+LW, K)
     &                         + AMSKV(IJ, K) * AMSKV(IJ+LW, K))
            CYE(IJ, K) = AMSKV(IJ, K) * AMSKV(IJ+LS, K) *
     &                   (3.D0 - AMSKV(IJ+LW, K) - AMSKV(IJ+LSW, K)
     &                         + AMSKV(IJ+LW, K) * AMSKV(IJ+LSW, K))
            CYW(IJ, K) = AMSKV(IJ+LW, K) * AMSKV(IJ+LSW, K) *
     &                   (3.D0 - AMSKV(IJ, K) - AMSKV(IJ+LS, K)
     &                         + AMSKV(IJ, K) * AMSKV(IJ+LS, K))
            CNE(IJ, K) = AMSKV(IJ, K) * AMSKV(IJ+LSW, K) *
     &                   (3.D0 - AMSKV(IJ+LW, K) - AMSKV(IJ+LS, K))
            CSE(IJ, K) = AMSKV(IJ+LW, K) * AMSKV(IJ+LS, K) *
     &                   (3.D0 - AMSKV(IJ, K) - AMSKV(IJ+LSW, K))
 40         CONTINUE
 41      CONTINUE
C#ifdef OPT_BBL
C         CALL ADMKV1
C#endif
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            FUX (IJ, K) = 0.D0
            FUY (IJ, K) = 0.D0
            FUNE(IJ, K) = 0.D0
            FUSE(IJ, K) = 0.D0
            FVX (IJ, K) = 0.D0
            FVY (IJ, K) = 0.D0
            FVNE(IJ, K) = 0.D0
            FVSE(IJ, K) = 0.D0
            FUZ (IJ, K) = 0.D0
            FUZU(IJ, K) = 0.D0
            FUZD(IJ, K) = 0.D0
            FVZ (IJ, K) = 0.D0
            FVZU(IJ, K) = 0.D0
            FVZD(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO IJ = IJVSTR, IJVEND
         HVBOT(IJ) = (  (HY(IJ)    + HY(IJ+LE) ) * DY(IJ)
     &                + (HY(IJ+LN) + HY(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
         HVBOT(IJ) = 1.D0 / HVBOT(IJ)
      END DO

      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND+NXDIM+1
            FVY(IJ, K) = (  CYE(IJ, K) * VADV(IJ, K) * HXT(IJ)
     &                    + CYW(IJ+LE, K) * VADV(IJ+LE, K) * HXT(IJ+LE))
     &                   / 6.D0
            FVX(IJ, K) = (  CXN(IJ, K) * UADV(IJ, K) * HYT(IJ)
     &                    + CXS(IJ+LN, K) * UADV(IJ+LN, K) * HYT(IJ+LN))
     &                   / 6.D0
            FVNE(IJ, K) = CNE(IJ, K) / 6.D0 *
     &                    (  UADV(IJ, K) * DY(IJ) * HYT(IJ)
     &                     + VADV(IJ, K) * DX * HXT(IJ))
            FVSE(IJ, K) = CSE(IJ, K) / 6.D0 *
     &                    (  UADV(IJ, K) * DY(IJ) * HYT(IJ)
     &                     - VADV(IJ, K) * DX * HXT(IJ))
         END DO
      END DO

      DO K = KSTR+1, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            FVZ(IJ, K) = WADV(IJ, K, 1)
            FVZU(IJ, K) = WADV(IJ, K, 2)
     &                  + WADV(IJ, K, 3)
     &                  + WADV(IJ, K, 4)
     &                  + WADV(IJ, K, 5)
     &                  + WADV(IJ, K, 6)
     &                  + WADV(IJ, K, 7)
     &                  + WADV(IJ, K, 8)
     &                  + WADV(IJ, K, 9)
            FVZD(IJ, K) = WADV(IJ+LN , K+1, 6)
     &                  + WADV(IJ+LNE, K+1, 7)
     &                  + WADV(IJ+LE , K+1, 8)
     &                  + WADV(IJ+LSE, K+1, 9)
     &                  + WADV(IJ+LS , K+1, 2)
     &                  + WADV(IJ+LSW, K+1, 3)
     &                  + WADV(IJ+LW , K+1, 4)
     &                  + WADV(IJ+LNW, K+1, 5)
         END DO
      END DO
      DO IJ = IJVSTR, IJVEND
         FVZ(IJ, KSTR+KZ) = WADV(IJ, KSTR+KZ, 1) * HVBOT(IJ)
         FVZD(IJ, KSTR) = WADV(IJ+LN , KSTR+1, 6)
     &                  + WADV(IJ+LNE, KSTR+1, 7)
     &                  + WADV(IJ+LE , KSTR+1, 8)
     &                  + WADV(IJ+LSE, KSTR+1, 9)
     &                  + WADV(IJ+LS , KSTR+1, 2)
     &                  + WADV(IJ+LSW, KSTR+1, 3)
     &                  + WADV(IJ+LW , KSTR+1, 4)
     &                  + WADV(IJ+LNW, KSTR+1, 5)
         FVZD(IJ, KSTR+KZ-1) = (  WADV(IJ+LN , KSTR+KZ, 6)
     &                          + WADV(IJ+LNE, KSTR+KZ, 7)
     &                          + WADV(IJ+LE , KSTR+KZ, 8)
     &                          + WADV(IJ+LSE, KSTR+KZ, 9)
     &                          + WADV(IJ+LS , KSTR+KZ, 2)
     &                          + WADV(IJ+LSW, KSTR+KZ, 3)
     &                          + WADV(IJ+LW , KSTR+KZ, 4)
     &                          + WADV(IJ+LNW, KSTR+KZ, 5)) * HVBOT(IJ)
      END DO

      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            DIV(IJ, K) = (  (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                       + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                      RXU(IJ) * RYU(IJ)
     &                    + (  FVNE(IJ+LNE, K) - FVNE(IJ, K)
     &                       + FVSE(IJ+LE, K) - FVSE(IJ+LN, K)) *
     &                      RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                    + (FVZ(IJ, K) - FVZ(IJ, K+1)) * RS(K)
     &                    + (FVZU(IJ, K) - FVZD(IJ, K)) * RS(K)) *
     &                   AMSKV(IJ, K)
         END DO
      END DO

      DO 1100 IJ = IJVSTR, IJVEND

         FUZD(IJ, KSTR) = - 0.5D0 *
     &   (  WADV(IJ+LN , KSTR+1, 6) *
     &      (UY(IJ+LN , KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LNE, KSTR+1, 7) *
     &      (UY(IJ+LNE, KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LE , KSTR+1, 8) *
     &      (UY(IJ+LE , KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LSE, KSTR+1, 9) *
     &      (UY(IJ+LSE, KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LS , KSTR+1, 2) *
     &      (UY(IJ+LS , KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LSW, KSTR+1, 3) *
     &      (UY(IJ+LSW, KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LW , KSTR+1, 4) *
     &      (UY(IJ+LW , KSTR+1) + UY(IJ, KSTR))
     &    + WADV(IJ+LNW, KSTR+1, 5) *
     &      (UY(IJ+LNW, KSTR+1) + UY(IJ, KSTR)))

         FVZD(IJ, KSTR) = - 0.5D0 *
     &   (  WADV(IJ+LN , KSTR+1, 6) *
     &      (VY(IJ+LN , KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LNE, KSTR+1, 7) *
     &      (VY(IJ+LNE, KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LE , KSTR+1, 8) *
     &      (VY(IJ+LE , KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LSE, KSTR+1, 9) *
     &      (VY(IJ+LSE, KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LS , KSTR+1, 2) *
     &      (VY(IJ+LS , KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LSW, KSTR+1, 3) *
     &      (VY(IJ+LSW, KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LW , KSTR+1, 4) *
     &      (VY(IJ+LW , KSTR+1) + VY(IJ, KSTR))
     &    + WADV(IJ+LNW, KSTR+1, 5) *
     &      (VY(IJ+LNW, KSTR+1) + VY(IJ, KSTR)))
 1100 CONTINUE

      DO 1210 K = KSTR+1, KEND
         DO 1200 IJ = IJVSTR, IJVEND
            FUZ(IJ, K) = 
     &                 - WADV(IJ, K, 1) * 0.5D0 *
     &                   (UY(IJ, K-1) + UY(IJ, K))
            FVZ(IJ, K) = 
     &                 - WADV(IJ, K, 1) * 0.5D0 *
     &                   (VY(IJ, K-1) + VY(IJ, K))

            FUZU(IJ, K) = - 0.5D0 *
     &      (  WADV(IJ, K, 2) * (UY(IJ, K) + UY(IJ+LN , K-1))
     &       + WADV(IJ, K, 3) * (UY(IJ, K) + UY(IJ+LNE, K-1))
     &       + WADV(IJ, K, 4) * (UY(IJ, K) + UY(IJ+LE , K-1))
     &       + WADV(IJ, K, 5) * (UY(IJ, K) + UY(IJ+LSE, K-1))
     &       + WADV(IJ, K, 6) * (UY(IJ, K) + UY(IJ+LS , K-1))
     &       + WADV(IJ, K, 7) * (UY(IJ, K) + UY(IJ+LSW, K-1))
     &       + WADV(IJ, K, 8) * (UY(IJ, K) + UY(IJ+LW , K-1))
     &       + WADV(IJ, K, 9) * (UY(IJ, K) + UY(IJ+LNW, K-1)))
            FUZD(IJ, K) = - 0.5D0 *
     &      (  WADV(IJ+LN , K+1, 6) * (UY(IJ+LN , K+1) + UY(IJ, K))
     &       + WADV(IJ+LNE, K+1, 7) * (UY(IJ+LNE, K+1) + UY(IJ, K))
     &       + WADV(IJ+LE , K+1, 8) * (UY(IJ+LE , K+1) + UY(IJ, K))
     &       + WADV(IJ+LSE, K+1, 9) * (UY(IJ+LSE, K+1) + UY(IJ, K))
     &       + WADV(IJ+LS , K+1, 2) * (UY(IJ+LS , K+1) + UY(IJ, K))
     &       + WADV(IJ+LSW, K+1, 3) * (UY(IJ+LSW, K+1) + UY(IJ, K))
     &       + WADV(IJ+LW , K+1, 4) * (UY(IJ+LW , K+1) + UY(IJ, K))
     &       + WADV(IJ+LNW, K+1, 5) * (UY(IJ+LNW, K+1) + UY(IJ, K)))

            FVZU(IJ, K) = - 0.5D0 *
     &      (  WADV(IJ, K, 2) * (VY(IJ, K) + VY(IJ+LN , K-1))
     &       + WADV(IJ, K, 3) * (VY(IJ, K) + VY(IJ+LNE, K-1))
     &       + WADV(IJ, K, 4) * (VY(IJ, K) + VY(IJ+LE , K-1))
     &       + WADV(IJ, K, 5) * (VY(IJ, K) + VY(IJ+LSE, K-1))
     &       + WADV(IJ, K, 6) * (VY(IJ, K) + VY(IJ+LS , K-1))
     &       + WADV(IJ, K, 7) * (VY(IJ, K) + VY(IJ+LSW, K-1))
     &       + WADV(IJ, K, 8) * (VY(IJ, K) + VY(IJ+LW , K-1))
     &       + WADV(IJ, K, 9) * (VY(IJ, K) + VY(IJ+LNW, K-1)))
            FVZD(IJ, K) = - 0.5D0 *
     &      (  WADV(IJ+LN , K+1, 6) * (VY(IJ+LN , K+1) + VY(IJ, K))
     &       + WADV(IJ+LNE, K+1, 7) * (VY(IJ+LNE, K+1) + VY(IJ, K))
     &       + WADV(IJ+LE , K+1, 8) * (VY(IJ+LE , K+1) + VY(IJ, K))
     &       + WADV(IJ+LSE, K+1, 9) * (VY(IJ+LSE, K+1) + VY(IJ, K))
     &       + WADV(IJ+LS , K+1, 2) * (VY(IJ+LS , K+1) + VY(IJ, K))
     &       + WADV(IJ+LSW, K+1, 3) * (VY(IJ+LSW, K+1) + VY(IJ, K))
     &       + WADV(IJ+LW , K+1, 4) * (VY(IJ+LW , K+1) + VY(IJ, K))
     &       + WADV(IJ+LNW, K+1, 5) * (VY(IJ+LNW, K+1) + VY(IJ, K)))
 1200    CONTINUE
 1210 CONTINUE

#ifdef OPT_BBL
      DO IJ = IJVSTR, IJVEND
         FUZD(IJ,KEND-1) = 0.D0
         FVZD(IJ,KEND-1) = 0.D0
         K = NBOTV(IJ) - AMSKVB(IJ)
            FUZD(IJ, K) = FUZD(IJ, K) - 0.5D0 *
     &      (  WADV(IJ+LN , KEND, 6) * (UY(IJ+LN , KEND) + UY(IJ, K))
     &       + WADV(IJ+LNE, KEND, 7) * (UY(IJ+LNE, KEND) + UY(IJ, K))
     &       + WADV(IJ+LE , KEND, 8) * (UY(IJ+LE , KEND) + UY(IJ, K))
     &       + WADV(IJ+LSE, KEND, 9) * (UY(IJ+LSE, KEND) + UY(IJ, K))
     &       + WADV(IJ+LS , KEND, 2) * (UY(IJ+LS , KEND) + UY(IJ, K))
     &       + WADV(IJ+LSW, KEND, 3) * (UY(IJ+LSW, KEND) + UY(IJ, K))
     &       + WADV(IJ+LW , KEND, 4) * (UY(IJ+LW , KEND) + UY(IJ, K))
     &       + WADV(IJ+LNW, KEND, 5) * (UY(IJ+LNW, KEND) + UY(IJ, K)))
            FVZD(IJ, K) = FVZD(IJ, K) - 0.5D0 *
     &      (  WADV(IJ+LN , KEND, 6) * (VY(IJ+LN , KEND) + VY(IJ, K))
     &       + WADV(IJ+LNE, KEND, 7) * (VY(IJ+LNE, KEND) + VY(IJ, K))
     &       + WADV(IJ+LE , KEND, 8) * (VY(IJ+LE , KEND) + VY(IJ, K))
     &       + WADV(IJ+LSE, KEND, 9) * (VY(IJ+LSE, KEND) + VY(IJ, K))
     &       + WADV(IJ+LS , KEND, 2) * (VY(IJ+LS , KEND) + VY(IJ, K))
     &       + WADV(IJ+LSW, KEND, 3) * (VY(IJ+LSW, KEND) + VY(IJ, K))
     &       + WADV(IJ+LW , KEND, 4) * (VY(IJ+LW , KEND) + VY(IJ, K))
     &       + WADV(IJ+LNW, KEND, 5) * (VY(IJ+LNW, KEND) + VY(IJ, K)))
      ENDDO
#endif

      DO 1310 K = KSTR, KEND
         DO 1300 IJ = IJVSTR, IJVEND+NXDIM
            FUY(IJ, K) =
     &                 - (  CYE(IJ, K) * VADV(IJ, K) * HXT(IJ)
     &                    + CYW(IJ+LE, K) * VADV(IJ+LE, K) * HXT(IJ+LE))
     &                   / 6.D0 * (UY(IJ, K) + UY(IJ+LS, K)) * 0.5D0
            FVY(IJ, K) =
     &                 - (  CYE(IJ, K) * VADV(IJ, K) * HXT(IJ)
     &                    + CYW(IJ+LE, K) * VADV(IJ+LE, K) * HXT(IJ+LE))
     &                   / 6.D0 * (VY(IJ, K) + VY(IJ+LS, K)) * 0.5D0
 1300    CONTINUE
 1310 CONTINUE

      DO 1410 K = KSTR, KEND
         DO 1400 IJ = IJVSTR, IJVEND+1
            FUX(IJ, K) =
     &                 - (  CXN(IJ, K) * UADV(IJ, K) * HYT(IJ)
     &                    + CXS(IJ+LN, K) * UADV(IJ+LN, K) * HYT(IJ+LN))
     &                   / 6.D0 * (UY(IJ, K) + UY(IJ+LW, K)) * 0.5D0
            FVX(IJ, K) =
     &                 - (  CXN(IJ, K) * UADV(IJ, K) * HYT(IJ)
     &                    + CXS(IJ+LN, K) * UADV(IJ+LN, K) * HYT(IJ+LN))
     &                   / 6.D0 * (VY(IJ, K) + VY(IJ+LW, K)) * 0.5D0
 1400    CONTINUE
 1410 CONTINUE

      DO 1510 K = KSTR, KEND
         DO 1500 IJ = IJVSTR, IJVEND+NXDIM+1
            FUNE(IJ, K) =
     &                  - CNE(IJ, K) / 6.D0 *
     &                    (  UADV(IJ, K) * DY(IJ) * HYT(IJ)
     &                     + VADV(IJ, K) * DX * HXT(IJ)) *
     &                    (UY(IJ, K) + UY(IJ+LSW, K)) * 0.5D0
            FVNE(IJ, K) =
     &                  - CNE(IJ, K) / 6.D0 *
     &                    (  UADV(IJ, K) * DY(IJ) * HYT(IJ)
     &                     + VADV(IJ, K) * DX * HXT(IJ)) *
     &                    (VY(IJ, K) + VY(IJ+LSW, K)) * 0.5D0
 1500    CONTINUE
 1510 CONTINUE

      DO 1610 K = KSTR, KEND
         DO 1600 IJ = IJVSTR, IJVEND+NXDIM
            FUSE(IJ, K) =
     &                  - CSE(IJ, K) / 6.D0 *
     &                    (  UADV(IJ, K) * DY(IJ) * HYT(IJ)
     &                     - VADV(IJ, K) * DX * HXT(IJ)) *
     &                    (UY(IJ+LW, K) + UY(IJ+LS, K)) * 0.5D0
            FVSE(IJ, K) =
     &                  - CSE(IJ, K) / 6.D0 *
     &                    (  UADV(IJ, K) * DY(IJ) * HYT(IJ)
     &                     - VADV(IJ, K) * DX * HXT(IJ)) *
     &                    (VY(IJ+LW, K) + VY(IJ+LS, K)) * 0.5D0
 1600    CONTINUE
 1610 CONTINUE

      DO 2010 K = KSTR, KSTR+KZ-2
         DO 2000 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FUNE(IJ+LNE, K) - FUNE(IJ   , K)
     &                      + FUSE(IJ+LE , K) - FUSE(IJ+LN, K)) *
     &                     RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                   + (FUZ(IJ, K) - FUZ(IJ, K+1)) * RS(K)
     &                   + (FUZU(IJ, K) - FUZD(IJ, K)) * RS(K)
     &                   + UY(IJ, K) * DIV(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FVNE(IJ+LNE, K) - FVNE(IJ   , K)
     &                      + FVSE(IJ+LE , K) - FVSE(IJ+LN, K)) *
     &                     RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                   + (FVZ(IJ, K) - FVZ(IJ, K+1)) * RS(K)
     &                   + (FVZU(IJ, K) - FVZD(IJ, K)) * RS(K)
     &                   + VY(IJ, K) * DIV(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 2000    CONTINUE
 2010 CONTINUE
      K = KSTR+KZ-1
         DO 2020 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FUNE(IJ+LNE, K) - FUNE(IJ   , K)
     &                      + FUSE(IJ+LE , K) - FUSE(IJ+LN, K)) *
     &                     RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                   + (  FUZ(IJ, K)
     &                      - FUZ(IJ, K+1) * HVBOT(IJ)) * RS(K)
     &                   + (  FUZU(IJ, K)
     &                      - FUZD(IJ, K) * HVBOT(IJ)) * RS(K)
     &                   + UY(IJ, K) * DIV(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FVNE(IJ+LNE, K) - FVNE(IJ   , K)
     &                      + FVSE(IJ+LE , K) - FVSE(IJ+LN, K)) *
     &                     RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                   + (  FVZ(IJ, K)
     &                      - FVZ(IJ, K+1) * HVBOT(IJ)) * RS(K)
     &                   + (  FVZU(IJ, K)
     &                      - FVZD(IJ, K) * HVBOT(IJ)) * RS(K)
     &                   + VY(IJ, K) * DIV(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 2020    CONTINUE
      DO 2040 K = KSTR+KZ, KEND
         DO 2030 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (  (  FUX(IJ+LE, K)
     &                            - FUX(IJ, K)) * RX
     &                         + (  FUY(IJ+LN, K)
     &                            - FUY(IJ, K)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + (  FUNE(IJ+LNE, K) - FUNE(IJ   , K)
     &                         + FUSE(IJ+LE , K) - FUSE(IJ+LN, K)) *
     &                        RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                      + FUZ(IJ, K) - FUZ(IJ, K+1)
     &                      + FUZU(IJ, K) - FUZD(IJ, K)) * RZ(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (   (  (  FVX(IJ+LE, K)
     &                             - FVX(IJ, K)) * RX
     &                          + (  FVY(IJ+LN, K)
     &                             - FVY(IJ, K)) * RYM(IJ)) *
     &                         RXU(IJ) * RYU(IJ)
     &                      + (  FVNE(IJ+LNE, K) - FVNE(IJ   , K)
     &                         + FVSE(IJ+LE , K) - FVSE(IJ+LN, K)) *
     &                        RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                      + FVZ(IJ, K) - FVZ(IJ, K+1)
     &                      + FVZU(IJ, K) - FVZD(IJ, K)) * RZ(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 2030    CONTINUE
 2040 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE ADVVLB(
     M                FUX,    FUY,   FUNE,   FUSE,
     M                FVX,    FVY,   FVNE,   FVSE,
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I               UADV,   VADV,   WADV)
      IMPLICIT NONE

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM, 9)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8    FUNE(NXYDIM, NZDIM),   FVNE(NXYDIM, NZDIM)
      REAL*8    FUSE(NXYDIM, NZDIM),   FVSE(NXYDIM, NZDIM)

      REAL*8      RZ(NXYDIM),    RZM(NXYDIM)
      SAVE RZ, RZM

      REAL*8     CXN(NXYDIM),    CXS(NXYDIM)
      REAL*8     CYE(NXYDIM),    CYW(NXYDIM)
      REAL*8     CNE(NXYDIM),    CSE(NXYDIM)
      SAVE CXN, CXS, CYE, CYW, CNE, CSE

      REAL*8     FUZ(NXYDIM, NZDIM),    FVZ(NXYDIM, NZDIM)
      REAL*8    FUZU(NXYDIM, NZDIM),   FVZU(NXYDIM, NZDIM)
      REAL*8    FUZD(NXYDIM, NZDIM),   FVZD(NXYDIM, NZDIM)
      COMMON /WORK/ FUZ, FUZU, FUZD, FVZ, FVZU, FVZD

      INTEGER     IJ,      K,    KUP

      LOGICAL OFIRST
      SAVE OFIRST

      DATA OFIRST / .TRUE. /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.

         DO IJ = 1, NXYDIM
            RZ   (IJ) = 1.D0 / DZV(IJ, KEND)
            RZM  (IJ) = 1.D0 / DZM(IJ, KEND)
         END DO

         DO 30 IJ = 1, NXYDIM
            CXN(IJ) = 0.D0
            CXS(IJ) = 0.D0
            CYE(IJ) = 0.D0
            CYW(IJ) = 0.D0
            CNE(IJ) = 0.D0
            CSE(IJ) = 0.D0
 30      CONTINUE
         DO 40 IJ = IJVSTR-NXDIM-1, IJVEND+NXDIM+1
         CXN(IJ) = AMSKVB(IJ) * AMSKVB(IJ+LW) *
     &             (3.D0 - AMSKVB(IJ+LS) - AMSKVB(IJ+LSW)
     &                   + AMSKVB(IJ+LS) * AMSKVB(IJ+LSW))
         CXS(IJ) = AMSKVB(IJ+LS) * AMSKVB(IJ+LSW) *
     &             (3.D0 - AMSKVB(IJ) - AMSKVB(IJ+LW)
     &                   + AMSKVB(IJ) * AMSKVB(IJ+LW))
         CYE(IJ) = AMSKVB(IJ) * AMSKVB(IJ+LS) *
     &             (3.D0 - AMSKVB(IJ+LW) - AMSKVB(IJ+LSW)
     &                   + AMSKVB(IJ+LW) * AMSKVB(IJ+LSW))
         CYW(IJ) = AMSKVB(IJ+LW) * AMSKVB(IJ+LSW) *
     &             (3.D0 - AMSKVB(IJ) - AMSKVB(IJ+LS)
     &                   + AMSKVB(IJ) * AMSKVB(IJ+LS))
         CNE(IJ) = AMSKVB(IJ) * AMSKVB(IJ+LSW) *
     &             (3.D0 - AMSKVB(IJ+LW) - AMSKVB(IJ+LS))
         CSE(IJ) = AMSKVB(IJ+LW) * AMSKVB(IJ+LS) *
     &             (3.D0 - AMSKVB(IJ) - AMSKVB(IJ+LSW))
 40      CONTINUE
      END IF

      DO IJ = IJVSTR, IJVEND
         KUP = MAX(NBOTV(IJ)-1, 1)
         FUZ(IJ, KEND) = 
     &                 - WADV(IJ, KEND, 1) * 0.5D0 *
     &                   (UY(IJ, KUP) + UY(IJ, KEND))
         FVZ(IJ, KEND) = 
     &                 - WADV(IJ, KEND, 1) * 0.5D0 *
     &                   (VY(IJ, KUP) + VY(IJ, KEND))

          K = NBOTV(IJ+LN ) - AMSKVB(IJ+LN )
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 2) * (UY(IJ, KEND) + UY(IJ+LN , K)))
          K = NBOTV(IJ+LNE) - AMSKVB(IJ+LNE)
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 3) * (UY(IJ, KEND) + UY(IJ+LNE, K)))
          K = NBOTV(IJ+LE ) - AMSKVB(IJ+LE )
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 4) * (UY(IJ, KEND) + UY(IJ+LE , K)))
          K = NBOTV(IJ+LSE) - AMSKVB(IJ+LSE)
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 5) * (UY(IJ, KEND) + UY(IJ+LSE, K)))
          K = NBOTV(IJ+LS ) - AMSKVB(IJ+LS )
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 6) * (UY(IJ, KEND) + UY(IJ+LS , K)))
          K = NBOTV(IJ+LSW) - AMSKVB(IJ+LSW)
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 7) * (UY(IJ, KEND) + UY(IJ+LSW, K)))
          K = NBOTV(IJ+LW ) - AMSKVB(IJ+LW )
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 8) * (UY(IJ, KEND) + UY(IJ+LW , K)))
          K = NBOTV(IJ+LNW ) - AMSKVB(IJ+LNW )
          FUZ(IJ, KEND) = FUZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 9) * (UY(IJ, KEND) + UY(IJ+LNW, K)))

          K = NBOTV(IJ+LN ) - AMSKVB(IJ+LN )
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 2) * (VY(IJ, KEND) + VY(IJ+LN , K)))
          K = NBOTV(IJ+LNE) - AMSKVB(IJ+LNE)
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 3) * (VY(IJ, KEND) + VY(IJ+LNE, K)))
          K = NBOTV(IJ+LE ) - AMSKVB(IJ+LE )
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 4) * (VY(IJ, KEND) + VY(IJ+LE , K)))
          K = NBOTV(IJ+LSE) - AMSKVB(IJ+LSE)
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 5) * (VY(IJ, KEND) + VY(IJ+LSE, K)))
          K = NBOTV(IJ+LS ) - AMSKVB(IJ+LS )
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 6) * (VY(IJ, KEND) + VY(IJ+LS , K)))
          K = NBOTV(IJ+LSW) - AMSKVB(IJ+LSW)
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 7) * (VY(IJ, KEND) + VY(IJ+LSW, K)))
          K = NBOTV(IJ+LW ) - AMSKVB(IJ+LW )
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 8) * (VY(IJ, KEND) + VY(IJ+LW , K)))
          K = NBOTV(IJ+LNW ) - AMSKVB(IJ+LNW )
          FVZ(IJ, KEND) = FVZ(IJ,KEND) - 0.5D0 *
     &      ( WADV(IJ, KEND, 9) * (VY(IJ, KEND) + VY(IJ+LNW, K)))
      END DO

      DO IJ = IJVSTR, IJVEND+NXDIM
         FUY(IJ, KEND) =
     &                 - (  CYE(IJ) * VADV(IJ, KEND) * HXT(IJ)
     &                    + CYW(IJ+LE) * VADV(IJ+LE, KEND) * HXT(IJ+LE))
     &                    / 6.D0 *
     &                   (UY(IJ, KEND) + UY(IJ+LS, KEND)) * 0.5D0
         FVY(IJ, KEND) =
     &                 - (  CYE(IJ) * VADV(IJ, KEND) * HXT(IJ)
     &                    + CYW(IJ+LE) * VADV(IJ+LE, KEND) * HXT(IJ+LE))
     &                   / 6.D0 *
     &                   (VY(IJ, KEND) + VY(IJ+LS, KEND)) * 0.5D0
      END DO

      DO IJ = IJVSTR, IJVEND+1
         FUX(IJ, KEND) =
     &                 - (  CXN(IJ) * UADV(IJ, KEND) * HYT(IJ)
     &                    + CXS(IJ+LN) * UADV(IJ+LN, KEND) * HYT(IJ+LN))
     &                   / 6.D0 *
     &                   (UY(IJ, KEND) + UY(IJ+LW, KEND)) * 0.5D0
         FVX(IJ, KEND) =
     &                 - (  CXN(IJ) * UADV(IJ, KEND) * HYT(IJ)
     &                    + CXS(IJ+LN) * UADV(IJ+LN, KEND) * HYT(IJ+LN))
     &                   / 6.D0 *
     &                   (VY(IJ, KEND) + VY(IJ+LW, KEND)) * 0.5D0
      END DO

      DO IJ = IJVSTR, IJVEND+NXDIM+1
         FUNE(IJ, KEND) =
     &                  - CNE(IJ) / 6.D0 *
     &                    (  UADV(IJ, KEND) * DY(IJ) * HYT(IJ)
     &                     + VADV(IJ, KEND) * DX * HXT(IJ)) *
     &                    (UY(IJ, KEND) + UY(IJ+LSW, KEND)) * 0.5D0
         FVNE(IJ, KEND) =
     &                  - CNE(IJ) / 6.D0 *
     &                    (  UADV(IJ, KEND) * DY(IJ) * HYT(IJ)
     &                     + VADV(IJ, KEND) * DX * HXT(IJ)) *
     &                    (VY(IJ, KEND) + VY(IJ+LSW, KEND)) * 0.5D0
      END DO

      DO IJ = IJVSTR, IJVEND+NXDIM
         FUSE(IJ, KEND) =
     &                  - CSE(IJ) / 6.D0 *
     &                    (  UADV(IJ, KEND) * DY(IJ) * HYT(IJ)
     &                     - VADV(IJ, KEND) * DX * HXT(IJ)) *
     &                    (UY(IJ+LW, KEND) + UY(IJ+LS, KEND)) * 0.5D0
         FVSE(IJ, KEND) =
     &                  - CSE(IJ) / 6.D0 *
     &                    (  UADV(IJ, KEND) * DY(IJ) * HYT(IJ)
     &                     - VADV(IJ, KEND) * DX * HXT(IJ)) *
     &                    (VY(IJ+LW, KEND) + VY(IJ+LS, KEND)) * 0.5D0
      END DO

      DO IJ = IJVSTR, IJVEND
         GX(IJ, KEND) = (  GX(IJ, KEND)
     &                   + (  (  (  FUX(IJ+LE, KEND)
     &                            - FUX(IJ, KEND)) * RX
     &                         + (  FUY(IJ+LN, KEND)
     &                            - FUY(IJ, KEND)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + (  FUNE(IJ+LNE, KEND)
     &                         - FUNE(IJ    , KEND)
     &                         + FUSE(IJ+LE , KEND)
     &                         - FUSE(IJ+LN , KEND)) *
     &                        RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                      + FUZ(IJ, KEND)) * RZ(IJ)
     &                   + VY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)
     &                   - UY(IJ, KEND) * VY(IJ, KEND) * HXYU(IJ)) *
     &                  AMSKVB(IJ)
         GY(IJ, KEND) = (  GY(IJ, KEND)
     &                   + (   (  (  FVX(IJ+LE, KEND)
     &                             - FVX(IJ, KEND)) * RX
     &                          + (  FVY(IJ+LN, KEND)
     &                             - FVY(IJ, KEND)) * RYM(IJ)) *
     &                         RXU(IJ) * RYU(IJ)
     &                      + (  FVNE(IJ+LNE, KEND)
     &                         - FVNE(IJ    , KEND)
     &                         + FVSE(IJ+LE , KEND)
     &                         - FVSE(IJ+LN , KEND)) *
     &                        RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
     &                      + FVZ(IJ, KEND)) * RZ(IJ)
     &                   + UY(IJ, KEND) * UY(IJ, KEND) * HXYU(IJ)
     &                   - UY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)) *
     &                  AMSKVB(IJ)
      END DO

      RETURN
      END
#endif
