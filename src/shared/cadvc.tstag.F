      SUBROUTINE ADVVEL(
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
c     '03.04.22  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '10.04.14  M.Kurogi: staggered time stepping
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8 C1,C2,C3
      PARAMETER (C1=23.D0/12.D0, C2=-16.D0/12.D0, C3=5.D0/12.D0)
      REAL*8      XX1(NXYDIM, NZDIM),     YY1(NXYDIM, NZDIM)
      REAL*8      XX2(NXYDIM, NZDIM),     YY2(NXYDIM, NZDIM)
      REAL*8      XX3(NXYDIM, NZDIM),     YY3(NXYDIM, NZDIM)
      SAVE XX1,XX2,XX3, YY1,YY2,YY3

      INTEGER NCALL
      DATA NCALL /0/
      SAVE NCALL
      LOGICAL OEOF
      SAVE OEOF


      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HY(NXYDIM)       ,     HX(NXYDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8     FUZ(NXYDIM, NZDIM),    FVZ(NXYDIM, NZDIM)
      REAL*8      RZ(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ, RZ, RZM
*** RZ is taken over from the previous subroutine VSCVEL ***

      REAL*8     DIV(NXYDIM, KSTR:KSTR+KZ-1)
      REAL*8   HVBOT(NXYDIM)
      INTEGER     IJ,      K

      IF (OINIT) THEN
         DO K=1,NZDIM
         DO IJ=1,NXYDIM
            XX1(IJ,K)=0.D0
            XX2(IJ,K)=0.D0
            XX3(IJ,K)=0.D0
            YY1(IJ,K)=0.D0
            YY2(IJ,K)=0.D0
            YY3(IJ,K)=0.D0
         END DO
         END DO
#ifdef OPT_TRIPOLE
         CALL RSTADD(XX2, OEOF, NXDIM, NYDIM, NZDIM, 'XX2', 'OCN',
     &                                            -1.D0,  -1,  -1)
         CALL RSTADD(XX3, OEOF, NXDIM, NYDIM, NZDIM, 'XX3', 'OCN',
     &                                            -1.D0,  -1,  -1)
         CALL RSTADD(YY2, OEOF, NXDIM, NYDIM, NZDIM, 'YY2', 'OCN',
     &                                            -1.D0,  -1,  -1)
         CALL RSTADD(YY3, OEOF, NXDIM, NYDIM, NZDIM, 'YY3', 'OCN',
     &                                            -1.D0,  -1,  -1)
#else
         CALL RSTADD(XX2, OEOF, NXDIM, NYDIM, NZDIM, 'XX2', 'OCN')
         CALL RSTADD(XX3, OEOF, NXDIM, NYDIM, NZDIM, 'XX3', 'OCN')
         CALL RSTADD(YY2, OEOF, NXDIM, NYDIM, NZDIM, 'YY2', 'OCN')
         CALL RSTADD(YY3, OEOF, NXDIM, NYDIM, NZDIM, 'YY3', 'OCN')
#endif
         RETURN
      END IF

      IF (OFINAL) THEN
         CALL FINADD(XX2, NXDIM, NYDIM, NZDIM, 'XX2', 'OCN')
         CALL FINADD(XX3, NXDIM, NYDIM, NZDIM, 'XX3', 'OCN')
         CALL FINADD(YY2, NXDIM, NYDIM, NZDIM, 'YY2', 'OCN')
         CALL FINADD(YY3, NXDIM, NYDIM, NZDIM, 'YY3', 'OCN')
         RETURN
      END IF


      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            FUX(IJ, K) = 0.D0
            FVX(IJ, K) = 0.D0
            FUY(IJ, K) = 0.D0
            FVY(IJ, K) = 0.D0
            FUZ(IJ, K) = 0.D0
            FVZ(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO IJ = IJVSTR, IJVEND
         HVBOT(IJ) = (  (HY(IJ)    + HY(IJ+LE) ) * DY(IJ)
     &                + (HY(IJ+LN) + HY(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
         HVBOT(IJ) = 1.D0 / HVBOT(IJ)
      END DO

      DO K = KSTR, KSTR+KZ-2
         DO IJ = IJVSTR, IJVEND
            DIV(IJ, K) = (  (  (  UADV(IJ+LE, K) *
     &                            (HYU(IJ+LE) + HYU(IJ))
     &                          - UADV(IJ, K) *
     &                            (HYU(IJ) + HYU(IJ+LW))) *
     &                         0.5D0 * RX
     &                       + (  VADV(IJ+LN, K) *
     &                            (HXU(IJ+LN) + HXU(IJ))
     &                          - VADV(IJ, K) *
     &                            (HXU(IJ) + HXU(IJ+LS))) *
     &                         0.5D0 * RYM(IJ)) * RXU(IJ) * RYU(IJ)
     &                    + (WADV(IJ, K) - WADV(IJ, K+1)) * RS(K)) *
     &                   AMSKV(IJ, K)
         END DO
      END DO
      K = KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            DIV(IJ, K) = (  (  (  UADV(IJ+LE, K) *
     &                            (HYU(IJ+LE) + HYU(IJ))
     &                          - UADV(IJ, K) *
     &                            (HYU(IJ) + HYU(IJ+LW))) *
     &                         0.5D0 * RX
     &                       + (  VADV(IJ+LN, K) *
     &                            (HXU(IJ+LN) + HXU(IJ))
     &                          - VADV(IJ, K) *
     &                            (HXU(IJ) + HXU(IJ+LS))) *
     &                         0.5D0 * RYM(IJ)) * RXU(IJ) * RYU(IJ)
     &                    + (  WADV(IJ, K)
     &                       - WADV(IJ, K+1) * HVBOT(IJ)) * RS(K)) *
     &                   AMSKV(IJ, K)
         END DO

      DO IJ = IJVSTR, IJVEND
         FUZ(IJ, KSTR) = 0.D0
         FVZ(IJ, KSTR) = 0.D0
      END DO
      DO 1210 K = KSTR+1, KEND
         DO 1200 IJ = IJVSTR, IJVEND
            FUZ(IJ, K) =
     &                 - WADV(IJ, K) * 0.5D0 *
     &                   (UY(IJ, K-1) + UY(IJ, K))
            FVZ(IJ, K) =
     &                 - WADV(IJ, K) * 0.5D0 *
     &                   (VY(IJ, K-1) + VY(IJ, K))
 1200    CONTINUE
 1210 CONTINUE

      DO 1310 K = KSTR, KEND
         DO 1300 IJ = IJVSTR, IJVEND+NXDIM
            FUY(IJ, K) =
     &                 - VADV(IJ, K) * 0.5D0 *
     &                   (  UY(IJ, K) * HXU(IJ)
     &                    + UY(IJ+LS, K) * HXU(IJ+LS))
            FVY(IJ, K) =
     &                 - VADV(IJ, K) * 0.5D0 *
     &                   (  VY(IJ, K) * HXU(IJ)
     &                    + VY(IJ+LS, K) * HXU(IJ+LS))
 1300    CONTINUE
 1310 CONTINUE

      DO 1410 K = KSTR, KEND
         DO 1400 IJ = IJVSTR, IJVEND+1
            FUX(IJ, K) =
     &                 - UADV(IJ, K) * 0.5D0 *
     &                   (  UY(IJ, K) * HYU(IJ)
     &                    + UY(IJ+LW, K) * HYU(IJ+LW))
            FVX(IJ, K) =
     &                 - UADV(IJ, K) * 0.5D0 *
     &                   (  VY(IJ, K) * HYU(IJ)
     &                    + VY(IJ+LW, K) * HYU(IJ+LW))
 1400    CONTINUE
 1410 CONTINUE

      DO 1510 K = KSTR, KSTR+KZ-2
         DO 1500 IJ = IJVSTR, IJVEND
CX            GX(IJ, K) = (  GX(IJ, K)
            XX1(IJ, K) =  (
     &                     (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (FUZ(IJ, K) - FUZ(IJ, K+1)) * RS(K)
     &                   + UY(IJ, K) * DIV(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
CX            GY(IJ, K) = (  GY(IJ, K)
            YY1(IJ, K) = ( 
     &                     (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (FVZ(IJ, K) - FVZ(IJ, K+1)) * RS(K)
     &                   + VY(IJ, K) * DIV(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 1500    CONTINUE
 1510 CONTINUE
      K = KSTR+KZ-1
      DO 1520 IJ = IJVSTR, IJVEND
CX            GX(IJ, K) = (  GX(IJ, K) +
            XX1(IJ, K) = (
     &                     (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FUZ(IJ, K)
     &                      - FUZ(IJ, K+1) * HVBOT(IJ)) * RS(K)
     &                   + UY(IJ, K) * DIV(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
CX            GY(IJ, K) = (  GY(IJ, K)+
            YY1(IJ, K) = (
     &                     (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FVZ(IJ, K)
     &                      - FVZ(IJ, K+1) * HVBOT(IJ)) * RS(K)
     &                   + VY(IJ, K) * DIV(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 1520 CONTINUE
      DO 1540 K = KSTR+KZ, KEND
         DO 1530 IJ = IJVSTR, IJVEND
CX            GX(IJ, K) = (  GX(IJ, K)+
            XX1(IJ, K) = (
     &                     (  (  (  FUX(IJ+LE, K)
     &                            - FUX(IJ   , K)) * RX
     &                         + (  FUY(IJ+LN, K)
     &                            - FUY(IJ   , K)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FUZ(IJ, K) - FUZ(IJ, K+1)) * RZ(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
CX            GY(IJ, K) = (  GY(IJ, K)+
            YY1(IJ, K) = (
     &                     (  (  (  FVX(IJ+LE, K)
     &                            - FVX(IJ   , K)) * RX
     &                         + (  FVY(IJ+LN, K)
     &                            - FVY(IJ   , K)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FVZ(IJ, K) - FVZ(IJ, K+1)) * RZ(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 1530    CONTINUE
 1540 CONTINUE


CCCCCCCCC Adams-Bashforth scheme
      NCALL=NCALL +1
      IF( (NCALL .GE. 3) .OR. (.NOT. OEOF)) THEN
        NCALL=3
        DO K = KSTR, KEND
        DO IJ = IJVSTR, IJVEND
            GX(IJ, K) =  GX(IJ, K) 
     & +  C1 *XX1(IJ, K)  + C2* XX2(IJ, K) + C3*XX3(IJ, K) 

            GY(IJ, K) =  GY(IJ, K) 
     & +  C1 *YY1(IJ, K)  + C2* YY2(IJ, K) + C3*YY3(IJ, K) 
        END DO
        END DO
      ELSE 
         IF(NCALL .EQ. 1) THEN ! forward
            DO K = KSTR, KEND
            DO IJ = IJVSTR, IJVEND
            GX(IJ, K) =  GX(IJ, K)  +  XX1(IJ, K)
            GY(IJ, K) =  GY(IJ, K)  +  YY1(IJ, K)
            END DO
            END DO
         ELSE IF(NCALL .EQ. 2) THEN !2nd order AB
            DO K = KSTR, KEND
            DO IJ = IJVSTR, IJVEND
               GX(IJ, K) =  GX(IJ, K) 
     &       +  1.5D0 *XX1(IJ, K)  - 0.5D0* XX2(IJ, K)

               GY(IJ, K) =  GY(IJ, K) 
     &       +  1.5D0 *YY1(IJ, K)  - 0.5D0* YY2(IJ, K)
            END DO
            END DO
         END IF
      END IF


      DO K = KSTR, KEND
      DO IJ = IJVSTR, IJVEND
          XX3(IJ, K)=XX2(IJ, K) !N-1 > N-2
          YY3(IJ, K)=YY2(IJ, K)

          XX2(IJ, K)=XX1(IJ, K) !N> N-1
          YY2(IJ, K)=YY1(IJ, K)
      END DO
      END DO

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE ADVVLB(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I               UADV,   VADV,   WADV)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Advection term of the BBL momentum eqs.
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8 C1,C2,C3
      PARAMETER (C1=23.D0/12.D0, C2=-16.D0/12.D0, C3=5.D0/12.D0)
      REAL*8      XX1(NXYDIM),     YY1(NXYDIM)
      REAL*8      XX2(NXYDIM),     YY2(NXYDIM)
      REAL*8      XX3(NXYDIM),     YY3(NXYDIM)
      SAVE XX1,XX2,XX3, YY1,YY2,YY3
      INTEGER NCALL
      DATA NCALL /0/
      SAVE NCALL
      LOGICAL OEOF
      SAVE OEOF

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8      RZ(NXYDIM)
      SAVE RZ

      REAL*8     FUX(NXYDIM),    FVX(NXYDIM)
      REAL*8     FUY(NXYDIM),    FVY(NXYDIM)
      REAL*8     FUZ(NXYDIM),    FVZ(NXYDIM)
      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ

      INTEGER     IJ,    KUP

      LOGICAL OFIRST
      SAVE OFIRST


      IF (OINIT) THEN
         DO IJ=1,NXYDIM
            XX1(IJ)=0.D0
            XX2(IJ)=0.D0
            XX3(IJ)=0.D0
            YY1(IJ)=0.D0
            YY2(IJ)=0.D0
            YY3(IJ)=0.D0
         END DO
#ifdef OPT_TRIPOLE
         CALL RSTADD(XX2, OEOF, NXDIM, NYDIM, 1, 'XX2', 'SFC',
     &                                        -1.D0,  -1,  -1)
         CALL RSTADD(XX3, OEOF, NXDIM, NYDIM, 1, 'XX3', 'SFC',
     &                                        -1.D0,  -1,  -1)
         CALL RSTADD(YY2, OEOF, NXDIM, NYDIM, 1, 'YY2', 'SFC',
     &                                        -1.D0,  -1,  -1)
         CALL RSTADD(YY3, OEOF, NXDIM, NYDIM, 1, 'YY3', 'SFC',
     &                                        -1.D0,  -1,  -1)
#else
         CALL RSTADD(XX2, OEOF, NXDIM, NYDIM, 1, 'XX2', 'SFC')
         CALL RSTADD(XX3, OEOF, NXDIM, NYDIM, 1, 'XX3', 'SFC')
         CALL RSTADD(YY2, OEOF, NXDIM, NYDIM, 1, 'YY2', 'SFC')
         CALL RSTADD(YY3, OEOF, NXDIM, NYDIM, 1, 'YY3', 'SFC')
#endif
         RETURN
      END IF

      IF (OFINAL) THEN
         CALL FINADD(XX2, NXDIM, NYDIM, 1, 'XX2', 'SFC')
         CALL FINADD(XX3, NXDIM, NYDIM, 1, 'XX3', 'SFC')
         CALL FINADD(YY2, NXDIM, NYDIM, 1, 'YY2', 'SFC')
         CALL FINADD(YY3, NXDIM, NYDIM, 1, 'YY3', 'SFC')
         RETURN
      END IF


      IF (OFIRST) THEN
         OFIRST = .FALSE.
         DO IJ = 1, NXYDIM
            RZ   (IJ) = 1.D0 / DZV(IJ, KEND)
         END DO
      END IF

      DO IJ = 1, NXYDIM
         FUX(IJ) = 0.D0
         FVX(IJ) = 0.D0
         FUY(IJ) = 0.D0
         FVY(IJ) = 0.D0
      END DO

      DO IJ = IJVSTR, IJVEND
         KUP = MAX(NBOTV(IJ)-1, 1)
         FUZ(IJ) = - WADV(IJ, KEND) * 0.5D0 *
     &               (UY(IJ, KUP) + UY(IJ, KEND))
         FVZ(IJ) = - WADV(IJ, KEND) * 0.5D0 *
     &               (VY(IJ, KUP) + VY(IJ, KEND))
      END DO

      DO IJ = IJVSTR, IJVEND+NXDIM
         FUY(IJ) = - VADV(IJ, KEND) * 0.5D0 *
     &               (  UY(IJ, KEND) * HXU(IJ)
     &                + UY(IJ+LS, KEND) * HXU(IJ+LS))
         FVY(IJ) = - VADV(IJ, KEND) * 0.5D0 *
     &               (  VY(IJ, KEND) * HXU(IJ)
     &                + VY(IJ+LS, KEND) * HXU(IJ+LS))
      END DO

      DO IJ = IJVSTR, IJVEND+1
         FUX(IJ) = - UADV(IJ, KEND) * 0.5D0 *
     &               (  UY(IJ, KEND) * HYU(IJ)
     &                + UY(IJ+LW, KEND) * HYU(IJ+LW))
         FVX(IJ) = - UADV(IJ, KEND) * 0.5D0 *
     &               (  VY(IJ, KEND) * HYU(IJ)
     &                + VY(IJ+LW, KEND) * HYU(IJ+LW))
      END DO

      DO IJ = IJVSTR, IJVEND
CX         GX(IJ, KEND) = (  GX(IJ, KEND)+
         XX1(IJ) = (
     &                     (  (  (FUX(IJ+LE) - FUX(IJ)) * RX
     &                         + (FUY(IJ+LN) - FUY(IJ)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FUZ(IJ)) * RZ(IJ)
     &                   + VY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)
     &                   - UY(IJ, KEND) * VY(IJ, KEND) * HXYU(IJ)) *
     &                  AMSKVB(IJ)
CX         GY(IJ, KEND) = (  GY(IJ, KEND)+
         YY1(IJ) = (
     &                     (  (  (FVX(IJ+LE) - FVX(IJ)) * RX
     &                         + (FVY(IJ+LN) - FVY(IJ)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FVZ(IJ)) * RZ(IJ)
     &                   + UY(IJ, KEND) * UY(IJ, KEND) * HXYU(IJ)
     &                   - UY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)) *
     &                  AMSKVB(IJ)
      END DO


CCCCCCCCC Adams-Bashforth scheme
      NCALL=NCALL +1

      IF( (NCALL .GE. 3) .OR. (.NOT. OEOF)) THEN
      NCALL=3
        DO IJ = IJVSTR, IJVEND
            GX(IJ, KEND) =  GX(IJ, KEND) 
     & +  C1 *XX1(IJ)  + C2* XX2(IJ) + C3*XX3(IJ) 

            GY(IJ, KEND) =  GY(IJ, KEND) 
     & +  C1 *YY1(IJ)  + C2* YY2(IJ) + C3*YY3(IJ) 
        END DO
      ELSE
        IF(NCALL .EQ. 1) THEN ! forward
           DO IJ = IJVSTR, IJVEND
            GX(IJ, KEND) =  GX(IJ, KEND) + XX1(IJ)
            GY(IJ, KEND) =  GY(IJ, KEND) + YY1(IJ)
           END DO
        ELSE IF(NCALL .EQ. 2) THEN !2nd order AB       
          DO IJ = IJVSTR, IJVEND
            GX(IJ, KEND) =  GX(IJ, KEND) 
     &    +  1.5D0 *XX1(IJ)  -0.5D0* XX2(IJ)

            GY(IJ, KEND) =  GY(IJ, KEND) 
     &    +  1.5D0 *YY1(IJ)  -0.5D0* YY2(IJ)
          END DO
        END IF
      END IF




      DO IJ = IJVSTR, IJVEND
          XX3(IJ)=XX2(IJ) !N-1 > N-2
          YY3(IJ)=YY2(IJ)

          XX2(IJ)=XX1(IJ) !N> N-1
          YY2(IJ)=YY1(IJ)
      END DO


      RETURN
      END
#endif
