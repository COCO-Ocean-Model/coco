      SUBROUTINE TNDUVD(
     O                GXX,    GYY,
     M                 GX,     GY,
     I                 XX,     YY,
     I                 UY,     VY,     HY,      R,
     I                 UX,     VX,   UADV,   VADV)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the pressure gradient and Coriolis terms of the equation
c of motion.
c
c  HISTORY
c     '03.04.22  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: for McDougall et al. (2003) eq. of state
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.01.15  T.Suzuki: bug fix vertical integration thanks to Hasumi
c     '10.04.14  M.Kurogi: staggered time stepping
c
c ---------------------------------------------------------------------

      REAL*8       R(NXYDIM, NZDIM),     HY(NXYDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),    GXX(NXYDIM)
      REAL*8      GY(NXYDIM, NZDIM),    GYY(NXYDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8     COF
      SAVE COF

      REAL*8      PX(NXYDIM, NZDIM),    PXM(NXYDIM, NZDIM)
      REAL*8      PY(NXYDIM, NZDIM),    PYM(NXYDIM, NZDIM)
      REAL*8      UU(NXYDIM)       ,     UV(NXYDIM)
      REAL*8      VU(NXYDIM)       ,     VV(NXYDIM)
      REAL*8   HVBOT(NXYDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM)
c      COMMON /WORK/ PX, PXM,  PY, PYM, UU, UV, VU, VV,
c     &              HVBOT, DZSIG
      COMMON /WORK2/ PX, PXM,  PY, PYM

      INTEGER     IJ,      K

      LOGICAL OFIRST
      SAVE OFIRST
      DATA OFIRST / .TRUE. /

      REAL*8 CX1,CX2,CX3
      PARAMETER (CX1=23.D0/12.D0, CX2=-16.D0/12.D0, CX3=5.D0/12.D0)

      REAL*8      XX1(NXYDIM),     YY1(NXYDIM)
      REAL*8      XX2(NXYDIM),     YY2(NXYDIM)
      REAL*8      XX3(NXYDIM),     YY3(NXYDIM)
      SAVE XX1,XX2,XX3, YY1,YY2,YY3
      INTEGER NCALL
      DATA NCALL /0/
      SAVE NCALL

      LOGICAL OEOF
      SAVE OEOF

      IF (OINIT) THEN
         DO IJ=1,NXYDIM
            XX2(IJ)=0.D0
            XX3(IJ)=0.D0
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
         COF = -0.25D0 * gravit / rhoo * 1.D-3
      END IF

      DO IJ = IJVSTR, IJVEND+NXDIM
         HVBOT(IJ) = (  (HY(IJ) + HY(IJ+LE)) * DY(IJ)
     &                + (HY(IJ+LN) + HY(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
      END DO
      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND+NXDIM
            DZSIG(IJ, K) = DS(K) * HVBOT(IJ)
         END DO
      END DO

      DO 100 IJ = IJVSTR, IJVEND
         PXM(IJ, KSTR) = (  R(IJ+LNE, KSTR) + R(IJ+LE, KSTR)
     &                    - R(IJ+LN , KSTR) - R(IJ   , KSTR)) *
     &                   COF * RX * RXU(IJ) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
         PX (IJ, KSTR) = PXM(IJ, KSTR)
         GX (IJ, KSTR) = GX(IJ, KSTR)
     &                 + (COR(IJ) * VX(IJ, KSTR) + PX(IJ, KSTR)) *
     &                   AMSKV(IJ, KSTR)
         XX (IJ, KSTR) = XX(IJ, KSTR)
     &                 + PX(IJ, KSTR) * AMSKV(IJ, KSTR)
         GXX(IJ)       = XX(IJ, KSTR) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
         PYM(IJ, KSTR) = (  R(IJ+LNE, KSTR) + R(IJ+LN, KSTR)
     &                    - R(IJ+LE , KSTR) - R(IJ   , KSTR)) *
     &                   COF * RYM(IJ) * RYU(IJ) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
         PY (IJ, KSTR) = PYM(IJ, KSTR)
         GY (IJ, KSTR) = GY(IJ, KSTR)
     &                 - (COR(IJ) * UX(IJ, KSTR) - PY(IJ, KSTR)) *
     &                   AMSKV(IJ, KSTR)
         YY (IJ, KSTR) = YY(IJ, KSTR)
     &                 + PY(IJ, KSTR) * AMSKV(IJ, KSTR)
         GYY(IJ)       = YY(IJ, KSTR) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
  100 CONTINUE

      DO 210 K = KSTR+1, KSTR+KZ-1
         DO 200 IJ = IJVSTR, IJVEND
            PXM(IJ, K) = (  R(IJ+LNE, K) + R(IJ+LE, K)
     &                    - R(IJ+LN , K) - R(IJ   , K)) *
     &                   COF * RX * RXU(IJ) * DZSIG(IJ, K) *
     &                   AMSKV(IJ, K)
            PX (IJ, K) = (PX(IJ, K-1) + PXM(IJ, K-1) + PXM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GX (IJ, K) = GX(IJ, K)
     &                 + (COR(IJ) * VX(IJ, K) + PX(IJ, K)) *
     &                   AMSKV(IJ, K)
            XX (IJ, K) = XX(IJ, K)
     &                 + PX(IJ, K) * AMSKV(IJ, K)
            GXX(IJ)    = GXX(IJ)
     &                 + XX(IJ, K) * DZSIG(IJ, K) * AMSKV(IJ, K)
            PYM(IJ, K) = (  R(IJ+LNE, K) + R(IJ+LN, K)
     &                    - R(IJ+LE , K) - R(IJ   , K)) *
     &                   COF * RYM(IJ) * RYU(IJ) * DZSIG(IJ, K) *
     &                   AMSKV(IJ, K)
            PY (IJ, K) = (PY(IJ, K-1) + PYM(IJ, K-1) + PYM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GY (IJ, K) = GY(IJ, K)
     &                 - (COR(IJ) * UX(IJ, K) - PY(IJ, K)) *
     &                   AMSKV(IJ, K)
            YY (IJ, K) = YY(IJ, K)
     &                 + PY(IJ, K) * AMSKV(IJ, K)
            GYY(IJ)    = GYY(IJ) + YY(IJ, K) * DZSIG(IJ, K) *
     &                   AMSKV(IJ, K)
  200    CONTINUE
  210 CONTINUE

      DO 240 K = KSTR+KZ, KEND
         DO 230 IJ = IJVSTR, IJVEND
            PXM(IJ, K) = (  R(IJ+LNE, K) - R(IJ+LN, K)
     &                    + R(IJ+LE , K) - R(IJ   , K)) *
     &                   COF * RX * RXU(IJ) * DZV(IJ, K) *
     &                   AMSKV(IJ, K)
            PX (IJ, K) = (PX(IJ, K-1) + PXM(IJ, K-1) + PXM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GX (IJ, K) = GX(IJ, K)
     &                 + (COR(IJ) * VX(IJ, K) + PX(IJ, K)) *
     &                   AMSKV(IJ, K)
            XX (IJ, K) = XX(IJ, K)
     &                 + PX(IJ, K) * AMSKV(IJ, K)
            GXX(IJ)    = GXX(IJ) + XX(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            PYM(IJ, K) = (  R(IJ+LNE, K) + R(IJ+LN, K)
     &                    - R(IJ+LE , K) - R(IJ   , K)) *
     &                   COF * RYM(IJ) * RYU(IJ) * DZV(IJ, K) *
     &                   AMSKV(IJ, K)
            PY (IJ, K) = (PY(IJ, K-1) + PYM(IJ, K-1) + PYM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GY (IJ, K) = GY(IJ, K)
     &                 - (COR(IJ) * UX(IJ, K) - PY(IJ, K)) *
     &                   AMSKV(IJ, K)
            YY (IJ, K) = YY(IJ, K)
     &                 + PY(IJ, K) * AMSKV(IJ, K)
            GYY(IJ)    = GYY(IJ) + YY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
  230    CONTINUE
  240 CONTINUE

C -----------------------------------------------------------
      DO IJ=1,NXYDIM
            XX1(IJ)=0.D0
            YY1(IJ)=0.D0
      END DO

      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND+NXDIM
            UU(IJ) = UADV(IJ, K) * 0.5D0 *
     &               (  UY(IJ, K) * HYU(IJ)
     &                + UY(IJ+LW, K) * HYU(IJ+LW))
            VU(IJ) = VADV(IJ, K) * 0.5D0 *
     &               (  UY(IJ, K) * HXU(IJ)
     &                + UY(IJ+LS, K) * HXU(IJ+LS))
            UV(IJ) = UADV(IJ, K) * 0.5D0 *
     &               (  VY(IJ, K) * HYU(IJ)
     &                + VY(IJ+LW, K) * HYU(IJ+LW))
            VV(IJ) = VADV(IJ, K) * 0.5D0 *
     &               (  VY(IJ, K) * HXU(IJ)
     &                + VY(IJ+LS, K) * HXU(IJ+LS))
         END DO
         DO IJ = IJVSTR, IJVEND
cx            GXX(IJ) = GXX(IJ)
            XX1(IJ) = XX1(IJ)
     &              - (  (UU(IJ+LE) - UU(IJ)) * RX
     &                 + (VU(IJ+LN) - VU(IJ)) * RYM(IJ)) *
     &                RXU(IJ) * RYU(IJ) * DZSIG(IJ,K)
cx            GYY(IJ) = GYY(IJ)
            YY1(IJ) = YY1(IJ)
     &              - (  (UV(IJ+LE) - UV(IJ)) * RX
     &                 + (VV(IJ+LN) - VV(IJ)) * RYM(IJ)) *
     &                RXU(IJ) * RYU(IJ) * DZSIG(IJ,K)
         END DO
      END DO
C -----------------------------------------------------------
C      DO IJ = 1, NXYDIM
C         UU(IJ) = 0.D0
C         UV(IJ) = 0.D0
C         VU(IJ) = 0.D0
C         VV(IJ) = 0.D0
C      END DO
C
C      DO 320 K = KSTR, KSTR+KZ-1
C         DO 310 IJ = IJVSTR, IJVEND+NXDIM
C            UU(IJ) = UU(IJ)
C     &             + UADV(IJ, K) * 0.5D0 *
C     &               (  UY(IJ, K) * HYU(IJ)
C     &                + UY(IJ+LW, K) * HYU(IJ+LW)) * DS(K)
C            VU(IJ) = VU(IJ)
C     &             + VADV(IJ, K) * 0.5D0 *
C     &               (  UY(IJ, K) * HXU(IJ)
C     &                + UY(IJ+LS, K) * HXU(IJ+LS)) * DS(K)
C            UV(IJ) = UV(IJ)
C     &             + UADV(IJ, K) * 0.5D0 *
C     &               (  VY(IJ, K) * HYU(IJ)
C     &                + VY(IJ+LW, K) * HYU(IJ+LW)) * DS(K)
C            VV(IJ) = VV(IJ)
C     &             + VADV(IJ, K) * 0.5D0 *
C     &               (  VY(IJ, K) * HXU(IJ)
C     &                + VY(IJ+LS, K) * HXU(IJ+LS)) * DS(K)
C 310     CONTINUE
C 320  CONTINUE
C
C      DO IJ = IJVSTR, IJVEND
C         GXX(IJ) = GXX(IJ)
C     &           - (  (UU(IJ+LE) - UU(IJ)) * RX
C     &              + (VU(IJ+LN) - VU(IJ)) * RYM(IJ)) *
C     &             RXU(IJ) * RYU(IJ) * HVBOT(IJ)
C         GYY(IJ) = GYY(IJ)
C     &           - (  (UV(IJ+LE) - UV(IJ)) * RX
C     &              + (VV(IJ+LN) - VV(IJ)) * RYM(IJ)) *
C     &             RXU(IJ) * RYU(IJ) * HVBOT(IJ)
C      ENDDO
C ------------------------------------------------------------
      DO IJ = 1, NXYDIM
         UU(IJ) = 0.D0
         UV(IJ) = 0.D0
         VU(IJ) = 0.D0
         VV(IJ) = 0.D0
      END DO
C
      DO 410 K = KSTR+KZ, KEND
         DO 400 IJ = IJVSTR, IJVEND+NXDIM
            UU(IJ) = UU(IJ)
     &             + UADV(IJ, K) * 0.5D0 *
     &               (  UY(IJ, K) * HYU(IJ)
     &                + UY(IJ+LW, K) * HYU(IJ+LW))
            VU(IJ) = VU(IJ)
     &             + VADV(IJ, K) * 0.5D0 *
     &               (  UY(IJ, K) * HXU(IJ)
     &                + UY(IJ+LS, K) * HXU(IJ+LS))
            UV(IJ) = UV(IJ)
     &             + UADV(IJ, K) * 0.5D0 *
     &               (  VY(IJ, K) * HYU(IJ)
     &                + VY(IJ+LW, K) * HYU(IJ+LW))
            VV(IJ) = VV(IJ)
     &             + VADV(IJ, K) * 0.5D0 *
     &               (  VY(IJ, K) * HXU(IJ)
     &                + VY(IJ+LS, K) * HXU(IJ+LS))
  400    CONTINUE
  410 CONTINUE

      DO 500 IJ = IJVSTR, IJVEND
cx         GXX(IJ) = GXX(IJ)
         XX1(IJ) = XX1(IJ)
     &           - (  (UU(IJ+LE) - UU(IJ)) * RX
     &              + (VU(IJ+LN) - VU(IJ)) * RYM(IJ)) *
     &             RXU(IJ) * RYU(IJ)
cx         GYY(IJ) = GYY(IJ)
         YY1(IJ) = YY1(IJ)
     &           - (  (UV(IJ+LE) - UV(IJ)) * RX
     &              + (VV(IJ+LN) - VV(IJ)) * RYM(IJ)) *
     &             RXU(IJ) * RYU(IJ)
  500 CONTINUE

      DO 600 IJ = IJVSTR, IJVEND
         UU(IJ) = UY(IJ, KSTR) * UY(IJ, KSTR) *
     &            DZSIG(IJ, KSTR) * AMSKV(IJ, KSTR)
         UV(IJ) = UY(IJ, KSTR) * VY(IJ, KSTR) *
     &            DZSIG(IJ, KSTR) * AMSKV(IJ, KSTR)
         VV(IJ) = VY(IJ, KSTR) * VY(IJ, KSTR) *
     &            DZSIG(IJ, KSTR) * AMSKV(IJ, KSTR)
  600 CONTINUE

      DO 620 K = KSTR+1, KSTR+KZ-1
         DO 610 IJ = IJVSTR, IJVEND
            UU(IJ) = UU(IJ)
     &             + UY(IJ, K) * UY(IJ, K) *
     &               DZSIG(IJ, K) * AMSKV(IJ, K)
            UV(IJ) = UV(IJ)
     &             + UY(IJ, K) * VY(IJ, K) *
     &               DZSIG(IJ, K) * AMSKV(IJ, K)
            VV(IJ) = VV(IJ)
     &             + VY(IJ, K) * VY(IJ, K) *
     &               DZSIG(IJ, K) * AMSKV(IJ, K)
  610    CONTINUE
  620 CONTINUE

      DO 710 K = KSTR+KZ, KEND
         DO 700 IJ = IJVSTR, IJVEND
            UU(IJ) = UU(IJ)
     &             + UY(IJ, K) * UY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            UV(IJ) = UV(IJ)
     &             + UY(IJ, K) * VY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            VV(IJ) = VV(IJ)
     &             + VY(IJ, K) * VY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
  700    CONTINUE
  710 CONTINUE

      DO 800 IJ = IJVSTR, IJVEND
cx         GXX(IJ) = GXX(IJ) + VV(IJ) * HYXU(IJ) - UV(IJ) * HXYU(IJ)
cx         GYY(IJ) = GYY(IJ) + UU(IJ) * HXYU(IJ) - UV(IJ) * HYXU(IJ)
         XX1(IJ) = XX1(IJ) + VV(IJ) * HYXU(IJ) - UV(IJ) * HXYU(IJ)
         YY1(IJ) = YY1(IJ) + UU(IJ) * HXYU(IJ) - UV(IJ) * HYXU(IJ)
         XX1(IJ) = XX1(IJ) * AMSKV(IJ,KSTR)
         YY1(IJ) = YY1(IJ) * AMSKV(IJ,KSTR)
  800 CONTINUE


CCCCCCCCC Adams-Bashforth scheme
      NCALL=NCALL +1

      IF( (NCALL .GE. 3) .OR. (.NOT. OEOF)) THEN
      NCALL=3
        DO IJ = IJVSTR, IJVEND
            GXX(IJ) =  GXX(IJ) 
     & +  CX1 *XX1(IJ)  + CX2* XX2(IJ) + CX3*XX3(IJ) 

            GYY(IJ) =  GYY(IJ) 
     & +  CX1 *YY1(IJ)  + CX2* YY2(IJ) + CX3*YY3(IJ) 
        END DO
      ELSE
        IF(NCALL .EQ. 1) THEN ! forward
           DO IJ = IJVSTR, IJVEND
            GXX(IJ) =  GXX(IJ) + XX1(IJ)
            GYY(IJ) =  GYY(IJ) + YY1(IJ)
           END DO
        ELSE IF(NCALL .EQ. 2) THEN !2nd order AB       
          DO IJ = IJVSTR, IJVEND
            GXX(IJ) =  GXX(IJ) 
     &    +  1.5D0 *XX1(IJ)  -0.5D0* XX2(IJ)

            GYY(IJ) =  GYY(IJ) 
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

#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE TNDUVB(
     M               GXX,    GYY,
     M                GX,     GY,
     I                XX,     YY,
     I                UY,     VY,     TY,  
     I                UX,     VX)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Pressure gradient and Coriolis terms for the BBL momentum eqs. See
c Nakano and Suginohara (2002, JPO).
c
c ---------------------------------------------------------------------

      REAL*8     GXX(NXYDIM)       ,    GYY(NXYDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)     

      REAL*8      PX(NXYDIM, NZDIM),    PXM(NXYDIM, NZDIM)
      REAL*8      PY(NXYDIM, NZDIM),    PYM(NXYDIM, NZDIM)
      COMMON /WORK2/ PX, PXM,  PY, PYM

      REAL*8  C0B(NXYDIM), C1B(NXYDIM), C2B(NXYDIM), C3B(NXYDIM)
      REAL*8  C4B(NXYDIM), C5B(NXYDIM), C6B(NXYDIM)
      REAL*8  D0B(NXYDIM), D1B(NXYDIM), D2B(NXYDIM), D3B(NXYDIM)
      REAL*8  D4B(NXYDIM), D5B(NXYDIM), D6B(NXYDIM), D7B(NXYDIM)
      REAL*8  D8B(NXYDIM), D9B(NXYDIM)
      SAVE C0B, C1B, C2B, C3B, C4B, C5B, C6B,
     &     D0B, D1B, D2B, D3B, D4B, D5B, D6B, D7B, D8B, D9B

      REAL*8  DEPT0(NXYDIM)
      REAL*8    COF
      LOGICAL OFIRST
      SAVE DEPT0, COF, OFIRST

      REAL*8    RBBL,     TB,     SB
      REAL*8  C0, C1, C2, C3, C4, C5, C6
      REAL*8  D0, D1, D2, D3, D4, D5, D6, D7, D8, D9
      REAL*8      TL,     SL,   RTMP
      REAL*8       P,     PE,     PN,    PNE
      INTEGER     IJ,      K,     KK
      INTEGER   IJLE,   IJLN,   IJNE

      DATA OFIRST / .TRUE. /

      REAL*8 CX1,CX2,CX3
      PARAMETER (CX1=23.D0/12.D0, CX2=-16.D0/12.D0, CX3=5.D0/12.D0)

      REAL*8      XX1(NXYDIM),     YY1(NXYDIM)
      REAL*8      XX2(NXYDIM),     YY2(NXYDIM)
      REAL*8      XX3(NXYDIM),     YY3(NXYDIM)
      SAVE XX1,XX2,XX3, YY1,YY2,YY3
      INTEGER NCALL
      DATA NCALL /0/
      SAVE NCALL

      LOGICAL OEOF
      SAVE OEOF

      RBBL(TB, SB, C0, C1, C2, C3, C4, C5, C6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9) = 
     &    (C0 + (C1 + (C2 + C3 * TB) * TB) * TB
     &        + (C4 + C5 * TB + C6 * SB) * SB)
     &  / (D0 + (D1 + (D2 + (D3 + D4 * TB) * TB) * TB) * TB
     &        + (D5 + (D6 + D7 * TB * TB) * TB
     &              + (D8 + D9 * TB * TB) * SQRT(SB)) * SB)
     &  - 1.D3

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
         CALL SECFBB(
     O               C0B,    C1B,    C2B,    C3B,
     O               C4B,    C5B,    C6B,
     O               D0B,    D1B,    D2B,    D3B,    D4B,
     O               D5B,    D6B,    D7B,    D8B,    D9B)

         COF = -0.5D0 * gravit / rhoo * 1.D-3
         DO IJ = 1, NXYDIM
            DEPT0(IJ) = 0.D0
            DO K = KSTR, NBOT(IJ)-1
               DEPT0(IJ) = DEPT0(IJ) + DZ0(K)
            END DO
         END DO
      END IF

      DO IJ = IJVSTR, IJVEND

         K = MAX(NBOTV(IJ)-1, 1)
         IJLE = IJ + LE
         IJLN = IJ + LN
         IJNE = IJ + LNE

         TL = TY(IJ, K, 1)
         SL = TY(IJ, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         P = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJ)-1
            TL = TY(IJ, KK, 1)
            SL = TY(IJ, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            P = P + RTMP * DZ0(KK)
         END DO
         TL = TY(IJ, KEND, 1)
         SL = TY(IJ, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         P = P + RTMP * DZ0(KEND) * 0.5D0

         TL = TY(IJLE, K, 1)
         SL = TY(IJLE, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PE = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJLE)-1
            TL = TY(IJLE, KK, 1)
            SL = TY(IJLE, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            PE = PE + RTMP * DZ0(KK)
         END DO
         TL = TY(IJLE, KEND, 1)
         SL = TY(IJLE, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PE = PE + RTMP * DZ0(KEND) * 0.5D0

         TL = TY(IJLN, K, 1)
         SL = TY(IJLN, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PN = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJLN)-1
            TL = TY(IJLN, KK, 1)
            SL = TY(IJLN, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            PN = PN + RTMP * DZ0(KK)
         END DO
         TL = TY(IJLN, KEND, 1)
         SL = TY(IJLN, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PN = PN + RTMP * DZ0(KEND) * 0.5D0

         TL = TY(IJNE, K, 1)
         SL = TY(IJNE, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PNE = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJNE)-1
            TL = TY(IJNE, KK, 1)
            SL = TY(IJNE, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            PNE = PNE + RTMP * DZ0(KK)
         END DO
         TL = TY(IJNE, KEND, 1)
         SL = TY(IJNE, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PNE = PNE + RTMP * DZ0(KEND) * 0.5D0

         TL = 0.25D0 * (  TY(IJ  , KEND, 1) + TY(IJLE, KEND, 1)
     &                  + TY(IJLN, KEND, 1) + TY(IJNE, KEND, 1))
         SL = 0.25D0 * (  TY(IJ  , KEND, 2) + TY(IJLE, KEND, 2)
     &                  + TY(IJLN, KEND, 2) + TY(IJNE, KEND, 2))
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))

         PX(IJ, KEND) = PX(IJ, K)
     &                + (  PNE + PE - PN - P
     &                   - RTMP * (  DEPT0(IJNE) + DEPT0(IJLE)
     &                             - DEPT0(IJLN) - DEPT0(IJ  ))) *
     &                  COF * RX * RXU(IJ) * AMSKVB(IJ)
         GX(IJ, KEND) = GX(IJ, KEND)
     &                + (COR(IJ) * VX(IJ, KEND) + PX(IJ, KEND)) *
     &                  AMSKVB(IJ)
         XX(IJ, KEND) = XX(IJ, KEND) + PX(IJ, KEND) * AMSKVB(IJ)
         GXX(IJ)      = GXX(IJ)
     &                + XX(IJ, KEND) * DZV(IJ, KEND) * AMSKVB(IJ)

         PY(IJ, KEND) = PY(IJ, K)
     &                + (  PNE + PN - PE - P
     &                   - RTMP * (  DEPT0(IJNE) + DEPT0(IJLN)
     &                             - DEPT0(IJLE) - DEPT0(IJ  ))) *
     &                  COF * RYM(IJ) * RYU(IJ) * AMSKVB(IJ)
         GY(IJ, KEND) = GY(IJ, KEND)
     &                - (COR(IJ) * UX(IJ, KEND) - PY(IJ, KEND)) *
     &                  AMSKVB(IJ)
         YY(IJ, KEND) = YY(IJ, KEND) + PY(IJ, KEND) * AMSKVB(IJ)
         GYY(IJ)      = GYY(IJ)
     &                + YY(IJ, KEND) * DZV(IJ, KEND) * AMSKVB(IJ)

      END DO

      DO IJ = IJVSTR, IJVEND
cx         GXX(IJ) = GXX(IJ) +
         XX1(IJ) =
     &             (  VY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)
     &              - UY(IJ, KEND) * VY(IJ, KEND) * HXYU(IJ)) *
     &             DZV(IJ, KEND) * AMSKVB(IJ)
cx         GYY(IJ) = GYY(IJ) +
         YY1(IJ) = 
     &             (  UY(IJ, KEND) * UY(IJ, KEND) * HXYU(IJ)
     &              - UY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)) *
     &             DZV(IJ, KEND) * AMSKVB(IJ)
      END DO

CCCCCCCCC Adams-Bashforth scheme
      NCALL=NCALL +1

      IF( (NCALL .GE. 3) .OR. (.NOT. OEOF)) THEN
      NCALL=3
        DO IJ = IJVSTR, IJVEND
            GXX(IJ) =  GXX(IJ) 
     & +  CX1 *XX1(IJ)  + CX2* XX2(IJ) + CX3*XX3(IJ) 

            GYY(IJ) =  GYY(IJ) 
     & +  CX1 *YY1(IJ)  + CX2* YY2(IJ) + CX3*YY3(IJ) 
        END DO
      ELSE
        IF(NCALL .EQ. 1) THEN ! forward
           DO IJ = IJVSTR, IJVEND
            GXX(IJ) =  GXX(IJ) + XX1(IJ)
            GYY(IJ) =  GYY(IJ) + YY1(IJ)
           END DO
        ELSE IF(NCALL .EQ. 2) THEN !2nd order AB       
          DO IJ = IJVSTR, IJVEND
            GXX(IJ) =  GXX(IJ) 
     &    +  1.5D0 *XX1(IJ)  -0.5D0* XX2(IJ)

            GYY(IJ) =  GYY(IJ) 
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
