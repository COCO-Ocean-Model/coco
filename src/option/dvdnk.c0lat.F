      SUBROUTINE VDIFF(
     O               AMV,    AHV,
     I                UY,     VY,      R,   TAUX,   TAUY,
     I                TY,     HY)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Sea surface mixed layer model of Noh and Kim (1999, JGR).
c
c  HISTORY
c     '02.08.19  H.Hasumi: from COCO3.4
c     '03.10.14  H.Hasumi: bug fix (typographical error in the
c                          original paper)
c     '03.10.22  A.Oka: bug fix ( calculation of AMV, DML )
c     '03.11.14  T.Suzuki: reduce background mixing on equator
c     '03.11.27  A.Oka: adjusted semi-implicit time integration
c                       time filter to TKE in strong damping phase (ALSC)
c                       time filter to AHV (AFLT)
c     '04.01.21  A.Oka: ALPHC, ALPH, BETA, BETAC. Pr depends on Ri.
c     '04.01.27  A.Oka: two way for calculation of DML.
c     '04.02.23  H.Hasumi: imported from MIROC3.2
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: McDougall et al. (2003) eq. of state
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.01.19  T.Suzuki: bug fix
c                          (TKE not initialized when read from RSTO)
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       R(NXYDIM, NZDIM)
      REAL*8     AMV(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)
      REAL*8    TAUX(NXYDIM)       ,   TAUY(NXYDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      HY(NXYDIM)

      REAL*8     TKE(NXYDIM, NZDIM)
      REAL*8  AHVBAK(NXYDIM, NZDIM)
      SAVE TKE, AHVBAK
      REAL*8  AHV03D(NXYDIM,NZDIM)
      REAL*8  ALPLAT(NXYDIM), C0LAT(NXYDIM)
      SAVE AHV03D, ALPLAT, C0LAT

      REAL*8    DRDZ(NXYDIM, NZDIM),  DUVDZ(NXYDIM, NZDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM), RZMSIG(NXYDIM, NZDIM)
      REAL*8   DEPTH(NXYDIM, NZDIM),     PR(NXYDIM, NZDIM)
      REAL*8       C(NXYDIM, NZDIM),   CDMP(NXYDIM, NZDIM)
      REAL*8  ADEFWD(NXYDIM, NZDIM),    RIT(NXYDIM, NZDIM)
      REAL*8     FEZ(NXYDIM, NZDIM),    TLS(NXYDIM, NZDIM)
      REAL*8  AMVTMP(NXYDIM, NZDIM),  DIFFZ(NXYDIM, NZDIM)
      REAL*8      AA(NXYDIM, NZDIM),     AB(NXYDIM, NZDIM)
      REAL*8      AC(NXYDIM, NZDIM),    ADE(NXYDIM, NZDIM)
      REAL*8     DML(NXYDIM)
c      COMMON /WORK/ DZSIG, RZMSIG, DEPTH, DRDZ, DUVDZ, TLS, C, FEZ
c      EQUIVALENCE(DZSIG(1, 1), ADE(1, 1))
c      EQUIVALENCE(DEPTH(1, 1), AA(1, 1))
c      EQUIVALENCE(FEZ(1, 1), AB(1, 1))
c      EQUIVALENCE(DUVDZ(1, 1), AC(1, 1))
c      EQUIVALENCE(C(1, 1), AMVTMP(1, 1))

      REAL*8 CC0(NZDIM), CC1(NZDIM), CC2(NZDIM)
      REAL*8 CC3(NZDIM), CC4(NZDIM), CC5(NZDIM), CC6(NZDIM)
      REAL*8 D0(NZDIM), D1(NZDIM), D2(NZDIM), D3(NZDIM), D4(NZDIM)
      REAL*8 D5(NZDIM), D6(NZDIM), D7(NZDIM), D8(NZDIM), D9(NZDIM)
      REAL*8 R0
      SAVE CC0, CC1, CC2, CC3, CC4, CC5, CC6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, R0

      REAL*8  P, ALPS, Q
      LOGICAL OFIRST, ODML(NXYDIM), OEOF
      SAVE OFIRST, OEOF

      REAL*8    CORT,  COR30,     PI,  OMEGA,    BFQ
      REAL*8   SLATS,  SLATN,  NLATS,  NLATN,    LAT
      REAL*8     QTL,   ZQTL, DMLTMP
      REAL*8      P1,     P2
      REAL*8      TL,     SL
      REAL*8      RL,    RLU
      REAL*8    DUDZ,   DVDZ
      REAL*8    FKLS,     RI
      REAL*8   AVRTX,  AVRTY
      REAL*8      FC
      INTEGER     IJ,      K,   IITR
      INTEGER     KT
      INTEGER  IFPAR,  JFPAR

      REAL*8  SM0, C0, PR0, SG
      REAL*8    RAHV(NZ)
      REAL*8    AMV0(NZ),   AHV0(NZ)
      REAL*8     EPS,     Z0,   ALPH,  CFTKE
      REAL*8   ALPHC, BETA,  BETAC,  PR1, PRMAX, C0EQ
      REAL*8    AHVB, AMVMAX, AFLT, ALSC, RITC
      REAL*8   LATAL, LATEQ, LATC0, LATCEQ, AL, ALEQ
      INTEGER     MZ,   NITR
      LOGICAL   OALLAT,  OC0LAT
      SAVE      AMV0, AHV0, EPS, Z0, ALPH, CFTKE, MZ, NITR
      SAVE      AMVMAX, AHVB, AFLT, ALSC, RITC, ALPHC, BETA, BETAC
      SAVE      SM0, C0, C0EQ, PR0, SG, PR1, PRMAX, OALLAT, OC0LAT
      SAVE      RAHV
      NAMELIST /NMVISV/ AMV0
      NAMELIST /NMDIFV/ AHV0
      NAMELIST /NMDFRE/ RAHV
      NAMELIST /NMDVNK/ EPS, Z0, ALPH, ALPHC, BETA, BETAC, 
     &                  CFTKE, MZ, NITR, PR0, PR1, PRMAX, 
     &                  SM0, C0, C0EQ, SG, AHVB, AMVMAX, AFLT,
     &                  ALSC, RITC, OALLAT, AL, ALEQ, LATAL, LATEQ,
     &                  OC0LAT, LATC0, LATCEQ
      DATA AMV0 / NZ*0.D0 /
      DATA AHV0 / NZ*0.D0 /
      DATA RAHV / NZ*0.D0 /

      DATA EPS / 1.D-5 /
      DATA Z0, CFTKE / 1.D2, 1.D2 /
      DATA ALPH, ALPHC / 3.0D0, -999.0D0 /
      DATA BETA, BETAC / 0.0D0, -999.0D0 /
      DATA MZ, NITR / NZ, 1 /
      DATA SM0 / 0.39D0 /
      DATA C0  / 0.06D0 /
      DATA PR0, PR1, PRMAX / 0.8D0, 0.5D0, 20.0D0 /
      DATA SG  / 1.95D0 /
      DATA AHVB / 0.1D0 /
      DATA AMVMAX / 1000.0D0 /
      DATA AFLT / 1.0D0 /
      DATA ALSC / 999.0D9 /
      DATA RITC / 1.0D0 /

      DATA OALLAT / .FALSE. /
      DATA LATEQ, LATAL  /  5.0D0, 30.0D0 /
      DATA  ALEQ, AL     /  0.1D0,  1.0D0 /

      DATA OC0LAT / .FALSE. /
      DATA LATCEQ, LATC0  /  5.0D0, 15.0D0 /
      DATA C0, C0EQ  / 0.06D0, -999.0D0 /

      DATA OFIRST / .TRUE. /

      IF (OINIT) THEN
         DO K = 1, NZDIM
            DO IJ = 1, NXYDIM
               TKE(IJ, K) = EPS
            END DO
         END DO
#ifdef OPT_TRIPOLE
         CALL RSTADD(TKE, OEOF, NXDIM, NYDIM, NZDIM, 'TKE', 'OCN',
     &                                              1.D0,  0,  0 )
#else
         CALL RSTADD(TKE, OEOF, NXDIM, NYDIM, NZDIM, 'TKE', 'OCN')
#endif
         RETURN
      END IF

      IF (OFINAL) THEN
         CALL FINADD(TKE, NXDIM, NYDIM, NZDIM, 'TKE', 'OCN')
         RETURN
      END IF

      IF (OFIRST) THEN
ccc         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMVISV, END=99)
 99      WRITE(JFPAR, NMVISV)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMDIFV, END=98)
 98      WRITE(JFPAR, NMDIFV)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMDVNK, END=97)
 97      WRITE(JFPAR, NMDVNK)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMDFRE, END=96)
 96      WRITE(JFPAR, NMDFRE)

         IF (ALPHC .LT. 0.0D0) THEN
            ALPHC = ALPH
            WRITE(JFPAR, *) '   ALPHC --> ', ALPHC
         ENDIF
         IF (BETAC .LT. 0.0D0) THEN
            BETAC = BETA
            WRITE(JFPAR, *) '   BETAC --> ', BETAC
         ENDIF

         CALL SECOEF(
     O            CC0(KSTR), CC1(KSTR), CC2(KSTR), CC3(KSTR),
     O            CC4(KSTR), CC5(KSTR), CC6(KSTR),
     O            D0(KSTR), D1(KSTR), D2(KSTR), D3(KSTR),  D4(KSTR),
     O            D5(KSTR), D6(KSTR), D7(KSTR), D8(KSTR),  D9(KSTR))

         R0 = rhoo * 1.D3

!---- reduce the background diffusion around EQ. ---
         PI = ATAN( 1.D0 )*4.D0
         OMEGA = 2.D0 * PI / 86400.D0
         COR30 = 2.D0 * OMEGA * SIN( PI*30.D0/180.D0 )
         BFQ   = 5.24D-3        ! /sec

         DO K = KSTR, KEND
            DO IJ = IJSTR, IJEND
               CORT=(COR(IJ)+COR(IJ+LW)+COR(IJ+LSW)+COR(IJ+LS))*0.25D0
               CORT=ABS(CORT)
               IF(CORT.GT.COR30) THEN
                  AHV03D(IJ,K)=AHV0(K-KSTR+1)
               ELSE
c               AHV03D(IJ,K)=AHV0(K-KSTR+1)
c     &              *CORT*ACOSH(BFQ/CORT)/COR30/ACOSH(BFQ/COR30)
                  AHV03D(IJ, K) = AHV0(K-KSTR+1) * CORT / COR30 *
     &                            LOG(  BFQ / CORT
     &                                + SQRT(BFQ**2 / CORT**2 - 1.D0))
     &                          / LOG(  BFQ / COR30
     &                                + SQRT(BFQ**2 / COR30**2 - 1.D0))
                  AHV03D(IJ,K)= MAX(AHV03D(IJ, K), 1.D-2)
                  AHV03D(IJ,K)= AHV03D(IJ,  K) * RAHV(K-KSTR+1) 
     &                        +AHV0(K-KSTR+1) * (1.0 - RAHV(K-KSTR+1)) 
               ENDIF
            END DO
         END DO

         IF (OALLAT) THEN
            WRITE(JFPAR, *)
     &        '   : DML is calculated following MY ',
     &        'with AL latitude dependency'
            SLATN = -LATEQ
            SLATS = -LATAL
            NLATS = +LATEQ
            NLATN = +LATAL
            DO IJ=IJSTR,IJEND
               CORT=(COR(IJ)+COR(IJ+LW)+COR(IJ+LSW)+COR(IJ+LS))*0.25D0
               IF ( CORT/(2.D0*OMEGA) .GT.  1.D0 ) THEN
                  CORT =   2.D0 * OMEGA
               END IF
               IF ( CORT/(2.D0*OMEGA) .LT. -1.D0 ) THEN
                  CORT = - 2.D0 * OMEGA
               END IF

               LAT = ASIN( CORT/2.D0/OMEGA ) * 180.D0 / PI
               IF(LAT.GE.NLATN.or.LAT.LE.SLATS) THEN
                  ALPLAT(IJ) = AL
               ELSEIF(LAT.GE.SLATS.and.LAT.LT.SLATN) THEN
                  ALPLAT(IJ)=(AL*(SLATN-LAT)
     &                 +ALEQ*(LAT-SLATS))/(SLATN-SLATS)
               ELSEIF(LAT.GT.NLATS.and.LAT.LE.NLATN) THEN
                  ALPLAT(IJ)=(AL*(LAT-NLATS)
     &                 +ALEQ*(NLATN-LAT))/(NLATN-NLATS)
               ELSE
                  ALPLAT(IJ) = ALEQ
               ENDIF
            ENDDO
         ELSE
            WRITE(JFPAR, *) 
     &        '   : DML is calculated following Noh et al.(2002)'
            DO IJ=IJSTR, IJEND
               ALPLAT(IJ) = AL  !! DUMMY
            ENDDO
         ENDIF

!---- increasing TKE dissipation coefficient around EQ.
         IF (OC0LAT) THEN
            WRITE(JFPAR, *)
     &        '   : C0 around the equator is enhanced'
            SLATN = -LATCEQ
            SLATS = -LATC0
            NLATS = +LATCEQ
            NLATN = +LATC0
            DO IJ=IJSTR,IJEND
               CORT=(COR(IJ)+COR(IJ+LW)+COR(IJ+LSW)+COR(IJ+LS))*0.25D0
               IF ( CORT/(2.D0*OMEGA) .GT.  1.D0 ) CORT =   2.D0 * OMEGA
               IF ( CORT/(2.D0*OMEGA) .LT. -1.D0 ) CORT = - 2.D0 * OMEGA
               LAT = ASIN( CORT/2.0/OMEGA ) * 180.D0 / PI
               IF(LAT.GE.NLATN.or.LAT.LE.SLATS) THEN
                  C0LAT(IJ)  = C0
               ELSEIF(LAT.GE.SLATS.and.LAT.LT.SLATN) THEN
                  C0LAT(IJ)=(C0*(SLATN-LAT)
     &                 +C0EQ*(LAT-SLATS))/(SLATN-SLATS)
               ELSEIF(LAT.GT.NLATS.and.LAT.LE.NLATN) THEN
                  C0LAT(IJ)=(C0*(LAT-NLATS)
     &                 +C0EQ*(NLATN-LAT))/(NLATN-NLATS)
               ELSE
                  C0LAT(IJ)  = C0EQ
               ENDIF
            ENDDO
         ELSE
            DO IJ=IJSTR, IJEND
               C0LAT(IJ) = C0
            ENDDO         
         ENDIF

      END IF

C     -- second step of Euler-Backward sheme --
      IF ( ( ITST .EQ. 1 ).AND.( IEULER .EQ. 2 ) ) THEN
         RETURN
      ENDIF

      IF (OFIRST) THEN
         DO K = 1, NZDIM
            DO IJ = 1, NXYDIM
               AHVBAK(IJ, K) = AHV(IJ, K)
            END DO
         END DO
      END IF

      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DEPTH (IJ, K) = 0.D0
            FEZ   (IJ, K) = 0.D0
         END DO
      END DO
      DO K = KSTR, KSTR+KZ-1
         DO IJ = 1, NXYDIM
            DZSIG (IJ, K) = (HY(IJ) + ZBOT) * DS(K)
            RZMSIG(IJ, K) = 1.D0 / (HY(IJ) + ZBOT) / DSM(K)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = 1, NXYDIM
            DZSIG (IJ, K) = DZ(IJ, K)
            RZMSIG(IJ, K) = 1.D0 / DZM(IJ, K)
         END DO
      END DO
      DO K = KSTR+1, KEND
         DO IJ = 1, NXYDIM
            DEPTH(IJ, K) = DEPTH(IJ, K-1) + DZSIG(IJ, K-1)
         END DO
      END DO

      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            TL = TY(IJ, K, 1)
            SL = TY(IJ, K, 2)
            P1 = CC0(K)
     &         + (CC1(K) + (CC2(K) + CC3(K) * TL) * TL) * TL
     &         + (CC4(K) + CC5(K) * TL + CC6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RL = P1 / P2

            TL = TY(IJ, K-1, 1)
            SL = TY(IJ, K-1, 2)
            P1 = CC0(K)
     &         + (CC1(K) + (CC2(K) + CC3(K) * TL) * TL) * TL
     &         + (CC4(K) + CC5(K) * TL + CC6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLU = P1 / P2

            DRDZ(IJ, K) = gravit / R0 * (RL - RLU) * RZMSIG(IJ, K)
            DRDZ(IJ, K) = MAX(DRDZ(IJ, K), 0.D0)
            DUDZ        = (  UY(IJ   , K-1) + UY(IJ+LW , K-1)
     &                     + UY(IJ+LS, K-1) + UY(IJ+LSW, K-1)
     &                     - UY(IJ   , K  ) - UY(IJ+LW , K  )
     &                     - UY(IJ+LS, K  ) - UY(IJ+LSW, K  )) *
     &                    0.25D0 * RZMSIG(IJ, K)
            DVDZ        = (  VY(IJ   , K-1) + VY(IJ+LW , K-1)
     &                     + VY(IJ+LS, K-1) + VY(IJ+LSW, K-1)
     &                     - VY(IJ   , K  ) - VY(IJ+LW , K  )
     &                     - VY(IJ+LS, K  ) - VY(IJ+LSW, K  )) *
     &                    0.25D0 * RZMSIG(IJ, K)
            DUVDZ(IJ, K) = DUDZ * DUDZ + DVDZ * DVDZ
         END DO
      END DO

      IF (OALLAT) THEN
c        -- DML is calculated following MY with latitude dependency --
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            QTL   = 0.0D0
            ZQTL  = 0.0D0
            DO K = KSTR, KSTR+MZ-1
               QTL  = QTL  + SQRT ( 2.0D0 * TKE(IJ, K) ) * DZSIG(IJ, K)
               ZQTL = ZQTL + DEPTH(IJ, K) * 
     &                    SQRT ( 2.0D0 * TKE(IJ, K) ) * DZSIG(IJ, K)
            END DO
            DMLTMP = ALPLAT(IJ) * ZQTL / MAX( QTL, EPS )
            KT = MIN( KSTR+MZ-1, MAX(KSTR+1, NBOT(IJ)))
            DML(IJ) = MAX( MIN( DMLTMP, DEPTH( IJ, KT ) ), EPS )
         END DO
      ELSE
c        -- DML is calculated following Noh et al.(2002) --
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            KT = MIN( KSTR+MZ-1, MAX(KSTR+1, NBOT(IJ)))
            DML(IJ) = DEPTH( IJ, KT )
            ODML(IJ)=.TRUE.
         ENDDO
         DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
               IF(ODML(IJ)) THEN
               IF ( AHV(IJ, K) .LE. MAX(AHVB, AHV0(K-KSTR+1)) ) THEN
                  DML(IJ) = DEPTH(IJ, K)
                  ODML(IJ)=.FALSE.
               END IF
               END IF
         END DO
         END DO
      ENDIF

C     -- initialize TKE -- 
      IF (OFIRST .AND. OEOF) THEN
         DO K = KSTR+1, KSTR+MZ-1
            DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
C     ---- using definition of Kondoh et. (1978) for Pr first guess,
C          since the method of Noh et al. (2005) requires TKE. ---
               PR(IJ, K) =  MIN ( PRMAX, 
     &                      PR0 + 7.D0 * DRDZ(IJ, K)
     &                            / MAX( DUVDZ(IJ, K), EPS ) )
               FKLS = ckarm * (DEPTH(IJ, K) + Z0)
               TLS(IJ, K) = FKLS / (1.D0 + FKLS / DML(IJ))
               TKE(IJ, K) =
     &              (   AHV(IJ, K) * PR0
     &                + SQRT( AHV(IJ, K)**2 * PR0**2
     &                        + 4.0D0 * ALPH * SM0 * SM0 *
     &                          DRDZ(IJ, K) * TLS(IJ, K)**4 )
     &              ) * AHV(IJ, K) * PR0
     &                / ( 4.0D0 * SM0 * SM0 * TLS(IJ, K) * TLS(IJ, K) )
               TKE(IJ, K) = MAX( TKE(IJ, K) * AMFTZ(IJ, K), EPS)
            ENDDO
         ENDDO
      END IF

      DO IITR = 1, NITR

      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            FKLS = ckarm * (DEPTH(IJ, K) + Z0)
            TLS(IJ, K) = FKLS / (1.D0 + FKLS / DML(IJ))
            RIT(IJ, K) = DRDZ(IJ, K) * TLS(IJ, K) * TLS(IJ, K)
     &                   / TKE(IJ,K) * 0.5D0
C           ---- Prantle number by definition of Noh et al. (2005, GRL)
            PR(IJ, K) =  MIN ( PRMAX, 
     &                         PR0 * SQRT( 1.D0 + PR1 * RIT(IJ, K) )  )
c            PR(IJ, K) =  MIN ( PRMAX, 
c     &                         PR0 + PR1 * DRDZ(IJ, K)
c     &                               / MAX( DUVDZ(IJ, K), EPS ) )
         ENDDO
      ENDDO
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            RI = RIT(IJ, K) * EXP( - MIN(300.0D0, 
     &                                  BETA / MAX( RIT(IJ, K), EPS )))
            Q  = SQRT( 2.D0 * TKE(IJ, K) )
            FC = TLS(IJ, K) * SM0 / SQRT( 1.0D0 + ALPH * RI )
            AMV(IJ, K) = MIN( AMVMAX, Q * FC )

C            --- Forward for shear and buoyancy ---
C            C(IJ, K)      =   4.D0 * C0 * TKE(IJ, K) / TLS(IJ, K) / RI
C            ADEFWD(IJ, K) = - 2.D0 * SM0 * TLS(IJ, K) / RI *
C     &                      ( DRDZ(IJ, K) / PR - DUVDZ(IJ, K)) *
C     &                      TKE(IJ, K) * AMFTZ(IJ, K)
C            --- Semi-implicit for shear and buoyancy term ---
C            CDMP(IJ, K) = ( 4.D0 * C0 * TKE(IJ, K) / TLS(IJ, K) / RI
C     &                    ) * AMFTZ(IJ, K)
            CDMP(IJ, K) = 
     &           ( C0LAT(IJ) * SQRT( 1.0D0 + ALPHC * RI )
     &                * 2.0D0 * Q / TLS(IJ, K) 
     &           ) * AMFTZ(IJ, K)
            C(IJ, K)    =
     &           ( FC * 2.0D0 / Q
     &                * (DRDZ(IJ, K) / PR(IJ, K) - DUVDZ(IJ, K)) 
     &           ) * AMFTZ(IJ, K)
            ADEFWD(IJ, K) =  0.0D0
         END DO
      END DO

C     -- semi-implicit : adjusted ALPS in each time step --
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            P = C(IJ, K) * DT
C           -- foward after initial-growth phase --
            IF ( (( RIT(IJ, K).LT.RITC)
     &             .AND.( C(IJ,K).LT.0.0D0) )
     &           .OR.( P .LT. -300.0D0 ) ) THEN
               ALPS  =  1.0D0
C           -- implicit for strong damping phase --
C            ELSEIF ( P .GE. 1.0D0 ) THEN
C               ALPS  =  0.0D0
C           -- semi-implicit for initial-growth and damping phase --
            ELSE
               P = SIGN(1.0D0, P) * MAX(1.0D-6, ABS(P))
               ALPS  =  ( -1.0D0 + (1.0D0+P) * EXP(-P) )
     &                  / ( P * ( EXP(-P) - 1.0D0 ) )
            ENDIF
            ADEFWD(IJ, K) = ADEFWD(IJ, K) 
     &                      - ALPS * C(IJ, K) * TKE(IJ, K)
            C(IJ, K)      = (1.0D0 - ALPS) * C(IJ, K)
     &                      + CDMP(IJ, K)
         END DO
      END DO

      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DIFFZ(IJ, K) = 0.D0
         END DO
      END DO
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         AVRTX = (  TAUX(IJ)    + TAUX(IJ+LW)
     &            + TAUX(IJ+LS) + TAUX(IJ+LSW)) * 0.25D0
         AVRTY = (  TAUY(IJ)    + TAUY(IJ+LW)
     &            + TAUY(IJ+LS) + TAUY(IJ+LSW)) * 0.25D0
         FEZ(IJ, KSTR+1) = ((AVRTX * AVRTX + AVRTY * AVRTY)
     &                      / rhoo / rhoo) ** 0.75D0 * CFTKE
      END DO
      DO K = KSTR+2, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            DIFFZ(IJ, K) = (AMV(IJ, K-1) + AMV(IJ, K)) *
     &                     0.5D0 / SG / DZSIG(IJ, K-1) *
     &                     AMSKT(IJ, K-1)
            FEZ(IJ, K) = DIFFZ(IJ, K) * (TKE(IJ, K-1) - TKE(IJ, K))
         END DO
      END DO
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            ADE(IJ, K) = (FEZ(IJ, K) - FEZ(IJ, K+1)) * RZMSIG(IJ, K)
     &                  - C(IJ, K) * TKE(IJ, K) + ADEFWD(IJ, K)
         END DO
      END DO

      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            AA(IJ, K) = - DT * DIFFZ(IJ, K) * RZMSIG(IJ, K)
            AC(IJ, K) = - DT * DIFFZ(IJ, K+1) * RZMSIG(IJ, K)
            AB(IJ, K) = 1.D0 - AA(IJ, K) - AC(IJ, K)
     &                 + DT * C(IJ, K)
         END DO
      END DO
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         AC(IJ, KSTR+1) = AC(IJ, KSTR+1) / AB(IJ, KSTR+1)
         ADE(IJ, KSTR+1) = ADE(IJ, KSTR+1) / AB(IJ, KSTR+1)
      END DO
      DO K = KSTR+2, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            FC = 1.D0 / (AB(IJ, K) - AA(IJ, K) * AC(IJ, K-1))
            AC(IJ, K) = AC(IJ, K) * FC
            ADE(IJ, K) = (ADE(IJ, K) - AA(IJ, K) * ADE(IJ, K-1)) * FC
         END DO
      END DO
      DO K = KSTR+MZ-2, KSTR+1, -1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            ADE(IJ, K) = ADE(IJ, K) - AC(IJ, K) * ADE(IJ, K+1)
         END DO
      END DO

C     -- time filter to TKE in case of strong damping --
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
c            IF (CDMP(IJ, K)*DT .GE. ALSC) THEN
c               TKE(IJ, K) = (  0.5D0 * TKE(IJ, K)
c     &                       + 0.5D0 * (TKE(IJ, K) + DT * ADE(IJ, K))
c     &                      ) * AMFTZ(IJ, K)
c            ELSE
               TKE(IJ, K) = ( TKE(IJ, K) + DT * ADE(IJ, K) 
     &                      ) * AMFTZ(IJ, K)
c            ENDIF
         END DO
      END DO
      DO K = KSTR, KEND
         DO IJ = 1, NXYDIM
            TKE(IJ, K) = MAX(TKE(IJ, K), EPS)
         END DO
      END DO

      END DO

      IF (OFIRST) THEN
         NITR = 1
         OFIRST = .FALSE.
      END IF


C     -- smoothing --
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND
            AMVTMP(IJ, K) = (  AMV(IJ   , K) + AMV(IJ+LE , K)
     &                       + AMV(IJ+LN, K) + AMV(IJ+LNE, K)
     &                      ) * 0.25D0 * AMFVZ(IJ, K)
         END DO
         DO IJ = IJSTR, IJEND
            AMV(IJ, K) = AMVTMP(IJ, K)
            AHV(IJ, K) = (  AMVTMP(IJ   , K) + AMVTMP(IJ+LW , K)
     &                    + AMVTMP(IJ+LS, K) + AMVTMP(IJ+LSW, K)) *
     &                   0.25D0 / PR(IJ, K)
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJSTR, IJEND
            AMV(IJ, K) = MIN(AMVMAX, MAX(AMV(IJ, K), AMV0(K-KSTR+1)))
            AHV(IJ, K) = MAX(AHV(IJ, K), AHV03D(IJ,K))
         END DO
      END DO

C     -- TIME FILTER --
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR, IJEND
            AHV(IJ, K) = (  AFLT * AHV(IJ, K) 
     &                     + ( 1.0D0 - AFLT ) * AHVBAK(IJ, K)
     &                   ) * AMFTZ(IJ, K)
         END DO
      END DO

C     -- MEMORY PREVIOUS AMV --
      DO K = KSTR+1, KSTR+MZ-1
         DO IJ = IJSTR, IJEND
            AHVBAK(IJ, K) = AHV(IJ, K)
         END DO
      END DO

#ifdef OPT_TRIPOLE
         CALL SHIFT1(   TKE,
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,      0,      0 )
#endif

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE VDIFFB(
     M                AMV,    AHV)
      use zocdim
      use zocmsk

c --- information -----------------------------------------------------
c
c  Vertical viscosity and diffusion coefficients for the bottom
c boundary layer.
c
c  HISTORY
c     '01.02.08  H.Hasumi
c
c ---------------------------------------------------------------------

      REAL*8     AMV(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      INTEGER     IJ
      INTEGER     KT,     KV
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  AMVBBL, AHVBBL
      SAVE AMVBBL, AHVBBL
      NAMELIST /NMBBDV/ AMVBBL, AHVBBL

      DATA OFIRST / .TRUE. /
      DATA AMVBBL, AHVBBL / 0.D0, 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMBBDV, END=99)
 99      WRITE(JFPAR, NMBBDV)
      END IF

      DO IJ = IJSTR, IJEND
         KT = NBOT(IJ)
         KV = NBOTV(IJ)
         AMV(IJ, KV) = MAX(AMVBBL, AMV(IJ, KV)) * AMSKVB(IJ)
     &               + AMV(IJ, KV) * (1.D0 - AMSKVB(IJ))
         AHV(IJ, KT) = MAX(AHVBBL, AHV(IJ, KT)) * AMSKTB(IJ)
     &               + AHV(IJ, KT) * (1.D0 - AMSKTB(IJ))
      END DO

      DO IJ = IJSTR, IJEND
         KT = NBOT(IJ)
         KV = NBOTV(IJ)
         AMV(IJ, KEND) = AMV(IJ, KV)
         AHV(IJ, KEND) = AHV(IJ, KT)
      END DO


      RETURN
      END
#endif
