      SUBROUTINE VDIFF(
     O               AMV,    AHV,
     I                UY,     VY,      R,   TAUX,   TAUY,
     I                TY,     HY)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Sea surface mixed layer model of Umlauf and Burchard (2003, JMR).
c  The implementation is based on Warner et al. (2005, OM).
c  Effects of surface wave breaking is also introduced 
c    after Carniel et al. (2009, OM).
c
c  HISTORY
c     '09.10.23  Y.Komuro: for COCO4.5
c     '10.04.09  Y.Komuro: many bug fixes + surface wave breaking
c     '10.04.13  Y.Komuro: bug fix (for non-BBL use etc.)
c     '12.01.19  Y.Komuro: bug fix (GH in unstable situation)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       R(NXYDIM, NZDIM)
      REAL*8     AMV(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)
      REAL*8    TAUX(NXYDIM)       ,   TAUY(NXYDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      HY(NXYDIM)

      REAL*8     TKE(NXYDIM, NZDIM),    PSI(NXYDIM, NZDIM)
      REAL*8  AHV03D(NXYDIM, NZDIM)
      REAL*8   CSAMV(NXYDIM, NZDIM),  CSAHV(NXYDIM, NZDIM)
      REAL*8  TAUAOX(NXYDIM), TAUAOY(NXYDIM)
      SAVE TKE, PSI, AHV03D, CSAMV, CSAHV, TAUAOX, TAUAOY

      REAL*8    AMVT(NXYDIM, NZDIM)
      REAL*8    DRDZ(NXYDIM, NZDIM),  DUVDZ(NXYDIM, NZDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM), RZMSIG(NXYDIM, NZDIM)
      REAL*8   DEPTH(NXYDIM, NZDIM), DEPTHM(NXYDIM, NZDIM)
      REAL*8   EPSIL(NXYDIM, NZDIM),    TLS(NXYDIM, NZDIM)
      REAL*8      SH(NXYDIM, NZDIM),     SM(NXYDIM, NZDIM)
      REAL*8   FWALL(NXYDIM, NZDIM)
      REAL*8     FEZ(NXYDIM, NZDIM),    RIT(NXYDIM, NZDIM)
      REAL*8      AA(NXYDIM, NZDIM),     AB(NXYDIM, NZDIM)
      REAL*8      AC(NXYDIM, NZDIM),    ADE(NXYDIM, NZDIM)
      REAL*8  ADEFWD(NXYDIM, NZDIM),   TKE0(NXYDIM, NZDIM)
      REAL*8  CDKDTL(NXYDIM, NZDIM),  DIFFZ(NXYDIM, NZDIM)
      REAL*8    CDMP(NXYDIM, NZDIM)
      REAL*8  UFRC2S(NXYDIM),        UFRC2B(NXYDIM)
      REAL*8  UFRC2O(NXYDIM),        UFRC2I(NXYDIM)
      REAL*8  UFRC3S(NXYDIM),        UFRC3O(NXYDIM)
      REAL*8  TAUBTM(NXYDIM)
      REAL*8      GH(NXYDIM, NZDIM),   GHUL(NXYDIM, NZDIM)
      REAL*8  Z0SF2D(NXYDIM),        PREPMX(NXYDIM)
      REAL*8   DPSI0(NXYDIM), SCNP3D(NXYDIM, NZDIM)
      REAL*8    TAOX(NXYDIM),   TAOY(NXYDIM)
      REAL*8  CTKEMN(NXYDIM, NZDIM), CTFILT(NXYDIM, NZDIM)
      COMMON /WORK/ DZSIG, RZMSIG, DRDZ, DUVDZ,
     &              DEPTH, DEPTHM, SCNP3D, FWALL

      REAL*8 CC0(NZDIM), CC1(NZDIM), CC2(NZDIM)
      REAL*8 CC3(NZDIM), CC4(NZDIM), CC5(NZDIM), CC6(NZDIM)
      REAL*8 D0(NZDIM), D1(NZDIM), D2(NZDIM), D3(NZDIM), D4(NZDIM)
      REAL*8 D5(NZDIM), D6(NZDIM), D7(NZDIM), D8(NZDIM), D9(NZDIM)
      REAL*8 R0
      SAVE CC0, CC1, CC2, CC3, CC4, CC5, CC6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, R0

      INTEGER NITR
      SAVE NITR

      LOGICAL OFIRST, OEOF
      SAVE OFIRST, OEOF

      REAL*8 CSFSH1, CSFSH2, CSFSM1, CSFSM2, CSFSM3
      REAL*8 CEPS, EEPST, EEPSP, CTLS
      REAL*8 CPSI30, CPSI3S, CPSIFE, CPSIBS, CPSIBB
      REAL*8 CPSILM, EPSILM
      REAL*8 SCNPSS, CSCNR, CSFTKW, CPSIFW
      REAL*8 CTAUBT
      SAVE CSFSH1, CSFSH2, CSFSM1, CSFSM2, CSFSM3
      SAVE CEPS, EEPST, EEPSP, CTLS
      SAVE CPSI30, CPSI3S, CPSIFE, CPSIBS, CPSIBB
      SAVE CPSILM, EPSILM
      SAVE SCNPSS, CSCNR, CSFTKW, CPSIFW
      SAVE CTAUBT

      REAL*8    CORT,  COR30,     PI,  OMEGA,    BFQ
      REAL*8    DUDZ,   DVDZ
      REAL*8      P1,     P2
      REAL*8      TL,     SL
      REAL*8      RL,    RLU
      REAL*8       Q,      P,   ALPS
      REAL*8      FC,  AVRTX,  AVRTY, AVRTOX, AVRTOY,  PSILM
      REAL*8  TKETMP, DSTWAL,  DZTMP,  RSCNP
      REAL*8   PREPS, PREPS0, ALPHCL
      REAL*8    CAIC,   CAIS
      INTEGER     IJ,      K,   IITR
      INTEGER  IFPAR,  JFPAR

      REAL*8    RAHV(NZ)
      REAL*8    AMV0(NZ),   AHV0(NZ)
      REAL*8  BTMFRC
      REAL*8   ESTRM, ESTRN, ESTRP
      REAL*8   CSFA1, CSFA2, CSFB1, CSFB2, CSFC2, CSFC3
      REAL*8  SFGHMN, SFGH0, SFGHCR
      REAL*8    CSFE, CMU0, CPSI1, CPSI2, CPSI3P, CPSI3N
      REAL*8  AMVMAX, AHVMAX, SCNTKE, SCNPSI
      REAL*8  TKEMIN, PSIMIN, RITC, Z0SFC, Z0BTM, EPSCMP
      REAL*8  CW, ALPHCH, Z0SFMN, ALPHCI, ATFILT
      INTEGER   NITR0, MZ
      LOGICAL OSFCWV, OSWNOI, OBTKEI
      SAVE      AMV0, AHV0, RAHV, BTMFRC,
     &          ESTRM, ESTRN, ESTRP,
     &          CSFA1, CSFA2, CSFB1, CSFB2, CSFC2, CSFC3,
     &          SFGHMN, SFGH0, SFGHCR,
     &          CSFE, CMU0, CPSI1, CPSI2, CPSI3P, CPSI3N,
     &          AMVMAX, AHVMAX, SCNTKE, SCNPSI,
     &          TKEMIN, PSIMIN, RITC, Z0SFC, Z0BTM, EPSCMP,
     &          CW, ALPHCH, Z0SFMN, ALPHCI, ATFILT,
     &          NITR0, MZ,
     &          OSFCWV, OSWNOI, OBTKEI
      NAMELIST /NMVISV/ AMV0
      NAMELIST /NMDIFV/ AHV0
      NAMELIST /NMDFRE/ RAHV
      NAMELIST /NMBTMF/ BTMFRC
      NAMELIST /NMDVGL/ ESTRM, ESTRN, ESTRP,
     &                  CSFA1, CSFA2, CSFB1, CSFB2, CSFC2, CSFC3,
     &                  SFGHMN, SFGH0, SFGHCR, CSFE,
     &                  CMU0, CPSI1, CPSI2, CPSI3P, CPSI3N,
     &                  AMVMAX, AHVMAX, SCNTKE, SCNPSI,
     &                  TKEMIN, PSIMIN, RITC, Z0SFC, Z0BTM,
     &                  NITR0, EPSCMP, MZ,
     &                  OSFCWV, CW, Z0SFMN, ALPHCH, OSWNOI, ALPHCI,
     &                  OBTKEI, ATFILT
      DATA AMV0 / NZ*0.D0 /
      DATA AHV0 / NZ*0.D0 /
      DATA RAHV / NZ*1.D0 /
      DATA BTMFRC / 0.D0 /
      DATA ESTRM, ESTRN, ESTRP / 1.0D0, -0.6666666666666666D0, 2.0D0 /
      DATA CSFA1, CSFA2 / 0.92D0, 0.74D0 /
      DATA CSFB1, CSFB2 / 16.6D0, 10.1D0 /
      DATA CSFC2, CSFC3 / 0.7D0, 0.2D0 /
      DATA SFGHMN, SFGH0, SFGHCR / -0.28D0, 0.0233D0, 0.01D0 /
      DATA CSFE / 1.0D0 /
      DATA CMU0 / 0.5544D0 /
      DATA CPSI1, CPSI2 / 1.0D0, 1.22D0 /
      DATA CPSI3P, CPSI3N / 1.0D0, 0.1D0 /
      DATA AMVMAX, AHVMAX / 1000.0D0, 1000.0D0 /
      DATA SCNTKE, SCNPSI / 0.8D0, 1.07D0 /
      DATA TKEMIN, PSIMIN / 7.6D-2, 1.0D-12 /
      DATA RITC / 1.0D0 /
      DATA Z0SFC, Z0BTM / 1.0D2, 1.0D2 /
      DATA NITR0 / 1 /
      DATA EPSCMP / 1.0D-15 /
      DATA MZ / NZ /
      DATA ATFILT / 0.0D0 /
      DATA OSFCWV, OSWNOI, OBTKEI / .FALSE., .FALSE., .FALSE. /
      DATA CW, Z0SFMN / 100.0D0, 1.0D2 /
      DATA ALPHCH, ALPHCI / 5.6D4, 1.4D3 /

      DATA OFIRST / .TRUE. /

      IF (OINIT) THEN
         DO K = 1, NZDIM
            DO IJ = 1, NXYDIM
               TKE(IJ, K) = TKEMIN
               PSI(IJ, K) = PSIMIN
            END DO
         END DO
#ifdef OPT_TRIPOLE
         CALL RSTADD(TKE, OEOF, NXDIM, NYDIM, NZDIM, 'TKE', 'OCN',
     &                                              1.D0,  0,  0 )
         CALL RSTADD(PSI, OEOF, NXDIM, NYDIM, NZDIM, 'PSI', 'OCN',
     &                                              1.D0,  0,  0 )
#else
         CALL RSTADD(TKE, OEOF, NXDIM, NYDIM, NZDIM, 'TKE', 'OCN')
         CALL RSTADD(PSI, OEOF, NXDIM, NYDIM, NZDIM, 'PSI', 'OCN')
#endif
         RETURN
      END IF

      IF (OFINAL) THEN
         CALL FINADD(TKE, NXDIM, NYDIM, NZDIM, 'TKE', 'OCN')
         CALL FINADD(PSI, NXDIM, NYDIM, NZDIM, 'PSI', 'OCN')
         RETURN
      END IF

      IF (OFIRST) THEN
CCC         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMVISV, END=99)
 99      WRITE(JFPAR, NMVISV)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMDIFV, END=98)
 98      WRITE(JFPAR, NMDIFV)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMBTMF, END=97)
 97      WRITE(JFPAR, NMBTMF)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMDVGL, END=96)
 96      WRITE(JFPAR, NMDVGL)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMDFRE, END=95)
 95      WRITE(JFPAR, NMDFRE)

         CALL SECOEF(
     O            CC0(KSTR), CC1(KSTR), CC2(KSTR), CC3(KSTR),
     O            CC4(KSTR), CC5(KSTR), CC6(KSTR),
     O            D0(KSTR), D1(KSTR), D2(KSTR), D3(KSTR),  D4(KSTR),
     O            D5(KSTR), D6(KSTR), D7(KSTR), D8(KSTR),  D9(KSTR))

         R0 = rhoo * 1.D3

         IF (.NOT. OSFCWV) THEN
            CW = 0.0D0
         END IF

         CSFSH1 = CSFA2 * (1.0D0 - 6.0D0 * CSFA1 / CSFB1)
         CSFSH2 = 3.0D0 * CSFA2 *
     &            (6.0D0 * CSFA1 + CSFB2 * (1.0D0 - CSFC3))
         CSFSM1 = CSFB1**(-1.0D0/3.0D0)
         CSFSM2 = 18.0D0 * CSFA1 * CSFA1
     &            + 9.0D0 * CSFA1 * CSFA2 * (1.0D0 - CSFC2)
         CSFSM3 = 9.0D0 * CSFA1 * CSFA2
         CEPS = CMU0**(3.0D0 + ESTRP/ESTRN)
         EEPST = 1.5D0 + ESTRM/ESTRN
         EEPSP = (-1.0D0) / ESTRN
         CTLS = CMU0**3.0D0
         CPSI30 = (CPSI3P + CPSI3N) * 0.5D0
         CPSI3S = (CPSI3P - CPSI3N) * 0.5D0
         CPSIFE = ESTRN * (CMU0**ESTRP)
         CPSIBS = CMU0**(ESTRP-2.0D0*ESTRM) 
         CPSIBB = CMU0**(ESTRP-2.0D0*ESTRM) * (ckarm*Z0BTM)**(ESTRN)
         CPSILM = SQRT(0.56D0) * (CMU0**(ESTRP/ESTRN))
         EPSILM = ESTRM/ESTRN + 0.5D0
         CSCNR = SQRT(1.5D0*SCNTKE)*CMU0/ckarm
         SCNPSS = ckarm*ckarm/(CPSI2*CMU0*CMU0) * 
     &          ( ESTRN*ESTRN - 4.0D0/3.0D0*CSCNR*ESTRN*ESTRM
     &          - 1.0D0/3.0D0*CSCNR*ESTRN
     &          + 2.0D0/9.0D0*ESTRM*CSCNR*CSCNR
     &          + 4.0D0/9.0D0*CSCNR*CSCNR*ESTRM*ESTRM )
         CSFTKW = SQRT(1.5D0*SCNTKE)*CMU0
         CPSIFW = ESTRM * SCNTKE * (CMU0**ESTRP)
         IF (OBTKEI) THEN
            CTAUBT = BTMFRC
         ELSE
            CTAUBT = 0.0D0
         END IF

         NITR = NITR0

         DO K = KSTR+1, KEND
            DO IJ = IJSTR-NXDIM-1, IJEND
               CSAMV(IJ, K) = (  AMFTZ(IJ   , K) + AMFTZ(IJ+LE , K)
     &                         + AMFTZ(IJ+LN, K) + AMFTZ(IJ+LNE, K)
     &                        ) * AMFVZ(IJ, K)
               IF (CSAMV(IJ, K) .EQ. 0.0D0) THEN
                  CSAMV(IJ, K) = 0.0D0
               ELSE
                  CSAMV(IJ, K) = 1.0D0 / CSAMV(IJ, K)
               END IF
            END DO
            DO IJ = IJSTR, IJEND
               CSAHV(IJ, K) = (  AMFVZ(IJ   , K) + AMFVZ(IJ+LW , K)
     &                         + AMFVZ(IJ+LS, K) + AMFVZ(IJ+LSW, K)
     &                        ) * AMFTZ(IJ, K)
               IF (CSAHV(IJ, K) .EQ. 0.0D0) THEN
                  CSAHV(IJ, K) = 0.0D0
               ELSE
                  CSAHV(IJ, K) = 1.0D0 / CSAHV(IJ, K)
               END IF
            END DO
         END DO

C --- reduce the background diffusion around EQ. ---
C
         PI = ATAN( 1.D0 )*4.D0
         OMEGA = 2.D0 * PI / 86400.D0
         COR30 =   2.D0 * OMEGA * SIN( PI*30.D0/180.D0 )
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
                  AHV03D(IJ,K)= MAX(AHV03D(IJ, K), 1.0D-2)
                  AHV03D(IJ,K)= AHV03D(IJ,  K) * RAHV(K-KSTR+1) 
     &                        +AHV0(K-KSTR+1) * (1.0 - RAHV(K-KSTR+1)) 
               ENDIF
            END DO
         END DO
      END IF

C     -- second step of Euler-Backward sheme --
      IF ( ( ITST .EQ. 1 ).AND.( IEULER .EQ. 2 ) ) THEN
         RETURN
      ENDIF

#ifdef OPT_BBL
         CALL RMMSKT
         CALL ADMKT1
#endif

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
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DEPTH (IJ, K) = 0.D0
            DEPTHM(IJ, K) = 0.0D0
         END DO
      END DO
      DO K = KSTR+1, KEND+1
         DO IJ = 1, NXYDIM
            DEPTH(IJ, K) = DEPTH(IJ, K-1) + DZSIG(IJ, K-1)
         END DO
      END DO
      DO K = KSTR, KSTR+KZ-1
         DO IJ = 1, NXYDIM
            DEPTHM(IJ, K) = DEPTHM(IJ, K-1) + (HY(IJ) + ZBOT) * DSM(K)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = 1, NXYDIM
            DEPTHM(IJ, K) = DEPTHM(IJ, K-1) + DZM(IJ, K)
         END DO
      END DO

      DO K = KSTR+1, KEND
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
c            DRDZ(IJ, K) = MAX(DRDZ(IJ, K), 0.D0)
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

C     Initial guess for TKE and PSI
c      (assume minimum TKE, neutral stratification(i.e., GH=0),
c        and shear production = disspiation)
c      NITR, N. of iteration, at the first step is set to 20.
      IF (OFIRST .AND. OEOF) THEN
c         DO K = KSTR+1, KEND
c            DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
c               TKE(IJ, K) = TKEMIN
c               TLS(IJ, K) = SQRT( (CMU0**3.0D0) * TKE(IJ, K)
c     &                    / ( SQRT(2.0D0)*CSFSM1*CSFE*DUVDZ(IJ, K) ) )
c               TLS(IJ, K) = MIN(TLS(IJ, K),
c     &            SQRT(0.56D0 * TKE(IJ, K) / MAX(DRDZ(IJ, K), EPSCMP)) )
c               PSI(IJ, K) =
c     &            (CMU0**ESTRP)*(TKE(IJ, K)**ESTRM)*(TLS(IJ, K)**ESTRN)
c            END DO
c         END DO
c
c         DO K = KSTR, KEND+1
c            DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
c               TKE(IJ, K) = MAX(TKEMIN, TKE(IJ, K)*AMFTZ(IJ, K))
c               PSI(IJ, K) = MAX(PSIMIN, PSI(IJ, K)*AMFTZ(IJ, K))
c            END DO
c         END DO
c
c         DO IJ = 1, NXYDIM
c            TAUBTM(IJ) = 0.0D0
c         END DO
cc         DO K = KSTR, KEND
cc            DO IJ = IJSTR-NXDIM-1-NXDIM-1, IJEND+NXDIM+1
cc               TAUBTM(IJ) = TAUBTM(IJ) + BTMFRC * rhoo *
cc     &              (UY(IJ, K) * UY(IJ, K) + VY(IJ, K) * VY(IJ, K))
cc     &              * AMSKB(IJ, K)
cc            END DO
cc         END DO
c         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
c            AVRTX = (  TAUX(IJ)    + TAUX(IJ+LW)
c     &           + TAUX(IJ+LS) + TAUX(IJ+LSW)) * 0.25D0
c            AVRTY = (  TAUY(IJ)    + TAUY(IJ+LW)
c     &           + TAUY(IJ+LS) + TAUY(IJ+LSW)) * 0.25D0
c            UFRC2S(IJ) = SQRT(AVRTX * AVRTX + AVRTY * AVRTY) / rhoo
c     &           * AMSKT(IJ, KSTR)
c            TKE(IJ, KSTR) = MAX( UFRC2S(IJ) / (CMU0*CMU0), TKEMIN)
cc            UFRC2B(IJ) =
cc     &           ( TAUBTM(IJ) + TAUBTM(IJ+LW)
cc     &           + TAUBTM(IJ+LS) + TAUBTM(IJ+LSW)) * 0.25D0
cc     &           / rhoo
cc            TKE(IJ, NBOT(IJ)+1) = MAX( UFRC2B(IJ) / (CMU0*CMU0), TKEMIN )
c         END DO

         NITR = 20
      END IF

      DO IITR = 1, NITR

c     TKE backup
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            TKE0(IJ, K) = TKE(IJ, K)
         END DO
      END DO

c     Dissipation EPSIL and turbulent length scale TLS.
c     An upper limit for TLS is also introduced (eq.(42))
c     Turbultent Richardson number RIT is also calculated.
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            EPSIL(IJ, K) = CEPS
     &                     * (TKE(IJ, K)**EEPST) * (PSI(IJ, K)**EEPSP)
            TLS(IJ, K) = CTLS * (TKE(IJ, K)**1.5D0) / EPSIL(IJ, K)
            IF (DRDZ(IJ, K) .GT. 0.0D0) THEN
               TLS(IJ, K) = MIN(TLS(IJ, K),
     &           SQRT( 0.56D0 * TKE(IJ, K) / MAX(DRDZ(IJ, K), EPSCMP)) )
            END IF
            RIT(IJ, K) = DRDZ(IJ, K) * TLS(IJ, K) * TLS(IJ, K)
     &                   / TKE(IJ,K) * 0.5D0
         END DO
      END DO
      
c     Kantha and Clayson (1994) quasi-equilibrium stability function
c      Eq.(33) of Warner et al. (2005, OM) has TYPOGRAPHICAL ERRORs;
c       - It should be Gh = Ghul - ( (Ghul-Ghcr)**2 / (Ghul+Gh0-2Ghcr) )
c       - SFGHCR should be 0.01 instead of 0.02.
c      See eq. (19) of Burchard et al. (1999, JMS) for the correct formula
c      (note that the definitions of the buoyancy parameter are different
c       between the two papers).
c
c     '12.01.19: bug fix (due to the typograpical error dscribed above)
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            GHUL(IJ, K) = (-0.5D0) *
     &               DRDZ(IJ, K) * TLS(IJ, K) * TLS(IJ, K) / TKE(IJ, K)
            IF (GHUL(IJ, K) .GT. SFGHCR) THEN
               GH(IJ, K) = GHUL(IJ, K) - ((GHUL(IJ, K) - SFGHCR)**2.0D0)
     &           / (GHUL(IJ, K) + SFGH0 - 2.0D0 * SFGHCR)
            ELSE
               GH(IJ, K) = MAX(SFGHMN, GHUL(IJ, K))
            ENDIF
            SH(IJ, K) = CSFSH1 / (1.0D0 - CSFSH2 * GH(IJ, K))
            SM(IJ, K) = (CSFSM1 + CSFSM2 * SH(IJ, K) * GH(IJ, K))
     &                  / (1.0D0 - CSFSM3 * GH(IJ, K))
         END DO
      END DO

c     Vertical eddy viscosity at T-grid AMVT and diffusivity AHV
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            Q = SQRT( 2.0D0 * TKE(IJ, K) )
            AMVT(IJ, K) = MIN( AMVMAX, 
     &                 CSFE * Q * TLS(IJ, K) * SM(IJ, K) )
            AHV(IJ, K) = MIN( AHVMAX,
     &                 CSFE * Q * TLS(IJ, K) * SH(IJ, K) )
         END DO
      END DO

c     Wall function FWALL is just a dummy
      DO K = KSTR-1, KEND+1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            FWALL(IJ, K) = 1.0D0
         END DO
      END DO
      
c     frictional velocity and surface roughness
      DO IJ = IJSTR-NXDIM-1-NXDIM-1, IJEND+NXDIM+1
            TAUBTM(IJ) = 0.0D0
      END DO
      DO K = KSTR, KEND
         DO IJ = IJSTR-NXDIM-1-NXDIM-1, IJEND+NXDIM+1
            TAUBTM(IJ) = TAUBTM(IJ) + CTAUBT * rhoo *
     &           (UY(IJ, K) * UY(IJ, K) + VY(IJ, K) * VY(IJ, K))
     &           * AMSKB(IJ, K)
         END DO
      END DO
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         AVRTX = (  TAUX(IJ)    + TAUX(IJ+LW)
     &            + TAUX(IJ+LS) + TAUX(IJ+LSW)) * 0.25D0
         AVRTY = (  TAUY(IJ)    + TAUY(IJ+LW)
     &            + TAUY(IJ+LS) + TAUY(IJ+LSW)) * 0.25D0
         UFRC2S(IJ) = SQRT(AVRTX * AVRTX + AVRTY * AVRTY) / rhoo
     &                * AMSKT(IJ, KSTR)
         IF (OSWNOI) THEN
            AVRTOX = (TAUAOX(IJ)    + TAUAOX(IJ+LW)
     &               + TAUAOX(IJ+LS) + TAUAOX(IJ+LSW)) * 0.25D0
            AVRTOY = (  TAUAOY(IJ)    + TAUAOY(IJ+LW)
     &               + TAUAOY(IJ+LS) + TAUAOY(IJ+LSW)) * 0.25D0
            UFRC2O(IJ) = SQRT(AVRTOX * AVRTOX + AVRTOY * AVRTOY) / rhoo
     &                * AMSKT(IJ, KSTR)
            UFRC2I(IJ) = SQRT( (AVRTX - AVRTOX)**2.0D0
     &                    + (AVRTY - AVRTOY)**2.0D0 ) / rhoo
     &                * AMSKT(IJ, KSTR)
         ELSE
            UFRC2O(IJ) = UFRC2S(IJ)
            UFRC2I(IJ) = 0.0D0
         END IF
         UFRC3S(IJ) = UFRC2S(IJ) * SQRT(UFRC2S(IJ))
         UFRC3O(IJ) = UFRC2O(IJ) * SQRT(UFRC2O(IJ))
         UFRC2B(IJ) =
     &        ( TAUBTM(IJ) + TAUBTM(IJ+LW)
     &        + TAUBTM(IJ+LS) + TAUBTM(IJ+LSW)) * 0.25D0
     &          / rhoo
         IF (OSFCWV) THEN
            Z0SF2D(IJ) = MAX(
     &           ALPHCH/gravit*UFRC2O(IJ)+ALPHCI/gravit*UFRC2I(IJ),
     &           Z0SFMN )
         ELSE
            Z0SF2D(IJ) = Z0SFC
         ENDIF
      END DO

c     Schmidt number for PSI
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            SCNP3D(IJ, K) = SCNPSI
         END DO
      END DO
      IF (OSFCWV) THEN
         DO IJ = 1, NXYDIM
            PREPMX(IJ) = EPSCMP
         END DO
         
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            PREPS0 = 0.0D0
            DO K = KSTR+1, NBOT(IJ)
               IF (PREPMX(IJ) .LT. 1.0D0) THEN
                  PREPS =
     &                 AMVT(IJ, K) * DUVDZ(IJ, K)
     &                 / MAX(EPSIL(IJ, K), EPSCMP)
                  PREPMX(IJ) = MAX(PREPMX(IJ), PREPS)
                  IF (PREPS .GE. 1.0D0) THEN
                     DPSI0(IJ) = DEPTH(IJ, K-1) + 
     &                    DZSIG(IJ, K-1) * (1.0D0-PREPS0)/(PREPS-PREPS0)
                  END IF
                  PREPS0 = PREPS
               END IF
            END DO
            IF (PREPMX(IJ) .LT. 1.0D0) THEN
               DPSI0(IJ) = DEPTH(IJ, NBOT(IJ)+1)
     &                     / MAX(PREPMX(IJ), EPSCMP)
            END IF
         END DO

         DO K = KSTR, KEND
            DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
               RSCNP = MIN(DEPTHM(IJ, K-1)/MAX(DPSI0(IJ),EPSCMP), 1.0D0)
               SCNP3D(IJ, K) = 
     &           (1.0D0-RSCNP)*SCNPSS*2.0D0/(FWALL(IJ,K-1)+FWALL(IJ,K))
     &           + RSCNP*SCNPSI
            END DO
         END DO
      ELSE
         DO K = KSTR, KEND+1
            DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
               SCNP3D(IJ, K) = SCNPSI
            END DO
         END DO
      END IF

c     TKE equation
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            CDKDTL(IJ, K) =   ! gamma = -(P+B)/k (temporary)
     &           ( AHV(IJ, K) * DRDZ(IJ, K)
     &           - AMVT(IJ, K) * DUVDZ(IJ, K) )  
     &             / TKE(IJ, K) * AMFTZ(IJ, K)
            CDMP(IJ, K) =     ! delta = epsil/k
     &           EPSIL(IJ, K) / TKE(IJ, K) * AMFTZ(IJ, K)
            ADEFWD(IJ, K) = 
     &         - ( AHV(IJ, K) * DRDZ(IJ, K) - AMVT(IJ, K) * DUVDZ(IJ, K) 
     &             + EPSIL(IJ, K)) * AMFTZ(IJ, K)
         END DO
      END DO

C     -- semi-implicit : adjusted ALPS in each time step --
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            P = CDKDTL(IJ, K) * DT   ! gamma * dt
C           -- foward after initial-growth phase --
            IF ( (( RIT(IJ, K).LT.RITC)
     &             .AND.( CDKDTL(IJ,K).LT.0.0D0) )
     &           .OR.( P .LT. -300.0D0 ) ) THEN
               ALPS  =  1.0D0
cc           -- implicit for strong damping phase --
c            ELSEIF ( P .GE. 1.0D0 ) THEN
c               ALPS  =  0.0D0
cc           -- semi-implicit for initial-growth and damping phase --
            ELSE
               P = SIGN(1.0D0, P) * MAX(1.0D-6, ABS(P))
               ALPS  =  ( -1.0D0 + (1.0D0+P) * EXP(-P) )
     &                  / ( P * ( EXP(-P) - 1.0D0 ) )
            ENDIF
            CDKDTL(IJ, K) = (1.0D0 - ALPS) * CDKDTL(IJ, K)
     &                    + CDMP(IJ, K)    ! (1-mu)*gamma + delta
         END DO
      END DO

c     --- boundary condition for TKE ---
c      "DIFFZ(IJ, KSTR+1)=0"  & "DIFFZ(IJ, NBOT(IJ)+1)=0"
c       ==>  "AA(IJ, KSTR+1)=0" & "AC(IJ, NBOT(IJ))=0".
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DIFFZ(IJ, K) = 0.D0
            FEZ(IJ, K) = 0.D0
         END DO
      END DO
      DO K = KSTR+2, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            DIFFZ(IJ, K) = (AMVT(IJ, K-1) + AMVT(IJ, K)) *
     &                     0.5D0 / SCNTKE / DZSIG(IJ, K-1) *
     &                     AMSKT(IJ, K)
            FEZ(IJ, K) = DIFFZ(IJ, K) * (TKE(IJ, K-1) - TKE(IJ, K))
         END DO
      END DO

c     TKE flux by surface wave breaking, after Carniel et al.(2009)
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         FEZ(IJ, KSTR+1) = CW * UFRC3O(IJ)
      END DO
      
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            ADE(IJ, K) = (FEZ(IJ, K) - FEZ(IJ, K+1)) * RZMSIG(IJ, K)
     &                  + ADEFWD(IJ, K)
         END DO
      END DO

      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            AA(IJ, K) = - DT * DIFFZ(IJ, K) * RZMSIG(IJ, K)
            AC(IJ, K) = - DT * DIFFZ(IJ, K+1) * RZMSIG(IJ, K)
            AB(IJ, K) = 1.D0 - AA(IJ, K) - AC(IJ, K)
     &                 + DT * CDKDTL(IJ, K)
         END DO
      END DO
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         AC(IJ, KSTR+1) = AC(IJ, KSTR+1) / AB(IJ, KSTR+1)
         ADE(IJ, KSTR+1) = ADE(IJ, KSTR+1) / AB(IJ, KSTR+1)
      END DO
      DO K = KSTR+2, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            FC = 1.D0 / (AB(IJ, K) - AA(IJ, K) * AC(IJ, K-1))
            AC(IJ, K) = AC(IJ, K) * FC
            ADE(IJ, K) = (ADE(IJ, K) - AA(IJ, K) * ADE(IJ, K-1)) * FC
         END DO
      END DO
      DO K = KEND-1, KSTR+1, -1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            ADE(IJ, K) = ADE(IJ, K) - AC(IJ, K) * ADE(IJ, K+1)
         END DO
      END DO

      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            TKE(IJ, K) = (  ATFILT * TKE(IJ, K)
     &                    + (1.0D0 - ATFILT) *
     &                      (TKE(IJ, K) + DT * ADE(IJ, K))
     &                   ) * AMFTZ(IJ, K)
         END DO
      END DO

c     surface and bottom TKE estimation
c     for surface TKE, surface wave breaking effect is accounted
c      (after Carniel et al. (2009)).
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         TKE(IJ, KSTR) = 1.0D0 / (CMU0*CMU0) *
     &        (UFRC3S(IJ) + UFRC3O(IJ) * CSFTKW * CW)**(2.0D0/3.0D0)
         TKE(IJ, NBOT(IJ)+1) = UFRC2B(IJ) / (CMU0*CMU0)
      END DO

      DO K = KSTR, KEND+1
         DO IJ = 1, NXYDIM
            IF (TKE(IJ, K) .LT. TKEMIN) THEN
               CTKEMN(IJ, K) = 1.0D0
            ELSE
               CTKEMN(IJ, K) = 0.0D0
            END IF
            TKE(IJ, K) = MAX(TKE(IJ, K), TKEMIN)
         END DO
      END DO

c     GLS quantity PSI equation
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            CDKDTL(IJ, K) =   ! gamma = -(c1*P+c3*B)/k (temporary)
     &           ( ( CPSI30 - CPSI3S * SIGN(1.0D0, DRDZ(IJ, K)) )
     &               * AHV(IJ, K) * DRDZ(IJ, K)
     &             - CPSI1 * AMVT(IJ, K) * DUVDZ(IJ, K) )  
     &               / TKE0(IJ, K) * AMFTZ(IJ, K)
            CDMP(IJ, K) =     ! delta = c2*epsil*Fwall/k
     &           CPSI2 * FWALL(IJ, K) * EPSIL(IJ, K) / TKE0(IJ, K)
     &           * AMFTZ(IJ, K)
            ADEFWD(IJ, K) =
     &           - ( CDKDTL(IJ, K) + CDMP(IJ, K) ) * PSI(IJ, K)
     &           * AMFTZ(IJ, K)
         END DO
      END DO

C     -- semi-implicit : adjusted ALPS in each time step --
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            P = CDKDTL(IJ, K) * DT   ! gamma * dt
C           -- foward after initial-growth phase --
            IF ( (( RIT(IJ, K).LT.RITC)
     &             .AND.( CDKDTL(IJ, K).LT.0.0D0) )
     &           .OR.( P .LT. -300.0D0 ) ) THEN
               ALPS  =  1.0D0
cc           -- implicit for strong damping phase --
c            ELSEIF ( P .GE. 1.0D0 ) THEN
c               ALPS  =  0.0D0
cc           -- semi-implicit for initial-growth and damping phase --
            ELSE
               P = SIGN(1.0D0, P) * MAX(1.0D-6, ABS(P))
               ALPS  =  ( -1.0D0 + (1.0D0+P) * EXP(-P) )
     &                  / ( P * ( EXP(-P) - 1.0D0 ) )
            ENDIF
            CDKDTL(IJ, K) = (1.0D0 - ALPS) * CDKDTL(IJ, K)
     &                    + CDMP(IJ, K)  ! (1-mu)*gamma + delta
         END DO
      END DO

c     --- boundary condition for PSI ---
c     "DIFFZ(IJ, KSTR+1)=0"  & "DIFFZ(IJ, NBOT(IJ)+1)=0"
c       ==> "AA(IJ, KSTR+1)=0" & "AC(IJ, NBOT(IJ))=0"
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DIFFZ(IJ, K) = 0.D0
            FEZ(IJ, K) = 0.D0
         END DO
      END DO
      DO K = KSTR+2, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            DIFFZ(IJ, K) = (AMVT(IJ, K-1) + AMVT(IJ, K)) *
     &                     0.5D0 / SCNP3D(IJ, K-1) / DZSIG(IJ, K-1) *
     &                     AMSKT(IJ, K)
            FEZ(IJ, K) = DIFFZ(IJ, K) * (PSI(IJ, K-1) - PSI(IJ, K))
         END DO
      END DO

c      FEZ(IJ, KSTR+1) & FEZ(IJ, NBOT(IJ)+1)
c      Eq.(54) of Warner et al. (2005, OM) has TYPOGRAPHICAL ERROR;
c        k^{n} must be \kappa^{n}.
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         K = KSTR+1
         TKETMP = (TKE0(IJ, K-1) + TKE0(IJ, K))*0.5D0
         DSTWAL = Z0SF2D(IJ) + DZSIG(IJ, K-1)*0.5D0
         FEZ(IJ, K) = 
     &        - CPSIFE / SCNP3D(IJ, K-1)
     &          * AMVT(IJ, K)
     &          * ( TKETMP**ESTRM )
     &          * ( ckarm**ESTRN )
     &          * ( DSTWAL**(ESTRN-1.0D0) )
     &        - CPSIFW / SCNP3D(IJ, K-1)
     &          * ( TKETMP**(ESTRM-1.0D0) )
     &          * ( ckarm**ESTRN )
     &          * ( DSTWAL**ESTRN )
     &          * CW * UFRC3O(IJ)
         K = NBOT(IJ)+1
         TKETMP = (TKE0(IJ, K-1) + TKE0(IJ, K))*0.5D0
         DSTWAL = Z0BTM + DZSIG(IJ, K-1)*0.5D0
         FEZ(IJ, K) = 
     &        CPSIFE / SCNP3D(IJ, K-1)
     &          * AMVT(IJ, K-1)
     &          * ( TKETMP**ESTRM )
     &          * ( ckarm**ESTRN )
     &          * ( DSTWAL**(ESTRN-1.0D0) )
      END DO
      
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            ADE(IJ, K) = (FEZ(IJ, K) - FEZ(IJ, K+1)) * RZMSIG(IJ, K)
     &                  + ADEFWD(IJ, K)
         END DO
      END DO

      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            AA(IJ, K) = - DT * DIFFZ(IJ, K) * RZMSIG(IJ, K)
            AC(IJ, K) = - DT * DIFFZ(IJ, K+1) * RZMSIG(IJ, K)
            AB(IJ, K) = 1.D0 - AA(IJ, K) - AC(IJ, K)
     &                + DT * CDKDTL(IJ, K)
         END DO
      END DO
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         AC(IJ, KSTR+1) = AC(IJ, KSTR+1) / AB(IJ, KSTR+1)
         ADE(IJ, KSTR+1) = ADE(IJ, KSTR+1) / AB(IJ, KSTR+1)
      END DO
      DO K = KSTR+2, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            FC = 1.D0 / (AB(IJ, K) - AA(IJ, K) * AC(IJ, K-1))
            AC(IJ, K) = AC(IJ, K) * FC
            ADE(IJ, K) = (ADE(IJ, K) - AA(IJ, K) * ADE(IJ, K-1)) * FC
         END DO
      END DO
      DO K = KEND-1, KSTR+1, -1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            ADE(IJ, K) = ADE(IJ, K) - AC(IJ, K) * ADE(IJ, K+1)
         END DO
      END DO

      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            PSI(IJ, K) = (  ATFILT * PSI(IJ, K)
     &                    + (1.0D0 - ATFILT) *
     &                      (PSI(IJ, K) + DT * ADE(IJ, K))
     &                   ) * AMFTZ(IJ, K)
         END DO
      END DO

c     Surface and bottom PSI estimation
c     (This part is only for output and will not be referred here,
c       thus can be removed).
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         PSI(IJ, KSTR) = CPSIBS * (ckarm*Z0SF2D(IJ))**(ESTRN)
     &                          * (UFRC2S(IJ)**(2.0D0*ESTRM))
         PSI(IJ, NBOT(IJ)+1) = CPSIBB * (UFRC2S(IJ)**(2.0D0*ESTRM))
      END DO

c     Limit on PSI (PSIMIN, and eq.(43) if stable stratification)
      DO K = KSTR, KEND+1
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            IF (DRDZ(IJ, K) .GT. 0.0D0) THEN
               PSILM = ( CPSILM * (TKE(IJ, K)**EPSILM)
     &              / SQRT(MAX(DRDZ(IJ, K), EPSCMP)) ) ** ESTRN
               IF (ESTRN .GT. 0.0D0) THEN
                  PSI(IJ, K) = MAX(MIN(PSI(IJ, K), PSILM), PSIMIN)
               ELSE
                  PSI(IJ, K) = MAX(MAX(PSI(IJ, K), PSILM), PSIMIN)
               ENDIF
            ELSE
               PSI(IJ, K) = MAX(PSI(IJ, K), PSIMIN)
            END IF
         END DO
      END DO

      END DO

      IF (OFIRST) THEN
         NITR = NITR0
         OFIRST = .FALSE.
      END IF

      DO K = KSTR+1, KEND
         DO IJ = IJSTR, IJEND
            AMV(IJ, K) = ( AMVT(IJ    , K) * AMFTZ(IJ    , K)
     &                   + AMVT(IJ+LE , K) * AMFTZ(IJ+LE , K)
     &                   + AMVT(IJ+LN , K) * AMFTZ(IJ+LN , K)
     &                   + AMVT(IJ+LNE, K) * AMFTZ(IJ+LNE, K)
     &                   ) * CSAMV(IJ, K)
         END DO
      END DO

C     -- smoothing --
      DO K = KSTR+1, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND
            AMVT(IJ, K) = ( AHV(IJ    , K) * AMFTZ(IJ    , K)
     &                    + AHV(IJ+LE , K) * AMFTZ(IJ+LE , K)
     &                    + AHV(IJ+LN , K) * AMFTZ(IJ+LN , K)
     &                    + AHV(IJ+LNE, K) * AMFTZ(IJ+LNE, K)
     &                   ) * CSAMV(IJ, K)
         END DO
         DO IJ = IJSTR, IJEND
            AHV(IJ, K) = ( AMVT(IJ    , K) * AMFVZ(IJ    , K)
     &                   + AMVT(IJ+LW , K) * AMFVZ(IJ+LW , K)
     &                   + AMVT(IJ+LS , K) * AMFVZ(IJ+LS , K)
     &                   + AMVT(IJ+LSW, K) * AMFVZ(IJ+LSW, K)
     &                   ) * CSAHV(IJ, K)
         END DO
      END DO

      DO K = KSTR, KSTR+MZ-1
         DO IJ = IJSTR-NXDIM-1, IJEND
            AMV(IJ, K) = MIN(AMVMAX, MAX(AMV(IJ, K), AMV0(K-KSTR+1)))
            AHV(IJ, K) = MAX(AHV(IJ, K), AHV03D(IJ, K))
         END DO
      END DO
      DO K = KSTR+MZ, KEND
         DO IJ = IJSTR, IJEND
            AMV(IJ, K) = AMV0(K-KSTR+1)
            AHV(IJ, K) = AHV03D(IJ, K)
         END DO
      END DO

#ifdef OPT_BBL
         CALL RMMSKT
         CALL ADMKTB
#endif

c      CALL CHEKIN(TKE  ,'TKE',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(PSI  ,'PSI',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(SM   ,'SM',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(EPSIL,'EPSIL',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(TLS  ,'TLS',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(GHUL ,'GHUL',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(GH   ,'GH',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(DUVDZ,'DUVDZ',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(DRDZ, 'DRDZ',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(TY(1,1,1), 'TODV',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(CTKEMN, 'CTKEMN',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')

#ifdef OPT_TRIPOLE
         CALL SHIFT2(   TKE,    PSI,
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,      0,      0 )
#endif

      RETURN

c =====================================================================

      ENTRY PUTTAO(
     I               TAOX,   TAOY,
     I               CAIC,   CAIS)

      DO 800 IJ = IJVSTR, IJVEND
         TAUAOX(IJ) = ( TAOX(IJ) * CAIC
     &              - TAOY(IJ) * SIGN(CAIS, COR(IJ))) *
     &              AMSKV(IJ, KSTR)
         TAUAOY(IJ) = ( TAOY(IJ) * CAIC
     &              + TAOX(IJ) * SIGN(CAIS, COR(IJ))) *
     &              AMSKV(IJ, KSTR)
  800 CONTINUE

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M            TAUAOX, TAUAOY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,     -1,     -1 )
#else
      CALL SHIFT2(
     M              TAUAOX,   TAUAOY,
     D               NXDIM,    NYDIM,     1)
#endif
#else
      CALL STBCTU(
     M              TAUAOX,   TAUAOY)
#endif

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE VDIFFB(
     M                AMV,    AHV)

c --- information -----------------------------------------------------
c
c  Vertical viscosity and diffusion coefficients for the bottom
c boundary layer.
c
c  HISTORY
c     '01.02.08  H.Hasumi
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocmsk.F"

      REAL*8     AMV(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      REAL*8    AMVB(NXYDIM),   AHVB(NXYDIM)
      REAL*8     KTB(NXYDIM),    KVB(NXYDIM)

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
         KVB(IJ) = KV
         AMVB(IJ) = AMV(IJ, KV)
         AHV(IJ, KT) = MAX(AHVBBL, AHV(IJ, KT)) * AMSKTB(IJ)
     &               + AHV(IJ, KT) * (1.D0 - AMSKTB(IJ))
         KTB(IJ) = KT
         AHVB(IJ) = AHV(IJ, KT)
      END DO

      DO IJ = IJSTR, IJEND
         KT = NBOT(IJ)
         KV = NBOTV(IJ)
         AMV(IJ, KEND) = AMV(IJ, KV)
         AHV(IJ, KEND) = AHV(IJ, KT)
      END DO

c      CALL CHEKIN(AHV  ,'AHV3',
c     &            NX,     NY,     NZ, NXYZDM, 'OCN')      
c      CALL CHEKIN(AHVB  ,'AHVB',
c     &            NX,     NY,      1, NXYDIM, 'SFC')      
c      CALL CHEKIN(AMVB  ,'AMVB',
c     &            NX,     NY,      1, NXYDIM, 'SFC')      
c      CALL CHEKIN(KVB   ,'KVB',
c     &            NX,     NY,      1, NXYDIM, 'SFC')      
c      CALL CHEKIN(KTB   ,'KTB',
c     &            NX,     NY,      1, NXYDIM, 'SFC')      
      RETURN
      END
#endif
