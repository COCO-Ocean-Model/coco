      SUBROUTINE OVTSET(
     O                       R,
     I                       T)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Convective adjustment for the unstable water column
c
c  HISTORY
c     '99.08.16  H.Hasumi: from CCSR2-MASK
c     '00.05.30  H.Hasumi: change in the initialization of R
c     '00.12.14  H.Hasumi: for hybrid vertical coordinate
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '01.09.17  H.Hasumi: GAMMA for BBL
c     '01.12.07  H.Hasumi
c     '02.05.29  H.Nakano: tracer dimension
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: McDougall et al. (2003) eq. of state
c     '07.09.25  H.Hasumi: arguments of CHEKIN
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.02.23  Y.Komuro: add DDENST (only diagnosing R)
c
c ---------------------------------------------------------------------

      REAL*8       R(NXYDIM, NZDIM)
      REAL*8       T(NXYDIM, NZDIM, NTDIM)
      REAL*8       H(NXYDIM)

      REAL*8     TTL(NXYDIM, NTDIM),     W2(NXYDIM)
      REAL*8    CONV(NXYDIM, NZDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM)
      INTEGER    LUP(NXYDIM)
c      COMMON /WORK/ TTL, W2, CONV, DZSIG, LUP

      REAL*8 C0(NZDIM), C1(NZDIM), C2(NZDIM)
      REAL*8 C3(NZDIM), C4(NZDIM), C5(NZDIM), C6(NZDIM)
      REAL*8 D0(NZDIM), D1(NZDIM), D2(NZDIM), D3(NZDIM), D4(NZDIM)
      REAL*8 D5(NZDIM), D6(NZDIM), D7(NZDIM), D8(NZDIM), D9(NZDIM)
      REAL*8 ZT(NXYDIM, NZDIM)
      SAVE ZT
      SAVE C0, C1, C2, C3, C4, C5, C6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9

      REAL*8      P1,     P2
      REAL*8      TU,     TL,     SU,     SL
      REAL*8      RU,     RL
      INTEGER     IJ,      K,     KK,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8   GAMMA(NZ)
      NAMELIST /NMACCT/ GAMMA
      SAVE GAMMA
      DATA GAMMA / NZ*1.D0 /

      CALL REWNML(IFPAR, JFPAR)
      WRITE(JFPAR, *) '*** OVTSET ***'
      READ(IFPAR, NMACCT, END=99)
   99 CONTINUE

      CALL SECOEF(
     O            C0(KSTR), C1(KSTR), C2(KSTR), C3(KSTR),
     O            C4(KSTR), C5(KSTR), C6(KSTR),
     O            D0(KSTR), D1(KSTR), D2(KSTR), D3(KSTR),  D4(KSTR),
     O            D5(KSTR), D6(KSTR), D7(KSTR), D8(KSTR),  D9(KSTR))
      WRITE(JFPAR, *) '*** the equation of state determined ***'

      DO 31 K = 1, NZDIM
         DO 30 IJ = 1, NXYDIM
            R(IJ, K) = 0.D0
   30    CONTINUE
   31 CONTINUE
      DO 41 K = KSTR, KEND
         DO 40 IJ = IJTSTR, IJTEND
            TL = T(IJ, K, 1)
            SL = T(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            R(IJ, K) = P1 / P2 - 1.D3
   40    CONTINUE
   41 CONTINUE

      RETURN

c =====================================================================

      ENTRY OVTURN(
     O                  R,
     M                  T,    
     I                  H)

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      DO 51 K = 1, NZDIM
         DO 50 IJ = 1, NXYDIM
            CONV(IJ, K) = 0.D0
   50    CONTINUE
   51 CONTINUE

      DO K = KSTR, KSTR+KZ-1
         DO IJ = 1, NXYDIM
            DZSIG(IJ, K) = DS(K) * (H(IJ) + ZBOT) * GAMMA(K-KSTR+1)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = 1, NXYDIM
            DZSIG(IJ, K) = DZ(IJ, K) * GAMMA(K-KSTR+1)
         END DO
      END DO
#ifdef OPT_BBL
      DO IJ = IJTSTR, IJTEND
         K = NBOT(IJ)
         DZSIG(IJ, K) = DZ(IJ, K) * GAMMA(NZ) * AMSKTB(IJ)
     &                + DZSIG(IJ, K) * (1.D0 - AMSKTB(IJ))
      END DO
#endif

      DO IJ = 1, NXYDIM
         ZT(IJ, KSTR) = 0.D0
      END DO
      DO K = KSTR, KEND
         DO IJ = 1, NXYDIM
            ZT(IJ, K+1) = ZT(IJ, K) + DZSIG(IJ, K)
         END DO
      END DO

      DO N = 1, NTDIM
         DO 100 IJ = IJTSTR, IJTEND
            TTL(IJ, N) = T(IJ, KSTR, N) * DZSIG(IJ, KSTR)
            LUP(IJ) = KSTR
  100    CONTINUE
      END DO

      DO 1000 IJ = IJTSTR, IJTEND
         DO 200 K = KSTR+1, NBOT(IJ)
            TU = T(IJ, K-1, 1)
            SU = T(IJ, K-1, 2)
            TL = T(IJ, K, 1)
            SL = T(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TU) * TU) * TU
     &         + (C4(K) + C5(K) * TU + C6(K) * SU) * SU
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TU) * TU) * TU) * TU
     &         + (D5(K) + (D6(K) + D7(K) * TU * TU) * TU
     &                  + (D8(K) + D9(K) * TU * TU) * SQRT(SU)) * SU
            RU = P1 / P2
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RL = P1 / P2
            IF (RU .GT. RL) THEN
               W2 (IJ) = 1.D0 / (ZT(IJ, K+1) - ZT(IJ, LUP(IJ)))
               CONV(IJ, K) = 1.D0
               DO N = 1, NTDIM
                  TTL(IJ, N) = TTL(IJ, N)
     &                 + T(IJ, K, N) * DZSIG(IJ, K)
                  DO 410 KK = KSTR, K
                     IF (KK .GE. LUP(IJ)) THEN
                        T(IJ, KK, N) = TTL(IJ, N) * W2(IJ)
                     END IF
  410             CONTINUE
               END DO
            ELSE
               DO N = 1, NTDIM  
                  TTL(IJ, N) = T(IJ, K, N) * DZSIG(IJ, K)
                  LUP(IJ) = K
               END DO
            END IF
  200    CONTINUE
 1000 CONTINUE

      DO 2010 K = KSTR, KEND
         DO 2000 IJ = IJTSTR, IJTEND
            TL = T(IJ, K, 1)
            SL = T(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            R(IJ, K) = P1 / P2 - 1.D3
 2000   CONTINUE
 2010 CONTINUE

      CALL CHEKIN(  CONV, 'CONV',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')

      RETURN

c =====================================================================

      ENTRY DDENST(
     O                  R,
     I                  T)

      DO 3010 K = KSTR, KEND
         DO 3000 IJ = IJTSTR, IJTEND
            TL = T(IJ, K, 1)
            SL = T(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            R(IJ, K) = P1 / P2 - 1.D3
 3000   CONTINUE
 3010 CONTINUE

      RETURN

      END
