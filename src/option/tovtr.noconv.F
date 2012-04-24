      SUBROUTINE OVTSET(
     O                       R,
     I                       T)
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
c     '10.04.14  Y.Komuro: dummy OVTURN (just diagnosing R)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8       R(NXYDIM, NZDIM)
      REAL*8       T(NXYDIM, NZDIM, NTDIM)
      REAL*8       H(NXYDIM)

      REAL*8     TTL(NXYDIM, NTDIM),     W2(NXYDIM)
      REAL*8    CONV(NXYDIM, NZDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM)
      INTEGER    LUP(NXYDIM)
      COMMON /WORK/ TTL, W2, CONV, DZSIG, LUP

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
