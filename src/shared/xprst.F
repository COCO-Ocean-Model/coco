      SUBROUTINE SECOEF(
     O                      C0,     C1,     C2,     C3,
     O                      C4,     C5,     C6,
     O                      D0,     D1,     D2,     D3,     D4,
     O                      D5,     D6,     D7,     D8,     D9)
      use zocdim
      use zocgrd
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Coefficients for the approximated equation of state by McDougall et
c al. (2003, JAOT).
c
c  HISTORY
c     '07.05.01  H.Hasumi
c
c ---------------------------------------------------------------------

      REAL*8 P10, P1T, P1TT, P1TTT, P1S, P1ST, P1SS,
     &       P1P, P1PTT,P1PS, P1PP, P1PPTT
      PARAMETER(P10    =  9.99843699D+2,
     &          P1T    =  7.35212840D+0,
     &          P1TT   = -5.45928211D-2,
     &          P1TTT  =  3.98476704D-4,
     &          P1S    =  2.96938239D+0,
     &          P1ST   = -7.23268813D-3,
     &          P1SS   =  2.12382341D-3,
     &          P1P    =  1.04004591D-2,
     &          P1PTT  =  1.03970529D-7,
     &          P1PS   =  5.18761880D-6,
     &          P1PP   = -3.24041825D-8,
     &          P1PPTT = -1.23869360D-11)
      REAL*8 P20, P2T, P2TT, P2TTT, P2TTTT, P2S, P2ST, P2STTT,
     &       P2SS, P2SSTT, P2P, P2PPTTT, P2PPPT
      PARAMETER(P20     =  1.D0,
     &          P2T     =  7.28606739D-3,
     &          P2TT    = -4.60835542D-5,
     &          P2TTT   =  3.68390573D-7,
     &          P2TTTT  =  1.80809186D-10,
     &          P2S     =  2.14691708D-3,
     &          P2ST    = -9.27062484D-6,
     &          P2STTT  = -1.78343643D-10,
     &          P2SS    =  4.76534122D-6,
     &          P2SSTT  =  1.63410736D-9,
     &          P2P     =  5.30848875D-6,
     &          P2PPTTT = -3.03175128D-16,
     &          P2PPPT  = -1.27934137D-17)

      REAL*8  C0(NZ), C1(NZ), C2(NZ), C3(NZ), C4(NZ), C5(NZ), C6(NZ)
      REAL*8  D0(NZ), D1(NZ), D2(NZ), D3(NZ), D4(NZ)
      REAL*8  D5(NZ), D6(NZ), D7(NZ), D8(NZ), D9(NZ)

      REAL*8  Z(NZ)
      INTEGER K

      Z(1) = 0.5D0 * DZ0(KSTR) / 1.0D2
      DO K = 2, NZ
         Z(K) = Z(K-1)
     &        + 0.5D0 * (DZ0(K+KSTR-1) + DZ0(K+KSTR-2)) / 1.D2
      END DO

      DO K = 1, NZ
         C0(K) = P10 + (P1P + P1PP * Z(K)) * Z(K)
         C1(K) = P1T
         C2(K) = P1TT + (P1PTT + P1PPTT * Z(K)) * Z(K)
         C3(K) = P1TTT
         C4(K) = P1S + P1PS * Z(K)
         C5(K) = P1ST
         C6(K) = P1SS

         D0(K) = P20 + P2P * Z(K)
         D1(K) = P2T + P2PPPT * Z(K) * Z(K) * Z(K)
         D2(K) = P2TT
         D3(K) = P2TTT + P2PPTTT * Z(K) * Z(K)
         D4(K) = P2TTTT
         D5(K) = P2S
         D6(K) = P2ST
         D7(K) = P2STTT
         D8(K) = P2SS
         D9(K) = P2SSTT
      END DO

      RETURN
      END

c **********************************************************************

      SUBROUTINE SECFBB(
     O                      C0,     C1,     C2,     C3,
     O                      C4,     C5,     C6,
     O                      D0,     D1,     D2,     D3,     D4,
     O                      D5,     D6,     D7,     D8,     D9)
      use zocdim
      use zocgrd
      IMPLICIT NONE

      REAL*8 P10, P1T, P1TT, P1TTT, P1S, P1ST, P1SS,
     &       P1P, P1PTT,P1PS, P1PP, P1PPTT
      PARAMETER(P10    =  9.99843699D+2,
     &          P1T    =  7.35212840D+0,
     &          P1TT   = -5.45928211D-2,
     &          P1TTT  =  3.98476704D-4,
     &          P1S    =  2.96938239D+0,
     &          P1ST   = -7.23268813D-3,
     &          P1SS   =  2.12382341D-3,
     &          P1P    =  1.04004591D-2,
     &          P1PTT  =  1.03970529D-7,
     &          P1PS   =  5.18761880D-6,
     &          P1PP   = -3.24041825D-8,
     &          P1PPTT = -1.23869360D-11)
      REAL*8 P20, P2T, P2TT, P2TTT, P2TTTT, P2S, P2ST, P2STTT,
     &       P2SS, P2SSTT, P2P, P2PPTTT, P2PPPT
      PARAMETER(P20     =  1.D0,
     &          P2T     =  7.28606739D-3,
     &          P2TT    = -4.60835542D-5,
     &          P2TTT   =  3.68390573D-7,
     &          P2TTTT  =  1.80809186D-10,
     &          P2S     =  2.14691708D-3,
     &          P2ST    = -9.27062484D-6,
     &          P2STTT  = -1.78343643D-10,
     &          P2SS    =  4.76534122D-6,
     &          P2SSTT  =  1.63410736D-9,
     &          P2P     =  5.30848875D-6,
     &          P2PPTTT = -3.03175128D-16,
     &          P2PPPT  = -1.27934137D-17)

      REAL*8  C0(NXYDIM), C1(NXYDIM), C2(NXYDIM), C3(NXYDIM)
      REAL*8  C4(NXYDIM), C5(NXYDIM), C6(NXYDIM)
      REAL*8  D0(NXYDIM), D1(NXYDIM), D2(NXYDIM), D3(NXYDIM)
      REAL*8  D4(NXYDIM), D5(NXYDIM), D6(NXYDIM), D7(NXYDIM)
      REAL*8  D8(NXYDIM), D9(NXYDIM)

      REAL*8  Z
      INTEGER K, IJ

      DO IJ = IJVSTR, IJVEND
         Z = 0.25D0 * 1.D-2 *
     &       (DEPT(IJ) + DEPT(IJ+LE) + DEPT(IJ+LN) + DEPT(IJ+LNE))

         C0(IJ) = P10 + (P1P + P1PP * Z) * Z
         C1(IJ) = P1T
         C2(IJ) = P1TT + (P1PTT + P1PPTT * Z) * Z
         C3(IJ) = P1TTT
         C4(IJ) = P1S + P1PS * Z
         C5(IJ) = P1ST
         C6(IJ) = P1SS

         D0(IJ) = P20 + P2P * Z
         D1(IJ) = P2T + P2PPPT * Z * Z * Z
         D2(IJ) = P2TT
         D3(IJ) = P2TTT + P2PPTTT * Z * Z
         D4(IJ) = P2TTTT
         D5(IJ) = P2S
         D6(IJ) = P2ST
         D7(IJ) = P2STTT
         D8(IJ) = P2SS
         D9(IJ) = P2SSTT
      END DO

      RETURN
      END

