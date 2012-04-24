      SUBROUTINE VELVAD(
     O               UADV,   VADV,   WADV,
     I                  U,      V,      W)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Velocity for momentum advection
c
c  HISTORY
c     '99.08.17  H.Hasumi
c     '00.03.17  H.Hasumi: for Ishizaki scheme
c        In this version, UADV and VADV are not advection velocity
c        for V-box, but velocity at T-points. WADV is nine-
c        direction vertical velocity: the direction is designated
c        by the last argument (1 for purely vertical, 2 for north-
c        ward, and rounds clockwise to 9 for northwestward).
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '01.02.20  H.Hasumi: bug fix (weight for UADV, VADV)
c                          for incorporating BBL model
c     '01.05.02  H.Hasumi
c     '01.05.17  H.Hasumi: bug fix in VELVAB
c     '02.10.17  T.Suzuki: diagonal momentum advection from BBL
c     '02.11.??  T.Suzuki: remove cross diagonal momentum advection
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.02.23  Y.Komuro: bug fix by Dr. Kurogi
c                          (loop length extended)
c     '10.04.14  M.kurogi
c     '10.04.14  M.Kurogi: for tripolar grid
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM, 9)
      REAL*8       U(NXYDIM, NZDIM),      V(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM)

      REAL*8      RN(NXYDIM, NZDIM),    RNN(NXYDIM, NZDIM)
      SAVE RN, RNN
      REAL*8     PLO(NXYDIM, NZDIM),    RCR(NXYDIM, NZDIM)
      SAVE PLO,RCR

      INTEGER     IJ,      K,      L
      REAL*8      AN

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

         DO 21 K = 1, NZDIM
            DO 20 IJ = 1, NXYDIM
               RN (IJ, K) = 0.D0
               RNN(IJ, K) = 0.D0
               PLO(IJ, K) = 0.D0
 20         CONTINUE
 21      CONTINUE

         DO 31 K = KSTR, KEND
            DO 30 IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
               AN = AMSKV(IJ   , K) + AMSKV(IJ+LW , K)
     &            + AMSKV(IJ+LS, K) + AMSKV(IJ+LSW, K)
               IF (AN .EQ. 0.D0) THEN
                  RN(IJ, K) = 0.D0
               ELSE
                  RN(IJ, K) = 1.D0 / AN
               END IF
 30         CONTINUE
 31      CONTINUE

         DO 41 K = KSTR+1, KEND
            DO 40 IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
               IF (RN(IJ, K) .EQ. RN(IJ, K-1)) THEN
                  RNN(IJ, K) = 0.D0
               ELSE IF (RN(IJ, K-1) .EQ. 0.D0) THEN
                  RNN(IJ, K) = 0.D0
               ELSE
                  RNN(IJ, K) = RN(IJ, K) * RN(IJ, K-1)
               END IF
 40         CONTINUE
 41      CONTINUE

         DO K = KSTR+1, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
               IF (      (RN(IJ,K) .EQ. 0.5D0)
     &             .AND. (RN(IJ,K-1) .EQ. 0.25D0)) THEN
                  RNN(IJ, K) = 0.25D0
                  RCR(IJ, K) = 0.0D0
               ELSE
                  RCR(IJ, K) = 1.0D0
               END IF
            END DO
         END DO

         DO K = KSTR+1, KEND
            DO IJ = IJVSTR-NXDIM-1, IJVEND+NXDIM+1
               IF (      (AMSKV(IJ, K) .EQ. 0.D0)
     &             .AND. (AMSKV(IJ, K-1) .NE. 0.D0)) THEN
                  PLO(IJ, K) = 1.D0
               END IF
            END DO
         END DO

#ifdef OPT_BBL
C         CALL ADMKV1
         CALL ADMSKV
#endif
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            UADV(IJ, K) = 0.D0
            VADV(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE
      DO 170 L = 1, 9
         DO 160 K = 1, NZDIM
            DO 150 IJ = 1, NXYDIM
               WADV(IJ, K, L) = 0.D0
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE

      DO 210 K = KSTR, KSTR+KZ-1
         DO 200 IJ = IJTSTR, IJTEND
            UADV(IJ, K) = 0.25D0 * (  U(IJ+LSW, K) + U(IJ+LW , K)
     &                              + U(IJ+LS , K) + U(IJ    , K))
            VADV(IJ, K) = 0.25D0 * (  V(IJ+LSW, K) + V(IJ+LW , K)
     &                              + V(IJ+LS , K) + V(IJ    , K))
  200    CONTINUE
  210 CONTINUE
      DO 230 K = KSTR+KZ, KEND
         DO 220 IJ = IJTSTR, IJTEND
            UADV(IJ, K) = 0.25D0 * (  U(IJ+LSW, K) * DZV(IJ+LSW, K)
     &                              + U(IJ+LW , K) * DZV(IJ+LW , K)
     &                              + U(IJ+LS , K) * DZV(IJ+LS , K)
     &                              + U(IJ    , K) * DZV(IJ    , K)) *
     &                    AMSKT(IJ, K)
            VADV(IJ, K) = 0.25D0 * (  V(IJ+LSW, K) * DZV(IJ+LSW, K)
     &                              + V(IJ+LW , K) * DZV(IJ+LW , K)
     &                              + V(IJ+LS , K) * DZV(IJ+LS , K)
     &                              + V(IJ    , K) * DZV(IJ    , K)) *
     &                    AMSKT(IJ, K)
  220    CONTINUE
  230 CONTINUE

      DO 310 K = KSTR, KEND
         DO 300 IJ = IJVSTR, IJVEND
            WADV(IJ, K, 1) = AMSKV(IJ, K) *
     &                       (  W(IJ+LNE, K) * RN(IJ+LNE, K-1)
     &                        + W(IJ+LN , K) * RN(IJ+LN , K-1)
     &                        + W(IJ+LE , K) * RN(IJ+LE , K-1)
     &                        + W(IJ    , K) * RN(IJ    , K-1))
  300    CONTINUE
  310 CONTINUE

      DO 410 K = KSTR+1, KEND
         DO 400 IJ = IJVSTR, IJVEND
            WADV(IJ, K, 2) = AMSKV(IJ, K) * PLO(IJ+LN, K) *
     &                       (  W(IJ+LN , K) * RNN(IJ+LN , K)
     &                        + W(IJ+LNE, K) * RNN(IJ+LNE, K))
            WADV(IJ, K, 4) = AMSKV(IJ, K) * PLO(IJ+LE, K) *
     &                       (  W(IJ+LNE, K) * RNN(IJ+LNE, K)
     &                        + W(IJ+LE , K) * RNN(IJ+LE , K))
            WADV(IJ, K, 6) = AMSKV(IJ, K) * PLO(IJ+LS, K) *
     &                       (  W(IJ+LE , K) * RNN(IJ+LE , K)
     &                        + W(IJ    , K) * RNN(IJ    , K))
            WADV(IJ, K, 8) = AMSKV(IJ, K) * PLO(IJ+LW, K) *
     &                       (  W(IJ    , K) * RNN(IJ    , K)
     &                        + W(IJ+LN , K) * RNN(IJ+LN , K))
            WADV(IJ, K, 3) = AMSKV(IJ, K) * PLO(IJ+LNE, K) *
     &                       W(IJ+LNE, K) * RNN(IJ+LNE, K) *
     &                       RCR(IJ+LNE, K)
            WADV(IJ, K, 5) = AMSKV(IJ, K) * PLO(IJ+LSE, K) *
     &                       W(IJ+LE , K) * RNN(IJ+LE , K) *
     &                       RCR(IJ+LE, K)
            WADV(IJ, K, 7) = AMSKV(IJ, K) * PLO(IJ+LSW, K) *
     &                       W(IJ    , K) * RNN(IJ    , K) *
     &                       RCR(IJ, K)
            WADV(IJ, K, 9) = AMSKV(IJ, K) * PLO(IJ+LNW, K) *
     &                       W(IJ+LN , K) * RNN(IJ+LN , K) *
     &                       RCR(IJ+LN, K)
  400    CONTINUE
  410 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE VELVAB(
     M               WADV,
     I                  W)
      IMPLICIT NONE

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8    WADV(NXYDIM, NZDIM, 9)
      REAL*8       W(NXYDIM, NZDIM)

      REAL*8      AN
      REAL*8      RB(NXYDIM)
      REAL*8     RNB(NXYDIM, NZDIM)
      REAL*8    RNNB(NXYDIM, NZDIM)
      REAL*8    RCRB(NXYDIM)
      INTEGER     IJ,      K,    KUP
      SAVE  RB, RNB, RNNB

      LOGICAL OFIRST
      SAVE OFIRST

      DATA OFIRST / .TRUE. /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.

         DO IJ = 1, NXYDIM
            RB(IJ) = 0.D0
         END DO
         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
            AN = AMSKVB(IJ)    + AMSKVB(IJ+LW)
     &         + AMSKVB(IJ+LS) + AMSKVB(IJ+LSW)
            IF (AN .EQ. 0.D0) THEN
               RB(IJ) = 0.D0
            ELSE
               RB(IJ) = 1.D0 / AN
            END IF
         END DO

         DO K = 1, NZDIM
            DO IJ = 1, NXYDIM
               RNB (IJ, K) = 0.D0
               RNNB(IJ, K) = 0.D0
            END DO
         END DO
         DO IJ = 1, NXYDIM
            RCRB(IJ) = 1.D0
         END DO

         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1
            K = NBOT(IJ) - AMSKTB(IJ)
            RNB(IJ, K) = 1.D0
            AN = AMSKV(IJ,     K)  + AMSKV(IJ+LW,  K)
     &         + AMSKV(IJ+LS,  K)  + AMSKV(IJ+LSW, K)
            IF (AN .EQ. 0.D0) THEN
               RNNB(IJ, K) = 0.D0
            ELSE
               RNNB(IJ, K) = RB(IJ) / AN
            END IF
            IF (AN .EQ. 2.D0) THEN
               RCRB(IJ) = 0.D0
               RNNB(IJ, K) = RB(IJ) 
            ELSE
               RCRB(IJ) = 1.D0
            END IF
         END DO
      END IF

      DO IJ = IJVSTR, IJVEND
            KUP = NBOTV(IJ) - AMSKVB(IJ)
            WADV(IJ, KEND, 1) = AMSKVB(IJ) *
     &              (  W(IJ+LNE, KEND) * RB(IJ+LNE) * RNB(IJ+LNE, KUP)
     &               + W(IJ+LN , KEND) * RB(IJ+LN ) * RNB(IJ+LN,  KUP)
     &               + W(IJ+LE , KEND) * RB(IJ+LE ) * RNB(IJ+LE,  KUP)
     &               + W(IJ    , KEND) * RB(IJ    ) * RNB(IJ,     KUP))

            K = NBOTV(IJ+LN) - AMSKVB(IJ+LN)
      WADV(IJ, KEND, 2) = AMSKVB(IJ) * AMSKV(IJ+LN, K) *
     & (  W(IJ+LN , KEND) * RNNB(IJ+LN , K) * (1.D0 - RNB(IJ+LN,  KUP))
     &  + W(IJ+LNE, KEND) * RNNB(IJ+LNE, K) * (1.D0 - RNB(IJ+LNE, KUP)))
            K = NBOTV(IJ+LE) - AMSKVB(IJ+LE)
      WADV(IJ, KEND, 4) = AMSKVB(IJ) * AMSKV(IJ+LE, K) *
     & (  W(IJ+LNE, KEND) * RNNB(IJ+LNE, K) * (1.D0 - RNB(IJ+LNE, KUP))
     &  + W(IJ+LE , KEND) * RNNB(IJ+LE , K) * (1.D0 - RNB(IJ+LE,  KUP)))
            K = NBOTV(IJ+LS) - AMSKVB(IJ+LS)
      WADV(IJ, KEND, 6) = AMSKVB(IJ) * AMSKV(IJ+LS, K) *
     & (  W(IJ+LE , KEND) * RNNB(IJ+LE , K) * (1.D0 - RNB(IJ+LE,  KUP))
     &  + W(IJ    , KEND) * RNNB(IJ    , K) * (1.D0 - RNB(IJ   ,  KUP)))
            K = NBOTV(IJ+LW) - AMSKVB(IJ+LW)
      WADV(IJ, KEND, 8) = AMSKVB(IJ) * AMSKV(IJ+LW, K) *
     & (  W(IJ    , KEND) * RNNB(IJ    , K) * (1.D0 - RNB(IJ   ,  KUP))
     &  + W(IJ+LN , KEND) * RNNB(IJ+LN , K) * (1.D0 - RNB(IJ+LN,  KUP)))

            K = NBOTV(IJ+LNE) - AMSKVB(IJ+LNE)
      WADV(IJ, KEND, 3) = AMSKVB(IJ) * AMSKV(IJ+LNE, K)* RCRB(IJ+LNE) *
     &    W(IJ+LNE, KEND) * RNNB(IJ+LNE, K) * (1.D0 - RNB(IJ+LNE, KUP)) 
            K = NBOTV(IJ+LSE) - AMSKVB(IJ+LSE)
      WADV(IJ, KEND, 5) = AMSKVB(IJ) * AMSKV(IJ+LSE, K)* RCRB(IJ+LE ) *
     &    W(IJ+LE , KEND) * RNNB(IJ+LE , K) * (1.D0 - RNB(IJ+LE,  KUP))
            K = NBOTV(IJ+LSW) - AMSKVB(IJ+LSW)
      WADV(IJ, KEND, 7) = AMSKVB(IJ) * AMSKV(IJ+LSW, K)* RCRB(IJ    ) *
     &    W(IJ    , KEND) * RNNB(IJ    , K) * (1.D0 - RNB(IJ,     KUP))
            K = NBOTV(IJ+LNW) - AMSKVB(IJ+LNW)
      WADV(IJ, KEND, 9) = AMSKVB(IJ) * AMSKV(IJ+LNW, K)* RCRB(IJ+LN ) *
     &    W(IJ+LN , KEND) * RNNB(IJ+LN , K) * (1.D0 - RNB(IJ+LN,  KUP))
      END DO

      DO IJ = IJVSTR, IJVEND
         K = NBOTV(IJ)
         WADV(IJ, K, 1) = WADV(IJ, K, 1) * (1.D0 - AMSKVB(IJ))
     &                  + WADV(IJ, KEND, 1) * AMSKVB(IJ)
      END DO

      RETURN
      END
#endif


#ifdef OPT_TRIPOLE
      SUBROUTINE EXCNGW(WADV)
      IMPLICIT NONE
#include "zocdim.F"
      REAL*8    WADV(NXDIM, NYDIM, NZDIM, 9)
      REAL*8   WADV2(NXDIM, NYDIM, NZDIM, 9)

      INTEGER I, J, K, N, ISTV, JSTV

      IF(JRANK .NE. JNODES-1) RETURN

      DO N=1,9
      DO K=1,NZDIM
      DO J=1,NYDIM
      DO I=1,NXDIM
         WADV2(I,J,K,N)=WADV(I,J,K,N)
      END DO
      END DO
      END DO
      END DO

      IF(INODES .EQ. 1) THEN
        DO K=1, NZDIM
        DO J=JEND+1, NYDIM
        DO I=1,NXDIM
          WADV(I,J,K,2)=WADV2(I,J,K,6)
          WADV(I,J,K,3)=WADV2(I,J,K,7)
          WADV(I,J,K,4)=WADV2(I,J,K,8)
          WADV(I,J,K,5)=WADV2(I,J,K,9)
          WADV(I,J,K,6)=WADV2(I,J,K,2)
          WADV(I,J,K,7)=WADV2(I,J,K,3)
          WADV(I,J,K,8)=WADV2(I,J,K,4)
          WADV(I,J,K,9)=WADV2(I,J,K,5)
        END DO
        END DO
        END DO

        J=JEND
        DO K=1, NZDIM
        DO I=NXDIM/2+1,NXDIM
          WADV(I,J,K,2)=WADV2(I,J,K,6)
          WADV(I,J,K,3)=WADV2(I,J,K,7)
          WADV(I,J,K,4)=WADV2(I,J,K,8)
          WADV(I,J,K,5)=WADV2(I,J,K,9)
          WADV(I,J,K,6)=WADV2(I,J,K,2)
          WADV(I,J,K,7)=WADV2(I,J,K,3)
          WADV(I,J,K,8)=WADV2(I,J,K,4)
          WADV(I,J,K,9)=WADV2(I,J,K,5)
        END DO
        END DO

      ELSE IF(IRANK .LT. INODES/2) THEN
        DO K=1, NZDIM
        DO J=JEND+1, NYDIM
        DO I=1,NXDIM
          WADV(I,J,K,2)=WADV2(I,J,K,6)
          WADV(I,J,K,3)=WADV2(I,J,K,7)
          WADV(I,J,K,4)=WADV2(I,J,K,8)
          WADV(I,J,K,5)=WADV2(I,J,K,9)
          WADV(I,J,K,6)=WADV2(I,J,K,2)
          WADV(I,J,K,7)=WADV2(I,J,K,3)
          WADV(I,J,K,8)=WADV2(I,J,K,4)
          WADV(I,J,K,9)=WADV2(I,J,K,5)
        END DO
        END DO
        END DO
      ELSE
        DO K=1, NZDIM
        DO J=JEND, NYDIM
        DO I=1,NXDIM
          WADV(I,J,K,2)=WADV2(I,J,K,6)
          WADV(I,J,K,3)=WADV2(I,J,K,7)
          WADV(I,J,K,4)=WADV2(I,J,K,8)
          WADV(I,J,K,5)=WADV2(I,J,K,9)
          WADV(I,J,K,6)=WADV2(I,J,K,2)
          WADV(I,J,K,7)=WADV2(I,J,K,3)
          WADV(I,J,K,8)=WADV2(I,J,K,4)
          WADV(I,J,K,9)=WADV2(I,J,K,5)
        END DO
        END DO
        END DO
      END IF
    
      RETURN
      END
#endif
