c --- information -----------------------------------------------------
c
c  Set boudary conditions
c
c  HISTORY
c     '02.10.10  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.03.25  H.Hasumi: STBBT2 restored (for SOM)
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
      SUBROUTINE STBCTR(
     M                       T,      R)
      use zocdim
      use zocmsk
      IMPLICIT NONE

#include "mpif.h"

      REAL*8       T(NXDIM, NYDIM, NZDIM, NTDIM)
      REAL*8       R(NXDIM, NYDIM, NZDIM)

      INTEGER      I,      J,      K,     IJ,      N 

      DO 230 N = 1, NTDIM
         DO 220 K = 1, KSTR-1
            DO 210 J = 1, NYDIM
               DO 200 I = 1, NXDIM
                  T(I, J, K, N) = T(I, J, KSTR, N)
  200          CONTINUE
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE

      DO 330 N = 1, NTDIM
         DO 320 J = 1, NYDIM
            DO 310 I = 1, NXDIM
               IJ = I + (J - 1) * NXDIM
               T(I, J, NBOT(IJ)+1, N) = T(I, J, NBOT(IJ), N)
  310       CONTINUE      
  320    CONTINUE
  330 CONTINUE

      RETURN
      END
#else
c *********************************************************************

      SUBROUTINE STBCSH(
     M                      SH)
      use zocdim
      IMPLICIT NONE

      REAL*8      SH(NXDIM, NYDIM)

      INTEGER      I,      J

      DO 100 J = JSTR, JEND
         DO 110 I = 1, ISTR-1
            SH(I, J) = SH(IEND-ISTR+1+I, J)
  110    CONTINUE
         DO 120 I = IEND+1, NXDIM
            SH(I, J) = SH(ISTR-1+I-IEND, J)
  120    CONTINUE
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCBT(
     M                     UBT,    VBT)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8     UBT(NXDIM, NYDIM),    VBT(NXDIM, NYDIM)

      INTEGER      I,      J

      DO 100 J = JSTR, JEND
         DO 110 I = 1, ISTR-1
            UBT(I, J) = UBT(IEND-ISTR+1+I, J)
            VBT(I, J) = VBT(IEND-ISTR+1+I, J)
  110    CONTINUE
         DO 120 I = IEND+1, NXDIM
            UBT(I, J) = UBT(ISTR-1+I-IEND, J)
            VBT(I, J) = VBT(ISTR-1+I-IEND, J)
  120    CONTINUE
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCTR(
     M                       T,      R)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8       T(NXDIM, NYDIM, NZDIM, NTDIM)
      REAL*8       R(NXDIM, NYDIM, NZDIM)

      INTEGER      I,      J,      K,     IJ,      N

      DO 130 K = KSTR, KEND
         DO 120 J = JSTR, JEND
            DO 100 I = 1, ISTR-1
               R(I, J, K) = R(IEND-ISTR+1+I, J, K)
  100       CONTINUE
            DO 110 I = IEND+1, NXDIM
               R(I, J, K) = R(ISTR-1+I-IEND, J, K)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
      
      DO 180 N = 1, NTDIM
         DO 170 K = KSTR, KEND
            DO 160 J = JSTR, JEND
               DO 140 I = 1, ISTR-1
                  T(I, J, K, N) = T(IEND-ISTR+1+I, J, K, N)
  140          CONTINUE
               DO 150 I = IEND+1, NXDIM
                  T(I, J, K, N) = T(ISTR-1+I-IEND, J, K, N)
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE

      DO 230 N = 1, NTDIM
         DO 220 K = 1, KSTR-1
            DO 210 J = 1, NYDIM
               DO 200 I = 1, NXDIM
                  T(I, J, K, N) = T(I, J, KSTR, N)
  200          CONTINUE
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE

      DO 330 N = 1, NTDIM
         DO 320 J = 1, NYDIM
            DO 310 I = 1, NXDIM
               IJ = I + (J - 1) * NXDIM
               T(I, J, NBOT(IJ)+1, N) = T(I, J, NBOT(IJ), N)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCUV(
     M                       U,      V)
      use zocdim
      IMPLICIT NONE

      REAL*8       U(NXDIM, NYDIM, NZDIM)
      REAL*8       V(NXDIM, NYDIM, NZDIM)

      INTEGER      I,      J,      K

      DO 100 K = KSTR, KEND
         DO 110 J = JSTR, JEND
            DO 120 I = 1, ISTR-1
               U(I, J, K) = U(IEND-ISTR+1+I, J, K)
               V(I, J, K) = V(IEND-ISTR+1+I, J, K)
  120       CONTINUE
            DO 130 I = IEND+1, NXDIM
               U(I, J, K) = U(ISTR-1+I-IEND, J, K)
               V(I, J, K) = V(ISTR-1+I-IEND, J, K)
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCGG(
     M                     GXX,    GYY)
      use zocdim
      IMPLICIT NONE

      REAL*8     GXX(NXDIM, NYDIM),    GYY(NXDIM, NYDIM)

      INTEGER      I,      J

      DO 100 J = JSTR, JEND
         DO 110 I = 1, ISTR-1
            GXX(I, J) = GXX(IEND-ISTR+1+I, J)
            GYY(I, J) = GYY(IEND-ISTR+1+I, J)
  110    CONTINUE
         DO 120 I = IEND+1, NXDIM
            GXX(I, J) = GXX(ISTR-1+I-IEND, J)
            GYY(I, J) = GYY(ISTR-1+I-IEND, J)
  120    CONTINUE
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCAD(
     M                    UADV,   VADV,   WADV)
      use zocdim
      IMPLICIT NONE

      REAL*8    UADV(NXDIM, NYDIM, NZDIM)
      REAL*8    VADV(NXDIM, NYDIM, NZDIM)
      REAL*8    WADV(NXDIM, NYDIM, NZDIM, 9)

      INTEGER      I,      J,      K,      L

      DO 130 K = KSTR, KEND
         DO 100 J = JSTR, JEND
            DO 110 I = 1, ISTR-1
               UADV(I, J, K) = UADV(IEND-ISTR+1+I, J, K)
               VADV(I, J, K) = VADV(IEND-ISTR+1+I, J, K)
  110       CONTINUE
            DO 120 I = IEND+1, NXDIM
               UADV(I, J, K) = UADV(ISTR-1+I-IEND, J, K)
               VADV(I, J, K) = VADV(ISTR-1+I-IEND, J, K)
  120       CONTINUE
  100    CONTINUE
  130 CONTINUE

      DO 240 L = 1, 9
         DO 230 K = KSTR, KEND
            DO 200 J = JSTR, JEND
               DO 210 I = I, ISTR-1
                  WADV(I, J, K, L) = WADV(IEND-ISTR+1+I, J, K, L)
  210          CONTINUE
               DO 220 I = IEND+1, NXDIM
                  WADV(I, J, K, L) = WADV(ISTR-1+I-IEND, J, K, L)
  220          CONTINUE
  200       CONTINUE
  230    CONTINUE
  240 CONTINUE


      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCAV(
     M                     AMV,    AHV)
      use zocdim
      IMPLICIT NONE

      REAL*8     AMV(NXDIM, NYDIM, NZDIM)
      REAL*8     AHV(NXDIM, NYDIM, NZDIM)

      INTEGER      I,      J,      K

      DO 100 K = KSTR, KEND
         DO 110 J = JSTR, JEND
            DO 120 I = 1, ISTR-1
               AMV(I, J, K) = AMV(IEND-ISTR+1+I, J, K)
               AHV(I, J, K) = AHV(IEND-ISTR+1+I, J, K)
  120       CONTINUE
            DO 130 I = IEND+1, NXDIM
               AMV(I, J, K) = AMV(ISTR-1+I-IEND, J, K)
               AHV(I, J, K) = AHV(ISTR-1+I-IEND, J, K)
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCTU(
     M                    TAUX,   TAUY)
      use zocdim
      IMPLICIT NONE

      REAL*8    TAUX(NXDIM, NYDIM),   TAUY(NXDIM, NYDIM)

      INTEGER      I,      J

      DO 100 J = JSTR, JEND
         DO 110 I = 1, ISTR-1
            TAUX(I, J) = TAUX(IEND-ISTR+1+I, J)
            TAUY(I, J) = TAUY(IEND-ISTR+1+I, J)
  110    CONTINUE
         DO 120 I = IEND+1, NXDIM
            TAUX(I, J) = TAUX(ISTR-1+I-IEND, J)
            TAUY(I, J) = TAUY(ISTR-1+I-IEND, J)
  120    CONTINUE
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCSI(
     M                       A,     HI,     EI,     TI,     HS)

      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8       A(NXDIM, NYDIM, 0:NIC)
      REAL*8      HI(NXDIM, NYDIM, 0:NIC)
      REAL*8      EI(NXDIM, NYDIM, 0:NIC)
      REAL*8      TI(NXDIM, NYDIM, 0:NIC)
      REAL*8      HS(NXDIM, NYDIM, 0:NIC)

      INTEGER      I,      J,      K

      DO K = 0, NIC
         DO J = JSTR, JEND
            DO I = 1, ISTR-1
               A (I, J, K) = A (IEND-ISTR+1+I, J, K)
               HI(I, J, K) = HI(IEND-ISTR+1+I, J, K)
               HS(I, J, K) = HS(IEND-ISTR+1+I, J, K)
               EI(I, J, K) = EI(IEND-ISTR+1+I, J, K)
               TI(I, J, K) = TI(IEND-ISTR+1+I, J, K)
            END DO
            DO I = IEND+1, NXDIM
               A (I, J, K) = A (ISTR-1+I-IEND, J, K)
               HI(I, J, K) = HI(ISTR-1+I-IEND, J, K)
               HS(I, J, K) = HS(ISTR-1+I-IEND, J, K)
               EI(I, J, K) = EI(ISTR-1+I-IEND, J, K)
               TI(I, J, K) = TI(ISTR-1+I-IEND, J, K)
            END DO
         END DO
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBCVI(
     M                      UI,     VI)

      use zocdim
      IMPLICIT NONE

      REAL*8      UI(NXDIM, NYDIM)
      REAL*8      VI(NXDIM, NYDIM)

      INTEGER      I,      J

      DO 110 J = JSTR, JEND
         DO 120 I = 1, ISTR-1
            UI(I, J) = UI(IEND-ISTR+1+I, J)
            VI(I, J) = VI(IEND-ISTR+1+I, J)
  120    CONTINUE
         DO 130 I = IEND+1, NXDIM
            UI(I, J) = UI(ISTR-1+I-IEND, J)
            VI(I, J) = VI(ISTR-1+I-IEND, J)
  130    CONTINUE
  110 CONTINUE

      RETURN
      END
#endif
#ifdef OPT_BBL
c *********************************************************************
c --- information -----------------------------------------------------
c
c  BBL variables are defined at two levels: one at KEND and the other
c at the NBOT (or NBOTV). This routine keeps two values consistent.
c
c  HISTORY
c     '01.02.09  H.Hasumi
c
c ---------------------------------------------------------------------

      SUBROUTINE STBBGV(
     M                      GX,     GY)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)

      INTEGER     IJ,      K

      DO IJ = IJVSTR, IJVEND
         K = NBOTV(IJ)
         GX(IJ, K) = GX(IJ, K) * (1.D0 - AMSKVB(IJ))
     &             + GX(IJ, KEND) * AMSKVB(IJ)
         GY(IJ, K) = GY(IJ, K) * (1.D0 - AMSKVB(IJ))
     &             + GY(IJ, KEND) * AMSKVB(IJ)
      END DO

      DO K = KSTR, KEND
         DO IJ = IJVSTR, IJVEND
            GX(IJ, K) = GX(IJ, K) * AMSKV(IJ, K)
            GY(IJ, K) = GY(IJ, K) * AMSKV(IJ, K)
         END DO
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBBUV(
     M                       U,      V)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8       U(NXYDIM, NZDIM),      V(NXYDIM, NZDIM)

      INTEGER     IJ,      K

      DO IJ = IJVSTR, IJVEND
         K = NBOTV(IJ)
         U(IJ, KEND) = U(IJ, K) * AMSKVB(IJ)
         V(IJ, KEND) = V(IJ, K) * AMSKVB(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBBVT(
     M                       U,      V)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8       U(NXYDIM, NZDIM),      V(NXYDIM, NZDIM)

      INTEGER     IJ,      K

      DO IJ = IJVSTR, IJVEND
         K = NBOTV(IJ)
         U(IJ, K) = U(IJ, KEND) * AMSKVB(IJ)
     &            + U(IJ, K) * (1.D0 - AMSKVB(IJ))
         V(IJ, K) = V(IJ, KEND) * AMSKVB(IJ)
     &            + V(IJ, K) * (1.D0 - AMSKVB(IJ))
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBBGT(
     M                     ADT,  DIFFZ)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)

      INTEGER     IJ,      K,      N

      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND
            K = NBOT(IJ)
            ADT(IJ, K, N) = ADT(IJ, K, N) * (1.D0 - AMSKTB(IJ))
     &                    + ADT(IJ, KEND, N) * AMSKTB(IJ)
         END DO

         DO K = KSTR, KEND
            DO IJ = IJTSTR, IJTEND
               ADT(IJ, K, N) = ADT(IJ, K, N) * AMSKT(IJ, K)
               DIFFZ(IJ, K) = DIFFZ(IJ, K) * AMSKT(IJ, K)
            END DO
         END DO
      ENDDO

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBBTR(
     M                       T)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8       T(NXYDIM, NZDIM, NTDIM)

      INTEGER     IJ,      K,      N

      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND
            K = NBOT(IJ)
            T(IJ, KEND, N) = T(IJ, K, N)
         END DO
      ENDDO

      RETURN
      END

c *********************************************************************

      SUBROUTINE STBBT2(
     M                       T)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8       T(NXYDIM, NZDIM, NTDIM)      

      INTEGER     IJ,      K,      N

      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND
            K = NBOT(IJ)
            T(IJ, K, N) = T(IJ, KEND, N) * AMSKTB(IJ)
     &                  + T(IJ, K, N) * (1.D0 - AMSKTB(IJ))
         END DO
      ENDDO

      RETURN
      END
#endif
