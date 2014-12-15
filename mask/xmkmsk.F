      PROGRAM MKMASK
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Generating masking arrays from the geometry file
c
c  HISTORY
c     '03.04.24  H.Hasumi: from COCO3.4
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '07.12.20  H.Hasumi: output file format is modified
c     '12.04.18  Y.Komuro: bug fix (argument inconsistency)
c
c ---------------------------------------------------------------------

#include "dimocn.F"

      INTEGER   ISTR, JSTR, KSTR
      PARAMETER(ISTR = 3, JSTR = 3, KSTR = 2)
      INTEGER   IEND, JEND, KEND
      PARAMETER(IEND = ISTR+NX-1, JEND = JSTR+NY-1, KEND = KSTR+NZ-1)
      INTEGER   NXDIM, NYDIM
      PARAMETER(NXDIM = NX+2*(ISTR-1), NYDIM = NY+2*(JSTR-1))
      INTEGER   NZDIM
      PARAMETER(NZDIM = NZ+KSTR)
      INTEGER   NXYDIM, NXYZDM
      PARAMETER(NXYDIM = NXDIM*NYDIM, NXYZDM = NXYDIM*NZDIM)
      INTEGER   ISSTR, JSSTR
      PARAMETER(ISSTR = 3, JSSTR = 3)

      REAL*8     ALD(NXYZDM),    BLD(NXYZDM)
      REAL*8      DY(NYDIM) ,     DZ(NZDIM)
      REAL*8     HIC(NIC+1)
      REAL*8     DZZ(NXYZDM),    DZV(NXYZDM)
      REAL*8    DGEO(NXYDIM)
      INTEGER   NBOT(NXYDIM),  NBOTV(NXYDIM)
      INTEGER   IBBL,  IMERC, ITHIRD
      REAL*8      DX
      REAL*8   ALONW,  ALONE,  ALATS,  ALATN
      REAL*8    PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT

      INTEGER INODES, JNODES

      REAL*8  MSK1D1(NZDIM ), MSK1D2(NZDIM )
      REAL*8  MSK1D3(NZDIM ), MSK1D4(NZDIM )
      REAL*8  MSK2D1(NXYDIM), MSK2D2(NXYDIM), MSK2D3(NXYDIM)
      REAL*8  MSK2D4(NXYDIM), MSK2D5(NXYDIM), MSK2D6(NXYDIM)
      REAL*8  MSK2D7(NXYDIM), MSK2D8(NXYDIM), MSK2D9(NXYDIM)
      REAL*8  MSK2DA(NXYDIM), MSK2DB(NXYDIM), MSK2DC(NXYDIM)
      REAL*8  MSK2DD(NXYDIM)
      REAL*8  MSK3D1(NXYZDM), MSK3D2(NXYZDM), MSK3D3(NXYZDM)
      REAL*8  MSK3D4(NXYZDM)
      INTEGER IMS2D1(NXYDIM)

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      PI = ATAN(1.D0) * 4.D0
      RAD = PI / 1.8D2
      RADIUS = 6370.D5
      OMEGA = 2.D0 * PI / 8.64D4
      CUNITI = (0.D0, 1.D0)

      CALL GETGEO(
     O               ALD,    BLD,   DGEO,   NBOT,
     O             IMERC, ITHIRD,   IBBL,
     O                DX,     DY,     DZ,    HIC,
     O             ALONW,  ALONE,  ALATS,  ALATN,
     O              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     O            INODES, JNODES,
     D                NX,     NY,     NZ,     KZ,    NIC,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W            MSK2D1, MSK2D2, IMS2D1)

      CALL GRDPRM(
     O               DZZ,    DZV,
     I               ALD,    BLD,   DGEO,   NBOT,  NBOTV,
     I                DX,     DY,     DZ,    HIC,
     I             IMERC, ITHIRD,   IBBL,
     I             ALONW,  ALONE,  ALATS,  ALATN,
     I              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     I            INODES, JNODES,
     D                NX,     NY,     NZ,     KZ,    NIC,
     D             NXDIM,  NYDIM,  NZDIM,  NTDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     C             ISSTR,  JSSTR,
     W            MSK1D1, MSK1D2, MSK1D3, MSK1D4,
     W            MSK2D1, MSK2D2, MSK2D3, MSK2D4, MSK2D5,
     W            MSK2D6, MSK2D7, MSK2D8, MSK2D9,
     W            MSK2DA, MSK2DB, MSK2DC, MSK2DD,
     W            MSK3D1)

      CALL TMASK (
     I               ALD,    BLD,    DZZ,   NBOT,   IBBL,
     D             NXDIM,  NYDIM,  NZDIM,     KZ,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W            MSK3D1, MSK3D2, MSK3D3, MSK3D4,
     W            MSK2D1, MSK2D2, MSK2D3)

      CALL VMASK (
     I               ALD,    BLD,    DZV,   NBOT,  NBOTV,   IBBL,
     D             NXDIM,  NYDIM,  NZDIM,     KZ,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W            MSK3D1, MSK3D2, MSK3D3, MSK3D4,
     W            MSK2D1, MSK2D2, MSK2D3)

      CALL BMASK (
     I               ALD,    BLD,   IBBL,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W            MSK3D1)

      WRITE(*, *)
      WRITE(*, '(2X,A)')
     & 'logical, save :: oinit, ofinal'

      WRITE(*, *)
      WRITE(*, '(A)') 'end module zocdim'

      END

c =====================================================================

      SUBROUTINE BMASK (
     I               ALD,    BLD,   IBBL,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W             AMSKB)

      INTEGER   JFILE
      PARAMETER(JFILE = 22)

      INTEGER  NXDIM,  NYDIM,  NZDIM
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND
      INTEGER   IBBL
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8   AMSKB(NXDIM, NYDIM, NZDIM)

      INTEGER      I,      J,      K
      REAL*8      AL,     AD

      DO 100 K = 1, NZDIM
         DO 100 J = 1, NYDIM
            DO 100 I = 1, NXDIM
               AMSKB(I, J, K) = 0.D0
  100 CONTINUE
      DO 110 K = KSTR, KEND
         DO 110 J = JSTR, JEND
            DO 110 I = ISTR, IEND
               AL = ALD(I, J  , K) * ALD(I+1, J  , K) *
     &              ALD(I, J+1, K) * ALD(I+1, J+1, K)
               AD = (ALD(I, J, K+1) + BLD(I, J, K+1)) *
     &              (ALD(I+1, J, K+1) + BLD(I+1, J, K+1)) *
     &              (ALD(I, J+1, K+1) + BLD(I, J+1, K+1)) *
     &              (ALD(I+1, J+1, K+1) + BLD(I+1, J+1, K+1))
               AMSKB(I, J, K) = AL * (1.D0 - AD)
  110 CONTINUE

      WRITE(JFILE) AMSKB

      RETURN
      END

c =====================================================================

      SUBROUTINE VMASK(
     I               ALD,    BLD,    DZV,   NBOT,  NBOTV,   IBBL,
     D             NXDIM,  NYDIM,  NZDIM,     KZ,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W             AMSKV,  AMFVX,  AMFVY,  AMFVZ,
     W            AMSKVB, AMSKV0, AMSKV1)

      INTEGER   JFILE
      PARAMETER(JFILE = 22)

      INTEGER  NXDIM,  NYDIM,  NZDIM,     KZ
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND
      INTEGER   IBBL
      INTEGER   NBOT(NXDIM, NYDIM),  NBOTV(NXDIM, NYDIM)
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8     DZV(NXDIM, NYDIM, NZDIM)
      REAL*8   AMSKV(NXDIM, NYDIM, NZDIM)
      REAL*8   AMFVX(NXDIM, NYDIM, NZDIM),  AMFVY(NXDIM, NYDIM, NZDIM)
      REAL*8   AMFVZ(NXDIM, NYDIM, NZDIM)
      REAL*8  AMSKVB(NXDIM, NYDIM)
      REAL*8  AMSKV0(NXDIM, NYDIM), AMSKV1(NXDIM, NYDIM)

      INTEGER      I,      J,      K
      INTEGER     KE,     KN,    KNE
      INTEGER IJVSTR, IJVEND
      REAL*8      AU,     AL

      DO 100 K = 1, NZDIM
         DO 100 J = 1, NYDIM
            DO 100 I = 1, NXDIM
               AMSKV(I, J, K) = 0.D0
  100 CONTINUE
      DO 110 K = KSTR, KEND
         DO 110 J = JSTR-1, JEND+1
            DO 110 I = ISTR-1, IEND+1
               AMSKV(I, J, K) = ALD(I, J, K) * ALD(I+1, J, K) *
     &                          ALD(I, J+1, K) * ALD(I+1, J+1, K)
  110 CONTINUE

      DO 200 K = 1, NZDIM
         DO 200 J = 1, NYDIM
            DO 200 I = 1, NXDIM
               AMFVX(I, J, K) = 0.D0
  200 CONTINUE
      DO 210 K = KSTR, KSTR+KZ-1
         DO 210 J = JSTR, JEND
            DO 210 I = ISTR, IEND+1
               AMFVX(I, J, K) = ALD(I, J, K) * ALD(I, J+1, K)
  210 CONTINUE
      DO 220 K = KSTR+KZ, KEND
         DO 220 J = JSTR, JEND
            DO 220 I = ISTR, IEND+1
               AMFVX(I, J, K) = ALD(I, J, K) * ALD(I, J+1, K) *
     &                          MIN(DZV(I, J, K), DZV(I-1, J, K))
  220 CONTINUE

      DO 300 K = 1, NZDIM
         DO 300 J = 1, NYDIM
            DO 300 I = 1, NXDIM
               AMFVY(I, J, K) = 0.D0
  300 CONTINUE
      DO 310 K = KSTR, KSTR+KZ-1
         DO 310 J = JSTR, JEND+1
            DO 310 I = ISTR, IEND
               AMFVY(I, J, K) = ALD(I, J, K) * ALD(I+1, J, K)
  310 CONTINUE
      DO 320 K = KSTR+KZ, KEND
         DO 320 J = JSTR, JEND+1
            DO 320 I = ISTR, IEND
               AMFVY(I, J, K) = ALD(I, J, K) * ALD(I+1, J, K) *
     &                          MIN(DZV(I, J, K), DZV(I, J-1, K))
  320 CONTINUE

      DO 400 K = 1, NZDIM
         DO 400 J = 1, NYDIM
            DO 400 I = 1, NXDIM
               AMFVZ(I, J, K) = 0.D0
  400 CONTINUE
      DO 410 K = KSTR, KEND
         DO 410 J = JSTR, JEND
            DO 410 I = ISTR, IEND
               AU = ALD(I, J, K-1) * ALD(I+1, J, K-1) *
     &              ALD(I, J+1, K-1) * ALD(I+1, J+1, K-1)
               AL = (1.D0 - ALD(I, J, K)) *
     &              (1.D0 - ALD(I+1, J, K)) *
     &              (1.D0 - ALD(I, J+1, K)) *
     &              (1.D0 - ALD(I+1, J+1, K))
               AMFVZ(I, J, K) = AU * (1.D0 - AL)
  410 CONTINUE

      DO 500 J = JSTR, JEND
         DO 500 I = ISTR, IEND
            IJVSTR = I + (J - 1) * NXDIM
            IF (AMSKV(I, J, KSTR) .EQ. 1.D0) GO TO 599
  500 CONTINUE
  599 CONTINUE

      DO 600 J = JEND, JSTR, -1
         DO 600 I = IEND, ISTR, -1
            IJVEND = I + (J - 1) * NXDIM
            IF (AMSKV(I, J, KSTR) .EQ. 1.D0) GO TO 699
  600 CONTINUE
  699 CONTINUE

      IF (IBBL .NE. 0) THEN
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               AMSKVB(I, J) = 0.D0
               AMSKV0(I, J) = 0.D0
               AMSKV1(I, J) = 0.D0
            END DO
         END DO

         DO J = JSTR-1, JEND+1
            DO I = ISTR-1, IEND+1
               K = NBOTV(I, J)
               AMSKV0(I, J) = AMSKV(I, J, K)
            END DO
         END DO

         DO J = JSTR-1, JEND+1
            DO I = ISTR-1, IEND+1
               K = NBOT(I, J)
               KE = NBOT(I+1, J)
               KN = NBOT(I, J+1)
               KNE = NBOT(I+1, J+1)
               AMSKVB(I, J) = BLD(I, J, K) * BLD(I+1, J, KE) *
     &                        BLD(I, J+1, KN) * BLD(I+1, J+1, KNE)
            END DO
         END DO

         DO J = JSTR-1, JEND+1
            DO I = ISTR-1, IEND+1
               K = NBOTV(I, J)
               AMSKV1(I, J) = AMSKV0(I, J) * (1.D0 - AMSKVB(I, J))
     &                      + AMSKVB(I, J)
            END DO
         END DO

         DO J = JSTR-1, JEND+1
            DO I = ISTR-1, IEND+1
               K = NBOTV(I, J)
               AMSKV(I, J, K) = AMSKV1(I, J)
            END DO
         END DO
         DO K = KSTR, KEND
            DO J = JSTR, JEND
               DO I = ISTR, IEND
                  AU = ALD(I, J, K-1) * ALD(I+1, J, K-1) *
     &                 ALD(I, J+1, K-1) * ALD(I+1, J+1, K-1)
                  AL = (1.D0 - ALD(I, J, K) - BLD(I, J, K)) *
     &                 (1.D0 - ALD(I+1, J, K) - BLD(I+1, J, K)) *
     &                 (1.D0 - ALD(I, J+1, K) - BLD(I, J+1, K)) *
     &                 (1.D0 - ALD(I+1, J+1, K) - BLD(I+1, J+1, K))
                  AMFVZ(I, J, K) = AU * (1.D0 - AL)
               END DO
            END DO
         END DO

         DO J = JSTR, JEND
            DO I = ISTR, IEND+1
               K = NBOT(I, J)
               KN = NBOT(I, J+1)
               AMFVX(I, J, KEND) = BLD(I, J, K) * BLD(I, J+1, KN) *
     &                             MIN(DZV(I, J, KEND),
     &                                 DZV(I-1, J, KEND))
            END DO
         END DO
         DO J = JSTR, JEND+1
            DO I = ISTR, IEND
               K = NBOT(I, J)
               KE = NBOT(I+1, J)
               AMFVY(I, J, KEND) = BLD(I, J, K) * BLD(I+1, J, KE) *
     &                             MIN(DZV(I, J, KEND),
     &                                 DZV(I, J-1, KEND))
            END DO
         END DO
      END IF

      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'ijvgst =', IJVSTR, ', ',
     & 'ijvged =', IJVEND

c      WRITE(JFILE) AMSKV, AMFVX, AMFVY, AMFVZ
      WRITE(JFILE) AMSKV
      WRITE(JFILE) AMFVX
      WRITE(JFILE) AMFVY
      WRITE(JFILE) AMFVZ
      IF (IBBL .NE. 0) THEN
c         WRITE(JFILE) AMSKVB, AMSKV0, AMSKV1
         WRITE(JFILE) AMSKVB
         WRITE(JFILE) AMSKV0
         WRITE(JFILE) AMSKV1
      END IF

      RETURN
      END

c =====================================================================

      SUBROUTINE TMASK(
     I               ALD,    BLD,    DZZ,   NBOT,   IBBL,
     D             NXDIM,  NYDIM,  NZDIM,     KZ,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W             AMSKT,  AMFTX,  AMFTY,  AMFTZ,
     W            AMSKTB, AMSKT0, AMSKT1)

      INTEGER   JFILE
      PARAMETER(JFILE = 22)

      INTEGER  NXDIM,  NYDIM,  NZDIM,     KZ
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND
      INTEGER   IBBL
      INTEGER   NBOT(NXDIM, NYDIM)
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8     DZZ(NXDIM, NYDIM, NZDIM)
      REAL*8   AMSKT(NXDIM, NYDIM, NZDIM)
      REAL*8   AMFTX(NXDIM, NYDIM, NZDIM),  AMFTY(NXDIM, NYDIM, NZDIM)
      REAL*8   AMFTZ(NXDIM, NYDIM, NZDIM)
      REAL*8  AMSKTB(NXDIM, NYDIM)
      REAL*8  AMSKT0(NXDIM, NYDIM), AMSKT1(NXDIM, NYDIM)

      INTEGER      I,      J,      K
      INTEGER IJTSTR, IJTEND

      DO 100 K = 1, NZDIM
         DO 100 J = 1, NYDIM
            DO 100 I = 1, NXDIM
               AMSKT(I, J, K) = 0.D0
  100 CONTINUE
      DO 110 K = KSTR, KEND
         DO 110 J = 1, NYDIM
            DO 110 I = 1, NXDIM
               AMSKT(I, J, K) = ALD(I, J, K)
  110 CONTINUE

      DO 200 K = 1, NZDIM
         DO 200 J = 1, NYDIM
            DO 200 I = 1, NXDIM
               AMFTX(I, J, K) = 0.D0
  200 CONTINUE
      DO 210 K = KSTR, KEND
         DO 210 J = JSTR-1, JEND+1
            DO 210 I = ISTR-1, IEND+2
               AMFTX(I, J, K) = ALD(I, J, K) * ALD(I-1, J, K) *
     &                          MIN(DZZ(I, J, K), DZZ(I-1, J, K))
  210 CONTINUE

      DO 300 K = 1, NZDIM
         DO 300 J = 1, NYDIM
            DO 300 I = 1, NXDIM
               AMFTY(I, J, K) = 0.D0
  300 CONTINUE
      DO 310 K = KSTR, KEND
         DO 310 J = JSTR-1, JEND+2
            DO 310 I = ISTR-1, IEND+1
               AMFTY(I, J, K) = ALD(I, J, K) * ALD(I, J-1, K) *
     &                          MIN(DZZ(I, J, K), DZZ(I, J-1, K))
  310 CONTINUE

      DO 400 K = 1, NZDIM
         DO 400 J = 1, NYDIM
            DO 400 I = 1, NXDIM
               AMFTZ(I, J, K) = 0.D0
  400 CONTINUE
      DO 410 K = KSTR, KEND
         DO 410 J = JSTR-1, JEND+1
            DO 410 I = ISTR-1, IEND+1
               AMFTZ(I, J, K) = ALD(I, J, K) * ALD(I, J, K-1)
  410 CONTINUE

      DO 500 J = JSTR, JEND
         DO 500 I = ISTR, IEND
            IJTSTR = I + (J - 1) * NXDIM
            IF (ALD(I, J, KSTR) .EQ. 1.D0) GO TO 599
  500 CONTINUE
  599 CONTINUE

      DO 600 J = JEND, JSTR, -1
         DO 600 I = IEND, ISTR, -1
            IJTEND = I + (J - 1) * NXDIM
            IF (ALD(I, J, KSTR) .EQ. 1.D0) GO TO 699
  600 CONTINUE
  699 CONTINUE

      IF (IBBL .NE. 0) THEN
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               K = NBOT(I, J)
               AMSKT0(I, J) = AMSKT(I, J, K)
               AMSKT1(I, J) = AMSKT(I, J, K) + BLD(I, J, K)
               AMSKTB(I, J) = BLD(I, J, K)
            END DO
         END DO

         DO J = 1, NYDIM
            DO I = 1, NXDIM
               K = NBOT(I, J)
               AMSKT(I, J, K) = AMSKT(I, J, K) + BLD(I, J, K)
            END DO
         END DO
         DO K = KSTR, KEND
            DO J = JSTR-1, JEND+1
               DO I = ISTR-1, IEND+1
                  AMFTZ(I, J, K) = AMSKT(I, J, K) * AMSKT(I, J, K-1)
               END DO
            END DO
         END DO
         DO J = JSTR-1, JEND+1
            DO I = ISTR-1, IEND+2
               AMFTX(I, J, KEND) = AMSKTB(I, J) * AMSKTB(I-1, J) *
     &                             MIN(DZZ(I, J, KEND),
     &                                 DZZ(I-1, J, KEND))
            END DO
         END DO
         DO J = JSTR-1, JEND+2
            DO I = ISTR-1, IEND+1
               AMFTY(I, J, KEND) = AMSKTB(I, J) * AMSKTB(I, J-1) *
     &                             MIN(DZZ(I, J, KEND),
     &                                 DZZ(I, J-1, KEND))
            END DO
         END DO
      END IF

      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'ijtgst =', IJTSTR, ', ',
     & 'ijtged =', IJTEND

c      WRITE(JFILE) AMSKT, AMFTX, AMFTY, AMFTZ
      WRITE(JFILE) AMSKT
      WRITE(JFILE) AMFTX
      WRITE(JFILE) AMFTY
      WRITE(JFILE) AMFTZ
      IF (IBBL .NE. 0) THEN
c         WRITE(JFILE) AMSKTB, AMSKT0, AMSKT1
         WRITE(JFILE) AMSKTB
         WRITE(JFILE) AMSKT0
         WRITE(JFILE) AMSKT1
      END IF

      RETURN
      END

c =====================================================================

      SUBROUTINE GRDPRM(
     O               DZZ,    DZV,
     I               ALD,    BLD,   DGEO,   NBOT,  NBOTV,
     I                DX,     DY,     DZ,    HIC,
     I             IMERC, ITHIRD,   IBBL,
     I             ALONW,  ALONE,  ALATS,  ALATN,
     I              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     I            INODES, JNODES,
     D                NX,     NY,     NZ,     KZ,    NIC,
     D             NXDIM,  NYDIM,  NZDIM,  NTDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     C             ISSTR,  JSSTR,
     W               DZM,  DEPTH,     DS,    DSM,
     W               COR,   DEPT,  RDEPV,    DYY,   DYYM,
     W               HXT,    HXU,    HYT,    HYU,
     W              HXYT,   HXYU,   HYXT,   HYXU,
     W              DZZM)
      IMPLICIT NONE

      INTEGER   JFILE
      PARAMETER(JFILE = 22)

      INTEGER     NX,     NY,     NZ,     KZ,    NIC
      INTEGER  NXDIM,  NYDIM,  NZDIM,  NTDIM
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND
      INTEGER  ISSTR,  JSSTR
      INTEGER  IMERC, ITHIRD,   IBBL
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8    DGEO(NXDIM, NYDIM)
      INTEGER   NBOT(NXDIM, NYDIM),  NBOTV(NXDIM, NYDIM)
      REAL*8      DX,     DY(NYDIM),     DZ(NZDIM)
      REAL*8     HIC(NIC+1)
      REAL*8   ALONW,  ALONE,  ALATS,  ALATN
      REAL*8    PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT
      INTEGER INODES, JNODES
      REAL*8     DZM(NZDIM),  DEPTH(NZDIM)
      REAL*8      DS(NZDIM),    DSM(NZDIM)
      REAL*8     DYY(NXDIM, NYDIM),   DYYM(NXDIM, NYDIM)
      REAL*8     DZZ(NXDIM, NYDIM, NZDIM)
      REAL*8     DZV(NXDIM, NYDIM, NZDIM)
      REAL*8    DZZM(NXDIM, NYDIM, NZDIM)

      REAL*8    ZBOT
      REAL*8     COR(NXDIM, NYDIM)
      REAL*8    DEPT(NXDIM, NYDIM),  RDEPV(NXDIM, NYDIM)
      REAL*8     HXT(NXDIM, NYDIM),    HXU(NXDIM, NYDIM)
      REAL*8     HYT(NXDIM, NYDIM),    HYU(NXDIM, NYDIM)
      REAL*8    HXYT(NXDIM, NYDIM),   HXYU(NXDIM, NYDIM)
      REAL*8    HYXT(NXDIM, NYDIM),   HYXU(NXDIM, NYDIM)

      REAL*8      CX,     CY,     CZ
      COMPLEX*16    ZPT,    ZPV,    ZQT,    ZQV
      COMPLEX*16    ZP1,    ZP2,    ZQ1,    ZQ2
      COMPLEX*16     Z1,     Z2,     Z3,     W3
      REAL*8    DELX,   DELY
      REAL*8   RLONT,  RLATT,  RLONV,  RLATV,  GLATV
      REAL*8   RLON1,  RLON2,  RLAT1,  RLAT2
      REAL*8     HX1,    HX2,    HY1,    HY2
      REAL*8      YS,     YN, DYMERC,     Y1,     Y2

      REAL*8     DEP
      REAL*8  FLGBBL
      INTEGER      I,      J,      K

      INTEGER    NXS,    NYS, NXSDIM, NYSDIM
      INTEGER  ISEND,  JSEND, IJSSTR, IJSEND


      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      REAL*8  ZTOLON, ZTOLAT
      REAL*8   FMTRC
      COMPLEX*16   CONF,  DCONF,   STOZ
      COMPLEX*16      Z,     ZP,     ZQ
      COMPLEX*16     ZA,     ZB,     ZC,     WC
      REAL*8    ALON,   ALAT
      STOZ(ALON, ALAT) = TAN(PI * 0.25D0 - ALAT * 0.5D0) *
     &                   EXP(ALON * CUNITI)
      ZTOLON(Z) = ATAN2(DIMAG(Z), DREAL(Z))
      ZTOLAT(Z) = PI * 0.5D0 - 2.D0 * ATAN(ABS(Z))
      CONF(Z, ZA, ZB, ZC, WC) = 
     &    (- ZB * (ZC - ZA) * Z / WC + ZA * (ZC - ZB))
     &    / (- (ZC - ZA) * Z / WC + ZC - ZB)
      DCONF(Z, ZA, ZB, ZC, WC) =
     &    (ZC - ZA) * (ZC - ZB) * (ZA - ZB) / WC
     &    / (ZC - ZB - (ZC - ZA) * Z / WC)
     &    / (ZC - ZB - (ZC - ZA) * Z / WC)
c      FMTRC(ZP, ZQ, ZA, ZB, ZC, WC) =
c     &     (1.D0 + ABS(ZP) * ABS(ZP))
c     &     / (1.D0 + ABS(ZQ) * ABS(ZQ))
c     &     / ABS(DCONF(ZQ, ZA, ZB, ZC, WC))
      FMTRC(ZP, ZQ, ZA, ZB, ZC, WC) =
     &     (1.D0 + ABS(ZQ) * ABS(ZQ)) *
     &     ABS(DCONF(ZQ, ZA, ZB, ZC, WC))
     &     / (1.D0 + ABS(ZP) * ABS(ZP))

      PI = 4.D0 * ATAN(1.D0)
      RAD = PI / 1.8D2
      RADIUS = 6370.D5
      OMEGA = 2.D0 * PI / 8.64D4

      ALONW = ALONW * RAD
      ALONE = ALONE * RAD
      ALATS = ALATS * RAD
      ALATN = ALATN * RAD
      PLON = PLON * RAD
      PLAT = PLAT * RAD
      QLON = QLON * RAD
      QLAT = QLAT * RAD
      ELON = ELON * RAD
      ELAT = ELAT * RAD

      DX = (ALONE - ALONW) / DBLE(NX)

      IF (IMERC .EQ. 2) THEN
         DO J = JSTR, JEND
            DY(J) = DY(J) * RAD
         END DO
      ELSE IF (IMERC .EQ. 1) THEN
         YS = LOG(TAN(ALATS * 0.5D0 + PI * 0.25D0))
         YN = LOG(TAN(ALATN * 0.5D0 + PI * 0.25D0))
         DYMERC = (YN - YS) / DBLE(NY)
         DO J = JSTR, JEND
            Y1 = YS + DYMERC * DBLE(J - JSTR)
            Y2 = YS + DYMERC * DBLE(J - JSTR + 1)
            DY(J) = 2.D0 * (ATAN(EXP(Y2)) - ATAN(EXP(Y1)))
         END DO
      ELSE
         DO J = JSTR, JEND
            DY(J) = (ALATN - ALATS) / DBLE(NY)
         END DO
      END IF
      DO J = 1, JSTR-1
         DY(J) = DY(JSTR)
      END DO
      DO J = JEND+1, NYDIM
         DY(J) = DY(JEND)
      END DO

      DO 100 J = 1, NYDIM
         DO 100 I = 1, NXDIM
            DYY(I, J) = DY(J)
  100 CONTINUE
      DO 110 J = 1, NYDIM-1
         DO 110 I = 1, NXDIM
            DYYM(I, J) = 0.5D0 * (DYY(I, J) + DYY(I, J+1))
  110 CONTINUE
      DO 120 I = 1, NXDIM
         DYYM(I, NYDIM) = 0.5D0 * DYY(I, NYDIM)
  120 CONTINUE

      DO 200 K = 1, NZDIM
         DZ(K) = DZ(K) * 1.D2
  200 CONTINUE
      DO 210 K = 1, KSTR-1
         DZM(K) = 0.5D0 * DZ(K)
         DEPTH(K) = 0.D0
  210 CONTINUE
      DO 220 K = KSTR, NZDIM
         DZM(K) = 0.5D0 * (DZ(K-1) + DZ(K))
         DEPTH(K) = DEPTH(K-1) + DZ(K)
  220 CONTINUE

      DO K = 1, NZDIM
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               DZZ (I, J, K) = DZ (K)
               DZZM(I, J, K) = DZM(K)
            END DO
         END DO
      END DO
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            K = NBOT(I, J) + KSTR - 1
            IF (K .GE. KSTR) THEN
               DZZ (I, J, K) = DZ(K) * DGEO(I, J)
               DZZM(I, J, K) = 0.5D0 * (DZ(K-1) + DZ(K) * DGEO(I, J))
            END IF
         END DO
      END DO
      IF (IBBL .NE. 0) THEN
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               DO K = NBOT(I, J)+KSTR, KEND-1
                  DZZ(I, J, K) = DZ(KEND)
               END DO
               K = NBOT(I, J)+KSTR
               DZZM(I, J, K) = 0.5D0 * (DZ(K-1) * DGEO(I, J) + DZ(KEND))
               DZZM(I, J, KEND) = DZZM(I, J, K)
            END DO
         END DO
      END IF

      DO K = 1, NZDIM
         DO J = JSTR, JEND
            DO I = ISTR, IEND
               DZV(I, J, K) = MIN(DZZ(I  , J  , K),
     &                            DZZ(I+1, J  , K),
     &                            DZZ(I  , J+1, K),
     &                            DZZ(I+1, J+1, K))
               IF (IBBL .NE. 0) THEN
                  FLGBBL = (1.D0 - BLD(I  , J  , K)) *
     &                     (1.D0 - BLD(I+1, J  , K)) *
     &                     (1.D0 - BLD(I  , J+1, K)) *
     &                     (1.D0 - BLD(I+1, J+1, K))
                  IF (FLGBBL .EQ. 0.D0) THEN
                     DZV(I, J, K) = MIN(DZZ(I  , J  , KEND),
     &                                  DZZ(I+1, J  , KEND),
     &                                  DZZ(I  , J+1, KEND),
     &                                  DZZ(I+1, J+1, KEND))
                  END IF
               END IF
            END DO
            DO I = 1, ISTR-1
               DZV(I, J, K) = DZV(IEND-ISTR+1+I, J, K)
            END DO
            DO I = IEND+1, NXDIM
               DZV(I, J, K) = DZV(ISTR-1+I-IEND, J, K)
            END DO
         END DO
         DO J = 1, JSTR-1
            DO I = 1, NXDIM
               DZV(I, J, K) = DZV(I, JSTR, K)
            END DO
         END DO
         DO J = JEND+1, NYDIM
            DO I = 1, NXDIM
               DZV(I, J, K) = DZV(I, JEND, K)
            END DO
         END DO
      END DO

      ZBOT = DEPTH(KSTR+KZ-1)
      DO K = KSTR, KSTR+KZ-1
         DS(K) = DZ(K) / ZBOT
      END DO
      DO K = 1, KSTR-1
         DS(K) = DS(KSTR)
      END DO
      DO K = KSTR+KZ, NZDIM
         DS(K) = DS(KSTR+KZ-1)
      END DO
      DO K = KSTR, KSTR+KZ-1
         DSM(K) = DZM(K) / ZBOT
      END DO
      DO K = 1, KSTR-1
         DSM(K) = DSM(KSTR)
      END DO
      DO K = KSTR+KZ, NZDIM
         DSM(K) = DS(KSTR+KZ-1)
      END DO

      WRITE(*, '(A)') 'module zocdim'
      WRITE(*, *)

      WRITE(*, '(2X,A)') 'use zocnod'
      WRITE(*, *)

      WRITE(*, '(2X,A)') 'implicit none'
      WRITE(*, *)

      WRITE(*, '(2X,A,2(A8,I7,A),A8,I7)')
     & 'integer, parameter :: ',
     & 'nxg    =', NX, ', ',
     & 'nyg    =', NY, ', ',
     & 'nz     =', NZ
      WRITE(*, '(2X,A,I7)')
     & 'integer, parameter :: kz     =', KZ
      WRITE(*, '(2X,A,I7)')
     & 'integer, parameter :: nic    =', NIC
      WRITE(*, '(2X,A,2(A8,I7,A),A8,I7)')
     & 'integer, parameter :: ',
     & 'nxgdim =', NXDIM, ', ',
     & 'nygdim =', NYDIM, ', ',
     & 'nzdim  =', NZDIM
      WRITE(*, '(2X,A,2(A8,I7,A),A8,I7)')
     & 'integer, parameter :: ',
     & 'igstr  =', ISTR, ', ',
     & 'jgstr  =', JSTR, ', ',
     & 'kstr   =', KSTR
      WRITE(*, '(2X,A,2(A8,I7,A),A8,I7)')
     & 'integer, parameter :: ',
     & 'igend  =', IEND, ', ',
     & 'jgend  =', JEND, ', ',
     & 'kend   =', KEND
      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'ijgstr =', ISTR+(JSTR-1)*NXDIM, ', ',
     & 'ijgend =', IEND+(JEND-1)*NXDIM
      WRITE(*, *)

      WRITE(*, '(2X,A)')
     & 'integer, parameter :: nxyg = nxg*nyg, nxyzg = nxyg*nz'
      WRITE(*, '(2X,A,A)')
     & 'integer, parameter :: nxygdm = nxgdim*nygdim, ',
     & 'nxyzgd = nxygdm*nzdim'
      WRITE(*, *)

      NXS = NX / INODES
      NYS = NY / JNODES
      NXSDIM = NXS + (ISSTR - 1) * 2
      NYSDIM = NYS + (JSSTR - 1) * 2
      ISEND = ISSTR + NXS - 1
      JSEND = JSSTR + NYS - 1
      IJSSTR = NXSDIM * (JSSTR - 1) + ISSTR
      IJSEND = NXSDIM * (JSEND - 1) + ISEND

      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'nx     =', NXS, ', ',
     & 'ny     =', NYS
      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'nxdim  =', NXSDIM, ', ',
     & 'nydim  =', NYSDIM
      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'istr   =', ISSTR, ', ',
     & 'jstr   =', JSSTR
      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'iend   =', ISEND, ', ',
     & 'jend   =', JSEND
      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'ijstr  =', IJSSTR, ', ',
     & 'ijend  =', IJSEND
      WRITE(*, *)

      WRITE(*, '(2X,A,I7)')
     & 'integer, parameter :: ntdim =', NTDIM
      WRITE(*, *)

      WRITE(*, '(2X,A)')
     & 'integer, parameter :: nxy = nx*ny, nxyz = nxy*nz'
      WRITE(*, '(2X,A,A)')
     & 'integer, parameter :: nxydim = nxdim*nydim, ',
     & 'nxyzdm = nxydim*nzdim'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: nztdim = nzdim*ntdim'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: nxyidm = nxydim*(nic+1)'
      WRITE(*, *)

      WRITE(*, '(2X,A)')
     & 'integer, parameter :: nwrk = nxyzdm*(ntdim*3+2)'
      WRITE(*, *)

      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lng = nxgdim, lsg = -nxgdim'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lneg = nxgdim+1, lnwg = nxgdim-1'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lseg = -nxgdim+1, lswg = -nxgdim-1'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lssg = -2 * nxgdim'
      WRITE(*, *)

      WRITE(*, '(2X,A)')
     & 'integer, parameter :: le = 1, lw = -1, ln = nxdim, ls = -nxdim'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lne = nxdim+1, lnw = nxdim-1'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lse = -nxdim+1, lsw = -nxdim-1'
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: lww = -2, lss = -2 * nxdim'
      WRITE(*, *)

      WRITE(*, '(2X,A,A8,I7,A,A8,I7)')
     & 'integer, parameter :: ',
     & 'inodes =', INODES, ', ',
     & 'jnodes =', JNODES
      WRITE(*, '(2X,A)')
     & 'integer, parameter :: icomm = istr-1, jcomm = jstr-1'
      WRITE(*, *)

      IF (ITHIRD .EQ. 0) THEN
         W3 = (-1.D0, 0.D0)
         CX = COS(PLON) * COS(PLAT) + COS(QLON) * COS(QLAT)
         CY = SIN(PLON) * COS(PLAT) + SIN(QLON) * COS(QLAT)
         CZ = SIN(PLAT) + SIN(QLAT)
         ELON = ATAN2(CY, CX)
         CX = SQRT(CX * CX + CY * CY)
         ELAT = ATAN2(CZ, CX)
      ELSE
         W3 = EXP(PLON * CUNITI)
      END IF
      Z1 = TAN(PI * 0.25D0 - PLAT * 0.5D0) * EXP(PLON * CUNITI)
      Z2 = TAN(PI * 0.25D0 - QLAT * 0.5D0) * EXP(QLON * CUNITI)
      Z3 = TAN(PI * 0.25D0 - ELAT * 0.5D0) * EXP(ELON * CUNITI)

      RLATV = ALATS
      DELX = DX * 0.1D0
      DO J = JSTR, JEND
         RLATT = RLATV + DY(J) * 0.5D0
         RLATV = RLATV + DY(J)
         DELY = DY(J) * 0.1D0
ccc         DO I = ISTR, IEND
         DO I = 1, NXDIM
            RLONT = ALONW + DX * (DBLE(I - ISTR) + 0.5D0)
            RLONV = ALONW + DX * DBLE(I - ISTR + 1)
            ZQT = STOZ(RLONT, RLATT)
            ZPT = CONF(ZQT, Z1, Z2, Z3, W3)
            ZQV = STOZ(RLONV, RLATV)
            ZPV = CONF(ZQV, Z1, Z2, Z3, W3)
            GLATV = ZTOLAT(ZPV)

            COR(I, J) = 2.D0 * OMEGA * SIN(GLATV)
            HYT(I, J) = RADIUS * FMTRC(ZPT, ZQT, Z1, Z2, Z3, W3)
            HXT(I, J) = HYT(I, J) * COS(RLATT)
            HYU(I, J) = RADIUS * FMTRC(ZPV, ZQV, Z1, Z2, Z3, W3)
            HXU(I, J) = HYU(I, J) * COS(RLATV)

            RLON1 = RLONT + DELX * 0.5D0
            RLON2 = RLONT - DELX * 0.5D0
            ZQ1 = STOZ(RLON1, RLATT)
            ZQ2 = STOZ(RLON2, RLATT)
            ZP1 = CONF(ZQ1, Z1, Z2, Z3, W3)
            ZP2 = CONF(ZQ2, Z1, Z2, Z3, W3)
            HY1 = FMTRC(ZP1, ZQ1, Z1, Z2, Z3, W3)
            HY2 = FMTRC(ZP2, ZQ2, Z1, Z2, Z3, W3)
            HYXT(I, J) = RADIUS * (HY1 - HY2) / DELX
     &                   / HXT(I, J) / HYT(I, J)
            
            
            RLON1 = RLONV + DELX * 0.5D0
            RLON2 = RLONV - DELX * 0.5D0
            ZQ1 = STOZ(RLON1, RLATV)
            ZQ2 = STOZ(RLON2, RLATV)
            ZP1 = CONF(ZQ1, Z1, Z2, Z3, W3)
            ZP2 = CONF(ZQ2, Z1, Z2, Z3, W3)
            HY1 = FMTRC(ZP1, ZQ1, Z1, Z2, Z3, W3)
            HY2 = FMTRC(ZP2, ZQ2, Z1, Z2, Z3, W3)
            HYXU(I, J) = RADIUS * (HY1 - HY2) / DELX
     &                   / HXU(I, J) / HYU(I, J)

            RLAT1 = RLATT + DELY * 0.5D0
            RLAT2 = RLATT - DELY * 0.5D0
            ZQ1 = STOZ(RLONT, RLAT1)
            ZQ2 = STOZ(RLONT, RLAT2)
            ZP1 = CONF(ZQ1, Z1, Z2, Z3, W3)
            ZP2 = CONF(ZQ2, Z1, Z2, Z3, W3)
            HX1 = FMTRC(ZP1, ZQ1, Z1, Z2, Z3, W3) * COS(RLAT1)
            HX2 = FMTRC(ZP2, ZQ2, Z1, Z2, Z3, W3) * COS(RLAT2)
            HXYT(I, J) = RADIUS * (HX1 - HX2) / DELY
     &                   / HXT(I, J) / HYT(I, J)

            RLAT1 = RLATV + DELY * 0.5D0
            RLAT2 = RLATV - DELY * 0.5D0
            ZQ1 = STOZ(RLONV, RLAT1)
            ZQ2 = STOZ(RLONV, RLAT2)
            ZP1 = CONF(ZQ1, Z1, Z2, Z3, W3)
            ZP2 = CONF(ZQ2, Z1, Z2, Z3, W3)
            HX1 = FMTRC(ZP1, ZQ1, Z1, Z2, Z3, W3) * COS(RLAT1)
            HX2 = FMTRC(ZP2, ZQ2, Z1, Z2, Z3, W3) * COS(RLAT2)
            HXYU(I, J) = RADIUS * (HX1 - HX2) / DELY
     &                   / HXU(I, J) / HYU(I, J)
            
         END DO
      END DO

      DO J = 1, JSTR-1
         DO I = 1, NXDIM
            COR(I, J) = COR(I, JSTR)
            HXT(I, J) = HXT(I, JSTR)
            HXU(I, J) = HXU(I, JSTR)
            HYT(I, J) = HYT(I, JSTR)
            HYU(I, J) = HYU(I, JSTR)
            HXYT(I, J) = HXYT(I, JSTR)
            HXYU(I, J) = HXYU(I, JSTR)
            HYXT(I, J) = HYXT(I, JSTR)
            HYXU(I, J) = HYXU(I, JSTR)
         END DO
      END DO
      DO J = JEND+1, NYDIM
         DO I = 1, NXDIM
            COR(I, J) = COR(I, JEND)
            HXT(I, J) = HXT(I, JEND)
            HXU(I, J) = HXU(I, JEND)
            HYT(I, J) = HYT(I, JEND)
            HYU(I, J) = HYU(I, JEND)
            HXYT(I, J) = HXYT(I, JEND)
            HXYU(I, J) = HXYU(I, JEND)
            HYXT(I, J) = HYXT(I, JEND)
            HYXU(I, J) = HYXU(I, JEND)
         END DO
      END DO

      DO 400 J = 1, NYDIM
         DO 400 I = 1, NXDIM
            RDEPV(I, J) = -999.D0
  400 CONTINUE
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            DEP = 0.D0
            DO K = KSTR, NBOT(I, J)+KSTR-1
               DEP = DEP + DZZ(I, J, K)
            END DO
            IF (IBBL .NE. 0) THEN
               K = NBOT(I, J) + KSTR
               DEP = DEP + DZZ(I, J, K) * BLD(I, J, K)
            END IF
            DEPT(I, J) = DEP
         END DO
      END DO

      DO J = 1, NYDIM
         DO I = 1, NXDIM
            NBOTV(I, J) = 0
         END DO
      END DO
      DO J = JSTR-1, JEND+1
         DO I = ISTR-1, IEND+1
            IF (      (NBOT(I  , J  ) .GT. 0)
     &          .AND. (NBOT(I+1, J  ) .GT. 0)
     &          .AND. (NBOT(I  , J+1) .GT. 0)
     &          .AND. (NBOT(I+1, J+1) .GT. 0)) THEN
               NBOTV(I, J) = MIN(NBOT(I, J), NBOT(I+1, J),
     &                           NBOT(I, J+1), NBOT(I+1, J+1))
               DEP = 0.D0
               DO K = KSTR, NBOTV(I, J)+KSTR-1
                  DEP = DEP + DZV(I, J, K)
               END DO
               IF (IBBL .NE. 0) THEN
                  FLGBBL = BLD(I, J, NBOT(I, J)+KSTR) *
     &                     BLD(I+1, J, NBOT(I+1, J)+KSTR) *
     &                     BLD(I, J+1, NBOT(I, J+1)+KSTR) *
     &                     BLD(I+1, J+1, NBOT(I+1, J+1)+KSTR)
                  IF (FLGBBL .GT. 0.D0) THEN
                     DEP = DEP + DZV(I, J, KEND)
                     NBOTV(I, J) = NBOTV(I, J) + 1
                  END IF
               END IF
               RDEPV(I, J) = 1.D0 / DEP
            END IF
         END DO
      END DO

c      WRITE(JFILE) DX, DYY, DYYM, DZZ, DZZM, DZV, DZ, DS, DSM, HIC
      WRITE(JFILE) DX
      WRITE(JFILE) DYY
      WRITE(JFILE) DYYM
      WRITE(JFILE) DZZ
      WRITE(JFILE) DZZM
      WRITE(JFILE) DZV
      WRITE(JFILE) DZ
      WRITE(JFILE) DS
      WRITE(JFILE) DSM
      WRITE(JFILE) HIC
c      WRITE(JFILE) RADIUS, ZBOT
      WRITE(JFILE) RADIUS
      WRITE(JFILE) ZBOT
c      WRITE(JFILE) COR, DEPT, RDEPV
      WRITE(JFILE) COR
      WRITE(JFILE) DEPT
      WRITE(JFILE) RDEPV
c      WRITE(JFILE) HXT, HXU, HYT, HYU, HXYT, HXYU, HYXT, HYXU
      WRITE(JFILE) HXT
      WRITE(JFILE) HXU
      WRITE(JFILE) HYT
      WRITE(JFILE) HYU
      WRITE(JFILE) HXYT
      WRITE(JFILE) HXYU
      WRITE(JFILE) HYXT
      WRITE(JFILE) HYXU

      IF (IBBL .NE. 0) THEN
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               K = NBOT(I, J) + KSTR
               NBOT(I, J) = NBOT(I, J) + NINT(BLD(I, J, K))
            END DO
         END DO
      END IF
      DO 500 J = 1, NYDIM
         DO 500 I = 1, NXDIM
            NBOT (I, J) = NBOT (I, J) + KSTR - 1
            NBOTV(I, J) = NBOTV(I, J) + KSTR - 1
  500 CONTINUE

      IF (IBBL .EQ. 0) THEN
         WRITE(JFILE) NBOT
      ELSE
c         WRITE(JFILE) NBOT, NBOTV
         WRITE(JFILE) NBOT
         WRITE(JFILE) NBOTV
      END IF

      RETURN
      END

c =====================================================================

      SUBROUTINE GETGEO(
     O               ALD,    BLD,   DGEO,   NBOT,
     O             IMERC, ITHIRD,   IBBL,
     O                DX,     DY,     DZ,    HIC,
     O             ALONW,  ALONE,  ALATS,  ALATN,
     O              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     O            INODES, JNODES,
     D                NX,     NY,     NZ,     KZ,    NIC,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     W              AGEO,   BGEO,   IGEO)

      INTEGER   IGFILE
      PARAMETER(IGFILE = 5)
      INTEGER   NROUND
      PARAMETER(NROUND = 1000)

      INTEGER     NX,     NY,     NZ,     KZ,    NIC
      INTEGER  NXDIM,  NYDIM,  NZDIM
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8    AGEO(NXDIM, NYDIM)
      REAL*8    BGEO(NXDIM, NYDIM)
      REAL*8    DGEO(NXDIM, NYDIM)
      INTEGER   IGEO(NXDIM, NYDIM),   NBOT(NXDIM, NYDIM)
      REAL*8      DX,     DY(NYDIM),     DZ(NZDIM)
      REAL*8     HIC(NIC+1)
      REAL*8   ALONW,  ALONE,  ALATS,  ALATN
      REAL*8    PLON,   PLAT,   QLON,   QLAT
      REAL*8    ELON,   ELAT
      INTEGER INODES, JNODES
      INTEGER  IMERC, ITHIRD,   IBBL
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND

      INTEGER     LX,     LY,     LZ,     JZ,    LIC
      INTEGER      I,      J,      K

      IBBL = 0

      READ(IGFILE, *) LX
      READ(IGFILE, *) LY
      READ(IGFILE, *) LZ
      READ(IGFILE, *) JZ
      READ(IGFILE, *) LIC
      IF (     (LX .NE. NX) .OR. (LY .NE. NY)
     &    .OR. (LZ .NE. NZ) .OR. (JZ .NE. KZ)
     &    .OR. (LIC .NE. NIC)) THEN
         WRITE(*, *) ' #### BAD GEOMETRY FILE #### '
         WRITE(*, *) ' #### DIMENSION ERROR   #### '
         STOP
      END IF
      READ(IGFILE, *) INODES, JNODES

      IF (     (MOD(NX, INODES) .NE. 0)
     &    .OR. (MOD(NY, JNODES) .NE. 0)) THEN
         WRITE(*, *) ' #### NODE NUMBER ERROR #### '
         WRITE(*, *) ' The size of sub-regions must be identical. '
         STOP
      END IF

      READ(IGFILE, *) ALONW, ALONE
      READ(IGFILE, *) ALATS, ALATN
      READ(IGFILE, *) IMERC
      IF (IMERC .EQ. 2) THEN
         READ(IGFILE, *) (DY(J), J = JSTR, JEND)
      END IF
      READ(IGFILE, *) (DZ(K), K = KSTR, KEND)
      READ(IGFILE, *) (HIC(K), K = 1, NIC)
      HIC(NIC+1) = 1.D+30
      READ(IGFILE, *) PLON, PLAT
      READ(IGFILE, *) QLON, QLAT
      READ(IGFILE, *) ITHIRD
      IF (ITHIRD .EQ. 1) THEN
         READ(IGFILE, *) ELON, ELAT
      END IF

      DO J = JEND, JSTR, -1
         READ(IGFILE, *) (AGEO(I, J), I = ISTR, IEND)
      END DO
      READ(IGFILE, *, END=999) IBBL
      IF (IBBL .NE. 0) THEN
         DO J = JEND, JSTR, -1
            READ(IGFILE, *) (BGEO(I, J), I = ISTR, IEND)
         END DO
      END IF
 999  CONTINUE
      DO J = JSTR, JEND
         DO I = ISTR, IEND
            IGEO(I, J) = NINT(AGEO(I, J) * DBLE(NROUND)) / NROUND
            DGEO(I, J) = DBLE(  NINT(AGEO(I, J) * DBLE(NROUND))
     &                        - IGEO(I, J) * NROUND)
     &                   / DBLE(NROUND)
         END DO
      END DO
      DO J = JSTR, JEND
         DO I = ISTR, IEND
            IF ((IGEO(I, J) .GT. 0) .AND. (IGEO(I, J) .LT. KZ)) THEN
               WRITE(*, *) ' #### BAD GEOMETRY FILE         #### '
               WRITE(*, *) ' #### BOTTOM LEVEL INCONSISTENT #### '
               STOP
            END IF
         END DO
      END DO

      DO 110 K = 1, KSTR-1
         DZ(K) = DZ(KSTR)
  110 CONTINUE
      DO 120 K = KEND+1, NZDIM
         DZ(K) = DZ(KEND)
  120 CONTINUE

      DO 200 J = JSTR, JEND
         DO 210 I = 1, ISTR-1
            IGEO(I, J) = IGEO(IEND-ISTR+1+I, J)
            DGEO(I, J) = DGEO(IEND-ISTR+1+I, J)
            BGEO(I, J) = BGEO(IEND-ISTR+1+I, J)
  210    CONTINUE
         DO 220 I = IEND+1, NXDIM
            IGEO(I, J) = IGEO(ISTR-1+I-IEND, J)
            DGEO(I, J) = DGEO(ISTR-1+I-IEND, J)
            BGEO(I, J) = BGEO(ISTR-1+I-IEND, J)
  220    CONTINUE
  200 CONTINUE
      DO 230 J = 1, JSTR-1
         DO 240 I = 1, NXDIM
            IGEO(I, J) = 0
            DGEO(I, J) = 0.D0
            BGEO(I, J) = 0.D0
  240    CONTINUE
  230 CONTINUE
      DO 250 J = JEND+1, NYDIM
         DO 260 I = 1, NXDIM
            IGEO(I, J) = 0
            DGEO(I, J) = 0.D0
            BGEO(I, J) = 0.D0
  260    CONTINUE
  250 CONTINUE

      DO 300 K = KSTR, KEND
         DO 300 J = 1, NYDIM
            DO 300 I = 1, NXDIM
               IF (IGEO(I, J) .GE. K-KSTR+1) THEN
                  ALD(I, J, K) = 1.D0
               ELSE IF (      (IGEO(I, J) .EQ. K-KSTR)
     &                  .AND. (DGEO(I, J) .GT. 0.D0)) THEN
                  ALD(I, J, K) = 1.D0
               ELSE
                  ALD(I, J, K) = 0.D0
               END IF
  300 CONTINUE
      DO 310 K = 1, KSTR-1
         DO 310 J = 1, NYDIM
            DO 310 I = 1, NXDIM
               ALD(I, J, K) = ALD(I, J, KSTR)
  310 CONTINUE
      DO 320 K = KEND+1, NZDIM
         DO 320 J = 1, NYDIM
            DO 320 I = 1, NXDIM
               ALD(I, J, K) = 0.D0
  320 CONTINUE

      DO 400 J = 1, NYDIM
         DO 400 I = 1, NXDIM
            IF (IGEO(I, J) .GE. 0) THEN
               IF (DGEO(I, J) .EQ. 0.D0) THEN
                  NBOT(I, J) = IGEO(I, J)
               ELSE
                  NBOT(I, J) = IGEO(I, J) + 1
               END IF
            ELSE
               NBOT(I, J) = 0
            END IF
  400 CONTINUE

      DO J = 1, NYDIM
         DO I = 1, NXDIM
            IF (DGEO(I, J) .EQ. 0.D0) THEN
               DGEO(I, J) = 1.D0
            END IF
         END DO
      END DO

      IF (IBBL .NE. 0) THEN
         DO K = 1, NZDIM
            DO J = 1, NYDIM
               DO I = 1, NXDIM
                  BLD(I, J, K) = 0.D0
               END DO
            END DO
         END DO
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               K = NBOT(I, J) + KSTR
               IF (      (IGEO(I, J) .GT. 0)
     &             .AND. (BGEO(I, J) .GT. 0.D0)) THEN
                  BLD(I, J, K) = 1.D0
               END IF
            END DO
         END DO
      END IF

      RETURN
      END
