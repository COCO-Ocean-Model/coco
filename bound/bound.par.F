      PROGRAM BOUND

c --- information -----------------------------------------------------
c
c  Generate simple boundary condition file
c
c  HISTORY
c     '97.03.24  H.Hasumi
c     '01.05.10  H.Hasumi
c
c ---------------------------------------------------------------------

#include "zocdim.F"

      INTEGER NFTSFC, NFTDMP
      INTEGER NFSSFC, NFSDMP
      INTEGER NFTAUX, NFTAUY
      PARAMETER(NFTSFC = 31, NFTDMP = 32)
      PARAMETER(NFSSFC = 33, NFSDMP = 34)
      PARAMETER(NFTAUX = 35, NFTAUY = 36)

      INTEGER IDATE(6), NTSTRT
      REAL*8  DIM2(NXG, NYG)
      REAL*8  TSFC(NYG), TDMP(NYG), SSFC(NYG), SDMP(NYG)
      REAL*8  TAUX(NYG), TAUY(NYG)
      CHARACTER CDATE*16, CHEAD(64)*16

      INTEGER I, J

      DATA NTSTRT / 0 /
      DATA IDATE / 0, 1, 1, 0, 0, 0 /
      DATA CHEAD / 64*'                ' /

      WRITE(CDATE, '(I6.6,5I2.2)') IDATE
      CHEAD(50) = CDATE
      WRITE(CDATE, '(I4.4,2I2.2,A1,3I2.2,A1)')
     &     IDATE(1), IDATE(2), IDATE(3), ' ',
     &     IDATE(4), IDATE(5), IDATE(6), ' '
      CHEAD(27) = CDATE
      CHEAD(38) = 'REAL8'
      WRITE(CHEAD(64), '(I16)') NXYG

      DO 100 J = 1, NYG
         TSFC(J) = 2.5D1 *
     &             (1.0D0 - (DBLE(J-1) - 0.5D0) / DBLE(NYG-2))
         TDMP(J) = 1.D0 / 3.D1 / 8.64D4
         SSFC(J) = 3.5D1
         SDMP(J) = 0.D0
         TAUX(J) = 0.D0
         TAUY(J) = 0.D0
  100 CONTINUE

      DO 200 J = 1, NYG
         DO 200 I = 1, NXG
            DIM2(I, J) = TSFC(J)
  200 CONTINUE
      CHEAD(3) = 'TSFC'
      WRITE(NFTSFC) CHEAD
      WRITE(NFTSFC) DIM2

      DO 300 J = 1, NYG
         DO 300 I = 1, NXG
            DIM2(I, J) = TDMP(J)
  300 CONTINUE
      CHEAD(3) = 'TDMP'
      WRITE(NFTDMP) CHEAD
      WRITE(NFTDMP) DIM2

      DO 400 J = 1, NYG
         DO 400 I = 1, NXG
            DIM2(I, J) = SSFC(J)
  400 CONTINUE
      CHEAD(3) = 'SSFC'
      WRITE(NFSSFC) CHEAD
      WRITE(NFSSFC) DIM2

      DO 500 J = 1, NYG
         DO 500 I = 1, NXG
            DIM2(I, J) = SDMP(J)
  500 CONTINUE
      CHEAD(3) = 'SDMP'
      WRITE(NFSDMP) CHEAD
      WRITE(NFSDMP) DIM2

      DO 600 J = 1, NYG
         DO 600 I = 1, NXG
            DIM2(I, J) = TAUX(J)
  600 CONTINUE
      CHEAD(3) = 'TAUX'
      WRITE(NFTAUX) CHEAD
      WRITE(NFTAUX) DIM2

      DO 700 J = 1, NYG
         DO 700 I = 1, NXG
            DIM2(I, J) = TAUY(J)
  700 CONTINUE
      CHEAD(3) = 'TAUY'
      WRITE(NFTAUY) CHEAD
      WRITE(NFTAUY) DIM2

      END
