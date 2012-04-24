      PROGRAM INITIAL

c --- information -----------------------------------------------------
c
c  Generate initial condition file (horizontally uniform, no currents)
c
c  HISTORY
c     '97.03.19  H.Hasumi
c     '01.05.10  H.Hasumi
c
c ---------------------------------------------------------------------

#include "zocdim.F"

      INTEGER NFINIT
      PARAMETER(NFINIT = 41)

      INTEGER IDATE(6)
      REAL*8  T(NX, NY, NZ), S(NX, NY, NZ)
      REAL*8  DIM3(NX, NY, NZ), DIM2(NX, NY)
      CHARACTER  CDATE*16, CHEAD(64)*16

      INTEGER I, J, K
      REAL*8  TINIT(NZ), SINIT(NZ)

      DATA CHEAD / 64*'                ' /
      DATA IDATE / 0, 1, 1, 0, 0, 0 /
      DATA TINIT / NZ * 5.D0 /
      DATA SINIT / NZ * 34.7D0 /

      DO 100 K = 1, NZ
         DO 100 J = 1, NY
            DO 100 I = 1, NX
               DIM3(I, J, K) = 0.0
               T(I, J, K) = TINIT(K)
               S(I, J, K) = SINIT(K)
  100 CONTINUE
      DO 200 J = 1, NY
         DO 200 I = 1, NX
            DIM2(I, J) = 0.0
  200 CONTINUE

      WRITE(CDATE, '(I6.6,5I2.2)') IDATE
      CHEAD(50) = CDATE
      WRITE(CDATE, '(I4.4,2I2.2,A1,3I2.2,A1)')
     &     IDATE(1), IDATE(2), IDATE(3), ' ',
     &     IDATE(4), IDATE(5), IDATE(6), ' '
      CHEAD(27) = CDATE
      CHEAD(38) = 'REAL8'

      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM3
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM3
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) T
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) S
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM2
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM2
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM2
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM3
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM3
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM3
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM2
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFINIT) CHEAD
      WRITE(NFINIT) DIM2

      END
