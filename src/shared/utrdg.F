      SUBROUTINE THMASC(
     M                  DR,     DI,
     C                   A,     BR,     BI,      C)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Solve the linear equations expressed by tri-diagonal matrix
c
c  HISTORY
c     '99.04.12  H.Hasumi: from CCSR2
c     '01.02.02  H.Hasumi: solve T and S at once
c     '02.05.29  H.Nakano: tracer dimension
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

      REAL*8       A(NXYDIM, NZDIM),      C(NXYDIM, NZDIM)
      REAL*8      BR(NXYDIM, NZDIM),     BI(NXYDIM, NZDIM)
      REAL*8      DR(NXYDIM, NZDIM),     DI(NXYDIM, NZDIM)

      REAL*8      R1,    PR1,    PI1,    QR1,    QI1
      REAL*8      RR,     RI,      R,     SR,     SI
      INTEGER     IJ,      K

      DO 100 IJ = IJVSTR, IJVEND
         R1 = BR(IJ, KSTR) * BR(IJ, KSTR)
     &      + BI(IJ, KSTR) * BI(IJ, KSTR)
         PR1 =   C(IJ, KSTR) * BR(IJ, KSTR) / R1
         PI1 = - C(IJ, KSTR) * BI(IJ, KSTR) / R1
         QR1 = (  DR(IJ, KSTR) * BR(IJ, KSTR)
     &          + DI(IJ, KSTR) * BI(IJ, KSTR)) / R1
         QI1 = (  DI(IJ, KSTR) * BR(IJ, KSTR)
     &          - DR(IJ, KSTR) * BI(IJ, KSTR)) / R1
         BR(IJ, KSTR) = PR1
         BI(IJ, KSTR) = PI1
         DR(IJ, KSTR) = QR1
         DI(IJ, KSTR) = QI1
  100 CONTINUE

      DO 210 K = KSTR+1, KEND
         DO 200 IJ = IJVSTR, IJVEND
            RR = BR(IJ, K) - A(IJ, K) * BR(IJ, K-1)
            RI = BI(IJ, K) - A(IJ, K) * BI(IJ, K-1)
            SR = DR(IJ, K) - A(IJ, K) * DR(IJ, K-1)
            SI = DI(IJ, K) - A(IJ, K) * DI(IJ, K-1)
            R = RR * RR + RI * RI
            BR(IJ, K) =   C(IJ, K) * RR / R
            BI(IJ, K) = - C(IJ, K) * RI / R
            DR(IJ, K) = (SR * RR + SI * RI) / R
            DI(IJ, K) = (SI * RR - SR * RI) / R
  200    CONTINUE
  210 CONTINUE

      DO 310 K = KEND-1, KSTR, -1
         DO 300 IJ = IJVSTR, IJVEND
            DR(IJ, K) = DR(IJ, K) - BR(IJ, K) * DR(IJ, K+1)
     &                            + BI(IJ, K) * DI(IJ, K+1)
            DI(IJ, K) = DI(IJ, K) - BI(IJ, K) * DR(IJ, K+1)
     &                            - BR(IJ, K) * DI(IJ, K+1)
  300    CONTINUE
  310 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE THOMAS(
     M                 ADT,     AC,
     I                  AA,     AB)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

      REAL*8      AA(NXYDIM, NZDIM),     AB(NXYDIM, NZDIM)
      REAL*8      AC(NXYDIM, NZDIM)
      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)   

      REAL*8      FC
      INTEGER     IJ,       K,       N

      DO 100 IJ = IJTSTR, IJTEND
         AC(IJ, KSTR) = AC(IJ, KSTR) / AB(IJ, KSTR)
  100 CONTINUE
      DO N = 1, NTDIM
         DO 110 IJ = IJTSTR, IJTEND
            ADT(IJ, KSTR, N) = ADT(IJ, KSTR, N) / AB(IJ, KSTR)
 110     CONTINUE
      END DO


      DO 210 K = KSTR+1, KEND
         DO 200 IJ = IJTSTR, IJTEND
            FC = 1.D0 / (AB(IJ, K) - AA(IJ, K) * AC(IJ, K-1))
            AC(IJ, K) = AC(IJ, K) * FC
            DO N = 1, NTDIM
               ADT(IJ, K, N) = (  ADT(IJ, K, N)
     &                          - AA(IJ, K) * ADT(IJ, K-1, N)) * FC
            END DO
  200    CONTINUE
  210 CONTINUE

      DO 330 N = 1, NTDIM
         DO 310 K = KEND-1, KSTR, -1
            DO 300 IJ = IJTSTR, IJTEND
               ADT(IJ, K, N) = ADT(IJ, K, N) 
     &                       - AC(IJ, K) * ADT(IJ, K+1, N)
  300       CONTINUE
  310    CONTINUE
  330 CONTINUE

      RETURN
      END
