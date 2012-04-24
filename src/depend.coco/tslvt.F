      SUBROUTINE SVTSET
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Time integration of the tracer equations.
c
c  HISTORY
c     '01.02.02  H.Hasumi: from COCO3.2
c     '01.07.23  H.Hasumi: bug fix (HZ should not be changed here)
c     '01.09.17  H.Hasumi: GAMMA for BBL
c     '01.12.07  H.Hasumi
c     '02.05.29  H.Nakano: tracer dimension
c     '02.06.04  H.Hasumi: 
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: arguments of CHEKIN
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8      HX(NXYDIM),     HZ(NXYDIM)
      REAL*8      FT(NXYDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)

      REAL*8  RGAMMA(NZDIM)
      SAVE RGAMMA

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      REAL*8      AA(NXYDIM, NZDIM),     AB(NXYDIM, NZDIM)
      REAL*8      AC(NXYDIM, NZDIM)
      COMMON /WORK/ FTX, FTY, FTZ
      EQUIVALENCE(FTX(1, 1, 1), AA(1, 1))
      EQUIVALENCE(FTY(1, 1, 1), AB(1, 1))
      EQUIVALENCE(FTZ(1, 1, 1), AC(1, 1))

      REAL*8      DH(NXYDIM),  HZBOT(NXYDIM),  HXBOT(NXYDIM)
      REAL*8     DZB,    DZT
      INTEGER     IJ,      K,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8   GAMMA(NZ)
      NAMELIST /NMACCT/ GAMMA

      DATA GAMMA / NZ*1.D0 /

      CALL REWNML(IFPAR, JFPAR)
      WRITE(JFPAR, *) '*** SVTSET ***'
      READ(IFPAR, NMACCT, END=99)
   99 CONTINUE
      WRITE(JFPAR, *) ' LEVEL     GAMMA'
      DO 10 K = 1, NZ
         WRITE(JFPAR, '(I6,F10.4)') K, GAMMA(K)
   10 CONTINUE

      DO 25 K = KSTR, KEND
         RGAMMA(K) = 1.D0 / GAMMA(K-KSTR+1)
   25 CONTINUE

      RETURN

c =====================================================================

      ENTRY SLVTRC(
     M                 TX,     HX,
     I                ADT,  DIFFZ,
     I                 FT,     HZ)

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      CALL CHEKIN(   FTX,  'FTX',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTY,  'FTY',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTZ,  'FTZ',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTX(1, 1, 2),  'FSX',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTY(1, 1, 2),  'FSY',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTZ(1, 1, 2),  'FSZ',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')

      DO IJ = IJTSTR, IJTEND
        HXBOT(IJ) = HX(IJ) + ZBOT
      ENDDO

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            AA(IJ, K) = 0.D0
            AB(IJ, K) = 0.D0
            AC(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO IJ = 1, NXYDIM
         DH(IJ) = HX(IJ) - HZ(IJ)
      END DO

      DO 810 K = KSTR, KSTR+KZ-1
         DO 800 IJ = IJTSTR, IJTEND
            AA(IJ, K) = - TS * DIFFZ(IJ, K) / DZ(IJ, K)
            AC(IJ, K) = - TS * DIFFZ(IJ, K+1) / DZ(IJ, K)
            AB(IJ, K) = (HX(IJ) + ZBOT) / ZBOT - AA(IJ, K) - AC(IJ, K)
  800    CONTINUE
  810 CONTINUE
      DO N = 1, NTDIM
         DO K = KSTR, KSTR+KZ-1
            DO IJ = IJTSTR, IJTEND
               ADT(IJ, K, N) = ADT(IJ, K, N)
     &                       - DH(IJ) / ZBOT / TS * TX(IJ, K, N)
            END DO
         END DO
      END DO
      DO 830 K = KSTR+KZ, KEND
         DO 820 IJ = IJTSTR, IJTEND
            AA(IJ, K) = - TS * DIFFZ(IJ, K) / DZ(IJ, K)
            AC(IJ, K) = - TS * DIFFZ(IJ, K+1) / DZ(IJ, K)
            AB(IJ, K) = 1.D0 - AA(IJ, K) - AC(IJ, K)
  820    CONTINUE
  830 CONTINUE

      DO IJ = IJSTR, IJEND
         ADT(IJ, KSTR, 1) = ADT(IJ, KSTR, 1)
     &                     + TX(IJ, KSTR, 1) * FT(IJ, 2) / ZBOT
     &                                       * AMSKT(IJ, KSTR)
      ENDDO

      CALL THOMAS(
     M               ADT,     AC,
     I                AA,     AB)

      DO N = 1, NTDIM
         DO K = KSTR, KEND
            DO IJ = IJTSTR, IJTEND
               ADT(IJ, K, N) = ADT(IJ, K, N) * RGAMMA(K) * AMSKT(IJ, K)
            END DO
         END DO
      END DO
#ifdef OPT_BBL
      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND
            K = MAX(NBOT(IJ), KSTR)
            ADT(IJ, K, N) = ADT(IJ, K, N) *
     &                (  1.D0
     &                 + AMSKTB(IJ) * (  GAMMA(K-KSTR+1) * RGAMMA(KEND)
     &                                 - 1.D0))
         END DO
      END DO
#endif

      DO 920 N = 1, NTDIM
         DO 910 K = KSTR, KEND
            DO 900 IJ = IJTSTR, IJTEND
               TX(IJ, K, N) = TX(IJ, K, N) + TS * ADT(IJ, K, N)
  900       CONTINUE
  910    CONTINUE
  920 CONTINUE

c     '12.01.30: removed
c      DO IJ = IJTSTR, IJTEND
c         HZBOT(IJ) = HX(IJ) + ZBOT
c         DH(IJ) = - TS * FT(IJ, 2) * AMSKT(IJ, KSTR)
c         HX(IJ) = HX(IJ) + DH(IJ)
c         HXBOT(IJ) = HX(IJ) + ZBOT
c         IF (DH(IJ) .GE. 0.D0) THEN
c            DZB = 0.D0
c            DZT = DH(IJ) * DS(KSTR+KZ-1)
c            DO K = KSTR+KZ-1, KSTR+1, -1
c               DO N = 1, NTDIM
c                  TX(IJ, K, N) = (  TX(IJ, K  , N) * DS(K) * HZBOT(IJ)
c     &                            - TX(IJ, K  , N) * DZB
c     &                            + TX(IJ, K-1, N) * DZT)
c     &                           / HXBOT(IJ) / DS(K)
c               END DO
c               DZB = DZT
c               DZT = DZB + DH(IJ) * DS(K-1)
c            END DO
c            TX(IJ, KSTR, 2) = (  TX(IJ, KSTR, 2) * DS(KSTR) * HZBOT(IJ)
c     &                         - TX(IJ, KSTR, 2) * DZB)
c     &                        / HXBOT(IJ) / DS(KSTR)
c            TX(IJ, KSTR, 1) = (  TX(IJ, KSTR, 1) * DS(KSTR) * HZBOT(IJ)
c     &                         - TX(IJ, KSTR, 1) * DZB
c     &                         + TX(IJ, KSTR, 1) * DZT
c     &                         + TS * FT(IJ, 1))
c     &                        / HXBOT(IJ) / DS(KSTR)
c            DO N = 3, NTDIM
c               TX(IJ, KSTR, N) = (  TX(IJ, KSTR, N)*DS(KSTR)*HZBOT(IJ)
c     &                            - TX(IJ, KSTR, N) * DZB
c     &                            + TX(IJ, KSTR, N) * DZT
c     &                            + TS * FT(IJ, N))
c     &                           / HXBOT(IJ) / DS(KSTR)
c            END DO
c         ELSE
c            DZT = - DH(IJ)
c            DZB = DZT + DH(IJ) * DS(KSTR)
c            TX(IJ, KSTR, 2) = (  TX(IJ, KSTR, 2) * DS(KSTR) * HZBOT(IJ)
c     &                         + TX(IJ, KSTR+1, 2) * DZB)
c     &                        / HXBOT(IJ) / DS(KSTR)
c            TX(IJ, KSTR, 1) = (  TX(IJ, KSTR, 1) * DS(KSTR) * HZBOT(IJ)
c     &                         + TX(IJ, KSTR+1, 1) * DZB
c     &                         - TX(IJ, KSTR, 1) * DZT
c     &                         + TS * FT(IJ, 1))
c     &                        / HXBOT(IJ) / DS(KSTR)
c            DO N = 3, NTDIM
c               TX(IJ, KSTR, N) = (  TX(IJ, KSTR, N)*DS(KSTR)*HZBOT(IJ)
c     &                            + TX(IJ, KSTR+1, N) * DZB
c     &                            - TX(IJ, KSTR  , N) * DZT
c     &                            + TS * FT(IJ, N))
c     &                           / HXBOT(IJ) / DS(KSTR)  
c            END DO            
c            DO K = KSTR+1, KSTR+KZ-2
c               DZT = DZB
c               DZB = DZT + DH(IJ) * DS(K)
c               DO N = 1, NTDIM
c                  TX(IJ, K, N) = (  TX(IJ, K, N) * DS(K) * HZBOT(IJ)
c     &                            + TX(IJ, K+1, N) * DZB
c     &                            - TX(IJ, K  , N) * DZT)
c     &                           / HXBOT(IJ) / DS(K)
c               END DO
c            END DO
c         END IF
c      END DO

      DO IJ = IJTSTR, IJTEND
         TX(IJ, KSTR, 1) = TX(IJ, KSTR, 1)
     &                  + TS * FT(IJ, 1) /  
     &                    HXBOT(IJ) / DS(KSTR) 
      ENDDO

      DO N = 3, NTDIM
         DO IJ = IJTSTR, IJTEND
            TX(IJ, KSTR, N) = TX(IJ, KSTR, N)
     &                + TS * FT(IJ, N) / HXBOT(IJ) / DS(KSTR)
         END DO
      END DO

      RETURN
      END
