      SUBROUTINE SVTSET
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Time integration of the tracer equations.
c
c  HISTORY
c     '02.10.10  H.Hasumi: from MIROC3.1-OMIP
c     '07.10.24  H.Hasumi: for SOM on COCO4.3
c     '07.11.28  Y.Komuro: tunnel diffusion from IcedCOCO4.1
c                          (originally Dr. Oka's)
c     '08.05.30  Y.Komuro: FTX2 etc. output code by Dr. Tatebe merged
c     '08.06.17  Y.Komuro: for COCO4.4; initial/final processing added
c     '08.07.25  Y.Komuro: initial/final processing added
c     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8      HX(NXYDIM),     HZ(NXYDIM)
      REAL*8      FT(NXYDIM, NTDIM)
      REAL*8   SWABS(NXYDIM),     FS(NXYDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8    SSFC(NXYDIM)
      REAL*8      AX(NXYDIM, 0:NIC)

      REAL*8  RGAMMA(NZDIM)
      REAL*8  SWCONV(NXYDIM, NZDIM)
      SAVE RGAMMA, SWCONV

c      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
c      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
c      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
c      REAL*8     WZC(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)

c      REAL*8     FTX2(NXYDIM, NZDIM, NTDIM)
c      REAL*8     FTY2(NXYDIM, NZDIM, NTDIM)
c      REAL*8     FTZ2(NXYDIM, NZDIM, NTDIM)
!      COMMON /WORK/ FTX, FTY, FTZ, WZC, RZM
c      COMMON /WORK/ FTX, FTY, FTZ, FTX2, FTY2, FTZ2, WZC, RZM
      REAL*8      AA(NXYDIM, NZDIM),     AB(NXYDIM, NZDIM)
      REAL*8      AC(NXYDIM, NZDIM)
c      EQUIVALENCE(FTX(1, 1, 1), AA(1, 1))
c      EQUIVALENCE(FTY(1, 1, 1), AB(1, 1))
c      EQUIVALENCE(FTZ(1, 1, 1), AC(1, 1))

      REAL*8      DH(NXYDIM),  HZBOT(NXYDIM),  HXBOT(NXYDIM)
      REAL*8     DZB,    DZT
      REAL*8   RADUP,  RADDN,  DEPTH
      REAL*8  SWCNV1(NZDIM), TSWCNV
      INTEGER     IJ,      K,      N,      L
      INTEGER  IFPAR,  JFPAR
      CHARACTER CTXNAM(NTDIM)*16, CTYNAM(NTDIM)*16, CTZNAM(NTDIM)*16
      SAVE CTXNAM, CTYNAM, CTZNAM

      REAL*8   GAMMA(NZ)
      REAL*8     RRR,  ZETA1,  ZETA2
      REAL*8    SDMP
      LOGICAL OSRSTR, OSRSTI
      SAVE SDMP, OSRSTR, OSRSTI
      NAMELIST /NMACCT/ GAMMA
      NAMELIST /NMSWAB/ RRR, ZETA1, ZETA2
      NAMELIST /NMSRST/ SDMP, OSRSTR, OSRSTI

      DATA GAMMA / NZ*1.D0 /
      DATA RRR, ZETA1, ZETA2 / 5.8D-1, 3.5D+1, 2.3D+3 /
      DATA SDMP / 0.D0 /
      DATA OSRSTR, OSRSTI / .FALSE., .FALSE. /

      CALL REWNML(IFPAR, JFPAR)
      WRITE(JFPAR, *) '*** SVTSET ***'
      READ(IFPAR, NMACCT, END=99)
   99 CONTINUE
      WRITE(JFPAR, *) ' LEVEL     GAMMA'
      DO 10 K = 1, NZ
         WRITE(JFPAR, '(I6,F10.4)') K, GAMMA(K)
   10 CONTINUE
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMSWAB, END=98)
 98   WRITE(JFPAR, NMSWAB)
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMSRST, END=97)
 97   WRITE(JFPAR, NMSRST)

      IF (OSRSTR) THEN
         SDMP = 1.D0 / 8.64D4 / SDMP
      END IF

      DO 25 K = KSTR, KEND
         RGAMMA(K) = 1.D0 / GAMMA(K-KSTR+1)
   25 CONTINUE

      DEPTH = 0.D0
      RADUP = 1.D0
      DO 50 K = KSTR, KEND
         DEPTH = DEPTH + DZ0(K)
         RADDN = RRR * EXP(- DEPTH / ZETA1)
     &         + (1.D0 - RRR) * EXP(- DEPTH / ZETA2)
         SWCNV1(K) = RADUP - RADDN
         RADUP = RADDN
   50 CONTINUE
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            SWCONV(IJ, K) = 0.D0
         END DO
      END DO
      DO 72 IJ = IJTSTR, IJTEND
         TSWCNV = 0.D0
         DO 70 K = KSTR, NBOT(IJ)-1
            SWCONV(IJ, K) = SWCNV1(K)
            TSWCNV = TSWCNV + SWCNV1(K)
   70    CONTINUE
         IF (NBOT(IJ) .GE. KSTR) THEN
            SWCONV(IJ, NBOT(IJ)) = 1.D0 - TSWCNV
         END IF
   72 CONTINUE
      CALL PUTSWC(
     I            SWCONV(1, KSTR))
      DO 90 IJ = 1, NXYDIM
         SWCONV(IJ, KSTR) = 0.D0
   90 CONTINUE
      DO K = KSTR+1, KEND
         DO IJ = IJTSTR, IJTEND
            SWCONV(IJ, K) = SWCONV(IJ, K) / rhoo / cpo
         END DO
      END DO

      CALL TTSSET

      DO L = 3, NTDIM
         WRITE(CTXNAM(L), '(A5,I2.2)') 'FTRCX', L
         WRITE(CTYNAM(L), '(A5,I2.2)') 'FTRCY', L
         WRITE(CTZNAM(L), '(A5,I2.2)') 'FTRCZ', L
      END DO

      RETURN

c =====================================================================

      ENTRY SLVTRC(
     M                 TX,     HX,
     I                ADT,  DIFFZ,
     I                 FT,  SWABS,     FS,     HZ,   SSFC,
     I                 AX)

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

c      moved to tflxt.F
c      CALL CHEKIN(   FTX,  'FTX',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTY,  'FTY',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTZ,  'FTZ',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTX2,  'FTX2',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTY2,  'FTY2',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTZ2,  'FTZ2',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTX(1, 1, 2),  'FSX',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTY(1, 1, 2),  'FSY',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTZ(1, 1, 2),  'FSZ',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTX2(1, 1, 2),  'FSX2',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTY2(1, 1, 2),  'FSY2',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(   FTZ2(1, 1, 2),  'FSZ2',
c     &                NX,     NY,     NZ, NXYZDM, 'OCN')
c      DO L=3, NTDIM
c         CALL CHEKIN(   FTX(1, 1, L), CTXNAM(L),
c     &                   NX,     NY,     NZ, NXYZDM,  'OCN')
c         CALL CHEKIN(   FTY(1, 1, L), CTYNAM(L),
c     &                   NX,     NY,     NZ, NXYZDM,  'OCN')
c         CALL CHEKIN(   FTZ(1, 1, L), CTZNAM(L),
c     &                   NX,     NY,     NZ, NXYZDM,  'OCN')
c      END DO

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

      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJTSTR, IJTEND
            TX(IJ, K, 1) = TX(IJ, K, 1)
     &                   + TS * SWCONV(IJ, K) * SWABS(IJ)
     &                     / HXBOT(IJ) / DS(K)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = IJTSTR, IJTEND
            TX(IJ, K, 1) = TX(IJ, K, 1)
     &                   + TS * SWCONV(IJ, K) * SWABS(IJ) / DZ(IJ, K)
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         TX(IJ, KSTR, 2) = TX(IJ, KSTR, 2)
     &                   - TS * FS(IJ) / HXBOT(IJ) / DS(KSTR)
      END DO

      DO N = 3, NTDIM
         DO IJ = IJTSTR, IJTEND
            TX(IJ, KSTR, N) = TX(IJ, KSTR, N)
     &                + TS * FT(IJ, N) / HXBOT(IJ) / DS(KSTR)
         END DO
      END DO

#ifdef OPT_SRST
      IF (OSRSTR) THEN
c         CALL TMINTP(SSFC, 10)
         IF (OSRSTI) THEN
            DO IJ = IJTSTR, IJTEND
               TX(IJ, KSTR, 2)
     &       = TX(IJ, KSTR, 2)
     &       + TS * SDMP * (SSFC(IJ) - TX(IJ, KSTR, 2)) *
     &         AMSKT(IJ, KSTR)
            END DO
         ELSE
            DO IJ = IJTSTR, IJTEND
               IF (AX(IJ, 0) .EQ. 1.D0) THEN
                  TX(IJ, KSTR, 2)
     &          = TX(IJ, KSTR, 2)
     &          + TS * SDMP * (SSFC(IJ) - TX(IJ, KSTR, 2)) *
     &            AMSKT(IJ, KSTR)
               END IF
            END DO
         END IF
      END IF
#endif

      RETURN
      END

#ifdef OPT_PARALLEL
c *********************************************************************

      SUBROUTINE TTSSET

      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      implicit none

c --- information -----------------------------------------------------
c
c  HISTORY
c     '04.04.13  A.Oka: from MIROC2.1
c     '07.11.28  Y.Komuro: for COCO4.3
c     '08.03.06  Y.Komuro: parallel code merged (from MIROC3.2)
c                          bug fix (LTNGG incorrect)
c
c ---------------------------------------------------------------------

#include "mpif.h"

      INTEGER   NTDIM1
      PARAMETER(NTDIM1 = NTDIM + 1)

      INTEGER NTUNNEL,NTNGR,NTIGR
      PARAMETER (NTUNNEL = 20)
      PARAMETER (NTNGR=NTUNNEL,NTIGR=NTUNNEL)
      INTEGER NTNI,NTN
      INTEGER LTNIG(NTUNNEL),LTNIGD(NTUNNEL),LTNIGG(NTUNNEL)
      INTEGER LTNG(NTUNNEL*NZDIM),LTNGD(NTUNNEL*NZDIM)
      INTEGER LTNGG(NTUNNEL*NZDIM)
      SAVE NTNI,LTNIG,LTNIGD,LTNIGG
      SAVE NTN,LTNG,LTNGD,LTNGG

      INTEGER  LTNGN(NTUNNEL*NZDIM),  LTNGR(NTUNNEL*NZDIM)
      INTEGER LTNGDN(NTUNNEL*NZDIM), LTNGDR(NTUNNEL*NZDIM)
      SAVE LTNGN, LTNGR, LTNGDN, LTNGDR

      REAL*8  BUF2D(NXYDIM)
      REAL*8  SNDBFI(NTDIM, NTUNNEL*2+1), RCVBFI(NTDIM, NTUNNEL*2+1)
      REAL*8  SNDBUF(NTDIM1, NTUNNEL*NZDIM*2+1)
      REAL*8  RCVBUF(NTDIM1, NTUNNEL*NZDIM*2+1)
      INTEGER ISTMPI(MPI_STATUS_SIZE)
      INTEGER  IRANKS,   NBUF
      INTEGER, ALLOCATABLE :: ISREQ(:), IRREQ(:)
      SAVE ISREQ, IRREQ

      INTEGER    LTNS(NTUNNEL*NZDIM*2+1), LTNSD(NTUNNEL*NZDIM*2+1)
      INTEGER   LTNSG(NTUNNEL*NZDIM*2+1)
      INTEGER    NTNS, NTNSST(0:NXYG)
      INTEGER   LTNIS(NTUNNEL*2+1), LTNISD(NTUNNEL*2+1)
      INTEGER  LTNISG(NTUNNEL*2+1)
      INTEGER   NTNIS, NTNIST(0:NXYG)
      SAVE LTNS, LTNSD, LTNSG, NTNS, NTNSST
      SAVE LTNIS, LTNISD, LTNISG, NTNIS, NTNIST
      REAL*8     BOTG(NXYGDM)

      REAL*8      TX(NXYZDM, NTDIM)
      REAL*8      HX(NXYDIM)

      REAL*8     VLM(NXYZDM),   AREA(NXYDIM)
      REAL*8    CRST(NTNGR),  CRSTI(NTIGR)
      SAVE VLM, AREA, CRST, CRSTI

      REAL*8     TSH(NXYDIM, NZDIM, NTDIM)
      REAL*8  RHXBOT(NXYDIM)
      REAL*8      AD(NTUNNEL*2+1),     HD(NTUNNEL*2+1)
      REAL*8    TSHD(NTUNNEL*2+1, NTDIM)
      REAL*8      VD(NTUNNEL*NZDIM*2+1),    TD(NTUNNEL*NZDIM*2+1, NTDIM)

      INTEGER     IJ,    IJK,      K,      N
      INTEGER      L,     LA,     LG,   ITRC
      INTEGER  IFPAR,  JFPAR
      INTEGER      I,      J,     IG,     JG,     M

      DATA CRST / NTNGR*0.D0 /
      DATA CRSTI / NTIGR*0.D0 /

      INTEGER NTUN,LTUN1(NTUNNEL),LTUN2(NTUNNEL),KTUN(NTUNNEL)

      INTEGER LEVTUN(NTUNNEL)
      REAL*8  DMPTUN(NTUNNEL),DSHTUN(NTUNNEL)
      INTEGER ITUN1(NTUNNEL),JTUN1(NTUNNEL)
      INTEGER ITUN2(NTUNNEL),JTUN2(NTUNNEL)
      NAMELIST /NMTUNL/ ITUN1,JTUN1,ITUN2,JTUN2,LEVTUN,DMPTUN,DSHTUN
      DATA LEVTUN / NTUNNEL*0 /
      DATA DMPTUN / NTUNNEL*30.0D0 /
      DATA DSHTUN / NTUNNEL*30.0D0 /
      DATA ITUN1 / NTUNNEL*0 /
      DATA JTUN1 / NTUNNEL*0 /
      DATA ITUN2 / NTUNNEL*0 /
      DATA JTUN2 / NTUNNEL*0 /

      ALLOCATE(ISREQ(0:NPROCS-1))
      ALLOCATE(IRREQ(0:NPROCS-1))

      DO IJ=1,NXYDIM
         BUF2D(IJ) = NBOT(IJ)
      ENDDO
      CALL GATHER_2D(BOTG, BUF2D)

C     --TUN DIF--
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMTUNL, END=95)
 95   CONTINUE
      WRITE(JFPAR,*) '*** TTSSET (TUNNEL DIFFUSION) ***'

      IF (MYRANK .EQ. IROOT) THEN
         NTUN = 0
         DO L=1,NTUNNEL
            IF (LEVTUN(L).EQ.0) GOTO 1999
            LTUN1(L) = IGSTR+ITUN1(L)-1 + NXGDIM*(JGSTR+JTUN1(L)-2)
            LTUN2(L) = IGSTR+ITUN2(L)-1 + NXGDIM*(JGSTR+JTUN2(L)-2)
            KTUN(L) = MIN( MIN( LEVTUN(L), INT(BOTG(LTUN1(L))) ), 
     &           INT(BOTG(LTUN2(L))) )
            NTUN = L
         ENDDO
 1999    CONTINUE
         NTNI = NTUN
         DO L=1,NTUN
            LTNIG(L) = LTUN1(L)
            LTNIGD(L) = LTUN2(L)
            LTNIGG(L) = L
         ENDDO
         NTN = 0
         DO L=1,NTUN
            DO K=KSTR+KZ,KSTR+KTUN(L)-1
               NTN = NTN+1
               LTNG(NTN) = LTUN1(L) + NXYGDM*(K-1)
               LTNGD(NTN) = LTUN2(L) + NXYGDM*(K-1)
               LTNGG(NTN) = L
            ENDDO
         ENDDO
         WRITE(JFPAR,*) '  TOTAL TUNNEL NUMBER = ',NTUN
         WRITE(JFPAR,*) '( I1, J1 ) .. ( I2, J2) KLEV DMP DMPSH'
         DO L=1,NTUN
            WRITE(JFPAR,*) '(',ITUN1(L),',',JTUN1(L),') .. (',
     &        ITUN2(L),',',JTUN2(L),') ',KTUN(L),DMPTUN(L),DSHTUN(L)
         ENDDO
      END IF
      N=1
      CALL MPI_BCAST(NTNI  ,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      IF (NTNI .GT. 0) THEN
      N=NTNI
      CALL MPI_BCAST(LTNIG ,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LTNIGD,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LTNIGG,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      END IF
      N=1
      CALL MPI_BCAST(NTN   ,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      IF (NTN .GT. 0) THEN
      N=NTN
      CALL MPI_BCAST(LTNG  ,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LTNGD ,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LTNGG ,N,MPI_INTEGER,IROOT,MPI_COMM_WORLD,IERR)
      END IF

      DO N = 1, NTN
         K  = (LTNG(N) - 1) / NXYGDM + 1
         JG = (LTNG(N) - (K - 1) * NXYGDM - 1) / NXGDIM + 1
         IG = LTNG(N) - (K - 1) * NXYGDM - (JG - 1) * NXGDIM
         J = MOD(JG-JGSTR, NY) + JSTR
         I = MOD(IG-IGSTR, NX) + ISTR
         LTNGN(N) = (K - 1) * NXYDIM + (J - 1) * NXDIM + I
         LTNGR(N) = (JG - JGSTR) / NY * INODES
     &            + (IG - IGSTR) / NX
      END DO
      DO N = 1, NTN
         K  = (LTNGD(N) - 1) / NXYGDM + 1
         JG = (LTNGD(N) - (K - 1) * NXYGDM - 1) / NXGDIM + 1
         IG = LTNGD(N) - (K - 1) * NXYGDM - (JG - 1) * NXGDIM
         J = MOD(JG-JGSTR, NY) + JSTR
         I = MOD(IG-IGSTR, NX) + ISTR
         LTNGDN(N) = (K - 1) * NXYDIM + (J - 1) * NXDIM + I
         LTNGDR(N) = (JG - JGSTR) / NY * INODES
     &             + (IG - IGSTR) / NX
      END DO
      IF (MYRANK .LT. IJNODE) THEN
         NTNS = 0
         NTNSST(0) = 1
         DO IRANKS = 0, IJNODE-1
            IF (IRANKS .EQ. MYRANK) THEN
               DO N = 1, NTN
                  IF (      (LTNGR(N)  .EQ. IRANKS)
     &                .AND. (LTNGDR(N) .EQ. IRANKS)) THEN
                     NTNS = NTNS + 1
                     LTNS (NTNS) = LTNGN(N)
                     LTNSD(NTNS) = LTNGDN(N)
                     LTNSG(NTNS) = LTNGG(N)
                     NTNS = NTNS + 1
                     LTNS (NTNS) = LTNGDN(N)
                     LTNSD(NTNS) = LTNGN(N)
                     LTNSG(NTNS) = LTNGG(N)
                  END IF
               END DO
            ELSE
               DO N = 1, NTN
                  IF (      (LTNGR(N)  .EQ. MYRANK)
     &                .AND. (LTNGDR(N) .EQ. IRANKS)) THEN
                     NTNS = NTNS + 1
                     LTNS (NTNS) = LTNGN(N)
                     LTNSD(NTNS) = -1
                     LTNSG(NTNS) = LTNGG(N)
                  ELSE IF (      (LTNGR(N)  .EQ. IRANKS)
     &                     .AND. (LTNGDR(N) .EQ. MYRANK)) THEN
                     NTNS = NTNS + 1
                     LTNS (NTNS) = LTNGDN(N)
                     LTNSD(NTNS) = -1
                     LTNSG(NTNS) = LTNGG(N)
                  END IF
               END DO
            END IF
            NTNSST(IRANKS+1) = NTNS + 1
         END DO
      END IF

      DO N = 1, NTNI
         JG = (LTNIG(N) - 1) / NXGDIM + 1
         IG = LTNIG(N) - (JG - 1) * NXGDIM
         J = MOD(JG-JGSTR, NY) + JSTR
         I = MOD(IG-IGSTR, NX) + ISTR
         LTNGN(N) = (J - 1) * NXDIM + I
         LTNGR(N) = (JG - JGSTR) / NY * INODES
     &            + (IG - IGSTR) / NX
      END DO
      DO N = 1, NTNI
         JG = (LTNIGD(N) - 1) / NXGDIM + 1
         IG = LTNIGD(N) - (JG - 1) * NXGDIM
         J = MOD(JG-JGSTR, NY) + JSTR
         I = MOD(IG-IGSTR, NX) + ISTR
         LTNGDN(N) = (J - 1) * NXDIM + I
         LTNGDR(N) = (JG - JGSTR) / NY * INODES
     &             + (IG - IGSTR) / NX
      END DO
      IF (MYRANK .LT. IJNODE) THEN
         NTNIS = 0
         NTNIST(0) = 1
         DO IRANKS = 0, IJNODE-1
            IF (IRANKS .EQ. MYRANK) THEN
               DO N = 1, NTNI
                  IF (      (LTNGR(N)  .EQ. IRANKS)
     &                .AND. (LTNGDR(N) .EQ. IRANKS)) THEN
                     NTNIS = NTNIS + 1
                     LTNIS (NTNIS) = LTNGN(N)
                     LTNISD(NTNIS) = LTNGDN(N)
                     LTNISG(NTNIS) = LTNIGG(N)
                     NTNIS = NTNIS + 1
                     LTNIS (NTNIS) = LTNGDN(N)
                     LTNISD(NTNIS) = LTNGN(N)
                     LTNISG(NTNIS) = LTNIGG(N)
                  END IF
               END DO
            ELSE
               DO N = 1, NTNI
                  IF (      (LTNGR(N)  .EQ. MYRANK)
     &                .AND. (LTNGDR(N) .EQ. IRANKS)) THEN
                     NTNIS = NTNIS + 1
                     LTNIS (NTNIS) = LTNGN(N)
                     LTNISD(NTNIS) = -1
                     LTNISG(NTNIS) = LTNIGG(N)
                  ELSE IF (      (LTNGR(N)  .EQ. IRANKS)
     &                     .AND. (LTNGDR(N) .EQ. MYRANK)) THEN
                     NTNIS = NTNIS + 1
                     LTNIS (NTNIS) = LTNGDN(N)
                     LTNISD(NTNIS) = -1
                     LTNISG(NTNIS) = LTNIGG(N)
                  END IF
               END DO
            END IF
            NTNIST(IRANKS+1) = NTNIS + 1
         END DO
      END IF

      DO N=1,NTNGR
         CRST(N) = DMPTUN(N)
      ENDDO
      DO N=1,NTIGR
         CRSTI(N) = DSHTUN(N)
      ENDDO

      DO 30 K = KSTR+KZ, KEND
         DO 32 IJ = IJSTR, IJEND
            IJK = (K - 1) * NXYDIM + IJ
            VLM(IJK) = HXT(IJ) * HYT(IJ) * DZ(IJ, K)
   32    CONTINUE
   30 CONTINUE

      DO IJ = IJSTR, IJEND
         AREA(IJ) = HXT(IJ) * HYT(IJ)
      END DO
      DO N = 1, NTNGR
         CRST(N) = 1.D0 / CRST(N) / 86400.D0
      END DO
      DO N = 1, NTIGR
         CRSTI(N) = 1.D0 / CRSTI(N) / 86400.D0
      END DO

      RETURN

c =====================================================================

      ENTRY TUNDIF(
     M                   TX,     HX)

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      DO N = 1, NTDIM
         DO K = KSTR, KSTR+KZ-1
            DO IJ = IJSTR, IJEND
               IJK = (K - 1) * NXYDIM + IJ
               TSH(IJ, K, N) = TX(IJK, N) * (HX(IJ) + ZBOT)
            END DO
         END DO
      END DO

      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .NE. MYRANK) THEN
            IF (NTNIST(IRANKS+1) .GT. NTNIST(IRANKS)) THEN
               NBUF = NTNIST(IRANKS+1) - NTNIST(IRANKS)
               DO N = NTNIST(IRANKS), NTNIST(IRANKS+1)-1
                  SNDBFI(1, N) = HX(LTNIS(N))
                  SNDBFI(2, N) = AREA(LTNIS(N))
               END DO
               CALL MPI_ISEND(
     &                        SNDBFI(1, NTNIST(IRANKS)),
c_bugfix_20060829
c     &                        2*NBUF, MPI_DOUBLE_PRECISION,
     &                        NTDIM*NBUF, MPI_DOUBLE_PRECISION,
c_bugfix_20060829
     &                         IRANKS,      1, MPI_COMM_WORLD,
     &                         ISREQ(IRANKS),   IERR)
               CALL MPI_IRECV(
     &                        RCVBFI(1, NTNIST(IRANKS)),
c_bugfix_20060829
c     &                        2*NBUF, MPI_DOUBLE_PRECISION,
     &                        NTDIM*NBUF, MPI_DOUBLE_PRECISION,
c_bugfix_20060829
     &                         IRANKS,      1, MPI_COMM_WORLD,
     &                         IRREQ(IRANKS),   IERR)
            END IF
         END IF
      END DO
      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .NE. MYRANK) THEN
            IF (NTNIST(IRANKS+1) .GT. NTNIST(IRANKS)) THEN
               CALL MPI_WAIT(ISREQ(IRANKS), ISTMPI, IERR)
               CALL MPI_WAIT(IRREQ(IRANKS), ISTMPI, IERR)
            END IF
         END IF
      END DO
      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .EQ. MYRANK) THEN
            DO N = NTNIST(IRANKS), NTNIST(IRANKS+1)-1
               HD(N) = HX(LTNISD(N))
               AD(N) = AREA(LTNISD(N))
            END DO
         ELSE
            DO N = NTNIST(IRANKS), NTNIST(IRANKS+1)-1
               HD(N) = RCVBFI(1, N)
               AD(N) = RCVBFI(2, N)
            END DO
         END IF
      END DO
*POPTION INDEP(HX)
      DO N = 1, NTNIS
         L  = LTNIS (N)
         LG = LTNISG(N)
         HX(L) = HX(L)
     &         - (HX(L) - HD(N)) * TS * CRSTI(LG) *
     &           (AREA(L) + AD(N)) * 0.5D0 / AREA(L)
      END DO

      DO IJ = IJSTR, IJEND
         RHXBOT(IJ) = 1.D0 / (HX(IJ) + ZBOT)
      END DO

      DO K = KSTR, KSTR+KZ-1

      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .NE. MYRANK) THEN
            IF (NTNIST(IRANKS+1) .GT. NTNIST(IRANKS)) THEN
               NBUF = NTNIST(IRANKS+1) - NTNIST(IRANKS)
               DO N = NTNIST(IRANKS), NTNIST(IRANKS+1)-1
                  DO L = 1, NTDIM
                     SNDBFI(L, N) = TSH(LTNIS(N), K, L)
                  END DO
               END DO
               CALL MPI_ISEND(
     &                        SNDBFI(1, NTNIST(IRANKS)),
     &                        NTDIM*NBUF, MPI_DOUBLE_PRECISION,
     &                         IRANKS,      1, MPI_COMM_WORLD,
     &                         ISREQ(IRANKS),   IERR)
               CALL MPI_IRECV(
     &                        RCVBFI(1, NTNIST(IRANKS)),
     &                        NTDIM*NBUF, MPI_DOUBLE_PRECISION,
     &                         IRANKS,      1, MPI_COMM_WORLD,
     &                         IRREQ(IRANKS),   IERR)
            END IF
         END IF
      END DO
      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .NE. MYRANK) THEN
            IF (NTNIST(IRANKS+1) .GT. NTNIST(IRANKS)) THEN
               CALL MPI_WAIT(ISREQ(IRANKS), ISTMPI, IERR)
               CALL MPI_WAIT(IRREQ(IRANKS), ISTMPI, IERR)
            END IF
         END IF
      END DO
      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .EQ. MYRANK) THEN
            DO L = 1, NTDIM
               DO N = NTNIST(IRANKS), NTNIST(IRANKS+1)-1
                  TSHD(N, L) = TSH(LTNISD(N), K, L)
               END DO
            END DO
         ELSE
            DO N = NTNIST(IRANKS), NTNIST(IRANKS+1)-1
               DO L = 1, NTDIM
                  TSHD(N, L) = RCVBFI(L, N)
               END DO
            END DO
         END IF
      END DO
      DO M = 1, NTDIM
*POPTION INDEP(TSH)
         DO N = 1, NTNIS
            L  = LTNIS (N)
            LG = LTNISG(N)
            TSH(L, K, M) = TSH(L, K, M)
     &                   - (TSH(L, K, M) - TSHD(N, M)) * TS *
     &                     CRSTI(LG) * (AREA(L) + AD(N)) * 0.5D0
     &                     / AREA(L)
         END DO
         DO IJ = IJSTR, IJEND
            IJK = (K - 1) * NXYDIM + IJ
            TX(IJK, M) = TSH(IJ, K, M) * RHXBOT(IJ)
         END DO
      END DO

      END DO

      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .NE. MYRANK) THEN
            IF (NTNSST(IRANKS+1) .GT. NTNSST(IRANKS)) THEN
               NBUF = NTNSST(IRANKS+1) - NTNSST(IRANKS)
               DO N = NTNSST(IRANKS), NTNSST(IRANKS+1)-1
                  DO L = 1, NTDIM
                     SNDBUF(L, N) = TX(LTNS(N), L)
                  END DO
                  SNDBUF(NTDIM1, N) = VLM(LTNS(N))
               END DO
               CALL MPI_ISEND(
     &                        SNDBUF(1, NTNSST(IRANKS)),
     &                        NTDIM1*NBUF, MPI_DOUBLE_PRECISION,
     &                         IRANKS,      1, MPI_COMM_WORLD,
     &                         ISREQ(IRANKS),   IERR)
               CALL MPI_IRECV(
     &                        RCVBUF(1, NTNSST(IRANKS)),
     &                        NTDIM1*NBUF, MPI_DOUBLE_PRECISION,
     &                         IRANKS,      1, MPI_COMM_WORLD,
     &                         IRREQ(IRANKS),   IERR)
            END IF
         END IF
      END DO
      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .NE. MYRANK) THEN
            IF (NTNSST(IRANKS+1) .GT. NTNSST(IRANKS)) THEN
               CALL MPI_WAIT(ISREQ(IRANKS), ISTMPI, IERR)
               CALL MPI_WAIT(IRREQ(IRANKS), ISTMPI, IERR)
            END IF
         END IF
      END DO

      DO IRANKS = 0, IJNODE-1
         IF (IRANKS .EQ. MYRANK) THEN
            DO N = NTNSST(IRANKS), NTNSST(IRANKS+1)-1
               DO L = 1, NTDIM
                  TD(N, L) = TX(LTNSD(N), L)
               END DO
               VD(N) = VLM(LTNSD(N))
            END DO
         ELSE
            DO N = NTNSST(IRANKS), NTNSST(IRANKS+1)-1
               DO L = 1, NTDIM
                  TD(N, L) = RCVBUF(L, N)
               END DO
               VD(N) = RCVBUF(NTDIM1, N)
            END DO
         END IF
      END DO

      DO M = 1, NTDIM
*POPTION INDEP(TX)
         DO 1000 N = 1, NTNS
            L  = LTNS (N)
            LG = LTNSG(N)
            TX(L, M) = TX(L, M)
     &               - (TX(L, M) - TD(N, M)) * TS * CRST(LG) *
     &                 (VLM(L) + VD(N)) * 0.5D0 / VLM(L)
 1000    CONTINUE
      END DO

      RETURN
      END
#else
c *********************************************************************

      SUBROUTINE TTSSET

      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      implicit none
c --- information -----------------------------------------------------
c
c  HISTORY
c     '04.04.13  A.Oka: from MIROC2.1
c     '07.11.28  Y.Komuro: for COCO4.3
c     '08.03.06  Y.Komuro: non-parallel only
c                          bug fix (LTNGG incorrect)
c
c ---------------------------------------------------------------------

      INTEGER NTUNNEL,NTNGR,NTIGR
      PARAMETER (NTUNNEL = 20)
      PARAMETER (NTNGR=NTUNNEL,NTIGR=NTUNNEL)
      INTEGER NTNI,NTN
      INTEGER LTNIG(NTUNNEL),LTNIGD(NTUNNEL),LTNIGG(NTUNNEL)
      INTEGER LTNG(NTUNNEL*NZDIM),LTNGD(NTUNNEL*NZDIM)
      INTEGER LTNGG(NTUNNEL*NZDIM)
      SAVE NTNI,LTNIG,LTNIGD,LTNIGG
      SAVE NTN,LTNG,LTNGD,LTNGG

      REAL*8      TXG(NXYZDM, NTDIM)
      REAL*8      HXG(NXYDIM),     DZG(NXYDIM,NZDIM)

      REAL*8      TX(NXYZDM, NTDIM)
      REAL*8      HX(NXYDIM)

      REAL*8     VLM(NXYZDM),   AREA(NXYDIM)
      REAL*8    CRST(NTNGR),  CRSTI(NTIGR)
      SAVE VLM, AREA, CRST

      REAL*8     DIFFSH(NXYDIM),  DIFFT(NXYZDM, NTDIM)
      REAL*8     TSH(NXYDIM, NZDIM, NTDIM)
      REAL*8  RHXBOT(NXYDIM)
C      COMMON /WORK/ DIFFT, DIFFSH, TSH, RHXBOT

      INTEGER     IJ,    IJK,      K,      N
      INTEGER      L,     LA,     LG,   ITRC
      INTEGER  IFPAR,  JFPAR

      NAMELIST /NMOTNG/ CRST
      NAMELIST /NMITNG/ CRSTI
      DATA CRST / NTNGR*0.D0 /
      DATA CRSTI / NTIGR*0.D0 /

      INTEGER NTUN,LTUN1(NTUNNEL),LTUN2(NTUNNEL),KTUN(NTUNNEL)

      INTEGER LEVTUN(NTUNNEL)
      REAL*8  DMPTUN(NTUNNEL),DSHTUN(NTUNNEL)
      INTEGER ITUN1(NTUNNEL),JTUN1(NTUNNEL)
      INTEGER ITUN2(NTUNNEL),JTUN2(NTUNNEL)
      NAMELIST /NMTUNL/ ITUN1,JTUN1,ITUN2,JTUN2,LEVTUN,DMPTUN,DSHTUN
      DATA LEVTUN / NTUNNEL*0 /
      DATA DMPTUN / NTUNNEL*30.0D0 /
      DATA DSHTUN / NTUNNEL*30.0D0 /
      DATA ITUN1 / NTUNNEL*0 /
      DATA JTUN1 / NTUNNEL*0 /
      DATA ITUN2 / NTUNNEL*0 /
      DATA JTUN2 / NTUNNEL*0 /

C     --TUN DIF--
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMTUNL, END=95)
 95   CONTINUE
      NTUN = 0
      DO L=1,NTUNNEL
         IF (LEVTUN(L).EQ.0) GOTO 1999
         LTUN1(L) = ISTR+ITUN1(L)-1 + NXDIM*(JSTR+JTUN1(L)-2)
         LTUN2(L) = ISTR+ITUN2(L)-1 + NXDIM*(JSTR+JTUN2(L)-2)
         KTUN(L) = MIN( MIN( LEVTUN(L), INT(NBOT(LTUN1(L))) ), 
     &        INT(NBOT(LTUN2(L))) )
         NTUN = L
      ENDDO
 1999 CONTINUE
      NTNI = NTUN
      DO L=1,NTUN
         LTNIG(L) = LTUN1(L)
         LTNIGD(L) = LTUN2(L)
         LTNIGG(L) = L
      ENDDO
      NTN = 0
      DO L=1,NTUN
         DO K=KSTR+KZ,KSTR+KTUN(L)-1
            NTN = NTN+1
            LTNG(NTN) = LTUN1(L) + NXYDIM*(K-1)
            LTNGD(NTN) = LTUN2(L) + NXYDIM*(K-1)
            LTNGG(NTN) = L
         ENDDO
      ENDDO
      WRITE(JFPAR,*) '*** TTSSET (TUNNEL DIFFUSION) ***'
      WRITE(JFPAR,*) '  TOTAL TUNNEL NUMBER = ',NTUN
      WRITE(JFPAR,*) '( I1, J1 ) .. ( I2, J2) KLEV DMP DMPSH'
      DO L=1,NTUN
         WRITE(JFPAR,*) '(',ITUN1(L),',',JTUN1(L),') .. (',
     &        ITUN2(L),',',JTUN2(L),') ',KTUN(L),DMPTUN(L),DSHTUN(L)
      ENDDO
      DO N=1,NTNGR
         CRST(N) = DMPTUN(N)
      ENDDO
      DO N=1,NTIGR
         CRSTI(N) = DSHTUN(N)
      ENDDO

      DO 30 K = KSTR+KZ, KEND
         DO 32 IJ = IJSTR, IJEND
            IJK = (K - 1) * NXYDIM + IJ
            VLM(IJK) = HXT(IJ) * HYT(IJ) * DZ(IJ, K)
   32    CONTINUE
   30 CONTINUE

      DO IJ = IJSTR, IJEND
         AREA(IJ) = HXT(IJ) * HYT(IJ)
      END DO
      DO N = 1, NTNGR
         CRST(N) = 1.D0 / CRST(N) / 86400.D0
      END DO
      DO N = 1, NTIGR
         CRSTI(N) = 1.D0 / CRSTI(N) / 86400.D0
      END DO

      RETURN

c =====================================================================

      ENTRY TUNDIF(
     M                   TX,     HX)

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF
c      IF (NTUN.EQ.0) GOTO 999

      DO ITRC = 1, NTDIM
      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJSTR, IJEND
            IJK = (K - 1) * NXYDIM + IJ
            TSH(IJ, K, ITRC) = TX(IJK, ITRC) * (HX(IJ) + ZBOT)
         END DO
      END DO
      END DO

*POPTION INDEP(DIFFSH)
      DO N = 1, NTNI
         L  = LTNIG (N)
         LA = LTNIGD(N)
         LG = LTNIGG(N)
         DIFFSH(L) = (HX(L) - HX(LA)) * TS * CRSTI(LG)
      END DO
*POPTION INDEP(HX)
      DO N = 1, NTNI
         L  = LTNIG (N)
         LA = LTNIGD(N)
         HX(L) = HX(L)
     &         - DIFFSH(L) * (AREA(L) + AREA(LA)) * 0.5D0 / AREA(L)
      END DO
*POPTION INDEP(HX)
      DO N = 1, NTNI
         L  = LTNIG (N)
         LA = LTNIGD(N)
         HX(LA) = HX(LA)
     &          + DIFFSH(L) * (AREA(L) + AREA(LA)) * 0.5D0 / AREA(LA)
      END DO

      DO IJ = IJSTR, IJEND
         RHXBOT(IJ) = 1.D0 / (HX(IJ) + ZBOT)
      END DO

      DO K = KSTR, KSTR+KZ-1
*POPTION INDEP(DIFFT)
         DO N = 1, NTNI
            L  = LTNIG (N)
            LA = LTNIGD(N)
            LG = LTNIGG(N)
            DO ITRC = 1, NTDIM
               DIFFT(L, ITRC) = (TSH(L, K, ITRC)
     &                      - TSH(LA, K, ITRC)) * TS * CRSTI(LG)
            ENDDO
         END DO
*POPTION INDEP(TSH)
         DO N = 1, NTNI
            L  = LTNIG (N)
            LA = LTNIGD(N)
            DO ITRC = 1, NTDIM
               TSH(L, K, ITRC) = TSH(L, K, ITRC)
     &                - DIFFT(L, ITRC) * (AREA(L) + AREA(LA)) *
     &                  0.5D0 / AREA(L)
            ENDDO
         END DO
*POPTION INDEP(TSH)
         DO N = 1, NTNI
            L  = LTNIG (N)
            LA = LTNIGD(N)
            DO ITRC = 1, NTDIM
               TSH(LA, K, ITRC) = TSH(LA, K, ITRC)
     &                 + DIFFT(L, ITRC) * (AREA(L) + AREA(LA)) *
     &                   0.5D0 / AREA(LA)
            ENDDO
         END DO

         DO IJ = IJSTR, IJEND
            IJK = (K - 1) * NXYDIM + IJ
            DO ITRC = 1, NTDIM
               TX(IJK, ITRC) = TSH(IJ, K, ITRC) * RHXBOT(IJ)
            ENDDO
         END DO
      END DO

*POPTION INDEP(DIFFT)
      DO 1100 N = 1, NTN
         L  = LTNG (N)
         LA = LTNGD(N)
         LG = LTNGG(N)
         DO ITRC = 1, NTDIM
            DIFFT(L, ITRC) = (TX(L, ITRC)
     &           - TX(LA, ITRC)) * TS * CRST(LG)
         ENDDO
 1100 CONTINUE
*POPTION INDEP(TX)
      DO 1200 N = 1, NTN
         L  = LTNG (N)
         LA = LTNGD(N)
         DO ITRC= 1, NTDIM
            TX(L, ITRC) = TX(L, ITRC)
     &         - (VLM(L) + VLM(LA)) * 0.5D0 / VLM(L) * DIFFT(L, ITRC)
         ENDDO
 1200 CONTINUE
*POPTION INDEP(TX)
      DO 1300 N = 1, NTN
         L  = LTNG (N)
         LA = LTNGD(N)
         DO ITRC = 1, NTDIM
            TX(LA, ITRC) = TX(LA, ITRC)
     &          + (VLM(L) + VLM(LA)) * 0.5D0 
     &            / VLM(LA) * DIFFT(L, ITRC)
         ENDDO
 1300 CONTINUE

 999  CONTINUE
      RETURN
      END
#endif
