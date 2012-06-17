      SUBROUTINE VELRDS(
     M                 UX,     VX,
     I                 GX,     GY,    AMV,
     I                 HX)
      use zocdim
      use zocgrd
      use zocmsk
      use utrdg
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Time integration of the equation of motion with semi-implicit
c treatment of the Coriolis term
c
c  HISTORY
c     '99.08.13  H.Hasumi: from CCSR2-MASK
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '01.05.02  H.Hasumi
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8     AMV(NXYDIM, NZDIM)
      REAL*8      HX(NXYDIM)

      REAL*8      AA(NXYDIM, NZDIM),     AC(NXYDIM, NZDIM)
      REAL*8     ABR(NXYDIM, NZDIM),    ABI(NXYDIM, NZDIM)
      REAL*8     ADR(NXYDIM, NZDIM),    ADI(NXYDIM, NZDIM)
      REAL*8  TDIFFZ(NXYDIM, NZDIM)
      REAL*8   HVBOT(NXYDIM)
c      COMMON /WORK/ AA, ABR, ABI, AC, ADR, ADI, TDIFFZ, HVBOT

      REAL*8     CF1,    CF2
      LOGICAL OFIRST
      SAVE CF1, CF2, OFIRST

      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      REAL*8     ACC,   AIMP
      NAMELIST /NMACCV/ ACC
      NAMELIST /NMIMPL/ AIMP

      DATA OFIRST / .TRUE. /
      DATA ACC / 1.D0 /
      DATA AIMP / 0.5D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMACCV, END=99)
 99      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMIMPL, END=98)
 98      CONTINUE
         WRITE(JFPAR, NMACCV)
         WRITE(JFPAR, NMIMPL)

         CF1   = AIMP / ACC
         CF2   = 1.0D+0 / ACC
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            AA    (IJ, K) = 0.D0
            ABR   (IJ, K) = 0.D0
            ABI   (IJ, K) = 0.D0
            AC    (IJ, K) = 0.D0
            ADR   (IJ, K) = 0.D0
            ADI   (IJ, K) = 0.D0
            TDIFFZ(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO IJ = IJVSTR, IJVEND
         HVBOT(IJ) = (  (HX(IJ)    + HX(IJ+LE) ) * DY(IJ)
     &                + (HX(IJ+LN) + HX(IJ+LNE)) * DY(IJ+LN))
     &               / DYM(IJ) * 0.25D0
     &             + ZBOT
         HVBOT(IJ) = 1.D0 / HVBOT(IJ)
      END DO

      DO 210 K = KSTR+1, KSTR+KZ-1
         DO 200 IJ = IJVSTR, IJVEND
            TDIFFZ(IJ, K) = TS * CF2 * AMV(IJ, K) * HVBOT(IJ) * RSM(K) *
     &                      AMFVZ(IJ, K)
  200    CONTINUE
  210 CONTINUE
      DO 230 K = KSTR+KZ, KEND
         DO 220 IJ = IJVSTR, IJVEND
            TDIFFZ(IJ, K) = TS * CF2 * AMV(IJ, K) / DZM(IJ, K) *
     &                      AMFVZ(IJ, K)
  220    CONTINUE
  230 CONTINUE

      DO 310 K = KSTR, KSTR+KZ-1
         DO 300 IJ = IJVSTR, IJVEND
            AA (IJ, K) = - TDIFFZ(IJ, K) * HVBOT(IJ) * RS(K)
            AC (IJ, K) = - TDIFFZ(IJ, K+1) * HVBOT(IJ) * RS(K)
            ABR(IJ, K) = 1.D0 - AA(IJ, K) - AC(IJ, K)
            ABI(IJ, K) = TS * CF1 * COR(IJ)
            ADR(IJ, K) = GX(IJ, K)
            ADI(IJ, K) = GY(IJ, K)
  300    CONTINUE
  310 CONTINUE
      DO 330 K = KSTR+KZ, KEND
         DO 320 IJ = IJVSTR, IJVEND
            AA (IJ, K) = - TDIFFZ(IJ, K) / DZV(IJ, K)
            AC (IJ, K) = - TDIFFZ(IJ, K+1) / DZV(IJ, K)
            ABR(IJ, K) = 1.D0 - AA(IJ, K) - AC(IJ, K)
            ABI(IJ, K) = TS * CF1 * COR(IJ)
            ADR(IJ, K) = GX(IJ, K)
            ADI(IJ, K) = GY(IJ, K)
  320    CONTINUE
  330 CONTINUE

      CALL THMASC(
     M               ADR,    ADI,
     I                AA,    ABR,    ABI,     AC)

      DO 410 K = KSTR, KEND
         DO 400 IJ = IJVSTR, IJVEND
            UX(IJ, K) = (  UX(IJ, K)
     &                   + TS * CF2 * ADR(IJ, K)) * AMSKV(IJ, K)
            VX(IJ, K) = (  VX(IJ, K)
     &                   + TS * CF2 * ADI(IJ, K)) * AMSKV(IJ, K)
  400    CONTINUE
  410 CONTINUE

      RETURN
      END
