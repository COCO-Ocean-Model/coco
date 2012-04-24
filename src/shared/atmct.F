      SUBROUTINE TMSTUP(
     O             TSTRT,   TEND,    DTT)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Controling the scheme of time integration and the timing of output
c
c  HISTORY
c     '00.05.30  H.Hasumi: parallelized COCO3
c     '00.12.07  H.Hasumi: combine par and non-par routines
c     '01.12.07  H.Hasumi
c     '03.06.04  H.Hasumi: allow to specify the time step by the unit
c                          other than 'hour'
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

#include "zocfil.F"
#ifdef OPT_PARALLEL
#include "mpif.h"
#include "zocnod.F"
#endif

      INTEGER   NTSCYC
      PARAMETER(NTSCYC = 100)

      REAL*8   TSTRT,   TEND,    DTT,     TT

      REAL*8      TS,    TSS
      INTEGER     NT,    ITS,   ITST,   NTSS
      LOGICAL OFLOUT(NFOMAX), OFLSTK(NFOMAX)
      LOGICAL ORSOUT, ORSRWD

      REAL*8  DT
      INTEGER I, IITEM
      INTEGER ITSCYC, ITSTT(0:NTSCYC-1)
      INTEGER IONEXT(6, NFOMAX)
      INTEGER IRNEXT(6)
      INTEGER IFLOUT(NFOMAX)
      REAL*8  TOSTRT(NFOMAX),  TOEND(NFOMAX), TONEXT(NFOMAX)
      REAL*8  TRNEXT
      CHARACTER CUNIT(6)*6
      INTEGER IOHITM
      CHARACTER CITEM(NFOMAX)*16

      INTEGER IEB, ILEAP
      INTEGER ITSTRT(6), ITEND(6)
      INTEGER ICALN
      INTEGER IRINTV, IURINT, IRSRWD
      INTEGER IOSTRT(6, NFOMAX), IOEND(6, NFOMAX)
      INTEGER IOINTV(NFOMAX), IUINTV(NFOMAX)
      INTEGER NTSPLT
      INTEGER IUTSTP, IDATE(6)
      REAL*8  TMSTP
      INTEGER NOHITM

      SAVE    DT, ITSTT, ITSCYC
      SAVE    IFLOUT
      SAVE    IOSTRT, IOEND, IOINTV, IUINTV, IONEXT
      SAVE    IRINTV, IURINT, IRNEXT, IRSRWD
      SAVE    TOSTRT, TOEND, TONEXT
      SAVE    TRNEXT
      SAVE    NTSPLT
      SAVE    NOHITM

      INTEGER IODSTR(6), IODEND(6), IODINT, IUDINT, IODAVR, IODSNG
      INTEGER IOHSTR(6), IOHEND(6), IOHINT, IUHINT, IOHAVR, IOHSNG
      INTEGER IOXSTR, IOXEND, IOYSTR, IOYEND, IOZSTR, IOZEND
      CHARACTER COHFIL*(NCF), COHITM*16
      INTEGER  IFPAR,  JFPAR

      NAMELIST /NMSCHM/ IEB, ILEAP
      NAMELIST /NMTIME/ ITSTRT, ITEND, TMSTP, IUTSTP, NTSPLT
      NAMELIST /NMCALN/ ICALN
      NAMELIST /NMRSTR/ IRINTV, IURINT, IRSRWD
      NAMELIST /NMDOUT/ IODSTR, IODEND, IODINT, IUDINT, IODAVR, IODSNG
      NAMELIST /NMHIST/ COHITM, COHFIL,
     &                  IOHSTR, IOHEND, IOHINT, IUHINT, IOHAVR, IOHSNG,
     &                  IOXSTR, IOXEND, IOYSTR, IOYEND, IOZSTR, IOZEND

      DATA IEB, ILEAP / 9, 10 /
      DATA ITSTRT / 0, 0, 0, 0, 0, 0 /
      DATA ITEND  / 0, 0, 0, 0, 0, 0 /
      DATA IUTSTP / 0 /
      DATA TMSTP  / 1.D10 /
      DATA NTSPLT / 0 /
      DATA ICALN  / 0 /
      DATA IRINTV, IURINT, IRSRWD / 1, 1, 1 /
      DATA IODSTR / 0, 0, 0, 0, 0, 0 /
      DATA IODEND / 0, 0, 0, 0, 0, 0 /
      DATA IODINT, IUDINT / 1, 1 /
      DATA IODAVR, IODSNG / 1, 1 /
      DATA IFLOUT / NFOMAX*0 /
      DATA COHITM / 'not-specified' /

      DATA CUNIT / 'year', 'month', 'day', 'hour', 'minute', 'second' /
      DATA IDATE / 0, 1, 1, 0, 0, 0 /

      ILEAP  = NTSCYC - IEB - 1

      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMSCHM, END=99)
   99 CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMTIME, END=98)
   98 CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMCALN, END=97)
   97 CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMRSTR, END=96)
   96 CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMDOUT, END=95)
 95   CONTINUE
      DO IITEM = 1, NFOMAX
         DO I = 1, 6
            IOSTRT(I, IITEM) = IODSTR(I)
            IOEND (I, IITEM) = IODEND(I)
         END DO
         IOINTV(IITEM) = IODINT
         IUINTV(IITEM) = IUDINT
      END DO
      CALL REWNML(IFPAR, JFPAR)
      NOHITM = 0
 90   CONTINUE
      COHITM = 'not-specified'
      IOHSTR(1) = -1
      IOHEND(1) = -1
      IOHINT = -1
      IUHINT = -1
      READ(IFPAR, NMHIST, END=94)
      IF (COHITM(1:13) .NE. 'not-specified') THEN
         NOHITM = NOHITM + 1
         IOHITM = NOHITM
         IFLOUT(IOHITM) = 1
         CITEM(IOHITM) = COHITM
         IF (IOHSTR(1) .GE. 0) THEN
            DO I = 1, 6
               IOSTRT(I, IOHITM) = IOHSTR(I)
            END DO
         END IF
         IF (IOHEND(1) .GE. 0) THEN
            DO I = 1, 6
               IOEND(I, IOHITM) = IOHEND(I)
            END DO
         END IF
         IF (IOHINT .GT. 0) THEN
            IOINTV(IOHITM) = IOHINT
         END IF
         IF (IUHINT .GT. 0) THEN
            IUINTV(IOHITM) = IUHINT
         END IF
      END IF
      GO TO 90
 94   CONTINUE

      ITSCYC= IEB + 1 + ILEAP
      IF ( ITSCYC .GT. NTSCYC ) THEN
         WRITE(JFPAR, *) '### Workspace Over : NTSCYC (TMSTUP) ###'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
      END IF

c     ***** Making scheme table *****
c
c     ITSTT = 1 : Euler-backward
c             2 : Forward
c             3 : Leapfrog B --> A
c             4 : Leapfrog A --> B

      ITSTT(0) = 4
      DO 10 I = 1, IEB
         ITSTT(I) = 1
   10 CONTINUE
      ITSTT(IEB+1) = 2
      DO 20 I = IEB+2, ITSCYC-2, 2
         ITSTT(I  ) = 3
         ITSTT(I+1) = 4
   20 CONTINUE
      ITSTT(ITSCYC-1) = 3

c     ***** Reading time control parameters *****

      CALL CALNDR(
     I             ICALN)
      CALL CYH2SS(
     O             TSTRT,
     I            ITSTRT)
      CALL CYH2SS(
     O              TEND,
     I             ITEND)
      DO I = 1, 6
         IRNEXT(I) = ITSTRT(I)
      END DO
      DO IITEM = 1, NOHITM
         CALL CYH2SS(
     O            TOSTRT(IITEM),
     I            IOSTRT(1, IITEM))
         CALL CYH2SS(
     O             TOEND(IITEM),
     I             IOEND(1, IITEM))
         IF (TOSTRT(IITEM) .LT. TSTRT) THEN
            TOSTRT(IITEM) = TSTRT
            CALL CSS2YH(
     O               IOSTRT(1, IITEM),
     I               TOSTRT(IITEM))
         END IF
         DO I = 1, 6
            IONEXT(I, IITEM) = IOSTRT(I, IITEM)
         END DO
      END DO

      IF (IUTSTP .EQ. 0) THEN
         DT = ANINT(TMSTP * 3.6D+3)
         DTT = DT
      ELSE IF ((IUTSTP .GE. 1) .AND. (IUTSTP .LE. 6)) THEN
         IDATE(IUTSTP) = IDATE(IUTSTP) + NINT(TMSTP)
         CALL CYH2SS(
     O                   DT,
     I                IDATE)
         DTT = DT
      ELSE
         WRITE(JFPAR, *) '*** TMSTUP: NO SUCH UNIT OF TIME ***'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
      END IF

      DO IITEM = 1, NOHITM
         IF ((IUINTV(IITEM) .GE. 1) .AND. (IUINTV(IITEM) .LE. 6)) THEN
            IONEXT(IUINTV(IITEM), IITEM)
     &    = IONEXT(IUINTV(IITEM), IITEM) + IOINTV(IITEM)
            CALL CYH2SS(
     O               TONEXT(IITEM),
     I               IONEXT(1, IITEM))
            CALL CSS2YH(
     O               IONEXT(1, IITEM),
     I               TONEXT(IITEM))
         ELSE IF (IFLOUT(IITEM) .EQ. 1) THEN
            WRITE(JFPAR, *) '*** TMSTUP: NO SUCH UNIT OF TIME ***'
#ifdef OPT_PARALLEL
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
            STOP
#endif
         END IF
      END DO

      IF ((IURINT .GE. 1) .AND. (IURINT .LE. 6)) THEN
         IRNEXT(IURINT) = IRNEXT(IURINT) + IRINTV
         CALL CYH2SS(
     O               TRNEXT,
     I               IRNEXT)
         CALL CSS2YH(
     O               IRNEXT,
     I               TRNEXT)
      ELSE
         WRITE(JFPAR, *) '*** TMSTUP: NO SUCH UNIT OF TIME ***'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
      END IF

      WRITE(JFPAR, *) '***** Time Control Parameters *****'
      WRITE(JFPAR, *)
      WRITE(JFPAR, *) '     Start time  :', ITSTRT
      WRITE(JFPAR, *) '       End time  :', ITEND
      IF (IUTSTP .EQ. 0) IUTSTP = 4
      WRITE(JFPAR, *) '      Time step  :', TMSTP, ' ', CUNIT(IUTSTP)
      WRITE(JFPAR, *) '     Time Split  :', NTSPLT
      WRITE(JFPAR, *)
      WRITE(JFPAR, *) 'Restart interval :', IRINTV, ' ', CUNIT(IURINT)
      WRITE(JFPAR, *)
      WRITE(JFPAR, *) ' Euler-Backward :',IEB   ,' step'
      WRITE(JFPAR, *) ' Leap-frog      :',ILEAP ,' step'
      WRITE(JFPAR, *) ' 1 cycle        :',ITSCYC,' step'
      DO IITEM = 1, NOHITM
         IF (IFLOUT(IITEM) .EQ. 1) THEN
            WRITE(JFPAR, *)
            WRITE(JFPAR, *) 'Output item      :', CITEM(IITEM)
            WRITE(JFPAR, *) 'Output starts at :',
     &                       (IOSTRT(I, IITEM), I = 1, 6)
            WRITE(JFPAR, *) 'Output  ends  at :',
     &                       (IOEND(I, IITEM), I = 1, 6)
            WRITE(JFPAR, *) 'Output interval  :', IOINTV(IITEM),
     &                       ' ', CUNIT(IUINTV(IITEM))
         END IF
      END DO

      RETURN

c =====================================================================

      ENTRY TMSTPC (
     O                ITST,     TS,    ITS,
     O                NTSS,    TSS,
     G              OFLOUT, OFLSTK, ORSOUT, ORSRWD,
     I                  NT,     TT)

c     ***** Selection of the time integration scheme *****

      ITST = ITSTT( MOD (NT,ITSCYC) )
      IF ( ITST .GE. 3 ) THEN
         TS  = 2.0D+0 * DT
         ITS = 2
         NTSS = NTSPLT
         TSS  = TS / DBLE(NTSS)
      ELSE
         TS  = DT
         ITS = 1
         NTSS = NTSPLT / 2
         TSS  = TS / DBLE(NTSS)
      END IF

      DO IITEM = 1, NOHITM
         IF (IFLOUT(IITEM) .EQ. 1) THEN
            IF (      (TT .GT. TOSTRT(IITEM))
     &          .AND. (TT .LT. TOEND(IITEM)+DT)) THEN
               OFLSTK(IITEM) = .TRUE.
               IF (TT .GE. TONEXT(IITEM)) THEN
                  OFLOUT(IITEM) = .TRUE.
                  IONEXT(IUINTV(IITEM), IITEM)
     &          = IONEXT(IUINTV(IITEM), IITEM) + IOINTV(IITEM)
                  CALL CYH2SS(
     O                  TONEXT(IITEM),
     I                  IONEXT(1, IITEM))
                  CALL CSS2YH(
     O                  IONEXT(1, IITEM),
     I                  TONEXT(IITEM))
               ELSE
                  OFLOUT(IITEM) = .FALSE.
               END IF
            ELSE
               OFLSTK(IITEM) = .FALSE.
               OFLOUT(IITEM) = .FALSE.
            END IF
         ELSE
            OFLSTK(IITEM) = .FALSE.
            OFLOUT(IITEM) = .FALSE.
         END IF
      END DO

      IF (TT .GE. TRNEXT) THEN
         ORSOUT = .TRUE.
         IF (IRSRWD .EQ. 1) THEN
            ORSRWD = .TRUE.
         ELSE
            ORSRWD = .FALSE.
         END IF
         IRNEXT(IURINT) = IRNEXT(IURINT) + IRINTV
         CALL CYH2SS(
     O               TRNEXT,
     I               IRNEXT)
         CALL CSS2YH(
     O               IRNEXT,
     I               TRNEXT)
      ELSE
         ORSOUT = .FALSE.
      END IF

      RETURN
      END
