      SUBROUTINE TMINTP(
     O                   DITEM,
     I                   IITEM)
      use zocdim
      use zocgrd
      use zocfil
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Time interpolation of the boundary condition data
c
c  HISTORY
c     '02.10.10  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.07.08  H.Hasumi
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

      INTEGER   NITEM
      PARAMETER(NITEM = (NTDIM-2)*2+10)
      REAL*8   DITEM(NXDIM, NYDIM)
      INTEGER  IITEM

      INTEGER IDATE1(6, NITEM), IDATE2(6, NITEM)
      REAL*8   TIME1(NITEM),  TIME2(NITEM)
      REAL*8   DATA1(NX, NY, NITEM),   DATA2(NX, NY, NITEM)
      LOGICAL OSNGLD(NITEM), OFIRST(NITEM),     OF
      LOGICAL   OEOF
      INTEGER IYSKIQ(NITEM)
      SAVE IDATE1, IDATE2, TIME1, TIME2, DATA1, DATA2,
     &     OSNGLD, OFIRST, IYSKIQ

#ifdef OPT_PARALLEL
      REAL*8   DATAG(NXG, NYG)
#endif
      INTEGER IDATES(6), IDATET(6)
      INTEGER   IYTT, IYDIFF
      REAL*8   TIMES,  TIMET
      REAL*8   TINTV,  TINTP,  TINTQ
      INTEGER      I,      J
      CHARACTER CDATE*16, CHEAD(64)*16, CFITEM*(NCF)
      INTEGER  IFPAR,  JFPAR

      INTEGER IYSKIP(NITEM)
      NAMELIST /NMSKIP/ IYSKIP

      INTEGER NFITEM(NITEM)
      CHARACTER CFTAUX*(NCF), CFTAUY*(NCF)
      CHARACTER CFUSFC*(NCF), CFTSFC*(NCF)
      CHARACTER CFQSFC*(NCF), CFWFLX*(NCF)
      CHARACTER CFSWNT*(NCF), CFDWLW*(NCF)
      CHARACTER CFPSFC*(NCF), CFSSFC*(NCF)
      CHARACTER CFTREF(NTDIM)*(NCF), CFTDMP(NTDIM)*(NCF)
      SAVE NFITEM
      NAMELIST /NMSFBC/ CFTAUX, CFTAUY, CFUSFC, CFTSFC, CFQSFC,
     &                  CFWFLX, CFSWNT, CFDWLW, CFPSFC, CFSSFC,
     &                  CFTREF, CFTDMP

      DATA CFTAUX, CFTAUY / 'not-specified', 'not-specified' /
      DATA CFUSFC, CFTSFC / 'not-specified', 'not-specified' /
      DATA CFQSFC, CFWFLX / 'not-specified', 'not-specified' /
      DATA CFSWNT, CFDWLW / 'not-specified', 'not-specified' /
      DATA CFPSFC, CFSSFC / 'not-specified', 'not-specified' /
      DATA CFTREF / NTDIM*'not-specified' /
      DATA CFTDMP / NTDIM*'not-specified' /

      DATA OSNGLD / NITEM*.FALSE. /
      DATA OFIRST / NITEM*.TRUE. /
      DATA OF / .TRUE. /
      DATA IYSKIP / NITEM*1 /

      IF (OF) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSFBC, END=132)
 132     CONTINUE
         OF = .FALSE.
      END IF

      IF (IITEM .EQ. 1) THEN
         CFITEM = CFTAUX
      ELSE IF (IITEM .EQ. 2) THEN
         CFITEM = CFTAUY
      ELSE IF (IITEM .EQ. 3) THEN
         CFITEM = CFUSFC
      ELSE IF (IITEM .EQ. 4) THEN
         CFITEM = CFTSFC
      ELSE IF (IITEM .EQ. 5) THEN
         CFITEM = CFQSFC
      ELSE IF (IITEM .EQ. 6) THEN
         CFITEM = CFWFLX
      ELSE IF (IITEM .EQ. 7) THEN
         CFITEM = CFSWNT
      ELSE IF (IITEM .EQ. 8) THEN
         CFITEM = CFDWLW
      ELSE IF (IITEM .EQ. 9) THEN
         CFITEM = CFPSFC
      ELSE IF (IITEM .EQ. 10) THEN
         CFITEM = CFSSFC
      ELSE IF (IITEM. GT. NITEM) THEN
         WRITE(JFPAR, *) '*** TMINTP: NO SUCH ITEM ***'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
      ELSE IF (MOD(IITEM, 2) .EQ. 1) THEN
         I = (IITEM - 10) / 2 + 1
         CFITEM = CFTREF(I)
      ELSE
         I = (IITEM - 10) / 2
         CFITEM = CFTDMP(I)
      END IF

      IF (OFIRST(IITEM)) THEN
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            CALL FILOPN(NFITEM(IITEM), CFITEM, 'READ')
         END IF
#else
         CALL FILOPN(NFITEM(IITEM), CFITEM, 'READ')
#endif
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSKIP, END=123)
  123    CONTINUE
         IF (IYSKIP(IITEM) .LT. 1) THEN
            IYSKIQ(IITEM) = 1
         ELSE
            IYSKIQ(IITEM) = IYSKIP(IITEM)
         END IF
         CALL CSS2YH(
     O               IDATET,
     I                   TT)
         IYTT = IDATET(1)

#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            REWIND(NFITEM(IITEM))
            READ(NFITEM(IITEM)) CHEAD
            READ(NFITEM(IITEM))
     &      ((DATAG(I, J), I = 1, NXG), J = 1, NYG)
         END IF
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_SFC(DATA1(1, 1, IITEM), DATAG)
#else
         REWIND(NFITEM(IITEM))
         READ(NFITEM(IITEM)) CHEAD
         READ(NFITEM(IITEM))
     &   ((DATA1(I, J, IITEM), I = 1, NX), J = 1, NY)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE1(I, IITEM), I = 1, 6)
         IDATET(1) = 0
         IDATES(1) = 0
         DO 10 I = 2, 6
            IDATES(I) = IDATE1(I, IITEM)
   10    CONTINUE
         CALL CYH2SS(
     O                TIMET,
     I               IDATET)
         CALL CYH2SS(
     O                TIMES,
     I               IDATES)
         IF (TIMET .GE. TIMES) THEN
            IDATES(1) = IYTT
         ELSE
            IDATES(1) = IYTT-1
         END IF
         CALL CYH2SS(
     O                TIME1(IITEM),
     I               IDATES)

#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            OEOF = .FALSE.
            READ(NFITEM(IITEM), END=198) CHEAD
            READ(NFITEM(IITEM))
     &      ((DATAG(I, J), I = 1, NXG), J = 1, NYG)
            GO TO 298
 198        OEOF = .TRUE.
 298        CONTINUE
         END IF
         CALL MPI_BCAST
     &        (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
         IF (OEOF) GO TO 98
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_SFC(DATA2(1, 1, IITEM), DATAG)
#else
         READ(NFITEM(IITEM), END=98) CHEAD
         READ(NFITEM(IITEM))
     &   ((DATA2(I, J, IITEM), I = 1, NX), J = 1, NY)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
         IYDIFF = IDATE2(1, IITEM) - IDATE1(1, IITEM)
         IDATET(1) = IDATES(1) + IYDIFF
         DO 20 I = 2, 6
            IDATET(I) = IDATE2(I, IITEM)
   20    CONTINUE
         CALL CYH2SS(
     O                TIME2(IITEM),
     I               IDATET)
         IF (      (TT .GE. TIME1(IITEM))
     &       .AND. (TT .LE. TIME2(IITEM))) GO TO 99

   90    CONTINUE
         TIME1(IITEM) = TIME2(IITEM)
         DO 30 I = 1, 6
            IDATE1(I, IITEM) = IDATE2(I, IITEM)
   30    CONTINUE
         DO 40 J = 1, NY
            DO 40 I = 1, NX
               DATA1(I, J, IITEM) = DATA2(I, J, IITEM)
   40    CONTINUE
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            OEOF = .FALSE.
            READ(NFITEM(IITEM), END=197) CHEAD
            READ(NFITEM(IITEM))
     &      ((DATAG(I, J), I = 1, NXG), J = 1, NYG)
            GO TO 297
 197        OEOF = .TRUE.
 297        CONTINUE
         END IF
         CALL MPI_BCAST
     &        (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
         IF (OEOF) GO TO 97
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_SFC(DATA2(1, 1, IITEM), DATAG)
#else
         READ(NFITEM(IITEM), END=97) CHEAD
         READ(NFITEM(IITEM))
     &   ((DATA2(I, J, IITEM), I = 1, NX), J = 1, NY)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
         IYDIFF = IDATE2(1, IITEM) - IDATE1(1, IITEM)
         IDATET(1) = IDATES(1) + IYDIFF
         DO 50 I = 2, 6
            IDATET(I) = IDATE2(I, IITEM)
   50    CONTINUE
         CALL CYH2SS(
     O                   TIME2(IITEM),
     I                  IDATET)
         IF (      (TT .GE. TIME1(IITEM))
     &       .AND. (TT .LE. TIME2(IITEM))) GO TO 99
         GO TO 90

 97      CONTINUE
         TIME1(IITEM) = TIME2(IITEM)
         DO 60 I = 1, 6
            IDATE1(I, IITEM) = IDATE2(I, IITEM)
   60    CONTINUE
         DO 70 J = 1, NY
            DO 70 I = 1, NX
               DATA1(I, J, IITEM) = DATA2(I, J, IITEM)
   70    CONTINUE
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            REWIND(NFITEM(IITEM))
            READ(NFITEM(IITEM)) CHEAD
            READ(NFITEM(IITEM))
     &      ((DATAG(I, J), I = 1, NXG), J = 1, NYG)
         END IF
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_SFC(DATA2(1, 1, IITEM), DATAG)
#else
         REWIND(NFITEM(IITEM))
         READ(NFITEM(IITEM)) CHEAD
         READ(NFITEM(IITEM))
     &   ((DATA2(I, J, IITEM), I = 1, NX), J = 1, NY)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
         CALL CSS2YH(
     O               IDATES,
     I                TIME1(IITEM))
         IDATET(1) = IDATES(1)
         DO 80 I = 2, 6
            IDATET(I) = IDATE2(I, IITEM)
   80    CONTINUE
         CALL CYH2SS(
     O                TIMES,
     I               IDATES)
         CALL CYH2SS(
     O                TIMET,
     I               IDATET)
         IF (TIMET .LE. TIMES) THEN
            IYDIFF = IYSKIQ(IITEM)
         ELSE
            IYDIFF = IYSKIQ(IITEM) - 1
         END IF
         IDATET(1) = IDATES(1) + IYDIFF
         CALL CYH2SS(
     O                TIME2(IITEM),
     I                IDATET)
         IF (      (TT .GE. TIME1(IITEM))
     &       .AND. (TT .LE. TIME2(IITEM))) GO TO 99

   96    WRITE(JFPAR, *)
     &   '*** TMINTP: UNEXPECTED ERROR ***'
         WRITE(JFPAR, *)
     &   '-> Please inform the developer of the situation.'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
   98    OSNGLD(IITEM) = .TRUE.
   99    CONTINUE
         OFIRST(IITEM) = .FALSE.
      END IF

      IF (OSNGLD(IITEM)) THEN
         DO 100 J = 1, NY
            DO 100 I = 1, NX
               DITEM(ISTR+I-1, JSTR+J-1) = DATA1(I, J, IITEM)
  100    CONTINUE
      ELSE
         IF (TT .GT. TIME2(IITEM)) THEN
            TIME1(IITEM) = TIME2(IITEM)
            DO 300 I = 1, 6
               IDATE1(I, IITEM) = IDATE2(I, IITEM)
  300       CONTINUE
            DO 310 J = 1, NY
               DO 310 I = 1, NX
                  DATA1(I, J, IITEM) = DATA2(I, J, IITEM)
  310       CONTINUE
#ifdef OPT_PARALLEL
            IF (MYRANK .EQ. IROOT) THEN
               OEOF = .FALSE.
               READ(NFITEM(IITEM), END=897) CHEAD
               READ(NFITEM(IITEM))
     &         ((DATAG(I, J), I = 1, NXG), J = 1, NYG)
               GO TO 797
 897           OEOF = .TRUE.
 797           CONTINUE
            END IF
            CALL MPI_BCAST
     &           (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
            IF (OEOF) GO TO 997
            CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                     IROOT, MPI_COMM_WORLD, IERR)
            CALL SCATTER_SFC(DATA2(1, 1, IITEM), DATAG)
#else
            READ(NFITEM(IITEM), END=997) CHEAD
            READ(NFITEM(IITEM))
     &      ((DATA2(I, J, IITEM), I = 1, NX), J = 1, NY)
#endif
            CDATE = CHEAD(50)
            READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
            IYDIFF = IDATE2(1, IITEM) - IDATE1(1, IITEM)
            CALL CSS2YH(
     O                  IDATES,
     I                   TIME1(IITEM))
            IDATET(1) = IDATES(1) + IYDIFF
            DO 320 I = 2, 6
               IDATET(I) = IDATE2(I, IITEM)
  320       CONTINUE
            CALL CYH2SS(
     O                   TIME2(IITEM),
     I                  IDATET)
            GO TO 999

 997        CONTINUE
#ifdef OPT_PARALLEL
            IF (MYRANK .EQ. IROOT) THEN
               OEOF = .FALSE.
               REWIND(NFITEM(IITEM))
               READ(NFITEM(IITEM), END=697) CHEAD
               READ(NFITEM(IITEM))
     &         ((DATAG(I, J), I = 1, NXG), J = 1, NYG)
               GO TO 597
 697           OEOF = .TRUE.
 597           CONTINUE
            END IF
            CALL MPI_BCAST
     &           (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
            IF (OEOF) GO TO 997
            CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                     IROOT, MPI_COMM_WORLD, IERR)
            CALL SCATTER_SFC(DATA2(1, 1, IITEM), DATAG)
#else
            REWIND(NFITEM(IITEM))
            READ(NFITEM(IITEM), END=997) CHEAD
            READ(NFITEM(IITEM))
     &      ((DATA2(I, J, IITEM), I = 1, NX), J = 1, NY)
#endif
            CDATE = CHEAD(50)
            READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
            CALL CSS2YH(
     O                  IDATES,
     I                   TIME1(IITEM))
            IDATET(1) = IDATES(1)
            DO 400 I = 2, 6
               IDATET(I) = IDATE2(I, IITEM)
  400       CONTINUE
            CALL CYH2SS(
     O                   TIMES,
     I                  IDATES)
            CALL CYH2SS(
     O                   TIMET,
     I                  IDATET)
            IF (TIMET .LE. TIMES) THEN
               IYDIFF = IYSKIQ(IITEM)
            ELSE
               IYDIFF = IYSKIQ(IITEM) - 1
            END IF
            IDATET(1) = IDATES(1) + IYDIFF
            CALL CYH2SS(
     O                   TIME2(IITEM),
     I                  IDATET)
            
  999       CONTINUE
         END IF
         TINTV = TIME2(IITEM) - TIME1(IITEM)
         TINTP = (TT - TIME1(IITEM)) / TINTV
         TINTQ = (TIME2(IITEM) - TT) / TINTV
         DO 200 J = 1, NY
            DO 200 I = 1, NX
               DITEM(ISTR+I-1, JSTR+J-1)
     &       = TINTQ * DATA1(I, J, IITEM)
     &       + TINTP * DATA2(I, J, IITEM)
  200    CONTINUE
      END IF

      RETURN
      END
#ifdef OPT_BODY
c *********************************************************************

      SUBROUTINE TMINTB(
     O                   DITEM,
     I                   IITEM)
      use zocdim
      use zocgrd
      use zocfil
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Time interpolation of the boundary condition data
c
c  HISTORY
c     '99.10.04  H.Hasumi
c     '00.12.07  H.Hasumi: combine par and non-par routines
c     '01.05.10  H.Hasumi
c     '01.12.07  H.Hasumi
c     '02.06.02  H.Hasumi: tracer dimension
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

      INTEGER   NITEM
c      PARAMETER(NITEM = NTDIM*2)
      PARAMETER(NITEM = 4)
      REAL*8   DITEM(NXDIM, NYDIM, NZDIM)
      INTEGER  IITEM

      INTEGER IDATE1(6, NITEM), IDATE2(6, NITEM)
      REAL*8   TIME1(NITEM),  TIME2(NITEM)
      REAL*8   DATA1(NX, NY, NZ, NITEM),   DATA2(NX, NY, NZ, NITEM)
      LOGICAL OSNGLD(NITEM), OFIRST(NITEM),     OF
      LOGICAL   OEOF
      INTEGER IYSKIQ(NITEM)
      SAVE IDATE1, IDATE2, TIME1, TIME2, DATA1, DATA2,
     &     OSNGLD, OFIRST, IYSKIQ

#ifdef OPT_PARALLEL
      REAL*8   DATAG(NXG, NYG, NZ)
#endif
      INTEGER IDATES(6), IDATET(6)
      INTEGER   IYTT, IYDIFF
      REAL*8   TIMES,  TIMET
      REAL*8   TINTV,  TINTP,  TINTQ
      INTEGER      I,      J,      K
      CHARACTER CDATE*16, CHEAD(64)*16, CFITEM*(NCF)
      INTEGER  IFPAR,  JFPAR

      INTEGER IYSKIB(NITEM)
      NAMELIST /NMSKIB/ IYSKIB

      INTEGER NFITEM(NITEM)
      CHARACTER CFTBDY(NTDIM)*(NCF), CFTDMB(NTDIM)*(NCF)
      SAVE NFITEM
      NAMELIST /NMBODY/ CFTBDY, CFTDMB

      DATA CFTBDY / NTDIM*'not-specified' /
      DATA CFTDMB / NTDIM*'not-specified' /

      DATA OSNGLD / NITEM*.FALSE. /
      DATA OFIRST / NITEM*.TRUE. /
      DATA OF / .TRUE. /
      DATA IYSKIB / NITEM*1 /

      IF (OF) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBODY, END=132)
 132     CONTINUE
         OF = .FALSE.
      END IF

c      IF (IITEM .LE. NTDIM) THEN
c         CFITEM = CFTBDY(IITEM)
c      ELSE IF (IITEM .LE. 2*NTDIM) THEN
c         CFITEM = CFTDMB(IITEM-NTDIM)
      IF (IITEM .EQ. 1) THEN
         CFITEM = CFTBDY(1)
      ELSE IF (IITEM .EQ. 2) THEN
         CFITEM = CFTBDY(2)
      ELSE IF (IITEM .EQ. 3) THEN
         CFITEM = CFTDMB(1)
      ELSE IF (IITEM .EQ. 4) THEN
         CFITEM = CFTDMB(2)
      ELSE
         WRITE(JFPAR, *) '*** TMINTB: NO SUCH ITEM ***'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
      END IF

      IF (OFIRST(IITEM)) THEN
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            CALL FILOPN(NFITEM(IITEM), CFITEM, 'READ')
         END IF
#else
         CALL FILOPN(NFITEM(IITEM), CFITEM, 'READ')
#endif
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSKIB, END=123)
  123    CONTINUE
         IF (IYSKIB(IITEM) .LT. 1) THEN
            IYSKIQ(IITEM) = 1
         ELSE
            IYSKIQ(IITEM) = IYSKIB(IITEM)
         END IF
         CALL CSS2YH(
     O               IDATET,
     I                   TT)
         IYTT = IDATET(1)

#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            REWIND(NFITEM(IITEM))
            READ(NFITEM(IITEM)) CHEAD
            READ(NFITEM(IITEM)) (((DATAG(I, J, K), I = 1, NXG),
     &                                             J = 1, NYG),
     &                                             K = 1, NZ)
         END IF
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_BDY(DATA1(1, 1, 1, IITEM), DATAG)
#else
         REWIND(NFITEM(IITEM))
         READ(NFITEM(IITEM)) CHEAD
         READ(NFITEM(IITEM)) (((DATA1(I, J, K, IITEM), I = 1, NX),
     &                                                 J = 1, NY),
     &                                                 K = 1, NZ)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE1(I, IITEM), I = 1, 6)
         IDATET(1) = 0
         IDATES(1) = 0
         DO 10 I = 2, 6
            IDATES(I) = IDATE1(I, IITEM)
   10    CONTINUE
         CALL CYH2SS(
     O                TIMET,
     I               IDATET)
         CALL CYH2SS(
     O                TIMES,
     I               IDATES)
         IF (TIMET .GE. TIMES) THEN
            IDATES(1) = IYTT
         ELSE
            IDATES(1) = IYTT-1
         END IF
         CALL CYH2SS(
     O                TIME1(IITEM),
     I               IDATES)

#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            OEOF = .FALSE.
            READ(NFITEM(IITEM), END=198) CHEAD
            READ(NFITEM(IITEM)) (((DATAG(I, J, K), I = 1, NXG),
     &                                             J = 1, NYG),
     &                                             K = 1, NZ)
            GO TO 298
 198        OEOF = .TRUE.
 298        CONTINUE
         END IF
         CALL MPI_BCAST
     &        (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
         IF (OEOF) GO TO 98
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_BDY(DATA2(1, 1, 1, IITEM), DATAG)
#else
         READ(NFITEM(IITEM), END=98) CHEAD
         READ(NFITEM(IITEM)) (((DATA2(I, J, K, IITEM), I = 1, NX),
     &                                                 J = 1, NY),
     &                                                 K = 1, NZ)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
         IYDIFF = IDATE2(1, IITEM) - IDATE1(1, IITEM)
         IDATET(1) = IDATES(1) + IYDIFF
         DO 20 I = 2, 6
            IDATET(I) = IDATE2(I, IITEM)
   20    CONTINUE
         CALL CYH2SS(
     O                TIME2(IITEM),
     I               IDATET)
         IF (      (TT .GE. TIME1(IITEM))
     &       .AND. (TT .LE. TIME2(IITEM))) GO TO 99

   90    CONTINUE
         TIME1(IITEM) = TIME2(IITEM)
         DO 30 I = 1, 6
            IDATE1(I, IITEM) = IDATE2(I, IITEM)
   30    CONTINUE
         DO 40 K = 1, NZ
            DO 40 J = 1, NY
               DO 40 I = 1, NX
                  DATA1(I, J, K, IITEM) = DATA2(I, J, K, IITEM)
   40    CONTINUE
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            OEOF = .FALSE.
            READ(NFITEM(IITEM), END=197) CHEAD
            READ(NFITEM(IITEM)) (((DATAG(I, J, K), I = 1, NXG),
     &                                             J = 1, NYG),
     &                                             K = 1, NZ)
            GO TO 297
 197        OEOF = .TRUE.
 297        CONTINUE
         END IF
         CALL MPI_BCAST
     &        (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
         IF (OEOF) GO TO 97
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_BDY(DATA2(1, 1, 1, IITEM), DATAG)
#else
         READ(NFITEM(IITEM), END=97) CHEAD
         READ(NFITEM(IITEM)) (((DATA2(I, J, K, IITEM), I = 1, NX),
     &                                                 J = 1, NY),
     &                                                 K = 1, NZ)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
         IYDIFF = IDATE2(1, IITEM) - IDATE1(1, IITEM)
         IDATET(1) = IDATES(1) + IYDIFF
         DO 50 I = 2, 6
            IDATET(I) = IDATE2(I, IITEM)
   50    CONTINUE
         CALL CYH2SS(
     O                   TIME2(IITEM),
     I                  IDATET)
         IF (      (TT .GE. TIME1(IITEM))
     &       .AND. (TT .LE. TIME2(IITEM))) GO TO 99
         GO TO 90

 97      CONTINUE
         TIME1(IITEM) = TIME2(IITEM)
         DO 60 I = 1, 6
            IDATE1(I, IITEM) = IDATE2(I, IITEM)
   60    CONTINUE
         DO 70 K = 1, NZ
            DO 70 J = 1, NY
               DO 70 I = 1, NX
                  DATA1(I, J, K, IITEM) = DATA2(I, J, K, IITEM)
   70    CONTINUE
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            REWIND(NFITEM(IITEM))
            READ(NFITEM(IITEM)) CHEAD
            READ(NFITEM(IITEM)) (((DATAG(I, J, K), I = 1, NXG),
     &                                             J = 1, NYG),
     &                                             K = 1, NZ)
         END IF
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL SCATTER_BDY(DATA2(1, 1, 1, IITEM), DATAG)
#else
         REWIND(NFITEM(IITEM))
         READ(NFITEM(IITEM)) CHEAD
         READ(NFITEM(IITEM)) (((DATA2(I, J, K, IITEM), I = 1, NX),
     &                                                 J = 1, NY),
     &                                                 K = 1, NZ)
#endif
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
         CALL CSS2YH(
     O               IDATES,
     I                TIME1(IITEM))
         IDATET(1) = IDATES(1)
         DO 80 I = 2, 6
            IDATET(I) = IDATE2(I, IITEM)
   80    CONTINUE
         CALL CYH2SS(
     O                TIMES,
     I               IDATES)
         CALL CYH2SS(
     O                TIMET,
     I               IDATET)
         IF (TIMET .LE. TIMES) THEN
            IYDIFF = IYSKIQ(IITEM)
         ELSE
            IYDIFF = IYSKIQ(IITEM) - 1
         END IF
         IDATET(1) = IDATES(1) + IYDIFF
         CALL CYH2SS(
     O                TIME2(IITEM),
     I                IDATET)
         IF (      (TT .GE. TIME1(IITEM))
     &       .AND. (TT .LE. TIME2(IITEM))) GO TO 99

   96    WRITE(JFPAR, *)
     &   '*** TMINTP: UNEXPECTED ERROR ***'
         WRITE(JFPAR, *)
     &   '-> Please inform the developer of the situation.'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
   98    OSNGLD(IITEM) = .TRUE.
   99    CONTINUE
         OFIRST(IITEM) = .FALSE.
      END IF

      IF (OSNGLD(IITEM)) THEN
         DO 100 K = 1, NZ
            DO 100 J = 1, NY
               DO 100 I = 1, NX
                  DITEM(ISTR+I-1, JSTR+J-1, KSTR+K-1)
     &          = DATA1(I, J, K, IITEM)
  100    CONTINUE
      ELSE
         IF (TT .GT. TIME2(IITEM)) THEN
            TIME1(IITEM) = TIME2(IITEM)
            DO 300 I = 1, 6
               IDATE1(I, IITEM) = IDATE2(I, IITEM)
  300       CONTINUE
            DO 310 K = 1, NZ
               DO 310 J = 1, NY
                  DO 310 I = 1, NX
                     DATA1(I, J, K, IITEM) = DATA2(I, J, K, IITEM)
  310       CONTINUE
#ifdef OPT_PARALLEL
            IF (MYRANK .EQ. IROOT) THEN
               OEOF = .FALSE.
               READ(NFITEM(IITEM), END=897) CHEAD
               READ(NFITEM(IITEM)) (((DATAG(I, J, K), I = 1, NXG),
     &                                                J = 1, NYG),
     &                                                K = 1, NZ)
               GO TO 797
 897           OEOF = .TRUE.
 797           CONTINUE
            END IF
            CALL MPI_BCAST
     &           (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
            IF (OEOF) GO TO 997
            CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                     IROOT, MPI_COMM_WORLD, IERR)
            CALL SCATTER_BDY(DATA2(1, 1, 1, IITEM), DATAG)
#else
            READ(NFITEM(IITEM), END=997) CHEAD
            READ(NFITEM(IITEM)) (((DATA2(I, J, K, IITEM), I = 1, NX),
     &                                                    J = 1, NY),
     &                                                    K = 1, NZ)
#endif
            CDATE = CHEAD(50)
            READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
            IYDIFF = IDATE2(1, IITEM) - IDATE1(1, IITEM)
            CALL CSS2YH(
     O                  IDATES,
     I                   TIME1(IITEM))
            IDATET(1) = IDATES(1) + IYDIFF
            DO 320 I = 2, 6
               IDATET(I) = IDATE2(I, IITEM)
  320       CONTINUE
            CALL CYH2SS(
     O                   TIME2(IITEM),
     I                  IDATET)
            GO TO 999

 997        CONTINUE
#ifdef OPT_PARALLEL
            IF (MYRANK .EQ. IROOT) THEN
               OEOF = .FALSE.
               REWIND(NFITEM(IITEM))
               READ(NFITEM(IITEM), END=697) CHEAD
               READ(NFITEM(IITEM)) (((DATAG(I, J, K), I = 1, NXG),
     &                                                J = 1, NYG),
     &                                                K = 1, NZ)
               GO TO 597
 697           OEOF = .TRUE.
 597           CONTINUE
            END IF
            CALL MPI_BCAST
     &           (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
            IF (OEOF) GO TO 997
            CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,
     &                     IROOT, MPI_COMM_WORLD, IERR)
            CALL SCATTER_BDY(DATA2(1, 1, 1, IITEM), DATAG)
#else
            REWIND(NFITEM(IITEM))
            READ(NFITEM(IITEM), END=997) CHEAD
            READ(NFITEM(IITEM)) (((DATA2(I, J, K, IITEM), I = 1, NX),
     &                                                    J = 1, NY),
     &                                                    K = 1, NZ)
#endif
            CDATE = CHEAD(50)
            READ(CDATE, '(I6.6,5I2.2)') (IDATE2(I, IITEM), I = 1, 6)
            CALL CSS2YH(
     O                  IDATES,
     I                   TIME1(IITEM))
            IDATET(1) = IDATES(1)
            DO 400 I = 2, 6
               IDATET(I) = IDATE2(I, IITEM)
  400       CONTINUE
            CALL CYH2SS(
     O                   TIMES,
     I                  IDATES)
            CALL CYH2SS(
     O                   TIMET,
     I                  IDATET)
            IF (TIMET .LE. TIMES) THEN
               IYDIFF = IYSKIQ(IITEM)
            ELSE
               IYDIFF = IYSKIQ(IITEM) - 1
            END IF
            IDATET(1) = IDATES(1) + IYDIFF
            CALL CYH2SS(
     O                   TIME2(IITEM),
     I                  IDATET)
            
  999       CONTINUE
         END IF
         TINTV = TIME2(IITEM) - TIME1(IITEM)
         TINTP = (TT - TIME1(IITEM)) / TINTV
         TINTQ = (TIME2(IITEM) - TT) / TINTV
         DO 200 K = 1, NZ
            DO 200 J = 1, NY
               DO 200 I = 1, NX
                  DITEM(ISTR+I-1, JSTR+J-1, KSTR+K-1)
     &          = TINTQ * DATA1(I, J, K, IITEM)
     &          + TINTP * DATA2(I, J, K, IITEM)
  200    CONTINUE
      END IF

      RETURN
      END
#endif
