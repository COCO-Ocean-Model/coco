      SUBROUTINE CHKSET
      use zocdim
      use zocgrd
      use zocfil
      use zocout
      use ufile
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Output the data
c
c  HISTORY
c     '99.04.14  H.Hasumi: from CCSR2
c     '01.05.09  H.Hasumi
c     '01.12.07  H.Hasumi
c     '02.05.29  H.Hasumi: combine parallel and nonparallel
c     '03.06.03  H.Hasumi: consistency with the change in atmct.F
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '08.??.??  Y.Komuro: SNGOUG/DBLOUG are allocated in all the nodes                           
c
c ---------------------------------------------------------------------

      LOGICAL OFLOUT(NFOMAX), OFLSTK(NFOMAX)
      INTEGER NXITEM, NYITEM, NZITEM, NDITEM
      CHARACTER CCITEM*(*), CLAS*(*)
      REAL*8   DITEM(NDITEM)

      INTEGER  NXITM(NFOMAX),  NYITM(NFOMAX),  NZITM(NFOMAX)
      INTEGER NFUNIT(NFOMAX)
      INTEGER ITOPAD(NFOMAX), KTOPAD(NFOMAX), NSZITM(NFOMAX)
      REAL*8  FCTAVR(NFOMAX), CNTAVR(NFOMAX),  TTOLD(NFOMAX)
      LOGICAL OSINGL(NFOMAX),   OADD(NFOMAX), OFIRST(NFOMAX)
      INTEGER IFLOUT(NFOMAX), IOAVRG(NFOMAX), ISINGL(NFOMAX)
      INTEGER  IXSTR(NFOMAX),  IXEND(NFOMAX)
      INTEGER  JYSTR(NFOMAX),  JYEND(NFOMAX)
      INTEGER  KZSTR(NFOMAX),  KZEND(NFOMAX)
      INTEGER NOHITM
      CHARACTER CITEM(NFOMAX)*16, CFITEM(NFOMAX)*(NCF)
      SAVE NXITM, NYITM, NZITM,
     &     NFUNIT, ITOPAD, KTOPAD, NSZITM,
     &     FCTAVR, CNTAVR, TTOLD, OSINGL, OADD, OFIRST,
     &     IFLOUT, IOAVRG, ISINGL,
     &     IXSTR, IXEND, JYSTR, JYEND, KZSTR, KZEND,
     &     CITEM, CFITEM, NOHITM

      INTEGER   ITEM
      INTEGER  IITEM, JTOPAD
      INTEGER      I,      J,      K
      INTEGER    IJK,   IJKM
      INTEGER  IDATE(6)
      INTEGER IOHITM
      REAL*8   TSTRT, TOSTRT,   TOUT
      CHARACTER  CHEAD(64)*16, CTITEM*16
      INTEGER  IFPAR, JFPAR
      SAVE IFPAR, JFPAR

      INTEGER ITSTRT(6), ITEND(6)
      INTEGER IUTSTP
      REAL*8  TMSTP
      INTEGER NTSPLT
      INTEGER IOSTRT(6, NFOMAX),  IOEND(6, NFOMAX)
      INTEGER IOINTV(NFOMAX), IUINTV(NFOMAX)
      INTEGER IODSTR(6), IODEND(6), IODINT, IUDINT, IODAVR, IODSNG
      INTEGER IOHSTR(6), IOHEND(6), IOHINT, IUHINT, IOHAVR, IOHSNG
      INTEGER IOXSTR, IOXEND, IOYSTR, IOYEND, IOZSTR, IOZEND
      INTEGER  IXDIM,  JYDIM,  KZDIM
      CHARACTER COHFIL*(NCF), COHITM*16

      NAMELIST /NMTIME/ ITSTRT, ITEND, TMSTP, IUTSTP, NTSPLT
      NAMELIST /NMDOUT/ IODSTR, IODEND, IODINT, IUDINT, IODAVR, IODSNG
      NAMELIST /NMHIST/ COHITM, COHFIL,
     &                  IOHSTR, IOHEND, IOHINT, IUHINT, IOHAVR, IOHSNG,
     &                  IOXSTR, IOXEND, IOYSTR, IOYEND, IOZSTR, IOZEND

      DATA ITSTRT / 0, 0, 0, 0, 0, 0 /
      DATA ITEND  / 0, 0, 0, 0, 0, 0 /
      DATA IUTSTP / 4 /
      DATA TMSTP  / 1.D10 /
      DATA NTSPLT / 0 /
      DATA IODSTR / 0, 0, 0, 0, 0, 0 /
      DATA IODEND / 0, 0, 0, 0, 0, 0 /
      DATA IODINT, IUDINT / 1, 1 /
      DATA IODAVR, IODSNG / 0, 1 /

      DATA CHEAD  / 64*'                ' /
      DATA CITEM  / NFOMAX*'                ' /
      DATA JTOPAD / 1 /
      DATA CNTAVR / NFOMAX*0.D0 /
      DATA IFLOUT / NFOMAX*0 /
      DATA IOAVRG / NFOMAX*0 /
      DATA ISINGL / NFOMAX*1 /
      DATA OFIRST / NFOMAX*.TRUE. /

      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMTIME, END=99)
 99   CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMDOUT, END=98)
 98   CONTINUE
      DO IITEM = 1, NFOMAX
         DO I = 1, 6
            IOSTRT(I, IITEM) = IODSTR(I)
            IOEND (I, IITEM) = IODEND(I)
         END DO
         IOINTV(IITEM) = IODINT
         IUINTV(IITEM) = IUDINT
         IOAVRG(IITEM) = IODAVR
         ISINGL(IITEM) = IODSNG
         IXSTR(IITEM) = 1
         JYSTR(IITEM) = 1
#ifdef OPT_PARALLEL
         IXEND(IITEM) = NXG
         JYEND(IITEM) = NYG
#else
         IXEND(IITEM) = NX
         JYEND(IITEM) = NY
#endif
         KZSTR(IITEM) = 1
         KZEND(IITEM) = NZ
      END DO
            
      CALL REWNML(IFPAR, JFPAR)
      NOHITM = 0
 90   CONTINUE
      COHITM = 'not-specified'
      COHFIL = 'not-specified'
      IOHSTR(1) = -1
      IOHEND(1) = -1
      IOHINT = -1
      IUHINT = -1
      IOHAVR = -1
      IOHSNG = -1
      IOXSTR = -1
      IOXEND = -1
      IOYSTR = -1
      IOYEND = -1
      IOZSTR = -1
      IOZEND = -1
      READ(IFPAR, NMHIST, END=97)
      IF (COHITM(1:13) .NE. 'not-specified') THEN
         NOHITM = NOHITM + 1
         IOHITM = NOHITM
         IFLOUT(IOHITM) = 1
         CITEM(IOHITM) = COHITM
         CFITEM(IOHITM) = COHFIL
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            CALL FILOPN(NFUNIT(IOHITM), COHFIL, 'WRITE')
         END IF
#else
         CALL FILOPN(NFUNIT(IOHITM), COHFIL, 'WRITE')
#endif
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
         IF (IOHAVR .GE. 0) THEN
            IOAVRG(IOHITM) = IOHAVR
         END IF
         IF (IOHSNG .GE. 0) THEN
            ISINGL(IOHITM) = IOHSNG
         END IF
#ifdef OPT_PARALLEL
         IF ((IOXSTR .GT. 0) .AND. (IOXSTR .LE. NXG)) THEN
#else
         IF ((IOXSTR .GT. 0) .AND. (IOXSTR .LE. NX)) THEN
#endif
            IXSTR(IOHITM) = IOXSTR
         END IF
#ifdef OPT_PARALLEL
         IF ((IOXEND .GE. IXSTR(IOHITM)) .AND. (IOXEND .LE. NXG)) THEN
#else
         IF ((IOXEND .GE. IXSTR(IOHITM)) .AND. (IOXEND .LE. NX)) THEN
#endif
            IXEND(IOHITM) = IOXEND
         END IF
#ifdef OPT_PARALLEL
         IF ((IOYSTR .GT. 0) .AND. (IOYSTR .LE. NYG)) THEN
#else
         IF ((IOYSTR .GT. 0) .AND. (IOYSTR .LE. NY)) THEN
#endif
            JYSTR(IOHITM) = IOYSTR
         END IF
#ifdef OPT_PARALLEL
         IF ((IOYEND .GE. JYSTR(IOHITM)) .AND. (IOYEND .LE. NYG)) THEN
#else
         IF ((IOYEND .GE. JYSTR(IOHITM)) .AND. (IOYEND .LE. NY)) THEN
#endif
            JYEND(IOHITM) = IOYEND
         END IF
         IF ((IOZSTR .GT. 0) .AND. (IOZSTR .LE. NZ)) THEN
            KZSTR(IOHITM) = IOZSTR
         END IF
         IF ((IOZEND .GE. KZSTR(IOHITM)) .AND. (IOZEND .LE. NZ)) THEN
            KZEND(IOHITM) = IOZEND
         END IF
      END IF
      GO TO 90
 97   CONTINUE

#ifdef OPT_PARALLEL
      ALLOCATE(DBLOUG(1, 1, 1))
      ALLOCATE(SNGOUG(1, 1, 1))
      IF (MYRANK .EQ. IROOT) THEN
         DEALLOCATE(DBLOUG, SNGOUG)
         ALLOCATE(DBLOUG(NXG, NYG, NZ))
         ALLOCATE(SNGOUG(NXG, NYG, NZ))
      END IF
#endif
      DO 10 I = 1, NWORK
         WRKOUT(I) = 0.D0
   10 CONTINUE
      DO 30 IITEM = 1, NOHITM
         IF (IOAVRG(IITEM) .EQ. 1) THEN
            FCTAVR(IITEM) = 1.D0
         ELSE
            FCTAVR(IITEM) = 0.D0
         END IF
         IF (ISINGL(IITEM) .EQ. 1) THEN
            OSINGL(IITEM) = .TRUE.
         ELSE
            OSINGL(IITEM) = .FALSE.
         END IF
   30 CONTINUE
      CALL CYH2SS(
     O             TSTRT,
     I            ITSTRT)
      DO IITEM = 1, NOHITM
         CALL CYH2SS(
     O            TOSTRT,
     I            IOSTRT(1, IITEM))
         TTOLD(IITEM) = MAX(TOSTRT,  TSTRT)
      END DO

      RETURN

c =====================================================================

      ENTRY CHKSTK(
     I             OFLSTK)

      DO IITEM = 1, NOHITM
         IF (OFLSTK(IITEM)) THEN
            CNTAVR(IITEM) = CNTAVR(IITEM) * FCTAVR(IITEM) + 1.D0
            OADD(IITEM) = .TRUE.
         ELSE
            OADD(IITEM) = .FALSE.
         END IF
      END DO

      RETURN

c =====================================================================

      ENTRY CHKOUT(
     I             OFLOUT)

      DO IITEM = 1, NOHITM
         IF (OFLOUT(IITEM)) THEN
            IF (IFLOUT(IITEM) .EQ. 1) THEN

               IF (IOAVRG(IITEM) .EQ. 1) THEN
                  TOUT = (TT + TTOLD(IITEM)) * 0.5D0
               ELSE
                  TOUT = TT
               END IF
               CALL CSS2YH(
     O                       IDATE,
     I                        TOUT)
               IXDIM = IXEND(IITEM) - IXSTR(IITEM) + 1
               JYDIM = JYEND(IITEM) - JYSTR(IITEM) + 1
               KZDIM = KZEND(IITEM) - KZSTR(IITEM) + 1
               CHEAD(3) = CITEM(IITEM)
               WRITE(CHEAD(27), '(I4.4,2I2.2,A1,3I2.2,A1)')
     &              IDATE(1), IDATE(2), IDATE(3), ' ',
     &              IDATE(4), IDATE(5), IDATE(6), ' '
               WRITE(CHEAD(50), '(I6.6,5I2.2)') IDATE
               WRITE(CHEAD(29), '(I15,A1)') NXITM(IITEM), 'X'
               WRITE(CHEAD(30), '(I16)') IXSTR(IITEM)
               WRITE(CHEAD(31), '(I16)') IXEND(IITEM)
               WRITE(CHEAD(32), '(I15,A1)') NYITM(IITEM), 'Y'
               WRITE(CHEAD(33), '(I16)') JYSTR(IITEM)
               WRITE(CHEAD(34), '(I16)') JYEND(IITEM)
               WRITE(CHEAD(35), '(I15,A1)') NZITM(IITEM), 'Z'
               WRITE(CHEAD(36), '(I16)') KZSTR(IITEM)
               WRITE(CHEAD(37), '(I16)') KZEND(IITEM)
               WRITE(CHEAD(64), '(I16)') IXDIM*JYDIM*KZDIM

#ifdef OPT_PARALLEL
               IF (MYRANK .LT. IJNODE) THEN
                  DO K = KZSTR(IITEM), KZEND(IITEM)
                     DO J = 1, NY
                        DO I = 1, NX
                           IJK = (K - KZSTR(IITEM)) * NXY
     &                         + (J - 1) * NX + I
                           IJKM = (K - KZSTR(IITEM)) * NXYDIM
     &                          + (JSTR + J - 2) * NXDIM
     &                          + ISTR + I - 1
                           DBLEOU(IJK)
     &                   = WRKOUT(ITOPAD(IITEM) + IJKM - 1)
     &                     / CNTAVR(IITEM)
                        END DO
                     END DO
                  END DO
               END IF
               CALL GATHER_CHK(DBLOUG, DBLEOU)
               IF (MYRANK .EQ. IROOT) THEN
                  IF (OSINGL(IITEM)) THEN
                     DO K = 1, KZDIM
                        DO J = 1, NYG
                           DO I = 1, NXG
                              SNGOUG(I, J, K) = DBLOUG(I, J, K)
                           END DO
                        END DO
                     END DO
                     CHEAD(38) = 'REAL4'
                     WRITE(NFUNIT(IITEM)) CHEAD
                     WRITE(NFUNIT(IITEM))
     &                    (((SNGOUG(I, J, K),
     &                       I = IXSTR(IITEM), IXEND(IITEM)),
     &                       J = JYSTR(IITEM), JYEND(IITEM)),
     &                       K = 1, KZDIM)
                  ELSE
                     CHEAD(38) = 'REAL8'
                     WRITE(NFUNIT(IITEM)) CHEAD
                     WRITE(NFUNIT(IITEM))
     &                    (((DBLOUG(I, J, K),
     &                       I = IXSTR(IITEM), IXEND(IITEM)),
     &                       J = JYSTR(IITEM), JYEND(IITEM)),
     &                       K = 1, KZDIM)
                  END IF
               END IF

#else
               IF (OSINGL(IITEM)) THEN
                  DO K = KZSTR(IITEM), KZEND(IITEM)
                     DO J = JYSTR(IITEM), JYEND(IITEM)
                        DO I = IXSTR(IITEM), IXEND(IITEM)
                           IJK = (K - KZSTR(IITEM)) * IXDIM * JYDIM
     &                         + (J - JYSTR(IITEM)) * IXDIM
     &                         + I - IXSTR(IITEM) + 1
                           IJKM = (K - KZSTR(IITEM)) * NXYDIM
     &                          + (JSTR + J - 2) * NXDIM
     &                          + ISTR + I - 1
                           SNGLOU(IJK)
     &                   = WRKOUT(ITOPAD(IITEM) + IJKM - 1)
     &                     / CNTAVR(IITEM)
                        END DO
                     END DO
                  END DO
                  CHEAD(38) = 'REAL4'
                  WRITE(NFUNIT(IITEM)) CHEAD
                  WRITE(NFUNIT(IITEM)) (SNGLOU(IJK),
     &                                  IJK = 1, IXDIM*JYDIM*KZDIM)
               ELSE
                  DO K = KZSTR(IITEM), KZEND(IITEM)
                     DO J = JYSTR(IITEM), JYEND(IITEM)
                        DO I = IXSTR(IITEM), IXEND(IITEM)
                           IJK = (K - KZSTR(IITEM)) * IXDIM * JYDIM
     &                         + (J - JYSTR(IITEM)) * IXDIM
     &                         + I - IXSTR(IITEM) + 1
                           IJKM = (K - KZSTR(IITEM)) * NXYDIM
     &                          + (JSTR + J - 2) * NXDIM
     &                          + ISTR + I - 1
                           DBLEOU(IJK)
     &                   = WRKOUT(ITOPAD(IITEM) + IJKM - 1)
     &                     / CNTAVR(IITEM)
                        END DO
                     END DO
                  END DO
                  CHEAD(38) = 'REAL8'
                  WRITE(NFUNIT(IITEM)) CHEAD
                  WRITE(NFUNIT(IITEM)) (DBLEOU(IJK),
     &                                  IJK = 1, IXDIM*JYDIM*KZDIM)
               END IF
#endif

#ifdef OPT_PARALLEL
               IF (MYRANK .LT. IJNODE) THEN
#endif
               DO 310 I = 1, NSZITM(IITEM)
                  WRKOUT(ITOPAD(IITEM) + I - 1) = 0.D0
 310           CONTINUE
#ifdef OPT_PARALLEL
               END IF
#endif

            END IF
            CNTAVR(IITEM) = 0.D0
            TTOLD(IITEM) = TT
#ifdef OPT_PARALLEL
            IF (MYRANK .EQ. IROOT) THEN
#endif
            WRITE(JFPAR, *) '*** File output ***'
            WRITE(JFPAR, *) ' Item :', CITEM(IITEM),
     &                       '  Time :', IDATE, '  Step :', NT
#ifdef OPT_PARALLEL
            END IF
#endif
         END IF
      END DO

      RETURN

c =====================================================================

      ENTRY CHEKIN(
     I              DITEM, CCITEM,
     I             NXITEM, NYITEM, NZITEM, NDITEM,   CLAS)

      DO IOHITM = 1, NOHITM
         CTITEM = CCITEM
         IF (CITEM(IOHITM) .EQ. CCITEM) THEN
            ITEM = IOHITM
            GO TO 123
         END IF
      END DO
      GO TO 456
 123  CONTINUE

      IF (OFIRST(ITEM)) THEN
         CALL REWNML(IFPAR, JFPAR)
         IF (CLAS(1:3) .EQ. 'SFC') THEN
            KZSTR(ITEM) = 1
            KZEND(ITEM) = 1
         ELSE IF (CLAS(1:3) .EQ. 'ICE') THEN
            KZSTR(ITEM) = 1
            KZEND(ITEM) = NIC
         END IF
#ifdef OPT_PARALLEL
         NXITM(ITEM) = NXG
         NYITM(ITEM) = NYG
#else
         NXITM(ITEM) = NXITEM
         NYITM(ITEM) = NYITEM
#endif
         NZITM(ITEM) = NZITEM
         NSZITM(ITEM) = NXYDIM * (KZEND(ITEM) - KZSTR(ITEM) + 1)

         IF (IFLOUT(ITEM) .EQ. 1) THEN
            ITOPAD(ITEM) = JTOPAD
            JTOPAD = JTOPAD + NSZITM(ITEM)
            IF (JTOPAD .GT. NWORK+1) THEN
               WRITE(JFPAR, *) '### WORK AREA SHORTAGE ###'
               WRITE(JFPAR, *) '### THE ITEM NUMBERED', ITEM,
     &                          'WILL NOT BE OUTPUT ###'
               ITOPAD(ITEM) = 0
               IFLOUT(ITEM) = 0
               JTOPAD = JTOPAD - NSZITM(ITEM)
            END IF
            IF (CLAS(1:3) .EQ. 'SFC') THEN
               KTOPAD(ITEM) = 1
            ELSE IF (CLAS(1:3) .EQ. 'ICE') THEN
               KTOPAD(ITEM) = NXYDIM + 1
            ELSE IF (CLAS(1:3) .EQ. 'OCN') THEN
               KTOPAD(ITEM) = (KSTR + KZSTR(ITEM) - 2) * NXYDIM + 1
            ELSE
               WRITE(NFSTDO, *) '### NO SUCH OUTPUT CLASS AS',
     &                          CLAS, '###'
               KTOPAD(ITEM) = 1
            END IF
         ELSE
            ITOPAD(ITEM) = 0
         END IF

         OFIRST(ITEM) = .FALSE.
      END IF

      IF (OADD(ITEM)) THEN
         IF (IFLOUT(ITEM) .EQ. 1) THEN
            DO 100 I = 1, NSZITM(ITEM)
               WRKOUT(ITOPAD(ITEM) + I - 1) =
     &         WRKOUT(ITOPAD(ITEM) + I - 1) * FCTAVR(ITEM)
     &       + DITEM(KTOPAD(ITEM) + I - 1)
  100       CONTINUE
         END IF
      END IF

 456  CONTINUE

      RETURN
      END
