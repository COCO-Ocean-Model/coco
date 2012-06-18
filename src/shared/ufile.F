      SUBROUTINE FLOPCL
      use zocfil
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Open and close files.
c
c  HISTORY
c     '01.12.06  H.Hasumi
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

      CHARACTER     CF*(NCF),   CACT*(*)
      INTEGER       NF
      LOGICAL      OEX
      INTEGER        I,   IERR

      LOGICAL     OPN(NFMAX)
      SAVE OPN

      DATA OPN / NFMAX*.FALSE. /

c =====================================================================

      ENTRY FILOPN(
     O                 NF,
     I                 CF,   CACT)

      IF (CACT .EQ. 'READ') THEN
         INQUIRE(FILE=CF, EXIST=OEX)
         IF (.NOT. OEX) THEN
            I = INDEX(CF, ' ') - 1
            WRITE(NFSTDO, *) '### FILE "', CF(1:I),
     &                       '" DOES NOT EXIST ###'
#ifdef OPT_PARALLEL
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
            STOP
#endif
         END IF
      ELSE IF (CF(1:13) .eq. 'not-specified') THEN
         WRITE(NFSTDO, *) '### NAME NOT SPECIFIED FOR OUTPUT FILE ###'
#ifdef OPT_PARALLEL
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
         STOP
#endif
      END IF

      DO I = NFMIN, NFMAX
         IF (.NOT. OPN(I)) THEN
            NF = I
            OPN(NF) = .TRUE.
            GO TO 123
         END IF
      END DO
      WRITE(NFSTDO, *) '### MAXIMUM FILE NUMBER EXCEEDED ###'
#ifdef OPT_PARALLEL
      CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
#else
      STOP
#endif
 123  CONTINUE
      OPEN(UNIT=NF, FILE=CF,
     &     FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
      REWIND(UNIT=NF)
      I = INDEX(CF, ' ') - 1
      WRITE(NFSTDO, *) '*** FILE "', CF(1:I),
     &                 '" OPENED FOR UNIT', NF, '***'

      RETURN

c =====================================================================

      ENTRY FILCLS(
     I                 NF)

      CLOSE(UNIT=NF)
      OPN(NF) = .FALSE.
      WRITE(NFSTDO, *) '*** FILE UNIT', NF, 'CLOSED ***'

      RETURN
      END

c *********************************************************************

      SUBROUTINE REWNML(
     O                   IFILE,  JFILE)
      use zocfil

      INTEGER  IFILE,  JFILE

      IFILE = NFPARM
      JFILE = NFSTDO
      REWIND(IFILE, ERR=999)
      RETURN

 999  WRITE(JFILE, *) '*** ERROR IN REWINDING THE NAMELIST FILE ***'
      RETURN
      END

c *********************************************************************

      subroutine cstnml(
     I                    jfile,  cnsbr,  cnnml,  istat)
      use zocfil

      integer  jfile, istat 
      character*(*)  cnsbr,  cnnml

      if (istat > 0) then ! an error occured while reading
         write(jfile, *)
     &      '*** ERROR OCCURS while reading namelist:', cnnml,
     &      ' at subroutine:', cnsbr, ' ***'
         stop
      end if

      return
      end subroutine cstnml
