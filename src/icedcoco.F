      PROGRAM ICEDCOCO

      use zocdim
      use zocfil
      use ucloc
      use atmct
      use brstt
      use ufile
      use bfrch
      implicit none

c --- information -----------------------------------------------------
c
c  COCO3.4 coupled with sea ice
c
c  HISTORY
c     '02.10.10  H.Hasumi: from COCO3.4
c     '06.01.30  H.Hasumi: change the timing of calling MPI_INIT
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '07.09.27  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.09.02  Y.Komuro: STDOUT filename format changed
c     '10.04.14  M.kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

      REAL*8     TA(NXDIM, NYDIM, NZDIM, NTDIM)
      REAL*8     TB(NXDIM, NYDIM, NZDIM, NTDIM)
      REAL*8     UA(NXDIM, NYDIM, NZDIM),     UB(NXDIM, NYDIM, NZDIM)
      REAL*8     VA(NXDIM, NYDIM, NZDIM),     VB(NXDIM, NYDIM, NZDIM)
      REAL*8     HA(NXDIM, NYDIM),            HB(NXDIM, NYDIM)
      REAL*8   UBTA(NXDIM, NYDIM),          UBTB(NXDIM, NYDIM)
      REAL*8   VBTA(NXDIM, NYDIM),          VBTB(NXDIM, NYDIM)
      REAL*8      W(NXDIM, NYDIM, NZDIM),      R(NXDIM, NYDIM, NZDIM)
      REAL*8    AMV(NXDIM, NYDIM, NZDIM),    AHV(NXDIM, NYDIM, NZDIM)

      REAL*8     AA(NXDIM, NYDIM, 0:NIC),     AB(NXDIM, NYDIM, 0:NIC)
      REAL*8    HIA(NXDIM, NYDIM, 0:NIC),    HIB(NXDIM, NYDIM, 0:NIC)
      REAL*8    UIA(NXDIM, NYDIM),    UIB(NXDIM, NYDIM)
      REAL*8    VIA(NXDIM, NYDIM),    VIB(NXDIM, NYDIM)
      REAL*8    TIA(NXDIM, NYDIM, 0:NIC),    TIB(NXDIM, NYDIM, 0:NIC)
      REAL*8    HSA(NXDIM, NYDIM, 0:NIC),    HSB(NXDIM, NYDIM, 0:NIC)
      REAL*8    TSI(NXDIM, NYDIM, 0:NIC)

      REAL*8     FT(NXDIM, NYDIM, NTDIM)
      REAL*8  SWABS(NXDIM, NYDIM),     FS(NXDIM, NYDIM)
      REAL*8   TAUX(NXDIM, NYDIM),   TAUY(NXDIM, NYDIM)
      REAL*8   PTOP(NXDIM, NYDIM)

c      REAL*8  WRK(NWRK)
c      COMMON /WORK/ WRK

      REAL*8      DT,     TT,     TS,    TSS
      INTEGER     NT,   ITST,    ITS,   NTSS
      REAL*8   TSTRT,   TEND
      LOGICAL OFLOUT(NFOMAX), OFLSTK(NFOMAX), ORSOUT, ORSRWD
      INTEGER    IJK
      INTEGER  IFPAR,  JFPAR
#ifdef OPT_PARALLEL
      INTEGER LENSTD
#endif

      CHARACTER CRUN*(NCF), CSTDO*(NCF)

      NAMELIST /NMRUN/ CRUN
      NAMELIST /NMSTDO/ CSTDO

      DATA CRUN / '(RUN NAME WAS NOT SET)' /
      DATA CSTDO / 'STDOUT' /

c *** Initial setup ***

#ifdef OPT_PARALLEL
      CALL MPI_INIT(IERR)
#endif
      CALL CLCSTR('SETUP')
      CALL REWNML(IFPAR, JFPAR)
      OPEN(UNIT=IFPAR, FILE='PARAMET', STATUS='OLD',
     &     ACCESS='SEQUENTIAL', FORM='FORMATTED')
#ifdef OPT_PARALLEL
      CALL PARSET
#endif
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMSTDO, END=98)
 98   CONTINUE
#ifdef OPT_PARALLEL
      LENSTD = INDEX(CSTDO, ' ')
      WRITE(CSTDO(LENSTD:LENSTD+3), '(A1,I3.3)') '.', MYRANK
#endif
      CALL REWNML(IFPAR, JFPAR)
      OPEN(UNIT=JFPAR, FILE=CSTDO,
     &     ACCESS='SEQUENTIAL', FORM='FORMATTED')
#ifdef OPT_PARALLEL
      WRITE(NFSTDO, *) 'MESSAGE OUTPUT FOR RANK', MYRANK
#endif
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMRUN, END=99)
   99 CONTINUE
      WRITE(NFSTDO, *) 'Run name :'//CRUN
c      DO IJK = 1, NWRK
c         WRK(IJK) = 0.D0
c      END DO

      OFINAL = .FALSE.
      OINIT = .FALSE.
      CALL TMSTUP(
     O             TSTRT,   TEND,     DT)
      CALL RESTRT(
     I             TSTRT,
     O                UB,     VB,     TB,     
     O                HB,   UBTB,   VBTB,
     O                 W,    AMV,    AHV,
     O                AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     O               TSI,
     O                FT,  SWABS,     FS,
     O              TAUX,   TAUY,   PTOP)
      CALL OCSTUP(
     O                UB,     VB,     TB,  
     O                HB,   UBTB,   VBTB,      W,      R,
     O               AMV,    AHV,
     O                FT,   PTOP,
     I                DT)
      OINIT = .TRUE.
      CALL OCEAN (
     M                UA,     VA,     TA,
     M                HA,   UBTA,   VBTA,
     M                AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     M                UB,     VB,     TB,
     M                HB,   UBTB,   VBTB,
     M                AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     M                 W,      R,    AMV,    AHV,
     M                FT,  SWABS,     FS,
     M              TAUX,   TAUY,   PTOP,    TSI,
     I                NT,     TT,   ITST,     TS,    ITS,
     I              NTSS,    TSS,
     I            OFLOUT, OFLSTK)
      OINIT = .FALSE.
      CALL FORSTO(
     O                UA,     VA,     TA,
     O                HA,   UBTA,   VBTA,
     I                UB,     VB,     TB,
     I                HB,   UBTB,   VBTB)
      CALL FORSTI(
     O                AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                AB,    HIB,    UIB,    VIB,    TIB,    HSB)
      CALL CLCEND('SETUP')

c *** Main loop ***

      CALL CLCSTR('MAIN')
c      DO IJK = 1, NWRK
c         WRK(IJK) = 0.D0
c      END DO
      TT = TSTRT
      NT = 0
  100 CONTINUE
         TT = TT + DT
         NT = NT + 1
         CALL TMSTPC(
     O                 ITST,     TS,    ITS,
     O                 NTSS,    TSS,
     O               OFLOUT, OFLSTK, ORSOUT, ORSRWD,
     I                   NT,     TT)
         CALL OCEAN (
     M                   UA,     VA,     TA,
     M                   HA,   UBTA,   VBTA,
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     M                   UB,     VB,     TB,
     M                   HB,   UBTB,   VBTB,
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     M                    W,      R,    AMV,    AHV,
     M                   FT,  SWABS,     FS,
     M                 TAUX,   TAUY,   PTOP,    TSI,
     I                   NT,     TT,   ITST,     TS,    ITS,
     I                 NTSS,    TSS,
     I               OFLOUT, OFLSTK)
         IF (TT .GE. TEND) GO TO 999
         IF (ITST .EQ. 3) THEN
            CALL FINOUT(
     I                      TT,     NT,
     I                      UA,     VA,     TA,
     I                      HA,   UBTA,   VBTA,
     I                       W,    AMV,    AHV,
     I                      AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                     TSI,
     I                      FT,  SWABS,     FS,
     I                    TAUX,   TAUY,   PTOP,
     I                  ORSOUT, ORSRWD)
         ELSE
            CALL FINOUT(
     I                      TT,     NT,
     I                      UB,     VB,     TB,
     I                      HB,   UBTB,   VBTB,
     I                       W,    AMV,    AHV,
     I                      AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     I                     TSI,
     I                      FT,  SWABS,     FS,
     I                    TAUX,   TAUY,   PTOP,
     I                  ORSOUT, ORSRWD)
         END IF
         IF (ORSOUT) THEN
            OFINAL = .TRUE.
            CALL OCEAN (
     M                   UA,     VA,     TA,
     M                   HA,   UBTA,   VBTA,
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     M                   UB,     VB,     TB,
     M                   HB,   UBTB,   VBTB,
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     M                    W,      R,    AMV,    AHV,
     M                   FT,  SWABS,     FS,
     M                 TAUX,   TAUY,   PTOP,    TSI,
     I                   NT,     TT,   ITST,     TS,    ITS,
     I                 NTSS,    TSS,
     I               OFLOUT, OFLSTK)
            OFINAL = .FALSE.
         END IF
      GO TO 100
  999 CONTINUE
      CALL CLCEND('MAIN')

c *** Handling for termination ***

      CALL CLCSTR('FINOUT')
      OFINAL = .TRUE.
      ORSOUT = .TRUE.
      IF (ITST .EQ. 3) THEN
         CALL FINOUT(
     I                      TT,     NT,
     I                      UA,     VA,     TA,
     I                      HA,   UBTA,   VBTA,
     I                       W,    AMV,    AHV,
     I                      AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                     TSI,
     I                      FT,  SWABS,     FS,
     I                    TAUX,   TAUY,   PTOP,
     I                  ORSOUT, ORSRWD)
      ELSE
         CALL FINOUT(
     I                      TT,     NT,
     I                      UB,     VB,     TB,
     I                      HB,   UBTB,   VBTB,
     I                       W,    AMV,    AHV,
     I                      AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     I                     TSI,
     I                      FT,  SWABS,     FS,
     I                    TAUX,   TAUY,   PTOP,
     I                  ORSOUT, ORSRWD)
      END IF
      CALL OCEAN (
     M                   UA,     VA,     TA,
     M                   HA,   UBTA,   VBTA,
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     M                   UB,     VB,     TB,
     M                   HB,   UBTB,   VBTB,
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     M                    W,      R,    AMV,    AHV,
     M                   FT,  SWABS,     FS,
     M                 TAUX,   TAUY,   PTOP,    TSI,
     I                   NT,     TT,   ITST,     TS,    ITS,
     I                 NTSS,    TSS,
     I               OFLOUT, OFLSTK)
      CALL CLCEND('FINOUT')
      CALL CLCOUT
#ifdef OPT_PARALLEL
      CALL PARFIN
#endif

      END
#ifdef OPT_PARALLEL
c *********************************************************************

      SUBROUTINE PARSET
      use zocdim
      use ufile

c --- information -----------------------------------------------------
c
c  Initialization for inter-node communication using MPI
c
c  HISTORY
c     '00.05.30  H.Hasumi: from COCO2
c
c ---------------------------------------------------------------------

#include "mpif.h"

      INTEGER  IFPAR,  JFPAR
#ifdef OPT_TRIPOLE
      INTEGER  II
#endif
      INTEGER NDROOT
      NAMELIST /NMROOT/ NDROOT
      DATA NDROOT / 0 /

ccc      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMROOT, END=99)
   99 CONTINUE

      IJNODE = INODES * JNODES
      IROOT  = NDROOT

      IF (IJNODE .GT. NPROCS) THEN
         WRITE(JFPAR, *)
     &   '### NODE NUMBER ERROR                   ###'
         WRITE(JFPAR, *)
     &   '### INODES*JNODES must be = or < NPROCS ###'
         WRITE(JFPAR, *) 'INODES:', INODES
         WRITE(JFPAR, *) 'JNODES:', JNODES
         WRITE(JFPAR, *) 'NPROCS:', NPROCS
         CALL MPI_FINALIZE(IERR)
         STOP
      END IF

      IF ((NDROOT .LT. 0) .OR. (NDROOT .GT. NPROCS-1)) THEN
         WRITE(JFPAR, *)
     &   '### ERROR in the choice of NDROOT          ###'
         WRITE(JFPAR, *) 'NDROOT:', NDROOT
         WRITE(JFPAR, *) 'NPROCS:', NPROCS
         WRITE(JFPAR, *)
     &   '### NDROOT must lie between 0 and NPROCS-1 ###'
         CALL MPI_FINALIZE(IERR)
         STOP
      END IF

      IRANK = MOD(MYRANK, INODES)
      JRANK = MYRANK / INODES

      IF (IRANK .EQ. 0) THEN
         IDOWN = INODES * (JRANK + 1) - 1
      ELSE
         IDOWN = INODES * JRANK + IRANK - 1
      END IF
      IF (IRANK .EQ. INODES-1) THEN
         IUP = INODES * JRANK
      ELSE
         IUP = INODES * JRANK + IRANK + 1
      END IF

      IF (JRANK .EQ. 0) THEN
         JDOWN = MPI_PROC_NULL
      ELSE
         JDOWN = INODES * (JRANK - 1) + IRANK
      END IF

#ifdef OPT_TRIPOLE 
      IF (JRANK .EQ. JNODES-1) THEN

      IF (INODES.GE.2) THEN

         II = MOD(INODES,2)
         IF(II.NE.0) THEN
            WRITE(JFPAR, *) 'INODES:', INODES
            WRITE(JFPAR, *)
     &      '### INODES must be 1 or multiples of 2 ###'
            CALL MPI_FINALIZE(IERR)
            STOP
         ENDIF
         
         II = INODES / 2
         IF(IRANK.LT.II) THEN
            JUPE = INODES * JNODES - 1 - IRANK
            JUPW = MPI_PROC_NULL
         ELSE
            JUPE = MPI_PROC_NULL
            JUPW = INODES * JNODES  -1 - IRANK  
         ENDIF

      ELSE

         JUPE = MYRANK
         JUPW = MYRANK

      ENDIF

      ELSE

         JUPE = MPI_PROC_NULL
         JUPW = MPI_PROC_NULL

      ENDIF
#endif
      IF (JRANK .EQ. JNODES-1) THEN
         JUP = MPI_PROC_NULL
      ELSE
         JUP = INODES * (JRANK + 1) + IRANK
      END IF

      IF (MYRANK .GE. INODES*JNODES) THEN
         IUP = MPI_PROC_NULL
         IDOWN = MPI_PROC_NULL
         JUP = MPI_PROC_NULL
         JDOWN = MPI_PROC_NULL
      END IF

      CALL GS2DST
      CALL GS3DST
      CALL GSIDST

      RETURN
      END

c *********************************************************************

      SUBROUTINE PARFIN
      use zocnod
c --- information -----------------------------------------------------
c
c  Finalization for parallelization
c
c  HISTORY
c     '00.05.30  H.Hasumi: from COCO3
c
c ---------------------------------------------------------------------

#include "mpif.h"

      CALL MPI_FINALIZE(IERR)

      RETURN
      END
#endif
