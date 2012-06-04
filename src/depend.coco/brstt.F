      SUBROUTINE RESTRT(
     I               TSTRT,
     O                  UB,     VB,     TB,    
     O                  HB,   UBTB,   VBTB,
     O                   W,    AMV,    AHV,
     O                  FT)
      use zocdim
      use zocfil
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Reading and writing restart file
c
c  HISTORY
c     '99.08.17  H.Hasumi: from CCSR2-MASK
c     '01.05.10  H.Hasumi
c     '01.12.07  H.Hasumi
c     '02.05.29  H.Nakano: tracer dimension
c     '02.06.02  H.Hasumi: modification to the above
c     '07.12.20  H.Hasumi
c     '09.09.07  Y.Komuro: save allocatable variables (G3D etc.)
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

      REAL*8      UB(NXDIM, NYDIM, NZDIM),     VB(NXDIM, NYDIM, NZDIM)
      REAL*8      TB(NXDIM, NYDIM, NZDIM, NTDIM)
      REAL*8      HB(NXDIM, NYDIM)
      REAL*8    UBTB(NXDIM, NYDIM)       ,   VBTB(NXDIM, NYDIM)
      REAL*8       W(NXDIM, NYDIM, NZDIM)
      REAL*8     AMV(NXDIM, NYDIM, NZDIM),    AHV(NXDIM, NYDIM, NZDIM)
      REAL*8      FT(NXDIM, NYDIM, NTDIM)

      INTEGER   IXDIM,  JYDIM,  KZDIM
      LOGICAL    OEOF
      REAL*8     ADDITM(IXDIM, JYDIM, KZDIM)
      CHARACTER  CCITEM*(*), CLAS*(*)

#ifdef OPT_PARALLEL
      REAL(8),ALLOCATABLE :: BUF2(:,:), BUF3(:,:,:)
      REAL(8),ALLOCATABLE :: G2D(:,:), G3D(:,:,:)
      SAVE BUF2, BUF3, G2D, G3D
#else
      REAL*8    BUF2(NX, NY),   BUF3(NX, NY, NZ)
c      COMMON /WORK/ BUF2, BUF3
#endif
#ifdef OPT_TRIPOLE
      REAL*8    FACT  ! set 1 for scalar, set -1 for vector
      INTEGER   IOFF, JOFF !(IOFF,JOFF)=(0,0)for T-cell, (-1,-1)for V-cell 
#endif
      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8   TSTRT,     TT
      INTEGER NTSTEP
      LOGICAL ORSOUT, ORSRWD

      REAL*8    TTT
      INTEGER IDATE(6)
      INTEGER      I,      J,      K,      L
      CHARACTER  CHEAD(64)*16,  CDATE*16
      INTEGER  IFPAR,  JFPAR
      SAVE IFPAR, JFPAR

      INTEGER   NFINIT, IRSTRT
      CHARACTER CFINIT*(NCF)
      INTEGER   NFREST
      CHARACTER CFREST*(NCF)
      SAVE NFINIT, NFREST
      NAMELIST /NMFINI/ CFINIT, IRSTRT
      NAMELIST /NMFRST/ CFREST

      DATA OFIRST / .TRUE. /
      DATA CFINIT, CFREST / 'not-specified', 'not-specified' /
      DATA IRSTRT / 0 /

      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMFINI, END=99)
   99 CONTINUE

#ifdef OPT_PARALLEL
      IF (MYRANK .EQ. IROOT) THEN
         CALL FILOPN(NFINIT, CFINIT, 'READ')
         WRITE(JFPAR, *) '*** Read initial condition file ***'
      END IF
#else
      CALL FILOPN(NFINIT, CFINIT, 'READ')
      WRITE(JFPAR, *) '*** Read initial condition file ***'
#endif

#ifdef OPT_PARALLEL
      TB = 0.D0
      FT = 0.D0

      IF(MYRANK .NE. IROOT) THEN
         ALLOCATE(BUF2(1,1))
         ALLOCATE(BUF3(1,1,1))
         ALLOCATE(G2D(1,1))
         ALLOCATE(G3D(1,1,1))
      ELSE
         ALLOCATE(BUF2(NXG,NYG))
         ALLOCATE(BUF3(NXG,NYG,NZ))
         ALLOCATE(G2D(NXGDIM,NYGDIM))
         ALLOCATE(G3D(NXGDIM,NYGDIM,NZDIM))
      ENDIF

      IF (MYRANK .EQ. IROOT) THEN
         REWIND(NFINIT)
         READ(NFINIT, END=109) CHEAD
         READ(NFINIT) BUF3
         DO 100 K = 1, NZ
            DO 100 J = 1, NYG
               DO 100 I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
 100     CONTINUE
 109     CONTINUE
      END IF
      CALL SCATTER_3D(UB, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=119) CHEAD
         READ(NFINIT) BUF3
         DO 110 K = 1, NZ
            DO 110 J = 1, NYG
               DO 110 I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
 110     CONTINUE
 119     CONTINUE
      END IF
      CALL SCATTER_3D(VB, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=129) CHEAD
         READ(NFINIT) BUF3
         DO 120 K = 1, NZ
            DO 120 J = 1, NYG
               DO 120 I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
 120     CONTINUE
 129     CONTINUE
      END IF
      CALL SCATTER_3D(TB, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=139) CHEAD
         READ(NFINIT) BUF3
         DO 130 K = 1, NZ
            DO 130 J = 1, NYG
               DO 130 I = 1, NXG
c                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) =
     &            ABS(BUF3(I, J, K))
 130     CONTINUE
 139     CONTINUE
      END IF
      CALL SCATTER_3D(TB(1, 1, 1, 2), G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=141) CHEAD
         READ(NFINIT) BUF2
         DO 140 J = 1, NYG
            DO 140 I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
 140     CONTINUE
 141     CONTINUE
      END IF
      CALL SCATTER_2D(HB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=143) CHEAD
         READ(NFINIT) BUF2
         DO 142 J = 1, NYG
            DO 142 I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
 142     CONTINUE
 143     CONTINUE
      END IF
      CALL SCATTER_2D(UBTB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=145) CHEAD
         READ(NFINIT) BUF2
         DO 144 J = 1, NYG
            DO 144 I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
 144     CONTINUE
 145     CONTINUE
      END IF
      CALL SCATTER_2D(VBTB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=159) CHEAD
         READ(NFINIT) BUF3
         DO 150 K = 1, NZ
            DO 150 J = 1, NYG
               DO 150 I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
 150     CONTINUE
 159     CONTINUE
      END IF
      CALL SCATTER_3D(W, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=169) CHEAD
         READ(NFINIT) BUF3
         DO 160 K = 1, NZ
            DO 160 J = 1, NYG
               DO 160 I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
 160     CONTINUE
 169     CONTINUE
      END IF
      CALL SCATTER_3D(AMV, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=179) CHEAD
         READ(NFINIT) BUF3
         DO 170 K = 1, NZ
            DO 170 J = 1, NYG
               DO 170 I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
 170     CONTINUE
 179     CONTINUE
      END IF
      CALL SCATTER_3D(AHV, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=189) CHEAD
         READ(NFINIT) BUF2
         DO 180 J = 1, NYG
            DO 180 I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
 180     CONTINUE
 189     CONTINUE
      END IF
      CALL SCATTER_2D(FT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=199) CHEAD
         READ(NFINIT) BUF2
         DO 190 J = 1, NYG
            DO 190 I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
 190     CONTINUE
 199     CONTINUE
      END IF
      CALL SCATTER_2D(FT(1, 1, 2), G2D)

      DO L = 3, NTDIM
         IF (MYRANK .EQ. IROOT) THEN
            READ(NFINIT, END=209) CHEAD
            READ(NFINIT) BUF3
            DO K = 1, NZ
               DO J = 1, NYG
                  DO I = 1, NXG
                     G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
                  END DO
               END DO
            END DO
 209        CONTINUE
         END IF
         CALL SCATTER_3D(TB(1, 1, 1, L), G3D)
      END DO

      DO L = 3, NTDIM
         IF (MYRANK .EQ. IROOT) THEN
            READ(NFINIT, END=219) CHEAD
            READ(NFINIT) BUF2
            DO J = 1, NYG
               DO I = 1, NXG
                  G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
               END DO
            END DO
 219        CONTINUE
         END IF
         CALL SCATTER_2D(FT(1, 1, L), G2D)
      END DO

#ifdef OPT_TRIPOLE
      CALL SHIFT2(    UB,      VB,
     D             NXDIM,   NYDIM,  NZDIM,
     I             -1.D0,      -1,     -1 )
      CALL SHIFT1(   TB,
     D             NXDIM,   NYDIM, NZTDIM,
     I              1.D0,       0,      0 )
      CALL SHIFT1(   AMV,
     D             NXDIM,   NYDIM,  NZDIM,
     I              1.D0,      -1,     -1 )
      CALL SHIFT2(     W,     AHV,
     D             NXDIM,   NYDIM,  NZDIM,
     I              1.D0,       0,      0 )
      CALL SHIFT1(    HB,
     D             NXDIM,   NYDIM,      1,
     I              1.D0,       0,      0 )
      CALL SHIFT2(  UBTB,    VBTB,
     D             NXDIM,   NYDIM,      1,
     I             -1.D0,      -1,     -1 )
      CALL SHIFT1(    FT,
     D             NXDIM,   NYDIM,  NTDIM,
     I              1.D0,       0,      0 )
#else
      CALL SHIFT3(
     M                UB,     VB,      W,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT1(
     M                TB, 
     D             NXDIM,  NYDIM, NZTDIM)
      CALL SHIFT3(
     M                HB,   UBTB,   VBTB,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT2(
     M               AMV,    AHV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT1(
     M                FT,       
     D             NXDIM,  NYDIM,  NTDIM)
#endif
      IF (MYRANK .EQ. IROOT) THEN
         CDATE = CHEAD(50)
         READ(CDATE, '(I6.6,5I2.2)') IDATE
         CALL CYH2SS(
     O               TTT,
     I             IDATE)
         WRITE(JFPAR, *) ' Start time :', IDATE
         IF (IRSTRT .NE. 0) THEN
            IF (TTT .NE. TSTRT) THEN
               WRITE(JFPAR, *) '*** Start time error ***'
               CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            END IF
         END IF
      END IF
#else
      DO L = 1, NTDIM
         DO K = 1, NZDIM
            DO J = 1, NYDIM
               DO I = 1, NXDIM
                  TB(I, J, K, L) = 0.D0
               END DO
            END DO
         END DO
      END DO
      DO L = 1, NTDIM
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               FT(I, J, L) = 0.D0
            END DO
         END DO
      END DO

      REWIND(NFINIT)
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 100 K = 1, NZ
         DO 100 J = 1, NY
            DO 100 I = 1, NX
               UB(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
  100 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 110 K = 1, NZ
         DO 110 J = 1, NY
            DO 110 I = 1, NX
               VB(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
  110 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 120 K = 1, NZ
         DO 120 J = 1, NY
            DO 120 I = 1, NX
               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 1) = BUF3(I, J, K)
  120 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 130 K = 1, NZ
         DO 130 J = 1, NY
            DO 130 I = 1, NX
c               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 2) = BUF3(I, J, K)
               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 2) = ABS(BUF3(I, J, K))
  130 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO 140 J = 1, NY
         DO 140 I = 1, NX
            HB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
  140 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO 142 J = 1, NY
         DO 142 I = 1, NX
            UBTB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
  142 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO 144 J = 1, NY
         DO 144 I = 1, NX
            VBTB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
  144 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 150 K = 1, NZ
         DO 150 J = 1, NY
            DO 150 I = 1, NX
               W(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
  150 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 160 K = 1, NZ
         DO 160 J = 1, NY
            DO 160 I = 1, NX
               AMV(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
  160 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO 170 K = 1, NZ
         DO 170 J = 1, NY
            DO 170 I = 1, NX
               AHV(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
  170 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO 180 J = 1, NY
         DO 180 I = 1, NX
            FT(ISTR+I-1, JSTR+J-1, 1) = BUF2(I, J)
  180 CONTINUE
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO 190 J = 1, NY
         DO 190 I = 1, NX
            FT(ISTR+I-1, JSTR+J-1, 2) = BUF2(I, J)
  190 CONTINUE

      DO L = 3, NTDIM
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, L) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      END DO
      DO L = 3, NTDIM
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            FT(ISTR+I-1, JSTR+J-1, L) = BUF2(I, J)
         END DO
      END DO
      END DO

  999 CDATE = CHEAD(50)
      READ(CDATE, '(I6.6,5I2.2)') IDATE

      CALL CYH2SS(
     O               TTT,
     I             IDATE)
      WRITE(JFPAR, *) ' Start time :', IDATE
      IF (IRSTRT .NE. 0) THEN
         IF (TTT .NE. TSTRT) THEN
            WRITE(JFPAR, *) '*** Start time error ***'
            STOP
         END IF
      END IF
#endif

ccc      CALL FILCLS(NFINIT)

      RETURN

c =====================================================================

      ENTRY RSTADD(
     O             ADDITM,   OEOF,
     I              IXDIM,  JYDIM,  KZDIM,
#ifndef OPT_TRIPOLE
     I             CCITEM,   CLAS )
#else
     I             CCITEM,   CLAS,
     I               FACT,   IOFF,   JOFF)
#endif


#ifdef OPT_PARALLEL
      IF (CLAS(1:3) .EQ. 'OCN') THEN
         IF (MYRANK .EQ. IROOT) THEN
            OEOF = .TRUE.
            READ(NFINIT, END=309) CHEAD
            READ(NFINIT) BUF3
            DO K = 1, NZ
               DO J = 1, NYG
                  DO I = 1, NXG
                     G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
                  END DO
               END DO
            END DO
            OEOF = .FALSE.
 309        CONTINUE
         END IF
         CALL MPI_BCAST
     &        (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
         IF (.NOT. OEOF) THEN
            CALL SCATTER_3D(ADDITM, G3D)
#ifdef OPT_TRIPOLE
            CALL SHIFT1(ADDITM, NXDIM, NYDIM, NZDIM,
     &                           FACT,  IOFF,  JOFF)
#else
            CALL SHIFT1(ADDITM, NXDIM, NYDIM, NZDIM)
#endif
         END IF
      ELSE
         IF (MYRANK .EQ. IROOT) THEN
            OEOF = .TRUE.
            READ(NFINIT, END=319) CHEAD
            READ(NFINIT) BUF2
            DO J = 1, NYG
               DO I = 1, NXG
                  G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
               END DO
            END DO
            OEOF = .FALSE.
 319        CONTINUE
         END IF
         CALL MPI_BCAST
     &        (OEOF, 1, MPI_LOGICAL, IROOT, MPI_COMM_WORLD, IERR)
         IF (.NOT. OEOF) THEN
            CALL SCATTER_2D(ADDITM, G2D)
#ifdef OPT_TRIPOLE
            CALL SHIFT1(ADDITM, NXDIM, NYDIM,     1,
     &                           FACT,  IOFF,  JOFF)
#else
            CALL SHIFT1(ADDITM, NXDIM, NYDIM, 1)
#endif
         END IF
      END IF
#else
      IF (CLAS(1:3) .EQ. 'OCN') THEN
         OEOF = .TRUE.
         READ(NFINIT, END=329) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  ADDITM(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
         OEOF = .FALSE.
 329     CONTINUE
      ELSE
         OEOF = .TRUE.
         READ(NFINIT, END=339) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NY
            DO I = 1, NX
               ADDITM(ISTR+I-1, JSTR+J-1, 1) = BUF2(I, J)
            END DO
         END DO
         OEOF = .FALSE.
 339     CONTINUE
      END IF
#endif

      RETURN

c =====================================================================

      ENTRY FINOUT(
     I                 TT, NTSTEP,
     I                 UB,     VB,     TB,     
     I                 HB,   UBTB,   VBTB,
     I                  W,    AMV,    AHV,
     I                 FT,    
     I             ORSOUT, ORSRWD)

      IF (OFIRST) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMFRST, END=98)
 98      CONTINUE
#ifdef OPT_PARALLEL
         IF (MYRANK .EQ. IROOT) THEN
            CALL FILOPN(NFREST, CFREST, 'WRITE')
         END IF
#else
         CALL FILOPN(NFREST, CFREST, 'WRITE')
#endif
         OFIRST = .FALSE.
      END IF

      IF (.NOT. ORSOUT) RETURN

#ifdef OPT_PARALLEL
      IF (MYRANK .EQ. IROOT) THEN
         IF (ORSRWD) THEN
            REWIND(NFREST)
         END IF
         CALL CSS2YH(
     O             IDATE,
     I                TT)
         WRITE(CHEAD(50), '(I6.6,5I2.2)') IDATE
         WRITE(CHEAD(27), '(I4.4,2I2.2,A1,3I2.2,A1)')
     &        IDATE(1), IDATE(2), IDATE(3), ' ',
     &        IDATE(4), IDATE(5), IDATE(6), ' '
         WRITE(CHEAD(29), '(I15,A1)') NXG, 'X'
         WRITE(CHEAD(30), '(I16)') 1
         WRITE(CHEAD(31), '(I16)') NXG
         WRITE(CHEAD(32), '(I15,A1)') NYG, 'Y'
         WRITE(CHEAD(33), '(I16)') 1
         WRITE(CHEAD(34), '(I16)') NYG
         WRITE(CHEAD(36), '(I16)') 1
         CHEAD(38) = 'REAL8'
      END IF
#else
      IF (ORSRWD) THEN
         REWIND(NFREST)
      END IF
      CALL CSS2YH(
     O             IDATE,
     I                TT)
      WRITE(CHEAD(27), '(I4.4,2I2.2,A1,3I2.2,A1)')
     &              IDATE(1), IDATE(2), IDATE(3), ' ',
     &              IDATE(4), IDATE(5), IDATE(6), ' '
      WRITE(CHEAD(50), '(I6.6,5I2.2)') IDATE
      WRITE(CHEAD(29), '(I15,A1)') NX, 'X'
      WRITE(CHEAD(30), '(I16)') 1
      WRITE(CHEAD(31), '(I16)') NX
      WRITE(CHEAD(32), '(I15,A1)') NY, 'Y'
      WRITE(CHEAD(33), '(I16)') 1
      WRITE(CHEAD(34), '(I16)') NY
      WRITE(CHEAD(36), '(I16)') 1
      CHEAD(38) = 'REAL8'
#endif

#ifdef OPT_PARALLEL
      CALL GATHER_3D(G3D, UB)
      IF (MYRANK .EQ. IROOT) THEN
         DO 200 K = 1, NZ
            DO 200 J = 1, NYG
               DO 200 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 200     CONTINUE
         CHEAD(3) = 'U'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, VB)
      IF (MYRANK .EQ. IROOT) THEN
         DO 210 K = 1, NZ
            DO 210 J = 1, NYG
               DO 210 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 210     CONTINUE
         CHEAD(3) = 'V'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, TB)
      IF (MYRANK .EQ. IROOT) THEN
         DO 220 K = 1, NZ
            DO 220 J = 1, NYG
               DO 220 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 220     CONTINUE
         CHEAD(3) = 'T'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, TB(1, 1, 1, 2))
      IF (MYRANK .EQ. IROOT) THEN
         DO 230 K = 1, NZ
            DO 230 J = 1, NYG
               DO 230 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 230     CONTINUE
         CHEAD(3) = 'S'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_2D(G2D, HB)
      IF (MYRANK .EQ. IROOT) THEN
         DO 240 J = 1, NYG
            DO 240 I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
 240     CONTINUE
         CHEAD(3) = 'SH'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, UBTB)
      IF (MYRANK .EQ. IROOT) THEN
         DO 242 J = 1, NYG
            DO 242 I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
 242     CONTINUE
         CHEAD(3) = 'UBT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, VBTB)
      IF (MYRANK .EQ. IROOT) THEN
         DO 244 J = 1, NYG
            DO 244 I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
 244     CONTINUE
         CHEAD(3) = 'VBT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_3D(G3D, W)
      IF (MYRANK .EQ. IROOT) THEN
         DO 250 K = 1, NZ
            DO 250 J = 1, NYG
               DO 250 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 250     CONTINUE
         CHEAD(3) = 'W'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, AMV)
      IF (MYRANK .EQ. IROOT) THEN
         DO 260 K = 1, NZ
            DO 260 J = 1, NYG
               DO 260 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 260     CONTINUE
         CHEAD(3) = 'AMV'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, AHV)
      IF (MYRANK .EQ. IROOT) THEN
         DO 270 K = 1, NZ
            DO 270 J = 1, NYG
               DO 270 I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
 270     CONTINUE
         CHEAD(3) = 'AHV'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_2D(G2D, FT)
      IF (MYRANK .EQ. IROOT) THEN
         DO 280 J = 1, NYG
            DO 280 I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
 280     CONTINUE
         CHEAD(3) = 'FT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, FT(1, 1, 2))
      IF (MYRANK .EQ. IROOT) THEN
         DO 290 J = 1, NYG
            DO 290 I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
 290     CONTINUE
         CHEAD(3) = 'FW'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      DO L = 3, NTDIM
         CALL GATHER_3D(G3D, TB(1, 1, 1, L))
         IF (MYRANK .EQ. IROOT) THEN
            DO K = 1, NZ
               DO J = 1, NYG
                  DO I = 1, NXG
                     BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
                  END DO
               END DO
            END DO
            WRITE(CHEAD(3), '(A6,I2.2)') 'TRACER', L
            WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
            WRITE(CHEAD(37), '(I16)') NZ
            WRITE(CHEAD(64), '(I16)') NXYZG
            WRITE(NFREST) CHEAD
            WRITE(NFREST) BUF3
         END IF
      END DO

      DO L = 3, NTDIM
         CALL GATHER_2D(G2D, FT(1, 1, L))
         IF (MYRANK .EQ. IROOT) THEN
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
               END DO
            END DO
            WRITE(CHEAD(3), '(A6,I2.2)') 'TRCFLX', L
            WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
            WRITE(CHEAD(37), '(I16)') 1
            WRITE(CHEAD(64), '(I16)') NXYG
            WRITE(NFREST) CHEAD
            WRITE(NFREST) BUF2
         END IF
      END DO

      IF (MYRANK .EQ. IROOT) THEN
         WRITE(JFPAR, *) '*** Write restart file ***'
         WRITE(JFPAR, *) ' Time :', IDATE
         WRITE(JFPAR, *) ' Step :', NTSTEP
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
#else
      DO 200 K = 1, NZ
         DO 200 J = 1, NY
            DO 200 I = 1, NX
               BUF3(I, J, K) = UB(ISTR+I-1, JSTR+J-1, KSTR+K-1)
  200 CONTINUE
      CHEAD(3) = 'U'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 210 K = 1, NZ
         DO 210 J = 1, NY
            DO 210 I = 1, NX
               BUF3(I, J, K) = VB(ISTR+I-1, JSTR+J-1, KSTR+K-1)
  210 CONTINUE
      CHEAD(3) = 'V'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 220 K = 1, NZ
         DO 220 J = 1, NY
            DO 220 I = 1, NX
               BUF3(I, J, K) = TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 1)
  220 CONTINUE
      CHEAD(3) = 'T'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 230 K = 1, NZ
         DO 230 J = 1, NY
            DO 230 I = 1, NX
               BUF3(I, J, K) = TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 2)
  230 CONTINUE
      CHEAD(3) = 'S'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 240 J = 1, NY
         DO 240 I = 1, NX
            BUF2(I, J) = HB(ISTR+I-1, JSTR+J-1)
  240 CONTINUE
      CHEAD(3) = 'SH'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO 242 J = 1, NY
         DO 242 I = 1, NX
            BUF2(I, J) = UBTB(ISTR+I-1, JSTR+J-1)
  242 CONTINUE
      CHEAD(3) = 'UBT'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO 244 J = 1, NY
         DO 244 I = 1, NX
            BUF2(I, J) = VBTB(ISTR+I-1, JSTR+J-1)
  244 CONTINUE
      CHEAD(3) = 'VBT'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO 250 K = 1, NZ
         DO 250 J = 1, NY
            DO 250 I = 1, NX
               BUF3(I, J, K) = W(ISTR+I-1, JSTR+J-1, KSTR+K-1)
  250 CONTINUE
      CHEAD(3) = 'W'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 260 K = 1, NZ
         DO 260 J = 1, NY
            DO 260 I = 1, NX
               BUF3(I, J, K) = AMV(ISTR+I-1, JSTR+J-1, KSTR+K-1)
  260 CONTINUE
      CHEAD(3) = 'AMV'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 270 K = 1, NZ
         DO 270 J = 1, NY
            DO 270 I = 1, NX
               BUF3(I, J, K) = AHV(ISTR+I-1, JSTR+J-1, KSTR+K-1)
  270 CONTINUE
      CHEAD(3) = 'AHV'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO 280 J = 1, NY
         DO 280 I = 1, NX
            BUF2(I, J) = FT(ISTR+I-1, JSTR+J-1, 1)
  280 CONTINUE
      CHEAD(3) = 'FT'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO 290 J = 1, NY
         DO 290 I = 1, NX
            BUF2(I, J) = FT(ISTR+I-1, JSTR+J-1, 2)
  290 CONTINUE
      CHEAD(3) = 'FW'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO L = 3, NTDIM
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  BUF3(I, J, K) = TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, L)
               END DO
            END DO
         END DO
         WRITE(CHEAD(3), '(A6,I2.2)') 'TRACER', L
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZ
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END DO

      DO L = 3, NTDIM
         DO J = 1, NY
            DO I = 1, NX
               BUF2(I, J) = FT(ISTR+I-1, JSTR+J-1, L)
            END DO
         END DO
         WRITE(CHEAD(3), '(A6,I2.2)') 'TRCFLX', L
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXY
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END DO

      WRITE(JFPAR, *) '*** Write restart file ***'
      WRITE(JFPAR, *) ' Time :', IDATE
      WRITE(JFPAR, *) ' Step :', NTSTEP
#endif

      RETURN

c =====================================================================

      ENTRY FINADD(
     I             ADDITM,  IXDIM,  JYDIM,  KZDIM,
     I             CCITEM,   CLAS)

#ifdef OPT_PARALLEL
      IF (CLAS(1:3) .EQ. 'OCN') THEN
         CALL GATHER_3D(G3D, ADDITM)
         IF (MYRANK .EQ. IROOT) THEN
            DO K = 1, NZ
               DO J = 1, NYG
                  DO I = 1, NXG
                     BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
                  END DO
               END DO
            END DO
            CHEAD(3) = CCITEM
            WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
            WRITE(CHEAD(37), '(I16)') NZ
            WRITE(CHEAD(64), '(I16)') NXYZG
            WRITE(NFREST) CHEAD
            WRITE(NFREST) BUF3
         END IF
      ELSE
         CALL GATHER_2D(G2D, ADDITM)
         IF (MYRANK .EQ. IROOT) THEN
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
               END DO
            END DO
            CHEAD(3) = CCITEM
            WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
            WRITE(CHEAD(37), '(I16)') 1
            WRITE(CHEAD(64), '(I16)') NXYG
            WRITE(NFREST) CHEAD
            WRITE(NFREST) BUF2
         END IF
      END IF
#else
      IF (CLAS(1:3) .EQ. 'OCN') THEN
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  BUF3(I, J, K) = ADDITM(ISTR+I-1, JSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = CCITEM
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZ
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      ELSE
         DO J = 1, NY
            DO I = 1, NX
               BUF2(I, J) = ADDITM(ISTR+I-1, JSTR+J-1, 1)
            END DO
         END DO
         CHEAD(3) = 'FT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXY
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF
#endif

      RETURN
      END
