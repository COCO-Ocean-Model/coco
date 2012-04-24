      SUBROUTINE RESTRT(
     I             TSTRT,
     O                UB,     VB,     TB,     
     O                HB,   UBTB,   VBTB,
     O                 W,    AMV,    AHV,
     O                AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     O               TSI,
     O                FT,  SWABS,     FS,
     O              TAUX,   TAUY,   PTOP)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Reading and writing restart file
c
c  HISTORY
c     '02.10.10  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: removal of negative salinity, which is
c                          sometime used for missing values
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '07.09.27  H.Hasumi: 1-layer sea ice thermodynamics
c     '07.11.06  T.Suzuki: for ES
c     '08.11.25  Y.Komuro: default value of TSI changed
c                          eliminate unphysical negative AB(IJ, 0)
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif
#include "zocdim.F"
#include "zocphy.F"
#include "zocfil.F"

      REAL*8      UB(NXDIM, NYDIM, NZDIM),     VB(NXDIM, NYDIM, NZDIM)
      REAL*8      TB(NXDIM, NYDIM, NZDIM, NTDIM)
      REAL*8      HB(NXDIM, NYDIM)
      REAL*8    UBTB(NXDIM, NYDIM)       ,   VBTB(NXDIM, NYDIM)
      REAL*8       W(NXDIM, NYDIM, NZDIM)
      REAL*8     AMV(NXDIM, NYDIM, NZDIM),    AHV(NXDIM, NYDIM, NZDIM)
      REAL*8      AB(NXDIM, NYDIM, 0:NIC),    HIB(NXDIM, NYDIM, 0:NIC)
      REAL*8     HSB(NXDIM, NYDIM, 0:NIC),    TSI(NXDIM, NYDIM, 0:NIC)
      REAL*8     UIB(NXDIM, NYDIM)       ,    VIB(NXDIM, NYDIM)
      REAL*8     TIB(NXDIM, NYDIM, 0:NIC)
      REAL*8      FT(NXDIM, NYDIM, NTDIM)
      REAL*8   SWABS(NXDIM, NYDIM)       ,     FS(NXDIM, NYDIM)
      REAL*8    TAUX(NXDIM, NYDIM)       ,   TAUY(NXDIM, NYDIM)
      REAL*8    PTOP(NXDIM, NYDIM)

      INTEGER   IXDIM,  JYDIM,  KZDIM
      LOGICAL    OEOF
      REAL*8     ADDITM(IXDIM, JYDIM, KZDIM)
      CHARACTER  CCITEM*(*), CLAS*(*)

#ifdef OPT_PARALLEL
C      REAL*8    BUF2(NXG, NYG)      ,   BUF3(NXG, NYG, NZ)
C      REAL*8    BUFI(NXG, NYG, NIC)
C      REAL*8     G2D(NXGDIM, NYGDIM),    G3D(NXGDIM, NYGDIM, NZDIM)
C      REAL*8     GID(NXGDIM, NYGDIM, 0:NIC)
C      COMMON /WORK/ BUF2, BUF3, BUFI, G2D, G3D, GID
      REAL(8),ALLOCATABLE :: BUF2(:,:), BUF3(:,:,:), BUFI(:,:,:)
      REAL(8),ALLOCATABLE :: G2D(:,:), G3D(:,:,:), GID(:,:,:)
      SAVE BUF2, BUF3, BUFI, G2D, G3D, GID
#else
      REAL*8    BUF2(NX, NY),   BUF3(NX, NY, NZ)
      REAL*8    BUFI(NX, NY, NIC)
      COMMON /WORK/ BUF2, BUF3, BUFI
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
      REAL*8    SI
      SAVE NFINIT, NFREST, SI
      NAMELIST /NMFINI/ CFINIT, IRSTRT
      NAMELIST /NMFRST/ CFREST
      NAMELIST /NMISLT/ SI
      DATA SI / 5.D0 /

      DATA OFIRST / .TRUE. /
      DATA CFINIT, CFREST / 'not-specified', 'not-specified' /
      DATA IRSTRT / 0 /

      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMFINI, END=99)
   99 CONTINUE
      CALL REWNML(IFPAR, JFPAR)
      READ(IFPAR, NMISLT, END=97)
   97 WRITE(JFPAR, NMISLT)

#ifdef OPT_PARALLEL
      IF(MYRANK .EQ. IROOT) THEN
      CALL FILOPN(NFINIT, CFINIT, 'READ')
      WRITE(JFPAR, *) '*** Read initial condition file ***'
      ENDIF
#else
      CALL FILOPN(NFINIT, CFINIT, 'READ')
      WRITE(JFPAR, *) '*** Read initial condition file ***'
#endif

#ifdef OPT_PARALLEL
      TB = 0.D0
      FT = 0.D0
      SWABS = 0.D0
      FS = 0.D0
      TAUX = 0.D0
      TAUY = 0.D0
      AMV = 0.D0
      AHV = 0.D0
      PTOP = 0.D0
      TSI = dtds * SI
      HIB = 0.D0
      HSB = 0.D0
      TIB = -0.1D0

      IF(MYRANK.NE.IROOT) THEN
         ALLOCATE(BUF2(1,1))
         ALLOCATE(BUF3(1,1,1))
         ALLOCATE(BUFI(1,1,1))
         ALLOCATE(G2D(1,1))
         ALLOCATE(G3D(1,1,1))
         ALLOCATE(GID(1,1,1))
      ELSE
         ALLOCATE(BUF2(NXG,NYG))
         ALLOCATE(BUF3(NXG,NYG,NZ))
         ALLOCATE(BUFI(NXG,NYG,NIC))
         ALLOCATE(G2D(NXGDIM,NYGDIM))
         ALLOCATE(G3D(NXGDIM,NYGDIM,NZDIM))
         ALLOCATE(GID(NXGDIM,NYGDIM,0:NIC))
      ENDIF

      IF (MYRANK .EQ. IROOT) THEN
         REWIND(NFINIT)
         READ(NFINIT, END=109) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
 109     CONTINUE
      END IF
      CALL SCATTER_3D(UB, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=119) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
 119     CONTINUE
      END IF
      CALL SCATTER_3D(VB, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=129) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
 129     CONTINUE
      END IF
      CALL SCATTER_3D(TB, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=139) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
c                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) =
     &            ABS(BUF3(I, J, K))
               END DO
            END DO
         END DO
 139     CONTINUE
      END IF
      CALL SCATTER_3D(TB(1, 1, 1, 2), G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=141) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 141     CONTINUE
      END IF
      CALL SCATTER_2D(HB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=143) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 143     CONTINUE
      END IF
      CALL SCATTER_2D(UBTB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=145) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 145     CONTINUE
      END IF
      CALL SCATTER_2D(VBTB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=159) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
 159     CONTINUE
      END IF
      CALL SCATTER_3D(W, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=301) CHEAD
         READ(NFINIT) BUFI
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  GID(IGSTR+I-1, JGSTR+J-1, K) = BUFI(I, J, K)
               END DO
            END DO
         END DO
 301     CONTINUE
      END IF
      CALL SCATTER_ID(AB, GID)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=311) CHEAD
         READ(NFINIT) BUFI
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  GID(IGSTR+I-1, JGSTR+J-1, K) = BUFI(I, J, K)
               END DO
            END DO
         END DO
 311     CONTINUE
      END IF
      CALL SCATTER_ID(HIB, GID)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=321) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 321     CONTINUE
      END IF
      CALL SCATTER_2D(UIB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=331) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 331     CONTINUE
      END IF
      CALL SCATTER_2D(VIB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=336) CHEAD
         READ(NFINIT) BUFI
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  GID(IGSTR+I-1, JGSTR+J-1, K) = BUFI(I, J, K)
               END DO
            END DO
         END DO
 336     CONTINUE
      END IF
      CALL SCATTER_ID(TIB, GID)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=341) CHEAD
         READ(NFINIT) BUFI
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  GID(IGSTR+I-1, JGSTR+J-1, K) = BUFI(I, J, K)
               END DO
            END DO
         END DO
 341     CONTINUE
      END IF
      CALL SCATTER_ID(HSB, GID)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=189) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 189     CONTINUE
      END IF
      CALL SCATTER_2D(FT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=401) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 401     CONTINUE
      END IF
      CALL SCATTER_2D(SWABS, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=199) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 199     CONTINUE
      END IF
      CALL SCATTER_2D(FT(1, 1, 2), G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=411) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 411     CONTINUE
      END IF
      CALL SCATTER_2D(FS, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=421) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 421     CONTINUE
      END IF
      CALL SCATTER_2D(TAUX, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=431) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 431     CONTINUE
      END IF
      CALL SCATTER_2D(TAUY, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=169) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
 169     CONTINUE
      END IF
      CALL SCATTER_3D(AMV, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=179) CHEAD
         READ(NFINIT) BUF3
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
               END DO
            END DO
         END DO
 179     CONTINUE
      END IF
      CALL SCATTER_3D(AHV, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=441) CHEAD
         READ(NFINIT) BUF2
         DO J = 1, NYG
            DO I = 1, NXG
               G2D(IGSTR+I-1, JGSTR+J-1) = BUF2(I, J)
            END DO
         END DO
 441     CONTINUE
      END IF
      CALL SCATTER_2D(PTOP, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFINIT, END=451) CHEAD
         READ(NFINIT) BUFI
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  GID(IGSTR+I-1, JGSTR+J-1, K) = BUFI(I, J, K)
               END DO
            END DO
         END DO
 451     CONTINUE
      END IF
      CALL SCATTER_ID(TSI, GID)

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
      CALL SHIFT3(    AB,     HIB,    HSB,
     D             NXDIM,   NYDIM,  NIC+1,
     I              1.D0,       0,      0 )
      CALL SHIFT2(   TSI,     TIB, 
     D             NXDIM,   NYDIM,  NIC+1,
     I              1.D0,       0,      0 )
      CALL SHIFT2(   UIB,     VIB,
     D             NXDIM,   NYDIM,      1,
     I             -1.D0,      -1,     -1 )
      CALL SHIFT1(    FT,
     D             NXDIM,   NYDIM,  NTDIM,
     I              1.D0,       0,      0 )
      CALL SHIFT3(    FS,   SWABS,   PTOP,
     D             NXDIM,   NYDIM,      1,
     I              1.D0,       0,      0 )
      CALL SHIFT2(  TAUX,    TAUY,
     D             NXDIM,   NYDIM,      1,
     I             -1.D0,      -1,     -1 )
#else
      CALL SHIFT3(
     M                UB,     VB,      W,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT1(
     M                TB, 
     D             NXDIM,  NYDIM, NZTDIM)
      CALL SHIFT2(
     M               AMV,    AHV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M                HB,   UBTB,   VBTB,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M                AB,    HIB,    HSB,
     D             NXDIM,  NYDIM,  NIC+1)
      CALL SHIFT2(
     M               TSI,    TIB,
     D             NXDIM,  NYDIM,  NIC+1)
      CALL SHIFT2(
     M               UIB,    VIB,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT1(
     M                FT,       
     D             NXDIM,  NYDIM,  NTDIM)
      CALL SHIFT2(
     M             SWABS,     FS,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M              TAUX,   TAUY,   PTOP,
     D             NXDIM,  NYDIM,      1)
#endif

      DO J = 1, NYDIM
         DO I = 1, NXDIM
            AB(I, J, 0) = 1.D0
         END DO
      END DO
      DO K = 1, NIC
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               AB(I, J, 0) = AB(I, J, 0) - AB(I, J, K)
            END DO
         END DO
      END DO

*** because of the finit precision, AB(IJ, 0) could be a small negative
*** value, which causes some problems
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            AB(I, J, 0) = MAX(0.D0, AB(I, J ,0))
         END DO
      END DO

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
      DO K = 1, NZDIM
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               AMV(I, J, K) = 0.D0
               AHV(I, J, K) = 0.D0
            END DO
         END DO
      END DO
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            SWABS(I, J) = 0.D0
            FS(I, J) = 0.D0
            TAUX(I, J) = 0.D0
            TAUY(I, J) = 0.D0
            PTOP(I, J) = 0.D0
         END DO
      END DO
      DO K = 0, NIC
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               HIB(I, J, K) = 0.D0
               HSB(I, J, K) = 0.D0
               TSI(I, J, K) = dtds * SI
               TIB(I, J, K) = -0.1D0
            END DO
         END DO
      END DO

      REWIND(NFINIT)
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               UB(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               VB(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 1) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
c               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 2) = BUF3(I, J, K)
               TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 2) = ABS(BUF3(I, J, K))
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            HB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            UBTB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            VBTB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               W(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUFI
      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               AB(ISTR+I-1, JSTR+J-1, K) = BUFI(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUFI
      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               HIB(ISTR+I-1, JSTR+J-1, K) = BUFI(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            UIB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            VIB(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUFI
      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               TIB(ISTR+I-1, JSTR+J-1, K) = BUFI(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUFI
      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               HSB(ISTR+I-1, JSTR+J-1, K) = BUFI(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            FT(ISTR+I-1, JSTR+J-1, 1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            SWABS(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            FT(ISTR+I-1, JSTR+J-1, 2) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            FS(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            TAUX(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            TAUY(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               AMV(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF3
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               AHV(ISTR+I-1, JSTR+J-1, KSTR+K-1) = BUF3(I, J, K)
            END DO
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUF2
      DO J = 1, NY
         DO I = 1, NX
            PTOP(ISTR+I-1, JSTR+J-1) = BUF2(I, J)
         END DO
      END DO
      READ(NFINIT, END=999) CHEAD
      READ(NFINIT) BUFI
      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               TSI(ISTR+I-1, JSTR+J-1, K) = BUFI(I, J, K)
            END DO
         END DO
      END DO

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

      CALL STBCSI(AB, HIB, TIB, TSI, HSB) ! TSI is dummy
      CALL STBCVI(UIB, VIB)
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            AB(I, J, 0) = 1.D0
         END DO
      END DO
      DO K = 1, NIC
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               AB(I, J, 0) = AB(I, J, 0) - AB(I, J, K)
            END DO
         END DO
      END DO

*** because of the finit precision, AB(IJ, 0) could be a small negative
*** value, which causes some problems
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            AB(I, J, 0) = MAX(0.D0, AB(I, J ,0))
         END DO
      END DO


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
     I                 AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     I                TSI,
     I                 FT,  SWABS,     FS,
     I               TAUX,   TAUY,   PTOP,
     I             ORSOUT, ORSRWD)

      IF (OFIRST) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMFRST, END=98)
   98    CONTINUE
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
         REWIND(UNIT=NFREST)
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
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = 'U'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, VB)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO   
         CHEAD(3) = 'V'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, TB)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = 'T'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, TB(1, 1, 1, 2))
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = 'S'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_2D(G2D, HB)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'SH'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, UBTB)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'UBT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, VBTB)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'VBT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_3D(G3D, W)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = 'W'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_ID(GID, AB)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  BUFI(I, J, K) = GID(IGSTR+I-1, JGSTR+J-1, K)
               END DO
            END DO
         END DO
         CHEAD(3) = 'AI'
         WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
         WRITE(CHEAD(37), '(I16)') NIC
         WRITE(CHEAD(64), '(I16)') NXYG*NIC
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUFI
      END IF

      CALL GATHER_ID(GID, HIB)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  BUFI(I, J, K) = GID(IGSTR+I-1, JGSTR+J-1, K)
               END DO
            END DO
         END DO
         CHEAD(3) = 'HI'
         WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
         WRITE(CHEAD(37), '(I16)') NIC
         WRITE(CHEAD(64), '(I16)') NXYG*NIC
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUFI
      END IF

      CALL GATHER_2D(G2D, UIB)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'UI'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, VIB)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'VI'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_ID(GID, TIB)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  BUFI(I, J, K) = GID(IGSTR+I-1, JGSTR+J-1, K)
               END DO
            END DO
         END DO
         CHEAD(3) = 'TI'
         WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
         WRITE(CHEAD(37), '(I16)') NIC
         WRITE(CHEAD(64), '(I16)') NXYG*NIC
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUFI
      END IF

      CALL GATHER_ID(GID, HSB)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  BUFI(I, J, K) = GID(IGSTR+I-1, JGSTR+J-1, K)
               END DO
            END DO
         END DO
         CHEAD(3) = 'HS'
         WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
         WRITE(CHEAD(37), '(I16)') NIC
         WRITE(CHEAD(64), '(I16)') NXYG*NIC
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUFI
      END IF

      CALL GATHER_2D(G2D, FT)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'FT'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, SWABS)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'SWABS'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, FT(1, 1, 2))
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'FW'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, FS)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'FS'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, TAUX)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'TAUX'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_2D(G2D, TAUY)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'TAUY'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_3D(G3D, AMV)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = 'AMV'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_3D(G3D, AHV)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NZ
            DO J = 1, NYG
               DO I = 1, NXG
                  BUF3(I, J, K) = G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1)
               END DO
            END DO
         END DO
         CHEAD(3) = 'AHV'
         WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
         WRITE(CHEAD(37), '(I16)') NZ
         WRITE(CHEAD(64), '(I16)') NXYZG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF3
      END IF

      CALL GATHER_2D(G2D, PTOP)
      IF (MYRANK .EQ. IROOT) THEN
         DO J = 1, NYG
            DO I = 1, NXG
               BUF2(I, J) = G2D(IGSTR+I-1, JGSTR+J-1)
            END DO
         END DO
         CHEAD(3) = 'PTOP'
         WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
         WRITE(CHEAD(37), '(I16)') 1
         WRITE(CHEAD(64), '(I16)') NXYG
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUF2
      END IF

      CALL GATHER_ID(GID, TSI)
      IF (MYRANK .EQ. IROOT) THEN
         DO K = 1, NIC
            DO J = 1, NYG
               DO I = 1, NXG
                  BUFI(I, J, K) = GID(IGSTR+I-1, JGSTR+J-1, K)
               END DO
            END DO
         END DO
         CHEAD(3) = 'TSI'
         WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
         WRITE(CHEAD(37), '(I16)') NIC
         WRITE(CHEAD(64), '(I16)') NXYG*NIC
         WRITE(NFREST) CHEAD
         WRITE(NFREST) BUFI
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
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = UB(ISTR+I-1, JSTR+J-1, KSTR+K-1)
            END DO
         END DO
      END DO
      CHEAD(3) = 'U'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = VB(ISTR+I-1, JSTR+J-1, KSTR+K-1)
            END DO
         END DO
      END DO
      CHEAD(3) = 'V'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 1)
            END DO
         END DO
      END DO
      CHEAD(3) = 'T'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = TB(ISTR+I-1, JSTR+J-1, KSTR+K-1, 2)
            END DO
         END DO
      END DO
      CHEAD(3) = 'S'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = HB(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'SH'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = UBTB(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'UBT'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = VBTB(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'VBT'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = W(ISTR+I-1, JSTR+J-1, KSTR+K-1)
            END DO
         END DO
      END DO
      CHEAD(3) = 'W'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               BUFI(I, J, K) = AB(ISTR+I-1, JSTR+J-1, K)
            END DO
         END DO
      END DO
      CHEAD(3) = 'AI'
      WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
      WRITE(CHEAD(37), '(I16)') NIC
      WRITE(CHEAD(64), '(I16)') NXY*NIC
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUFI

      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               BUFI(I, J, K) = HIB(ISTR+I-1, JSTR+J-1, K)
            END DO
         END DO
      END DO
      CHEAD(3) = 'HI'
      WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
      WRITE(CHEAD(37), '(I16)') NIC
      WRITE(CHEAD(64), '(I16)') NXY*NIC
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUFI

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = UIB(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'UI'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = VIB(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'VI'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               BUFI(I, J, K) = TIB(ISTR+I-1, JSTR+J-1, K)
            END DO
         END DO
      END DO
      CHEAD(3) = 'TI'
      WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
      WRITE(CHEAD(37), '(I16)') NIC
      WRITE(CHEAD(64), '(I16)') NXY*NIC
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUFI

      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               BUFI(I, J, K) = HSB(ISTR+I-1, JSTR+J-1, K)
            END DO
         END DO
      END DO
      CHEAD(3) = 'HS'
      WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
      WRITE(CHEAD(37), '(I16)') NIC
      WRITE(CHEAD(64), '(I16)') NXY*NIC
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUFI

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = FT(ISTR+I-1, JSTR+J-1, 1)
         END DO
      END DO
      CHEAD(3) = 'FT'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = SWABS(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'SWABS'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = FT(ISTR+I-1, JSTR+J-1, 2)
         END DO
      END DO
      CHEAD(3) = 'FW'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = FS(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'FS'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = TAUX(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'TAUX'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = TAUY(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'TAUY'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = AMV(ISTR+I-1, JSTR+J-1, KSTR+K-1)
            END DO
         END DO
      END DO
      CHEAD(3) = 'AMV'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               BUF3(I, J, K) = AHV(ISTR+I-1, JSTR+J-1, KSTR+K-1)
            END DO
         END DO
      END DO
      CHEAD(3) = 'AHV'
      WRITE(CHEAD(35), '(I15,A1)') NZ, 'Z'
      WRITE(CHEAD(37), '(I16)') NZ
      WRITE(CHEAD(64), '(I16)') NXYZ
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF3

      DO J = 1, NY
         DO I = 1, NX
            BUF2(I, J) = PTOP(ISTR+I-1, JSTR+J-1)
         END DO
      END DO
      CHEAD(3) = 'PTOP'
      WRITE(CHEAD(35), '(I15,A1)') 1, 'Z'
      WRITE(CHEAD(37), '(I16)') 1
      WRITE(CHEAD(64), '(I16)') NXY
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUF2

      DO K = 1, NIC
         DO J = 1, NY
            DO I = 1, NX
               BUFI(I, J, K) = TSI(ISTR+I-1, JSTR+J-1, K)
            END DO
         END DO
      END DO
      CHEAD(3) = 'TSI'
      WRITE(CHEAD(35), '(I15,A1)') NIC, 'Z'
      WRITE(CHEAD(37), '(I16)') NIC
      WRITE(CHEAD(64), '(I16)') NXY*NIC
      WRITE(NFREST) CHEAD
      WRITE(NFREST) BUFI

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
