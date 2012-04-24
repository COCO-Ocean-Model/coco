      SUBROUTINE RDGEO
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Reading grid information and masking array
c
c  HISTORY
c     '03.04.21  H.Hasumi: from COCO3.4
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '07.11.06  T.Suzuki: for ES
c     '07.12.20  H.Hasumi: MASK file format is modified
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif
#include "zocdim.F"
#include "zocmsk.F"
#include "zocgrd.F"
#include "zocfil.F"

#ifdef OPT_PARALLEL
      INTEGER,ALLOCATABLE :: G2DI(:)
      REAL(8),ALLOCATABLE :: G2D(:), G3D(:)
      INTEGER    N
#endif

      INTEGER     IJ,      K
      INTEGER      I,      J,    IJS,   IJE,   IJN,  IJW
      INTEGER  IFPAR,  JFPAR

      INTEGER   NFMASK
      CHARACTER CFMASK*(NCF)
      NAMELIST /NMMASK/ CFMASK
      DATA CFMASK / 'not-specified' /

      CALL REWNML(IFPAR, JFPAR)
      WRITE(JFPAR, *) '*** RDGEO ***'
      READ(IFPAR, NMMASK, END=99)
 99   CONTINUE

#ifdef OPT_PARALLEL
      IF (MYRANK .EQ. IROOT) THEN
        CALL FILOPN(NFMASK, CFMASK, 'READ')
      ENDIF
#else
        CALL FILOPN(NFMASK, CFMASK, 'READ')
#endif

#ifdef OPT_PARALLEL
      AMSKT = 0.D0
      AMFTX = 0.D0
      AMFTY = 0.D0
      AMFTZ = 0.D0
      AMSKV = 0.D0
      AMFVX = 0.D0
      AMFVY = 0.D0
      AMFVZ = 0.D0
      AMSKB = 0.D0
      DY = 1.D0
      DYM = 1.D0
      DZ  = 1.D0
      DZM = 1.D0
      DZV = 1.D0
      COR = 1.D0
      HXT = 1.D0
      HXU = 1.D0
      HYT = 1.D0
      HYU = 1.D0
      HXYT = 1.D0
      HXYU = 1.D0
      HYXT = 1.D0
      HYXU = 1.D0
      DEPT = 1.D0
      RDEPV = 1.D0
      NBOT = KSTR - 1
#ifdef OPT_BBL
      NBOTV = KSTR - 1
      AMSKTB = 0.D0
      AMSKT0 = 0.D0
      AMSKT1 = 0.D0
      AMSKVB = 0.D0
      AMSKV0 = 0.D0
      AMSKV1 = 0.D0
#endif
#endif

#ifdef OPT_PARALLEL
      IF (MYRANK .EQ. IROOT) THEN
         ALLOCATE(G2D(NXYGDM))
         ALLOCATE(G3D(NXYZGD))
         ALLOCATE(G2DI(NXYGDM))
      ELSE
         ALLOCATE(G2D(1))
         ALLOCATE(G3D(1))
         ALLOCATE(G2DI(1))
      END IF

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) DX
      END IF
      CALL MPI_BCAST(DX, 1, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(DY, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(DYM, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(DZ, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(DZM, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(DZV, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) DZ0
      END IF
      CALL MPI_BCAST(DZ0, NZDIM, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) DS
      END IF
      CALL MPI_BCAST(DS, NZDIM, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) DSM
      END IF
      CALL MPI_BCAST(DSM, NZDIM, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) HIC
      END IF
      CALL MPI_BCAST(HIC, NIC+1, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) REA
      END IF
      CALL MPI_BCAST(REA, 1, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) ZBOT
      END IF
      CALL MPI_BCAST(ZBOT, 1, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(COR, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(DEPT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(RDEPV, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HXT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HXU, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HYT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HYU, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HXYT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HXYU, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HYXT, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(HYXU, G2D)

#ifdef OPT_TRIPOLE
      CALL SHIFT2(    DY,     DEPT, 
     D             NXDIM,  NYDIM,      1,
     I              1.D0,      0,      0 )
      CALL SHIFT3(   DYM,    COR,  RDEPV,
     D             NXDIM,  NYDIM,      1,
     I              1.D0,     -1,     -1 )
      CALL SHIFT2(    DZ,    DZM, 
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,      0,      0 )
      CALL SHIFT1(   DZV, 
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,     -1,     -1 )
      CALL SHIFT2(   HXT,    HYT, 
     D             NXDIM,  NYDIM,      1,
     I              1.D0,      0,      0 )
      CALL SHIFT2(   HXU,    HYU, 
     D             NXDIM,  NYDIM,      1,
     I              1.D0,     -1,     -1 )
      CALL SHIFT2(  HXYT,   HYXT, 
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
      CALL SHIFT2(  HXYU,   HYXU, 
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,     -1,     -1 )
#else
      CALL SHIFT2(
     M                DY,    DYM,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M                DZ,    DZM,    DZV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M               COR,   DEPT,  RDEPV,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M               HXT,    HXU,    HYT,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M               HYU,   HXYT,   HXYU,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT2(
     M              HYXT,   HYXU,
     D             NXDIM,  NYDIM,      1)

      IF (JUP .EQ. MPI_PROC_NULL) THEN
         DO J = JEND+1, NYDIM
            DO I = 1, NXDIM
               IJE = (JEND - 1) * NXDIM + I
               IJ  = (J - 1) * NXDIM + I
               DY(IJ) = DY(IJE)
               DYM(IJ) = DYM(IJE)
               COR(IJ) = COR(IJE)
               DEPT(IJ) = DEPT(IJE)
               RDEPV(IJ) = RDEPV(IJE)
               HXT(IJ) = HXT(IJE)
               HXU(IJ) = HXU(IJE)
               HYT(IJ) = HYT(IJE)
               HYU(IJ) = HYU(IJE)
               HXYT(IJ) = HXYT(IJE)
               HXYU(IJ) = HXYU(IJE)
               HYXT(IJ) = HYXT(IJE)
               HYXU(IJ) = HYXU(IJE)
            END DO
         END DO
         DO K = 1, NZDIM
            DO J = JEND+1, NYDIM
               DO I = 1, NXDIM
                  IJE = (JEND - 1) * NXDIM + I
                  IJ  = (J - 1) * NXDIM + I
                  DZ(IJ, K) = DZ(IJE, K)
                  DZM(IJ, K) = DZM(IJE, K)
                  DZV(IJ, K) = DZV(IJE, K)
               END DO
            END DO
         END DO
      END IF
#endif

      IF (JDOWN .EQ. MPI_PROC_NULL) THEN
         DO J = 1, JSTR-1
            DO I = 1, NXDIM
               IJS = (JSTR - 1) * NXDIM + I
               IJ  = (J - 1) * NXDIM + I
               DY(IJ) = DY(IJS)
               DYM(IJ) = DYM(IJS)
               COR(IJ) = COR(IJS)
               DEPT(IJ) = DEPT(IJS)
               RDEPV(IJ) = RDEPV(IJS)
               HXT(IJ) = HXT(IJS)
               HXU(IJ) = HXU(IJS)
               HYT(IJ) = HYT(IJS)
               HYU(IJ) = HYU(IJS)
               HXYT(IJ) = HXYT(IJS)
               HXYU(IJ) = HXYU(IJS)
               HYXT(IJ) = HYXT(IJS)
               HYXU(IJ) = HYXU(IJS)
            END DO
         END DO
         DO K = 1, NZDIM
            DO J = 1, JSTR-1
               DO I = 1, NXDIM
                  IJS = (JSTR - 1) * NXDIM + I
                  IJ  = (J - 1) * NXDIM + I
                  DZ(IJ, K) = DZ(IJS, K)
                  DZM(IJ, K) = DZM(IJS, K)
                  DZV(IJ, K) = DZV(IJS, K)
               END DO
            END DO
         END DO
      END IF

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2DI
      END IF
      CALL SCATTER_2D_INT(NBOT, G2DI)

#ifdef OPT_BBL
      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2DI
      END IF
      CALL SCATTER_2D_INT(NBOTV, G2DI)
#endif

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMSKT, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMFTX, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMFTY, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMFTZ, G3D)

#ifdef OPT_BBL
      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(AMSKTB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(AMSKT0, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(AMSKT1, G2D)
#endif

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMSKV, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMFVX, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMFVY, G3D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G3D
      END IF
      CALL SCATTER_3D(AMFVZ, G3D)

#ifdef OPT_BBL
      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(AMSKVB, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(AMSKV0, G2D)

      IF (MYRANK .EQ. IROOT) THEN
         READ(NFMASK) G2D
      END IF
      CALL SCATTER_2D(AMSKV1, G2D)
#endif

      IF (MYRANK .EQ. IROOT) THEN      
         READ(NFMASK) G3D
      ENDIF
      CALL SCATTER_3D(AMSKB, G3D)

#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M             AMSKT,  AMFTZ, 
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,      0,      0 )
      CALL SHIFT2(
     M             AMSKV,  AMFVZ,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,     -1,     -1 )
C
      CALL SHFT1I(
     M              NBOT,
     D             NXDIM,  NYDIM)
      CALL SHFTINT(
     M              NBOT,
     D             NXDIM,  NYDIM)
C
      IF (JUPE .NE. MPI_PROC_NULL.or.JUPW .NE. MPI_PROC_NULL) THEN
C      DO J = JEND+1, NYDIM
C         DO I = 2, NXDIM
C            IJW = (J - 1) * NXDIM + I  - 1
C            IJ  = (J - 1) * NXDIM + I
C            DO K=1,NZDIM
C               AMFTX(IJ,K) = AMSKT(IJ,K)*AMSKT(IJW,K)
C     &                     * MIN(DZ(IJ,K),DZ(IJW,K))
C             ENDDO
C         ENDDO
C      ENDDO
      DO J = JEND+1, NYDIM
         DO I = 1, NXDIM
            IJS = (J - 2) * NXDIM + I
            IJ  = (J - 1) * NXDIM + I
            DO K=1,NZDIM
                AMFTY(IJ,K) = AMSKT(IJ,K)*AMSKT(IJS,K)
     &                      * MIN(DZ(IJ,K),DZ(IJS,K))
             ENDDO
         ENDDO
      ENDDO
C      DO J = JEND+1, NYDIM-1
C         DO I = 2, NXDIM
C            IJN = J * NXDIM + I
C            IJW = (J - 1) * NXDIM + I  - 1
C            IJ  = (J - 1) * NXDIM + I
C            DO K= 1, KSTR+KZ-1
C                AMFVX(IJ,K) = AMSKT(IJ,K)*AMSKT(IJN,K)
C            ENDDO
C            DO K=KSTR+KZ,NZDIM
C                AMFVX(IJ,K) = AMSKT(IJ,K)*AMSKT(IJN,K)
C     &                      * MIN(DZV(IJ,K),DZV(IJW,K))
C             ENDDO
C         ENDDO
C      ENDDO
C      DO J = JEND+1, NYDIM
C         DO I = 1, NXDIM-1
C            IJS = (J - 2) * NXDIM + I
C            IJE = (J - 1) * NXDIM + I  + 1
C            IJ  = (J - 1) * NXDIM + I
C            DO K= 1, KSTR+KZ-1
C                AMFVY(IJ,K) = AMSKT(IJ,K)*AMSKT(IJE,K)
C            ENDDO
C            DO K= KSTR+KZ,NZDIM
C                AMFVY(IJ,K) = AMSKT(IJ,K)*AMSKT(IJE,K)
C     &                      * MIN(DZV(IJ,K),DZV(IJS,K))
C            ENDDO
C         ENDDO
C      ENDDO
      ENDIF

      CALL SHIFT1(
     M             AMFTX, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,      1,      0)
      CALL SHIFT1(
     M             AMFTY, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,      0,      1)
      CALL SHIFT1(
     M             AMFVX, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,      0,     -1)
      CALL SHIFT1(
     M             AMFVY, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,     -1,      0)
#else
      CALL SHIFT2(
     M             AMSKT,  AMSKV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M             AMFTX,  AMFTY,  AMFTZ,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M             AMFVX,  AMFVY,  AMFVZ,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHFT1I(
     M              NBOT,
     D             NXDIM,  NYDIM)
#endif

#ifdef OPT_BBL
#ifdef OPT_TRIPOLE
      CALL SHIFT3(
     M            AMSKTB, AMSKT0, AMSKT1,
     D             NXDIM,  NYDIM,      1,
     I              1.D0,      0,      0 )
      CALL SHIFT3(
     M            AMSKVB, AMSKV0, AMSKV1,
     D             NXDIM,  NYDIM,      1,
     I              1.D0,     -1,     -1 )
      CALL SHFT1I(
     M             NBOTV,
     D             NXDIM,  NYDIM)
      CALL SHFTINV(
     M             NBOTV,
     D             NXDIM,  NYDIM)
      CALL SHFT1I(
     M             NBOTV,
     D             NXDIM,  NYDIM)
C
      IF (JUPE .NE. MPI_PROC_NULL.or.JUPW .NE. MPI_PROC_NULL) THEN
C      DO J = JEND+1, NYDIM
C         DO I = 2, NXDIM
C            IJS = (J - 2) * NXDIM + I
C            IJ  = (J - 1) * NXDIM + I
C            AMFTX(IJ,KEND) =  AMSKTB(IJ)*AMSKTB(IJW)
C     &                      * MIN(DZ(IJ,KEND),DZ(IJW,KEND))
C         ENDDO
C      ENDDO
      DO J = JEND+1, NYDIM
         DO I = 1, NXDIM
            IJS = (J - 2) * NXDIM + I
            IJ  = (J - 1) * NXDIM + I
            AMFTY(IJ,KEND) =  AMSKTB(IJ)*AMSKTB(IJS)
     &                      * MIN(DZ(IJ,KEND),DZ(IJS,KEND))
         ENDDO
      ENDDO
C      DO J = JEND+1, NYDIM-1
C         DO I = 2, NXDIM
C            IJN = J * NXDIM + I
C            IJW = (J - 1) * NXDIM + I  - 1
C            IJ  = (J - 1) * NXDIM + I
C            AMFVX(IJ,KEND) =  AMSKTB(IJ)*AMSKTB(IJN)
C     &                      * MIN(DZV(IJ,KEND),DZV(IJW,KEND))
C         ENDDO
C      ENDDO
C      DO J = JEND+1, NYDIM
C         DO I = 1, NXDIM-1
C            IJS = (J - 2) * NXDIM + I
C            IJE = (J - 1) * NXDIM + I  + 1
C            IJ  = (J - 1) * NXDIM + I
C            AMFVY(IJ,KEND) =  AMSKTB(IJ)*AMSKTB(IJE)
C     &                      * MIN(DZV(IJ,KEND),DZV(IJS,KEND))
C         ENDDO
C      ENDDO
      ENDIF
      CALL SHIFT1(
     M             AMFTX, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,      1,      0)
      CALL SHIFT1(
     M             AMFTY, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,      0,      1)
      CALL SHIFT1(
     M             AMFVX, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,      0,     -1)
      CALL SHIFT1(
     M             AMFVY, 
     D             NXDIM,  NYDIM,  NZDIM,
     I             1.D0 ,     -1,      0)
#else
      CALL SHFT1I(
     M             NBOTV,
     D             NXDIM,  NYDIM)
      CALL SHIFT3(
     M            AMSKTB, AMSKT0, AMSKT1,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M            AMSKVB, AMSKV0, AMSKV1,
     D             NXDIM,  NYDIM,      1)
#endif
#endif

      DO 200 IJ = IJSTR, IJEND
         IJTSTR = IJ
         IF (AMSKT(IJ, KSTR) .EQ. 1.D0) GO TO 209
  200 CONTINUE
  209 CONTINUE
      DO 210 IJ = IJEND, IJSTR, -1
         IJTEND = IJ
         IF (AMSKT(IJ, KSTR) .EQ. 1.D0) GO TO 219
  210 CONTINUE
  219 CONTINUE

      DO 300 IJ = IJSTR, IJEND
         IJVSTR = IJ
         IF (AMSKV(IJ, KSTR) .EQ. 1.D0) GO TO 309
  300 CONTINUE
  309 CONTINUE
      DO 310 IJ = IJEND, IJSTR, -1
         IJVEND = IJ
         IF (AMSKV(IJ, KSTR) .EQ. 1.D0) GO TO 319
  310 CONTINUE
  319 CONTINUE
#else
      READ(NFMASK) DX
      READ(NFMASK) DY
      READ(NFMASK) DYM
      READ(NFMASK) DZ
      READ(NFMASK) DZM
      READ(NFMASK) DZV
      READ(NFMASK) DZ0
      READ(NFMASK) DS
      READ(NFMASK) DSM
      READ(NFMASK) HIC
      READ(NFMASK) REA
      READ(NFMASK) ZBOT
      READ(NFMASK) COR
      READ(NFMASK) DEPT
      READ(NFMASK) RDEPV
      READ(NFMASK) HXT
      READ(NFMASK) HXU
      READ(NFMASK) HYT
      READ(NFMASK) HYU
      READ(NFMASK) HXYT
      READ(NFMASK) HXYU
      READ(NFMASK) HYXT
      READ(NFMASK) HYXU
      READ(NFMASK) NBOT
#ifdef OPT_BBL
      READ(NFMASK) NBOTV
#endif
      READ(NFMASK) AMSKT
      READ(NFMASK) AMFTX
      READ(NFMASK) AMFTY
      READ(NFMASK) AMFTZ
#ifdef OPT_BBL
      READ(NFMASK) AMSKTB
      READ(NFMASK) AMSKT0
      READ(NFMASK) AMSKT1
#endif
      READ(NFMASK) AMSKV
      READ(NFMASK) AMFVX
      READ(NFMASK) AMFVY
      READ(NFMASK) AMFVZ
#ifdef OPT_BBL
      READ(NFMASK) AMSKVB
      READ(NFMASK) AMSKV0
      READ(NFMASK) AMSKV1
#endif
      READ(NFMASK) AMSKB
#endif

#ifdef OPT_PARALLEL
      IF (MYRANK .EQ. IROOT) THEN
         CALL FILCLS(NFMASK)
      ENDIF
      DEALLOCATE(G2D)
      DEALLOCATE(G3D)
      DEALLOCATE(G2DI)
#else
      CALL FILCLS(NFMASK)
#endif

      RX = 1.D0 / DX
      DO IJ = 1, NXYDIM
         RY(IJ) = 1.D0 / DY(IJ)
         RYM(IJ) = 1.D0 / DYM(IJ)
         RXT(IJ) = 1.D0 / HXT(IJ)
         RXU(IJ) = 1.D0 / HXU(IJ)
         RYT(IJ) = 1.D0 / HYT(IJ)
         RYU(IJ) = 1.D0 / HYU(IJ)
      END DO
      DO K = 1, NZDIM
         RS(K) = 1.D0 / DS(K)
         RSM(K) = 1.D0 / DSM(K)
      END DO

      RETURN
      END
