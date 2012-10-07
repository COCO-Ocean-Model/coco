      SUBROUTINE INISET(
     O               UADV,   VADV,   WADV,      R,
     I                 UB,     VB,     TB,
     I                 HB,   UBTB,   VBTB,
     I                  W,    AMV,    AHV,
     I               TAUX,   TAUY,
     I                 FT)
      use zocdim
      use zocgrd
      use zocphy
      use zocmsk
      use bchmk

      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '00.05.30  H.Hasumi: parallelized COCO3
c     '00.12.07  H.Hasumi: combine par and non-par routines
c     '02.05.29  H.Nakano: tracer dimension
c     '07.04.23  H.Hasumi
c     '09.02.23  Y.Komuro: call DDENST(diagnose R) instead of OVTURN
c     '10.04.14  M.kurogi
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8      TB(NXYZDM, NTDIM)
      REAL*8      UB(NXYZDM),     VB(NXYZDM)
      REAL*8      HB(NXYDIM)
      REAL*8    UBTB(NXYDIM),   VBTB(NXYDIM)
      REAL*8    UADV(NXYZDM),   VADV(NXYZDM),   WADV(NXYZDM, 9)
      REAL*8       W(NXYZDM),      R(NXYZDM)
      REAL*8     AMV(NXYZDM),    AHV(NXYZDM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8      FT(NXYDIM, NTDIM)

      REAL*8      UZ(NXYDIM),     VZ(NXYDIM)
c      COMMON /WORK/ UZ, VZ

      INTEGER    IJK, IJ, K

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT1(    HB,
     D             NXDIM,  NYDIM,      1,
     I              1.D0,      0,      0 )
      CALL SHIFT2(  UBTB,   VBTB,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,     -1,     -1 )
      CALL SHIFT2(    UB,     VB,
     D             NXDIM,  NYDIM,  NZDIM,
     I             -1.D0,     -1,     -1 )
      CALL SHIFT1(     W,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,      0,      0 )
#else
      CALL SHIFT3(
     M                HB,   UBTB,   VBTB,
     D             NXDIM,  NYDIM,      1)
      CALL SHIFT3(
     M                UB,     VB,      W,
     D             NXDIM,  NYDIM,  NZDIM)
#endif
      DO 10 IJK = 1, NXYZDM
         R(IJK) = 0.0D+0
   10 CONTINUE
      CALL DDENST(
     O                    R,
     I                   TB)
#ifdef OPT_TRIPOLE
      CALL SHIFT1(     R,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,      0,      0 )
      CALL SHIFT1(    TB,
     D             NXDIM,  NYDIM, NZTDIM,
     I              1.D0,      0,      0 )
#else
      CALL SHIFT1(
     M                 R,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT1(
     M                TB,
     D             NXDIM,  NYDIM, NZTDIM)
#endif
      CALL STBCTR(
     M                TB,      R)
#ifdef OPT_BBL
         CALL RMMSKV
         CALL ADMKVB
         DO K=1, NZDIM
         DO IJ=1,NXYDIM
            IJK= IJ + (K-1)*NXYDIM
            UB(IJK)=UB(IJK)*AMSKV(IJ,K)
            VB(IJK)=VB(IJK)*AMSKV(IJ,K)
         END DO
         END DO
         CALL ADMSKV

         CALL RMMSKT
         CALL ADMKTB
#endif
      CALL VELVAD(
     O              UADV,   VADV,   WADV,
     I                UB,     VB,      W)
#ifdef OPT_BBL
         CALL VELVAB(
     M                 WADV,
     I                    W)
         CALL STBBVT(
     M                   UB,     VB)
#ifdef OPT_TRIPOLE
         CALL SHIFT2(    UB,     VB,
     D                NXDIM,  NYDIM,  NZDIM,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                   UB,     VB,
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#endif
#ifdef OPT_TRIPOLE
         CALL SHIFT2(
     M                 UADV,  VADV,
     D                NXDIM,  NYDIM,  NZDIM,
     I                -1.D0,      0,      0)
         CALL SHIFT3(
     M                 WADV(1, 1),   WADV(1, 2),   WADV(1, 3),
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,     -1,     -1)
         CALL SHIFT3(
     M                 WADV(1, 4),   WADV(1, 5),   WADV(1, 6),
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,     -1,     -1)
         CALL SHIFT3(
     M                 WADV(1, 7),   WADV(1, 8),   WADV(1, 9),
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,     -1,     -1)
         CALL EXCNGW(WADV)

      CALL SHIFT1(   AMV,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,     -1,     -1 )
      CALL SHIFT1(   AHV,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,      0,      0 )
      CALL SHIFT2(  TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I              1.D0,      0,      0 )
#else
      CALL SHIFT2(
     M              UADV,   VADV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M              WADV(1, 1),   WADV(1, 2),   WADV(1, 3),
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M              WADV(1, 4),   WADV(1, 5),   WADV(1, 6),
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT3(
     M              WADV(1, 7),   WADV(1, 8),   WADV(1, 9),
     D             NXDIM,  NYDIM,  NZDIM)

      CALL SHIFT2(
     M               AMV,    AHV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1)
#endif
#else

      CALL STBCSH(
     M                HB)
      CALL STBCBT(
     M              UBTB,   VBTB)
      DO IJK = 1, NXYZDM
         WADV(IJK, 1) = W(IJK)
      END DO
      CALL STBCAD(
     M                UB,     VB,   WADV)
      DO IJK = 1, NXYZDM
         W(IJK) = WADV(IJK, 1)
      END DO

      DO 10 IJK = 1, NXYZDM
         R(IJK) = 0.0D+0
   10 CONTINUE
      CALL DDENST(
     O                    R,
     I                   TB)
      CALL STBCTR(
     M                TB,      R)

#ifdef OPT_BBL
         CALL RMMSKV
         CALL ADMKVB
         DO K=1, NZDIM
         DO IJ=1,NXYDIM
            IJK= IJ + (K-1)*NXYDIM
            UB(IJK)=UB(IJK)*AMSKV(IJ,K)
            VB(IJK)=VB(IJK)*AMSKV(IJ,K)
         END DO
         END DO
         CALL ADMSKV
         CALL RMMSKT
         CALL ADMKTB
#endif
      CALL VELVAD(
     O              UADV,   VADV,   WADV,
     I                UB,     VB,      W)

#ifdef OPT_BBL
         CALL VELVAB(
     M                 WADV,
     I                    W)
         CALL STBBVT(
     M                   UB,     VB)
         CALL STBCUV(
     M                   UB,     VB)
#endif

      CALL STBCAD(
     M              UADV,   VADV,   WADV)

      CALL STBCAV(
     M               AMV,    AHV)
      CALL STBCTU(
     M              TAUX,   TAUY)
#endif

      RETURN
      END
