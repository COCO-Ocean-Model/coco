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
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '99.08.17  H.Hasumi: from CCSR2-MASK
c     '01.02.16  H.Hasumi: for incorporating BBL model
c     '02.05.29  H.Nakano: tracer dimension
c     '02.06.02  H.Hasumi: combine parallel and nonparallel
c     '07.04.23  H.Hasumi
c     '09.02.23  Y.Komuro: call DDENST(diagnose R) instead of OVTURN
c                          bug fix (SHIFT UX/VX after STBBVT)
c     '10.04.14  M.kurogi
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8      TB(NXYZDM, NTDIM)
      REAL*8      UB(NXYZDM),     VB(NXYZDM)
      REAL*8      HB(NXYDIM)
      REAL*8    UBTB(NXYDIM),   VBTB(NXYDIM)
      REAL*8    UADV(NXYZDM),   VADV(NXYZDM),   WADV(NXYZDM)
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
#else
      CALL STBCSH(
     M                HB)
      CALL STBCBT(
     M              UBTB,   VBTB)
      CALL STBCAD(
     M                UB,     VB,      W)
#endif

      DO 10 IJK = 1, NXYZDM
         R(IJK) = 0.0D+0
   10 CONTINUE
      CALL DDENST(
     O                    R,
     I                   TB)
#ifdef OPT_BBL
      CALL STBBTR(
     M                TB)
#endif
#ifdef OPT_PARALLEL
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

         CALL RMMSKT
         CALL ADMKTB
         CALL ADMSKV
#endif
         CALL VELVAD(
     O                 UADV,   VADV,   WADV,
     I                   UB,     VB,      W)
#ifdef OPT_BBL
         CALL STBBVT(
     M                   UB,     VB)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT2(    UB,     VB,
     D                NXDIM,  NYDIM,  NZDIM,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                   UB,     VB,
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#else
         CALL STBCUV(
     M                   UB,     VB)
#endif
         CALL RMMSKT
         CALL ADMKT1
         CALL RMMSKV
         CALL ADMKV1
#endif

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT1(  UADV,
     D             NXDIM,  NYDIM,  NZDIM,
     I             -1.D0,      0,     -1 )
      CALL SHIFT1(  VADV,
     D             NXDIM,  NYDIM,  NZDIM,
     I             -1.D0,     -1,      0 )
      CALL SHIFT1(  WADV,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,     -1,     -1 )
      CALL SHIFT1(   AMV,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,     -1,     -1 )
      CALL SHIFT1(   AHV,
     D             NXDIM,  NYDIM,  NZDIM,
     I              1.D0,      0,      0 )
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
#else
      CALL SHIFT3(
     M              UADV,   VADV,   WADV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT2(
     M               AMV,    AHV,
     D             NXDIM,  NYDIM,  NZDIM)
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1)
#endif
#else
      CALL STBCAD(
     M              UADV,   VADV,   WADV)
      CALL STBCAV(
     M               AMV,    AHV)
      CALL STBCTU(
     M              TAUX,   TAUY)
#endif

      RETURN
      END
