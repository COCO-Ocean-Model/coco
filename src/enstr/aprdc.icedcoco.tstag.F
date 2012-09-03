      SUBROUTINE PREDCO(
     M                    HX,   UBTX,   VBTX,      W,      R,
     M                    UX,     VX,     TX,
     M                   AMV,    AHV,
     I                  TAUX,   TAUY,   PTOP,
     I                    FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                    TQ,
#endif
     I                    UY,     VY,     TY,
     I                    HY,   UBTY,   VBTY,     AX,
     W                    GX,     GY,     XX,     YY,
     W                   GXX,    GYY,
     W                  UADV,   VADV,   WADV)
      use zocdim
      use zocgrd
      use dvdif
      use tslvt
      use tflxt
      use tovtr
      use fbtav
      use dvlta
      use dwdns
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '02.10.10  H.Hasumi: from MIROC3.1-OMIP
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.01.20  T.Suzuki: bug fix (SHIFT UX/VX after STBBVT)
c     '09.10.06  Y.Komuro: bug fix (skip barotro. part when INIT/FINAL)
c     '10.04.14  M.Kurogi: staggered time stepping
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c     '12.09.04  Y.Komuro: (change surface water flux by T. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8      TY(NXYZDM, NTDIM),     TX(NXYZDM, NTDIM)
      REAL*8      UY(NXYZDM),     UX(NXYZDM)
      REAL*8      VY(NXYZDM),     VX(NXYZDM)
      REAL*8      HY(NXYDIM),     HX(NXYDIM)
      REAL*8    UBTY(NXYDIM),   UBTX(NXYDIM)
      REAL*8    VBTY(NXYDIM),   VBTX(NXYDIM)
      REAL*8      AX(NXYDIM, 0:NIC)
      REAL*8    UADV(NXYZDM),   VADV(NXYZDM),   WADV(NXYZDM, 9)
      REAL*8       W(NXYZDM),      R(NXYZDM)
      REAL*8     AMV(NXYZDM),    AHV(NXYZDM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8      FT(NXYDIM, NTDIM)
      REAL*8   SWABS(NXYDIM),     FS(NXYDIM)
      REAL*8    SSFC(NXYDIM)
      REAL*8    PTOP(NXYDIM)
#ifdef OPT_BODY
      REAL*8      TQ(NXYZDM, NTDIM)
#endif
      REAL*8      GX(NXYZDM),     GY(NXYZDM)
      REAL*8      XX(NXYZDM),     YY(NXYZDM)
      REAL*8     GXX(NXYDIM),    GYY(NXYDIM)
      REAL*8     TMP(NXYZDM, NTDIM)

      REAL*8   UBTAV(NXYDIM),  VBTAV(NXYDIM)
      SAVE UBTAV, VBTAV

      REAL*8     FUX(NXYZDM),    FVX(NXYZDM)
      REAL*8     FUY(NXYZDM),    FVY(NXYZDM)
      REAL*8    FUNE(NXYZDM),   FVNE(NXYZDM)
      REAL*8    FUSE(NXYZDM),   FVSE(NXYZDM)
      REAL*8      HZ(NXYDIM)
      REAL*8    HTMP(NXYDIM),  UBTMP(NXYDIM),  VBTMP(NXYDIM)
      INTEGER ITSPLT,     IJ,    IJK,      N

      REAL*8  HXB(NXYDIM)
      REAL*8      H1(NXYDIM)
      SAVE H1
      LOGICAL OEOF

#ifdef OPT_PARALLEL
      IF (      (MYRANK .GE. IJNODE)
     &    .AND. (.NOT. OINIT) .AND. (.NOT. OFINAL)) RETURN
#endif

      IF(OINIT) THEN
#ifdef OPT_TRIPOLE
         CALL RSTADD(H1, OEOF, NXDIM, NYDIM, 1, 'H1', 'SFC',
     &                                        1.D0,  0,  0 )
#else
         CALL RSTADD(H1, OEOF, NXDIM, NYDIM, 1, 'H1', 'SFC')
#endif
         IF(OEOF) THEN
         DO IJ = 1, NXYDIM
            H1(IJ) = HX(IJ)
         END DO
         END IF
      END IF



      IF (OFINAL) THEN
         CALL FINADD(H1, NXDIM, NYDIM, 1, 'H1', 'SFC')
      END IF




c *** vertical diffusivity and viscosity ***
      CALL CLCSTR('COEFF')
         CALL VDIFF (
     M                  AMV,    AHV,
     I                   UY,     VY,      R,   TAUX,   TAUY,
     I                   TY,     HY)
#ifdef OPT_BBL
         CALL VDIFFB(
     M                  AMV,    AHV)
#endif
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT1(   AMV,
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,     -1,     -1 )
         CALL SHIFT1(   AHV,
     D                NXDIM,  NYDIM,  NZDIM,
     I                 1.D0,      0,      0 )
#else
         CALL SHIFT2(
     M                  AMV,    AHV,
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#else
         CALL STBCAV(
     M                  AMV,    AHV)
#endif
      CALL CLCEND('COEFF')

c *** baroclinic flow ***

      CALL CLCSTR('BRCLI')
#ifdef OPT_BBL
         CALL RMMSKV
#endif
         CALL SRCVEL(
     O                   GX,     GY,     XX,     YY,
     I                   UX,     VX)
         CALL VSCVEL(
     M                   GX,     GY,     XX,     YY,
     I                   UY,     UX,     VY,     VX,
     I                   H1,     HY,
C     I                   HX,     HY,
     I                  AMV,   TAUX,   TAUY)
         CALL ADVVEL(
     O                  FUX,    FUY,   FUNE,   FUSE,
     O                  FVX,    FVY,   FVNE,   FVSE,
     M                   GX,     GY,     XX,     YY,
     I                   UY,     UX,     VY,     VX,
     I                   H1,     HY,
C     I                   HX,     HY,
     I                 UADV,   VADV,   WADV)
#ifdef OPT_BBL
         CALL SRCVLB(
     M                   GX,     GY,     XX,     YY,
     I                   UX,     VX)
         CALL VSCVLB(
     M                   GX,     GY,     XX,     YY,
     I                   UY,     UX,     VY,     VX,
     I                  AMV)
         CALL ADVVLB(
     M                  FUX,    FUY,   FUNE,   FUSE,
     M                  FVX,    FVY,   FVNE,   FVSE,
     M                   GX,     GY,     XX,     YY,
     I                   UY,     UX,     VY,     VX,
     I                 UADV,   VADV,   WADV)
#endif
         CALL TNDUVD(
     O                  GXX,    GYY,
     M                   GX,     GY,
     I                   XX,     YY,
     I                   UY,     VY,     HY,      R,
     I                   UX,     VX,
     I                  FUX,    FUY,   FUNE,   FUSE,
     I                  FVX,    FVY,   FVNE,   FVSE)
#ifdef OPT_BBL
         CALL TNDUVB(
     M                  GXX,    GYY,
     M                   GX,     GY,
     I                   XX,     YY,
     I                   UY,     VY,     TY,    
     I                   UX,     VX,
     I                  FUX,    FUY,   FUNE,   FUSE,
     I                  FVX,    FVY,   FVNE,   FVSE)
         CALL ADMKV1
         CALL STBBGV(
     M                   GX,     GY)
#endif
         CALL VELRDS(
     M                   UX,     VX,
     I                   GX,     GY,    AMV,
     I                   H1)
C     I                   HX)
#ifdef OPT_BBL
         CALL STBBUV(
     M                   UX,     VX)
#endif
      CALL CLCEND('BRCLI')

c *** barotropic flow and surface elevation ***

      CALL CLCSTR('BRTRO')
      IF (OINIT .OR. OFINAL) THEN
         CALL MODGXY(
     M                  GXX,    GYY,
     I                 UBTX,   VBTX)
      ELSE
      IF (IEULER .EQ. 2) THEN
         DO IJ = 1, NXYDIM
            HZ   (IJ) = HX   (IJ)
            HTMP (IJ) = HX   (IJ)
            UBTMP(IJ) = UBTX (IJ)
            VBTMP(IJ) = VBTX (IJ)
            HX   (IJ) = HY   (IJ)
            UBTX (IJ) = UBTY (IJ)
            VBTX (IJ) = VBTY (IJ)
            HY   (IJ) = HTMP (IJ)
            UBTY (IJ) = UBTMP(IJ)
            VBTY (IJ) = VBTMP(IJ)
         END DO
      ELSE
         DO IJ = 1, NXYDIM
            HZ(IJ) = HX(IJ)
         END DO
         CALL MODGXY(
     M                  GXX,    GYY,
     I                 UBTX,   VBTX)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT2(   GXX,    GYY,
     D                NXDIM,  NYDIM,      1,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                  GXX,    GYY,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCGG(
     M                  GXX,    GYY)
#endif
         CALL BTAVST(
     O                UBTAV,  VBTAV)
         DO ITSPLT = 1, NTSS
            DO IJ = 1, NXYDIM
               HTMP(IJ) = HX(IJ)
               UBTMP(IJ) = UBTX(IJ)
               VBTMP(IJ) = VBTX(IJ)
            END DO
            CALL SHALOW(
     M                    HTMP,  UBTMP,  VBTMP,
     I                      HX,   UBTX,   VBTX,
     I                     GXX,    GYY,   PTOP,  FT(1,2))
            CALL BTAVAD(
     M                   UBTAV,  VBTAV,
     I                   UBTMP,  VBTMP)
            CALL SHALOW(
     M                      HX,   UBTX,   VBTX,
     I                    HTMP,  UBTMP,  VBTMP,
     I                     GXX,    GYY,   PTOP,  FT(1,2))
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
            CALL SHIFT2(  UBTX,   VBTX,
     D                   NXDIM,  NYDIM,      1,
     I                   -1.D0,     -1,     -1 )
            CALL SHIFT1(    HX,
     D                   NXDIM,  NYDIM,      1,
     I                    1.D0,      0,      0 )
#else
            CALL SHIFT3(
     M                      HX,   UBTX,   VBTX,
     D                   NXDIM,  NYDIM,      1)
#endif
#else
            CALL STBCSH(
     M                      HX)
            CALL STBCBT(
     M                    UBTX,   VBTX)
#endif

         IF (ITSPLT .EQ. NTSS/2) THEN
            DO IJ = 1, NXYDIM
               H1(IJ) = HX(IJ)
            END DO
         END IF


         END DO
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT2( UBTAV,  VBTAV,
     D                NXDIM,  NYDIM,      1,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                UBTAV,  VBTAV,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCBT(
     M                UBTAV,  VBTAV)
#endif
      END IF
         DO IJ = 1, NXYDIM
            HXB(IJ) = HX(IJ)
         END DO
      END IF
      CALL CLCEND('BRTRO')

c *** velocity for tracer advection ***

      CALL CLCSTR('TDIAG')
#ifdef OPT_BBL
         CALL RMMSKV
         CALL ADMKVB
#endif
         CALL VELTAD(
     O                 UADV,   VADV,
     I                UBTAV,  VBTAV,     UX,     VX)
C     I                UBTAV,  VBTAV,     UY,     VY)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      IF (.NOT.(OINIT .OR. OFINAL)) THEN
         CALL SHIFT2(   UADV,    VADV,
     D                 NXDIM,   NYDIM,  NZDIM,
     I                 -1.D0,      -1,     -1 )
      END IF
#else
         CALL SHIFT2(
     M                 UADV,   VADV,
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#else
         CALL STBCUV(
     M                 UADV,   VADV)
#endif
#ifdef OPT_BBL
         CALL RMMSKT
         CALL ADMKTB
#endif
         CALL WDENST(
     O                 WADV,
     I                 UADV,   VADV,
     I                   HX,     HZ)
      CALL CLCEND('TDIAG')

c *** tracer ***

      CALL CLCSTR('TRACE')
      IF (OINIT .OR. OFINAL) THEN
         CALL FLXTRC(
     O                  TMP,     XX,
     I                   TX,     HX,
     I                   TY,     HZ,
     I                 UADV,   VADV,   WADV,
     I                  AHV)
#ifdef OPT_BBL
         CALL FLXTRB(
     O                  TMP,     XX,
     I                   TX,
     I                   TY,
     I                 UADV,   VADV,   WADV,
     I                  AHV)
#endif
         CALL SLVTRC(
     M                   TX,     HX,
     I                  TMP,     XX,
     I                   FT,  SWABS,     FS,     HZ,   SSFC,
     I                   AX)
         CALL SHDIFF(
     M                   TX,     HX)
         CALL OVTURN(
     O                    R,
     M                   TX,
     I                   HX)
      ELSE
      IF (IEULER .EQ. 2) THEN
         DO N = 1, NTDIM
            DO IJK = 1, NXYZDM
               GX(IJK)    = TX(IJK, N)
               TX(IJK, N) = TY(IJK, N)
               TY(IJK, N) = GX(IJK)
            END DO
         END DO
      ELSE
         CALL FLXTRC(
     O                  TMP,     XX,
     I                   TX,     HX,
     I                   TY,     HZ,
     I                 UADV,   VADV,   WADV,
     I                  AHV)
#ifdef OPT_BBL
         CALL FLXTRB(
     O                  TMP,     XX,
     I                   TX,
     I                   TY,
     I                 UADV,   VADV,   WADV,
     I                  AHV)
         CALL RMMSKT
         CALL ADMKT1
         CALL STBBGT(
     M                  TMP,     XX)
#endif
#ifdef OPT_BODY
         DO N = 1, NTDIM
            DO IJK = 1, NXYZDM
               TMP(IJK, N) = TMP(IJK, N) + TQ(IJK, N)
            END DO
         END DO
#endif
         CALL SLVTRC(
     M                   TX,     HX,
     I                  TMP,     XX,
     I                   FT,  SWABS,     FS,     HZ,   SSFC,
     I                   AX)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
            CALL SHIFT1(     TX,
     D                    NXDIM,   NYDIM, NZTDIM,
     I                     1.D0,       0,      0 )
            CALL SHIFT1(     HX,
     D                    NXDIM,   NYDIM,      1,
     I                     1.D0,       0,      0 )
#else
         CALL SHIFT1(
     M                   TX,
     D                NXDIM,  NYDIM, NZTDIM)
         CALL SHIFT1(
     M                   HX,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCTR(
     M                   TX,      R)
         CALL STBCSH(
     M                   HX)
#endif
         CALL SHDIFF(
     M                   TX,     HX)
         CALL OVTURN(
     O                    R,
     M                   TX,
     I                   HX)
#ifdef OPT_BBL
         CALL STBBTR(
     M                   TX)
#endif
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
        CALL SHIFT1(    R,
     D                NXDIM,   NYDIM,  NZDIM,
     I                 1.D0,       0,      0 )
        CALL SHIFT1(   TX,
     D                NXDIM,   NYDIM, NZTDIM,
     I                 1.D0,       0,      0 )
        CALL SHIFT1(   HX,
     D                NXDIM,   NYDIM,      1,
     I                 1.D0,       0,      0 )
#else
         CALL SHIFT1(
     M                    R,
     D                NXDIM,  NYDIM,  NZDIM)
         CALL SHIFT1(
     M                   TX,  
     D                NXDIM,  NYDIM, NZTDIM)
         CALL SHIFT1(
     M                   HX,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCSH(
     M                   HX)
#endif
         CALL STBCTR(
     M                   TX,      R)
      END IF
      END IF
      CALL CLCEND('TRACE')

c *** velocity for momentum advection ***

      CALL CLCSTR('VDIAG')
         IF (.NOT. (OINIT .OR. OFINAL)) THEN
         DO IJK = 1, NXYZDM
            UADV(IJK) = UX(IJK)
            VADV(IJK) = VX(IJK)
         END DO
         END IF
         CALL VELTAD(
     O                   UX,     VX,
     I                 UBTAV,   VBTAV,   UADV,   VADV)
C     I                 UBTX,   VBTX,   UADV,   VADV)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT2(    UX,     VX,
     D                NXDIM,  NYDIM,  NZDIM,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                   UX,     VX,
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#else
         CALL STBCUV(
     M                   UX,     VX)
#endif
#ifdef OPT_BBL
         CALL RMMSKT
         CALL ADMKTB
#endif
         CALL WDENST(
     O                    W,
     I                   UX,     VX,
     I                   HXB,     HZ)
C     I                   HX,     HZ)
#ifdef OPT_BBL
         CALL ADMSKV
#endif
         CALL VELVAD(
     O                 UADV,   VADV,   WADV,
     I                   UX,     VX,      W)
#ifdef OPT_BBL
         CALL VELVAB(
     M                 WADV,
     I                    W)
         CALL STBBVT(
     M                   UX,     VX)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT2(    UX,     VX,
     D                NXDIM,  NYDIM,  NZDIM,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                   UX,     VX,
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#else
         CALL STBCUV(
     M                   UX,     VX)
#endif 
         CALL RMMSKV
         CALL ADMKV1
#endif
#ifdef OPT_PARALLEL
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
#else
         CALL SHIFT2(
     M                 UADV,   VADV,
     D                NXDIM,  NYDIM,  NZDIM)
         CALL SHIFT3(
     M                 WADV(1, 1),   WADV(1, 2),   WADV(1, 3),
     D                NXDIM,  NYDIM,  NZDIM)
         CALL SHIFT3(
     M                 WADV(1, 4),   WADV(1, 5),   WADV(1, 6),
     D                NXDIM,  NYDIM,  NZDIM)
         CALL SHIFT3(
     M                 WADV(1, 7),   WADV(1, 8),   WADV(1, 9),
     D                NXDIM,  NYDIM,  NZDIM)
#endif
#else
         CALL STBCAD(
     M                 UADV,   VADV,   WADV)
#endif
      CALL CLCEND('VDIAG')

      IF ( (.NOT. OINIT) .AND. (.NOT. OFINAL)) THEN
         CALL CHEKIN(  UBTAV,  'UBTAV',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(  VBTAV,  'VBTAV',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
      END IF

      RETURN
      END
