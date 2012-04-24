      program trans
      implicit none

      integer NX, NY, NZ, KZ
      parameter(NX = 128, NY = 120, NZ = 40, KZ = 5)
      integer ISTR, JSTR, KSTR
      parameter(ISTR = 3, JSTR = 3, KSTR = 2)
      INTEGER   IEND, JEND, KEND
      PARAMETER(IEND = ISTR+NX-1, JEND = JSTR+NY-1, KEND = KSTR+NZ-1)
      integer NXDIM, NYDIM, NZDIM
      parameter(NXDIM = NX+4, NYDIM = NY+4, NZDIM = NZ+2)

      integer NXG, NYG
      parameter(NXG = 320, NYG = 160)
      real*4 RMISS
      parameter(RMISS = -99.0)

      REAL*8 ALD(-1:NX+2, -1:NY+2, 0:NZ+1)
      REAL*8 BLD(-1:NX+2, -1:NY+2, 0:NZ+1)
      real*8 AGEO(-1:NX+2, -1:NY+2), BGEO(-1:NX+2, -1:NY+2)
      real*8 DGEO(-1:NX+2, -1:NY+2)
      integer NBOT(-1:NX+2, -1:NY+2), NBOTV(-1:NX+2, -1:NY+2)
      integer IGEO(-1:NX+2, -1:NY+2)
      real*8 DX, DY(-1:NY+2), DZ(0:NZ+1)
      real*8 DZZ(-1:NX+2, -1:NY+2, 0:NZ+1)
      real*8 DZV(-1:NX+2, -1:NY+2, 0:NZ+1)
      real*8 ALONW, ALONE, ALATS, ALATN
      real*8 PLON, PLAT, QLON, QLAT, ELON, ELAT
      integer IBBL, IMERC, ITHIRD
      integer INODES, JNODES
      real*8 DYM(-1:NY+2)
      real*8 HXT(-1:NX+2, -1:NY+2), HXU(-1:NX+2, -1:NY+2)
      real*8 HYT(-1:NX+2, -1:NY+2), HYU(-1:NX+2, -1:NY+2)

      real*4 datin(0:NXG+1, 0:NYG+1)
      real*4 OMASK(0:NXG+1, 0:NYG+1)
      real*8 GLONT(0:NXG+1), GLONV(0:NXG)
      real*8 GLATT(0:NYG+1), GLATV(0:NYG)
      real*8 GCST(0:NYG+1)
      integer GCOAST(NXG, NYG)
      character GEOIN*256, GRIDIN*256, FILEIN*256

      real*8 datou(NX, NY)
      real*8 RLONT(NX), RLONV(0:NX)
      real*8 RLATT(NY), RLATV(0:NY)
      integer RCOAST(NX, NY)
      character GEOOU*256, FILEOU*256

      character chead(64)*16, cdummy*1
      integer idate(6)
      integer i, j
      real*8 cst

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      data chead / 64*'                ' /
      data GEOIN / 'Geometry.Roske' /
      data GRIDIN / 'Grid.Roske' /
      data FILEIN / 'runoff.mon' /
      data GEOOU / 'GEO.std2_bbl' /
      data FILEOU / 'ROFF' /

      PI = ATAN(1.D0) * 4.D0
      RAD = PI / 1.8D2
      RADIUS = 6370.D5
      OMEGA = 2.D0 * PI / 8.64D4
      CUNITI = (0.D0, 1.D0)

      open(unit=11, file=GEOIN, status='OLD',
     &     form='FORMATTED', access='SEQUENTIAL')
      do j = NYG, 1, -1
         read(11, *) (OMASK(i, j), i = 1, NXG)
      end do
      close(unit=11)
      do j = 1, NYG
         OMASK(0, j) = OMASK(NXG, j)
         OMASK(NXG+1, j) = OMASK(1, j)
      end do
      do i = 0, NXG+1
         OMASK(i, 0) = OMASK(i, 1)
         OMASK(i, NYG+1) = OMASK(i, NYG)
      end do

      open(unit=11, file=GRIDIN, status='OLD',
     &     form='FORMATTED', access='SEQUENTIAL')
      do i = 1, NXG
         read(11, *) GLONT(i)
         GLONT(i) = GLONT(i) * RAD
      end do
      do i = 0, NXG
         read(11, *) GLONV(i)
         GLONV(i) = GLONV(i) * RAD
      end do
      do j = 1, NYG
         read(11, *) GLATT(j)
         GLATT(j) = GLATT(j) * RAD
      end do
      do j = 0, NYG
         read(11, *) GLATV(j)
         GLATV(j) = GLATV(j) * RAD
      end do
      GLONT(0) = GLONT(1) - (GLONT(2) - GLONT(1))
      GLONT(NXG+1) = GLONT(NXG) + (GLONT(NXG) - GLONT(NXG-1))
      GLATT(0) = GLATT(1) - (GLATT(2) - GLATT(1))
      GLATT(NYG+1) = GLATT(NYG) + (GLATT(NYG) - GLATT(NYG-1))
      do j = 1, NYG
c         GCST(j) = cos(GLATT(j))
         GCST(j) = (cos(GLATV(j-1)) + cos(GLATV(j))) * 0.5D0
      end do
      GCST(0) = GCST(1)
      GCST(NYG+1) = GCST(NYG)

      do j = 1, NYG
         do i = 1, NXG
            if (OMASK(i, j) .gt. 0.) then
               cst = OMASK(i+1, j) * OMASK(i-1, j) *
     &               OMASK(i, j+1) * OMASK(i, j-1) *
     &               OMASK(i-1, j-1) * OMASK(i-1, j+1) *
     &               OMASK(i+1, j-1) * OMASK(i+1, j+1)
               if (cst .eq. 0.) then
                  GCOAST(i, j) = 1
               else
                  GCOAST(i, j) = 0
               end if
            else
               GCOAST(i, j) = 0
            end if
         end do
      end do

      CALL GETGEO(
     O               ALD,    BLD,   DGEO,   NBOT,
     O             IMERC, ITHIRD,   IBBL,
     O                DX,     DY,     DZ,
     O             ALONW,  ALONE,  ALATS,  ALATN,
     O              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     O            INODES, JNODES,
     D                NX,     NY,     NZ,     KZ,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     O              AGEO,   BGEO,   IGEO,
     I             GEOOU)

      CALL GRDPRM(
     O               DZZ,    DZV,
     I               ALD,    BLD,   DGEO,   NBOT,  NBOTV,
     I                DX,     DY,     DZ,
     I             IMERC, ITHIRD,   IBBL,
     I             ALONW,  ALONE,  ALATS,  ALATN,
     I              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     D                NX,     NY,     NZ,     KZ,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     O               DYM,
     O               HXT,    HXU,    HYT,    HYU)

      do j = 1, NY
         do i = 1, NX
            if (AGEO(i, j) .gt. 0.) then
               cst = AGEO(i+1, j) * AGEO(i-1, j) *
     &               AGEO(i, j+1) * AGEO(i, j-1) *
     &               AGEO(i-1, j-1) * AGEO(i-1, j+1) *
     &               AGEO(i+1, j-1) * AGEO(i+1, j+1)
               if (cst .eq. 0.) then
                  RCOAST(i, j) = 1
               else
                  RCOAST(i, j) = 0
               end if
            else
               RCOAST(i, j) = 0
            end if
         end do
      end do

      RLONV(0) = ALONW
      do i = 1, NX
         RLONV(i) = RLONV(i-1) + DX
      end do
      RLONT(1) = ALONW + DX * 0.5D0
      do i = 2, NX
         RLONT(i) = RLONT(i-1) + DX
      end do
      RLATV(0) = ALATS
      do j = 1, NY
         RLATV(j) = RLATV(j-1) + DY(j)
      end do
      RLATT(1) = ALATS + DY(1) * 0.5D0
      do j = 2, NY
         RLATT(j) = (RLATV(j) + RLATV(j-1)) * 0.5D0
      end do

      open(unit=11, status='OLD', file=FILEIN,
     &     form='UNFORMATTED', access='SEQUENTIAL')
      open(unit=12, file=FILEOU,
     &     form='UNFORMATTED', access='SEQUENTIAL')

 10   continue

      read(11, end=999) chead
      write(*, *) chead(27)
      read(chead(27), '(I4.4,2I2.2,A1,3I2.2,A1)')
     &    idate(1), idate(2), idate(3), cdummy,
     &    idate(4), idate(5), idate(6)
      read(11) ((datin(i, j), i = 1, NXG), j = 1, NYG)

      call roff(
     o        datou,
     i        datin,
     i        RLONT, RLONV, RLATT, RLATV,
     i        HXT, HYT, RCOAST,
     i        GLONT, GLONV, GLATT, GLATV,
     i        GCST, GCOAST,
     c        ALONW, ALONE, ALATS, ALATN,
     d        NX, NY, NXG, NYG)

      chead(38) = 'UR8'
      write(chead(50), '(I6.6,5I2.2)') idate
      write(chead(64), '(I16)') NX*NY
      write(12) chead
      write(12) datou

      go to 10
 999  continue

      end

c **********************************************************************

      subroutine roff(
     o        datou,
     i        datin,
     i        RLONT, RLONV, RLATT, RLATV,
     i        HXT, HYT, RCOAST,
     i        GLONT, GLONV, GLATT, GLATV,
     i        GCST, GCOAST,
     c        ALONW, ALONE, ALATS, ALATN,
     d        NX, NY, NXG, NYG)
      implicit none

      integer NX, NY, NXG, NYG
      real*4  datin(0:NXG+1, 0:NYG+1)
      real*8  datou(NX, NY)
      real*8  RLONT(NX), RLONV(0:NX)
      real*8  RLATT(NY), RLATV(0:NY)
      real*8  HXT(-1:NX+2, -1:NY+2), HYT(-1:NX+2, -1:NY+2)
      integer RCOAST(NX, NY)
      real*8  GLONT(0:NXG+1), GLONV(0:NXG)
      real*8  GLATT(0:NYG+1), GLATV(0:NYG)
      real*8  GCST(0:NYG+1)
      integer GCOAST(NXG, NYG)
      real*8  ALONW, ALONE, ALATS, ALATN

      real*8 glon, glat, cosang, cosangmax
      integer i, j, ig, jg, im, jm

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      do j = 1, NY
         do i = 1, NX
            datou(i, j) = 0.D0
         end do
      end do

      do jg = 1, NYG
         do ig = 1, NXG
            if (      (GCOAST(ig, jg) .eq. 1)
     &          .and. (datin(ig, jg) .gt. 0.)) then
               cosangmax = -1.D0
               do j = 1, NY
                  do i = 1, NX
                     if (RCOAST(i, j) .eq. 1) then
                        call RTOG(glon, glat, RLONT(i), RLATT(j))
                        cosang = cos(glat) * cos(glon) *
     &                           cos(GLATT(jg)) * cos(GLONT(ig))
     &                         + cos(glat) * sin(glon) *
     &                           cos(GLATT(jg)) * sin(GLONT(ig))
     &                         + sin(glat) *
     &                           sin(GLATT(jg))
                        if (cosang .gt. cosangmax) then
                           im = i
                           jm = j
                           cosangmax = cosang
                        end if
                     end if
                  end do
               end do
               datou(im, jm) = datou(im, jm)
     &                       + datin(ig, jg) *
     &                         (GLONV(ig) - GLONV(ig-1)) * GCST(jg) *
     &                         (GLATV(jg) - GLATV(jg-1))
     &                          / (RLONV(im) - RLONV(im-1))
     &                          / (RLATV(jm) - RLATV(jm-1))
     &                          / HXT(im, jm) / HYT(im, jm)
            end if
         end do
      end do

      return
      end

c **********************************************************************

      SUBROUTINE GRDPRM(
     O               DZZ,    DZV,
     I               ALD,    BLD,   DGEO,   NBOT,  NBOTV,
     I                DX,     DY,     DZ,
     I             IMERC, ITHIRD,   IBBL,
     I             ALONW,  ALONE,  ALATS,  ALATN,
     I              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     D                NX,     NY,     NZ,     KZ,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     O               DYM,
     O               HXT,    HXU,    HYT,    HYU)
      IMPLICIT NONE

      INTEGER     NX,     NY,     NZ,     KZ
      INTEGER  NXDIM,  NYDIM,  NZDIM
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND
      INTEGER  IMERC, ITHIRD,   IBBL
      REAL*8      DX,     DY(NYDIM),     DZ(NZDIM)
      REAL*8   ALONW,  ALONE,  ALATS,  ALATN
      REAL*8    PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT
      REAL*8     DYM(NYDIM)
      REAL*8     HXT(NXDIM, NYDIM),    HXU(NXDIM, NYDIM)
      REAL*8     HYT(NXDIM, NYDIM),    HYU(NXDIM, NYDIM)
      REAL*8     DZZ(NXDIM, NYDIM, NZDIM)
      REAL*8     DZV(NXDIM, NYDIM, NZDIM)
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8    DGEO(NXDIM, NYDIM)
      INTEGER   NBOT(NXDIM, NYDIM),  NBOTV(NXDIM, NYDIM)

      REAL*8      CX,     CY,     CZ
      COMPLEX*16    ZPT,    ZPV,    ZQT,    ZQV,     DF
      COMPLEX*16     Z1,     Z2,     Z3,     W3
      REAL*8   RLONT,  RLATT,  RLONV,  RLATV,   RANG
      REAL*8      YS,     YN, DYMERC,     Y1,     Y2

      REAL*8  FLGBBL,    DEP
      INTEGER      I,      J,      K
      REAL*8    RLON,   RLAT,   GLON,   GLAT
      INTEGER     IV,     JV,     IT,     JT

      SAVE

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      REAL*8  ZTOLON, ZTOLAT
      REAL*8   FMTRC
      COMPLEX*16   CONF,  DCONF, CONFINV,  STOZ
      COMPLEX*16      Z,     ZP,     ZQ
      COMPLEX*16     ZA,     ZB,     ZC,     WC
      REAL*8    ALON,   ALAT
      STOZ(ALON, ALAT) = TAN(PI * 0.25D0 - ALAT * 0.5D0) *
     &                   EXP(ALON * CUNITI)
      ZTOLON(Z) = ATAN2(DIMAG(Z), DREAL(Z))
      ZTOLAT(Z) = PI * 0.5D0 - 2.D0 * ATAN(ABS(Z))
      CONF(Z, ZA, ZB, ZC, WC) = 
     &    (- ZB * (ZC - ZA) * Z / WC + ZA * (ZC - ZB))
     &    / (- (ZC - ZA) * Z / WC + ZC - ZB)
      DCONF(Z, ZA, ZB, ZC, WC) =
     &    (ZC - ZA) * (ZC - ZB) * (ZA - ZB) / WC
     &    / (ZC - ZB - (ZC - ZA) * Z / WC)
     &    / (ZC - ZB - (ZC - ZA) * Z / WC)
      CONFINV(Z, ZA, ZB, ZC, WC) =
     &    WC * (Z - ZA) * (ZC - ZB) / (Z - ZB) / (ZC - ZA)
c      FMTRC(ZP, ZQ, ZA, ZB, ZC, WC) =
c     &     (1.D0 + ABS(ZP) * ABS(ZP))
c     &     / (1.D0 + ABS(ZQ) * ABS(ZQ))
c     &     / ABS(DCONF(ZQ, ZA, ZB, ZC, WC))
      FMTRC(ZP, ZQ, ZA, ZB, ZC, WC) =
     &     (1.D0 + ABS(ZQ) * ABS(ZQ)) *
     &     ABS(DCONF(ZQ, ZA, ZB, ZC, WC))
     &     / (1.D0 + ABS(ZP) * ABS(ZP))

      ALONW = ALONW * RAD
      ALONE = ALONE * RAD
      ALATS = ALATS * RAD
      ALATN = ALATN * RAD
      PLON = PLON * RAD
      PLAT = PLAT * RAD
      QLON = QLON * RAD
      QLAT = QLAT * RAD
      ELON = ELON * RAD
      ELAT = ELAT * RAD

      DX = (ALONE - ALONW) / DBLE(NX)

       IF (IMERC .EQ. 2) THEN
         DO J = JSTR, JEND
            DY(J) = DY(J) * RAD
         END DO
      ELSE IF (IMERC .EQ. 1) THEN
         YS = LOG(TAN(ALATS * 0.5D0 + PI * 0.25D0))
         YN = LOG(TAN(ALATN * 0.5D0 + PI * 0.25D0))
         DYMERC = (YN - YS) / DBLE(NY)
         DO J = JSTR, JEND
            Y1 = YS + DYMERC * (J - JSTR)
            Y2 = YS + DYMERC * (J - JSTR + 1)
            DY(J) = 2.D0 * (ATAN(EXP(Y2)) - ATAN(EXP(Y1)))
         END DO
      ELSE
         DO J = JSTR, JEND
            DY(J) = (ALATN - ALATS) / DBLE(NY)
         END DO
      END IF
      DO J = 1, JSTR-1
         DY(J) = DY(JSTR)
      END DO
      DO J = JEND+1, NYDIM
         DY(J) = DY(JEND)
      END DO
      DO J = 1, NYDIM-1
         DYM(J) = 0.5D0 * (DY(J) + DY(J+1))
      END DO
      DYM(NYDIM) = 0.5D0 * DY(NYDIM)

      DO 200 K = 1, NZDIM
         DZ(K) = DZ(K) * 1.D2
  200 CONTINUE

      DO K = 1, NZDIM
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               DZZ (I, J, K) = DZ (K)
            END DO
         END DO
      END DO
      DO J = 1, NYDIM
         DO I = 1, NXDIM
            K = NBOT(I, J) + KSTR - 1
            IF (K .GE. KSTR) THEN
               DZZ (I, J, K) = DZ(K) * DGEO(I, J)
            END IF
         END DO
      END DO
      IF (IBBL .NE. 0) THEN
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               DO K = NBOT(I, J)+KSTR, KEND-1
                  DZZ(I, J, K) = DZ(KEND)
               END DO
            END DO
         END DO
      END IF

      DO K = 1, NZDIM
         DO J = JSTR, JEND
            DO I = ISTR, IEND
               DZV(I, J, K) = MIN(DZZ(I  , J  , K),
     &                            DZZ(I+1, J  , K),
     &                            DZZ(I  , J+1, K),
     &                            DZZ(I+1, J+1, K))
               IF (IBBL .NE. 0) THEN
                  FLGBBL = (1.D0 - BLD(I  , J  , K)) *
     &                     (1.D0 - BLD(I+1, J  , K)) *
     &                     (1.D0 - BLD(I  , J+1, K)) *
     &                     (1.D0 - BLD(I+1, J+1, K))
                  IF (FLGBBL .EQ. 0.D0) THEN
                     DZV(I, J, K) = MIN(DZZ(I  , J  , KEND),
     &                                  DZZ(I+1, J  , KEND),
     &                                  DZZ(I  , J+1, KEND),
     &                                  DZZ(I+1, J+1, KEND))
                  END IF
               END IF
            END DO
            DO I = 1, ISTR-1
               DZV(I, J, K) = DZV(IEND-ISTR+1+I, J, K)
            END DO
            DO I = IEND+1, NXDIM
               DZV(I, J, K) = DZV(ISTR-1+I-IEND, J, K)
            END DO
         END DO
         DO J = 1, JSTR-1
            DO I = 1, NXDIM
               DZV(I, J, K) = DZV(I, JSTR, K)
            END DO
         END DO
         DO J = JEND+1, NYDIM
            DO I = 1, NXDIM
               DZV(I, J, K) = DZV(I, JEND, K)
            END DO
         END DO
      END DO

      IF (ITHIRD .EQ. 0) THEN
         W3 = (-1.D0, 0.D0)
         CX = COS(PLON) * COS(PLAT) + COS(QLON) * COS(QLAT)
         CY = SIN(PLON) * COS(PLAT) + SIN(QLON) * COS(QLAT)
         CZ = SIN(PLAT) + SIN(QLAT)
         ELON = ATAN2(CY, CX)
         CX = SQRT(CX * CX + CY * CY)
         ELAT = ATAN2(CZ, CX)
      ELSE
         W3 = EXP(PLON * CUNITI)
      END IF
      Z1 = TAN(PI * 0.25D0 - PLAT * 0.5D0) * EXP(PLON * CUNITI)
      Z2 = TAN(PI * 0.25D0 - QLAT * 0.5D0) * EXP(QLON * CUNITI)
      Z3 = TAN(PI * 0.25D0 - ELAT * 0.5D0) * EXP(ELON * CUNITI)

      RLATV = ALATS
      DO J = JSTR, JEND
         RLATT = RLATV + DY(J) * 0.5D0
         RLATV = RLATV + DY(J)
ccc         DO I = ISTR, IEND
         DO I = 1, NXDIM
            RLONT = ALONW + DX * (DBLE(I - ISTR) + 0.5D0)
            RLONV = ALONW + DX * DBLE(I - ISTR + 1)
            ZQT = STOZ(RLONT, RLATT)
            ZPT = CONF(ZQT, Z1, Z2, Z3, W3)
            ZQV = STOZ(RLONV, RLATV)
            ZPV = CONF(ZQV, Z1, Z2, Z3, W3)
            HYT(I, J) = FMTRC(ZPT, ZQT, Z1, Z2, Z3, W3)
            HXT(I, J) = HYT(I, J) * COS(RLATT)
            HYU(I, J) = FMTRC(ZPV, ZQV, Z1, Z2, Z3, W3)
            HXU(I, J) = HYU(I, J) * COS(RLATV)
         END DO
      END DO

      DO J = 1, JSTR-1
         DO I = 1, NXDIM
            HXT(I, J) = HXT(I, JSTR)
            HXU(I, J) = HXU(I, JSTR)
            HYT(I, J) = HYT(I, JSTR)
            HYU(I, J) = HYU(I, JSTR)
         END DO
      END DO
      DO J = JEND+1, NYDIM
         DO I = 1, NXDIM
            HXT(I, J) = HXT(I, JEND)
            HXU(I, J) = HXU(I, JEND)
            HYT(I, J) = HYT(I, JEND)
            HYU(I, J) = HYU(I, JEND)
         END DO
      END DO

      DO J = 1, NYDIM
         DO I = 1, NXDIM
            NBOTV(I, J) = 0
         END DO
      END DO
      DO J = JSTR-1, JEND+1
         DO I = ISTR-1, IEND+1
            IF (      (NBOT(I  , J  ) .GT. 0)
     &          .AND. (NBOT(I+1, J  ) .GT. 0)
     &          .AND. (NBOT(I  , J+1) .GT. 0)
     &          .AND. (NBOT(I+1, J+1) .GT. 0)) THEN
               NBOTV(I, J) = MIN(NBOT(I, J), NBOT(I+1, J),
     &                           NBOT(I, J+1), NBOT(I+1, J+1))
               DEP = 0.D0
               DO K = KSTR, NBOTV(I, J)+KSTR-1
                  DEP = DEP + DZV(I, J, K)
               END DO
            END IF
         END DO
      END DO

      RETURN

c ----------------------------------------------------------------------

      ENTRY VROT(
     O             RANG,
     I             RLON,   RLAT)

      ZQT = STOZ(RLON, RLAT)
      ZPT = CONF(ZQT, Z1, Z2, Z3, W3)
      DF  = DCONF(ZQT, Z1, Z2, Z3, W3)
      RANG = ATAN2(DIMAG(ZPT), DREAL(ZPT))
     &     - ATAN2(DIMAG(ZQT), DREAL(ZQT))
     &     - ATAN2(DIMAG(DF) , DREAL(DF) )
      RANG = - RANG

      RETURN

c ----------------------------------------------------------------------

      ENTRY RTOG(
     O             GLON,   GLAT,
     I             RLON,   RLAT)

      ZQT = STOZ(RLON, RLAT)
      ZPT = CONF(ZQT, Z1, Z2, Z3, W3)
      GLON = ZTOLON(ZPT)
      GLAT = ZTOLAT(ZPT)

      RETURN

c ----------------------------------------------------------------------

      ENTRY GTOR(
     O             RLON,   RLAT,
     I             GLON,   GLAT)

      ZPT = STOZ(GLON, GLAT)
      ZQT = CONFINV(ZPT, Z1, Z2, Z3, W3)
      RLON = ZTOLON(ZQT)
      RLAT = ZTOLAT(ZQT)

      RETURN

c ----------------------------------------------------------------------

      ENTRY IJVTOR(
     O               RLON,   RLAT,
     I                 IV,     JV,
     I              ALONW,  ALATS,     DX,     DY,
     D                 NY,  NYDIM,   JSTR)

      RLON = ALONW + DX * DBLE(IV)
      RLAT = ALATS
      do j = 1, JV
         RLAT = RLAT + DY(j+JSTR-1)
      end do

      RETURN

c ----------------------------------------------------------------------

      ENTRY IJTTOR(
     O               RLON,   RLAT,
     I                 IT,     JT,
     I              ALONW,  ALATS,     DX,     DY,
     D                 NY,  NYDIM,   JSTR)

      RLON = ALONW + DX * (DBLE(IT-1) + 0.5D0)
      RLAT = ALATS
      do j = 1, JT-1
         RLAT = RLAT + DY(j+JSTR-1)
      end do
      RLAT = RLAT + DY(JT+JSTR-1) * 0.5D0

      RETURN

      END

c **********************************************************************

      SUBROUTINE GETGEO(
     O               ALD,    BLD,   DGEO,   NBOT,
     O             IMERC, ITHIRD,   IBBL,
     O                DX,     DY,     DZ,
     O             ALONW,  ALONE,  ALATS,  ALATN,
     O              PLON,   PLAT,   QLON,   QLAT,   ELON,   ELAT,
     O            INODES, JNODES,
     D                NX,     NY,     NZ,     KZ,
     D             NXDIM,  NYDIM,  NZDIM,
     C              ISTR,   JSTR,   KSTR,   IEND,   JEND,   KEND,
     O              AGEO,   BGEO,   IGEO,
     I            GEOFIL)

      INTEGER   IGFILE
      PARAMETER(IGFILE = 11)
      INTEGER   NROUND
      PARAMETER(NROUND = 1000)

      INTEGER     NX,     NY,     NZ,     KZ
      INTEGER  NXDIM,  NYDIM,  NZDIM
      REAL*8     ALD(NXDIM, NYDIM, NZDIM)
      REAL*8     BLD(NXDIM, NYDIM, NZDIM)
      REAL*8    AGEO(NXDIM, NYDIM), BGEO(NXDIM, NYDIM)
      REAL*8    DGEO(NXDIM, NYDIM)
      INTEGER   IGEO(NXDIM, NYDIM), NBOT(NXDIM, NYDIM)
      REAL*8      DX,     DY(NYDIM),     DZ(NZDIM)
      REAL*8   ALONW,  ALONE,  ALATS,  ALATN
      REAL*8    PLON,   PLAT,   QLON,   QLAT
      REAL*8    ELON,   ELAT
      INTEGER INODES, JNODES
      INTEGER  IMERC, ITHIRD,   IBBL
      INTEGER   ISTR,   JSTR,   KSTR
      INTEGER   IEND,   JEND,   KEND
      CHARACTER GEOFIL*(*)

      INTEGER     LX,     LY,     LZ,     JZ
      INTEGER      I,      J,      K

      IBBL = 0

      OPEN(UNIT=IGFILE, FILE=GEOFIL,
     &     ACCESS='SEQUENTIAL', FORM='FORMATTED')

      READ(IGFILE, *) LX
      READ(IGFILE, *) LY
      READ(IGFILE, *) LZ
      READ(IGFILE, *) JZ
      IF (     (LX .NE. NX) .OR. (LY .NE. NY)
     &    .OR. (LZ .NE. NZ) .OR. (JZ .NE. KZ)) THEN
         WRITE(*, *) ' #### BAD GEOMETRY FILE #### '
         WRITE(*, *) ' #### DIMENSION ERROR   #### '
         STOP
      END IF
      READ(IGFILE, *) INODES, JNODES
      READ(IGFILE, *) ALONW, ALONE
      READ(IGFILE, *) ALATS, ALATN
      READ(IGFILE, *) IMERC
      IF (IMERC .EQ. 2) THEN
         READ(IGFILE, *) (DY(J), J = JSTR, JEND)
      END IF
      READ(IGFILE, *) (DZ(K), K = KSTR, KEND)
      READ(IGFILE, *) PLON, PLAT
      READ(IGFILE, *) QLON, QLAT
      READ(IGFILE, *) ITHIRD
      IF (ITHIRD .EQ. 1) THEN
         READ(IGFILE, *) ELON, ELAT
      END IF

      DO J = JEND, JSTR, -1
         READ(IGFILE, *) (AGEO(I, J), I = ISTR, IEND)
      END DO
      READ(IGFILE, *, END=999) IBBL
      IF (IBBL .NE. 0) THEN
         DO J = JEND, JSTR, -1
            READ(IGFILE, *) (BGEO(I, J), I = ISTR, IEND)
         END DO
      END IF
 999  CONTINUE
      CLOSE(UNIT=IGFILE)
      DO J = JSTR, JEND
         DO I = ISTR, IEND
            IGEO(I, J) = NINT(AGEO(I, J) * DBLE(NROUND)) / NROUND
            DGEO(I, J) = DBLE(  NINT(AGEO(I, J) * DBLE(NROUND))
     &                        - IGEO(I, J) * NROUND)
     &                   / DBLE(NROUND)
         END DO
      END DO
      DO J = JSTR, JEND
         DO I = ISTR, IEND
            IF ((IGEO(I, J) .GT. 0) .AND. (IGEO(I, J) .LT. KZ)) THEN
               WRITE(*, *) ' #### BAD GEOMETRY FILE         #### '
               WRITE(*, *) ' #### BOTTOM LEVEL INCONSISTENT #### '
               STOP
            END IF
         END DO
      END DO

      DO 110 K = 1, KSTR-1
         DZ(K) = DZ(KSTR)
  110 CONTINUE
      DO 120 K = KEND+1, NZDIM
         DZ(K) = DZ(KEND)
  120 CONTINUE

      DO 200 J = JSTR, JEND
         DO 210 I = 1, ISTR-1
            IGEO(I, J) = IGEO(IEND-ISTR+1+I, J)
            DGEO(I, J) = DGEO(IEND-ISTR+1+I, J)
            BGEO(I, J) = BGEO(IEND-ISTR+1+I, J)
  210    CONTINUE
         DO 220 I = IEND+1, NXDIM
            IGEO(I, J) = IGEO(ISTR-1+I-IEND, J)
            DGEO(I, J) = DGEO(ISTR-1+I-IEND, J)
            BGEO(I, J) = BGEO(ISTR-1+I-IEND, J)
  220    CONTINUE
  200 CONTINUE
      DO 230 J = 1, JSTR-1
         DO 240 I = 1, NXDIM
            IGEO(I, J) = 0
            DGEO(I, J) = 0.D0
            BGEO(I, J) = 0.D0
  240    CONTINUE
  230 CONTINUE
      DO 250 J = JEND+1, NYDIM
         DO 260 I = 1, NXDIM
            IGEO(I, J) = 0
            DGEO(I, J) = 0.D0
            BGEO(I, J) = 0.D0
  260    CONTINUE
  250 CONTINUE

      DO 300 K = KSTR, KEND
         DO 300 J = 1, NYDIM
            DO 300 I = 1, NXDIM
               IF (IGEO(I, J) .GE. K-KSTR+1) THEN
                  ALD(I, J, K) = 1.D0
               ELSE IF (      (IGEO(I, J) .EQ. K-KSTR)
     &                  .AND. (DGEO(I, J) .GT. 0.D0)) THEN
                  ALD(I, J, K) = 1.D0
               ELSE
                  ALD(I, J, K) = 0.D0
               END IF
  300 CONTINUE
      DO 310 K = 1, KSTR-1
         DO 310 J = 1, NYDIM
            DO 310 I = 1, NXDIM
               ALD(I, J, K) = ALD(I, J, KSTR)
  310 CONTINUE
      DO 320 K = KEND+1, NZDIM
         DO 320 J = 1, NYDIM
            DO 320 I = 1, NXDIM
               ALD(I, J, K) = 0.D0
  320 CONTINUE

      DO 400 J = 1, NYDIM
         DO 400 I = 1, NXDIM
            IF (IGEO(I, J) .GE. 0) THEN
               IF (DGEO(I, J) .EQ. 0.D0) THEN
                  NBOT(I, J) = IGEO(I, J)
               ELSE
                  NBOT(I, J) = IGEO(I, J) + 1
               END IF
            ELSE
               NBOT(I, J) = 0
            END IF
  400 CONTINUE

      DO J = 1, NYDIM
         DO I = 1, NXDIM
            IF (DGEO(I, J) .EQ. 0.D0) THEN
               DGEO(I, J) = 1.D0
            END IF
         END DO
      END DO

      IF (IBBL .NE. 0) THEN
         DO K = 1, NZDIM
            DO J = 1, NYDIM
               DO I = 1, NXDIM
                  BLD(I, J, K) = 0.D0
               END DO
            END DO
         END DO
         DO J = 1, NYDIM
            DO I = 1, NXDIM
               K = NBOT(I, J) + KSTR
               IF (      (IGEO(I, J) .GT. 0)
     &             .AND. (BGEO(I, J) .GT. 0.D0)) THEN
                  BLD(I, J, K) = 1.D0
               END IF
            END DO
         END DO
      END IF

      RETURN
      END
