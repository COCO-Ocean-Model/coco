      integer NX, NY, NZ, KZ
      parameter(NX = 128, NY = 120, NZ = 40, KZ = 5)
      integer ISTR, JSTR, KSTR
      parameter(ISTR = 3, JSTR = 3, KSTR = 2)
      INTEGER   IEND, JEND, KEND
      PARAMETER(IEND = ISTR+NX-1, JEND = JSTR+NY-1, KEND = KSTR+NZ-1)
      integer NXDIM, NYDIM, NZDIM
      parameter(NXDIM = NX+4, NYDIM = NY+4, NZDIM = NZ+2)

      integer NYP
      parameter(NYP = 90)
      integer NXDIV, NYDIV
      parameter(NXDIV = 3, NYDIV = 3)
      real*4 UNDEF
      parameter(UNDEF = -999.)

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

      character CGEO(0:NX+1, 0:NY+1, 5)*1
c 1: global, 2: Atlantic, 3: Pacific, 4: Indian, 5: Indo-Pacific
      real*4 v(0:NX, 0:NY, NZ), w(NX, NY, NZ), ww(NX, NY, NZ)
      real*4 ps(0:NZ), p(0:NYP, 0:NZ)
      real*4 wint(NYP, 0:NZ)
      real*8 GLATP(0:NYP)
      real*4 platp(0:NYP), ddep(0:NZ)
      integer NBOTP(NYP)

      real*8 RLON, RLAT, GLON, GLAT
      real*8 RLON0, RLAT0, RLON1, RLAT1
      real*8 SLATS, SLATN

      character chead(64)*16
      character GEOFIL*256

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      data chead / 64*'                ' /
      data GEOFIL / 'GEO.std2_bbl' /

      PI = ATAN(1.D0) * 4.D0
      RAD = PI / 1.8D2
      RADIUS = 6370.D5
      OMEGA = 2.D0 * PI / 8.64D4
      CUNITI = (0.D0, 1.D0)

      do j = 0, NYP
         GLATP(j) = -90.D0 + 180.D0 * dble(j) / dble(NYP)
         platp(j) = sngl(GLATP(j))
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
     I            GEOFIL)

      ddep(0) = 0.
      do k = 1, NZ
         ddep(k) = ddep(k-1) + sngl(DZ(k))
      end do

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

      do ib = 1, 5
         if (ib .eq. 1) then
            open(unit=11, file='Ocean.std2_bbl',
     &           access='SEQUENTIAL', form='FORMATTED')
         else if (ib .eq. 2) then
            open(unit=11, file='Atlantic.std2_bbl',
     &           access='SEQUENTIAL', form='FORMATTED')
         else if (ib .eq. 3) then
            open(unit=11, file='Pacific.std2_bbl',
     &           access='SEQUENTIAL', form='FORMATTED')
         else if (ib .eq. 4) then
            open(unit=11, file='Indian.std2_bbl',
     &           access='SEQUENTIAL', form='FORMATTED')
         else
            open(unit=11, file='IndoPac.std2_bbl',
     &           access='SEQUENTIAL', form='FORMATTED')
         end if
         do j = NY, 1, -1
            read(11, '(128A1)') (CGEO(i, j, ib), i = 1, NX)
         end do
         close(unit=11)
         do i = 0, NX+1
            CGEO(i, 0, ib) = '0'
            CGEO(i, NY+1, ib) = '0'
         end do
         do j = 1, NY
            CGEO(0, j, ib) = CGEO(NX, j, ib)
            CGEO(NX+1, j, ib) = CGEO(1, j, ib)
         end do
      end do

      open(unit=21, file='moc.grd',
     &     access='DIRECT', form='UNFORMATTED', recl=4*(NYP+1))
      irec = 0

      open(unit=11, file='v',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=12, file='wzc',
     &     access='SEQUENTIAL', form='UNFORMATTED')

      read(11) chead
      read(11) (((v(i, j, k), i = 1, NX), j = 1, NY), k = 1, NZ)
      do k = 1, NZ
         do j = 1, NY
            v(0, j, k) = v(NX, j, k)
         end do
         do i = 0, NX
            v(i, 0, k) = 0.
         end do
      end do
      read(12) chead
      read(12) ww

      do ib = 1, 5

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               w(i, j, k) = ww(i, j, k)
            end do
         end do
      end do

      do j = 1, NY
         do i = 1, NX
            if (CGEO(i, j, ib) .ne. '.') then
               do k = 1, NZ
                  w(i, j, k) = 0.
               end do
            end if
         end do
      end do

      KYS = NY
      do j = NY, 1, -1
         do i = 1, NX
            if (CGEO(i, j, ib) .eq. '*') go to 99
         end do
         KYS = j
      end do
 99   continue

      SLATS = PI / 2.
      SLATN = - PI / 2.
      if (KYS .eq. 1) then
         do i = 1, NX
            CALL IJVTOR(RLON, RLAT, i-1, KYS-1,
     &                  ALONW, ALATS, DX, DY, NY, NYDIM, JSTR)
            CALL RTOG(GLON, GLAT, RLON, RLAT)
            SLATS = MIN(SLATS, GLAT)
            SLATN = MAX(SLATN, GLAT)
            CALL IJVTOR(RLON, RLAT, i, KYS-1,
     &                  ALONW, ALATS, DX, DY, NY, NYDIM, JSTR)
            CALL RTOG(GLON, GLAT, RLON, RLAT)
            SLATS = MIN(SLATS, GLAT)
            SLATN = MAX(SLATN, GLAT)
         end do
      else
         do i = 1, NX
            if (CGEO(i, KYS, ib) .eq. '.') then
               CALL IJVTOR(RLON, RLAT, i-1, KYS-1,
     &                     ALONW, ALATS, DX, DY, NY, NYDIM, JSTR)
               CALL RTOG(GLON, GLAT, RLON, RLAT)
               SLATS = MIN(SLATS, GLAT)
               SLATN = MAX(SLATN, GLAT)
               CALL IJVTOR(RLON, RLAT, i, KYS-1,
     &                     ALONW, ALATS, DX, DY, NY, NYDIM, JSTR)
               CALL RTOG(GLON, GLAT, RLON, RLAT)
               SLATS = MIN(SLATS, GLAT)
               SLATN = MAX(SLATN, GLAT)
            end if
         end do
      end if
      SLATS = SLATS / RAD
      SLATN = SLATN / RAD

      do k = 0, NZ
         ps(k) = 0.
      end do
      do k = NZ-1, 0, -1
         do i = 1, NX
            if (CGEO(i, KYS, ib) .eq. '.') then
               ps(k) = ps(k)
     &               + DZV(i-1, KYS-1, k+1) * v(i-1, KYS-1, k+1) *
     &                 DX * RADIUS * HXU(i-1, KYS-1) * 1.e-12 * 0.5
     &               + DZV(i, KYS-1, k+1) * v(i, KYS-1, k+1) *
     &                 DX * RADIUS * HXU(i, KYS-1) * 1.e-12 * 0.5
            end if
         end do
         ps(k) = ps(k+1) - ps(k)
      end do

      do k = 0, NZ
         do j = 0, NYP
            p(j, k) = 0.
         end do
      end do

      do j = 0, NYP
         if (GLATP(j) .le. SLATS) then
            do k = 0, NZ
               p(j, k) = ps(k)
            end do
         end if
      end do

      do k = 0, NZ
         do j = 1, NYP
            wint(j, k) = 0.
         end do
      end do

      do j = 1, NYP
         NBOTP(j) = 0
      end do

      do j = 1, NY
c         write(*, *) j
         do i = 1, NX
            if (CGEO(i, j, ib) .eq. '.') then
            CALL IJVTOR(RLON0, RLAT0, i-1, j-1,
     &                  ALONW, ALATS, DX, DY, NY, NYDIM, JSTR)
            CALL IJVTOR(RLON1, RLAT1, i, j,
     &                  ALONW, ALATS, DX, DY, NY, NYDIM, JSTR)
            do jdiv = 1, NYDIV
               RLAT = RLAT0
     &              + (RLAT1 - RLAT0) / dble(NYDIV) *
     &                (dble(jdiv - 1) + 0.5D0)
               do idiv = 1, NXDIV
                  RLON = RLON0
     &                 + (RLON1 - RLON0) / dble(NXDIV) *
     &                   (dble(idiv - 1) + 0.5D0)
                  CALL RTOG(GLON, GLAT, RLON, RLAT)
                  GLON = GLON / RAD
                  GLAT = GLAT / RAD
                  do jp = 1, NYP
                     if (      (GLAT .gt. GLATP(jp-1))
     &                   .and. (GLAT .le. GLATP(jp))) then
                        NBOTP(jp) = max(NBOTP(jp), NBOT(i, j))
                        do k = 0, NZ-1
                           wint(jp, k) = wint(jp, k)
     &                                 + w(i, j, k+1) *
     &                                   HXT(i, j) * DX * RADIUS *
     &                                   HYT(i, j) * DY(j) * RADIUS
     &                                   / dble(NXDIV * NYDIV) *
     &                                   1.e-12
                        end do
                     end if
                  end do
               end do
            end do
         end if
         end do
      end do

      do k = 0, NZ
         do j = 1, NYP
            p(j, k) = p(j-1, k) + wint(j, k)
         end do
      end do

      do j = 1, NYP-1
         nbotpv = min(NBOTP(j), NBOTP(j+1))
         if (IBBL .gt. 0) then
            p(j, nbotpv) = p(j, NZ-1)
         end if
         do k = nbotpv+1, NZ
            p(j, k) = UNDEF
c            p(j, k) = 0.
         end do
      end do

      if (KYS .gt. 1) then
         do j = 0, NYP
            if (GLATP(j) .lt. SLATN) then
               do k = 0, NZ
                  p(j, k) = UNDEF
               end do
            end if
         end do
      end if

      do k = 0, NZ
         irec = irec + 1
         write(21, rec=irec) (p(j, k), j = 0, NYP)
      end do

      end do ! ib loop

      END

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
      COMPLEX*16    ZPT,    ZPV,    ZQT,    ZQV
      COMPLEX*16     Z1,     Z2,     Z3,     W3
      REAL*8   RLONT,  RLATT,  RLONV,  RLATV
      REAL*8      YS,     YN, DYMERC,     Y1,     Y2

      REAL*8  FLGBBL
      INTEGER      I,      J,      K
      REAL*8    RLON,   RLAT,   GLON,   GLAT

      SAVE

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      REAL*8  ZTOLON, ZTOLAT
      REAL*8   FMTRC
      COMPLEX*16   CONF,  DCONF,   STOZ, CONFIN
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
      CONFIN(Z, ZA, ZB, ZC, WC) =
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
      ZQT = CONFIN(ZPT, Z1, Z2, Z3, W3)
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
