      implicit none

      integer NX, NY, NZ, KZ
      parameter(NX = 128, NY = 120, NZ = 40, KZ = 5)
      integer ISTR, JSTR, KSTR
      parameter(ISTR = 3, JSTR = 3, KSTR = 2)
      INTEGER   IEND, JEND, KEND
      PARAMETER(IEND = ISTR+NX-1, JEND = JSTR+NY-1, KEND = KSTR+NZ-1)
      integer NXDIM, NYDIM, NZDIM
      parameter(NXDIM = NX+4, NYDIM = NY+4, NZDIM = NZ+2)

      integer NXO, NYO
      parameter(NXO = 180, NYO = 90)
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

c      character CGEO(0:NX+1, 0:NY+1, 5)*1
c 1: global, 2: Atlantic, 3: Pacific, 4: Indian, 5: Indo-Pacific
      real*4 datin1(0:NX+1, 0:NY+1, NZ), datin2(0:NX+1, 0:NY+1, NZ)
      real*4 datou1(NXO, NYO, NZ), datou2(NXO, NYO, NZ)
      real*8 RLONT(0:NX+1), RLONV(0:NX)
      real*8 RLATT(0:NY+1), RLATV(0:NY)

      character chead(64)*16
      character GEOFIL*256
      integer i, j, k, l, irec, ifile

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

      RLONV(0) = ALONW
      RLONT(0) = ALONW - DX * 0.5D0
      do i = 1, NX
         RLONV(i) = RLONV(i-1) + DX
         RLONT(i) = RLONT(i-1) + DX
      end do
      RLONT(NX+1) = RLONT(NX) + DX

      RLATV(0) = ALATS
      RLATT(0) = ALATS - DY(1) * 0.5D0
      do j = 1, NY
         RLATV(j) = RLATV(j-1) + DY(j)
         RLATT(j) = (RLATV(j) + RLATV(j-1)) * 0.5D0
      end do
      RLATT(NY+1) = ALATN + DY(NY) * 0.5D0

      open(unit=51, file='ocn.grd',
     &     access='DIRECT', form='UNFORMATTED', recl=4*NXO*NYO)
      irec = 0

      open(unit=11, file='t',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=12, file='s',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=13, file='u',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=14, file='v',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=15, file='wzc',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=16, file='conv',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=21, file='sh',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=22, file='ai',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=23, file='hi',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=24, file='ui',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      open(unit=25, file='vi',
     &     access='SEQUENTIAL', form='UNFORMATTED')

      do ifile = 11, 12
         read(ifile) chead
         read(ifile) (((datin1(i, j, k), i = 1, NX),
     &                                   j = 1, NY),
     &                                   k = 1, NZ)
         do k = 1, NZ
            call intp(
     o        datou1(1, 1, k),
     i        datin1(0, 0, k),
     i        NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i        IBBL, k,
     c        ALONW, ALONE, ALATS, ALATN, UNDEF,
     d        NX, NY, NZ, NXO, NYO)
            irec = irec + 1
            write(51, rec=irec)
     &      ((datou1(i, j, k), i = 1, NXO), j = 1, NYO)
         end do
      end do

      read(13) chead
      read(13) (((datin1(i, j, k), i = 1, NX),
     &                             j = 1, NY),
     &                             k = 1, NZ)
      read(14) chead
      read(14) (((datin2(i, j, k), i = 1, NX),
     &                             j = 1, NY),
     &                             k = 1, NZ)
      do k = 1, NZ
         call intpv(
     o        datou1(1, 1, k), datou2(1, 1, k),
     i        datin1(0, 0, k), datin2(0, 0, k),
     i        NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i        IBBL, k,
     c        ALONW, ALONE, ALATS, ALATN, UNDEF,
     d        NX, NY, NZ, NXO, NYO)
      end do
      do k = 1, NZ
         irec = irec + 1
         write(51, rec=irec)
     &   ((datou1(i, j, k), i = 1, NXO), j = 1, NYO)
      end do
      do k = 1, NZ
         irec = irec + 1
         write(51, rec=irec)
     &   ((datou2(i, j, k), i = 1, NXO), j = 1, NYO)
      end do

      do ifile = 15, 16
         read(ifile) chead
         read(ifile) (((datin1(i, j, k), i = 1, NX),
     &                                   j = 1, NY),
     &                                   k = 1, NZ)
         do k = 1, NZ
            call intp(
     o        datou1(1, 1, k),
     i        datin1(0, 0, k),
     i        NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i        IBBL, k,
     c        ALONW, ALONE, ALATS, ALATN, UNDEF,
     d        NX, NY, NZ, NXO, NYO)
            irec = irec + 1
            write(51, rec=irec)
     &      ((datou1(i, j, k), i = 1, NXO), j = 1, NYO)
         end do
      end do

      do ifile = 21, 23
         read(ifile) chead
         read(ifile) ((datin1(i, j, 1), i = 1, NX), j = 1, NY)
         call intp(
     o        datou1(1, 1, 1),
     i        datin1(0, 0, 1),
     i        NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i        IBBL, 1,
     c        ALONW, ALONE, ALATS, ALATN, UNDEF,
     d        NX, NY, NZ, NXO, NYO)
         irec = irec + 1
         write(51, rec=irec)
     &   ((datou1(i, j, 1), i = 1, NXO), j = 1, NYO)
      end do

      read(24) chead
      read(24) ((datin1(i, j, 1), i = 1, NX), j = 1, NY)
      read(25) chead
      read(25) ((datin2(i, j, 1), i = 1, NX), j = 1, NY)
         call intpv(
     o        datou1(1, 1, 1), datou2(1, 1, 1),
     i        datin1(0, 0, 1), datin2(0, 0, 1),
     i        NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i        IBBL, 1,
     c        ALONW, ALONE, ALATS, ALATN, UNDEF,
     d        NX, NY, NZ, NXO, NYO)
         irec = irec + 1
         write(51, rec=irec)
     &   ((datou1(i, j, 1), i = 1, NXO), j = 1, NYO)
         irec = irec + 1
         write(51, rec=irec)
     &   ((datou2(i, j, 1), i = 1, NXO), j = 1, NYO)

      end

c **********************************************************************

      subroutine intpv(
     o                datou1, datou2,
     i                datin1, datin2,
     i                NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i                IBBL, k,
     c                ALONW, ALONE, ALATS, ALATN, UNDEF,
     d                NX, NY, NZ, NXO, NYO)
      implicit none

      integer NX, NY, NZ, NXO, NYO
      real*4  datin1(0:NX+1, 0:NY+1), datou1(NXO, NYO)
      real*4  datin2(0:NX+1, 0:NY+1), datou2(NXO, NYO)
      integer NBOT(-1:NX+2, -1:NY+2), IBBL, k
      real*8  BGEO(-1:NX+2, -1:NY+2)
      real*8  RLONT(0:NX+1), RLONV(0:NX), RLATT(0:NY+1), RLATV(0:NY)
      real*8  ALONW, ALONE, ALATS, ALATN
      real*4  UNDEF

      real*8 GLONO, GLATO, RLONO, RLATO, RANG
      integer i, j, io, jo, ilon, jlat
      real*4 d1, d2, d3, d4, d12, d34, ru, rv

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      do j = 1, NY
         datin1(0, j) = datin1(NX, j)
         datin2(0, j) = datin2(NX, j)
         datin1(NX+1, j) = datin1(1, j)
         datin2(NX+1, j) = datin2(1, j)
      end do
      do i = 0, NX+1
         datin1(i, 0) = datin1(i, 1)
         datin2(i, 0) = datin2(i, 1)
         datin1(i, NY+1) = datin1(i, NY)
         datin2(i, NY+1) = datin2(i, NY)
      end do

      do jo = 1, NYO
         GLATO = - 0.5D0 * PI
     &         + (dble(jo) - 0.5D0) * PI / dble(NYO)
         do io = 1, NXO
            GLONO = (dble(io) - 0.5D0) * 2.D0 * PI / dble(NXO)
            call GTOR(RLONO, RLATO, GLONO, GLATO)
            if (RLONO .lt. 0.D0) then
               RLONO = RLONO + 2.D0 * PI
            else if (RLONO .gt. 2.D0 * PI) then
               RLONO = RLONO - 2.D0 * PI
            end if
            call VROT(RANG, RLONO, RLATO)
            if ((RLONO .lt. ALONW) .or. (RLONO .gt. ALONE)) then
               ilon = 0
            else
               do i = 1, NX
                  if (RLONV(i) .ge. RLONO) then
                     ilon = i
                     go to 201
                  end if
               end do
 201           continue
            end if
            if ((RLATO .lt. ALATS) .or. (RLATO .gt. ALATN)) then
               jlat = 0
            else
               do j = 1, NY
                  if (RLATV(j) .ge. RLATO) then
                     jlat = j
                     go to 202
                  end if
               end do
 202           continue
            end if
            if ((ilon .eq. 0) .or. (jlat .eq. 0)) then
               datou1(io, jo) = UNDEF
               datou2(io, jo) = UNDEF
            else if ((IBBL .gt. 0) .and. (k .eq. NZ)) then
               if (BGEO(ilon, jlat) .eq. 0.D0) then
                  datou1(io, jo) = UNDEF
                  datou2(io, jo) = UNDEF
               else
                  d1 = datin1(ilon-1, jlat-1)
                  d2 = datin1(ilon, jlat-1)
                  d3 = datin1(ilon-1, jlat)
                  d4 = datin1(ilon, jlat)
                  d12 = (  d1 * (RLONV(ilon) - RLONO)
     &                   + d2 * (RLONO - RLONV(ilon-1)))
     &                  / (RLONV(ilon) - RLONV(ilon-1))
                  d34 = (  d3 * (RLONV(ilon) - RLONO)
     &                   + d4 * (RLONO - RLONV(ilon-1)))
     &                  / (RLONV(ilon) - RLONV(ilon-1))
                  ru = (  d12 * (RLATV(jlat) - RLATO)
     &                  + d34 * (RLATO - RLATV(jlat-1)))
     &                 / (RLATV(jlat) - RLATV(jlat-1))
                  d1 = datin2(ilon-1, jlat-1)
                  d2 = datin2(ilon, jlat-1)
                  d3 = datin2(ilon-1, jlat)
                  d4 = datin2(ilon, jlat)
                  d12 = (  d1 * (RLONV(ilon) - RLONO)
     &                   + d2 * (RLONO - RLONV(ilon-1)))
     &                  / (RLONV(ilon) - RLONV(ilon-1))
                  d34 = (  d3 * (RLONV(ilon) - RLONO)
     &                   + d4 * (RLONO - RLONV(ilon-1)))
     &                  / (RLONV(ilon) - RLONV(ilon-1))
                  rv = (  d12 * (RLATV(jlat) - RLATO)
     &                  + d34 * (RLATO - RLATV(jlat-1)))
     &                 / (RLATV(jlat) - RLATV(jlat-1))
                  datou1(io, jo) = ru * COS(RANG) - rv * SIN(RANG)
                  datou2(io, jo) = ru * SIN(RANG) + rv * COS(RANG)
               end if
            else if (NBOT(ilon, jlat) .lt. k) then
               datou1(io, jo) = UNDEF
               datou2(io, jo) = UNDEF
            else
               d1 = datin1(ilon-1, jlat-1)
               d2 = datin1(ilon, jlat-1)
               d3 = datin1(ilon-1, jlat)
               d4 = datin1(ilon, jlat)
               d12 = (  d1 * (RLONV(ilon) - RLONO)
     &                + d2 * (RLONO - RLONV(ilon-1)))
     &               / (RLONV(ilon) - RLONV(ilon-1))
               d34 = (  d3 * (RLONV(ilon) - RLONO)
     &                + d4 * (RLONO - RLONV(ilon-1)))
     &               / (RLONV(ilon) - RLONV(ilon-1))
               ru = (  d12 * (RLATV(jlat) - RLATO)
     &               + d34 * (RLATO - RLATV(jlat-1)))
     &              / (RLATV(jlat) - RLATV(jlat-1))
               d1 = datin2(ilon-1, jlat-1)
               d2 = datin2(ilon, jlat-1)
               d3 = datin2(ilon-1, jlat)
               d4 = datin2(ilon, jlat)
               d12 = (  d1 * (RLONV(ilon) - RLONO)
     &                + d2 * (RLONO - RLONV(ilon-1)))
     &               / (RLONV(ilon) - RLONV(ilon-1))
               d34 = (  d3 * (RLONV(ilon) - RLONO)
     &                + d4 * (RLONO - RLONV(ilon-1)))
     &               / (RLONV(ilon) - RLONV(ilon-1))
               rv = (  d12 * (RLATV(jlat) - RLATO)
     &               + d34 * (RLATO - RLATV(jlat-1)))
     &              / (RLATV(jlat) - RLATV(jlat-1))
               datou1(io, jo) = ru * COS(RANG) - rv * SIN(RANG)
               datou2(io, jo) = ru * SIN(RANG) + rv * COS(RANG)
            end if
         end do
      end do

      return
      end

c **********************************************************************

      subroutine intp(
     o                datou,
     i                datin,
     i                NBOT, BGEO, RLONT, RLONV, RLATT, RLATV,
     i                IBBL, k,
     c                ALONW, ALONE, ALATS, ALATN, UNDEF,
     d                NX, NY, NZ, NXO, NYO)
      implicit none

      integer NX, NY, NZ, NXO, NYO
      real*4  datin(0:NX+1, 0:NY+1), datou(NXO, NYO)
      integer NBOT(-1:NX+2, -1:NY+2), IBBL, k
      real*8  BGEO(-1:NX+2, -1:NY+2)
      real*8  RLONT(0:NX+1), RLONV(0:NX), RLATT(0:NY+1), RLATV(0:NY)
      real*8  ALONW, ALONE, ALATS, ALATN
      real*4  UNDEF

      real*8 GLONO, GLATO, RLONO, RLATO
      integer i, j, io, jo, ilon, jlat
      real*4 d1, d2, d3, d4, d12, d34

      REAL*8      PI,    RAD, RADIUS,  OMEGA
      COMPLEX*16 CUNITI
      COMMON /CONST/ PI, RAD, RADIUS, OMEGA, CUNITI

      do j = 1, NY
         datin(0, j) = datin(NX, j)
         datin(NX+1, j) = datin(1, j)
      end do
      do i = 0, NX+1
         datin(i, 0) = datin(i, 1)
         datin(i, NY+1) = datin(i, NY)
      end do

      do jo = 1, NYO
         GLATO = - 0.5D0 * PI
     &         + (dble(jo) - 0.5D0) * PI / dble(NYO)
         do io = 1, NXO
            GLONO = (dble(io) - 0.5D0) * 2.D0 * PI / dble(NXO)
            call GTOR(RLONO, RLATO, GLONO, GLATO)
            if (RLONO .lt. 0.D0) then
               RLONO = RLONO + 2.D0 * PI
            else if (RLONO .gt. 2.D0 * PI) then
               RLONO = RLONO - 2.D0 * PI
            end if
            if ((RLONO .lt. ALONW) .or. (RLONO .gt. ALONE)) then
               ilon = 0
            else
               do i = 1, NX+1
                  if (RLONT(i) .gt. RLONO) then
                     ilon = i
                     go to 101
                  end if
               end do
 101           continue
            end if
            if ((RLATO .lt. ALATS) .or. (RLATO .gt. ALATN)) then
               jlat = 0
            else
               do j = 1, NY+1
                  if (RLATT(j) .gt. RLATO) then
                     jlat = j
                     go to 102
                  end if
               end do
 102           continue
            end if
            if ((ilon .eq. 0) .or. (jlat .eq. 0)) then
               datou(io, jo) = UNDEF
            else if ((IBBL .gt. 0) .and. (k .eq. NZ)) then
               if (     (BGEO(ilon-1, jlat-1) .eq. 0.D0)
     &             .or. (BGEO(ilon, jlat-1) .eq. 0.D0)
     &             .or. (BGEO(ilon-1, jlat) .eq. 0.D0)
     &             .or. (BGEO(ilon, jlat) .eq. 0.D0)) then
                  datou(io, jo) = UNDEF
               else
                  d1 = datin(ilon-1, jlat-1)
                  d2 = datin(ilon, jlat-1)
                  d3 = datin(ilon-1, jlat)
                  d4 = datin(ilon, jlat)
                  d12 = (  d1 * (RLONT(ilon) - RLONO)
     &                + d2 * (RLONO - RLONT(ilon-1)))
     &               / (RLONT(ilon) - RLONT(ilon-1))
                  d34 = (  d3 * (RLONT(ilon) - RLONO)
     &                + d4 * (RLONO - RLONT(ilon-1)))
     &               / (RLONT(ilon) - RLONT(ilon-1))
                  datou(io, jo) = (  d12 * (RLATT(jlat) - RLATO)
     &                          + d34 * (RLATO - RLATT(jlat-1)))
     &                         / (RLATT(jlat) - RLATT(jlat-1))
               end if
            else if (     (NBOT(ilon-1, jlat-1) .lt. k)
     &               .or. (NBOT(ilon, jlat-1) .lt. k)
     &               .or. (NBOT(ilon-1, jlat) .lt. k)
     &               .or. (NBOT(ilon, jlat) .lt. k)) then
               datou(io, jo) = UNDEF
            else
               d1 = datin(ilon-1, jlat-1)
               d2 = datin(ilon, jlat-1)
               d3 = datin(ilon-1, jlat)
               d4 = datin(ilon, jlat)
               d12 = (  d1 * (RLONT(ilon) - RLONO)
     &                + d2 * (RLONO - RLONT(ilon-1)))
     &               / (RLONT(ilon) - RLONT(ilon-1))
               d34 = (  d3 * (RLONT(ilon) - RLONO)
     &                + d4 * (RLONO - RLONT(ilon-1)))
     &               / (RLONT(ilon) - RLONT(ilon-1))
               datou(io, jo) = (  d12 * (RLATT(jlat) - RLATO)
     &                          + d34 * (RLATO - RLATT(jlat-1)))
     &                         / (RLATT(jlat) - RLATT(jlat-1))
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
