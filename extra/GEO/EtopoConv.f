      implicit none

      integer NX, NY
      parameter(NX = 4320, NY = 2160)

      integer ITOPO(0:NX+1, 0:NY+1), ITOPOR(NX, NY)
      real*8 alon(0:NX+1), alat(0:NY+1)

      real*8 glon, glat
      complex*16 zp, zq

      real*8 pi, rad
      complex*16 cuniti

      real*8 plon, plat
      real*8 qlon, qlat
      real*8 elon, elat
      real*8 cx, cy, cz
      complex*16 z1, z2, z3, w3

      integer i, j

      data plon, plat / -40.D0, 75.D0 /
      data qlon, qlat / -40.D0, -90.D0 /

      pi = 4.D0 * atan(1.D0)
      rad = pi / 180.D0
      cuniti = (0.D0, 1.D0)

      plon = plon * rad
      plat = plat * rad
      qlon = qlon * rad
      qlat = qlat * rad
c      w3 = (1.D0, 0.D0)
      w3 = (-1.D0, 0.D0)
      cx = cos(plon) * cos(plat) + cos(qlon) * cos(qlat)
      cy = sin(plon) * cos(plat) + sin(qlon) * cos(qlat)
      cz = sin(plat) + sin(qlat)
      elon = atan2(cy, cx)
      cx = sqrt(cx*cx + cy*cy)
      elat = atan2(cz, cx)

      z1 = tan(pi * 0.25D0 - plat * 0.5D0) * exp(plon * cuniti)
      z2 = tan(pi * 0.25D0 - qlat * 0.5D0) * exp(qlon * cuniti)
      z3 = tan(pi * 0.25D0 - elat * 0.5D0) * exp(elon * cuniti)

      open(unit=11, file='ETOPO5.reformat.data',
     &     form='UNFORMATTED', access='SEQUENTIAL')
      do j = NY, 1, -1
         read(11) (ITOPO(i, j), i = 1, NX)
      end do
      do j = 1, NY
         ITOPO(0, j) = ITOPO(NX, j)
         ITOPO(NX+1, j) = ITOPO(1, j)
      end do
      do i = 0, NX+1
         ITOPO(i, 0) = ITOPO(i, 1)
         ITOPO(i, NY+1) = ITOPO(i, NY)
      end do

      do j = 0, NY+1
         alat(j) = (-90.D0 + (dble(j) - 0.5D0) / 12.D0) * rad
      end do
      do i = 0, NX+1
         alon(i) = (dble(i) - 0.5D0) / 12.D0 * rad
      end do

      do j = 1, NY
         write(*, *) j
         do i = 1, NX
            zq = tan(pi * 0.25D0 - alat(j) * 0.5D0) *
     &           exp(alon(i) * cuniti)
            zp = (-z2 * (z3 - z1) * zq / w3 + z1 * (z3 - z2))
     &           / (-(z3 - z1) * zq / w3 + z3 - z2)
            glon = atan2(dimag(zp), dreal(zp))
            if (glon .lt. 0.D0) then
               glon = glon + 2.D0 * pi
            end if
            glat = pi * 0.5D0 - 2.D0 * atan(abs(zp))
            call intp(
     o             ITOPOR(i, j),
     i             alon, alat, ITOPO, glon, glat,
     d             NX, NY)
         end do
      end do

      open(unit=12, file='ETOPO5-c',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      do j = NY, 1, -1
         write(12) (ITOPOR(i, j), i = 1, NX)
      end do

      end

c *********************************************************************

      subroutine intp(
     o                idept,
     i                alon, alat, idep, glon, glat,
     d                NX, NY, NZ)

      integer NX, NY, NZ
      real*8 alon(0:NX+1), alat(0:NY+1)
      real*8 glon, glat
      integer idept, idep(0:NX+1, 0:NY+1)

      real*4 d1, d2, d3, d4, d12, d34

      do i = 1, NX+1
         if (alon(i) .gt. glon) then
            ilon = i
            go to 101
         end if
      end do
 101  continue
      do j = 1, NY+1
         if (alat(j) .gt. glat) then
            jlat = j
            go to 102
         end if
      end do
 102  continue

      d1 = real(idep(ilon-1, jlat-1))
      d2 = real(idep(ilon, jlat-1))
      d3 = real(idep(ilon-1, jlat))
      d4 = real(idep(ilon, jlat))

      d12 = (  d1 * (alon(ilon) - glon)
     &       + d2 * (glon - alon(ilon-1)))
     &      / (alon(ilon) - alon(ilon-1))
      d34 = (  d3 * (alon(ilon) - glon)
     &       + d4 * (glon - alon(ilon-1)))
     &      / (alon(ilon) - alon(ilon-1))
      idept = (  d12 * (alat(jlat) - glat)
     &         + d34 * (glat - alat(jlat-1)))
     &        / (alat(jlat) - alat(jlat-1))

      return
      end
