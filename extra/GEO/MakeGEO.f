      integer NX, NY, NZ, KZ
      parameter(NX = 128, NY = 120, NZ = 40, KZ = 5)
      integer NX5, NY5
      parameter(NX5 = 4320, NY5 = 2160)

      integer ITOPO5(NX5, NY5)
      integer ISTR(NX), IEND(NX), JSTR(NY), JEND(NY)
      real*4  DDEPTH(0:NZ), DEPTH(NZ)
      real*4  AGEO(NX, NY), BGEO(NX, NY)
      integer IGEO(NX, NY)

      real*8 dz(NZ)
      real*8 alonw, alone, dx
      real*8 alats, alatn, ys, yn, dymerc, y1(NY), y2(NY)
      real*8 pi, rad
      real*8 plon, plat, qlon, qlat

      data plon, plat / -40.D0, 75.D0 /
      data qlon, qlat / -40.D0, 90.D0 /
      data alonw, alone / 0.D0, 360.D0 /
      data alats, alatn / -78.D0, 85.D0 /

      data DDEPTH /    0.0,
     &                 5.0,   10.0,   18.0,   30.0,   50.0,
     &                80.0,  120.0,  170.0,  230.0,  300.0,
     &               380.0,  470.0,  570.0,  680.0,  800.0,
     &               930.0, 1070.0, 1220.0, 1380.0, 1550.0,
     &              1730.0, 1915.0, 2105.0, 2300.0, 2500.0,
     &              2700.0, 2900.0, 3100.0, 3300.0, 3500.0,
     &              3700.0, 3900.0, 4100.0, 4300.0, 4500.0,
     &              4700.0, 4900.0, 5100.0, 5300.0, 5400.0 /

      pi = 4.D0 * atan(1.D0)
      rad = pi / 180.D0
      alats = alats * rad
      alatn = alatn * rad

      dx = (alone - alonw) / dble(NX)
      ys = log(tan(alats * 0.5D0 + pi * 0.25D0))
      yn = log(tan(alatn * 0.5D0 + pi * 0.25D0))
      dymerc = (yn - ys) / dble(NY)
      do j = 1, NY
         y1(j) = ys + dymerc * dble(j - 1)
         y2(j) = ys + dymerc * dble(j)
c         dy(j) = 2.D0 * (atan(exp(y2(j))) - atan(exp(y1(j))))
         write(97, *) j, 2.D0*(atan(exp(y2(j)))-atan(exp(y1(j))))
      end do
      do j = 1, NY
         y1(j) = (2.D0 * atan(exp(y1(j))) - pi * 0.5D0) / rad
         y2(j) = (2.D0 * atan(exp(y2(j))) - pi * 0.5D0) / rad
c         dy(j) = dy(j) / rad
      end do
      do k = 1, NZ
         dz(k) = DDEPTH(k) - DDEPTH(k-1)
      end do
      do k = 1, NZ
         DEPTH(k) = (DDEPTH(k-1) + DDEPTH(k)) * 0.5
      end do

      do i = 1, NX
         ISTR(i) = int((alonw + dx * dble(i - 1)) * 12.D0)
         IEND(i) = int((alonw + dx * dble(i)) * 12.D0)
      end do
      do j = 1, NY
c         JSTR(j) = int(y1(j) + 90.D0) * 12
c         JEND(j) = int(y2(j) + 90.D0) * 12
         JSTR(j) = int((y1(j) + 90.D0) * 12.D0)
         JEND(j) = int((y2(j) + 90.D0) * 12.D0)
      end do

      do j = 1, NY
         do i = 1, NX
            IGEO(i, j) = 0
         end do
      end do

      open(unit=11, file='./ETOPO5-c', status='OLD',
     &     access='SEQUENTIAL', form='UNFORMATTED')
      do j = NY5, 1, -1
         read(11) (ITOPO5(i, j), i = 1, NX5)
      end do
      close(unit=11)

      do j = 1, NY
         do i = 1, NX
            dtopo = 0.0
            icount = 0
            do jdiv = JSTR(j), JEND(j)
               do idiv = ISTR(i), IEND(i)
                  icount = icount + 1
                  dtopo = dtopo + real(ITOPO5(idiv, jdiv))
               end do
            end do
            dtopo = - dtopo / real(icount)
            if (dtopo .ge. DEPTH(NZ)) then
               IGEO(i, j) = NZ
               AGEO(i, j) = real(NZ)
            else
               do k = 1, NZ-1
                  if (      (dtopo .ge. DEPTH(k))
     &                .and. (dtopo .lt. DEPTH(k+1))) then
                     IGEO(i, j) = k
                     AGEO(i, j) = real(k)
     &                          + (dtopo - DEPTH(k))
     &                             / (DEPTH(k+1) - DEPTH(k))
                  end if
               end do
            end if
         end do
      end do

      do j = 1, NY
         do i = 1, NX
            if (AGEO(i, j) .gt. 0.0) then
               AGEO(i, j) = max(AGEO(i, j), real(KZ))
            end if
         end do
      end do

      do j = 1, NY
         do i = 1, NX
            BGEO(i, j) = AGEO(i, j)
         end do
      end do
      do j = 1, NY
         BGEO(0, j) = BGEO(NX, j)
         BGEO(NX+1, j) = BGEO(1, j)
      end do

      do j = 2, NY-1
         do i = 1, NX
            if (BGEO(i, j) .gt. 0.) then
               AGEO(i, j) = BGEO(i, j) * 0.5
               if (BGEO(i+1, j) .gt. 0) then
                  AGEO(i, j) = AGEO(i, j) + BGEO(i+1, j) * 0.125
               else
                  AGEO(i, j) = AGEO(i, j) + BGEO(i, j) * 0.125
               end if
               if (BGEO(i-1, j) .gt. 0) then
                  AGEO(i, j) = AGEO(i, j) + BGEO(i-1, j) * 0.125
               else
                  AGEO(i, j) = AGEO(i, j) + BGEO(i, j) * 0.125
               end if
               if (BGEO(i, j+1) .gt. 0) then
                  AGEO(i, j) = AGEO(i, j) + BGEO(i, j+1) * 0.125
               else
                  AGEO(i, j) = AGEO(i, j) + BGEO(i, j) * 0.125
               end if
               if (BGEO(i, j-1) .gt. 0) then
                  AGEO(i, j) = AGEO(i, j) + BGEO(i, j-1) * 0.125
               else
                  AGEO(i, j) = AGEO(i, j) + BGEO(i, j) * 0.125
               end if
            end if
         end do
      end do

      open(unit=12, file='GEO',
     &     access='SEQUENTIAL', form='FORMATTED')
      write(12, *) NX
      write(12, *) NY
      write(12, *) NZ
      write(12, *) KZ
      write(12, *) 1, 1
      write(12, *) alonw, alone
      write(12, *) alats/rad, alatn/rad
      write(12, *) 1
      write(12, *) dz
      write(12, *) plon, plat
      write(12, *) qlon, qlat
      write(12, *) 0
      do j = NY, 1, -1
         write(12, '(512F5.1)') (AGEO(i, j), i = 1, NX)
      end do
      close(unit=12)

      end
