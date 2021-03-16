      implicit none

      integer, parameter :: NNC = 1  !! N. of coorinates defined
      integer, parameter :: NSIGMX = 106  !! maximum number of layers

      character*16, allocatable :: cname(:)  !! Name of coordinate
      integer, allocatable :: nsig(:)        !! N. of layers of coordinate
      real*8, allocatable :: zref(:)         !! reference level [cm]

!     For the coordinate n, boundaries of lsig(k,n) are
!      lsigp(k,n) (lower) and lsigp(k+1,n) (upper).
!     dsig(k,n) = lsigp(k+1,n) - lsigp(k,n)
      real*8, allocatable :: lsig(:,:),      !! Sigma coordinate
     &                       lsigp(:,:),     !! Coordinate boundary
     &                       dsig(:,:)

      integer :: n,k

c     sigma2_l106 coordinate
      allocate(cname(NNC),zref(NNC),nsig(NNC))
      cname(1) = 'sigma2_l106     '
      nsig(1) = 106
      zref(1) = 2000.0D2

      allocate(lsig(NSIGMX,NNC),lsigp(1:NSIGMX+1,NNC),dsig(NSIGMX,NNC))
      lsigp(1,1) = -999.9D0
      lsig(1,1) = 27.875D0
      lsigp(2,1) = 28.0D0
      dsig(1,1) = lsigp(2,1)-lsigp(1,1)
      do n=2, 25  !! dsig=0.25 for 28.0-34.0
         dsig(n,1) = 0.25D0
         lsig(n,1) = 28.0D0 + (dble(n)-1.5D0) * 0.25D0
         lsigp(n+1,1) = 28.0D0 + (dble(n)-1.0D0) * 0.25D0
      end do
      do n=26, 50  !! dsig=0.1 for 34.0-36.5
         dsig(n,1) = 0.1D0
         lsig(n,1) = 34.0D0 + (dble(n)-25.5D0) * 0.1D0
         lsigp(n+1,1) = 34.0D0 + (dble(n)-25.0D0) * 0.1D0
      end do
      do n=51, 100  !! dsig=0.02 for 36.5-37.5
         dsig(n,1) = 0.02D0
         lsig(n,1) = 36.5D0 + (dble(n)-50.5D0) * 0.02D0
         lsigp(n+1,1) = 36.5D0 + (dble(n)-50.0D0) * 0.02D0
      end do
      do n=101, 105  !! dsig=0.1 for 37.5-38.0
         dsig(n,1) = 0.1D0
         lsig(n,1) = 37.5D0 + (dble(n)-100.5D0) * 0.1D0
         lsigp(n+1,1) = 37.5D0 + (dble(n)-100.0D0) * 0.1D0
      end do
      lsig(106,1) = 38.125D0
      lsigp(107,1) = 999.9D0
      dsig(106,1) = lsigp(107,1)-lsigp(106,1)

      open(21, file='SCOORD', form='unformatted')
      write(21) NNC
      write(21) cname
      write(21) zref
      write(21) nsig
      write(21) LSIG
      write(21) LSIGP
      write(21) DSIG
      close(21)

      stop
      end
