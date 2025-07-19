module tslvt

! --- information -----------------------------------------------------
!
!  Time integration of the tracer equations.
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1-OMIP
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
!     '12.08.15  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim, only :                                   &
       &   nxyzdm,  nxydim,   nzdim,  ntdim,           &
       &     nxyg,      nz, mpi_comm_ogcm
  use zocfil, only :                                   &
       &      ncf

  implicit none
  private

  real(8),        save  ::  rgamma(nzdim)
                  
  real(8),        save  ::      aa(nxydim, nzdim)
  real(8),        save  ::      ab(nxydim, nzdim)
  real(8),        save  ::      ac(nxydim, nzdim)
                  
  real(8),        save  ::      dh(nxydim)  
  real(8),        save  ::   hzbot(nxydim),   hxbot(nxydim)
  real(8),        save  ::       gamma(nz)
                  
  integer(4)            ::   ifpar,  jfpar,  istat

!---- for tunnel diffusion
  integer(4), parameter ::  ntdim1 = ntdim + 1
  integer(4), parameter ::  ntunnel = 20 
  integer(4), parameter ::  ntngr = ntunnel
  integer(4), parameter ::  ntigr = ntunnel
  real(8),                 save  ::   crst(ntngr),  crsti(ntigr)
  integer(4),              save  ::   ntni
  integer(4),              save  ::   ntns,         ntnis
  integer(4),              save  ::   ltns(ntunnel*nzdim*2+1)
  integer(4),              save  ::  ltnsd(ntunnel*nzdim*2+1)
  integer(4),              save  ::  ltnsg(ntunnel*nzdim*2+1)
  integer(4),              save  ::  ltnis(ntunnel*2+1)
  integer(4),              save  :: ltnisd(ntunnel*2+1)
  integer(4),              save  :: ltnisg(ntunnel*2+1)
  integer(4),              save  :: ntnist(0:nxyg), ntnsst(0:nxyg)
  real(8),                 save  ::    vlm(nxyzdm),  area(nxydim)
  integer(4), allocatable, save  ::  isreq(:),      irreq(:)
  data crst  / ntngr*0.d0 /
  data crsti / ntigr*0.d0 /

  character(16),  save  :: ctxnam(ntdim)
  character(16),  save  :: ctynam(ntdim)
  character(16),  save  :: ctznam(ntdim)

  public  ::  slvtrc     !   aprdc.F
  public  ::  svtset     !   aocea.F
  public  ::  tundif     !   aprdc.F

  namelist /nmacct/ gamma

contains

! **********************************************************************

  subroutine svtset

    use zocdim,  only :                                &
         &  nxydim,  nxdim,  nydim,  nzdim,  ntdim,    &
         &    istr,   iend,   jstr,   jend,            &
         &    kstr,   kend,     kz,                    &
         &   ijstr,  ijend, ijtstr,  ijtend,           &
#ifndef OPT_IO_COCOMPI
         &     nxg,    nyg, nxgdim,  nygdim,           &
         &   igstr,  jgstr,                            &
#endif
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dx,     dy,    dz0,     hxt,    hyt
    use zocmsk,  only :                                &
         &   amskt,   nbot
    use bshft
    use qckag
    use ufile
    use zocfil,  only :                                &
         &     ncf

    implicit none
#include "mpif.h"

    integer(4)        ::      ij,      k,      n,     l
    integer(4)        ::       i,      j
    logical,   save   ::  ofirst = .true.

    if ( ofirst ) then

       ofirst = .false.
       gamma(1:nz) = 1.d0
       call rewnml( ifpar, jfpar )
       write(jfpar, *) '*** svtset ***'
       read(ifpar, nmacct, iostat = istat )
       call cstnml( jfpar, 'svtset', 'nmacct', istat )
       write(jfpar, nmacct)
       write(jfpar, *) ' level     gamma'
       do k = 1, nz
          write(jfpar, '(i6,f10.4)') k, gamma(k)
       end do
       do k = kstr, kend
          rgamma(k) = 1.d0 / gamma(k-kstr+1)
       end do

!----- for tunnel diffusion
       call ttsset
       do l = 3, ntdim
          write(ctxnam(l), '(a5,i2.2)') 'ftrcx', l
          write(ctynam(l), '(a5,i2.2)') 'ftrcy', l
          write(ctznam(l), '(a5,i2.2)') 'ftrcz', l
       end do

    end if

  end subroutine svtset

! **********************************************************************

  subroutine slvtrc(    tx,   adt, diffz,    ha,    hb,    hc,    ft )

    use zocdim,  only :                                &
         &      nx,     ny, nxydim,  nzdim,  ntdim,    &
         &    kstr,   kend,     kz,                    &
         &   ijstr,  ijend, ijtstr,  ijtend,           &
         &     nic,                                    &
         &  inodes, jnodes,                            &
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dz,    dz0,     ds,   zbot,            &
         &      ts
    use zocmsk,  only :                                &
#ifdef OPT_BBL
         &   amsktb,                                   &
#endif
         &   amskt,  nbot
    use zocnod,  only :                                &
         &   iroot,  ierr,  myrank
    use utrdg
    use qckag
    use qckot

    implicit none

#include "mpif.h"

    real(8),   intent(inout)  ::     tx(nxydim, nzdim, ntdim)
    real(8),   intent(inout)  ::    adt(nxydim, nzdim, ntdim)
    real(8),   intent(in)     ::  diffz(nxydim, nzdim)
    real(8),   intent(in)     ::     ha(nxydim)
    real(8),   intent(in)     ::     hb(nxydim)
    real(8),   intent(in)     ::     hc(nxydim)
    real(8),   intent(in)     ::     ft(nxydim, ntdim)

    real(8)                   ::    dzb,  dzt
    integer(4)                ::     ij,    k,     n,      i,     j
    
    if ( oinit .or. ofinal ) then
       return
    end if

    do k = 1, nzdim
       do ij = 1, nxydim
          aa(ij, k) = 0.d0
          ab(ij, k) = 0.d0
          ac(ij, k) = 0.d0
       end do
    end do

    do ij = ijtstr, ijtend
       hxbot(ij) = hb(ij) + zbot
    end do

    do ij = 1, nxydim
       dh(ij) = hc(ij) - ha(ij)
    end do

    do k = kstr, kstr+kz-1
       do ij = ijtstr, ijtend
          aa(ij, k) = - ts * diffz(ij, k)   / dz(ij, k)
          ac(ij, k) = - ts * diffz(ij, k+1) / dz(ij, k)
          ab(ij, k) = (hc(ij) + zbot) / zbot - aa(ij, k) - ac(ij, k)
       end do
    end do

    do n = 3, ntdim
       do k = kstr, kstr+kz-1
          do ij = ijtstr, ijtend
             adt(ij, k, n) = adt(ij, k, n)                            &
    &                      - dh(ij) / zbot / ts * tx(ij, k, n)
          end do
       end do
    end do

    do k = kstr+kz, kend
       do ij = ijtstr, ijtend
          aa(ij, k) = - ts * diffz(ij, k  ) / dz(ij, k)
          ac(ij, k) = - ts * diffz(ij, k+1) / dz(ij, k)
          ab(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k)
       end do
    end do

    call thomas( adt, ac, aa, ab )

    do n = 3, ntdim
       do k = kstr, kend
          do ij = ijtstr, ijtend
             adt(ij, k, n) = adt(ij, k, n) * rgamma(k) * amskt(ij, k)
          end do
       end do
    end do
#ifdef OPT_BBL
    do n = 3, ntdim
       do ij = ijtstr, ijtend
          k = max(nbot(ij), kstr)
          adt(ij, k, n) = adt(ij, k, n) *                                  &
    &                (  1.d0                                               &
    &                 + amsktb(ij) * (  gamma(k-kstr+1) * rgamma(kend)     &
    &                                 - 1.d0) )
       end do
    end do
#endif

    do n = 3, ntdim
       do k = kstr, kend
          do ij = ijtstr, ijtend
             tx(ij, k, n) = tx(ij, k, n) + ts * adt(ij, k, n)
          end do
       end do
    end do

!---- tracer redistribution among sigma-layers (T. Suzuki) 
    do n = 3, ntdim
       do k = kstr, kstr+kz-1
          do ij = ijtstr, ijtend
             tx(ij, k, n) = tx(ij, k, n )                 &
    &                     * ( ( hc(ij) + zbot ) * ds(k) ) &
    &                     / ( ( hb(ij) + zbot ) * ds(k) )
          end do
       end do
    end do

!---- tracer redistribution among sigma-layers (COCO's default since F77 ver.)
!    do ij = ijtstr, ijtend
!       hzbot(ij) = hc(ij) + zbot
!       hxbot(ij) = hb(ij) + zbot
!       dh(ij) = (hb(ij) - hc(ij)) * amskt(ij, kstr)
!       if (dh(ij) .ge. 0.d0) then
!          dzb = 0.d0
!          dzt = dh(ij) * ds(kstr+kz-1)
!          do k = kstr+kz-1, kstr+1, -1
!             do n = 3, ntdim
!                tx(ij, k, n) = (  tx(ij, k  , n) * ds(k) * hzbot(ij) &
!     &                          - tx(ij, k  , n) * dzb               &
!     &                          + tx(ij, k-1, n) * dzt               &
!     &                         ) / hxbot(ij) / ds(k)
!             end do
!             dzb = dzt
!             dzt = dzb + dh(ij) * ds(k-1)
!          end do
!          do n = 3, ntdim
!             tx(ij, kstr, n) = ( tx(ij, kstr, n)*ds(kstr)*hzbot(ij) &
!     &                         - tx(ij, kstr, n) * dzb              &
!     &                         ) / hxbot(ij) / ds(kstr)
!          end do
!       else
!          dzt = - dh(ij)
!          dzb = dzt + dh(ij) * ds(kstr)
!          do n = 3, ntdim
!             tx(ij, kstr, n) = ( tx(ij, kstr, n)*ds(kstr)*hzbot(ij) &
!     &                         + tx(ij, kstr+1, n) * dzb            &
!     &                         ) / hxbot(ij) / ds(kstr)  
!          end do
!          do k = kstr+1, kstr+kz-2
!             dzt = dzb
!             dzb = dzt + dh(ij) * ds(k)
!             do n = 3, ntdim
!                tx(ij, k, n) = (  tx(ij, k, n) * ds(k) * hzbot(ij)  &
!     &                          + tx(ij, k+1, n) * dzb              &
!     &                          - tx(ij, k  , n) * dzt              &
!     &                         ) / hxbot(ij) / ds(k)
!             end do
!          end do
!       end if
!    end do
    
    do n = 3, ntdim
       do ij = ijtstr, ijtend
          tx(ij, kstr, n) = tx(ij, kstr, n)                          &
    &                     - ts * ft(ij, n) / hxbot(ij) / ds(kstr)
       end do
    end do

    do ij = 1, nxydim
       dh(ij) = - dh(ij) / ts
    end do
    call chekin(    dh,       'DIAGFW', &
         &         'diagnosed freshwater flux',   'cm/s', &
         &          nx,     ny,      1, nxydim, 'OCSFCT') 

  end subroutine slvtrc

! --- information -----------------------------------------------------
!
!  HISTORY
!     '04.04.13  A.Oka: from MIROC2.1
!     '07.11.28  Y.Komuro: for COCO4.3
!     '08.03.06  Y.Komuro: parallel code merged (from MIROC3.2)
!                          bug fix (LTNGG incorrect)
!
! ---------------------------------------------------------------------

  subroutine ttsset

    use zocdim,  only  :                                              &
         &    nxgdim, nygdim,                                         &
         &     igstr,  jgstr,                                         &
         &    nxydim, nxyzdm,  nxdim,  nzdim,  ntdim,                 &
         &    nxygdm,   nxyg,                                         &
         &        nx,     ny,                                         &
         &       istr,   jstr, ijstr,  ijend,                         &
         &    inodes, jnodes, ijnode,                                 &
         &      kstr,   kend,     kz
    use zocnod,  only  :  iroot,   nprocs,  myrank
    use zocgrd
    use zocmsk
    use zocphy
    use ufile
    use bgs2d

    implicit none

#include "mpif.h"

!    integer(4), save  ::  ntni,ntn
    integer(4), save  ::  ntn
    integer(4), save  ::  ltnig (ntunnel),       ltnigd(ntunnel)
    integer(4), save  ::  ltnigg(ntunnel)
    integer(4), save  ::    ltng(ntunnel*nzdim),  ltngd(ntunnel*nzdim)
    integer(4), save  ::   ltngg(ntunnel*nzdim)

    integer(4), save  ::   ltngn(ntunnel*nzdim),  ltngr(ntunnel*nzdim)
    integer(4), save  ::  ltngdn(ntunnel*nzdim), ltngdr(ntunnel*nzdim)

    real(8)           ::   buf2d (nxydim)
!    real(8)           ::   sndbfi(ntdim,  ntunnel*2+1)
!    real(8)           ::   rcvbfi(ntdim,  ntunnel*2+1)
!    real(8)           ::   sndbuf(ntdim1, ntunnel*nzdim*2+1)
!    real(8)           ::   rcvbuf(ntdim1, ntunnel*nzdim*2+1)
!    integer(4)        ::   istmpi(mpi_status_size)

!    integer(4), allocatable, save :: isreq(:), irreq(:)

!    integer(4), save  ::   ltns(ntunnel*nzdim*2+1)
!    integer(4), save  ::  ltnsd(ntunnel*nzdim*2+1)
!    integer(4), save  ::  ltnsg(ntunnel*nzdim*2+1)
!    integer(4), save  ::  ltnis(ntunnel*2+1)
!    integer(4), save  :: ltnisd(ntunnel*2+1)
!    integer(4), save  :: ltnisg(ntunnel*2+1)
!    integer(4), save  :: ntnist(0:nxyg), ntnsst(0:nxyg)
!    integer(4), save  ::   ntns,          ntnis

    real(8)           ::   botg(nxygdm)

!    real(8),    save  ::    vlm(nxyzdm),  area(nxydim)
!    real(8),    save  ::   crst(ntngr),  crsti(ntigr)

!    real(8)           ::     ad(ntunnel*2+1)
!    real(8)           ::     hd(ntunnel*2+1)
!    real(8)           ::   tshd(ntunnel*2+1, ntdim)
!    real(8)           ::     vd(ntunnel*nzdim*2+1)
!    real(8)           ::     td(ntunnel*nzdim*2+1, ntdim)
    
    integer(4)        ::     ij,    ijk,      k,      n
    integer(4)        ::      l,     la,     lg,   itrc
    integer(4)        ::  ifpar,  jfpar, iranks
    integer(4)        ::      i,      j,     ig,     jg,     m
    integer(4)        ::   ierr

    real(8),    save  :: dmptun(ntunnel), dshtun(ntunnel)
    integer(4), save  ::  ltun1(ntunnel),  ltun2(ntunnel),  ktun(ntunnel)
    integer(4), save  :: levtun(ntunnel)
    integer(4), save  ::  itun1(ntunnel),  jtun1(ntunnel)
    integer(4), save  ::  itun2(ntunnel),  jtun2(ntunnel)
    integer(4)        ::   ntun

    namelist /nmtunl/ itun1, jtun1, itun2, jtun2, levtun, dmptun, dshtun
    data levtun / ntunnel*0 /
    data dmptun / ntunnel*30.0d0 /
    data dshtun / ntunnel*30.0d0 /
    data  itun1 / ntunnel*0 /
    data  jtun1 / ntunnel*0 /
    data  itun2 / ntunnel*0 /
    data  jtun2 / ntunnel*0 /

    allocate(isreq(0:nprocs-1))
    allocate(irreq(0:nprocs-1))

    do ij=1,nxydim
       buf2d(ij) = nbot(ij)
    enddo
    call gather_2d(botg, buf2d)

!     --TUN DIF--
    call rewnml(ifpar, jfpar)
    read(ifpar, nmtunl, iostat = istat )
    call cstnml( jfpar, 'ttsset', 'nmtunl', istat )
    write(jfpar,*) '*** ttsset (tunnel diffusion) ***'

    if ( myrank == iroot ) then
       ntun = 0
       do l = 1, ntunnel
          if ( levtun(l) == 0 ) goto 1999
          ltun1(l) = igstr+itun1(l)-1 + nxgdim*(jgstr+jtun1(l)-2)
          ltun2(l) = igstr+itun2(l)-1 + nxgdim*(jgstr+jtun2(l)-2)
          ktun (l) = min( min( levtun(l), int(botg(ltun1(l))) ),      &
    &                     int(botg(ltun2(l))) )
          ntun = l
       end do
1999   continue
       ntni = ntun
       do l = 1, ntun
          ltnig (l) = ltun1(l)
          ltnigd(l) = ltun2(l)
          ltnigg(l) = l
       end do
       ntn = 0
       do l = 1, ntun
          do k = kstr+kz, kstr+ktun(l)-1
             ntn = ntn+1
             ltng (ntn) = ltun1(l) + nxygdm*(k-1)
             ltngd(ntn) = ltun2(l) + nxygdm*(k-1)
             ltngg(ntn) = l
          end do
       end do
       write(jfpar,*) '  TOTAL TUNNEL NUMBER = ',ntun
       if (ntun /= 0) then
          write(jfpar,*) '( I1, J1 ) .. ( I2, J2) KLEV DMP DMPSH'
          do l = 1, ntun
             write(jfpar,*) '(',itun1(l),',',jtun1(l),') .. (',          &
          &        itun2(l),',',jtun2(l),') ',ktun(l),dmptun(l),dshtun(l)
          enddo
       end if
    end if

    n = 1
    call mpi_bcast(ntni  ,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
    if ( ntni == 0 ) return

    if ( ntni > 0 ) then
       n = ntni
       call mpi_bcast(ltnig ,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
       call mpi_bcast(ltnigd,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
       call mpi_bcast(ltnigg,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
    end if
    n = 1
    call mpi_bcast(ntn   ,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
    if ( ntn > 0 ) then
       n = ntn
       call mpi_bcast(ltng  ,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
       call mpi_bcast(ltngd ,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
       call mpi_bcast(ltngg ,n,mpi_integer,iroot,mpi_comm_ogcm,ierr)
    end if

    do n = 1, ntn
       k  = (ltng(n) - 1) / nxygdm + 1
       jg = (ltng(n) - (k - 1) * nxygdm - 1) / nxgdim + 1
       ig = ltng(n) - (k - 1) * nxygdm - (jg - 1) * nxgdim
       j = mod(jg-jgstr, ny) + jstr
       i = mod(ig-igstr, nx) + istr
       ltngn(n) = (k - 1) * nxydim + (j - 1) * nxdim + i
       ltngr(n) = (jg - jgstr) / ny * inodes + (ig - igstr) / nx
    end do
    do n = 1, ntn
       k  = (ltngd(n) - 1) / nxygdm + 1
       jg = (ltngd(n) - (k - 1) * nxygdm - 1) / nxgdim + 1
       ig = ltngd(n) - (k - 1) * nxygdm - (jg - 1) * nxgdim
       j = mod(jg-jgstr, ny) + jstr
       i = mod(ig-igstr, nx) + istr
       ltngdn(n) = (k - 1) * nxydim + (j - 1) * nxdim + i
       ltngdr(n) = (jg - jgstr) / ny * inodes + (ig - igstr) / nx
    end do
    

    if ( myrank < ijnode ) then
       ntns = 0
       ntnsst(0) = 1
       do iranks = 0, ijnode-1
          if (iranks .eq. myrank) then
             do n = 1, ntn
                if (      (ltngr(n)  .eq. iranks)                     &
    &               .and. (ltngdr(n) .eq. iranks)) then
                   ntns = ntns + 1
                   ltns (ntns) = ltngn(n)
                   ltnsd(ntns) = ltngdn(n)
                   ltnsg(ntns) = ltngg(n)
                   ntns = ntns + 1
                   ltns (ntns) = ltngdn(n)
                   ltnsd(ntns) = ltngn(n)
                   ltnsg(ntns) = ltngg(n)
                end if
             end do
          else
             do n = 1, ntn
                if (      (ltngr(n)  .eq. myrank)                     &
    &               .and. (ltngdr(n) .eq. iranks)) then
                   ntns = ntns + 1
                   ltns (ntns) = ltngn(n)
                   ltnsd(ntns) = -1
                   ltnsg(ntns) = ltngg(n)
                else if (      (ltngr(n)  .eq. iranks)                &
    &                    .and. (ltngdr(n) .eq. myrank)) then
                   ntns = ntns + 1
                   ltns (ntns) = ltngdn(n)
                   ltnsd(ntns) = -1
                   ltnsg(ntns) = ltngg(n)
                end if
             end do
          end if
          ntnsst(iranks+1) = ntns + 1
       end do
    end if

    do n = 1, ntni
       jg = (ltnig(n) - 1) / nxgdim + 1
       ig = ltnig(n) - (jg - 1) * nxgdim
       j = mod(jg-jgstr, ny) + jstr
       i = mod(ig-igstr, nx) + istr
       ltngn(n) = (j - 1) * nxdim + i
       ltngr(n) = (jg - jgstr) / ny * inodes + (ig - igstr) / nx
    end do
    do n = 1, ntni
       jg = (ltnigd(n) - 1) / nxgdim + 1
       ig = ltnigd(n) - (jg - 1) * nxgdim
       j = mod(jg-jgstr, ny) + jstr
       i = mod(ig-igstr, nx) + istr
       ltngdn(n) = (j - 1) * nxdim + i
       ltngdr(n) = (jg - jgstr) / ny * inodes + (ig - igstr) / nx
    end do

    if ( myrank < ijnode ) then
       ntnis = 0
       ntnist(0) = 1
       do iranks = 0, ijnode-1
          if (iranks .eq. myrank) then
             do n = 1, ntni
                if (      (ltngr(n)  .eq. iranks)                     &
    &               .and. (ltngdr(n) .eq. iranks)) then
                   ntnis = ntnis + 1
                   ltnis (ntnis) = ltngn(n)
                   ltnisd(ntnis) = ltngdn(n)
                   ltnisg(ntnis) = ltnigg(n)
                   ntnis = ntnis + 1
                   ltnis (ntnis) = ltngdn(n)
                   ltnisd(ntnis) = ltngn(n)
                   ltnisg(ntnis) = ltnigg(n)
                end if
             end do
          else
             do n = 1, ntni
                if (      (ltngr(n)  .eq. myrank)                     &
    &               .and. (ltngdr(n) .eq. iranks)) then
                   ntnis = ntnis + 1
                   ltnis (ntnis) = ltngn(n)
                   ltnisd(ntnis) = -1
                   ltnisg(ntnis) = ltnigg(n)
                else if (      (ltngr(n)  .eq. iranks)                &
    &                    .and. (ltngdr(n) .eq. myrank)) then
                   ntnis = ntnis + 1
                   ltnis (ntnis) = ltngdn(n)
                   ltnisd(ntnis) = -1
                   ltnisg(ntnis) = ltnigg(n)
                end if
             end do
          end if
          ntnist(iranks+1) = ntnis + 1
       end do
    end if
    
    do n = 1, ntngr
       crst (n) = dmptun(n)
    end do
    do n =1, ntigr
       crsti(n) = dshtun(n)
    end do

    do k = kstr+kz, kend
       do ij = ijstr, ijend
          ijk = (k - 1) * nxydim + ij
          vlm(ijk) = hxt(ij) * hyt(ij) * dz(ij, k)
       end do
    end do

    do ij = ijstr, ijend
       area(ij) = hxt(ij) * hyt(ij)
    end do
    do n = 1, ntngr
       crst(n) = 1.d0 / crst(n) / 86400.d0
    end do
    do n = 1, ntigr
       crsti(n) = 1.d0 / crsti(n) / 86400.d0
    end do

  end subroutine ttsset

! =====================================================================

  subroutine tundif( tx, hx )

    use zocdim,  only :                                &
         &  nxydim,  nxyzdm,  nzdim,  ntdim,           &
         &    kstr,    kend,     kz,                   &
         &   ijstr,   ijend,                           &
         &  inodes,  jnodes, ijnode,                   &
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dz,    dz0,     ds,   zbot,            &
         &      ts
    use zocnod,  only  :  myrank

    implicit none

#include "mpif.h"

    real(8),   intent(inout)  ::     tx(nxyzdm, ntdim)
    real(8),   intent(inout)  ::     hx(nxydim)

!---- local variables
    real(8)                   :: sndbfi(ntdim,ntunnel*2+1)
    real(8)                   :: rcvbfi(ntdim,ntunnel*2+1)
    real(8)                   :: sndbuf(ntdim1,ntunnel*nzdim*2+1)
    real(8)                   :: rcvbuf(ntdim1,ntunnel*nzdim*2+1)
    integer(4)                :: istmpi(mpi_status_size)
    real(8)                   :: rhxbot(nxydim)
    real(8)                   ::   tshd(ntunnel*2+1,ntdim)
    real(8)                   ::     vd(ntunnel*nzdim*2+1)
    real(8)                   ::     td(ntunnel*nzdim*2+1, ntdim)
    real(8)                   ::     ad(ntunnel*2+1)
    real(8)                   ::     hd(ntunnel*2+1)
    real(8)                   ::    tsh(nxydim, nzdim, ntdim)
    integer(4)                ::     ij,     k,      n,     m
    integer(4)                ::     lg,     l
    integer(4)                ::    ijk
    integer(4)                :: iranks,  nbuf,   ierr

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ntni == 0 ) return

    do n = 3, ntdim
       do k = kstr, kstr+kz-1
          do ij = ijstr, ijend
             ijk = (k - 1) * nxydim + ij
             tsh(ij, k, n) = tx(ijk, n) * (hx(ij) + zbot)
          end do
       end do
    end do

    do iranks = 0, ijnode-1
       if ( iranks .ne. myrank ) then
          if (ntnist(iranks+1) .gt. ntnist(iranks)) then
             nbuf = ntnist(iranks+1) - ntnist(iranks)
             do n = ntnist(iranks), ntnist(iranks+1)-1
                sndbfi(1, n) = hx(ltnis(n))
                sndbfi(2, n) = area(ltnis(n))
             end do
             call mpi_isend(                                          &
    &                        sndbfi(1, ntnist(iranks)),               &
    &                        ntdim*nbuf, mpi_double_precision,        &
    &                         iranks,      1, mpi_comm_ogcm,         &
    &                         isreq(iranks),   ierr)
             call mpi_irecv(                                          &
    &                        rcvbfi(1, ntnist(iranks)),               &
    &                        ntdim*nbuf, mpi_double_precision,        &
    &                         iranks,      1, mpi_comm_ogcm,         &
    &                         irreq(iranks),   ierr)
          end if
       end if
    end do

    do iranks = 0, ijnode-1
       if (iranks .ne. myrank) then
          if (ntnist(iranks+1) .gt. ntnist(iranks)) then
             call mpi_wait(isreq(iranks), istmpi, ierr)
             call mpi_wait(irreq(iranks), istmpi, ierr)
          end if
       end if
    end do
    do iranks = 0, ijnode-1
       if (iranks .eq. myrank) then
          do n = ntnist(iranks), ntnist(iranks+1)-1
             hd(n) = hx(ltnisd(n))
             ad(n) = area(ltnisd(n))
          end do
       else
          do n = ntnist(iranks), ntnist(iranks+1)-1
             hd(n) = rcvbfi(1, n)
             ad(n) = rcvbfi(2, n)
          end do
       end if
    end do
!*POPTION INDEP(HX)
    do n = 1, ntnis
       l  = ltnis (n)
       lg = ltnisg(n)
       hx(l) = hx(l)                                                  &
    &         - (hx(l) - hd(n)) * ts * crsti(lg) *                    &
    &           (area(l) + ad(n)) * 0.5d0 / area(l)
    end do

    do ij = ijstr, ijend
       rhxbot(ij) = 1.d0 / (hx(ij) + zbot)
    end do

    do k = kstr, kstr+kz-1

       do iranks = 0, ijnode-1
          if (iranks .ne. myrank) then
             if (ntnist(iranks+1) .gt. ntnist(iranks)) then
                nbuf = ntnist(iranks+1) - ntnist(iranks)
                do n = ntnist(iranks), ntnist(iranks+1)-1
                   do l = 1, ntdim
                      sndbfi(l, n) = tsh(ltnis(n), k, l)
                   end do
                end do
                call mpi_isend(                                       &
    &                       sndbfi(1, ntnist(iranks)),                &
    &                       ntdim*nbuf, mpi_double_precision,         &
    &                        iranks,      1, mpi_comm_ogcm,          &
    &                        isreq(iranks),   ierr)
                call mpi_irecv(                                       &
    &                       rcvbfi(1, ntnist(iranks)),                &
    &                       ntdim*nbuf, mpi_double_precision,         &
    &                        iranks,      1, mpi_comm_ogcm,          &
    &                        irreq(iranks),   ierr)
             end if
          end if
       end do

       do iranks = 0, ijnode-1
          if (iranks .ne. myrank) then
             if (ntnist(iranks+1) .gt. ntnist(iranks)) then
                call mpi_wait(isreq(iranks), istmpi, ierr)
                call mpi_wait(irreq(iranks), istmpi, ierr)
             end if
          end if
       end do

       do iranks = 0, ijnode-1
          if (iranks .eq. myrank) then
             do l = 1, ntdim
                do n = ntnist(iranks), ntnist(iranks+1)-1
                   tshd(n, l) = tsh(ltnisd(n), k, l)
                end do
             end do
          else
             do n = ntnist(iranks), ntnist(iranks+1)-1
                do l = 1, ntdim
                   tshd(n, l) = rcvbfi(l, n)
                end do
             end do
          end if
       end do

       do m = 1, ntdim
!*POPTION INDEP(TSH)
          do n = 1, ntnis
             l  = ltnis (n)
             lg = ltnisg(n)
             tsh(l, k, m) = tsh(l, k, m)                              &
    &                     - (tsh(l, k, m) - tshd(n, m)) * ts *        &
    &                       crsti(lg) * (area(l) + ad(n)) * 0.5d0     &
    &                       / area(l)
          end do
          do ij = ijstr, ijend
             ijk = (k - 1) * nxydim + ij
             tx(ijk, m) = tsh(ij, k, m) * rhxbot(ij)
          end do
       end do
       
    end do

    do iranks = 0, ijnode-1
       if (iranks .ne. myrank) then
          if (ntnsst(iranks+1) .gt. ntnsst(iranks)) then
             nbuf = ntnsst(iranks+1) - ntnsst(iranks)
             do n = ntnsst(iranks), ntnsst(iranks+1)-1
                do l = 1, ntdim
                   sndbuf(l, n) = tx(ltns(n), l)
                end do
                sndbuf(ntdim1, n) = vlm(ltns(n))
             end do
             call mpi_isend(                                          &
    &                        sndbuf(1, ntnsst(iranks)),               &
    &                        ntdim1*nbuf, mpi_double_precision,       &
    &                         iranks,      1, mpi_comm_ogcm,         &
    &                         isreq(iranks),   ierr)
             call mpi_irecv(                                          &
    &                        rcvbuf(1, ntnsst(iranks)),               &
    &                        ntdim1*nbuf, mpi_double_precision,       &
    &                         iranks,      1, mpi_comm_ogcm,         &
    &                         irreq(iranks),   ierr)
          end if
       end if
    end do

    do iranks = 0, ijnode-1
       if (iranks .ne. myrank) then
          if (ntnsst(iranks+1) .gt. ntnsst(iranks)) then
             call mpi_wait(isreq(iranks), istmpi, ierr)
             call mpi_wait(irreq(iranks), istmpi, ierr)
          end if
       end if
    end do

    do iranks = 0, ijnode-1
       if (iranks .eq. myrank) then
          do n = ntnsst(iranks), ntnsst(iranks+1)-1
             do l = 1, ntdim
                td(n, l) = tx(ltnsd(n), l)
             end do
             vd(n) = vlm(ltnsd(n))
          end do
       else
          do n = ntnsst(iranks), ntnsst(iranks+1)-1
             do l = 1, ntdim
                td(n, l) = rcvbuf(l, n)
             end do
             vd(n) = rcvbuf(ntdim1, n)
          end do
       end if
    end do

    do m = 3, ntdim
!*POPTION INDEP(TX)
       do n = 1, ntns
          l  = ltns (n)
          lg = ltnsg(n)
          tx(l, m) = tx(l, m)                                         &
    &              - (tx(l, m) - td(n, m)) * ts * crst(lg) *          &
    &                (vlm(l) + vd(n)) * 0.5d0 / vlm(l)
       end do
    end do

  end subroutine tundif
  
end module tslvt
