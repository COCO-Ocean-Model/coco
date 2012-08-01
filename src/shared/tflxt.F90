module tflxt

! --- information -----------------------------------------------------
!
!  Estimate the advection and diffusion terms of the tracer equations.
!
!  HISTORY
!     '03.04.23  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.08.01  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim, nxydim, nxyzdm, nzdim, ntdim, &
    &   kstr,   kend,     kz, &
    &     nx,     ny,     nz, &
    & ijtstr, ijtend, &
    &     le,     lw,     ln,     ls,    lsw, &
    &  oinit, ofinal    
  use zocgrd, only: &
    &     dz,    dzm,    dzv,    dsm, &
    &     rx,     ry,    rym, &
    &   zbot, &
    &    hxt,    hxu,    hyt,    hyu,    rxt,    ryt
  use zocmsk, only: &
#ifdef OPT_BBL
    & amsktb, &
#endif
    &  amftx,  amfty,  amftz, &
    &   nbot

  implicit none

  private

  real(8) ::    ftx(nxydim, nzdim, ntdim)
  real(8) ::    fty(nxydim, nzdim, ntdim)
  real(8) ::    ftz(nxydim, nzdim, ntdim)

  public :: flxtrc, chkftx
#ifdef OPT_BBL
  public :: flxtrb
#endif

contains

subroutine flxtrc( &
  &                   adt,  diffz, &
  &                    tx,     hx, &
  &                    ty,     hz, &
  &                    uy,     vy, &
  &                     w,    ahv )

  real(8), intent(out) ::    adt(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  diffz(nxydim, nzdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     hx(nxydim),     hz(nxydim)
  real(8), intent(in)  ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::      w(nxydim, nzdim),    ahv(nxydim, nzdim)

  logical, save :: ofirst = .true.
  real(8), save ::    alp,   alpv

  real(8) ::  hzbot(nxydim)
  real(8) ::      u,      v,     wt
  integer ::     ij,      k,      n
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::  ahh = 0.0d0,  alpha = 0.5d0,  alphav = 0.5d0

  namelist /nmdifh/ ahh
  namelist /nmwupc/ alpha, alphav

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifh, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifh', istat)
     write(jfpar, nmdifh)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmwupc, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmwupc', istat)
     write(jfpar, nmwupc)

     alp  = (alpha - 0.5d+0) * 0.5d+0
     alpv = alphav - 0.5d+0
  end if

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           ftx(ij, k, n) = 0.d0
           fty(ij, k, n) = 0.d0
           ftz(ij, k, n) = 0.d0
           adt(ij, k, n) = 0.d0
        end do
     end do
  end do
  do k = 1, nzdim
     do ij = 1, nxydim
        diffz(ij, k) = 0.d0
     end do
  end do
  do ij = 1, nxydim
     hzbot(ij) = hz(ij) + zbot
  end do

  do k = kstr+1, kstr+kz-1
     do ij = ijtstr, ijtend
        diffz(ij, k) = ahv(ij, k) / dsm(k) / hzbot(ij) * &
          &            amftz(ij, k)
     end do
  end do
  do n = 1, ntdim
     do k = kstr+1, kstr+kz-1
        do ij = ijtstr, ijtend
           wt = sign(alpv, w(ij, k))
           ftz(ij, k, n) =  &
             &           (  diffz(ij, k) * (tx(ij, k-1, n) - tx(ij, k, n)) &
             &            - w(ij, k) * hzbot(ij) * &
             &              (  (0.5d0 - wt) * tx(ij, k-1, n) &
             &               + (0.5d0 + wt) * tx(ij, k  , n))) * &
             &           amftz(ij, k)
        end do
     end do
  end do

  do k = kstr+kz, kend
     do ij = ijtstr, ijtend
        diffz(ij, k) = ahv(ij, k) / dzm(ij, k) * amftz(ij, k)
     enddo
  enddo
  do n = 1, ntdim
     do k = kstr+kz, kend
        do ij = ijtstr, ijtend
           wt = sign(alpv, w(ij, k))
           ftz(ij, k, n) = &
             &           (  diffz(ij, k) * (tx(ij, k-1, n) - tx(ij, k, n)) &
             &            - w(ij, k) * (  (0.5d0 - wt) * tx(ij, k-1, n) &
             &                          + (0.5d0 + wt) * tx(ij, k,   n))) * &
             &           amftz(ij, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend

        do ij = ijtstr, ijtend+nxdim
           v   = vy(ij+ls , k) * dzv(ij+ls , k) * hxu(ij+ls) &
             & + vy(ij+lsw, k) * dzv(ij+lsw, k) * hxu(ij+lsw)
           wt = sign(alp, v)
           fty(ij, k, n) = &
             &     ahh * (tx(ij, k, n) - tx(ij+ls, k, n)) * rym(ij+ls) * &
             &     (hxu(ij+ls) + hxu(ij+lsw)) / (hyt(ij) + hyt(ij+ls)) * &
             &     amfty(ij, k) &
             &   - v * (  (0.25d0 - wt) * tx(ij, k,    n) &
             &          + (0.25d0 + wt) * tx(ij+ls, k, n)) 
        end do

        do ij = ijtstr, ijtend+1
           u   = uy(ij+lw , k) * dzv(ij+lw , k) * hyu(ij+lw) &
             & + uy(ij+lsw, k) * dzv(ij+lsw, k) * hyu(ij+lsw)
           wt = sign(alp, u)
           ftx(ij, k, n) = &
             &     ahh * (tx(ij, k, n) - tx(ij+lw, k, n)) * rx * &
             &     (hyu(ij+lw) + hyu(ij+lsw)) / (hxt(ij) + hxt(ij+lw)) * &
             &     amftx(ij, k) &
             &   - u * (  (0.25d0 - wt) * tx(ij, k, n) &
             &          + (0.25d0 + wt) * tx(ij+lw, k, n))
        end do

        do ij = ijtstr, ijtend
           adt(ij, k, n) = &
             &     (  (  (ftx(ij+le, k, n) - ftx(ij, k, n)) * rx &
             &         + (fty(ij+ln, k, n) - fty(ij, k, n)) * ry(ij)) * &
             &        rxt(ij) * ryt(ij) &
             &      + ftz(ij, k, n) - ftz(ij, k+1, n)) / dz(ij, k)
        end do

     end do
  end do

  return

end subroutine flxtrc

#ifdef OPT_BBL
! *********************************************************************

subroutine flxtrb( &
  &                   adt,  diffz, &
  &                    tx, &
  &                    ty, &
  &                    uy,     vy,      w, &
  &                   ahv )

! --- information -----------------------------------------------------
!
!  Estimate the advection and diffusion terms of the tracer equations
! for BBL.
!
!  HISTORY
!     '01.02.13  H.Hasumi
!     '01.05.02  H.Hasumi
!     '02.05.29  H.Nakano: tracer dimension
!     '12.08.01  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  real(8), intent(out) ::    adt(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  diffz(nxydim, nzdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::      w(nxydim, nzdim),    ahv(nxydim, nzdim)

  real(8), save ::    alp,   alpv
  logical, save :: ofirst = .true.

  real(8) ::      u,      v,     wt
  integer ::     ij,      k,    kup,      n
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::  ahhbbl = 0.0d0,  alpha = 0.5d0,  alphav = 0.5d0

  namelist /nmbbdh/ ahhbbl
  namelist /nmwupc/ alpha, alphav

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmbbdh, iostat=istat)
     call cstnml(jfpar, 'flxtrb', 'nmbbdh', istat)
     write(jfpar, nmbbdh)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmwupc, iostat=istat)
     call cstnml(jfpar, 'flxtrb', 'nmwupc', istat)
     write(jfpar, nmwupc)

     alp  = (alpha - 0.5d0) * 0.5d0
     alpv = alphav - 0.5d0
  end if

  do n = 1, ntdim

     do ij = ijtstr, ijtend
        k = nbot(ij)
        kup = max(k-1, 1)
        wt = sign(alpv, w(ij, kend))
        ftz(ij, kend, n) = &
          &      (  diffz(ij, k) * (tx(ij, kup, n) - tx(ij, kend, n)) &
          &       - w(ij, kend) * (  (0.5d0 - wt) * tx(ij, kup, n) &
          &                        + (0.5d0 + wt) * tx(ij, kend, n))) * &
          &      amsktb(ij)
     end do

     do ij = ijtstr, ijtend+nxdim
        v   = vy(ij+ls , kend) * dzv(ij+ls , kend) * hxu(ij+ls) &
          & + vy(ij+lsw, kend) * dzv(ij+lsw, kend) * hxu(ij+lsw)
        wt = sign(alp, v)
        fty(ij, kend, n) = &
          &        ahhbbl * (tx(ij, kend, n) - tx(ij+ls, kend, n)) * &
          &        rym(ij) * &
          &        (hxu(ij+ls) + hxu(ij+lsw)) / (hyt(ij) + hyt(ij+ls)) * &
          &        amfty(ij, kend) &
          &      - v * (  (0.25d0 - wt) * tx(ij, kend, n) &
          &             + (0.25d0 + wt) * tx(ij+ls, kend, n))
     end do

     do ij = ijtstr, ijtend+1
        u   = uy(ij+lw , kend) * dzv(ij+lw , kend) * hyu(ij+lw) &
          & + uy(ij+lsw, kend) * dzv(ij+lsw, kend) * hyu(ij+lsw)
        wt = sign(alp, u)
        ftx(ij, kend, n) = &
          &        ahhbbl * (tx(ij, kend, n) - tx(ij+lw, kend, n)) * rx * &
          &        (hyu(ij+lw) + hyu(ij+lsw)) / (hxt(ij) + hxt(ij+lw)) * &
          &        amftx(ij, kend) &
          &      - u * (  (0.25d0 - wt) * tx(ij, kend, n) &
          &             + (0.25d0 + wt) * tx(ij+lw, kend, n))
     end do

     do ij = ijtstr, ijtend
        adt(ij, kend, n) = &
          &      (  (  (  ftx(ij+le, kend, n) &
          &             - ftx(ij, kend, n)) * rx &
          &          + (  fty(ij+ln, kend, n) &
          &             - fty(ij, kend, n)) * ry(ij)) * rxt(ij) * ryt(ij) &
          &       + ftz(ij, kend, n)) / dz(ij, kend)
     end do

  end do

  return

end subroutine flxtrb

#endif
! *********************************************************************
subroutine chkftx

  if (oinit .or. ofinal) then
     return
  end if

  call chekin(   ftx,  'FTX', &
    &             nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   fty,  'FTY', &
    &             nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftz,  'FTZ', &
    &             nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftx(1, 1, 2),  'FSX', &
    &             nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   fty(1, 1, 2),  'FSY', &
    &             nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftz(1, 1, 2),  'FSZ', &
    &             nx,     ny,     nz, nxyzdm, 'OCN')

  return

end subroutine chkftx

end module tflxt
