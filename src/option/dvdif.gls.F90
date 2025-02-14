module dvdif

! --- information -----------------------------------------------------
!
!  Sea surface mixed layer model of Umlauf and Burchard (2003, JMR).
!  The implementation is based on Warner et al. (2005, OM).
!  Effects of surface wave breaking is also introduced 
!    after Carniel et al. (2009, OM).
!
!  HISTORY
!     '09.10.23  Y.Komuro: for COCO4.5
!     '10.04.09  Y.Komuro: many bug fixes + surface wave breaking
!     '10.04.13  Y.Komuro: bug fix (for non-BBL use etc.)
!     '12.01.19  Y.Komuro: bug fix (GH in unstable situation)
!     '12.07.11  Y.Komuro: for COCO5.0
!     '13.02.13  Y.Komuro: remove non-parallel code 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim,  ntdim, nxydim, nxyzdm, &
    &   kstr,   kend,     kz,     nz, &
    &  ijstr,  ijend, ijvstr, ijvend, &
    &     le,     lw,     ln,     ls,    lne,    lsw, &
    &  oinit, ofinal, &
#ifndef OPT_IO_COCOMPI
    &    nxg,    nyg, nxgdim, nygdim,  igstr,  jgstr, &
#endif
    &     nx,     ny
  use zocfil, only: &
    &    ncf
  use zocgrd, only: &
    &     dz,    dzm,     ds,    dsm,     dt, &
    &   zbot,    cor,   itst, ieuler
  use zocmsk,  only: &
#ifdef OPT_BBL
    &  amsktb, amskvb, nbotv, &
#endif
    &  amskt,  amskv,  amskb,  amftz,  amfvz, &
    &   nbot
  use zocphy, only: &
    & rhoo, gravit,  ckarm

  implicit none
  private
#include "mpif.h"

  real(8), save :: tauaox(nxydim), tauaoy(nxydim)

  public ::  vdiff, puttao
#ifdef OPT_BBL
  public :: vdiffb
#endif

contains

subroutine vdiff( &
  &                  amv,    ahv, &
  &                   uy,     vy,      r,   taux,   tauy, &
  &                   ty,     hy)

  use xprst
  use brstt
  use ufile
  use bchmk
  use qckot
  use bshft
#ifdef OPT_IO_COCOMPI
  use mpiio
#else
  use bgs2d
  use zocnod,  only :  iroot,  myrank
#endif

  real(8), intent(out) ::     amv(nxydim, nzdim),    ahv(nxydim, nzdim)
  real(8), intent(in)  ::      uy(nxydim, nzdim),     vy(nxydim, nzdim)

  real(8), intent(in)  ::       r(nxydim, nzdim)
  real(8), intent(in)  ::    taux(nxydim)       ,   tauy(nxydim)
  real(8), intent(in)  ::      ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::      hy(nxydim)

  real(8), save ::    tke(nxydim, nzdim),    psi(nxydim, nzdim)
  real(8), save :: ahv03d(nxydim, nzdim)
  real(8), save ::  csamv(nxydim, nzdim),  csahv(nxydim, nzdim)

  real(8) ::   amvt(nxydim, nzdim)
  real(8) ::   drdz(nxydim, nzdim),  duvdz(nxydim, nzdim)
  real(8) ::  dzsig(nxydim, nzdim), rzmsig(nxydim, nzdim)
  real(8) ::  depth(nxydim, nzdim), depthm(nxydim, nzdim)
  real(8) ::  epsil(nxydim, nzdim),    tls(nxydim, nzdim)
  real(8) ::     sh(nxydim, nzdim),     sm(nxydim, nzdim)
  real(8) ::  fwall(nxydim, nzdim)
  real(8) ::    fez(nxydim, nzdim),    rit(nxydim, nzdim)
  real(8) ::     aa(nxydim, nzdim),     ab(nxydim, nzdim)
  real(8) ::     ac(nxydim, nzdim),    ade(nxydim, nzdim)
  real(8) :: adefwd(nxydim, nzdim),   tke0(nxydim, nzdim)
  real(8) :: cdkdtl(nxydim, nzdim),  diffz(nxydim, nzdim)
  real(8) ::   cdmp(nxydim, nzdim)
  real(8) :: ufrc2s(nxydim),        ufrc2b(nxydim)
  real(8) :: ufrc2o(nxydim),        ufrc2i(nxydim)
  real(8) :: ufrc3s(nxydim),        ufrc3o(nxydim)
  real(8) :: taubtm(nxydim)
  real(8) ::     gh(nxydim, nzdim),   ghul(nxydim, nzdim)
  real(8) :: z0sf2d(nxydim),        prepmx(nxydim)
  real(8) ::  dpsi0(nxydim), scnp3d(nxydim, nzdim)
  real(8) :: ctkemn(nxydim, nzdim), ctfilt(nxydim, nzdim)
!      common /work/ dzsig, rzmsig, drdz, duvdz, &
!     &              depth, depthm, scnp3d, fwall

  real(8), save :: cc0(nzdim), cc1(nzdim), cc2(nzdim)
  real(8), save :: cc3(nzdim), cc4(nzdim), cc5(nzdim), cc6(nzdim)
  real(8), save :: d0(nzdim), d1(nzdim), d2(nzdim), d3(nzdim), d4(nzdim)
  real(8), save :: d5(nzdim), d6(nzdim), d7(nzdim), d8(nzdim), d9(nzdim)
  real(8), save :: r0

  logical, save :: ofirst = .true.,   oeof

  integer, save ::   nitr
  real(8), save :: csfsh1, csfsh2, csfsm1, csfsm2, csfsm3
  real(8), save ::   ceps,  eepst,  eepsp,   ctls
  real(8), save :: cpsi30, cpsi3s, cpsife, cpsibs, cpsibb
  real(8), save :: cpsilm, epsilm
  real(8), save :: scnpss,  cscnr, csftkw, cpsifw
  real(8), save :: ctaubt

  real(8) ::   cort,  cor30,     pi,  omega,    bfq
  real(8) ::   dudz,   dvdz
  real(8) ::     p1,     p2
  real(8) ::     tl,     sl
  real(8) ::     rl,    rlu
  real(8) ::      q,      p,   alps
  real(8) ::     fc,  avrtx,  avrty, avrtox, avrtoy,  psilm
  real(8) :: tketmp, dstwal,  dztmp,  rscnp
  real(8) ::  preps, preps0, alphcl
  integer ::     ij,      k,   iitr
  integer ::      i,      j
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::  amv0(nz),  ahv0(nz)
  real(8), save ::  rahv(nz)
  real(8), save ::  btmfrc = 0.0d0
  real(8), save ::  estrm = 1.0d0,  estrn = -0.6666666666666666d0
  real(8), save ::  estrp = 2.0d0
  real(8), save ::  csfa1 = 0.92d0,  csfa2 = 0.74d0
  real(8), save ::  csfb1 = 16.6d0,  csfb2 = 10.1d0
  real(8), save ::  csfc2 =  0.7d0,  csfc3 =  0.2d0
  real(8), save ::  sfghmn = -0.28d0,  sfgh0 = 0.0233d0,  sfghcr = 0.01d0
  real(8), save ::  csfe = 1.0d0,  cmu0 = 0.5544d0
  real(8), save ::  cpsi1 = 1.0d0,  cpsi2 = 1.22d0
  real(8), save ::  cpsi3p = 1.0d0,  cpsi3n = 0.1d0
  real(8), save ::  amvmax = 1000.0d0,  ahvmax = 1000.0d0
  real(8), save ::  amvmin = 0.0d0,  ahvmin = 0.0d0
  real(8), save ::  scntke = 0.8d0,  scnpsi = 1.07d0
  real(8), save ::  tkemin = 7.6d-2,  psimin = 1.0d-12, ritc = 1.0d0
  real(8), save ::  z0sfc = 1.0d2,  z0btm = 1.0d2,  epscmp = 1.0d-15
  real(8), save ::  cw = 100.0d0,  alphch = 5.6d4,  z0sfmn = 1.0d2
  real(8), save ::  alphci = 1.4d3,  atfilt = 0.0d0
  integer, save ::  nitr0 = 1,  mz = nz
  integer, save ::  nitr00 = 20
  logical, save ::  osfcwv = .false.,  oswnoi = .false.,  obtkei = .false.

  logical, save ::  ovdfao = .false.
  real(8), save ::  ahv0ao(nz) = 0.d0
  integer, save ::  mzao = 0
  real(8)       ::  corao
  
  real(8), save ::  tedn2d(nxydim)
  real(8), save ::  tedf2d(nxydim)
  real(8), save ::  tedn3d(nxydim, nzdim)
  real(8), save ::  tedf3d(nxydim, nzdim)
  real(8), save ::  cgamma = 0.2d0,  ahvemx = 1000.0d0,  epst = 1.d-20
  real(8), save ::    zeta = 500.d0 ! [m]
  logical, save ::  ofvcnt = .false., ofvpn = .false. ! vert. prof of far-field mixing
  ! ofvcnt = .true. : vertically constant
  ! ofvcnt = .false, ofvpn = .false. : prop. to N^2
  ! ofvcnt = .false, ofvpn = .true. : prop. to N
  real(8), save ::   rzeta ! [1/cm]
  integer, save ::  iamn = 0, iamf = 0
  character(len = ncf) ::  cftedn = 'not-specified'
  character(len = ncf) ::  cftedf = 'not-specified'
  character(len = 16)  ::  chead(1:64)
#ifdef OPT_IO_COCOMPI
  integer :: mpi_fh
  integer :: icread
  integer (kind = mpi_offset_kind) :: disp
#else
  integer :: nfted
  real(8), allocatable :: buf2(:, :),  g2d(:, :)
#endif
  real(8), save ::  depth0(nxydim, nzdim)
  real(8)       ::  dzmsig(nxydim, nzdim)
  real(8)       ::    gint(nxydim)
  real(8)       ::  ahvted(nxydim, nzdim),  tedr(nxydim, nzdim)
  real(8)       ::     dep ! [cm]
  
  namelist /nmvisv/ amv0
  namelist /nmdifv/ ahv0
  namelist /nmdfre/ rahv
  namelist /nmbtmf/ btmfrc
  namelist /nmdvgl/ estrm, estrn, estrp, &
    &               csfa1, csfa2, csfb1, csfb2, csfc2, csfc3, &
    &               sfghmn, sfgh0, sfghcr, csfe, &
    &               cmu0, cpsi1, cpsi2, cpsi3p, cpsi3n, &
    &               amvmax, ahvmax, scntke, scnpsi, &
    &               tkemin, psimin, ritc, z0sfc, z0btm, &
    &               nitr0, nitr00, epscmp, mz, &
    &               osfcwv, cw, z0sfmn, alphch, oswnoi, alphci, &
    &               obtkei, atfilt, amvmin, ahvmin
  namelist /nmdifvao/ ovdfao, ahv0ao, mzao
  namelist /nmdved/ iamn, iamf, cftedn, cftedf, cgamma, ahvemx, epst, zeta, ofvcnt, ofvpn
 
  if (oinit) then
     !$acc enter data create(tke,  psi)
     !$acc enter data create( ahv03d)
     !$acc enter data create(  csamv,  csahv)
     
     !$acc enter data create(   amvt)
     !$acc enter data create(   drdz,  duvdz)
     !$acc enter data create(  dzsig, rzmsig)
     !$acc enter data create(  depth, depthm)
     !$acc enter data create(  epsil,    tls)
     !$acc enter data create(     sh,     sm)
     !$acc enter data create(  fwall)
     !$acc enter data create(    fez,    rit)
     !$acc enter data create(     aa,     ab)
     !$acc enter data create(     ac,    ade)
     !$acc enter data create( adefwd,   tke0)
     !$acc enter data create( cdkdtl,  diffz)
     !$acc enter data create(   cdmp)
     !$acc enter data create( ufrc2s,        ufrc2b)
     !$acc enter data create( ufrc2o,        ufrc2i)
     !$acc enter data create( ufrc3s,        ufrc3o)
     !$acc enter data create( taubtm)
     !$acc enter data create(     gh,   ghul)
     !$acc enter data create( z0sf2d,        prepmx)
     !$acc enter data create(  dpsi0, scnp3d)
     !$acc enter data create( ctkemn, ctfilt)

     !$acc enter data create( cc0, cc1, cc2 )
     !$acc enter data create( cc3, cc4, cc5, cc6 )
     !$acc enter data create( d0, d1, d2, d3, d4 )
     !$acc enter data create( d5, d6, d7, d8, d9 )
     
     !$acc enter data create( amv0, ahv0)
     !$acc enter data create( rahv)
     
     !$acc enter data create( tedn2d)
     !$acc enter data create( tedf2d)
     !$acc enter data create( tedn3d)
     !$acc enter data create( tedf3d)
     
     !$acc enter data create( depth0)
     !$acc enter data create( dzmsig)
     !$acc enter data create(gint)
     !$acc enter data create(ahvted, tedr)
     
     do k = 1, nzdim
        do ij = 1, nxydim
           tke(ij, k) = tkemin
           psi(ij, k) = psimin
        end do
     end do
#ifdef OPT_TRIPOLE
     call rstadd(tke, oeof, nxdim, nydim, nzdim, 'TKE', 'OCN', &
       &                                          1.d0,  0,  0 )
     call rstadd(psi, oeof, nxdim, nydim, nzdim, 'PSI', 'OCN', &
       &                                          1.d0,  0,  0 )
#else
     call rstadd(tke, oeof, nxdim, nydim, nzdim, 'TKE', 'OCN')
     call rstadd(psi, oeof, nxdim, nydim, nzdim, 'PSI', 'OCN')
#endif
     !$acc update device(tke, psi)
     return
  end if

  if (ofinal) then
     !$acc update self(tke, psi)
     call finadd(tke, nxdim, nydim, nzdim, 'TKE', 'OCN')
     call finadd(psi, nxdim, nydim, nzdim, 'PSI', 'OCN')
     return
  end if

  if (ofirst) then
     !$acc kernels default(present)
     amv0(:)=0.d0
     ahv0(:)=0.d0
     rahv(:)=1.d0
     tedn3d(:,:)=0.d0
     tedf3d(:,:)=0.d0
     !$acc end kernels

!     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmvisv, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmvisv', istat)
     write(jfpar, nmvisv)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmdifv, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdifv', istat)
     write(jfpar, nmdifv)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmbtmf, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmbtmf', istat)
     write(jfpar, nmbtmf)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmdvgl, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdvgl', istat)
     write(jfpar, nmdvgl)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmdfre, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdfre', istat)
     write(jfpar, nmdfre)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmdifvao, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdifvao', istat)
     write(jfpar, nmdifvao)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmdved, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdved', istat)
     write(jfpar, nmdved)
     !$acc update device(amv0, ahv0)
     
     if (oeof) then
        !$acc kernels default(present)
        do k = 1, nzdim
           do ij = 1, nxydim
              tke(ij, k) = tkemin
              psi(ij, k) = psimin
           end do
        end do
        !$acc end kernels
     end if

     call secoef( &
       &          cc0(kstr), cc1(kstr), cc2(kstr), cc3(kstr), &
       &          cc4(kstr), cc5(kstr), cc6(kstr), &
       &          d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
       &          d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr))
     !$acc update device( cc0, cc1, cc2 )
     !$acc update device( cc3, cc4, cc5, cc6 )
     !$acc update device( d0, d1, d2, d3, d4 )
     !$acc update device( d5, d6, d7, d8, d9 )
     
     r0 = rhoo * 1.d3

     if (.not. osfcwv) then
        cw = 0.0d0
     end if

     csfsh1 = csfa2 * (1.0d0 - 6.0d0 * csfa1 / csfb1)
     csfsh2 = 3.0d0 * csfa2 * &
       &      (6.0d0 * csfa1 + csfb2 * (1.0d0 - csfc3))
     csfsm1 = csfb1**(-1.0d0/3.0d0)
     csfsm2 = 18.0d0 * csfa1 * csfa1 &
       &      + 9.0d0 * csfa1 * csfa2 * (1.0d0 - csfc2)
     csfsm3 = 9.0d0 * csfa1 * csfa2
     ceps = cmu0**(3.0d0 + estrp/estrn)
     eepst = 1.5d0 + estrm/estrn
     eepsp = (-1.0d0) / estrn
     ctls = cmu0**3.0d0
     cpsi30 = (cpsi3p + cpsi3n) * 0.5d0
     cpsi3s = (cpsi3p - cpsi3n) * 0.5d0
     cpsife = estrn * (cmu0**estrp)
     cpsibs = cmu0**(estrp-2.0d0*estrm) 
     cpsibb = cmu0**(estrp-2.0d0*estrm) * (ckarm*z0btm)**(estrn)
     cpsilm = sqrt(0.56d0) * (cmu0**(estrp/estrn))
     epsilm = estrm/estrn + 0.5d0
     cscnr = sqrt(1.5d0*scntke)*cmu0/ckarm
     scnpss = ckarm*ckarm/(cpsi2*cmu0*cmu0) * &
       &    ( estrn*estrn - 4.0d0/3.0d0*cscnr*estrn*estrm &
       &    - 1.0d0/3.0d0*cscnr*estrn &
       &    + 2.0d0/9.0d0*estrm*cscnr*cscnr &
       &    + 4.0d0/9.0d0*cscnr*cscnr*estrm*estrm )
     csftkw = sqrt(1.5d0*scntke)*cmu0
     cpsifw = estrm * scntke * (cmu0**estrp)
     if (obtkei) then
        ctaubt = btmfrc
     else
        ctaubt = 0.0d0
     end if

     nitr = nitr0
     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend
           csamv(ij, k) = (  amftz(ij   , k) + amftz(ij+le , k) &
             &             + amftz(ij+ln, k) + amftz(ij+lne, k) &
             &            ) * amfvz(ij, k)
           if (csamv(ij, k) .eq. 0.0d0) then
              csamv(ij, k) = 0.0d0
           else
              csamv(ij, k) = 1.0d0 / csamv(ij, k)
           end if
        end do
        do ij = ijstr, ijend
           csahv(ij, k) = (  amfvz(ij   , k) + amfvz(ij+lw , k) &
             &             + amfvz(ij+ls, k) + amfvz(ij+lsw, k) &
             &            ) * amftz(ij, k)
           if (csahv(ij, k) .eq. 0.0d0) then
              csahv(ij, k) = 0.0d0
           else
              csahv(ij, k) = 1.0d0 / csahv(ij, k)
           end if
        end do
     end do
     !$acc end kernels
     
! --- reduce the background diffusion around EQ. ---
!
     pi = atan( 1.d0 )*4.d0
     omega = 2.d0 * pi / 86400.d0
     cor30 =   2.d0 * omega * sin( pi*30.d0/180.d0 )
     bfq   = 5.24d-3        ! /sec

     do k = kstr, kend
        do ij = ijstr-nxdim-1, ijend
           cort=(cor(ij)+cor(ij+lw)+cor(ij+lsw)+cor(ij+ls))*0.25d0
           cort=abs(cort)
           if(cort.gt.cor30) then
              ahv03d(ij,k)=ahv0(k-kstr+1)
           else
!              ahv03d(ij,k)=ahv0(k-kstr+1) &
!                &         *cort*acosh(bfq/cort)/cor30/acosh(bfq/cor30)
              ahv03d(ij, k) = ahv0(k-kstr+1) * cort / cor30 * &
                &             log(  bfq / cort &
                &                 + sqrt(bfq**2 / cort**2 - 1.d0)) &
                &           / log(  bfq / cor30 &
                &                 + sqrt(bfq**2 / cor30**2 - 1.d0))
!              ahv03d(ij,k)= max(ahv03d(ij, k), 1.d-2)
              ahv03d(ij,k)= max(ahv03d(ij, k), 1.d-4)
              ahv03d(ij,k)= ahv03d(ij,  k) * rahv(k-kstr+1) &
                &         + ahv0(k-kstr+1) * (1.0 - rahv(k-kstr+1)) 
           endif
        end do
     end do
     
!----- reducing background vert. diffusivity in Arctic Ocean
     if ( ovdfao ) then
        corao = 2.D0 * omega * sin( pi * 65.D0 / 180.D0 )
        do k = kstr, kstr+mzao-1
           do ij = ijstr-nxdim-1, ijend
              cort = (  cor(ij)     + cor(ij+lw) &
                   &  + cor(ij+lsw) + cor(ij+ls) ) * 0.25d0
              if ( cort > corao ) then
                 ahv03d(ij,k) = ahv0ao(k-kstr+1)
              end if
           end do
        end do
     end if
     !$acc update device(ahv03d)
     
!---- applying turbulent energy dissipation rate
     !$acc kernels default(present)
     do ij = 1, nxydim
        tedn2d(ij) = 0.d0
        tedf2d(ij) = 0.d0
     end do
     !$acc end kernels
     if ( iamn == 0 ) then
        write(jfpar, *) ' Turbulent energy dissipation rate (near-field) is not used.'
     else
!---- reading file of near-field tidal energy dissipation rate
#ifdef OPT_IO_COCOMPI
        call mpi_filopn(mpi_fh, cftedn, 'READ')
        disp=0
        call mpi_read_chead(chead, mpi_fh, disp, icread)
        call mpi_read_2d(tedn2d, mpi_fh  , disp)
        call mpi_filcls(mpi_fh)
#else
        allocate ( buf2(1:nxg,1:nyg) )
        allocate ( g2d(1:nxgdim,1:nygdim) )
        if ( myrank == iroot ) then
           call filopn( nfted, cftedn, 'READ' )
           rewind( nfted )
           read( nfted ) chead
           read( nfted ) buf2
           call filcls( nfted )
           do j = 1, nyg
              do i = 1, nxg
                 g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
              end do
           end do
        end if
        call scatter_2d( tedn2d, g2d )
        deallocate( buf2, g2d )
        !$acc update device(tedn2d)
#endif
#ifdef OPT_TRIPOLE
        call shift1(tedn2d,                                      &
    &                nxdim,  nydim,      1,                      &
    &                 1.d0,      0,      0 )
#else
        call shift1(tedn2d,                                      &
    &                nxdim,  nydim,      1)
#endif
        rzeta = 1.d0 / zeta * 1.d-2
        !$acc kernels default(present)
        do ij = 1, nxydim
           depth0(ij, kstr) = 0.d0
        end do
        do k = kstr+1, kend+1
           do ij = 1, nxydim
              depth0(ij, k) = depth0(ij, k-1) + dz(ij, k-1)
           end do
        end do
        !$acc end kernels
     end if
     if ( iamf == 0 ) then
        write(jfpar, *) ' Turbulent energy dissipation rate (far-field) is not used.'
        !$acc kernels default(present)
        do ij = 1, nxydim
           tedf2d(ij) = 0.d0
        end do
        !$acc end kernels
     else
!---- reading file of far-field tidal energy dissipation rate
#ifdef OPT_IO_COCOMPI
        call mpi_filopn(mpi_fh, cftedf, 'READ')
        disp=0
        call mpi_read_chead(chead, mpi_fh, disp, icread)
        call mpi_read_2d(tedf2d, mpi_fh  , disp)
        call mpi_filcls(mpi_fh)
#else
        allocate ( buf2(1:nxg,1:nyg) )
        allocate ( g2d(1:nxgdim,1:nygdim) )
        if ( myrank == iroot ) then
           call filopn( nfted, cftedf, 'READ' )
           rewind( nfted )
           read( nfted ) chead
           read( nfted ) buf2
           call filcls( nfted )
           do j = 1, nyg
              do i = 1, nxg
                 g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
              end do
           end do
        end if
        call scatter_2d( tedf2d, g2d )
        deallocate( buf2, g2d )
        !$acc update device(tedf2d)
#endif
        if (ofvcnt) then
           write(jfpar, *) 'Vertically constant for far-field mixing'
        else
           if (ofvpn) then
              write(jfpar, *) 'Dissipation rate is prop. to N for far-field mixing'
           else
              write(jfpar, *) 'Dissipation rate is prop. to N^2 for far-field mixing'
           end if
        end if
#ifdef OPT_TRIPOLE
        call shift1(tedf2d,                                      &
    &                nxdim,  nydim,      1,                      &
    &                 1.d0,      0,      0 )
#else
        call shift1(tedf2d,                                      &
    &                nxdim,  nydim,      1)
#endif
     end if
!---
  end if

!     -- second step of Euler-Eackward sheme --
  if ( ( itst .eq. 1 ).and.( ieuler .eq. 2 ) ) then
     return
  endif

#ifdef OPT_BBL
  call rmmskt
  call admkt1
#endif

  !$acc kernels default(present)
  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        dzsig (ij, k) = (hy(ij) + zbot) * ds(k)
        rzmsig(ij, k) = 1.d0 / (hy(ij) + zbot) / dsm(k)
        dzmsig(ij, k) = (hy(ij) + zbot) * dsm(k)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        dzsig (ij, k) = dz(ij, k)
        rzmsig(ij, k) = 1.d0 / dzm(ij, k)
        dzmsig(ij, k) = dzm(ij, k)
     end do
  end do
  do k = 1, nzdim
     do ij = 1, nxydim
        depth (ij, k) = 0.d0
        depthm(ij, k) = 0.0d0
     end do
  end do
  do k = kstr+1, kend+1
     do ij = 1, nxydim
        depth(ij, k) = depth(ij, k-1) + dzsig(ij, k-1)
     end do
  end do
  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        depthm(ij, k) = depthm(ij, k-1) + (hy(ij) + zbot) * dsm(k)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        depthm(ij, k) = depthm(ij, k-1) + dzm(ij, k)
     end do
  end do

  do k = 1, nzdim
     do ij = 1, nxydim
        drdz(ij, k) = 0.d0
     end do
  end do

  do k = kstr+1, kend
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1   = cc0(k) &
          & + (cc1(k) + (cc2(k) + cc3(k) * tl) * tl) * tl &
          & + (cc4(k) + cc5(k) * tl + cc6(k) * sl) * sl
        p2   = d0(k) &
          & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2

        tl = ty(ij, k-1, 1)
        sl = ty(ij, k-1, 2)
        p1   = cc0(k) &
          & + (cc1(k) + (cc2(k) + cc3(k) * tl) * tl) * tl &
          & + (cc4(k) + cc5(k) * tl + cc6(k) * sl) * sl
        p2   = d0(k) &
          & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlu = p1 / p2

        drdz(ij, k) = gravit / r0 * (rl - rlu) * rzmsig(ij, k)
!        drdz(ij, k) = max(drdz(ij, k), 0.d0)
        dudz        = (  uy(ij   , k-1) + uy(ij+lw , k-1) &
          &            + uy(ij+ls, k-1) + uy(ij+lsw, k-1) &
          &            - uy(ij   , k  ) - uy(ij+lw , k  ) &
          &            - uy(ij+ls, k  ) - uy(ij+lsw, k  )) * &
          &           0.25d0 * rzmsig(ij, k)
        dvdz        = (  vy(ij   , k-1) + vy(ij+lw , k-1) &
          &            + vy(ij+ls, k-1) + vy(ij+lsw, k-1) &
          &            - vy(ij   , k  ) - vy(ij+lw , k  ) &
          &            - vy(ij+ls, k  ) - vy(ij+lsw, k  )) * &
          &           0.25d0 * rzmsig(ij, k)
        duvdz(ij, k) = dudz * dudz + dvdz * dvdz
     end do
  end do
  !$acc end kernels

!     Initial guess for TKE and PSI
!     NITR, N. of iteration, at the first step is set to 20.
  if (ofirst .and. oeof) then
!    (Assume minimum TKE, neutral stratification(i.e., GH=0),
!     and shear production = disspiation.  Not used now.)
!     do k = kstr+1, kend
!        do ij = ijstr-nxdim-1, ijend+nxdim+1
!           tke(ij, k) = tkemin
!           tls(ij, k) = sqrt( (cmu0**3.0d0) * tke(ij, k) &
!             &            / ( sqrt(2.0d0)*csfsm1*csfe*duvdz(ij, k) ) )
!           tls(ij, k) = min(tls(ij, k), &
!             &    sqrt(0.56d0 * tke(ij, k) / max(drdz(ij, k), epscmp)) )
!           psi(ij, k) = &
!             &    (cmu0**estrp)*(tke(ij, k)**estrm)*(tls(ij, k)**estrn)
!        end do
!     end do
!
!     do k = kstr, kend+1
!        do ij = ijstr-nxdim-1, ijend+nxdim+1
!           tke(ij, k) = max(tkemin, tke(ij, k)*amftz(ij, k))
!           psi(ij, k) = max(psimin, psi(ij, k)*amftz(ij, k))
!        end do
!     end do
!
!     do ij = 1, nxydim
!        taubtm(ij) = 0.0d0
!     end do
!!     do k = kstr, kend
!!        do ij = ijstr-nxdim-1-nxdim-1, ijend+nxdim+1
!!           taubtm(ij) = taubtm(ij) + btmfrc * rhoo * &
!!             &      (uy(ij, k) * uy(ij, k) + vy(ij, k) * vy(ij, k)) &
!!             &      * amskb(ij, k)
!!        end do
!!     end do
!     do ij = ijstr-nxdim-1, ijend+nxdim+1
!        avrtx = (  taux(ij)    + taux(ij+lw) &
!          &      + taux(ij+ls) + taux(ij+lsw)) * 0.25d0
!        avrty = (  tauy(ij)    + tauy(ij+lw) &
!          &      + tauy(ij+ls) + tauy(ij+lsw)) * 0.25d0
!        ufrc2s(ij) = sqrt(avrtx * avrtx + avrty * avrty) / rhoo &
!          &      * amskt(ij, kstr)
!        tke(ij, kstr) = max( ufrc2s(ij) / (cmu0*cmu0), tkemin)
!!        ufrc2b(ij) = &
!!          &      ( taubtm(ij) + taubtm(ij+lw) &
!!          &      + taubtm(ij+ls) + taubtm(ij+lsw)) * 0.25d0 &
!!          &      / rhoo
!!        tke(ij, nbot(ij)+1) = max( ufrc2b(ij) / (cmu0*cmu0), tkemin )
!     end do
!
     nitr = nitr00
  end if

  do iitr = 1, nitr
!    tke backup
     !$acc kernels default(present)
     do k = 1, nzdim
        do ij = 1, nxydim
           tke0(ij, k) = tke(ij, k)
        end do
     end do
     !$acc end kernels
     
!     dissipation epsil and turbulent length scale tls.
!     an upper limit for tls is also introduced (eq.(42))
!     turbultent richardson number rit is also calculated.

     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           epsil(ij, k) = ceps &
             &          * (tke(ij, k)**eepst) * (psi(ij, k)**eepsp)
           tls(ij, k) = ctls * (tke(ij, k)**1.5d0) / epsil(ij, k)
           if (drdz(ij, k) .gt. 0.0d0) then
              tls(ij, k) = min(tls(ij, k), &
                &              sqrt( 0.56d0 * tke(ij, k) / &
                &                    max(drdz(ij, k), epscmp)) )
           end if
           rit(ij, k) = drdz(ij, k) * tls(ij, k) * tls(ij, k) &
             &          / tke(ij,k) * 0.5d0
        end do
     end do
     !$acc end kernels
     
!    Kantha and Clayson (1994) quasi-equilibrium stability function
!     Eq.(33) of Warner et al. (2005, OM) has TYPOGRAPHICAL ERRORs;
!      - It should be Gh = Ghul - ( (Ghul-Ghcr)**2 / (Ghul+Gh0-2Ghcr) )
!      - SFGHCR should be 0.01 instead of 0.02.
!     See eq. (19) of Burchard et al. (1999, JMS) for the correct formula
!     (note that the definitions of the buoyancy parameter are different
!      between the two papers).
!
!    '12.01.19: bug fix (due to the typograpical error dscribed above)

     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ghul(ij, k) = (-0.5d0) * &
             &           drdz(ij, k) * tls(ij, k) * tls(ij, k) / tke(ij, k)
           if (ghul(ij, k) .gt. sfghcr) then
              gh(ij, k) = ghul(ij, k) - ((ghul(ij, k) - sfghcr)**2.0d0) &
                &              / (ghul(ij, k) + sfgh0 - 2.0d0 * sfghcr)
           else
              gh(ij, k) = max(sfghmn, ghul(ij, k))
           endif
           sh(ij, k) = csfsh1 / (1.0d0 - csfsh2 * gh(ij, k))
           sm(ij, k) = (csfsm1 + csfsm2 * sh(ij, k) * gh(ij, k)) &
             &         / (1.0d0 - csfsm3 * gh(ij, k))
        end do
     end do
     !$acc end kernels

!     Vertical eddy viscosity at T-grid amvt and diffusivity ahv
     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           q = sqrt( 2.0d0 * tke(ij, k) )
           amvt(ij, k) = max(min( amvmax, &
             &                csfe * q * tls(ij, k) * sm(ij, k) ),amvmin)
           ahv(ij, k) = max(min( ahvmax, &
             &               csfe * q * tls(ij, k) * sh(ij, k) ),ahvmin)
        end do
     end do
     !$acc end kernels
     
!     Wall function fwall is just a dummy
     !$acc kernels default(present)
     do k = kstr-1, kend+1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           fwall(ij, k) = 1.0d0
        end do
     end do
     !$acc end kernels
     
!     Frictional velocity and surface roughness
     !$acc kernels default(present)
     do ij = ijstr-nxdim-1-nxdim-1, ijend+nxdim+1
        taubtm(ij) = 0.0d0
     end do
     do k = kstr, kend
        do ij = ijstr-nxdim-1-nxdim-1, ijend+nxdim+1
           taubtm(ij) = taubtm(ij) + ctaubt * rhoo * &
             &        ( uy(ij, k) * uy(ij, k) + vy(ij, k) * vy(ij, k) ) &
             &        * amskb(ij, k)
        end do
     end do
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        avrtx = (  taux(ij)    + taux(ij+lw) &
          &      + taux(ij+ls) + taux(ij+lsw)) * 0.25d0
        avrty = (  tauy(ij)    + tauy(ij+lw) &
          &      + tauy(ij+ls) + tauy(ij+lsw)) * 0.25d0
        ufrc2s(ij) = sqrt(avrtx * avrtx + avrty * avrty) / rhoo &
          &          * amskt(ij, kstr)
        if (oswnoi) then
           avrtox = ( tauaox(ij)    + tauaox(ij+lw) &
             &      + tauaox(ij+ls) + tauaox(ij+lsw)) * 0.25d0
           avrtoy = ( tauaoy(ij)    + tauaoy(ij+lw) &
             &      + tauaoy(ij+ls) + tauaoy(ij+lsw)) * 0.25d0
           ufrc2o(ij) = sqrt(avrtox * avrtox + avrtoy * avrtoy) / rhoo &
             &        * amskt(ij, kstr)
           ufrc2i(ij) = sqrt( (avrtx - avrtox)**2.0d0 &
             &              + (avrty - avrtoy)**2.0d0 ) / rhoo &
             &        * amskt(ij, kstr)
        else
           ufrc2o(ij) = ufrc2s(ij)
           ufrc2i(ij) = 0.0d0
        end if
        ufrc3s(ij) = ufrc2s(ij) * sqrt(ufrc2s(ij))
        ufrc3o(ij) = ufrc2o(ij) * sqrt(ufrc2o(ij))
        ufrc2b(ij) = &
          &        ( taubtm(ij) + taubtm(ij+lw) &
          &        + taubtm(ij+ls) + taubtm(ij+lsw)) * 0.25d0 &
          &        / rhoo
        if (osfcwv) then
           z0sf2d(ij) = max( &
             &          alphch/gravit*ufrc2o(ij)+alphci/gravit*ufrc2i(ij), &
             &          z0sfmn )
        else
           z0sf2d(ij) = z0sfc
        endif
     end do
     !$acc end kernels

!    Schmidt number for psi
     !$acc kernels default(present)
     do k = 1, nzdim
        do ij = 1, nxydim
           scnp3d(ij, k) = scnpsi
        end do
     end do
     !$acc end kernels
     if (osfcwv) then
        !$acc kernels default(present)
        do ij = 1, nxydim
           prepmx(ij) = epscmp
        end do
         
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           preps0 = 0.0d0
           do k = kstr+1, nbot(ij)
              if (prepmx(ij) .lt. 1.0d0) then
                 preps = amvt(ij, k) * duvdz(ij, k) &
                   &   / max(epsil(ij, k), epscmp)
                 prepmx(ij) = max(prepmx(ij), preps) 
                 if (preps .ge. 1.0d0) then
                    dpsi0(ij) = depth(ij, k-1) + dzsig(ij, k-1) &
                      &         * (1.0d0-preps0)/(preps-preps0)
                 end if
                 preps0 = preps
              end if
           end do
           if (prepmx(ij) .lt. 1.0d0) then
              dpsi0(ij) = depth(ij, nbot(ij)+1) &
                &         / max(prepmx(ij), epscmp)
           end if
        end do

        do k = kstr, kend
           do ij = ijstr-nxdim-1, ijend+nxdim+1
              rscnp = min(depthm(ij, k-1)/max(dpsi0(ij),epscmp), 1.0d0)
              scnp3d(ij, k) = &
                &  (1.0d0-rscnp)*scnpss*2.0d0/(fwall(ij,k-1)+fwall(ij,k)) &
                &  + rscnp*scnpsi
           end do
        end do
        !$acc end kernels
     else
        !$acc kernels default(present)
        do k = kstr, kend+1
           do ij = ijstr-nxdim-1, ijend+nxdim+1
              scnp3d(ij, k) = scnpsi
           end do
        end do
        !$acc end kernels
     end if


!    tke equation
     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           cdkdtl(ij, k) = &  ! gamma = -(p+b)/k (temporary)
             &           (  ahv(ij, k) *  drdz(ij, k) &
             &           - amvt(ij, k) * duvdz(ij, k) ) &
             &           / tke(ij, k) * amftz(ij, k)
           cdmp(ij, k) = &    ! delta = epsil/k
             &   epsil(ij, k) / tke(ij, k) * amftz(ij, k)
           adefwd(ij, k) = &
             &  - ( ahv(ij, k) * drdz(ij, k) - amvt(ij, k) * duvdz(ij, k) &
             &    + epsil(ij, k)) * amftz(ij, k)
        end do
     end do
     !$acc end kernels
     
!    -- Semi-implicit : adjusted alps in each time step --
     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           p = cdkdtl(ij, k) * dt   ! gamma * dt
!           -- foward after initial-growth phase --
           if ( ( ( rit(ij, k).lt.ritc) &
             &    .and.( cdkdtl(ij,k).lt.0.0d0) ) &
             &  .or.( p .lt. -300.0d0 ) ) then
              alps  =  1.0d0
!!           -- implicit for strong damping phase --
!           elseif ( p .ge. 1.0d0 ) then
!              alps  =  0.0d0
           else
!           -- semi-implicit for initial-growth and damping phase --
              p = sign(1.0d0, p) * max(1.0d-6, abs(p))
              alps  =  ( -1.0d0 + (1.0d0+p) * exp(-p) ) &
                &      / ( p * ( exp(-p) - 1.0d0 ) )
           endif
           cdkdtl(ij, k) = (1.0d0 - alps) * cdkdtl(ij, k) &
                &        + cdmp(ij, k)    ! (1-mu)*gamma + delta
        end do
     end do
     !$acc end kernels
     
!    --- Boundary condition for tke ---
!     "diffz(ij, kstr+1)=0"  & "diffz(ij, nbot(ij)+1)=0"
!      ==>  "aa(ij, kstr+1)=0" & "ac(ij, nbot(ij))=0".
     !$acc kernels default(present)
     do k = 1, nzdim
        do ij = 1, nxydim
           diffz(ij, k) = 0.d0
           fez(ij, k) = 0.d0
        end do
     end do
     do k = kstr+2, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           diffz(ij, k) = ( amvt(ij, k-1) + amvt(ij, k) ) * &
             &            0.5d0 / scntke / dzsig(ij, k-1) * &
             &            amskt(ij, k)
           fez(ij, k) = diffz(ij, k) * (tke(ij, k-1) - tke(ij, k))
        end do
     end do
     !$acc end kernels
     
!    tke flux by surface wave breaking, after Carniel et al.(2009)
     !$acc kernels default(present)
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        fez(ij, kstr+1) = cw * ufrc3o(ij)
     end do
      
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ade(ij, k) = (fez(ij, k) - fez(ij, k+1)) * rzmsig(ij, k) &
             &        + adefwd(ij, k)
        end do
     end do

     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           aa(ij, k) = - dt * diffz(ij, k) * rzmsig(ij, k)
           ac(ij, k) = - dt * diffz(ij, k+1) * rzmsig(ij, k)
           ab(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k) &
             &       + dt * cdkdtl(ij, k)
        end do
     end do
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        ac(ij, kstr+1) = ac(ij, kstr+1) / ab(ij, kstr+1)
        ade(ij, kstr+1) = ade(ij, kstr+1) / ab(ij, kstr+1)
     end do
     do k = kstr+2, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           fc = 1.d0 / (ab(ij, k) - aa(ij, k) * ac(ij, k-1))
           ac(ij, k) = ac(ij, k) * fc
           ade(ij, k) = (ade(ij, k) - aa(ij, k) * ade(ij, k-1)) * fc
        end do
     end do
     do k = kend-1, kstr+1, -1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ade(ij, k) = ade(ij, k) - ac(ij, k) * ade(ij, k+1)
        end do
     end do

     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           tke(ij, k) = (  atfilt * tke(ij, k) &
             &           + (1.0d0 - atfilt) * &
             &             (tke(ij, k) + dt * ade(ij, k)) &
             &          ) * amftz(ij, k)
        end do
     end do
     !$acc end kernels
     
!    Surface and bottom tke estimation
!    for surface tke, surface wave breaking effect is accounted
!     (after Carniel et al. (2009)).

     !$acc kernels default(present)
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        tke(ij, kstr) = 1.0d0 / (cmu0*cmu0) * &
          &        (ufrc3s(ij) + ufrc3o(ij) * csftkw * cw)**(2.0d0/3.0d0)
        tke(ij, nbot(ij)+1) = ufrc2b(ij) / (cmu0*cmu0)
     end do

     do k = kstr, kend+1
        do ij = 1, nxydim
           if (tke(ij, k) .lt. tkemin) then
              ctkemn(ij, k) = 1.0d0
           else
              ctkemn(ij, k) = 0.0d0
           end if
           tke(ij, k) = max(tke(ij, k), tkemin)
        end do
     end do
     !$acc end kernels
     
!    GLS quantity psi equation
     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           cdkdtl(ij, k) = &  ! gamma = -(c1*p+c3*b)/k (temporary) 
             &  ( ( cpsi30 - cpsi3s * sign(1.0d0, drdz(ij, k)) ) &
             &      * ahv(ij, k) * drdz(ij, k) &
             &    - cpsi1 * amvt(ij, k) * duvdz(ij, k) ) &
             &      / tke0(ij, k) * amftz(ij, k)
           cdmp(ij, k) = &   ! delta = c2*epsil*fwall/k
             &  cpsi2 * fwall(ij, k) * epsil(ij, k) / tke0(ij, k) &
             &  * amftz(ij, k)
           adefwd(ij, k) = &
             &  - ( cdkdtl(ij, k) + cdmp(ij, k) ) * psi(ij, k) &
             &  * amftz(ij, k)
        end do
     end do
     !$acc end kernels
     
!    -- Semi-implicit : adjusted alps in each time step --

     !$acc kernels default(present)
     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           p = cdkdtl(ij, k) * dt   ! gamma * dt
!          -- foward after initial-growth phase --
           if ( ( ( rit(ij, k).lt.ritc) &
             &    .and.( cdkdtl(ij, k).lt.0.0d0) ) &
             &  .or.( p .lt. -300.0d0 ) ) then
              alps  =  1.0d0
!!           -- implicit for strong damping phase --
!            elseif ( p .ge. 1.0d0 ) then
!               alps  =  0.0d0
           else
!          -- semi-implicit for initial-growth and damping phase --
              p = sign(1.0d0, p) * max(1.0d-6, abs(p))
              alps  =  ( -1.0d0 + (1.0d0+p) * exp(-p) ) &
                &      / ( p * ( exp(-p) - 1.0d0 ) )
           endif
           cdkdtl(ij, k) = (1.0d0 - alps) * cdkdtl(ij, k) &
             &           + cdmp(ij, k)  ! (1-mu)*gamma + delta
        end do
     end do
     !$acc end kernels
     
!    --- Boundary condition for psi ---
!     "diffz(ij, kstr+1)=0"  & "diffz(ij, nbot(ij)+1)=0"
!       ==> "aa(ij, kstr+1)=0" & "ac(ij, nbot(ij))=0"

     !$acc kernels default(present)
     do k = 1, nzdim
        do ij = 1, nxydim
           diffz(ij, k) = 0.d0
           fez(ij, k) = 0.d0
        end do
     end do
     do k = kstr+2, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           diffz(ij, k) = ( amvt(ij, k-1) + amvt(ij, k) ) * &
             &            0.5d0 / scnp3d(ij, k-1) / dzsig(ij, k-1) * &
             &            amskt(ij, k)
           fez(ij, k) = diffz(ij, k) * (psi(ij, k-1) - psi(ij, k))
        end do
     end do
     !$acc end kernels
     
!    fez(ij, kstr+1) & fez(ij, nbot(ij)+1)
!    eq.(54) of Warner et al. (2005, om) has TYPOGRAPHICAL ERROR;
!      k^{n} must be \kappa^{n}.

     !$acc kernels default(present)
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        k = kstr+1
        tketmp = (tke0(ij, k-1) + tke0(ij, k))*0.5d0
        dstwal = z0sf2d(ij) + dzsig(ij, k-1)*0.5d0
        fez(ij, k) = &
          &   - cpsife / scnp3d(ij, k-1) &
          &     * amvt(ij, k) &
          &     * ( tketmp**estrm ) &
          &     * ( ckarm**estrn ) &
          &     * ( dstwal**(estrn-1.0d0) ) &
          &   - cpsifw / scnp3d(ij, k-1) &
          &     * ( tketmp**(estrm-1.0d0) ) &
          &     * ( ckarm**estrn ) &
          &     * ( dstwal**estrn ) &
          &     * cw * ufrc3o(ij)
        if (nbot(ij) > kstr) then
           k = nbot(ij)+1
           tketmp = (tke0(ij, k-1) + tke0(ij, k))*0.5d0
           dstwal = z0btm + dzsig(ij, k-1)*0.5d0
           fez(ij, k) = &
             &   cpsife / scnp3d(ij, k-1) &
             &     * amvt(ij, k-1) &
             &     * ( tketmp**estrm ) &
             &     * ( ckarm**estrn ) &
             &     * ( dstwal**(estrn-1.0d0) )
        end if
     end do

     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ade(ij, k) = ( fez(ij, k) - fez(ij, k+1) ) * rzmsig(ij, k) &
             &          + adefwd(ij, k)
        end do
     end do

     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           aa(ij, k) = - dt * diffz(ij, k) * rzmsig(ij, k)
           ac(ij, k) = - dt * diffz(ij, k+1) * rzmsig(ij, k)
           ab(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k) &
             &       + dt * cdkdtl(ij, k)
        end do
     end do
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        ac(ij, kstr+1) = ac(ij, kstr+1) / ab(ij, kstr+1)
        ade(ij, kstr+1) = ade(ij, kstr+1) / ab(ij, kstr+1)
     end do
     do k = kstr+2, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           fc = 1.d0 / (ab(ij, k) - aa(ij, k) * ac(ij, k-1))
           ac(ij, k) = ac(ij, k) * fc
           ade(ij, k) = (ade(ij, k) - aa(ij, k) * ade(ij, k-1)) * fc
        end do
     end do
     do k = kend-1, kstr+1, -1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ade(ij, k) = ade(ij, k) - ac(ij, k) * ade(ij, k+1)
        end do
     end do

     do k = kstr+1, kend
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           psi(ij, k) = (  atfilt * psi(ij, k) &
             &           + (1.0d0 - atfilt) * &
             &             (psi(ij, k) + dt * ade(ij, k) ) &
             &          ) * amftz(ij, k)
        end do
     end do
     !$acc end kernels
     
!    Surface and bottom psi estimation
!    (this part is only for output and will not be referred below,
!     thus can be removed).

     !$acc kernels default(present)
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        psi(ij, kstr) = cpsibs * (ckarm*z0sf2d(ij))**(estrn) &
          &                    * (ufrc2s(ij)**(2.0d0*estrm))
        psi(ij, nbot(ij)+1) = cpsibb * (ufrc2s(ij)**(2.0d0*estrm))
     end do
     !$acc end kernels
     
!    Limit on psi (psimin, and eq.(43) if stable stratification)

     !$acc kernels default(present)
     do k = kstr, kend+1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           if (drdz(ij, k) .gt. 0.0d0) then
              psilm = ( cpsilm * (tke(ij, k)**epsilm) &
                &   / sqrt(max(drdz(ij, k), epscmp)) ) ** estrn
              if (estrn .gt. 0.0d0) then
                 psi(ij, k) = max(min(psi(ij, k), psilm), psimin)
              else
                 psi(ij, k) = max(max(psi(ij, k), psilm), psimin)
              endif
           else
              psi(ij, k) = max(psi(ij, k), psimin)
           end if
        end do
     end do
     !$acc end kernels
  end do
  
  if (ofirst) then
     nitr = nitr0
     ofirst = .false.
  end if

  !$acc kernels default(present)
  do k = kstr+1, kend
     do ij = ijstr, ijend
        amv(ij, k) = ( amvt(ij    , k) * amftz(ij    , k) &
          &          + amvt(ij+le , k) * amftz(ij+le , k) &
          &          + amvt(ij+ln , k) * amftz(ij+ln , k) &
          &          + amvt(ij+lne, k) * amftz(ij+lne, k) &
          &          ) * csamv(ij, k)
     end do
  end do
  !$acc end kernels

  
! -- Smoothing --
  !$acc kernels default(present)
  do k = kstr+1, kend
     do ij = ijstr-nxdim-1, ijend
        amvt(ij, k) = ( ahv(ij    , k) * amftz(ij    , k) &
          &           + ahv(ij+le , k) * amftz(ij+le , k) &
          &           + ahv(ij+ln , k) * amftz(ij+ln , k) &
          &           + ahv(ij+lne, k) * amftz(ij+lne, k) &
          &           ) * csamv(ij, k)
     end do
     do ij = ijstr, ijend
        ahv(ij, k) = ( amvt(ij    , k) * amfvz(ij    , k) &
          &          + amvt(ij+lw , k) * amfvz(ij+lw , k) &
          &          + amvt(ij+ls , k) * amfvz(ij+ls , k) &
          &          + amvt(ij+lsw, k) * amfvz(ij+lsw, k) &
          &          ) * csahv(ij, k)
     end do
  end do
  do k = kstr, kstr+mz-1
     do ij = ijstr-nxdim-1, ijend
        amv(ij, k) = min(amvmax, max(amv(ij, k), amv0(k-kstr+1)))
        ahv(ij, k) = max(ahv(ij, k), ahv03d(ij, k))
     end do
  end do
  do k = kstr+mz, kend
     do ij = ijstr, ijend
        amv(ij, k) = amv0(k-kstr+1)
        ahv(ij, k) = ahv03d(ij, k)
     end do
  end do
  !$acc end kernels


#ifdef OPT_TRIPOLE
  call shift2(   tke,    psi, &
    &          nxdim,  nydim,  nzdim, &
    &           1.d0,      0,      0 )
#endif

!!--- sea-surface elevation is not considered
!!---   for vertical structure function of energy dissipation rate
!  if ( iamn /= 0 .or. iamf /= 0 ) then
!     dzmsig = dzm
!     depth = depth0
!  end if
!!---
  
!--- tidal turbulent energy dissipation rate
! near-field
  if ( iamn /= 0 ) then
     !$acc kernels default(present)
     do ij = ijstr, ijend
        gint(ij) = 0.d0
     end do
     do k = kstr+1, kend
        do ij = 1, nxydim
           dep = depth(ij, nbot(ij) + 1) ! depth of bottom
           gint(ij) = gint(ij) + &
                & dzmsig(ij, k) * exp((depth(ij, k) - dep) * rzeta) * amftz(ij, k)
        end do
     end do
     where(gint /= 0.d0) gint = 1.d0 / gint
     do k = kstr+1, kend
        do ij = 1, nxydim
           dep = depth(ij, nbot(ij) + 1) ! depth of bottom
           tedn3d(ij, k) = gint(ij) * tedn2d(ij) * exp((depth(ij, k) - dep) * rzeta) * amftz(ij, k)
        end do
     end do
     !$acc end kernels
  end if
! far-field
  if ( iamf /= 0 ) then
     !$acc kernels default(present)
     do ij = ijstr, ijend
        gint(ij) = 0.d0
     end do
     !$acc end kernels
     if (ofvcnt) then
        !$acc kernels default(present)
        do k = kstr+1, kend
           do ij = 1, nxydim
              gint(ij) = gint(ij) + &
                   & dzmsig(ij, k) * amftz(ij, k)
           end do
        end do
        do ij = ijstr, ijend
           if (gint(ij) /= 0.d0) then
              gint(ij) = 1.d0 / gint(ij)
           end if
        end do
        do k = kstr+1, kend
           do ij = 1, nxydim
              tedf3d(ij, k) = gint(ij) * tedf2d(ij) * amftz(ij, k)
           end do
        end do
        !$acc end kernels
     else
        if (ofvpn) then ! prop to N
           !$acc kernels default(present)
           do k = kstr+1, kend
              do ij = 1, nxydim
                 gint(ij) = gint(ij) + &
                      & dzmsig(ij, k) * sqrt(abs(drdz(ij, k))) * amftz(ij, k)
              end do
           end do
           do ij = ijstr, ijend
              if (gint(ij) /= 0.d0) then
                 gint(ij) = 1.d0 / gint(ij)
              end if
           end do
           do k = kstr+1, kend
              do ij = 1, nxydim
                 tedf3d(ij, k) = gint(ij) * tedf2d(ij) * sqrt(abs(drdz(ij, k))) * amftz(ij, k)
              end do
           end do
           !$acc end kernels
        else ! prop to N2
           !$acc kernels default(present)
           do k = kstr+1, kend
              do ij = 1, nxydim
                 gint(ij) = gint(ij) + &
                      & dzmsig(ij, k) * drdz(ij, k) * amftz(ij, k)
              end do
           end do
           do ij = ijstr, ijend
              if (gint(ij) /= 0.d0) then
                 gint(ij) = 1.d0 / gint(ij)
              end if
           end do
           do k = kstr+1, kend
              do ij = 1, nxydim
                 tedf3d(ij, k) = gint(ij) * tedf2d(ij) * drdz(ij, k) * amftz(ij, k)
              end do
           end do
           !$acc end kernels
        end if
     end if
  end if

  !$acc kernels default(present)
  do k = kstr+1, kend
     do ij = ijstr, ijend
        ahvted(ij, k) = cgamma * (tedn3d(ij, k) + tedf3d(ij, k)) &
             &              / max(drdz(ij, k), epst) * amftz(ij, k)
        ahvted(ij, k) = min(ahvemx, ahvted(ij, k))
        ahv(ij, k) = max(ahv(ij, k), ahvted(ij, k))
        tedr(ij, k) = ahv(ij, k) * drdz(ij, k) / cgamma &
             &                  * amftz(ij, k)
     end do
  end do
  !$acc end kernels
  
  call chekin(ahvted, 'AHVTED', &
       &          'ahv by ted', 'cm^2/s', &
       &          nx,       ny,       nz,     nxyzdm, 'OCLVMT')
  call chekin(  tedr,   'TEDR', &
       &         'realized energy dissipation rate', 'cm^2/s^3', &
       &          nx,       ny,       nz,     nxyzdm, 'OCLVMT')
  call chekin(tedn3d,   'TEDN', &
       &              'near-field tidal energy dissipation rate', 'cm^2/s^3', &
       &          nx,       ny,       nz,     nxyzdm, 'OCLVMT')
  call chekin(tedf3d,   'TEDF', &
       &               'far-field tidal energy dissipation rate', 'cm^2/s^3', &
       &          nx,       ny,       nz,     nxyzdm, 'OCLVMT')


#ifdef OPT_BBL
  call rmmskt
  call admktb
#endif

!  call chekin(tke  ,'TKE', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')
!  call chekin(psi  ,'PSI', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(sm   ,'SM', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(epsil,'EPSIL', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(tls  ,'TLS', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(ghul ,'GHUL', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(gh   ,'GH', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(duvdz,'DUVDZ', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(drdz, 'DRDZ', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(ty(1,1,1), 'TODV', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(ctkemn, 'CTKEMN', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')

  return

end subroutine vdiff
! =====================================================================
subroutine puttao( &
  &                 taox,   taoy, &
  &                 caic,   cais)

  use bshft

  real(8), intent(in) ::   taox(nxydim),   taoy(nxydim)
  real(8), intent(in) ::   caic,   cais

  integer ::     ij
  logical, save :: ofirst = .true.
  if (ofirst) then
     ofirst = .false.
     !$acc enter data create(tauaox, tauaoy)
  end if

  !$acc kernels default(present)
  do ij = 1, nxydim
     tauaox(ij) = ( taox(ij) * caic &
       &          - taoy(ij) * sign(cais, cor(ij))) * &
       &          amskv(ij, kstr)
     tauaoy(ij) = ( taoy(ij) * caic &
       &          + taox(ij) * sign(cais, cor(ij))) * &
       &          amskv(ij, kstr)
  end do
  !$acc end kernels

#ifdef OPT_TRIPOLE
  call shift2( &
    &          tauaox, tauaoy, &
    &           nxdim,  nydim,      1, &
    &           -1.d0,     -1,     -1 )
#else
  call shift2( &
    &          tauaox,   tauaoy, &
    &           nxdim,    nydim,     1)
#endif

  return

end subroutine puttao
#ifdef OPT_BBL
! *********************************************************************
subroutine vdiffb( &
  &                   amv,    ahv)

! --- information -----------------------------------------------------
!
!  Vertical viscosity and diffusion coefficients for the bottom
! boundary layer.
!
!  HISTORY
!     '01.02.08  H.Hasumi
!
! ---------------------------------------------------------------------
  use ufile
  use qckot

  real(8) ::    amv(nxydim, nzdim),    ahv(nxydim, nzdim)

  real(8) ::   amvb(nxydim),   ahvb(nxydim)
  real(8) ::    ktb(nxydim),    kvb(nxydim)

  integer ::     ij
  integer ::     kt,     kv
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save :: amvbbl = 0.d0, ahvbbl = 0.d0
  namelist /nmbbdv/ amvbbl, ahvbbl

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmbbdv, iostat=istat)
     call cstnml(jfpar, 'vdiffb', 'nmbbdv', istat)
     write(jfpar, nmbbdv)
  end if

  do ij = ijstr, ijend
     kt = nbot(ij)
     kv = nbotv(ij)
     amv(ij, kv) = max( amvbbl, amv(ij, kv) ) * amskvb(ij) &
       &         + amv(ij, kv) * ( 1.d0 - amskvb(ij) )
     kvb(ij) = kv
     amvb(ij) = amv(ij, kv)
     ahv(ij, kt) = max( ahvbbl, ahv(ij, kt) ) * amsktb(ij) &
       &         + ahv(ij, kt) * ( 1.d0 - amsktb(ij) )
     ktb(ij) = kt
     ahvb(ij) = ahv(ij, kt)
  end do

  do ij = ijstr, ijend
     kt = nbot(ij)
     kv = nbotv(ij)
     amv(ij, kend) = amv(ij, kv)
     ahv(ij, kend) = ahv(ij, kt)
  end do

!  call chekin(ahv   ,'AHV3', &
!    &         nx,     ny,     nz, nxyzdm, 'OCN')      
!  call chekin(ahvb  ,'AHVB', &
!    &         nx,     ny,      1, nxydim, 'SFC')      
!  call chekin(amvb  ,'AMVB', &
!    &         nx,     ny,      1, nxydim, 'SFC')      
!  call chekin(kvb   ,'KVB', &
!    &         nx,     ny,      1, nxydim, 'SFC')      
!  call chekin(ktb   ,'KTB', &
!    &         nx,     ny,      1, nxydim, 'SFC')      

  return
end subroutine vdiffb
#endif
end module dvdif
