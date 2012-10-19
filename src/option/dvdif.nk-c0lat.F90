module dvdif

! --- information -----------------------------------------------------
!
!  Sea surface mixed layer model of Noh and Kim (1999, JGR).
!
!  HISTORY
!     '02.08.19  H.Hasumi: from COCO3.4
!     '03.10.14  H.Hasumi: bug fix (typographical error in the
!                          original paper)
!     '03.10.22  A.Oka: bug fix ( calculation of AMV, DML )
!     '03.11.14  T.Suzuki: reduce background mixing on equator
!     '03.11.27  A.Oka: adjusted semi-implicit time integration
!                       time filter to TKE in strong damping phase (ALSC)
!                       time filter to AHV (AFLT)
!     '04.01.21  A.Oka: ALPHC, ALPH, BETA, BETAC. Pr depends on Ri.
!     '04.01.27  A.Oka: two way for calculation of DML.
!     '04.02.23  H.Hasumi: imported from MIROC3.2
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: McDougall et al. (2003) eq. of state
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.19  T.Suzuki: bug fix
!                          (TKE not initialized when read from RSTO)
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.06.29  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim,  ntdim, nxydim, nxyzdm, &
    &   kstr,   kend,     kz,     nz, &
    &  ijstr,  ijend, &
    &     le,     lw,     ln,     ls,    lne,    lsw, &
    &  oinit, ofinal
  use zocgrd,  only: &
    &     dz,    dzm,     ds,    dsm,     dt, &
    &   zbot,    cor,   itst, ieuler
  use zocmsk,  only: &
#ifdef OPT_BBL
    & amsktb, amskvb, &
#endif
    &  amskt,  amftz,  amfvz, &
    &   nbot,  nbotv
  use zocphy, only: &
    &   rhoo, gravit,  ckarm

  implicit none
  private

  public :: vdiff
#ifdef OPT_BBL
  public :: vdiffb
#endif

contains

subroutine vdiff( &
  &                  amv,    ahv, &
  &                   uy,     vy,      r,   taux,   tauy, &
  &                   ty,     hy )

  use xprst
  use brstt
  use ufile

  real(8), intent(out) ::     amv(nxydim, nzdim),    ahv(nxydim, nzdim)
  real(8), intent(in)  ::      uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::       r(nxydim, nzdim)
  real(8), intent(in)  ::    taux(nxydim)       ,   tauy(nxydim)
  real(8), intent(in)  ::      ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::      hy(nxydim)

  real(8), save ::     tke(nxydim, nzdim)
  real(8), save ::  ahvbak(nxydim, nzdim)
  real(8), save ::  ahv03d(nxydim, nzdim)
  real(8), save ::  alplat(nxydim),  c0lat(nxydim)

  real(8) ::    drdz(nxydim, nzdim),  duvdz(nxydim, nzdim)
  real(8) ::   dzsig(nxydim, nzdim), rzmsig(nxydim, nzdim)
  real(8) ::   depth(nxydim, nzdim),     pr(nxydim, nzdim)
  real(8) ::       c(nxydim, nzdim),   cdmp(nxydim, nzdim)
  real(8) ::  adefwd(nxydim, nzdim),    rit(nxydim, nzdim)
  real(8) ::     fez(nxydim, nzdim),    tls(nxydim, nzdim)
  real(8) ::  amvtmp(nxydim, nzdim),  diffz(nxydim, nzdim)
  real(8) ::      aa(nxydim, nzdim),     ab(nxydim, nzdim)
  real(8) ::      ac(nxydim, nzdim),    ade(nxydim, nzdim)
  real(8) ::     dml(nxydim)
!      common /work/ dzsig, rzmsig, depth, drdz, duvdz, tls, c, fez

  real(8), save :: cc0(nzdim), cc1(nzdim), cc2(nzdim)
  real(8), save :: cc3(nzdim), cc4(nzdim), cc5(nzdim), cc6(nzdim)
  real(8), save :: d0(nzdim), d1(nzdim), d2(nzdim), d3(nzdim), d4(nzdim)
  real(8), save :: d5(nzdim), d6(nzdim), d7(nzdim), d8(nzdim), d9(nzdim)
  real(8), save :: r0

  logical, save :: ofirst = .true.,   oeof

  real(8) ::      p,   alps,      q
  logical ::   odml(nxydim)

  real(8) ::   cort,  cor30,     pi,  omega,    bfq
  real(8) ::  slats,  slatn,  nlats,  nlatn,    lat
  real(8) ::    qtl,   zqtl, dmltmp
  real(8) ::     p1,     p2
  real(8) ::     tl,     sl
  real(8) ::     rl,    rlu
  real(8) ::   dudz,   dvdz
  real(8) ::   fkls,     ri
  real(8) ::  avrtx,  avrty
  real(8) ::     fc
  integer ::     ij,      k,   iitr
  integer ::     kt
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::  amv0(nz) = 0.d0,  ahv0(nz) = 0.d0
  real(8), save ::  rahv(nz) = 0.d0
  real(8), save ::  sm0 = 0.39d0,  c0 = 0.06d0,  c0eq = -999.0D0
  real(8), save ::  pr0 = 0.8d0,  sg = 1.95d0
  real(8), save ::  eps = 1.0d-5,  z0 = 1.0d2,  alph = 3.0d0,  cftke = 1.0d2
  real(8), save ::  alphc = -999.0d0,  beta = 0.0d0,  betac = -999.0d0
  real(8), save ::  pr1 = 0.5d0,  prmax = 20.0d0
  real(8), save ::  ahvb = 0.1d0,  amvmax = 1000.0d0,  aflt = 1.0d0
  real(8), save ::  alsc = 999.0d0,  ritc = 1.0d0
  real(8), save ::  latal = 30.0d0,  lateq = 5.0d0,  al = 1.0d0,  aleq =0.1d0
  real(8), save ::  latc0 = 15.0D0,  latceq = 5.0D0
  integer, save ::  mz = nz,  nitr = 1
  logical, save ::  oallat = .false.,  oc0lat = .false.

  namelist /nmvisv/ amv0
  namelist /nmdifv/ ahv0
  namelist /nmdfre/ rahv
  namelist /nmdvnk/ eps, z0, alph, alphc, beta, betac, &
    &               cftke, mz, nitr, pr0, pr1, prmax, &
    &               sm0, c0, c0eq, sg, ahvb, amvmax, aflt, &
    &               alsc, ritc, oallat, al, aleq, latal, lateq, &
    &               oc0lat, latc0, latceq

  if (oinit) then
     do k = 1, nzdim
        do ij = 1, nxydim
           tke(ij, k) = eps
        end do
     end do
#ifdef OPT_TRIPOLE
     call rstadd(tke, oeof, nxdim, nydim, nzdim, 'TKE', 'OCN', &
       &                                         1.d0,  0,  0 )
#else
     call rstadd(tke, oeof, nxdim, nydim, nzdim, 'TKE', 'OCN')
#endif
     return
  end if

  if (ofinal) then
     call finadd(tke, nxdim, nydim, nzdim, 'TKE', 'OCN')
     return
  end if

  if (ofirst) then
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
     read (ifpar, nmdvnk, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdvnk', istat)
     write(jfpar, nmdvnk)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmdfre, iostat=istat)
     call cstnml(jfpar, 'vdiff', 'nmdfre', istat)
     write(jfpar, nmdfre)

     if (alphc < 0.0d0) then
        alphc = alph
        write(jfpar, *) '   alphc --> ', alphc
     endif
     if (betac < 0.0d0) then
        betac = beta
        write(jfpar, *) '   betac --> ', betac
     endif

     call secoef( &
       &         cc0(kstr), cc1(kstr), cc2(kstr), cc3(kstr), &
       &         cc4(kstr), cc5(kstr), cc6(kstr), &
       &         d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
       &         d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr) )

     r0 = rhoo * 1.d3

! --- reduce the background diffusion around EQ. ---
!
     pi = atan( 1.d0 )*4.d0
     omega = 2.d0 * pi / 86400.d0
     cor30 = 2.d0 * omega * sin( pi*30.d0/180.d0 )
     bfq   = 5.24d-3        ! /sec

     do k = kstr, kend
        do ij = ijstr, ijend
           cort=(cor(ij)+cor(ij+lw)+cor(ij+lsw)+cor(ij+ls))*0.25d0
           cort=abs(cort)
           if(cort.gt.cor30) then
              ahv03d(ij,k)=ahv0(k-kstr+1)
           else
!              ahv03d(ij,k)=ahv0(k-kstr+1) &
!                   *cort*acosh(bfq/cort)/cor30/acosh(bfq/cor30)
              ahv03d(ij, k) = ahv0(k-kstr+1) * cort / cor30 * &
                &             log(  bfq / cort &
                &                 + sqrt(bfq**2 / cort**2 - 1.d0)) &
                &           / log(  bfq / cor30 &
                &                 + sqrt(bfq**2 / cor30**2 - 1.d0))
              ahv03d(ij,k)= max(ahv03d(ij, k), 1.d-2)
              ahv03d(ij,k)= ahv03d(ij,  k) * rahv(k-kstr+1) &
                &         + ahv0(k-kstr+1) * (1.0 - rahv(k-kstr+1)) 
           endif
        end do
     end do
!
     if (oallat) then
        write(jfpar, *) &
          &  '   : DML is calculated following MY ', &
          &  'with AL latitude dependency'
        slatn = -lateq
        slats = -latal
        nlats = +lateq
        nlatn = +latal
        do ij = ijstr, ijend
           cort = (cor(ij)+cor(ij+lw)+cor(ij+lsw)+cor(ij+ls)) * 0.25d0
           if ( cort/(2.d0*omega) .gt.  1.d0 ) then
              cort =   2.d0 * omega
           end if
           if ( cort/(2.d0*omega) .lt. -1.d0 ) then
              cort = - 2.d0 * omega
           end if

           lat = asin( cort/2.d0/omega ) * 180.d0 / pi
           if ( (lat >= nlatn) .or. (lat <= slats) ) then
              alplat(ij) = al
           elseif ( (lat >= slats) .and. (lat < slatn) ) then
              alplat(ij) = (al*(slatn-lat) &
                &        + aleq*(lat-slats)) / (slatn-slats)
           elseif ( (lat > nlats) .and. (lat <= nlatn) ) then
              alplat(ij) = (al*(lat-nlats) &
                &        + aleq*(nlatn-lat)) / (nlatn-nlats)
           else
              alplat(ij) = aleq
           endif
        enddo
     else
        write(jfpar, *) &
          &  '   : DML is calculated following Noh et al.(2002)'
        do ij=ijstr, ijend
           alplat(ij) = al  !! dummy
        enddo
     endif

!    ---- increasing TKE dissipation coefficient around EQ.
     if (oc0lat) then
        write(jfpar, *) &
          &  '   : C0 around the equator is enhanced'
        slatn = -latceq
        slats = -latc0
        nlats = +latceq
        nlatn = +latc0
        do ij = ijstr, ijend
           cort = (cor(ij)+cor(ij+lw)+cor(ij+lsw)+cor(ij+ls)) * 0.25d0
           if ( cort/(2.d0*omega) .gt.  1.d0 ) then
              cort =   2.d0 * omega
           end if
           if ( cort/(2.d0*omega) .lt. -1.d0 ) then
              cort = - 2.d0 * omega
           end if
           lat = asin( cort/2.0/omega ) * 180.d0 / pi
           if ( lat.ge.nlatn .or. lat.le.slats ) then
              c0lat(ij)  = c0
           elseif ( lat.ge.slats .and. lat.lt.slatn ) then
              c0lat(ij) = (c0*(slatn-lat) &
                &       + c0eq*(lat-slats))/(slatn-slats)
           elseif ( lat.gt.nlats .and. lat.le.nlatn ) then
              c0lat(ij) = (c0*(lat-nlats) &
                &       + c0eq*(nlatn-lat))/(nlatn-nlats)
           else
              c0lat(ij) = c0eq
           endif
        enddo
     else
        do ij = ijstr, ijend
           c0lat(ij) = c0
        enddo
     endif
  end if

! -- second step of Euler-Backward sheme --
  if ( ( itst == 1 ).and.( ieuler == 2 ) ) then
     return
  endif

  if (ofirst) then
     do k = 1, nzdim
        do ij = 1, nxydim
           ahvbak(ij, k) = ahv(ij, k)
        end do
     end do
  end if

  do k = 1, nzdim
     do ij = 1, nxydim
        depth (ij, k) = 0.d0
        fez   (ij, k) = 0.d0
     end do
  end do
  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        dzsig (ij, k) = (hy(ij) + zbot) * ds(k)
        rzmsig(ij, k) = 1.d0 / (hy(ij) + zbot) / dsm(k)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        dzsig (ij, k) = dz(ij, k)
        rzmsig(ij, k) = 1.d0 / dzm(ij, k)
     end do
  end do
  do k = kstr+1, kend
     do ij = 1, nxydim
        depth(ij, k) = depth(ij, k-1) + dzsig(ij, k-1)
     end do
  end do

  do k = kstr+1, kstr+mz-1
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 =       cc0(k) &
          &     + (cc1(k) + (cc2(k) + cc3(k) * tl) * tl) * tl &
          &     + (cc4(k) + cc5(k) * tl + cc6(k) * sl) * sl
        p2 =       d0(k) &
          &     + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &     + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &              + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2

        tl = ty(ij, k-1, 1)
        sl = ty(ij, k-1, 2)
        p1 =       cc0(k) &
          &     + (cc1(k) + (cc2(k) + cc3(k) * tl) * tl) * tl &
          &     + (cc4(k) + cc5(k) * tl + cc6(k) * sl) * sl 
        p2 =       d0(k) &
          &     + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &     + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &              + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlu = p1 / p2

        drdz(ij, k) = gravit / r0 * (rl - rlu) * rzmsig(ij, k)
        drdz(ij, k) = max(drdz(ij, k), 0.d0)
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

  if (oallat) then
!    -- DML is calculated following MY with latitude dependency --
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        qtl   = 0.0d0
        zqtl  = 0.0d0
        do k = kstr, kstr+mz-1
           qtl  = qtl  + sqrt ( 2.0d0 * tke(ij, k) ) * dzsig(ij, k)
           zqtl = zqtl + depth(ij, k) * &
             &           sqrt ( 2.0d0 * tke(ij, k) ) * dzsig(ij, k)
        end do
        dmltmp = alplat(ij) * zqtl / max( qtl, eps )
        kt = min( kstr+mz-1, max(kstr+1, nbot(ij)))
        dml(ij) = max( min( dmltmp, depth( ij, kt ) ), eps )
     end do
  else
!    -- DML is calculated following Noh et al.(2002) --
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        kt = min( kstr+mz-1, max(kstr+1, nbot(ij)))
        dml(ij) = depth( ij, kt )
        odml(ij) = .true.
     enddo
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           if (odml(ij)) then
              if ( ahv(ij, k) .le. max(ahvb, ahv0(k-kstr+1)) ) then
                 dml(ij) = depth(ij, k)
                 odml(ij)=.false.
              end if
           end if
        end do
     end do
  endif

!  -- initialize TKE -- 
  if (ofirst .and. oeof) then
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
!          ---- using definition of Kondoh et. (1978) for Pr first guess,
!               since the method of Noh et al. (2005) requires TKE. ---
           pr(ij, k) =  min ( prmax, &
             &                pr0 + 7.d0 * drdz(ij, k) &
             &                    / max( duvdz(ij, k), eps ) )
           fkls = ckarm * (depth(ij, k) + z0)
           tls(ij, k) = fkls / (1.d0 + fkls / dml(ij))
           tke(ij, k) = &
             &        (   ahv(ij, k) * pr0 &
             &          + sqrt( ahv(ij, k)**2 * pr0**2 &
             &                  + 4.0d0 * alph * sm0 * sm0 * &
             &                    drdz(ij, k) * tls(ij, k)**4 ) &
             &        ) * ahv(ij, k) * pr0 &
             &          / ( 4.0d0 * sm0 * sm0 * tls(ij, k) * tls(ij, k) )
           tke(ij, k) = max( tke(ij, k) * amftz(ij, k), eps)
        enddo
     enddo
  end if

  do iitr = 1, nitr

     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           fkls = ckarm * (depth(ij, k) + z0)
           tls(ij, k) = fkls / (1.d0 + fkls / dml(ij))
           rit(ij, k) = drdz(ij, k) * tls(ij, k) * tls(ij, k) &
             &        / tke(ij,k) * 0.5d0
!          ---- Prantle number by definition of Noh et al. (2005, GRL)
           pr(ij, k) =  min ( prmax, &
             &                pr0 * sqrt( 1.d0 + pr1 * rit(ij, k) )  )
!           pr(ij, k) =  min ( prmax, &
!             &                pr0 + pr1 * drdz(ij, k) &
!             &                      / max( duvdz(ij, k), eps ) )
        enddo
     enddo
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ri = rit(ij, k) * exp( - min(300.0d0, &
             &                          beta / max( rit(ij, k), eps )))
           q  = sqrt( 2.d0 * tke(ij, k) )
           fc = tls(ij, k) * sm0 / sqrt( 1.0d0 + alph * ri )
           amv(ij, k) = min( amvmax, q * fc )

!           --- Forward for shear and buoyancy ---
!           c(ij, k)      =   4.d0 * c0 * tke(ij, k) / tls(ij, k) / ri
!           adefwd(ij, k) = - 2.d0 * sm0 * tls(ij, k) / ri * &
!             &             ( drdz(ij, k) / pr - duvdz(ij, k)) * &
!             &             tke(ij, k) * amftz(ij, k)
!           --- Semi-implicit for shear and buoyancy term ---
!           cdmp(ij, k) = ( 4.d0 * c0 * tke(ij, k) / tls(ij, k) / ri &
!             &           ) * amftz(ij, k)
           cdmp(ij, k) = &
             &         ( c0lat(ij) * sqrt( 1.0d0 + alphc * ri ) &
             &              * 2.0d0 * q / tls(ij, k) &
             &         ) * amftz(ij, k)
           c(ij, k)    = &
             &         ( fc * 2.0d0 / q &
             &              * (drdz(ij, k) / pr(ij, k) - duvdz(ij, k)) &
             &         ) * amftz(ij, k)
           adefwd(ij, k) =  0.0d0
        end do
     end do

!    -- semi-implicit : adjusted ALPS in each time step --
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           p = c(ij, k) * dt
!          -- foward after initial-growth phase --
           if ( (( rit(ij, k) < ritc) &
             &  .and. ( c(ij,k) < 0.0d0) ) &
             &  .or. ( p < -300.0d0 ) ) then
              alps  =  1.0d0
!          -- implicit for strong damping phase --
!           elseif ( p .ge. 1.0d0 ) then
!              alps  =  0.0d0
!          -- semi-implicit for initial-growth and damping phase --
           else
              p = sign(1.0d0, p) * max(1.0d-6, abs(p))
              alps  =  ( -1.0d0 + (1.0d0+p) * exp(-p) ) &
                &      / ( p * ( exp(-p) - 1.0d0 ) )
           endif
           adefwd(ij, k) = adefwd(ij, k) &
             &             - alps * c(ij, k) * tke(ij, k)
           c(ij, k)      = (1.0d0 - alps) * c(ij, k) &
             &             + cdmp(ij, k)
        end do
     end do

     do k = 1, nzdim
        do ij = 1, nxydim
           diffz(ij, k) = 0.d0
        end do
     end do
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        avrtx = (  taux(ij)    + taux(ij+lw) &
          &      + taux(ij+ls) + taux(ij+lsw)) * 0.25d0
        avrty = (  tauy(ij)    + tauy(ij+lw) &
          &      + tauy(ij+ls) + tauy(ij+lsw)) * 0.25d0
        fez(ij, kstr+1) = ((avrtx * avrtx + avrty * avrty) &
          &               / rhoo / rhoo) ** 0.75d0 * cftke
     end do
     do k = kstr+2, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           diffz(ij, k) = (amv(ij, k-1) + amv(ij, k)) * &
             &            0.5d0 / sg / dzsig(ij, k-1) * &
             &            amskt(ij, k-1)
           fez(ij, k) = diffz(ij, k) * (tke(ij, k-1) - tke(ij, k))
        end do
     end do
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ade(ij, k) = (fez(ij, k) - fez(ij, k+1)) * rzmsig(ij, k) &
             &          - c(ij, k) * tke(ij, k) + adefwd(ij, k)
        end do
     end do
     
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           aa(ij, k) = - dt * diffz(ij, k) * rzmsig(ij, k)
           ac(ij, k) = - dt * diffz(ij, k+1) * rzmsig(ij, k)
           ab(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k) &
             &         + dt * c(ij, k)
        end do
     end do
     do ij = ijstr-nxdim-1, ijend+nxdim+1
        ac(ij, kstr+1) = ac(ij, kstr+1) / ab(ij, kstr+1)
        ade(ij, kstr+1) = ade(ij, kstr+1) / ab(ij, kstr+1)
     end do
     do k = kstr+2, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           fc = 1.d0 / (ab(ij, k) - aa(ij, k) * ac(ij, k-1))
           ac(ij, k) = ac(ij, k) * fc
           ade(ij, k) = (ade(ij, k) - aa(ij, k) * ade(ij, k-1)) * fc
        end do
     end do
     do k = kstr+mz-2, kstr+1, -1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
           ade(ij, k) = ade(ij, k) - ac(ij, k) * ade(ij, k+1)
        end do
     end do
     
!    -- time filter to TKE in case of strong damping --
     do k = kstr+1, kstr+mz-1
        do ij = ijstr-nxdim-1, ijend+nxdim+1
!           if (cdmp(ij, k)*dt .ge. alsc) then
!              tke(ij, k) = (  0.5d0 * tke(ij, k) &
!                &           + 0.5d0 * (tke(ij, k) + dt * ade(ij, k)) &
!                &          ) * amftz(ij, k)
!           else
           tke(ij, k) = ( tke(ij, k) + dt * ade(ij, k) ) * amftz(ij, k)
!           endif
        end do
     end do
     do k = kstr, kend
        do ij = 1, nxydim
           tke(ij, k) = max(tke(ij, k), eps)
        end do
     end do

  end do

  if (ofirst) then
     nitr = 1
     ofirst = .false.
  end if

! -- smoothing --
  do k = kstr+1, kstr+mz-1
     do ij = ijstr-nxdim-1, ijend
        amvtmp(ij, k) = (  amv(ij   , k) + amv(ij+le , k) &
          &              + amv(ij+ln, k) + amv(ij+lne, k) &
          &             ) * 0.25d0 * amfvz(ij, k)
     end do
     do ij = ijstr, ijend
        amv(ij, k) = amvtmp(ij, k)
        ahv(ij, k) = (  amvtmp(ij   , k) + amvtmp(ij+lw , k) &
          &           + amvtmp(ij+ls, k) + amvtmp(ij+lsw, k)) * &
          &          0.25d0 / pr(ij, k)
     end do
  end do
  
  do k = kstr, kend
     do ij = ijstr, ijend
        amv(ij, k) = min(amvmax, max(amv(ij, k), amv0(k-kstr+1)))
        ahv(ij, k) = max(ahv(ij, k), ahv03d(ij,k))
     end do
  end do

! -- TIME FILTER --
  do k = kstr+1, kstr+mz-1
     do ij = ijstr, ijend
        ahv(ij, k) = (  aflt * ahv(ij, k) &
          &               + ( 1.0d0 - aflt ) * ahvbak(ij, k) &
          &          ) * amftz(ij, k)
     end do
  end do

! -- MEMORY PREVIOUS AMV --
  do k = kstr+1, kstr+mz-1
     do ij = ijstr, ijend
        ahvbak(ij, k) = ahv(ij, k)
     end do
  end do
  
#ifdef OPT_TRIPOLE
  call shift1(   tke, &
    &          nxdim,  nydim,  nzdim, &
    &           1.d0,      0,      0 )
#endif
   
  return
end subroutine vdiff
#ifdef OPT_BBL
! *********************************************************************

subroutine vdiffb( &
  &                   amv,    ahv )

! --- information -----------------------------------------------------
!
!  Vertical viscosity and diffusion coefficients for the bottom
! boundary layer.
!
!  HISTORY
!     '01.02.08  H.Hasumi
!     '12.06.29  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------
  use ufile

  real(8), intent(inout) ::    amv(nxydim, nzdim),    ahv(nxydim, nzdim)

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
     amv(ij, kv) = max(amvbbl, amv(ij, kv)) * amskvb(ij) &
       &         + amv(ij, kv) * (1.d0 - amskvb(ij))
     ahv(ij, kt) = max(ahvbbl, ahv(ij, kt)) * amsktb(ij) &
       &         + ahv(ij, kt) * (1.d0 - amsktb(ij))
  end do

  do ij = ijstr, ijend
     kt = nbot(ij)
     kv = nbotv(ij)
     amv(ij, kend) = amv(ij, kv)
     ahv(ij, kend) = ahv(ij, kt)
  end do

  return
end subroutine vdiffb
#endif

end module dvdif
