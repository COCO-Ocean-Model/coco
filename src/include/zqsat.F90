! PACKAGE ZQSAT     !! Clausis-Clapeyron
!
!     QSAT:  saturation water vapour mixing ratio
!     DQSAT: d(QSAT)/d(T)
!
! The following use statement is needed before this
! use zocphy, only:   epsv,    es0,     el,  emelt,  tqice,   rvap,  tmelt
!
  real(8) ::    fqsat
  real(8) ::    fdqsat
#ifdef OPT_FQSATL
  real(8) ::    fqsatl
  real(8) ::    fqsati
  real(8) ::    fqsatr
#endif
  real(8) ::    t, p, qs
!
  fqsat ( t,p )  = epsv * es0 / p &
    &            * exp( (el+emelt/2.d0*(1.d0-sign(1.d0,t-tqice))) &
    &                   /rvap *( 1.d0/tmelt - 1.d0/t )           )
!
  fdqsat( t,qs ) = (el+emelt/2.d0*(1.d0-sign(1.d0,t-tmelt))) &
    &            * qs / ( rvap * t*t )
!
#ifdef OPT_FQSATL
  fqsatl( t,p )  = epsv * es0 / p &
    &            * exp( el/rvap *( 1.d0/tmelt - 1.d0/t ) )
!
  fqsati( t,p )  = epsv * es0 / p &
    &            * exp( (el+emelt)/rvap *( 1.d0/tmelt - 1.d0/t ) )
!
  fqsatr( t   )  = exp( -emelt/rvap *( 1.d0/tmelt - 1.d0/t ) )
#endif
