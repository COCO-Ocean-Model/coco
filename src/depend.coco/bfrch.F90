module bfrch

! --- information -----------------------------------------------------
!
!  HISTORY
!     '99.08.10  H.Hasumi
!     '02.05.29  H.Nakano: tracer dimension
!     '07.04.23  H.Hasumi
!     '12.11.23  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  implicit none
  private

  public  ::  forsto,  excngo

contains

  
  subroutine forsto(  ua,     va,     ta,                             &
      &               ha,   ubta,   vbta,                             &
      &               ub,     vb,     tb,                             &
      &               hb,   ubtb,   vbtb   )

    use zocdim,  only  :  nxyzdm,  nxydim,  ntdim

    implicit none

    real(8),        intent(inout)  ::    ua(1:nxyzdm),   va(1:nxyzdm)
    real(8),        intent(inout)  ::    ta(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    ha(1:nxydim)
    real(8),        intent(inout)  ::  ubta(1:nxydim), vbta(1:nxydim)
    real(8),        intent(in)     ::    ub(1:nxyzdm),   vb(1:nxyzdm)
    real(8),        intent(in)     ::    tb(1:nxyzdm,1:ntdim)
    real(8),        intent(in)     ::    hb(1:nxydim)
    real(8),        intent(in)     ::  ubtb(1:nxydim), vbtb(1:nxydim)

    integer(4)  ::     ij,      ijk,     n 

    do ijk = 1, nxyzdm
       ua(ijk) = ub(ijk)
       va(ijk) = vb(ijk)
    end do

    do n = 1, ntdim
       do ijk = 1, nxyzdm
          ta(ijk, n) = tb(ijk, n)
       end do
    end do

    do ij  = 1, nxydim
       ha(ij)   = hb(ij)
       ubta(ij) = ubtb(ij)
       vbta(ij) = vbtb(ij)
    end do
    
  end subroutine forsto

! *********************************************************************

  subroutine excngo(  ua,     va,     ta,                             &
      &               ha,   ubta,   vbta,                             &
      &               ub,     vb,     tb,                             &
      &               hb,   ubtb,   vbtb   )

    use zocdim,  only  :  nxyzdm,  nxydim,  ntdim

    implicit none

    real(8),        intent(inout)  ::    ua(1:nxyzdm),   va(1:nxyzdm)
    real(8),        intent(inout)  ::    ta(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    ha(1:nxydim)
    real(8),        intent(inout)  ::  ubta(1:nxydim), vbta(1:nxydim)
    real(8),        intent(inout)  ::    ub(1:nxyzdm),   vb(1:nxyzdm)
    real(8),        intent(inout)  ::    tb(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    hb(1:nxydim)
    real(8),        intent(inout)  ::  ubtb(1:nxydim), vbtb(1:nxydim)

    real(8)     ::      u,      v,      t
    real(8)     ::     sh,    ubt,    vbt
    integer(4)  ::     ij,     ijk,     n 

    do ijk = 1, nxyzdm
       u       = ub(ijk)
       v       = vb(ijk)
       ub(ijk) = ua(ijk)
       vb(ijk) = va(ijk)
       ua(ijk) = u
       va(ijk) = v
    end do

    do n = 1, ntdim
       do ijk = 1, nxyzdm
          t          = tb(ijk, n)
          tb(ijk, n) = ta(ijk, n)
          ta(ijk, n) = t
       end do
    end do

    do ij  = 1, nxydim
       sh       = hb(ij)
       ubt      = ubtb(ij)
       vbt      = vbtb(ij)
       hb(ij)   = ha(ij)
       ubtb(ij) = ubta(ij)
       vbtb(ij) = vbta(ij)
       ha(ij)   = sh
       ubta(ij) = ubt
       vbta(ij) = vbt
    end do

  end subroutine excngo

end module bfrch
