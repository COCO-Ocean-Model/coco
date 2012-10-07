module bfrch
! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Nakano: from COCO3.4
!     '07.04.23  H.Hasumi
!     '12.10.07  M.Kurogi: rewrite in F95 format
!
! ---------------------------------------------------------------------
contains
 subroutine forsto(                             &
  &                      ua,     va,     ta,      &
  &                      ha,   ubta,   vbta,      &
  &                      ub,     vb,     tb,      &
  &                      hb,   ubtb,   vbtb)
 use zocdim, only : nxydim, nxyzdm, ntdim
 implicit none
 real(8), intent(out) ::  ua(nxyzdm),   va(nxyzdm),   ta(nxyzdm, ntdim)
 real(8), intent(out) ::  ha(nxydim), ubta(nxydim), vbta(nxydim)

 real(8), intent(in ) ::  ub(nxyzdm),   vb(nxyzdm),   tb(nxyzdm, ntdim)
 real(8), intent(in)  ::  hb(nxydim), ubtb(nxydim), vbtb(nxydim)

 integer ::    ij,    ijk,      n

 do ijk = 1, nxyzdm
    ua(ijk) = ub(ijk)
    va(ijk) = vb(ijk)
 end do
 do  n = 1, ntdim
    do  ijk = 1, nxyzdm
       ta(ijk, n) = tb(ijk, n)
    end do
 end do
 do  ij  = 1, nxydim
    ha(ij)   = hb(ij)
    ubta(ij) = ubtb(ij)
    vbta(ij) = vbtb(ij)
 end do
 return
 end subroutine forsto

! *********************************************************************

 subroutine excngo(                           &
  &                      ua,     va,     ta,  &
  &                      ha,   ubta,   vbta,  &
  &                      ub,     vb,     tb,  &
  &                      hb,   ubtb,   vbtb)
 use zocdim, only : nxydim, nxyzdm, ntdim
 implicit none

 real(8), intent(inout) ::      ta(nxyzdm, ntdim),     tb(nxyzdm, ntdim)
 real(8), intent(inout) ::      ua(nxyzdm),     ub(nxyzdm)
 real(8), intent(inout) ::      va(nxyzdm),     vb(nxyzdm)
 real(8), intent(inout) ::      ha(nxydim),     hb(nxydim)
 real(8), intent(inout) ::    ubta(nxydim),   ubtb(nxydim)
 real(8), intent(inout) ::    vbta(nxydim),   vbtb(nxydim)

 real(8) ::       u,      v,      t
 real(8) ::      sh,    ubt,    vbt
 integer ::    ij,    ijk,      n

 do  ijk = 1, nxyzdm
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

 return
 end subroutine excngo

! *********************************************************************

 subroutine forsti(                                             &
  &                aa,    hia,    uia,    via,    tia,    hsa,  &
  &                ab,    hib,    uib,    vib,    tib,    hsb)
 use zocdim, only : nxydim, nic
 implicit none
 real(8), intent(out) ::  aa(nxydim, 0:nic),  hia(nxydim, 0:nic)
 real(8), intent(out) :: uia(nxydim),         via(nxydim)
 real(8), intent(out) :: tia(nxydim, 0:nic),  hsa(nxydim, 0:nic)    

 real(8), intent(in ) ::  ab(nxydim, 0:nic), hib(nxydim, 0:nic)
 real(8), intent(in ) :: uib(nxydim),        vib(nxydim)
 real(8), intent(in ) :: tib(nxydim, 0:nic), hsb(nxydim, 0:nic)

 integer     ij,      k

 do ij = 1, nxydim
    uia(ij) = uib(ij)
    via(ij) = vib(ij)
 end do

 do k = 0, nic
    do ij = 1, nxydim
       aa (ij, k) = ab (ij, k)
       hia(ij, k) = hib(ij, k)
       hsa(ij, k) = hsb(ij, k)
       tia(ij, k) = tib(ij, k)
    end do
 end do

 return
 end subroutine forsti

! *********************************************************************

 subroutine excngi(                                             &
  &                aa,    hia,    uia,    via,    tia,    hsa,  &
  &                ab,    hib,    uib,    vib,    tib,    hsb)
 use zocdim, only : nxydim, nic
 implicit none
      
 real(8), intent(inout) ::  aa(nxydim, 0:nic),     ab(nxydim, 0:nic)
 real(8), intent(inout) :: hia(nxydim, 0:nic),    hib(nxydim, 0:nic)
 real(8), intent(inout) :: uia(nxydim),           uib(nxydim)
 real(8), intent(inout) :: via(nxydim),           vib(nxydim)
 real(8), intent(inout) :: tia(nxydim, 0:nic),    tib(nxydim, 0:nic)
 real(8), intent(inout) :: hsa(nxydim, 0:nic),    hsb(nxydim, 0:nic)

 integer ::    ij,      k
 real(8) ::       a,     hi,     ui,     vi,     ti,     hs

 do ij = 1, nxydim
    ui      = uib(ij)
    vi      = vib(ij)
    uib(ij) = uia(ij)
    vib(ij) = via(ij)
    uia(ij) = ui
    via(ij) = vi
 end do
 do k = 0, nic
    do ij = 1, nxydim
       a          = ab (ij, k)
       hi         = hib(ij, k)
       hs         = hsb(ij, k)
       ti         = tib(ij, k)
       ab (ij, k) = aa (ij, k)
       hib(ij, k) = hia(ij, k)
       hsb(ij, k) = hsa(ij, k)
       tib(ij, k) = tia(ij, k)
       aa (ij, k) = a
       hia(ij, k) = hi
       hsa(ij, k) = hs
       tia(ij, k) = ti
    end do
 end do

 return
 end subroutine excngi
end module bfrch
