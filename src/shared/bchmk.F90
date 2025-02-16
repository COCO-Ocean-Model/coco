module bchmk
! --- information -----------------------------------------------------
!
!  Change the values of masking array for BBL grids.
!
!  HISTORY
!     '01.02.08  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.10.07  M.Kurogi: rewrite in F95 format
!
! ---------------------------------------------------------------------
contains
#ifdef OPT_BBL
 subroutine rmmskv
 use zocdim, only : nxydim, kend
 use zocmsk, only : amskv, amskv0, nbotv
 implicit none
 integer ::  ij

 do ij = 1, nxydim
    amskv(ij, kend) = 0.d0
    amskv(ij, nbotv(ij)) = amskv0(ij)
 end do

 return
 end subroutine rmmskv

! *********************************************************************

 subroutine admkv1
 use zocdim, only : nxydim
 use zocmsk, only : amskv, nbotv, amskv1
 implicit none
 integer ::  ij

 do ij = 1, nxydim
    amskv(ij, nbotv(ij)) = amskv1(ij)
 end do

 return
 end subroutine admkv1

! *********************************************************************

 subroutine admkvb
 use zocdim, only : nxydim, kend
 use zocmsk, only : amskv, amskvb
 implicit none
 integer ::  ij

 do ij = 1, nxydim
    amskv(ij, kend) = amskvb(ij)
 end do

 return
 end subroutine admkvb

! *********************************************************************

 subroutine admskv
 use zocdim, only : nxydim, kend
 use zocmsk, only : amskv, amskvb, amskv1, nbotv
 implicit none
 integer ::  ij

 do ij = 1, nxydim
    amskv(ij, kend) = amskvb(ij)
    amskv(ij, nbotv(ij)) = amskv1(ij)
 end do

 return
 end subroutine admskv

! *********************************************************************

 subroutine rmmskt
 use zocdim, only : nxydim, kend
 use zocmsk, only : amskt, amskt0, nbot
 implicit none
 integer ::  ij

#ifdef ACC_
!$acc data copyin(amskt0, nbot)
!$acc data copy(amskt)
!$acc kernels
#endif
 do ij = 1, nxydim
    amskt(ij, kend) = 0.d0
    amskt(ij, nbot(ij)) = amskt0(ij)
 end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copyin(amskt0)
!$acc end data ! copy(amskt)
#endif

 return
 end subroutine rmmskt

! *********************************************************************

 subroutine admktb
 use zocdim, only : nxydim, kend
 use zocmsk, only : amskt, amsktb
 implicit none
 integer ::  ij

#ifdef ACC_
!$acc data copyin(amsktb)
!$acc data copy(amskt)
!$acc kernels
#endif
 do ij = 1, nxydim
    amskt(ij, kend) = amsktb(ij)
 end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copyin(amsktb)
!$acc end data ! copy(amskt)
#endif

 return
 end subroutine admktb

! *********************************************************************

 subroutine admkt1
 use zocdim, only : nxydim
 use zocmsk, only : amskt, amskt1, nbot
 implicit none
 integer ::  ij

#ifdef ACC_
!$acc data copyin(amskt1)
!$acc data copyin(nbot)
!$acc data copy(amskt)
!$acc kernels
#endif
 do ij = 1, nxydim
    amskt(ij, nbot(ij)) = amskt1(ij)
 end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(amskt)
!$acc end data ! copyin(amskt1)
!$acc end data ! copyin(nbot)
#endif

 return
 end subroutine admkt1

! *********************************************************************

 subroutine admskt
 use zocdim, only : nxydim, kend
 use zocmsk, only : amskt, amsktb, amskt1, nbot
 implicit none
 integer ::  ij

 do ij = 1, nxydim
    amskt(ij, kend) = amsktb(ij)
    amskt(ij, nbot(ij)) = amskt1(ij)
 end do

 return 
 end subroutine admskt
#else
! *********************************************************************
 subroutine admskt
 implicit none
 return
 end subroutine admskt
#endif

end module bchmk
