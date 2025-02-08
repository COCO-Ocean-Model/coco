module dvlva
! --- information -----------------------------------------------------
!
!  Velocity for momentum advection
!
!  HISTORY
!     '99.08.17  H.Hasumi
!     '00.03.17  H.Hasumi: for Ishizaki scheme
!        In this version, UADV and VADV are not advection velocity
!        for V-box, but velocity at T-points. WADV is nine-
!        direction vertical velocity: the direction is designated
!        by the last argument (1 for purely vertical, 2 for north-
!        ward, and rounds clockwise to 9 for northwestward).
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '01.02.20  H.Hasumi: bug fix (weight for UADV, VADV)
!                          for incorporating BBL model
!     '01.05.02  H.Hasumi
!     '01.05.17  H.Hasumi: bug fix in VELVAB
!     '02.10.17  T.Suzuki: diagonal momentum advection from BBL
!     '02.11.??  T.Suzuki: remove cross diagonal momentum advection
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.02.23  Y.Komuro: bug fix by Dr. Kurogi
!                          (loop length extended)
!     '10.04.14  M.kurogi
!     '10.04.14  M.Kurogi: for tripolar grid
!     '12.10.22  T.Suzuki: for COCO 5.0
!
! ----------------------------------------------------------------------
  private
  public :: velvad, velvab, excngw

contains
  
  subroutine velvad(                                                        &
    &               uadv,   vadv,   wadv,                                   &
    &                  u,      v,      w )
   use zocdim, only :                                                       &
    &  nxydim, nxdim, nzdim, kstr, kend, kz, ijtstr, ijtend, ijvstr, ijvend,&
    &  le, ln, lnw, lne, lw, ls, lsw, lse,                                  &
    &  oinit, ofinal

   use zocgrd, only : dy, dym, rym, dzv
   use zocmsk, only : amskt, amskv
   use bchmk

   implicit none

   real(8), intent(out) :: uadv(nxydim, nzdim),   vadv(nxydim, nzdim)
   real(8), intent(out) :: wadv(nxydim, nzdim,  9)
   real(8), intent(in)  ::    u(nxydim, nzdim),      v(nxydim, nzdim)
   real(8), intent(in)  ::    w(nxydim, nzdim)
   
   integer ::    ij,      k,      l
   real(8) ::    an

   real(8), save ::   rn(nxydim, nzdim),    rnn(nxydim, nzdim)
   real(8), save ::  plo(nxydim, nzdim),    rcr(nxydim, nzdim)

   logical, save :: ofirst = .true.

   if (oinit .or. ofinal) then
      return
   end if
 
   if (ofirst) then
      !$acc enter data create(rn, rnn, plo, rcr)
      ofirst = .false.
#ifdef OPT_BBL
      call rmmskv
#endif

    !$acc kernels default(present)
    do k = 1, nzdim
       do ij = 1, nxydim
          rn (ij, k) = 0.d0
          rnn(ij, k) = 0.d0
          plo(ij, k) = 0.d0
       end do
    end do

    do k = kstr, kend
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          an = amskv(ij   , k) + amskv(ij+lw , k)             &   
    &        + amskv(ij+ls, k) + amskv(ij+lsw, k)
       if (an .eq. 0.d0) then
          rn(ij, k) = 0.d0
       else
          rn(ij, k) = 1.d0 / an
       end if
       end do
    end do

    do k = kstr+1, kend
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          if (rn(ij, k) .eq. rn(ij, k-1)) then
             rnn(ij, k) = 0.d0
          else if (rn(ij, k-1) .eq. 0.d0) then
             rnn(ij, k) = 0.d0
          else
             rnn(ij, k) = rn(ij, k) * rn(ij, k-1)
          end if
       end do
    end do

    do k = kstr+1, kend
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          if (      (rn(ij,k) .eq. 0.5d0)                     &
    &         .and. (rn(ij,k-1) .eq. 0.25d0)) then
             rnn(ij, k) = 0.25d0
             rcr(ij, k) = 0.0d0
          else
             rcr(ij, k) = 1.0d0
          end if
       end do
    end do

    do k = kstr+1, kend
       do ij = ijvstr-nxdim-1, ijvend+nxdim+1
          if (      (amskv(ij, k) .eq. 0.d0)                  &
    &         .and. (amskv(ij, k-1) .ne. 0.d0)) then
                  plo(ij, k) = 1.d0
          end if
       end do
    end do
    !$acc end  kernels
#ifdef OPT_BBL
!         call admkv1
         call admskv
#endif
    end if

    !$acc kernels default(present)
    do k = 1, nzdim
       do ij = 1, nxydim
          uadv(ij, k) = 0.d0
          vadv(ij, k) = 0.d0
       enddo
    enddo
    do l = 1, 9
       do k = 1, nzdim
          do ij = 1, nxydim
             wadv(ij, k, l) = 0.d0
          enddo
       enddo
    enddo

    do k = kstr, kstr+kz-1
       do ij = ijtstr, ijtend
           uadv(ij, k) = 0.25d0 * (  u(ij+lsw, k) + u(ij+lw , k)   &
    &                              + u(ij+ls , k) + u(ij    , k))  
           vadv(ij, k) = 0.25d0 * (  v(ij+lsw, k) + v(ij+lw , k)   &
    &                              + v(ij+ls , k) + v(ij    , k))
       enddo
    enddo

    do k = kstr+kz, kend
       do ij = ijtstr, ijtend
          uadv(ij, k) = 0.25d0 * (  u(ij+lsw, k) * dzv(ij+lsw, k)     &
    &                             + u(ij+lw , k) * dzv(ij+lw , k)     &
    &                             + u(ij+ls , k) * dzv(ij+ls , k)     &
    &                             + u(ij    , k) * dzv(ij    , k)) *  &
    &                   amskt(ij, k)
          vadv(ij, k) = 0.25d0 * (  v(ij+lsw, k) * dzv(ij+lsw, k)     &
    &                             + v(ij+lw , k) * dzv(ij+lw , k)     &
    &                             + v(ij+ls , k) * dzv(ij+ls , k)     &
    &                             + v(ij    , k) * dzv(ij    , k)) *  &
    &                   amskt(ij, k)
       enddo
    enddo

    do k = kstr, kend
       do ij = ijvstr, ijvend
          wadv(ij, k, 1) = amskv(ij, k) *                     &
    &                     (  w(ij+lne, k) * rn(ij+lne, k-1)   &
    &                      + w(ij+ln , k) * rn(ij+ln , k-1)   &
    &                      + w(ij+le , k) * rn(ij+le , k-1)   &
    &                      + w(ij    , k) * rn(ij    , k-1)) 
       enddo
    enddo

    do k = kstr+1, kend
       do ij = ijvstr, ijvend
          wadv(ij, k, 2) = amskv(ij, k) * plo(ij+ln, k) *     &
    &                    (  w(ij+ln , k) * rnn(ij+ln , k)     &
    &                     + w(ij+lne, k) * rnn(ij+lne, k))    
          wadv(ij, k, 4) = amskv(ij, k) * plo(ij+le, k) *     &
    &                    (  w(ij+lne, k) * rnn(ij+lne, k)     &
    &                     + w(ij+le , k) * rnn(ij+le , k))
          wadv(ij, k, 6) = amskv(ij, k) * plo(ij+ls, k) *     &
    &                    (  w(ij+le , k) * rnn(ij+le , k)     &
    &                     + w(ij    , k) * rnn(ij    , k))
          wadv(ij, k, 8) = amskv(ij, k) * plo(ij+lw, k) *     &
    &                    (  w(ij    , k) * rnn(ij    , k)     &
    &                     + w(ij+ln , k) * rnn(ij+ln , k))
          wadv(ij, k, 3) = amskv(ij, k) * plo(ij+lne, k) *    &
    &                    w(ij+lne, k) * rnn(ij+lne, k) *      &
    &                    rcr(ij+lne, k)
          wadv(ij, k, 5) = amskv(ij, k) * plo(ij+lse, k) *    &
    &                    w(ij+le , k) * rnn(ij+le , k) *      &
    &                    rcr(ij+le, k)
          wadv(ij, k, 7) = amskv(ij, k) * plo(ij+lsw, k) *    &
    &                    w(ij    , k) * rnn(ij    , k) *      &
    &                    rcr(ij, k)
          wadv(ij, k, 9) = amskv(ij, k) * plo(ij+lnw, k) *    &
    &                   w(ij+ln , k) * rnn(ij+ln , k) *       &
    &                   rcr(ij+ln, k)
       enddo
    enddo
    !$acc end  kernels
    
    return
 end subroutine velvad
#ifdef OPT_BBL
! *********************************************************************
  subroutine velvab(             &
     &               wadv,       &
     &                  w)

    use zocdim
    use zocgrd
    use zocmsk

    implicit none

    real(8), intent(inout) ::    wadv(nxydim, nzdim, 9)
    real(8), intent(in)    ::       w(nxydim, nzdim)

    real(8) ::   an
    real(8), save ::     rb(nxydim)
    real(8), save ::    rnb(nxydim, nzdim)
    real(8), save ::   rnnb(nxydim, nzdim)
    real(8), save ::   rcrb(nxydim)
    integer ::    ij,      k,    kup
 
    logical, save :: ofirst = .true.

    if (oinit .or. ofinal) then
       return
    end if

    if (ofirst) then
       ofirst = .false.

       do ij = 1, nxydim
          rb(ij) = 0.d0
       end do
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          an = amskvb(ij)    + amskvb(ij+lw)   &
       &     + amskvb(ij+ls) + amskvb(ij+lsw)
          if (an .eq. 0.d0) then
             rb(ij) = 0.d0
          else
             rb(ij) = 1.d0 / an
          end if
       end do

       do k = 1, nzdim
          do ij = 1, nxydim
             rnb (ij, k) = 0.d0
             rnnb(ij, k) = 0.d0
          end do
       end do
       do ij = 1, nxydim
          rcrb(ij) = 1.d0
       end do

       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          k = nbot(ij) - amsktb(ij)
          rnb(ij, k) = 1.d0
          an = amskv(ij,     k)  + amskv(ij+lw,  k)  &
       &     + amskv(ij+ls,  k)  + amskv(ij+lsw, k)
          if (an .eq. 0.d0) then
             rnnb(ij, k) = 0.d0
          else
             rnnb(ij, k) = rb(ij) / an
          end if
          if (an .eq. 2.d0) then
             rcrb(ij) = 0.d0
             rnnb(ij, k) = rb(ij) 
          else
             rcrb(ij) = 1.d0
          end if
       end do
    end if

    do ij = ijvstr, ijvend
          kup = nbotv(ij) - amskvb(ij)
          wadv(ij, kend, 1) = amskvb(ij) *  &
    &             (  w(ij+lne, kend) * rb(ij+lne) * rnb(ij+lne, kup)  &
    &              + w(ij+ln , kend) * rb(ij+ln ) * rnb(ij+ln,  kup)  &
    &              + w(ij+le , kend) * rb(ij+le ) * rnb(ij+le,  kup)  &
    &              + w(ij    , kend) * rb(ij    ) * rnb(ij,     kup))

          k = nbotv(ij+ln) - amskvb(ij+ln)
          wadv(ij, kend, 2) = amskvb(ij) * amskv(ij+ln, k) *               &
    &     (  w(ij+ln , kend) * rnnb(ij+ln , k) * (1.d0 - rnb(ij+ln,  kup)) &
    &      + w(ij+lne, kend) * rnnb(ij+lne, k) * (1.d0 - rnb(ij+lne, kup)))
          k = nbotv(ij+le) - amskvb(ij+le)
          wadv(ij, kend, 4) = amskvb(ij) * amskv(ij+le, k) *               &
    &     (  w(ij+lne, kend) * rnnb(ij+lne, k) * (1.d0 - rnb(ij+lne, kup)) &
    &     + w(ij+le , kend) * rnnb(ij+le , k) * (1.d0 - rnb(ij+le,  kup)))
          k = nbotv(ij+ls) - amskvb(ij+ls)
          wadv(ij, kend, 6) = amskvb(ij) * amskv(ij+ls, k) *               &
    &     (  w(ij+le , kend) * rnnb(ij+le , k) * (1.d0 - rnb(ij+le,  kup)) &
    &     + w(ij    , kend) * rnnb(ij    , k) * (1.d0 - rnb(ij   ,  kup)))
          k = nbotv(ij+lw) - amskvb(ij+lw)
          wadv(ij, kend, 8) = amskvb(ij) * amskv(ij+lw, k) *               &
    &     (  w(ij    , kend) * rnnb(ij    , k) * (1.d0 - rnb(ij   ,  kup)) &
    &      + w(ij+ln , kend) * rnnb(ij+ln , k) * (1.d0 - rnb(ij+ln,  kup)))

          k = nbotv(ij+lne) - amskvb(ij+lne)
          wadv(ij, kend, 3) = amskvb(ij) * amskv(ij+lne, k)* rcrb(ij+lne) * &
    &        w(ij+lne, kend) * rnnb(ij+lne, k) * (1.d0 - rnb(ij+lne, kup)) 
          k = nbotv(ij+lse) - amskvb(ij+lse)
          wadv(ij, kend, 5) = amskvb(ij) * amskv(ij+lse, k)* rcrb(ij+le ) * &
    &        w(ij+le , kend) * rnnb(ij+le , k) * (1.d0 - rnb(ij+le,  kup))
          k = nbotv(ij+lsw) - amskvb(ij+lsw)
          wadv(ij, kend, 7) = amskvb(ij) * amskv(ij+lsw, k)* rcrb(ij    ) * &
    &        w(ij    , kend) * rnnb(ij    , k) * (1.d0 - rnb(ij,     kup))
          k = nbotv(ij+lnw) - amskvb(ij+lnw)
          wadv(ij, kend, 9) = amskvb(ij) * amskv(ij+lnw, k)* rcrb(ij+ln ) * &
    &        w(ij+ln , kend) * rnnb(ij+ln , k) * (1.d0 - rnb(ij+ln,  kup))
    end do

    do ij = ijvstr, ijvend
       k = nbotv(ij)
       wadv(ij, k, 1) = wadv(ij, k, 1) * (1.d0 - amskvb(ij)) &
    &                 + wadv(ij, kend, 1) * amskvb(ij)
    end do

  return
  end subroutine velvab
#endif


#ifdef OPT_TRIPOLE
  subroutine excngw(wadv)
    use zocdim
    implicit none

    real(8), intent(inout) ::    wadv(nxdim, nydim, nzdim, 9)
    real(8) ::  wadv2(nxdim, nydim, nzdim, 9)

    integer i, j, k, n, istv, jstv
    logical, save :: ofirst = .true.

    if (ofirst) then
       !$acc enter data create(wadv2)
       ofirst = .false.
    end if    
    
    if(jrank .ne. jnodes-1) return

    !$acc kernels default(present)
    do n=1,9
    do k=1,nzdim
    do j=1,nydim
    do i=1,nxdim
       wadv2(i,j,k,n)=wadv(i,j,k,n)
    end do
    end do
    end do
    end do
    !$acc end kernels
    if(inodes .eq. 1) then
       !$acc kernels default(present)
      do k=1, nzdim
      do j=jend+1, nydim
      do i=1,nxdim
        wadv(i,j,k,2)=wadv2(i,j,k,6)
        wadv(i,j,k,3)=wadv2(i,j,k,7)
        wadv(i,j,k,4)=wadv2(i,j,k,8)
        wadv(i,j,k,5)=wadv2(i,j,k,9)
        wadv(i,j,k,6)=wadv2(i,j,k,2)
        wadv(i,j,k,7)=wadv2(i,j,k,3)
        wadv(i,j,k,8)=wadv2(i,j,k,4)
        wadv(i,j,k,9)=wadv2(i,j,k,5)
      end do
      end do
      end do

      j=jend
      do k=1, nzdim
      do i=nxdim/2+1,nxdim
        wadv(i,j,k,2)=wadv2(i,j,k,6)
        wadv(i,j,k,3)=wadv2(i,j,k,7)
        wadv(i,j,k,4)=wadv2(i,j,k,8)
        wadv(i,j,k,5)=wadv2(i,j,k,9)
        wadv(i,j,k,6)=wadv2(i,j,k,2)
        wadv(i,j,k,7)=wadv2(i,j,k,3)
        wadv(i,j,k,8)=wadv2(i,j,k,4)
        wadv(i,j,k,9)=wadv2(i,j,k,5)
      end do
      end do
      !$acc end kernels
    else if(irank .lt. inodes/2) then
      !$acc kernels default(present)
      do k=1, nzdim
      do j=jend+1, nydim
      do i=1,nxdim
        wadv(i,j,k,2)=wadv2(i,j,k,6)
        wadv(i,j,k,3)=wadv2(i,j,k,7)
        wadv(i,j,k,4)=wadv2(i,j,k,8)
        wadv(i,j,k,5)=wadv2(i,j,k,9)
        wadv(i,j,k,6)=wadv2(i,j,k,2)
        wadv(i,j,k,7)=wadv2(i,j,k,3)
        wadv(i,j,k,8)=wadv2(i,j,k,4)
        wadv(i,j,k,9)=wadv2(i,j,k,5)
      end do
      end do
      end do
      !$acc end kernels
    else
      !$acc kernels default(present)
      do k=1, nzdim
      do j=jend, nydim
      do i=1,nxdim
        wadv(i,j,k,2)=wadv2(i,j,k,6)
        wadv(i,j,k,3)=wadv2(i,j,k,7)
        wadv(i,j,k,4)=wadv2(i,j,k,8)
        wadv(i,j,k,5)=wadv2(i,j,k,9)
        wadv(i,j,k,6)=wadv2(i,j,k,2)
        wadv(i,j,k,7)=wadv2(i,j,k,3)
        wadv(i,j,k,8)=wadv2(i,j,k,4)
        wadv(i,j,k,9)=wadv2(i,j,k,5)
      end do
      end do
      end do
      !$acc end kernels
    end if
    
  return
  end subroutine excngw
#endif

end module dvlva
