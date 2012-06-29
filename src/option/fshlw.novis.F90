module fshlw

! ---- information ----------------------------------------------------
!
!  HISTORY
!     '03.04.23  H.Hasumi: from COCO3.4
!     '06.01.10  H.Hasumi: bug fix (calculation of SYY)
!     '07.04.23  H.Hasumi: Formulation changed
!     '07.07.30  H.Hasumi: Loop boundary for S?? calculation reduced
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
!     '12.06.29  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim,  only :  nxydim

  implicit none

  private
  public  ::  modgxy,  shalow

  real(8),     save  ::    fux(nxydim),     fuy(nxydim)
  real(8),     save  ::    fvx(nxydim),     fvy(nxydim)
  real(8),     save  ::    fhx(nxydim),     fhy(nxydim)
  real(8),     save  ::     gu(nxydim),      gv(nxydim)
  real(8),     save  ::     gh(nxydim)

  real(8),     save  ::   accb,   acc
  data accb, acc / -1.d0, 1.d0 /

  logical,     save  ::  ofirst
  character(len=64)  ::  chead(1:16)
  data ofirst / .true. /

contains

  subroutine modgxy(                                                  &
         &      gxx,    gyy,                                          &
         &      ubtx,   vbtx  )

    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nzdim,                                   &
         &   oinit, ofinal
    use zocphy,  only :  gravit
    use zocgrd,  only :  rdepv

    implicit none
    
    real(8),   intent(inout)  ::    gxx(nxydim),   gyy(nxydim)
    real(8),   intent(in)     ::   ubtx(nxydim),  vbtx(nxydim)

!----- local variables
    integer(4)  ::     ij,      i,      j
    integer(4)  ::   ijle,   ijln,   ijlw,   ijls
    integer(4)  ::   ijnw,   ijse,   ijsw
    integer(4)  ::  ifpar,  jfpar,   istat

    namelist /nmaccb/ accb
    namelist /nmaccv/ acc

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then

       ofirst = .false.
       call rewnml( ifpar, jfpar )
       read( ifpar, nmaccb, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmaccb', istat )
       write( jfpar, nmaccb )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmaccv, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmaccv', istat )
       write( jfpar, nmaccv )
       
       if ( accb <= 0.d0 ) then
          accb = acc
       end if
       
       do ij = 1, nxydim
          gh(ij) = gravit / rdepv(ij)
       end do
       
    end if

  end subroutine modgxy

! =====================================================================
  subroutine shalow(                                                  &
         &       hx,   ubtx,   vbtx,                                  &
         &       hy,   ubty,   vbty,                                  &
         &      gxx,    gyy,   ptop,  fw  )
    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nzdim,                                   &
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijtstr, ijtend, ijvstr, ijvend,           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lne,    lse,    lsw
    use zocphy,   only :  rhoo
    use zocgrd,   only :                                              &
         &     cor,  rdepv,    tss,                                   &  
         &      rx,     ry,    rym,                                   & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu,    rxt,    ryt
    use zocmsk,  only :  amskt,  amskv
 
    implicit none
 
    real(8),   intent(inout)  ::    hx(nxydim),  ubtx(nxydim),  vbtx(nxydim)
    real(8),   intent(in)     ::    hy(nxydim),  ubty(nxydim),  vbty(nxydim)
    real(8),   intent(in)     ::   gxx(nxydim),   gyy(nxydim)
    real(8),   intent(in)     ::  ptop(nxydim),    fw(nxydim)
 
!----- local variables
    real(8)     ::     cf
    integer(4)  ::     ij,      i,      j
    integer(4)  ::   ijle,   ijln,   ijlw,   ijls
    integer(4)  ::   ijnw,   ijse,   ijsw
    integer(4)  ::  ifpar,  jfpar,   istat
     
    do ij = 1, nxydim
       fux(ij) = 0.d0
       fuy(ij) = 0.d0
       fvx(ij) = 0.d0
       fvy(ij) = 0.d0
       fhx(ij) = 0.d0
       fhy(ij) = 0.d0
    end do
 
    do ij = ijtstr-nxdim-1, ijtend+nxdim+2
       fhx(ij) = - (  ubty(ij+lw)  * hyu(ij+lw)                       &
    &               + ubty(ij+lsw) * hyu(ij+lsw)) * 0.5d0
    end do
    do ij = ijtstr-nxdim-1, ijtend+nxdim+nxdim+1
       fhy(ij) = - (  vbty(ij+ls) * hxu(ij+ls)                        &
    &               + vbty(ij+lsw) * hxu(ij+lsw)) * 0.5d0
    end do
 
    do ij = ijtstr-nxdim-1, ijtend+nxdim+1
       hx(ij) = hx(ij)                                                &
    &         + tss * (  (fhx(ij+le) - fhx(ij)) * rx                  &
    &                  + (fhy(ij+ln) - fhy(ij)) * ry(ij)) *           &
    &           rxt(ij) * ryt(ij) * amskt(ij, kstr)                   &
    &         - tss * fw(ij) * amskt(ij, kstr)
    end do
 
    do ij = ijstr-nxdim-1, ijend
       gu(ij) = gxx(ij) + cor(ij) * vbtx(ij)                          &
    &         + (  (fux(ij+le) - fux(ij)) * rx * ryu(ij)              &
    &            + (fuy(ij+ln) - fuy(ij)) * rym(ij) * rxu(ij)) *      &
    &           rxu(ij) * ryu(ij)                                     &
    &         - (  (  hy(ij+lne) + hy(ij+le)                          &
    &               - hy(ij+ln ) - hy(ij   )) * gh(ij)                &
    &            + (  ptop(ij+lne) + ptop(ij+le)                      &
    &               - ptop(ij+ln ) - ptop(ij   ))                     &
    &              / rdepv(ij) / rhoo                                 &
    &           ) * 0.5d0 * rx * rxu(ij)                              
       gv(ij) = gyy(ij) - cor(ij) * ubtx(ij)                          &
    &         + (  (fvx(ij+le) - fvx(ij)) * rx * ryu(ij)              &
    &            + (fvy(ij+ln) - fvy(ij)) * rym(ij) * rxu(ij)) *      &
    &           rxu(ij) * ryu(ij)                                     &
    &         - (  (  hy(ij+lne) + hy(ij+ln)                          &
    &               - hy(ij+le ) - hy(ij   )) * gh(ij)                & 
    &            + (  ptop(ij+lne) + ptop(ij+ln)                      &
    &               - ptop(ij+le ) - ptop(ij   ))                     &
    &              / rdepv(ij) / rhoo                                 &
    &           ) * 0.5d0 * rym(ij) * ryu(ij)                         
    end do                                                            
                                                                      
    do ij = ijstr-nxdim-1, ijend                                      
       cf = cor(ij) * tss / accb * 0.5d0                              
       ubtx(ij) = ubtx(ij)                                            &
    &           + tss / accb / (1.d0 + cf * cf) *                     &
    &             (gu(ij) + cf * gv(ij)) * amskv(ij, kstr)            
       vbtx(ij) = vbtx(ij)                                            &
    &           + tss / accb / (1.d0 + cf * cf) *                     &
    &             (gv(ij) - cf * gu(ij)) * amskv(ij, kstr)
    end do

  end subroutine shalow

end module fshlw


