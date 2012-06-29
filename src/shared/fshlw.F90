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
  real(8),     save  ::    sxx(nxydim),     syy(nxydim)
  real(8),     save  ::    sxy(nxydim),     syx(nxydim)
  real(8),     save  ::     gh(nxydim)

  real(8),     save  ::    amh
  data amh / 0.d0 /

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
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijvstr, ijvend,                           &
         &   igstr,  jgstr,                                           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lne,    lse,    lsw,                           &
         &   oinit,  ofinal
    use zocphy,   only :  gravit
    use zocgrd,   only :                                              &
         &   rdepv,                                                   &  
         &      dy,    hxt,    rea,                                   &  
         &      rx,     ry,    rym,                                   & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
    use zocmsk,  only :  amskv,  amfvx,  amfvy

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
    namelist /nmvish/ amh

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
       call rewnml( ifpar, jfpar )
       read( ifpar, nmvish, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmvish', istat )
       write( jfpar, nmvish )
       
       if ( accb <= 0.d0 ) then
          accb = acc
       end if
       
       do ij = 1, nxydim
          gh(ij) = gravit / rdepv(ij)
       end do
       
    end if

    do ij = 1, nxydim
       fux(ij) = 0.d0
       fuy(ij) = 0.d0
       fvx(ij) = 0.d0
       fvy(ij) = 0.d0
    end do
       
    do ij = ijstr-nxdim-1, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijln = ij + ln
       ijle = ij + le
       ijnw = ij + lnw
       ijse = ij + lse
       ijsw = ij + lsw
       sxx(ij) = (ubtx(ij) - ubtx(ijlw)) * amh *                      &
    &             rx / (hxu(ij) + hxu(ijlw)) * 2.d0                   &
    &          + (vbtx(ij) + vbtx(ijlw)) * amh *                      &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0                     &
    &          - (  vbtx(ijln) + vbtx(ijnw)                           & 
    &             - vbtx(ijls) - vbtx(ijsw)) * amh                    &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0
       sxy(ij) = (ubtx(ij) - ubtx(ijls)) * amh *                      &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          - (ubtx(ij) + ubtx(ijls)) * amh *                      &
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     &
    &          + (  vbtx(ijle) + vbtx(ijse)                           &
    &             - vbtx(ijlw) - vbtx(ijsw)) * amh                    &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amh *                      &
    &           (hyxu(ij) + hyxu(ijls)) * 0.25d0
       syy(ij) = (vbtx(ij) - vbtx(ijls)) * amh *                      &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          + (ubtx(ij) + ubtx(ijls)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     &
    &          - (  ubtx(ijle) + ubtx(ijse)                           &
    &             - ubtx(ijlw) - ubtx(ijsw)) * amh                    &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amh *                      & 
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0
       syx(ij) = (vbtx(ij) - vbtx(ijlw)) * amh *                      &
    &            rx / (hxu(ij) + hxu(ijlw)) * 2.d0                    &
    &          - (vbtx(ij) + vbtx(ijlw)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     &
    &          + (  ubtx(ijln) + ubtx(ijnw)                           &
    &             - ubtx(ijls) - ubtx(ijsw)) * amh                    &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amh *                      &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0
    end do

    do ij = ijvstr-nxdim-1, ijvend+1
       ijlw = ij + lw
       fux(ij) = sxx(ij) * (hyu(ij) + hyu(ijlw)) *                   &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
       fvx(ij) = syx(ij) * (hyu(ij) + hyu(ijlw)) *                   &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
    end do
    do ij = ijvstr-nxdim-1, ijvend+nxdim
       ijls = ij + ls
       fuy(ij) = sxy(ij) * (hxu(ij) + hxu(ijls)) *                   &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
       fvy(ij) = syy(ij) * (hxu(ij) + hxu(ijls)) *                   &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
    end do

    do ij = ijvstr-nxdim-1, ijvend+nxdim+1
       gxx(ij) = gxx(ij)                                             &
    &          - (  (fux(ij+le) - fux(ij)) * rx * ryu(ij)            & 
    &             + (fuy(ij+ln) - fuy(ij)) * rym(ij) * rxu(ij)) *    &
    &            rxu(ij) * ryu(ij) * amskv(ij, kstr)
       gyy(ij) = gyy(ij)                                             &
    &          - (  (fvx(ij+le) - fvx(ij)) * rx * ryu(ij)            &
    &             + (fvy(ij+ln) - fvy(ij)) * rym(ij) * rxu(ij)) *    &
    &            rxu(ij) * ryu(ij) * amskv(ij, kstr)
    end do

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
         &      dy,    hxt,    rea,                                   &  
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
       fhx(ij) =  - (  ubty(ij+lw) * hyu(ij+lw)                       &
    &                + ubty(ij+lsw) * hyu(ij+lsw)) * 0.5d0
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
 
    do ij = ijstr-nxdim-1, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijln = ij + ln
       ijle = ij + le
       ijnw = ij + lnw
       ijse = ij + lse
       ijsw = ij + lsw
       sxx(ij) = (ubtx(ij) - ubtx(ijlw)) * amh *                      &
    &             rx / (hxu(ij) + hxu(ijlw)) * 2.d0                   &
    &          + (vbtx(ij) + vbtx(ijlw)) * amh *                      &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0                     &
    &          - (  vbtx(ijln) + vbtx(ijnw)                           &
    &             - vbtx(ijls) - vbtx(ijsw)) * amh                    &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     
       sxy(ij) = (ubtx(ij) - ubtx(ijls)) * amh *                      &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          - (ubtx(ij) + ubtx(ijls)) * amh *                      & 
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     &
    &          + (  vbtx(ijle) + vbtx(ijse)                           &
    &             - vbtx(ijlw) - vbtx(ijsw)) * amh                    &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     
       syy(ij) = (vbtx(ij) - vbtx(ijls)) * amh *                      &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          + (ubtx(ij) + ubtx(ijls)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     &
    &          - (  ubtx(ijle) + ubtx(ijse)                           &
    &             - ubtx(ijlw) - ubtx(ijsw)) * amh                    &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amh *                      &
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     
       syx(ij) = (vbtx(ij) - vbtx(ijlw)) * amh *                      &
    &            rx / (hxu(ij) + hxu(ijlw)) * 2.d0                    &
    &          - (vbtx(ij) + vbtx(ijlw)) * amh *                      &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     &
    &          + (  ubtx(ijln) + ubtx(ijnw)                           &
    &             - ubtx(ijls) - ubtx(ijsw)) * amh                    &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amh *                      &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0
    end do
 
    do ij = ijvstr-nxdim-1, ijvend+1
       ijlw = ij + lw
       fux(ij) = sxx(ij) * (hyu(ij) + hyu(ijlw)) *                    &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
       fvx(ij) = syx(ij) * (hyu(ij) + hyu(ijlw)) *                    &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
    end do
 
    do ij = ijvstr-nxdim-1, ijvend+nxdim
       ijls = ij + ls
       fuy(ij) = sxy(ij) * (hxu(ij) + hxu(ijls)) *                    & 
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
       fvy(ij) = syy(ij) * (hxu(ij) + hxu(ijls)) *                    &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
    end do
 
    do ij = ijstr-nxdim-1, ijend
       gu(ij) = gxx(ij) + cor(ij) * vbtx(ij)                         &
    &         + (  (fux(ij+le) - fux(ij)) * rx * ryu(ij)             &
    &            + (fuy(ij+ln) - fuy(ij)) * rym(ij) * rxu(ij)) *     &
    &           rxu(ij) * ryu(ij)                                    &
    &         - (  (  hy(ij+lne) + hy(ij+le)                         &
    &               - hy(ij+ln ) - hy(ij   )) * gh(ij)               &
    &            + (  ptop(ij+lne) + ptop(ij+le)                     &
    &               - ptop(ij+ln ) - ptop(ij   ))                    &
    &              / rdepv(ij) / rhoo                                &
    &           ) * 0.5d0 * rx * rxu(ij)
       gv(ij) = gyy(ij) - cor(ij) * ubtx(ij)                         &
    &         + (  (fvx(ij+le) - fvx(ij)) * rx * ryu(ij)             &
    &            + (fvy(ij+ln) - fvy(ij)) * rym(ij) * rxu(ij)) *     &
    &           rxu(ij) * ryu(ij)                                    &
    &         - (  (  hy(ij+lne) + hy(ij+ln)                         &
    &               - hy(ij+le ) - hy(ij   )) * gh(ij)               & 
    &            + (  ptop(ij+lne) + ptop(ij+ln)                     &
    &               - ptop(ij+le ) - ptop(ij   ))                    &
    &              / rdepv(ij) / rhoo                                &
    &           ) * 0.5d0 * rym(ij) * ryu(ij)
    end do
 
    do ij = ijstr-nxdim-1, ijend
       cf = cor(ij) * tss / accb * 0.5d0
       ubtx(ij) = ubtx(ij)                                           &
    &           + tss / accb / (1.d0 + cf * cf) *                    &
    &             (gu(ij) + cf * gv(ij)) * amskv(ij, kstr)
       vbtx(ij) = vbtx(ij)                                           &
    &           + tss / accb / (1.d0 + cf * cf) *                    &
    &             (gv(ij) - cf * gu(ij)) * amskv(ij, kstr)
    end do

  end subroutine shalow

end module fshlw


