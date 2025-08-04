module sfcng
  
  implicit none
  private
   
  public  ::  sfcflx,   ocnvar
#ifdef OPT_BODY
  public  ::  bdyflx
#endif

contains
  
  subroutine sfcflx(    ft,     t )

    use zocdim
    use zocgrd
    use zocmsk
    use zocphy
    
    implicit none
    
    real(8),    intent(in)   ::  t(nxydim, nzdim, ntdim)
    real(8),    intent(out)  :: ft(nxydim, ntdim)

!---- local
    integer(4)  ::    ij,     l      
    
! only for passive tracer fluxes
    do l = 3, ntdim
       do ij = ijstr, ijend
          ft(ij, l) = 0.d0
       end do
    end do
    
  end subroutine sfcflx

! ---------------------------------------------------------------------  
  subroutine ocnvar(     t,     u,     v,    hb,   ahv )
    
    use zocdim
    use zocgrd
    use utint
    
    implicit none

    real(8),    intent(out)  ::    t(nxydim, nzdim, ntdim)
    real(8),    intent(out)  ::    u(nxydim, nzdim)
    real(8),    intent(out)  ::    v(nxydim, nzdim)
    real(8),    intent(out)  ::   hb(nxydim)
    real(8),    intent(out)  ::  ahv(nxydim, nzdim)
      
    call tmintb( t(1, 1, 1), 1 )
    call tmintb( t(1, 1, 2), 2 )
    call tmintb(   u, 3 )
    call tmintb(   v, 4 )
    call tmintb( ahv, 5 )
      
    tt = tt + ts
    call tmintp(hb, 1) ! one-time-step advanced
    tt = tt - ts
    
  end subroutine ocnvar

#ifdef OPT_BODY
! ---------------------------------------------------------------------  
  subroutine bdyflx(    tq,     t )
     
    use zocdim
    use zocgrd
    use zocmsk
     
    implicit none
     
    real(8),    intent(in)    ::   t(nxdim,  nzdim,  ntdim)
    real(8),    intent(inout) ::  tq(nxdim,  nzdim,  ntdim)      
     
  end subroutine bdyflx
#endif
   
end module sfcng

 
