#ifdef OPT_TOUZA
subroutine ucaln_dummy
return
end subroutine ucaln_dummy
#else
module ucaln

! --- information -----------------------------------------------------
!
!  Conversion between calendar date and total second
!
!  HISTORY
!     '97.03.21  H.Hasumi: from AGCM5.4 developed by A.Numaguti
!     '12.11.27  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  implicit none

  private

  integer(4),    save  ::   monday(1:12,1:2)
  integer(4),    save  ::   idaymo,  imonyr,  isecmn,  iminhr,  ihrday
  logical,       save  ::   ogrego,  oideal

  data  monday / 31,28,31,30,31,30,31,31,30,31,30,31,                 &
                 31,29,31,30,31,30,31,31,30,31,30,31 /
  data  idaymo, imonyr / 30, 12 / 
  data  isecmn, iminhr, ihrday / 60, 60, 24 /
  data  ogrego, oideal / .true., .false. /

  public  ::  calndr,  cyh2ss,  css2yh

contains

  subroutine calndr(  icaln )

    implicit none
    
    integer(4),         intent(in)     ::  icaln

    if ( icaln == 1 ) then
       oideal = .false.
       ogrego = .false.
    else if ( icaln == 2 ) then
       oideal = .false.
       ogrego = .true.
    else
       oideal = .true.
       ogrego = .false.
    end if
    
  end subroutine calndr

!==================================================================
  subroutine cdaymo( ndaymo, iyear, imonth )

    implicit none

    integer(4),   intent(inout)  ::  ndaymo
    integer(4),   intent(in)     ::   iyear,  imonth
     
    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndaymo = monday( imonth,2 )
       else
          ndaymo = monday( imonth,1 )            
       end if
    else if ( .not. oideal ) then
       ndaymo = monday( imonth,1 )            
    else
       ndaymo = idaymo
    end if
    
  end subroutine cdaymo

!==================================================================
  subroutine cdayyr( ndayyr, iyear )

    implicit none

    integer(4),   intent(inout)  ::  ndayyr
    integer(4),   intent(in)     ::   iyear

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndayyr = 366
       else
          ndayyr = 365
       end if
    else if ( .not. oideal ) then
       ndayyr = 365
    else
       ndayyr = idaymo*imonyr
    end if
    
  end subroutine cdayyr
!==================================================================
  subroutine cmonyr( nmonyr, iyear )

    implicit none

    integer(4),   intent(inout)  ::  nmonyr
    integer(4),   intent(in)     ::   iyear

    nmonyr = imonyr

  end subroutine cmonyr
!==================================================================
  subroutine csecdy( nsecdy )

    implicit none

    integer(4),   intent(inout)  ::  nsecdy

    nsecdy = isecmn*iminhr*ihrday

  end subroutine csecdy
!==================================================================
  subroutine csecmi( nsecmi )

    implicit none

    integer(4),   intent(inout)  ::  nsecmi

    nsecmi = isecmn

  end subroutine csecmi
!==================================================================
  subroutine csechr( nsechr )

    implicit none
    integer(4),   intent(inout)  ::  nsechr

    nsechr = isecmn*iminhr
    
  end subroutine csechr
!==================================================================
  subroutine  css2ds( idays, rsec, dsec )

    implicit none

    integer(4),    intent(inout) ::  idays
    real(8),       intent(inout) ::  rsec
    real(8),       intent(in)    ::  dsec

    integer(4)  ::  isecdy  

    isecdy = isecmn*iminhr*ihrday
    idays  = int( dsec/dble (isecdy) ) + 1
    rsec   = dsec - dble (idays-1)*dble (isecdy)
    if ( nint( rsec ) >= isecdy ) then
       idays = idays + 1
       rsec  = rsec - dble(isecdy)
    end if

  end subroutine css2ds
!==================================================================
  subroutine cds2ss( dsec, idays, rsec )

    implicit none

    real(8),       intent(inout) ::  dsec
    integer(4),    intent(in)    ::  idays
    real(8),       intent(in)    ::  rsec

    integer(4)  ::  isecdy

    isecdy = isecmn*iminhr*ihrday
    dsec   = dble (idays-1)*dble (isecdy) + dble (rsec)

  end subroutine CDS2SS
!==================================================================
  subroutine crs2hm( ihour, imin, isec, rsec )

    implicit none

    integer(4),    intent(inout) ::  ihour,  imin,   isec
    real(8),       intent(in)    ::  rsec

    integer(4)  ::  isechr

    isechr = isecmn*iminhr
    ihour  = int ( rsec / dble(isechr ) )
    imin   = int ( ( rsec - dble(ihour*isechr) )/dble(isecmn) )
    isec   = nint( rsec - dble(ihour*isechr) - dble(imin*isecmn) )
    if ( isec >= isecmn ) then
       imin  = imin + 1
       isec  = isec - isecmn
    end if
    if ( imin .eq. iminhr ) then
       ihour = ihour + 1
       imin  = imin  - iminhr
    end if

  end subroutine crs2hm
!==================================================================
  subroutine chm2rs( rsec, ihour, imin, isec )

    implicit none

    real(8),       intent(inout) ::  rsec
    integer(4),    intent(in)    ::  ihour,  imin,   isec

    rsec = ihour*isecmn*iminhr + imin*isecmn + isec

  end subroutine chm2rs
!==================================================================
  subroutine cdd2ym( iyear, imonth, iday, idays )

    implicit none

    integer(4),    intent(inout) ::  iyear,  imonth, iday
    integer(4),    intent(in)    ::  idays

!---- internal work
    integer(4)   ::    jy4,  jcent,  jcent4,   jy
    integer(4)   :: idays0,    idy,  idayyr
    integer(4)   ::  ileap,     id
    integer(4)   ::      m

    if ( ogrego ) then
       jy = int(dble(idays)/365.24)
1100   continue 
       jy4    = (jy+3)/4
       jcent  = (jy+99)/100
       jcent4 = (jy+399)/400
       idays0 = jy*365+jy4-jcent+jcent4
       if ( idays <= idays0 ) then
          jy = jy -1 
          if ( jy >= 0 ) goto 1100
       end if
       iyear = jy
       idy   = idays - idays0
       if ( ocleap( iyear ) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    else if ( .not. oideal ) then
!       iyear = (idays-1)/365
       iyear = floor(dble(idays-1)/365.d0)
       idy   = idays - iyear*365
       ileap = 1
    end if

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, 12
          id = id + monday(m,ileap)
          if ( idy <= id ) then
             imonth = m
             iday   = idy-id+monday(m,ileap)
             goto 3190
          end if
       end do
3190   continue 
    else 
       idayyr = idaymo*imonyr
       iyear  = ( idays-1 ) / idayyr
       imonth = ( idays-1 - iyear*idayyr )/idaymo+1
       iday   = idays - iyear*idayyr - (imonth-1)*idaymo
    end if

  end subroutine cdd2ym
!======================================================================
  subroutine cym2dd( idays, iyear , imonth, iday )

    implicit none

    integer(4),    intent(inout) ::  idays
    integer(4),    intent(in)    ::  iyear,  imonth, iday

!---- local variables
    integer(4)   ::    jy4,  jcent,  jcent4
    integer(4)   ::  jyear, jmonth
    integer(4)   :: idays0
    integer(4)   ::  ileap,     id
    integer(4)   ::      m

    if ( ogrego .or. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       end if
    end if
    
    if ( ogrego ) then
       jy4    = (jyear+3)/4
       jcent  = (jyear+99)/100
       jcent4 = (jyear+399)/400
       idays0 = jyear*365+jy4-jcent+jcent4
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       end if
    else if ( .not. oideal ) then
       IDAYS0 = JYEAR*365
       ileap  = 1
    end if

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       end do
    else
       idays0 = iyear*idaymo*imonyr
       id     = (imonth-1)*idaymo
    end if

    idays = idays0 + id + iday

  end subroutine cym2dd
!======================================================================
  subroutine cym2yd( idaysy, iyear, imonth, iday )

    implicit none

    integer(4),    intent(inout) ::  idaysy
    integer(4),    intent(in)    ::  iyear,  imonth, iday

!---- local variables
    integer(4)   ::  jyear, jmonth
    integer(4)   ::  ileap,     id
    integer(4)   ::      m

    if ( ogrego .or. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       end if
    end if

    if ( ogrego ) then
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       end if
    else if ( .not. oideal ) then
       ileap  = 1
    end if

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       end do
    else
       id     = (imonth-1)*idaymo
    end if

    idaysy = id + iday

  end subroutine cym2yd

!*********************************************************************
  subroutine css2yh( idate, dsec )

    implicit none

    integer(4),    intent(inout) ::  idate(1:6)
    real(8),       intent(in)    ::  dsec

!---- internal work
    real(8)      ::   rsec
    integer(4)   ::  idays

    call css2ds( idays, rsec, dsec )
    call cdd2ym( idate(1), idate(2), idate(3), idays )
    call crs2hm( idate(4), idate(5), idate(6), rsec  )

  end subroutine css2yh

!===============================================================
  subroutine cyh2ss( dsec, idate )

    implicit none

    real(8),       intent(inout) ::  dsec
    integer(4),    intent(in)    ::  idate(1:6)

!---- internal work
    real(8)      ::   rsec
    integer(4)   ::  idays

    call cym2dd( idays, idate(1), idate(2), idate(3) )
    call chm2rs(  rsec, idate(4), idate(5), idate(6) )
    call cds2ss(  dsec, idays, rsec )

  end subroutine cyh2ss

!*********************************************************************
  function ocleap( iyear )

    implicit none

    integer(4),    intent(in)    ::  iyear
    logical    ocleap

!---- internal work
    integer(4)   ::  iy, iycen, icent

    iy     = mod( iyear,     4  )
    iycen  = mod( iyear,    100 )
    icent  = mod( iyear/100, 4  )

    if ( iy == 0 .and. ( iycen /= 0 .or. icent == 0 ) ) then
       ocleap = .true.
    else
       ocleap = .false.
    end if

  end function ocleap

end module ucaln
#endif

