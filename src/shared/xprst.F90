module xprst

! --- information -----------------------------------------------------
!
!  Coefficients for the approximated equation of state by McDougall et
! al. (2003, JAOT).
!
!  HISTORY
!     '07.05.01  H.Hasumi
!     '12.06.15  H.Tatebe: for COCO5.0 in F90
!
! ---------------------------------------------------------------------

  implicit none
  private

  real(8),     parameter  ::  P10    =  9.99843699D+2
  real(8),     parameter  ::  P1T    =  7.35212840D+0
  real(8),     parameter  ::  P1TT   = -5.45928211D-2
  real(8),     parameter  ::  P1TTT  =  3.98476704D-4
  real(8),     parameter  ::  P1S    =  2.96938239D+0
  real(8),     parameter  ::  P1ST   = -7.23268813D-3
  real(8),     parameter  ::  P1SS   =  2.12382341D-3
  real(8),     parameter  ::  P1P    =  1.04004591D-2
  real(8),     parameter  ::  P1PTT  =  1.03970529D-7
  real(8),     parameter  ::  P1PS   =  5.18761880D-6
  real(8),     parameter  ::  P1PP   = -3.24041825D-8
  real(8),     parameter  ::  P1PPTT = -1.23869360D-11

  real(8),     parameter  ::  P20     =  1.D0
  real(8),     parameter  ::  P2T     =  7.28606739D-3
  real(8),     parameter  ::  P2TT    = -4.60835542D-5
  real(8),     parameter  ::  P2TTT   =  3.68390573D-7
  real(8),     parameter  ::  P2TTTT  =  1.80809186D-10
  real(8),     parameter  ::  P2S     =  2.14691708D-3
  real(8),     parameter  ::  P2ST    = -9.27062484D-6
  real(8),     parameter  ::  P2STTT  = -1.78343643D-10
  real(8),     parameter  ::  P2SS    =  4.76534122D-6
  real(8),     parameter  ::  P2SSTT  =  1.63410736D-9
  real(8),     parameter  ::  P2P     =  5.30848875D-6
  real(8),     parameter  ::  P2PPTTT = -3.03175128D-16
  real(8),     parameter  ::  P2PPPT  = -1.27934137D-17

  public  ::  secoef,  secfbb  !  used in tflxt.F90, dvdif.F90, tovtr.F90, ctnuv.F

contains


  subroutine secoef(                                                  &
         &     c0,     c1,     c2,     c3,                            &
         &     c4,     c5,     c6,                                    &
         &     d0,     d1,     d2,     d3,     d4,                    &
         &     d5,     d6,     d7,     d8,     d9   )

    use zocdim,  only  :   kstr,     nz
    use zocgrd,  only  :    dz0

    implicit none

    real(8),  intent(out)  ::    c0(nz), c1(nz), c2(nz), c3(nz), c4(nz)
    real(8),  intent(out)  ::    c5(nz), c6(nz)
    real(8),  intent(out)  ::    d0(nz), d1(nz), d2(nz), d3(nz), d4(nz)
    real(8),  intent(out)  ::    d5(nz), d6(nz), d7(nz), d8(nz), d9(nz)

    real(8)                ::    z(nz)
    integer(4)             ::  k

    z(1) = 0.5d0 * dz0(kstr) / 1.0d2
    do k = 2, nz
       z(k) = z(k-1) + 0.5d0 * (dz0(k+kstr-1) + dz0(k+kstr-2)) / 1.d2
    end do
    
    do k = 1, nz
       c0(k) = p10 + (p1p + p1pp * z(k)) * z(k)
       c1(k) = p1t
       c2(k) = p1tt + (p1ptt + p1pptt * z(k)) * z(k)
       c3(k) = p1ttt
       c4(k) = p1s + p1ps * z(k)
       c5(k) = p1st
       c6(k) = p1ss

       d0(k) = p20 + p2p * z(k)
       d1(k) = p2t + p2pppt * z(k) * z(k) * z(k)
       d2(k) = p2tt
       d3(k) = p2ttt + p2ppttt * z(k) * z(k)
       d4(k) = p2tttt
       d5(k) = p2s
       d6(k) = p2st
       d7(k) = p2sttt
       d8(k) = p2ss
       d9(k) = p2sstt
    end do

  end subroutine secoef

! **********************************************************************

  subroutine secfbb(                                                  &
         &     c0,     c1,     c2,     c3,                            &
         &     c4,     c5,     c6,                                    &
         &     d0,     d1,     d2,     d3,     d4,                    &
         &     d5,     d6,     d7,     d8,     d9   )

    use zocdim,  only  :   kstr,     nz,  nxydim,  ijvstr,  ijvend,   &  
                             le,     ln,     lne
    use zocgrd,  only  :   dept,    dz0

    implicit none

    real(8),  intent(out)    ::    c0(nxydim), c1(nxydim), c2(nxydim)
    real(8),  intent(out)    ::    c3(nxydim), c4(nxydim), c5(nxydim)
    real(8),  intent(out)    ::    c6(nxydim)
    real(8),  intent(out)    ::    d0(nxydim), d1(nxydim), d2(nxydim)
    real(8),  intent(out)    ::    d3(nxydim), d4(nxydim), d5(nxydim)
    real(8),  intent(out)    ::    d6(nxydim), d7(nxydim), d8(nxydim)
    real(8),  intent(out)    ::    d9(nxydim)

    real(8)                  ::     z
    integer(4)               ::     k,         ij

    do ij = ijvstr, ijvend

       z = 0.25d0 * 1.d-2 * ( dept(ij)    + dept(ij+le)               &    
    &                       + dept(ij+ln) + dept(ij+lne) )

       c0(ij) = p10 + (p1p + p1pp * z) * z
       c1(ij) = p1t
       c2(ij) = p1tt + (p1ptt + p1pptt * z) * z
       c3(ij) = p1ttt
       c4(ij) = p1s + p1ps * z
       c5(ij) = p1st
       c6(ij) = p1ss
       
       d0(ij) = p20 + p2p * z
       d1(ij) = p2t + p2pppt * z * z * z
       d2(ij) = p2tt
       d3(ij) = p2ttt + p2ppttt * z * z
       d4(ij) = p2tttt
       d5(ij) = p2s
       d6(ij) = p2st
       d7(ij) = p2sttt
       d8(ij) = p2ss
       d9(ij) = p2sstt

    end do
    
  end subroutine secfbb

end module xprst


