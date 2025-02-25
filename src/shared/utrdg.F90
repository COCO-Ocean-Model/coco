module utrdg

  implicit none
  private

  public  ::  thmasc  !  used in 
  public  ::  thomas  !  used in 

contains

  subroutine thmasc(                                     &
         &      dr,     di,                              &
         &       a,     br,     bi,      c      )
  
  ! --- information -----------------------------------------------------
  !
  !  Solve the linear equations expressed by tri-diagonal matrix
  !
  !  HISTORY
  !     '99.04.12  H.Hasumi: from CCSR2
  !     '01.02.02  H.Hasumi: solve T and S at once
  !     '02.05.29  H.Nakano: tracer dimension
  !     '07.04.23  H.Hasumi
  !     '12.06.15  H.Tatebe: for COCO5.0 in F90
  !
  ! ---------------------------------------------------------------------
  
    use zocdim,  only :                                  &
         &  nxydim,  nzdim,                              &
         &    kstr,   kend,                              &
         &  ijvstr,  ijvend
  
    implicit none
  
    real(8),   intent(inout)  ::    dr(nxydim, nzdim),   di(nxydim, nzdim)
    real(8),   intent(inout)  ::    br(nxydim, nzdim),   bi(nxydim, nzdim)
    real(8),   intent(in)     ::     a(nxydim, nzdim),    c(nxydim, nzdim)
  
    real(8)                   ::    r1,    pr1,    pi1,    qr1,    qi1
    real(8)                   ::    rr,     ri,      r,     sr,     si
    integer(4)                ::    ij,      k

#ifdef ADF_
!$acc data copyin(a, c) 
!$acc data copy(dr, di, br, bi) 
#endif

#ifdef OMP_
!$omp parallel do private(ij, r1, pr1, pi1, qr1, qi1)
#elif  ACC_
!$acc data copyin(c) 
!$acc data copy(bi, br, di, dr)  
!$acc kernels 
!$acc loop independent
#endif
    do ij = ijvstr, ijvend
       r1 = br(ij, kstr) * br(ij, kstr)                                 &
    &     + bi(ij, kstr) * bi(ij, kstr)
       pr1 =   c(ij, kstr) * br(ij, kstr) / r1
       pi1 = - c(ij, kstr) * bi(ij, kstr) / r1
       qr1 = (  dr(ij, kstr) * br(ij, kstr)                             &
    &         + di(ij, kstr) * bi(ij, kstr)) / r1
       qi1 = (  di(ij, kstr) * br(ij, kstr)                             &
    &        - dr(ij, kstr) * bi(ij, kstr)) / r1
       br(ij, kstr) = pr1
       bi(ij, kstr) = pi1
       dr(ij, kstr) = qr1
       di(ij, kstr) = qi1
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(c)
!$acc end data ! copy(bi, br, di, dr)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copyin(c, a) 
!$acc data copy(bi, br, di, dr)  
!$acc kernels 
#endif
    do k = kstr+1, kend
#ifdef OMP_
!$omp do private(ij, rr, ri, sr, si, r)
#endif
       do ij = ijvstr, ijvend
          rr = br(ij, k) - a(ij, k) * br(ij, k-1)
          ri = bi(ij, k) - a(ij, k) * bi(ij, k-1)
          sr = dr(ij, k) - a(ij, k) * dr(ij, k-1)
          si = di(ij, k) - a(ij, k) * di(ij, k-1)
          r = rr * rr + ri * ri
          br(ij, k) =   c(ij, k) * rr / r
          bi(ij, k) = - c(ij, k) * ri / r
          dr(ij, k) = (sr * rr + si * ri) / r
          di(ij, k) = (si * rr - sr * ri) / r
       end do
#ifdef OMP_
!$omp end do
#endif
    end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(c, a)
!$acc end data ! copy(bi, br, di, dr)
#endif
 
#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copy(dr, di)   
!$acc data copyin(bi, br)  
!$acc kernels 
#endif
    do k = kend-1, kstr, -1
#ifdef OMP_
!$omp do private(ij)
#endif
       do ij = ijvstr, ijvend
          dr(ij, k) = dr(ij, k) - br(ij, k) * dr(ij, k+1)               &
    &                           + bi(ij, k) * di(ij, k+1)
          di(ij, k) = di(ij, k) - bi(ij, k) * dr(ij, k+1)               &
    &                           - br(ij, k) * di(ij, k+1)
       end do
#ifdef OMP_
!$omp end do
#endif
    end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(dr, di)
!$acc end data ! copyin(bi, br) 
#endif
 
#ifdef ADF_
!$acc end data ! copyin(a, c)
!$acc end data ! copy(dr, di, br, bi)
#endif

  end subroutine thmasc
  
  ! *********************************************************************
  
  subroutine thomas(                                     &
         &     adt,     ac,                              &
         &      aa,     ab    )
  
    use zocdim,  only :                                  &
         &  nxydim,  nzdim,  ntdim,                      &
         &    kstr,   kend,                              &
         &  ijtstr,  ijtend
  
    implicit none
  
    real(8),   intent(inout)  ::   adt(nxydim, nzdim, ntdim)
    real(8),   intent(inout)  ::    ac(nxydim, nzdim)
    real(8),   intent(in)     ::    aa(nxydim, nzdim),   ab(nxydim, nzdim)
  
    real(8)                   ::    fc
    integer(4)                ::    ij,       k,       n
  
#ifdef ADF_
!$acc data copy(adt, ac) 
!$acc data copyin(aa, ab) 
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(ab) 
!$acc data copy(ac) 
!$acc kernels 
#endif
    do ij = ijtstr, ijtend
       ac(ij, kstr) = ac(ij, kstr) / ab(ij, kstr)
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(ab)
!$acc end data ! copy(ac)
#endif
#ifdef OMP_
!$omp parallel do collapse(2) private(n, ij)
#elif  ACC_
!$acc data copyin(ab) 
!$acc data copy(adt) 
!$acc kernels
#endif
    do n = 1, ntdim
!$omp do
       do ij = ijtstr, ijtend
          adt(ij, kstr, n) = adt(ij, kstr, n) / ab(ij, kstr)
       end do
!$omp end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(ab)
!$acc end data ! copy(adt)
#endif
  
#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copy(adt) 
!$acc data copy(ac) 
!$acc data copyin(aa, ab)
!$acc kernels 
!$acc loop seq
#endif
    do k = kstr+1, kend
#ifdef OMP_
!$omp do private(ij, fc, n)
#elif  ACC_
!$acc loop independent
#endif
       do ij = ijtstr, ijtend
          fc = 1.d0 / ( ab(ij, k) - aa(ij, k) * ac(ij, k-1) )
          ac(ij, k) = ac(ij, k) * fc
#ifdef OMP_
#elif  ACC_
!$acc loop independent
#endif
          do n = 1, ntdim
             adt(ij, k, n) = ( adt(ij, k, n)                            &
    &                      - aa(ij, k) * adt(ij, k-1, n) ) * fc
          end do
       end do
#ifdef OMP_
!$omp end do
#endif
    end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(adt)
!$acc end data ! copy(ac)
!$acc end data ! copyin(aa, ab)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copy(adt) 
!$acc data copyin(ac)
!$acc kernels 
#endif
    do n = 1, ntdim
#ifdef OMP_
#elif  ACC_
!$acc loop seq
#endif
       do k = kend-1, kstr, -1
#ifdef OMP_
!$omp do private(ij)
#endif
          do ij = ijtstr, ijtend
             adt(ij, k, n) = adt(ij, k, n) - ac(ij, k) * adt(ij, k+1, n)
          end do
#ifdef OMP_
!$omp end do
#endif
       end do
    end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(adt)
!$acc end data ! copyin(ac)
#endif
  
#ifdef ADF_
!$acc end data ! copy(adt, ac)
!$acc end data ! copyin(aa, ab)
#endif
  end subroutine thomas
  
end module utrdg

