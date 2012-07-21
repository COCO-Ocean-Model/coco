      SUBROUTINE FLXTRC(
     O                ADT,  DIFFZ,
     I                 TX,     HX,
     I                 TY,     HZ,
     I                 UY,     VY,
     I                  W,    AHV)
      use zocdim
      use zocgrd
      use zocmsk
      use zocfil
      use bstbc
      IMPLICIT NONE

c --- information -----------------------------------------------------
c  Second order moment scheme for tracer advection of Prather (1988, JGR)
c
c  HISTORY
c     '07.07.10  H.Tatebe
c     '07.10.24  H.Hasumi: for COCO4.3
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.04.02  H.Tatebe: bug fix
c     '09.04.08  H.Tatebe: bug fix 2
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)    
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      HX(NXYDIM),     HZ(NXYDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      REAL*8     WZC(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
c      COMMON /WORK/ FTX, FTY, FTZ, WZC, RZM
      COMMON /WORKFT/ FTX, FTY, FTZ

      REAL*8   HZBOT(NXYDIM)

      REAL*8   XDZDX(NXYDIM, NZDIM),  YDZDY(NXYDIM, NZDIM)
      REAL*8   ZDZDX(NXYDIM, NZDIM),  ZDZDY(NXYDIM, NZDIM)
      REAL*8   XDTDZ(NXYDIM, NZDIM, NTDIM),  YDTDZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   ZDTDX(NXYDIM, NZDIM, NTDIM),  ZDTDY(NXYDIM, NZDIM, NTDIM)

!---- spatially varying isopycnal diffusion coefficient
      REAL*8    AHH3D(NXYDIM, NZDIM),  AHI3D(NXYDIM, NZDIM)
      REAL*8    AHG3D(NXYDIM, NZDIM) 
      SAVE AHI3D, AHG3D, AHH3D

!---- for second order moment
      REAL*8   S0 (NXYDIM, NZDIM, NTDIM)
      REAL*8   SM (NXYDIM, NZDIM, NTDIM)
      REAL*8   SX (NXYDIM, NZDIM, NTDIM), SXX(NXYDIM, NZDIM, NTDIM)
      REAL*8   SY (NXYDIM, NZDIM, NTDIM), SYY(NXYDIM, NZDIM, NTDIM)
      REAL*8   SZ (NXYDIM, NZDIM, NTDIM), SZZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   SXY(NXYDIM, NZDIM, NTDIM), SXZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   SYZ(NXYDIM, NZDIM, NTDIM)   

      REAL*8   F0 (NXYDIM, NZDIM)
      REAL*8   FM (NXYDIM, NZDIM)
      REAL*8   FX (NXYDIM, NZDIM), FXX(NXYDIM, NZDIM)
      REAL*8   FY (NXYDIM, NZDIM), FYY(NXYDIM, NZDIM) 
      REAL*8   FZ (NXYDIM, NZDIM), FZZ(NXYDIM, NZDIM) 
      REAL*8   FXY(NXYDIM, NZDIM), FXZ(NXYDIM, NZDIM)
      REAL*8   FYZ(NXYDIM, NZDIM) 

      REAL*8   VLMX(NXYDIM, NZDIM), VLMY(NXYDIM, NZDIM)
      REAL*8   VLMZ(NXYDIM)
      REAL*8   R(NXYZDM)
      REAL*8   ALF(NXYDIM, NZDIM), UV(NXYDIM, NZDIM)

!---- bug fix
      SAVE S0, SM, SX, SXX, SY, SYY, SZ, SZZ, SXY, SXZ, SYZ
      SAVE F0, FM, FX, FXX, FY, FYY, FZ, FZZ, FXY, FXZ, FYZ
      SAVE VLMX, VLMY, VLMZ, R, ALF, UV

      REAL*8   S0M,    S1M,    S0P,    SXP
      REAL*8   ALFQ,   ALF1,   ALF1Q
      REAL*8   U,      V,      TMP

      REAL*8   EPS,   SQ3,   CI3,   TSIV
      SAVE EPS,   SQ3,   CI3,   TSIV

      INTEGER     IJ,     K,      N,      I,     J
      INTEGER     IJLW,   IJLSW,  IJLE
      INTEGER     IJLN,   IJLS
      INTEGER     KUU,     KU,     KD
      INTEGER   IFPAR,  JFPAR

      LOGICAL OFIRST,   OEOF
      SAVE OFIRST, OEOF

      REAL*8     AHH,    AHI,    AHG
      SAVE       AHH,    AHI,    AHG
      NAMELIST /NMDIFH/ AHH
      NAMELIST /NMDIFI/ AHI
      NAMELIST /NMDIFG/ AHG

!---- 
#ifdef OPT_PARALLEL
      REAL*8   BUF3(NXG, NYG, NZ)
      REAL*8   G3D(NXGDIM, NYGDIM, NZDIM)
      EQUIVALENCE(FTX(1, 1, 1), G3D(1, 1, 1))
#else
      REAL*8   BUF3(NX, NY, NZ)
      EQUIVALENCE(FTX(1, 1, 1), BUF3(1, 1, 1))
#endif

!---- file name of isotropic diffusion and thickness diffusion
      CHARACTER  CFAHI*(NCF), CFAHG*(NCF)
      CHARACTER  CHEAD(64)*16 
      INTEGER  NFAHI, NFAHG, IAH
      NAMELIST /NMCAH/ CFAHI, CFAHG, IAH
      DATA CFAHI, CFAHG / 'not-specified', 'not-specified' /

#ifdef OPT_BBL
      REAL*8  AHHBBL
      SAVE AHHBBL
      NAMELIST /NMBBDH/ AHHBBL
      DATA AHHBBL / 0.D0 /
#endif

      DATA OFIRST / .TRUE. /
      DATA AHH / 0.D0 /
      DATA AHI / 0.D0 /
      DATA AHG / 0.D0 /
      DATA IAH / 0 /

      IF (OINIT) THEN
         DO N = 1, NTDIM
#ifdef OPT_TRIPOLE
            CALL RSTADD(SX(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SX', 'OCN',
     &                                    1.D0,  0,  0 )
            CALL RSTADD(SY(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SY', 'OCN',
     &                                   -1.D0,  0,  0 )
            CALL RSTADD(SZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SZ', 'OCN',
     &                                    1.D0,  0,  0 )
            CALL RSTADD(SXX(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SXX', 'OCN',
     &                                    1.D0,  0,  0 )
            CALL RSTADD(SYY(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SYY', 'OCN',
     &                                    1.D0,  0,  0 )
            CALL RSTADD(SZZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SZZ', 'OCN',
     &                                    1.D0,  0,  0 )
            CALL RSTADD(SXY(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SXY', 'OCN',
     &                                   -1.D0,  0,  0 )
            CALL RSTADD(SXZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SXZ', 'OCN',
     &                                    1.D0,  0,  0 )
            CALL RSTADD(SYZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SYZ', 'OCN',
     &                                   -1.D0,  0,  0 )
#else
            CALL RSTADD(SX(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SX', 'OCN')
            CALL RSTADD(SY(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SY', 'OCN')
            CALL RSTADD(SZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SZ', 'OCN')
            CALL RSTADD(SXX(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SXX', 'OCN')
            CALL RSTADD(SYY(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SYY', 'OCN')
            CALL RSTADD(SZZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SZZ', 'OCN')
            CALL RSTADD(SXY(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SXY', 'OCN')
            CALL RSTADD(SXZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SXZ', 'OCN')
            CALL RSTADD(SYZ(1, 1, N), OEOF,
     &                  NXDIM, NYDIM, NZDIM, 'SYZ', 'OCN')
#endif
         END DO
         RETURN
      END IF

      IF (OFINAL) THEN
         DO N = 1, NTDIM
            CALL FINADD( SX(1, 1, N), NXDIM, NYDIM, NZDIM, 'SX',  'OCN')
            CALL FINADD( SY(1, 1, N), NXDIM, NYDIM, NZDIM, 'SY',  'OCN')
            CALL FINADD( SZ(1, 1, N), NXDIM, NYDIM, NZDIM, 'SZ',  'OCN')
            CALL FINADD(SXX(1, 1, N), NXDIM, NYDIM, NZDIM, 'SXX', 'OCN')
            CALL FINADD(SYY(1, 1, N), NXDIM, NYDIM, NZDIM, 'SYY', 'OCN')
            CALL FINADD(SZZ(1, 1, N), NXDIM, NYDIM, NZDIM, 'SZZ', 'OCN')
            CALL FINADD(SXY(1, 1, N), NXDIM, NYDIM, NZDIM, 'SXY', 'OCN')
            CALL FINADD(SXZ(1, 1, N), NXDIM, NYDIM, NZDIM, 'SXZ', 'OCN')
            CALL FINADD(SYZ(1, 1, N), NXDIM, NYDIM, NZDIM, 'SYZ', 'OCN')
         END DO
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFH, END=97)
 97      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFI, END=96)
 96      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFG, END=95)
 95      CONTINUE

         EPS = 1.D-20
         SQ3 = SQRT( 3.D0 )
         CI3 = 1.D0 / 3.D0

!----- initialization
         DO N = 1, NTDIM
            DO K = 1, NZDIM
               DO IJ = 1, NXYDIM
                  F0 (IJ, K)    = 0.D0
                  FM (IJ, K)    = 0.D0
                  FX (IJ, K)    = 0.D0
                  FY (IJ, K)    = 0.D0
                  FZ (IJ, K)    = 0.D0
                  FXX(IJ, K)    = 0.D0
                  FYY(IJ, K)    = 0.D0
                  FZZ(IJ, K)    = 0.D0
                  FXY(IJ, K)    = 0.D0
                  FXZ(IJ, K)    = 0.D0
                  FYZ(IJ, K)    = 0.D0
                  S0 (IJ, K, N) = 0.D0
                  SM (IJ, K, N) = 0.D0
                  VLMX(IJ, K)   = 0.D0
                  VLMY(IJ, K)   = 0.D0
                  VLMZ(IJ)      = 0.D0
                  ALF(IJ, K)    = 0.D0
               END DO
            END DO
         END DO

         IF (OEOF) THEN
            DO N = 1, NTDIM
               DO K = 1, NZDIM
                  DO IJ = 1, NXYDIM
                     SX (IJ, K, N) = 0.D0
                     SY (IJ, K, N) = 0.D0
                     SZ (IJ, K, N) = 0.D0
                     SXX(IJ, K, N) = 0.D0
                     SYY(IJ, K, N) = 0.D0
                     SZZ(IJ, K, N) = 0.D0
                     SXY(IJ, K, N) = 0.D0
                     SXZ(IJ, K, N) = 0.D0
                     SYZ(IJ, K, N) = 0.D0
                  END DO
               END DO
            END DO
         END IF

!---- area normal to u defined on UV-grid
         DO K = KSTR, KSTR+KZ-1
            DO IJ = 1, NXYDIM

               VLMX(IJ, K) =  HXU(IJ) * ZBOT * 0.5d0 * DS(K) * DX
               VLMY(IJ, K) =  HYU(IJ) * ZBOT * 0.5d0 * DS(K)

            END DO
         END DO

         DO K = KSTR+KZ, KEND
            DO IJ = 1, NXYDIM
            
               VLMX(IJ, K) =  HXU(IJ) * DZV(IJ, K) * 0.5d0 * DX
               VLMY(IJ, K) =  HYU(IJ) * DZV(IJ, K) * 0.5d0

            END DO
         END DO

         DO IJ = 1, NXYDIM
            VLMX(IJ, KSTR-1) = VLMX(IJ, KSTR)
            VLMX(IJ, KEND+1) = VLMX(IJ, KEND)
            VLMY(IJ, KSTR-1) = VLMY(IJ, KSTR)
            VLMY(IJ, KEND+1) = VLMY(IJ, KEND)
         END DO

         DO IJ = 1, NXYDIM
         
            VLMZ(IJ) = HXT(IJ) * DX * HYT(IJ) * DY(IJ)

         END DO

!---- 
         REWIND(NFPARM)
         READ(NFPARM, NMCAH , END = 90 )
 90      CONTINUE

         IF ( IAH .EQ. 0 ) THEN

            WRITE(JFPAR, *) 'Background horizontal diffusion :', AHH
            WRITE(JFPAR, *) 'Isopycnal diffusion             :', AHI
            WRITE(JFPAR, *) 'G-M thickness diffusion         :', AHG

            DO K = 1, NZDIM
               DO IJ = 1, NXYDIM
                  AHI3D(IJ, K) = AHI
                  AHG3D(IJ, K) = AHG
                  AHH3D(IJ, K) = AHH
               END DO
            END DO

         END IF

         IF ( IAH .EQ. 0 ) GO TO 80

         DO K = 1, NZDIM
            DO IJ = 1, NXYDIM
               AHH3D(IJ, K) = AHH
            END DO
         END DO

         WRITE(JFPAR, *) '  file name of AHI: ', CFAHI
         WRITE(JFPAR, *) '  file name of AHG: ', CFAHG

!---- reading diffusion coefficient
#ifdef OPT_PARALLEL

!---- AHI
         IF ( MYRANK .EQ. IROOT ) THEN

            CALL FILOPN( NFAHI, CFAHI, 'READ' )
            REWIND( NFAHI )
            READ( NFAHI ) CHEAD
            READ( NFAHI ) BUF3

            DO K = 1, NZ
               DO J = 1, NYG
                  DO I = 1, NXG

                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K )

                  END DO
               END DO
            END DO
            CALL FILCLS( NFAHI )

         END IF
         CALL SCATTER_3D( AHI3D, G3D )

!---- AHG
         IF ( MYRANK .EQ. IROOT ) THEN

            CALL FILOPN( NFAHG, CFAHG, 'READ' )
            REWIND( NFAHG )
            READ( NFAHG ) CHEAD
            READ( NFAHG ) BUF3

            DO K = 1, NZ
               DO J = 1, NYG
                  DO I = 1, NXG

                  G3D(IGSTR+I-1, JGSTR+J-1, KSTR+K-1) = BUF3(I, J, K )

                  END DO
               END DO
            END DO
            CALL FILCLS( NFAHG )

         END IF
         CALL SCATTER_3D( AHG3D, G3D )

#else

!---- AHI
         CALL FILOPN( NFAHI, CFAHI, 'READ' ) 
         REWIND( NFAHI )
         READ( NFAHI ) CHEAD
         READ( NFAHI ) BUF3
         DO K = KSTR, KEND
            DO J = JSTR, JEND
               DO I = ISTR, IEND

                  IJ = ( J - 1 ) * NXDIM + I
                  AHI3D(IJ,K) = BUF3(I-ISTR+1, J-JSTR+1, K-KSTR+1)
            
               END DO
            END DO
         END DO
         CALL FILCLS( NFAHI )

!---- AHG
         CALL FILOPN( NFAHG, CFAHG, 'READ' ) 
         REWIND( NFAHG )
         READ( NFAHG ) CHEAD
         READ( NFAHG ) BUF3
         DO K = KSTR, KEND
            DO J = JSTR, JEND
               DO I = ISTR, IEND

                  IJ = ( J - 1 ) * NXDIM + I
                  AHG3D(IJ,K) = BUF3(I-ISTR+1, J-JSTR+1, K-KSTR+1)
            
               END DO
            END DO
         END DO
         CALL FILCLS( NFAHG )

#endif

 80      CONTINUE

#ifdef OPT_BBL
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBBDH, END=98)
 98      CONTINUE
         WRITE(JFPAR, NMBBDH)
#endif

      END IF

      TSIV = 1.D0 / TS

      CALL DNSGRD(
     O             XDZDX,  YDZDY,  ZDZDX,  ZDZDY,
     O             XDTDZ,  YDTDZ,  ZDTDX,  ZDTDY,
     I                TY,     TX,     
     W               FTX,    FTY,    FTZ)

      DO 120 N = 1, NTDIM
         DO 110 K = 1, NZDIM
            DO 100 IJ = 1, NXYDIM
               ADT (IJ, K, N) = 0.D0
               FTX (IJ, K, N) = 0.D0
               FTY (IJ, K, N) = 0.D0
               FTZ (IJ, K, N) = 0.D0
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE

      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DIFFZ(IJ, K) = 0.D0
         END DO
      END DO

      DO IJ = 1, NXYDIM
         HZBOT(IJ) = HZ(IJ) + ZBOT
      END DO

!---- vertical velocity on sigma coordinate
      DO K = KSTR, KSTR+KZ-1
         DO IJ = 1, NXYDIM
            WZC(IJ, K) = W(IJ, K) * HZBOT(IJ)
            RZM(IJ, K) = 1.D0 / DSM(K) / HZBOT(IJ)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = 1, NXYDIM
            WZC(IJ, K) = W(IJ, K)
            RZM(IJ, K) = 1.D0 / DZM(IJ, K)
         END DO
      END DO

      CALL CHEKIN(WZC, 'WZC', NX, NY, NZ, NXYZDM, 'OCN')

!======  GM  isopycnal and diapycnal diffusion  ======
!---- z diffusion flux of GM
      DO 230 N = 1, NTDIM
         DO 220 K = KSTR+1, KEND
            KUU = K - 2
            KU  = K - 1
            KD  = K + 1
            DO 210 IJ = IJTSTR, IJTEND

               DIFFZ(IJ, K) =
     &         (  AHV(IJ, K)
     &          + AHI3D(IJ, K) * (  ZDZDX(IJ, K) * ZDZDX(IJ, K)
     &                            + ZDZDY(IJ, K) * ZDZDY(IJ, K) ) ) *
     &         RZM(IJ, K) * AMFTZ(IJ, K)

               FTZ(IJ, K, N) = 
     &         (  DIFFZ(IJ, K) * (TX(IJ, KU, N) - TX(IJ, K, N))
     &          - ( AHI3D(IJ, K) + AHG3D(IJ, K) ) *
     &            (  ZDZDX(IJ, K) * ZDTDX(IJ, K, N)
     &             + ZDZDY(IJ, K) * ZDTDY(IJ, K, N) )
     &         ) * AMFTZ(IJ, K)

  210       CONTINUE
  220    CONTINUE
  230 CONTINUE

!---- y diffusion flux of GM
      DO 340 N = 1, NTDIM

         DO 330 K = KSTR, KEND
            DO 320 IJ = IJTSTR, IJTEND+NXDIM

               IJLS = IJ + LS

               FTY(IJ, K, N) =
     &           (  ( AHH3D(IJ, K) + AHI3D(IJ, K) ) * RYM(IJLS) 
     &            / ( HYT(IJ) + HYT(IJLS) ) *
     &              ( TX(IJ, K, N) - TX(IJLS, K, N) ) * 2.D0
     &            - ( AHI3D(IJ, K) - AHG3D(IJ, K) ) 
     &            * YDZDY(IJ, K) * YDTDZ(IJ, K, N) )
     &            * ( HXU(IJLS) + HXU(IJ+LSW) ) * 0.5D0 * AMFTY(IJ, K)

  320       CONTINUE
  330    CONTINUE

  340 CONTINUE

!---- x diffusion flux of GM
      DO 440 N = 1, NTDIM

         DO 430 K = KSTR, KEND
            DO 420 IJ = IJTSTR, IJTEND+1

               IJLW = IJ + LW
               
               FTX(IJ, K, N) =
     &           (  ( AHH3D(IJ, K) + AHI3D(IJ, K) ) * RX 
     &            / ( HXT(IJ) + HXT(IJLW) ) 
     &            * ( TX(IJ, K, N) - TX(IJLW, K, N) ) * 2.D0
     &            - ( AHI3D(IJ, K) - AHG3D(IJ, K) )  
     &            * XDZDX(IJ, K) * XDTDZ(IJ, K, N) ) 
     &            * ( HYU(IJLW) + HYU(IJ+LSW) ) * 0.5D0 * AMFTX(IJ, K)

  420       CONTINUE
  430    CONTINUE
         
  440 CONTINUE

!---- divergence of diffusion fluxes
      DO 630 N = 1, NTDIM
         DO 620 K = KSTR, KEND
            DO 610 IJ = IJTSTR, IJTEND

               ADT(IJ, K, N) =
     &         (  (  (FTX(IJ+LE, K, N) - FTX(IJ, K, N)) * RX
     &             + (FTY(IJ+LN, K, N) - FTY(IJ, K, N)) * RY(IJ)) *
     &            RXT(IJ) * RYT(IJ)
     &          + FTZ(IJ, K, N) - FTZ(IJ, K+1, N)) / DZ(IJ, K)

  610       CONTINUE
  620    CONTINUE
  630 CONTINUE

!---- diffusion in BBL
#ifdef OPT_BBL

      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND
            K = NBOT(IJ)
            FTZ(IJ, KEND, N) = FTZ(IJ, K, N)
         END DO
      END DO

      DO N = 1, NTDIM

         DO IJ = IJTSTR, IJTEND+NXDIM
            IJLS = IJ + LS
            FTY(IJ, KEND, N) =
     &        AHHBBL * (TX(IJ, KEND, N) - TX(IJLS, KEND, N)) *
     &        RYM(IJ) * 
     &        (HXU(IJLS) + HXU(IJ+LSW)) / (HYT(IJ) + HYT(IJLS)) *
     &        AMFTY(IJ, KEND)
         END DO

         DO IJ = IJTSTR, IJTEND+1
            IJLW = IJ + LW
            FTX(IJ, KEND, N) = 
     &        AHHBBL * (TX(IJ, KEND, N) - TX(IJLW, KEND, N)) * RX *
     &        (HYU(IJLW) + HYU(IJ+LSW)) / (HXT(IJ) + HXT(IJLW)) *
     &        AMFTX(IJ, KEND)
         END DO

         DO IJ = IJTSTR, IJTEND
            ADT(IJ, KEND, N) =
     &      (  (  (  FTX(IJ+LE, KEND, N) - FTX(IJ, KEND, N)) * RX
     &          + (  FTY(IJ+LN, KEND, N) - FTY(IJ, KEND, N)) * RY(IJ)) *
     &         RXT(IJ) * RYT(IJ)
     &       + FTZ(IJ, KEND, N)) / DZ(IJ, KEND)
         END DO

      END DO

#endif

!---- SOM 

!---- mass contained in a tracer grid
      DO N = 1, NTDIM

          DO K = KSTR, KSTR+KZ-1
            DO IJ = 1, NXYDIM

               SM(IJ, K, N) = VLMZ(IJ) * HZBOT(IJ) * DS(K)

            END DO
         END DO

         DO K = KSTR+KZ, KEND
            DO IJ = 1, NXYDIM

               SM(IJ, K, N) = VLMZ(IJ) * DZ(IJ, K)

            END DO
         END DO

         DO K = KSTR, KEND
            DO IJ = 1, NXYDIM

               S0(IJ, K, N) = SM(IJ, K, N) * TX(IJ, K, N)
               
            END DO
         END DO

         DO IJ = 1, NXYDIM
            SM(IJ, KSTR-1, N) = SM(IJ, KSTR, N)
            SM(IJ, KEND+1, N) = SM(IJ, KEND, N)
            S0(IJ, KSTR-1, N) = S0(IJ, KSTR, N)
            S0(IJ, KEND+1, N) = S0(IJ, KEND, N)
         END DO

      END DO

!---- in X-direction

      DO K = KSTR, KEND
         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

            IJLW  = IJ + LW
            IJLSW = IJ + LSW
               
            UV(IJ, K) = ( UY(IJLW , K) * VLMY(IJLW, K) 
     &                  + UY(IJLSW, K) * VLMY(IJLSW, K) ) 
     &                * AMSKT(IJ, K) * AMSKT(IJLW, K) * DY(IJLW)

         END DO
      END DO

      DO N = 1, NTDIM

!---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW = IJ + LW * NINT( AMSKT(IJ+LW, K) )
               IJLE = IJ + LE * NINT( AMSKT(IJ+LE, K) )

               S0M  = S0(IJ, K, N) 
     &              - MIN( S0(IJLW, K, N) / SM(IJLW, K, N), 
     &                     S0(IJ  , K, N) / SM(IJ  , K, N), 
     &                     S0(IJLE, K, N) / SM(IJLE, K, N)  ) 
     &              * SM(IJ, K, N)

               S1M  = SQ3 * S0M
               SX(IJ, K, N) = MIN( S1M, MAX( - S1M, SX(IJ, K, N) ) )

               S0P = S0M * S0M
               SXP = SQRT( MAX( S0P - SX(IJ, K, N) * SX(IJ, K, N) * CI3, 
     &                          EPS ) )
          
               IF ( ABS( SX(IJ, K, N) ) < 1.5D0 * S0M ) THEN
             
                  SXX(IJ, K, N) = MIN( S0M + SXP, 
     &                                 MAX( ABS( SX(IJ, K, N) ) - S0M,
     &                                      SXX(IJ, K, N) ) )

               ELSE
                  
                  SXX(IJ, K, N) = MIN( S0M + SXP,
     &                                 MAX( S0M - SXP, SXX(IJ, K, N) ) )
             
               END IF

               SXY(IJ, K, N) = MIN( S0M, MAX( - S0M, SXY(IJ, K, N) ) )
               SXZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SXZ(IJ, K, N) ) )
          
            END DO
         END DO

!---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW = IJ + LW * NINT( AMSKT(IJ+LW, K) )
               IJLE = IJ + LE * NINT( AMSKT(IJ+LE, K) )

               S0M = - S0(IJ, K, N)
     &               + MAX( S0(IJLW, K, N) / SM(IJLW, K, N), 
     &                      S0(IJ  , K, N) / SM(IJ  , K, N), 
     &                      S0(IJLE, K, N) / SM(IJLE, K, N)  ) 
     &               * SM(IJ, K, N)

               S1M  = SQ3 * S0M
               SX(IJ, K, N) = MIN( S1M, MAX( - S1M, SX(IJ, K, N) ) )

               S0P = S0M * S0M
               SXP = SQRT( MAX( S0P - SX(IJ, K, N) * SX(IJ, K, N) * CI3, 
     &                          EPS ) )
          
               IF ( ABS( SX(IJ, K, N) ) < 1.5D0 * S0M ) THEN
             
                  SXX(IJ, K, N) = MIN( S0M + SXP, 
     &                                 MAX( ABS( SX(IJ, K, N) ) - S0M,
     &                                      SXX(IJ, K, N) ) )

               ELSE
                  
                  SXX(IJ, K, N) = MIN( S0M + SXP,
     &                                 MAX( S0M - SXP, SXX(IJ, K, N) ) )
             
               END IF

               SXY(IJ, K, N) = MIN( S0M, MAX( - S0M, SXY(IJ, K, N) ) )
               SXZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SXZ(IJ, K, N) ) )
          
            END DO
         END DO

!---- bug fix 2
!#ifdef OPT_PARALLEL
!         CALL SHIFT1( SX (:,:,N), NXDIM, NYDIM, NZDIM )
!         CALL SHIFT1( SXX(:,:,N), NXDIM, NYDIM, NZDIM )
!         CALL SHIFT1( SXY(:,:,N), NXDIM, NYDIM, NZDIM )
!         CALL SHIFT1( SXZ(:,:,N), NXDIM, NYDIM, NZDIM )
!#endif

!---- calculating ALF and MASS between box (I-1,J,K) <---> (I,J,K)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW  = IJ + LW
               
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN
                  FM(IJLW, K) = UV(IJ, K) * TS 
                  ALF(IJ, K)  = FM(IJLW, K) / SM(IJLW, K, N)
               ELSE
                  FM(IJLW, K) = - UV(IJ, K) * TS
                  ALF(IJ, K)  = FM(IJLW, K) / SM(IJ  , K, N)
               END IF

            END DO
         END DO

!---- calculating flux and moments between box (I-1,J,K) <---> (I,J,K)                  
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW  = IJ + LW

!---- flux from (I-1) to (I),  when U > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1

!---- moments be transported
                  F0(IJLW, K) = 
     &              ALF(IJ, K) * ( S0(IJLW, K, N)                                         
     &                    + ALF1 * ( 
     &                         SX(IJLW, K, N)                             
     &                     + ( ALF1 - ALF(IJ, K) ) * SXX(IJLW, K, N) ) )

                  FX(IJLW, K) = 
     &              ALFQ * ( SX(IJLW, K, N) 
     &                     + 3.D0 * ALF1 * SXX(IJLW, K, N) ) 

                  FXX(IJLW, K) = ALF(IJ, K) * ALFQ * SXX(IJLW, K, N)

                  FY (IJLW, K) = ALF(IJ, K) * ( SY(IJLW, K, N) 
     &                                        + ALF1 * SXY(IJLW, K, N) )
                  FZ (IJLW, K) = ALF(IJ, K) * ( SZ(IJLW, K, N) 
     &                                        + ALF1 * SXZ(IJLW, K, N) )

                  FXY(IJLW, K) = ALFQ       * SXY(IJLW, K, N)
                  FXZ(IJLW, K) = ALFQ       * SXZ(IJLW, K, N)
                  FYY(IJLW, K) = ALF(IJ, K) * SYY(IJLW, K, N)
                  FZZ(IJLW, K) = ALF(IJ, K) * SZZ(IJLW, K, N)
                  FYZ(IJLW, K) = ALF(IJ, K) * SYZ(IJLW, K, N)

                  FTX (IJ, K, N) = FTX(IJ, K, N) 
     &                           - F0(IJLW, K) * TSIV * RY(IJLW)

               ELSE
             
!---- flux from (I) to (I-1),  when U < 0

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1
             
                  F0 (IJLW, K) = 
     &               ALF(IJ, K) * ( S0(IJ, K, N)                                       
     &                     - ALF1 * (  SX(IJ, K, N)                            
     &                                - ( ALF1 - ALF(IJ, K) ) 
     &                                * SXX(IJ, K, N) ) )
                  FX (IJLW, K) = ALFQ * ( 
     &               SX(IJ, K, N) - 3.D0 * ALF1 * SXX(IJ, K, N) )
                  FXX(IJLW, K) = ALF(IJ, K) * ALFQ * SXX(IJ, K, N)
                  
                  FY (IJLW, K) = ALF(IJ, K)  * ( SY(IJ, K, N) 
     &                                         - ALF1 * SXY(IJ, K, N) )
                  FZ (IJLW, K) = ALF(IJ, K)  * ( SZ(IJ, K, N) 
     &                                         - ALF1 * SXZ(IJ, K, N) )

                  FXY(IJLW, K) = ALFQ       * SXY(IJ, K, N)
                  FXZ(IJLW, K) = ALFQ       * SXZ(IJ, K, N)
                  FYY(IJLW, K) = ALF(IJ, K) * SYY(IJ, K, N)
                  FZZ(IJLW, K) = ALF(IJ, K) * SZZ(IJ, K, N)
                  FYZ(IJLW, K) = ALF(IJ, K) * SYZ(IJ, K, N)
             
                  FTX (IJ, K, N) = FTX(IJ, K, N) 
     &                           + F0(IJLW, K) * TSIV * RY(IJLW)

               END IF
          
            END DO
         END DO

!---- calculating flux and moments between box (I-1,J,K) <---> (I,J,K)                  
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW  = IJ + LW

!---- flux from (I-1) to (I),  when U > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1

!---- moments remaining 
                  SM (IJLW, K, N) = SM(IJLW, K, N) - FM(IJLW, K)
                  S0 (IJLW, K, N) = S0(IJLW, K, N) - F0(IJLW, K)

                  SX (IJLW, K, N) = 
     &                ALF1Q * ( SX(IJLW, K, N) 
     &                        - 3.D0 * ALF(IJ, K) * SXX(IJLW, K, N) )
                  SXX(IJLW, K, N) = ALF1 * ALF1Q * SXX(IJLW, K, N)

                  SY (IJLW, K, N) = SY (IJLW, K, N) - FY (IJLW, K)
                  SYY(IJLW, K, N) = SYY(IJLW, K, N) - FYY(IJLW, K)

                  SZ (IJLW, K, N) = SZ (IJLW, K, N) - FZ (IJLW, K)
                  SZZ(IJLW, K, N) = SZZ(IJLW, K, N) - FZZ(IJLW, K)

                  SXY(IJLW, K, N) = ALF1Q * SXY(IJLW, K, N)
                  SXZ(IJLW, K, N) = ALF1Q * SXZ(IJLW, K, N)
                  SYZ(IJLW, K, N) = SYZ(IJLW, K, N) - FYZ(IJLW, K)
             
               ELSE
             
!---- flux from (I) to (I-1),  when U < 0

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1
             
                  SM (IJ, K, N) = SM(IJ, K, N) - FM(IJLW, K)
                  S0 (IJ, K, N) = S0(IJ, K, N) - F0(IJLW, K)

                  SX (IJ, K, N) = ALF1Q * ( 
     &               SX(IJ, K, N) + 3.D0 * ALF(IJ, K) * SXX(IJ, K, N)  )
                  SXX(IJ, K, N) = ALF1 * ALF1Q * SXX(IJ, K, N)

                  SY (IJ, K, N) = SY (IJ, K, N) - FY (IJLW, K)
                  SYY(IJ, K, N) = SYY(IJ, K, N) - FYY(IJLW, K)

                  SZ (IJ, K, N) = SZ (IJ, K, N) - FZ (IJLW, K)
                  SZZ(IJ, K, N) = SZZ(IJ, K, N) - FZZ(IJLW, K)

                  SXY(IJ, K, N) = ALF1Q * SXY(IJ, K, N)
                  SXZ(IJ, K, N) = ALF1Q * SXZ(IJ, K, N)
                  SYZ(IJ, K, N) = SYZ(IJ, K, N) - FYZ(IJLW, K)
             
               END IF
          
            END DO
         END DO

!---- put the temporary moments (fi) into appropriate neighboring boxes
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW  = IJ + LW

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  SM(IJ  , K, N) = SM(IJ, K, N) + FM(IJLW, K)
                  ALF(IJ, K)     = FM(IJLW, K) / SM(IJ, K, N)

               ELSE

                  SM(IJLW, K, N) = SM(IJLW, K, N) + FM(IJLW, K)
                  ALF(IJ, K)     = FM(IJLW, K) / SM(IJLW, K, N)

               END IF

            END DO
         END DO

         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLW  = IJ + LW

!---- flux from (I-1) to (I),  when U > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALF1 = 1.D0 - ALF(IJ, K)

                  TMP = ALF(IJ, K) * S0(IJ, K, N) 
     &                - ALF1 * F0(IJLW, K)

                  S0(IJ, K, N) = S0(IJ, K, N) + F0(IJLW, K)
                  SXX(IJ, K, N) = 
     &               ALF(IJ, K)  * ALF(IJ, K) * FXX(IJLW, K) 
     &             + ALF1 * ALF1 * SXX(IJ, K, N)
     &             + 5.D0 * (  ALF(IJ, K) * ALF1 * ( 
     &                           SX(IJ, K, N) - FX(IJLW, K) )
     &                       - ( ALF1 - ALF(IJ, K) ) * TMP )
                  SX (IJ, K, N) = 
     &               ALF(IJ, K) * FX(IJLW, K) 
     &             + ALF1 * SX(IJ, K, N) + 3.D0 * TMP

                  SXY(IJ, K, N) = 
     &               ALF(IJ, K)  * FXY(IJLW, K) 
     &             + ALF1 * SXY(IJ, K, N)
     &             + 3.D0 * ( ALF(IJ, K) * SY(IJ, K, N) 
     &                      - ALF1 * FY(IJLW, K) )
                  SXZ(IJ, K, N) = 
     &               ALF(IJ, K)  * FXZ(IJLW, K) 
     &             + ALF1 * SXZ(IJ, K, N)
     &             + 3.D0 * ( ALF(IJ, K) * SZ(IJ, K, N) 
     &                      - ALF1 * FZ(IJLW, K) )

                  SY (IJ, K, N) = SY (IJ, K, N) + FY (IJLW, K)
                  SYY(IJ, K, N) = SYY(IJ, K, N) + FYY(IJLW, K)

                  SZ (IJ, K, N) = SZ (IJ, K, N) + FZ (IJLW, K)
                  SZZ(IJ, K, N) = SZZ(IJ, K, N) + FZZ(IJLW, K)

                  SYZ(IJ, K, N) = SYZ(IJ, K, N) + FYZ(IJLW, K)

!---- flux from (I) to (I-1),  when U < 0
               ELSE

                  ALF1 = 1.D0 - ALF(IJ, K)

                  TMP = - ALF(IJ, K) * S0(IJLW, K, N) 
     &                  + ALF1 * F0(IJLW, K)

                  S0(IJLW, K, N) = S0(IJLW, K, N) + F0(IJLW, K)
                  SXX(IJLW, K, N) = 
     &                 ALF(IJ, K) * ALF(IJ, K) * FXX(IJLW, K)
     &               + ALF1 * ALF1 * SXX(IJLW, K, N)
     &               + 5.D0 * (  ALF(IJ, K) * ALF1 * ( - SX(IJLW, K, N) 
     &                                          + FX(IJLW, K)    )
     &                         + ( ALF1 - ALF(IJ, K) ) * TMP  )
                  SX (IJLW, K, N) = 
     &                 ALF(IJ, K) * FX(IJLW, K) + ALF1 * SX(IJLW, K, N) 
     &               + 3.D0 * TMP

                  SXY(IJLW, K, N) = 
     &                 ALF(IJ, K) * FXY(IJLW, K) 
     &               + ALF1 * SXY(IJLW, K, N)
     &               + 3.D0 * ( ALF1 * FY(IJLW, K) 
     &                        - ALF(IJ, K) * SY(IJLW, K, N) )
                  SXZ(IJLW, K, N) = 
     &                 ALF(IJ, K) * FXZ(IJLW, K) 
     &               + ALF1 * SXZ(IJLW, K, N)
     &               + 3.D0 * ( ALF1 * FZ(IJLW, K) 
     &                        - ALF(IJ, K) * SZ(IJLW, K, N) )

                  SY (IJLW, K, N) = SY (IJLW, K, N) + FY (IJLW, K)
                  SZ (IJLW, K, N) = SZ (IJLW, K, N) + FZ (IJLW, K)

                  SYY(IJLW, K, N) = SYY(IJLW, K, N) + FYY(IJLW, K)
                  SZZ(IJLW, K, N) = SZZ(IJLW, K, N) + FZZ(IJLW, K)
                  SYZ(IJLW, K, N) = SYZ(IJLW, K, N) + FYZ(IJLW, K)

               END IF

            END DO
         END DO

      END DO

!---- Y-direction

      DO K = KSTR, KEND
         DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

            IJLS  = IJ + LS
            IJLSW = IJ + LSW

            UV(IJ, K) = ( VY(IJLS , K) * VLMX(IJLS , K)
     &                  + VY(IJLSW, K) * VLMX(IJLSW, K) )
     &                * AMSKT(IJ, K) * AMSKT(IJLS, K) 

         END DO
      END DO

      DO N = 1, NTDIM

!---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS * NINT( AMSKT(IJ+LS, K) )
               IJLN  = IJ + LN * NINT( AMSKT(IJ+LN, K) )
                
               S0M = S0(IJ, K, N) 
     &             - MIN( S0(IJLS, K, N) / SM(IJLS, K, N), 
     &                    S0(IJ  , K, N) / SM(IJ  , K, N), 
     &                    S0(IJLN, K, N) / SM(IJLN, K, N)  ) 
     &             * SM(IJ, K, N)

               S1M = SQ3 * S0M
               SY(IJ, K, N) = MIN( S1M, MAX( - S1M, SY(IJ, K, N) ) )

               S0P = S0M * S0M
               SXP = SQRT( MAX( S0P - SY(IJ, K, N) * SY(IJ, K, N) * CI3, 
     &                          EPS ) )
               
               IF ( ABS( SY(IJ, K, N) ) < 1.5D0 * S0M ) THEN

                  SYY(IJ, K, N) = MIN( S0M + SXP, 
     &                                 MAX( ABS( SY(IJ, K, N) ) - S0M,
     &                                      SYY(IJ, K, N) ) )

               ELSE
                  
                  SYY(IJ, K, N) = MIN( S0M + SXP,
     &                                 MAX( S0M - SXP, SYY(IJ, K, N) ) )

               END IF
             
               SXY(IJ, K, N) = MIN( S0M, MAX( - S0M, SXY(IJ, K, N) ) )
               SYZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SYZ(IJ, K, N) ) )

            END DO
         END DO

!---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS * NINT( AMSKT(IJ+LS, K) )
               IJLN  = IJ + LN * NINT( AMSKT(IJ+LN, K) )

               S0M = - S0(IJ, K, N) 
     &               + MAX( S0(IJLS, K, N) / SM(IJLS, K, N), 
     &                      S0(IJ  , K, N) / SM(IJ  , K, N), 
     &                      S0(IJLN, K, N) / SM(IJLN, K, N)  ) 
     &               * SM(IJ, K, N)

               S1M = SQ3 * S0M
               SY(IJ, K, N) = MIN( S1M, MAX( - S1M, SY(IJ, K, N) ) )

               S0P = S0M * S0M
               SXP = SQRT( MAX( S0P - SY(IJ, K, N) * SY(IJ, K, N) * CI3, 
     &                          EPS ) )
               
               IF ( ABS( SY(IJ, K, N) ) < 1.5D0 * S0M ) THEN

                  SYY(IJ, K, N) = MIN( S0M + SXP, 
     &                                 MAX( ABS( SY(IJ, K, N) ) - S0M,
     &                                      SYY(IJ, K, N) ) )

               ELSE
                  
                  SYY(IJ, K, N) = MIN( S0M + SXP,
     &                                 MAX( S0M - SXP, SYY(IJ, K, N) ) )

               END IF
             
               SXY(IJ, K, N) = MIN( S0M, MAX( - S0M, SXY(IJ, K, N) ) )
               SYZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SYZ(IJ, K, N) ) )

            END DO
         END DO

!---- bug fix 2
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT1( SY(:,:,N),
     D                NXDIM, NYDIM, NZDIM,
     I               -1.D0,      0,      0 )
         CALL SHIFT1( SYY(:,:,N),
     D                NXDIM, NYDIM, NZDIM,
     I                1.D0,      0,      0 )
         CALL SHIFT1( SXY(:,:,N),
     D                NXDIM, NYDIM, NZDIM,
     I               -1.D0,      0,      0 )
         CALL SHIFT1( SYZ(:,:,N),
     D                NXDIM, NYDIM, NZDIM,
     I               -1.D0,      0,      0 )
#else
      CALL SHIFT1( SY (:,:,N), NXDIM, NYDIM, NZDIM )
      CALL SHIFT1( SYY(:,:,N), NXDIM, NYDIM, NZDIM )
      CALL SHIFT1( SXY(:,:,N), NXDIM, NYDIM, NZDIM )
      CALL SHIFT1( SYZ(:,:,N), NXDIM, NYDIM, NZDIM )
#endif
#endif

!---- calculating ALF  and MASS between box (I,J-1,K) <---> (I,J,K)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN
                  FM(IJLS, K) = UV(IJ, K) * TS
                  ALF(IJ, K)  = FM(IJLS, K) / SM(IJLS, K, N)
               ELSE
                  FM(IJLS, K) = - UV(IJ, K) * TS
                  ALF(IJ, K)  = FM(IJLS, K) / SM(IJ, K, N)
               END IF

            END DO
         END DO

!---- calculating flux between box (I,J-1,K) <---> (I,J,K)
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS

!---- flux from (J-1) to (J),  when V > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1

!---- moments be transported
                  F0(IJLS, K) = 
     &              ALF(IJ, K) * ( S0(IJLS, K, N)                                         
     &                  + ALF1 * ( 
     &                     SY(IJLS, K, N)
     &                   + ( ALF1 - ALF(IJ, K) ) * SYY(IJLS, K, N) ) )

                  FY(IJLS, K) = 
     &              ALFQ * ( SY(IJLS, K, N) 
     &                     + 3.D0 * ALF1 * SYY(IJLS, K, N) ) 

                  FYY(IJLS, K) = ALF(IJ, K) * ALFQ * SYY(IJLS, K, N)

                  FX (IJLS, K) = ALF(IJ, K) * ( SX(IJLS, K, N) 
     &                                        + ALF1 * SXY(IJLS, K, N) )
                  FZ (IJLS, K) = ALF(IJ, K) * ( SZ(IJLS, K, N) 
     &                                        + ALF1 * SYZ(IJLS, K, N) )

                  FXY(IJLS, K) = ALFQ       * SXY(IJLS, K, N)
                  FYZ(IJLS, K) = ALFQ       * SYZ(IJLS, K, N)
                  FXX(IJLS, K) = ALF(IJ, K) * SXX(IJLS, K, N)
                  FZZ(IJLS, K) = ALF(IJ, K) * SZZ(IJLS, K, N)
                  FXZ(IJLS, K) = ALF(IJ, K) * SXZ(IJLS, K, N)
             
                  FTY (IJ, K, N) = FTY(IJ, K, N)
     &                           - F0(IJLS, K) * TSIV * RX

               ELSE
             
!---- flux from (J) to (J-1),  when V < 0

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1
             
                  F0 (IJLS, K) = 
     &               ALF(IJ, K) * ( S0(IJ, K, N)                                       
     &                     - ALF1 * (  SY(IJ, K, N)                            
     &                               - ( ALF1 - ALF(IJ, K) ) 
     &                                * SYY(IJ, K, N) ) )
                  FY (IJLS, K) = ALFQ * ( 
     &               SY(IJ, K, N) - 3.D0 * ALF1 * SYY(IJ, K, N) )
                  FYY(IJLS, K) = ALF(IJ, K) * ALFQ * SYY(IJ, K, N)
                  
                  FX (IJLS, K) = ALF(IJ, K) * ( SX(IJ, K, N) 
     &                                        - ALF1 * SXY(IJ, K, N) )
                  FZ (IJLS, K) = ALF(IJ, K) * ( SZ(IJ, K, N) 
     &                                        - ALF1 * SYZ(IJ, K, N) )

                  FXY(IJLS, K) = ALFQ       * SXY(IJ, K, N)
                  FYZ(IJLS, K) = ALFQ       * SYZ(IJ, K, N)
                  FXX(IJLS, K) = ALF(IJ, K) * SXX(IJ, K, N)
                  FZZ(IJLS, K) = ALF(IJ, K) * SZZ(IJ, K, N)
                  FXZ(IJLS, K) = ALF(IJ, K) * SXZ(IJ, K, N)
             
                  FTY (IJ, K, N) = FTY(IJ, K, N) 
     &                           + F0(IJLS, K) * TSIV * RX

               END IF
          
            END DO
         END DO

!---- calculating moments
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS

!---- flux from (J-1) to (J),  when V > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1

!---- moments remaining 
                  SM (IJLS, K, N) = SM(IJLS, K, N) - FM(IJLS, K)
                  S0 (IJLS, K, N) = S0(IJLS, K, N) - F0(IJLS, K)

                  SY (IJLS, K, N) = 
     &                ALF1Q * ( SY(IJLS, K, N) 
     &                        - 3.D0 * ALF(IJ, K) * SYY(IJLS, K, N) )
                  SYY(IJLS, K, N) = ALF1 * ALF1Q * SYY(IJLS, K, N)

                  SX (IJLS, K, N) = SX (IJLS, K, N) - FX (IJLS, K)
                  SXX(IJLS, K, N) = SXX(IJLS, K, N) - FXX(IJLS, K)

                  SZ (IJLS, K, N) = SZ (IJLS, K, N) - FZ (IJLS, K)
                  SZZ(IJLS, K, N) = SZZ(IJLS, K, N) - FZZ(IJLS, K)

                  SXY(IJLS, K, N) = ALF1Q * SXY(IJLS, K, N)
                  SYZ(IJLS, K, N) = ALF1Q * SYZ(IJLS, K, N)
                  SXZ(IJLS, K, N) = SXZ(IJLS, K, N) - FXZ(IJLS, K)
             
               ELSE
             
!---- flux from (J) to (J-1),  when V < 0

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1
             
                  SM (IJ, K, N) = SM(IJ, K, N) - FM(IJLS, K)
                  S0 (IJ, K, N) = S0(IJ, K, N) - F0(IJLS, K)

                  SY (IJ, K, N) = ALF1Q * ( 
     &                    SY(IJ, K, N) 
     &                  + 3.D0 * ALF(IJ, K) * SYY(IJ, K, N)  )
                  SYY(IJ, K, N) = ALF1 * ALF1Q * SYY(IJ, K, N)

                  SX (IJ, K, N) = SX (IJ, K, N) - FX (IJLS, K)
                  SXX(IJ, K, N) = SXX(IJ, K, N) - FXX(IJLS, K)

                  SZ (IJ, K, N) = SZ (IJ, K, N) - FZ (IJLS, K)
                  SZZ(IJ, K, N) = SZZ(IJ, K, N) - FZZ(IJLS, K)

                  SXY(IJ, K, N) = ALF1Q * SXY(IJ, K, N)
                  SYZ(IJ, K, N) = ALF1Q * SYZ(IJ, K, N)
                  SXZ(IJ, K, N) = SXZ(IJ, K, N) - FXZ(IJLS, K)
             
               END IF
          
            END DO
         END DO

!---- put the temporary moments (fi) into appropriate neighboring boxes
         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN
                  SM(IJ, K, N) = SM(IJ, K, N) + FM(IJLS, K)
                  ALF(IJ, K)   = FM(IJLS, K) / SM(IJ, K, N)
               ELSE
                  SM(IJLS, K, N) = SM(IJLS, K, N) + FM(IJLS, K)
                  ALF(IJ, K)     = FM(IJLS, K) / SM(IJLS, K, N)
               END IF

            END DO
         END DO

         DO K = KSTR, KEND
            DO IJ = IJTSTR-NXDIM-1, IJTEND+NXDIM+1

               IJLS  = IJ + LS

!---- flux from (J-1) to (J),  when V > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALF1 = 1.D0 - ALF(IJ, K)

                  TMP = ALF(IJ, K) * S0(IJ, K, N) - ALF1 * F0(IJLS, K)

                  S0(IJ, K, N) = S0(IJ, K, N) + F0(IJLS, K)
                  SYY(IJ, K, N) = 
     &               ALF(IJ, K) * ALF(IJ, K) * FYY(IJLS, K) 
     &             + ALF1 * ALF1 * SYY(IJ, K, N)
     &             + 5.D0 * (  ALF(IJ, K) * ALF1 * ( 
     &                           SY(IJ, K, N) - FY(IJLS, K) )
     &                       - ( ALF1 - ALF(IJ, K) ) * TMP )
                  SY (IJ, K, N) = 
     &               ALF(IJ, K) * FY(IJLS, K) 
     &             + ALF1 * SY(IJ, K, N) + 3.D0 * TMP

                  SXY(IJ, K, N) = 
     &               ALF(IJ, K) * FXY(IJLS, K) 
     &             + ALF1 * SXY(IJ, K, N)
     &             + 3.D0 * ( ALF(IJ, K) * SX(IJ, K, N) 
     &                      - ALF1 * FX(IJLS, K) )
                  SYZ(IJ, K, N) = 
     &               ALF(IJ, K) * FYZ(IJLS, K) 
     &             + ALF1 * SYZ(IJ, K, N)
     &             + 3.D0 * ( ALF(IJ, K) * SZ(IJ, K, N) 
     &                      - ALF1 * FZ(IJLS, K) )

                  SX (IJ, K, N) = SX (IJ, K, N) + FX (IJLS, K)
                  SXX(IJ, K, N) = SXX(IJ, K, N) + FXX(IJLS, K)

                  SZ (IJ, K, N) = SZ (IJ, K, N) + FZ (IJLS, K)
                  SZZ(IJ, K, N) = SZZ(IJ, K, N) + FZZ(IJLS, K)

                  SXZ(IJ, K, N) = SXZ(IJ, K, N) + FXZ(IJLS, K)

!---- flux from (J) to (J-1),  when V < 0
               ELSE

                  ALF1 = 1.D0 - ALF(IJ, K)

                  TMP = - ALF(IJ, K) * S0(IJLS, K, N) 
     &                  + ALF1 * F0(IJLS, K)

                  S0(IJLS, K, N) = S0(IJLS, K, N) + F0(IJLS, K)
                  SYY(IJLS, K, N) = 
     &                 ALF(IJ, K) * ALF(IJ, K) * FYY(IJLS, K)
     &               + ALF1 * ALF1 * SYY(IJLS, K, N)
     &               + 5.D0 * (  ALF(IJ, K) * ALF1 * ( 
     &                         - SY(IJLS, K, N) + FY(IJLS, K) )
     &                         + ( ALF1 - ALF(IJ, K) ) * TMP  )
                  SY (IJLS, K, N) = 
     &                 ALF(IJ, K) * FY(IJLS, K) + ALF1 * SY(IJLS, K, N) 
     &               + 3.D0 * TMP

                  SXY(IJLS, K, N) = 
     &                 ALF(IJ, K) * FXY(IJLS, K) 
     &               + ALF1 * SXY(IJLS, K, N)
     &               + 3.D0 * ( ALF1       * FX(IJLS, K) 
     &                        - ALF(IJ, K) * SX(IJLS, K, N) )
                  SYZ(IJLS, K, N) = 
     &                 ALF(IJ, K) * FYZ(IJLS, K) 
     &               + ALF1 * SYZ(IJLS, K, N)
     &               + 3.D0 * ( ALF1       * FZ(IJLS, K) 
     &                        - ALF(IJ, K) * SZ(IJLS, K, N) )

                  SX (IJLS, K, N) = SX (IJLS, K, N) + FX (IJLS, K)
                  SZ (IJLS, K, N) = SZ (IJLS, K, N) + FZ (IJLS, K)
                  SXX(IJLS, K, N) = SXX(IJLS, K, N) + FXX(IJLS, K)
                  SZZ(IJLS, K, N) = SZZ(IJLS, K, N) + FZZ(IJLS, K)
                  SXZ(IJLS, K, N) = SXZ(IJLS, K, N) + FXZ(IJLS, K)

               END IF

            END DO
         END DO

      END DO

#ifdef OPT_BBL
!---- keeping BBL variables consistent at two levels
      CALL STBBT2( S0 )
      CALL STBBT2( SM )
      CALL STBBT2( SX )
      CALL STBBT2( SY )
      CALL STBBT2( SZ )
      CALL STBBT2( SXX )
      CALL STBBT2( SYY )
      CALL STBBT2( SZZ )
      CALL STBBT2( SXY )
      CALL STBBT2( SXZ )
      CALL STBBT2( SYZ )
#endif

!---- Z-direction

      DO K = KSTR, KEND
         DO IJ = IJTSTR, IJTEND

!            UV(IJ, K) = - WZC(IJ, K) * VLMZ(IJ)
            UV(IJ, K) = - WZC(IJ, K) * VLMZ(IJ)
     &                * AMSKT(IJ, K) * AMSKT(IJ, K-1)
               
         END DO
      END DO

#ifdef OPT_BBL
      DO IJ = IJTSTR, IJTEND
         K = NBOT(IJ)
         UV(IJ, K) = - WZC(IJ, K) * VLMZ(IJ) * AMSKTB(IJ)
      END DO
#endif

      DO N = 1, NTDIM

!---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
         DO K = KSTR, KEND

            KU = MAX( K - 1, KSTR )
            KD = MIN( K + 1, KEND )
            
            DO IJ = IJTSTR, IJTEND

               S0M = S0(IJ, K, N) 
     &             - MIN( S0(IJ, KU, N) / SM(IJ, KU, N), 
     &                    S0(IJ, K , N) / SM(IJ, K , N),
     &                    S0(IJ, KD, N) / SM(IJ, KD, N)  ) 
     &             * SM(IJ, K, N)

               S1M = SQ3 * S0M
               SZ(IJ, K, N) = MIN( S1M, MAX( - S1M, SZ(IJ, K, N) ) )
               S0P = S0M * S0M
               SXP = SQRT( MAX( S0P - SZ(IJ, K, N) * SZ(IJ, K, N) * CI3, 
     &                          EPS ) )
               
               IF ( ABS( SZ(IJ, K, N) ) < 1.5D0 * S0M ) THEN

                  SZZ(IJ, K, N) = MIN( S0M + SXP, 
     &                                 MAX( ABS( SZ(IJ, K, N) ) - S0M,
     &                                      SZZ(IJ, K, N) ) )

               ELSE
                  
                  SZZ(IJ, K, N) = MIN( S0M + SXP,
     &                                 MAX( S0M - SXP, SZZ(IJ, K, N) ) )

               END IF
             
               SXZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SXZ(IJ, K, N) ) )
               SYZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SYZ(IJ, K, N) ) )

            END DO
         END DO

!---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
         DO K = KSTR, KEND

            KU = MAX( K - 1, KSTR )
            KD = MIN( K + 1, KEND )

            DO IJ = IJTSTR, IJTEND

               S0M = - S0(IJ, K, N) 
     &               + MAX( S0(IJ, KU, N) / SM(IJ, KU, N), 
     &                      S0(IJ, K,  N) / SM(IJ, K , N),
     &                      S0(IJ, KD, N) / SM(IJ, KD, N)  ) 
     &               * SM(IJ, K, N)

               S1M = SQ3 * S0M
               SZ(IJ, K, N) = MIN( S1M, MAX( - S1M, SZ(IJ, K, N) ) )
               S0P = S0M * S0M
               SXP = SQRT( MAX( S0P - SZ(IJ, K, N) * SZ(IJ, K, N) * CI3, 
     &                          EPS ) )
               
               IF ( ABS( SZ(IJ, K, N) ) < 1.5D0 * S0M ) THEN

                  SZZ(IJ, K, N) = MIN( S0M + SXP, 
     &                                 MAX( ABS( SZ(IJ, K, N) ) - S0M,
     &                                      SZZ(IJ, K, N) ) )

               ELSE
                  
                  SZZ(IJ, K, N) = MIN( S0M + SXP,
     &                                 MAX( S0M - SXP, SZZ(IJ, K, N) ) )

               END IF
             
               SXZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SXZ(IJ, K, N) ) )
               SYZ(IJ, K, N) = MIN( S0M, MAX( - S0M, SYZ(IJ, K, N) ) )

            END DO
         END DO

!---- calculating ALF
         DO K = KSTR, KEND

            KU = K - 1

            DO IJ = IJTSTR, IJTEND

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN
                  FM(IJ, KU) = UV(IJ, K) * TS
                  ALF(IJ, K) = FM(IJ, KU) / SM(IJ, KU, N)
               ELSE
                  FM(IJ, KU) = - UV(IJ, K) * TS
                  ALF(IJ, K) = FM(IJ, KU) / SM(IJ, K, N)
               END IF

            END DO
         END DO

!---- calculating flux between box (I,J,K-1) <---> (I,J,K)
         DO K = KSTR, KEND

            KU = K - 1

            DO IJ = IJTSTR, IJTEND

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

!---- flux from (K-1) to (K),  when W > 0
                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1

!---- moments be transported
                  F0(IJ, KU) = 
     &              ALF(IJ, K) * ( S0(IJ, KU, N)                                         
     &                    + ALF1 * ( 
     &                         SZ(IJ, KU, N)                             
     &                       + ( ALF1 - ALF(IJ, K) ) 
     &                       * SZZ(IJ, KU, N) ) )

                  FZ(IJ, KU) = 
     &              ALFQ * ( SZ(IJ, KU, N) 
     &                     + 3.D0 * ALF1 * SZZ(IJ, KU, N) ) 

                  FZZ(IJ, KU) = ALF(IJ, K) * ALFQ * SZZ(IJ, KU, N)

                  FX (IJ, KU) = ALF(IJ, K) * ( SX(IJ, KU, N) 
     &                                       + ALF1 * SXZ(IJ, KU, N) )
                  FY (IJ, KU) = ALF(IJ, K) * ( SY(IJ, KU, N) 
     &                                       + ALF1 * SYZ(IJ, KU, N) )

                  FXZ(IJ, KU) = ALFQ       * SXZ(IJ, KU, N)
                  FYZ(IJ, KU) = ALFQ       * SYZ(IJ, KU, N)
                  FXX(IJ, KU) = ALF(IJ, K) * SXX(IJ, KU, N)
                  FYY(IJ, KU) = ALF(IJ, K) * SYY(IJ, KU, N)
                  FXY(IJ, KU) = ALF(IJ, K) * SXY(IJ, KU, N)
                  
                  FTZ (IJ, K, N) = FTZ(IJ, K, N) 
     &                           + F0(IJ, KU) * TSIV / VLMZ(IJ)

               ELSE
             
!---- flux from (K) to (K-1),  when W < 0

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1
             
                  F0 (IJ, KU) = 
     &               ALF(IJ, K) * ( S0(IJ, K, N)                                       
     &                     - ALF1 * (  SZ(IJ, K, N)                            
     &                                - ( ALF1 - ALF(IJ, K) ) 
     &                                * SZZ(IJ, K, N) ) )
                  FZ (IJ, KU) = ALFQ * ( 
     &               SZ(IJ, K, N) - 3.D0 * ALF1 * SZZ(IJ, K, N) )
                  FZZ(IJ, KU) = ALF(IJ, K) * ALFQ * SZZ(IJ, K, N)
                  
                  FX (IJ, KU) = ALF(IJ, K)  * ( SX(IJ, K, N) 
     &                                        - ALF1 * SXZ(IJ, K, N) )
                  FY (IJ, KU) = ALF(IJ, K)  * ( SY(IJ, K, N) 
     &                                        - ALF1 * SYZ(IJ, K, N) )

                  FXZ(IJ, KU) = ALFQ       * SXZ(IJ, K, N)
                  FYZ(IJ, KU) = ALFQ       * SYZ(IJ, K, N)
                  FXX(IJ, KU) = ALF(IJ, K) * SXX(IJ, K, N)
                  FYY(IJ, KU) = ALF(IJ, K) * SYY(IJ, K, N)
                  FXY(IJ, KU) = ALF(IJ, K) * SXY(IJ, K, N)
             
                  FTZ (IJ, K, N) = FTZ(IJ, K, N) 
     &                           - F0(IJ, KU) * TSIV / VLMZ(IJ)

               END IF
          
            END DO
         END DO

!---- calculating flux between box (I,J,K-1) <---> (I,J,K)
         DO K = KSTR, KEND

            KU = K - 1

            DO IJ = IJTSTR, IJTEND

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

!---- flux from (K-1) to (K),  when W > 0
                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1

!---- moments remaining 
                  SM (IJ, KU, N) = SM(IJ, KU, N) - FM(IJ, KU)
                  S0 (IJ, KU, N) = S0(IJ, KU, N) - F0(IJ, KU)

                  SZ (IJ, KU, N) = 
     &                ALF1Q * ( SZ(IJ, KU, N) 
     &                        - 3.D0 * ALF(IJ, K) * SZZ(IJ, KU, N) )
                  SZZ(IJ, KU, N) = ALF1 * ALF1Q * SZZ(IJ, KU, N)

                  SX (IJ, KU, N) = SX (IJ, KU, N) - FX (IJ, KU)
                  SXX(IJ, KU, N) = SXX(IJ, KU, N) - FXX(IJ, KU)

                  SY (IJ, KU, N) = SY (IJ, KU, N) - FY (IJ, KU)
                  SYY(IJ, KU, N) = SYY(IJ, KU, N) - FYY(IJ, KU)

                  SXZ(IJ, KU, N) = ALF1Q * SXZ(IJ, KU, N)
                  SYZ(IJ, KU, N) = ALF1Q * SYZ(IJ, KU, N)
                  SXY(IJ, KU, N) = SXY(IJ, KU, N) - FXY(IJ, KU)
             
               ELSE
             
!---- flux from (K) to (K-1),  when W < 0

                  ALFQ  = ALF(IJ, K) * ALF(IJ, K)
                  ALF1  = 1.D0 - ALF(IJ, K)
                  ALF1Q = ALF1 * ALF1
             
                  SM (IJ, K, N) = SM(IJ, K, N) - FM(IJ, KU)
                  S0 (IJ, K, N) = S0(IJ, K, N) - F0(IJ, KU)

                  SZ (IJ, K, N) = ALF1Q * ( 
     &                SZ(IJ, K, N) 
     &              + 3.D0 * ALF(IJ, K) * SZZ(IJ, K, N)  )
                  SZZ(IJ, K, N) = ALF1 * ALF1Q * SZZ(IJ, K, N)

                  SX (IJ, K, N) = SX (IJ, K, N) - FX (IJ, KU)
                  SXX(IJ, K, N) = SXX(IJ, K, N) - FXX(IJ, KU)

                  SY (IJ, K, N) = SY (IJ, K, N) - FY (IJ, KU)
                  SYY(IJ, K, N) = SYY(IJ, K, N) - FYY(IJ, KU)

                  SXZ(IJ, K, N) = ALF1Q * SXZ(IJ, K, N)
                  SYZ(IJ, K, N) = ALF1Q * SYZ(IJ, K, N)
                  SXY(IJ, K, N) = SXY(IJ, K, N) - FXY(IJ, KU)
             
               END IF
          
            END DO
         END DO

!---- put the temporary moments (fi) into appropriate neighboring boxes
         DO K = KSTR, KEND

            KU = K - 1

            DO IJ = IJTSTR, IJTEND

               IF ( UV(IJ, K) .GT. 0.D0 ) THEN
                  SM(IJ, K, N) = SM(IJ, K, N) + FM(IJ, KU)
                  ALF(IJ, K)   = FM(IJ, KU) / SM(IJ, K, N)
               ELSE
                  SM(IJ, KU, N) = SM(IJ, KU, N) + FM(IJ, KU)
                  ALF(IJ, K)    = FM(IJ, KU) / SM(IJ, KU, N)
               END IF

            END DO
         END DO

         DO K = KSTR, KEND

            KU = K - 1

            DO IJ = IJTSTR, IJTEND

!---- flux from (K-1) to (K),  when W > 0
               IF ( UV(IJ, K) .GT. 0.D0 ) THEN

                  ALF1 = 1.D0 - ALF(IJ, K)

                  TMP = ALF(IJ, K) * S0(IJ, K, N) - ALF1 * F0(IJ, KU)

                  S0(IJ, K, N) = S0(IJ, K, N) + F0(IJ, KU)

                  SZZ(IJ, K, N) = 
     &               ALF(IJ, K) * ALF(IJ, K) * FZZ(IJ, KU) 
     &             + ALF1 * ALF1 * SZZ(IJ, K, N)
     &             + 5.D0 * (  ALF(IJ, K) * ALF1 * ( 
     &                           SZ(IJ, K, N) - FZ(IJ, KU) )
     &                       - ( ALF1 - ALF(IJ, K) ) * TMP )
                  SZ (IJ, K, N) = 
     &               ALF(IJ, K) * FZ(IJ, KU) 
     &             + ALF1       * SZ(IJ, K, N) + 3.D0 * TMP

                  SXZ(IJ, K, N) = 
     &               ALF(IJ, K) * FXZ(IJ, KU) 
     &             + ALF1       * SXZ(IJ, K, N)
     &             + 3.D0 * ( ALF(IJ, K) * SX(IJ, K, N) 
     &                      - ALF1 * FX(IJ, KU) )
                  SYZ(IJ, K, N) = 
     &               ALF(IJ, K) * FYZ(IJ, KU) 
     &             + ALF1       * SYZ(IJ, K, N)
     &             + 3.D0 * ( ALF(IJ, K) * SY(IJ, K, N) 
     &                      - ALF1 * FY(IJ, KU) )

                  SX (IJ, K, N) = SX (IJ, K, N) + FX (IJ, KU)
                  SXX(IJ, K, N) = SXX(IJ, K, N) + FXX(IJ, KU)

                  SY (IJ, K, N) = SY (IJ, K, N) + FY (IJ, KU)
                  SYY(IJ, K, N) = SYY(IJ, K, N) + FYY(IJ, KU)

                  SXY(IJ, K, N) = SXY(IJ, K, N) + FXY(IJ, KU)

!---- flux from (K) to (K-1),  when W < 0
               ELSE

                  ALF1 = 1.D0 - ALF(IJ, K)

                  TMP = - ALF(IJ, K) * S0(IJ, KU, N) + ALF1 * F0(IJ, KU)

                  S0(IJ, KU, N) = S0(IJ, KU, N) + F0(IJ, KU)

                  SZZ(IJ, KU, N) = 
     &                 ALF(IJ, K) * ALF(IJ, K) * FZZ(IJ, KU)
     &               + ALF1 * ALF1 * SZZ(IJ, KU, N)
     &               + 5.D0 * (  ALF(IJ, K) * ALF1 * ( - SZ(IJ, KU, N) 
     &                                                 + FZ(IJ, KU)    )
     &                         + ( ALF1 - ALF(IJ, K) ) * TMP  )
                  SZ (IJ, KU, N) = 
     &                 ALF(IJ, K) * FZ(IJ, KU) + ALF1 * SZ(IJ, KU, N) 
     &               + 3.D0 * TMP

                  SXZ(IJ, KU, N) = 
     &                 ALF(IJ, K) * FXZ(IJ, KU) + ALF1 * SXZ(IJ, KU, N)
     &               + 3.D0 * ( ALF1       * FX(IJ, KU) 
     &                        - ALF(IJ, K) * SX(IJ, KU, N) )
                  SYZ(IJ, KU, N) = 
     &                 ALF(IJ, K) * FYZ(IJ, KU) + ALF1 * SYZ(IJ, KU, N)
     &               + 3.D0 * ( ALF1       * FY(IJ, KU) 
     &                        - ALF(IJ, K) * SY(IJ, KU, N) )

                  SX (IJ, KU, N) = SX (IJ, KU, N) + FX (IJ, KU)
                  SY (IJ, KU, N) = SY (IJ, KU, N) + FY (IJ, KU)
                  SXX(IJ, KU, N) = SXX(IJ, KU, N) + FXX(IJ, KU)
                  SYY(IJ, KU, N) = SYY(IJ, KU, N) + FYY(IJ, KU)
                  SXY(IJ, KU, N) = SXY(IJ, KU, N) + FXY(IJ, KU)

               END IF

            END DO
         END DO

      END DO

!---- TX
      DO N = 1, NTDIM
         DO K = KSTR, KEND
            DO IJ = IJTSTR, IJTEND

               TX(IJ, K, N) = S0(IJ, K, N) / SM(IJ, K, N)

            END DO
         END DO
      END DO

#ifdef OPT_BBL
!---- keeping BBL variables consistent at two levels
      CALL STBBTR( S0 )
      CALL STBBTR( SM )
      CALL STBBTR( SX )
      CALL STBBTR( SY )
      CALL STBBTR( SZ )
      CALL STBBTR( SXX )
      CALL STBBTR( SYY )
      CALL STBBTR( SZZ )
      CALL STBBTR( SXY )
      CALL STBBTR( SXZ )
      CALL STBBTR( SYZ )
#endif

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT1(    SX,
     D             NXDIM, NYDIM, NZTDIM,
     I              1.D0,     0,      0 )
      CALL SHIFT1(    SY,
     D             NXDIM, NYDIM, NZTDIM,
     I             -1.D0,     0,      0 )
      CALL SHIFT1(    SZ,
     D             NXDIM, NYDIM, NZTDIM,
     I              1.D0,     0,      0 )
      CALL SHIFT1(   SXX,
     D             NXDIM, NYDIM, NZTDIM,
     I              1.D0,     0,      0 )
      CALL SHIFT1(   SYY,
     D             NXDIM, NYDIM, NZTDIM,
     I              1.D0,     0,      0 )
      CALL SHIFT1(   SZZ,
     D             NXDIM, NYDIM, NZTDIM,
     I              1.D0,     0,      0 )
      CALL SHIFT1(   SXY,
     D             NXDIM, NYDIM, NZTDIM,
     I             -1.D0,     0,      0 )
      CALL SHIFT1(   SXZ,
     D             NXDIM, NYDIM, NZTDIM,
     I              1.D0,     0,      0 )
      CALL SHIFT1(   SYZ,
     D             NXDIM, NYDIM, NZTDIM,
     I             -1.D0,     0,      0 )
#else
!      CALL SHIFT1( S0 , NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SX , NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SY , NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SZ , NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SXX, NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SYY, NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SZZ, NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SXY, NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SXZ, NXDIM, NYDIM, NZTDIM )
      CALL SHIFT1( SYZ, NXDIM, NYDIM, NZTDIM )
#endif
#else
!      CALL STBCTR( S0 , R )
      CALL STBCTR( SX , R )
      CALL STBCTR( SY , R )
      CALL STBCTR( SZ , R )
      CALL STBCTR( SXX, R )
      CALL STBCTR( SYY, R )
      CALL STBCTR( SZZ, R )
      CALL STBCTR( SXY, R )
      CALL STBCTR( SXZ, R )
      CALL STBCTR( SYZ, R )
#endif

      RETURN
      END

c *********************************************************************

      SUBROUTINE DNSGRD(
     O              XDZDX,  YDZDY,  ZDZDX,  ZDZDY,
     O              XDTDZ,  YDTDZ,  ZDTDX,  ZDTDY,
     I                 TY,     TX,     
     W               DTDX,   DTDY,   DTDZ)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

      REAL*8   XDZDX(NXYDIM, NZDIM),  YDZDY(NXYDIM, NZDIM)
      REAL*8   ZDZDX(NXYDIM, NZDIM),  ZDZDY(NXYDIM, NZDIM)
      REAL*8   XDTDZ(NXYDIM, NZDIM, NTDIM),  YDTDZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   ZDTDX(NXYDIM, NZDIM, NTDIM),  ZDTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8    DTDX(NXYDIM, NZDIM, NTDIM),   DTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8    DTDZ(NXYDIM, NZDIM, NTDIM)

      REAL*8 C0(NZDIM), C1(NZDIM), C2(NZDIM)
      REAL*8 C3(NZDIM), C4(NZDIM), C5(NZDIM), C6(NZDIM)
      REAL*8 D0(NZDIM), D1(NZDIM), D2(NZDIM), D3(NZDIM), D4(NZDIM)
      REAL*8 D5(NZDIM), D6(NZDIM), D7(NZDIM), D8(NZDIM), D9(NZDIM)
      REAL*8 EPS
      LOGICAL OFIRST
      SAVE C0, C1, C2, C3, C4, C5, C6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9
      SAVE EPS, OFIRST

      REAL*8      P1,     P2
      REAL*8      TL,     SL
      REAL*8      RL,    RLW,    RLS,    RLU
      REAL*8    DZDX,   DZDY
      INTEGER     IJ,      K,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8  SLPMAX
      REAL*8  SLPZ(NZDIM), DZM0(NZDIM)
      INTEGER KSLP
      SAVE SLPMAX, SLPZ, KSLP
      NAMELIST /NMSLPM/ SLPMAX, KSLP

      DATA OFIRST / .TRUE. /
      DATA SLPMAX, KSLP / 1.D-2, -999 /

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSLPM, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMSLPM)

         CALL SECOEF(
     O            C0(KSTR), C1(KSTR), C2(KSTR), C3(KSTR),
     O            C4(KSTR), C5(KSTR), C6(KSTR),
     O            D0(KSTR), D1(KSTR), D2(KSTR), D3(KSTR),  D4(KSTR),
     O            D5(KSTR), D6(KSTR), D7(KSTR), D8(KSTR),  D9(KSTR))

         EPS = 1.D-20

         IF ( KSLP < 0 ) THEN

            DO K = 1, NZDIM
               SLPZ(K) = SLPMAX
            END DO

         ELSE

!---- vertically varying SLPMAX (like OPA)
            DO K = 1, KSTR-1
               DZM0(K) = 0.5D0 * DZ0(K)
            END DO
            DO K = KSTR, NZDIM
               DZM0(K) = 0.5D0 * (DZ0(K-1) + DZ0(K))
            END DO

            P1 = 0.D0
            DO K = KSTR, KSTR+KSLP-1
               P1 = P1 + DZM0(K)
            END DO

            P2 = 0.D0
            DO K = KSTR, KSTR+KSLP-1
               P2 = P2 + DZM0(K)
!            SLPZ(K) = SLPMAX * EXP( - 10.D0 * ( P1 - P2 ) / P1 )
               SLPZ(K) = ( SLPMAX - 0.D0 ) / ( P1 - 0.D0 ) * P2
            END DO

            DO K = 1, KSTR-1
               SLPZ(K) = SLPZ(KSTR)
            END DO

            DO K = KSTR+KSLP, NZDIM
               SLPZ(K) = SLPZ(KSTR+KSLP-1)
            END DO

            WRITE(JFPAR, *) 'depth dependent Maximum slope:'
            DO K = 1, NZDIM
               WRITE(JFPAR, '(I6,F10.3)') K, SLPZ(K)
            END DO

         END IF

      END IF

      DO 115 N = 1, NTDIM
         DO 110 K = 1, NZDIM
            DO 100 IJ = 1, NXYDIM
               DTDX(IJ, K, N) = 0.D0
               DTDY(IJ, K, N) = 0.D0
               DTDZ(IJ, K, N) = 0.D0
  100       CONTINUE
  110    CONTINUE
  115 CONTINUE

      DO 200 K = KSTR, KEND
         DO 210 IJ = IJTSTR, IJTEND+NXDIM
            TL = TY(IJ, K, 1)
            SL = TY(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RL = P1 / P2

            TL = TY(IJ+LW, K, 1)
            SL = TY(IJ+LW, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLW = P1 / P2

            TL = TY(IJ+LS, K, 1)
            SL = TY(IJ+LS, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLS = P1 / P2

            DTDX(IJ, K, 1) = (RL - RLW) * RX * 2.D0
     &                       / (HXT(IJ) + HXT(IJ+LW)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LW, K)
            DTDY(IJ, K, 1) = (RL - RLS) * RYM(IJ+LS) * 2.D0
     &                       / (HYT(IJ) + HYT(IJ+LS)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LS, K)
  210    CONTINUE
  200 CONTINUE

      DO 230 K = KSTR+1, KEND
         DO 240 IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            TL = TY(IJ, K, 1)
            SL = TY(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RL = P1 / P2

            TL = TY(IJ, K-1, 1)
            SL = TY(IJ, K-1, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLU = P1 / P2

            DTDZ(IJ, K, 1) = MIN((RLU - RL) / DZM(IJ, K), 0.D0) *
     &                       AMFTZ(IJ, K)
  240    CONTINUE
  230 CONTINUE

      DO 300 K = KSTR, KEND
         DO 310 IJ = IJTSTR, IJTEND+NXDIM
            DZDX = DTDX(IJ, K, 1) * 4.D0
     &             / (  DTDZ(IJ, K  , 1) + DTDZ(IJ+LW, K ,  1)
     &                + DTDZ(IJ, K+1, 1) + DTDZ(IJ+LW, K+1, 1) - EPS)
            XDZDX(IJ, K) = MIN(SLPZ(K), MAX(-SLPZ(K), DZDX))
            DZDY = DTDY(IJ, K, 1) * 4.D0 * DYM(IJ+LS)
     &          / (  (DTDZ(IJ, K, 1) + DTDZ(IJ, K+1, 1)) * DY(IJ)
     &          + (DTDZ(IJ+LS, K, 1) + DTDZ(IJ+LS, K+1, 1)) * DY(IJ+LS)
     &                - EPS)
            YDZDY(IJ, K) = MIN(SLPZ(K), MAX(-SLPZ(K), DZDY))
  310    CONTINUE
  300 CONTINUE

      DO 330 K = KSTR+1, KEND
         DO 340 IJ = IJTSTR, IJTEND
            DZDX = ((DTDX(IJ, K-1, 1)+DTDX(IJ+LE, K-1, 1))*DZ(IJ, K-1)
     &            + (DTDX(IJ, K  , 1)+DTDX(IJ+LE, K  , 1))*DZ(IJ, K)) *
     &             0.25D0 / DZM(IJ, K) / (DTDZ(IJ, K, 1) - EPS)
            ZDZDX(IJ, K) = MIN(SLPZ(K), MAX(-SLPZ(K), DZDX)) *
     &                     AMFTZ(IJ, K)
            DZDY = 
     &        (  (DTDY(IJ, K-1, 1)+DTDY(IJ+LN, K-1, 1))*DZ(IJ, K-1)
     &         + (DTDY(IJ, K  , 1)+DTDY(IJ+LN, K  , 1))*DZ(IJ, K)) *
     &          0.25D0 / DZM(IJ, K) / (DTDZ(IJ, K, 1) - EPS)
            ZDZDY(IJ, K) = MIN(SLPZ(K), MAX(-SLPZ(K), DZDY)) *
     &                     AMFTZ(IJ, K)
  340    CONTINUE
  330 CONTINUE

      DO 165 N = 1, NTDIM
         DO 160 K = 1, NZDIM
            DO 150 IJ = 1, NXYDIM
               DTDX(IJ, K, N) = 0.D0
               DTDY(IJ, K, N) = 0.D0
               DTDZ(IJ, K, N) = 0.D0
  150      CONTINUE
  160   CONTINUE
  165 CONTINUE

      DO 420 N = 1, NTDIM
         DO 410 K = KSTR, KEND
            DO 400 IJ = IJTSTR, IJTEND+NXDIM
            DTDX(IJ, K, N) = (TX(IJ, K, N) - TX(IJ+LW, K, N)) * RX *
     &                       2.D0 / (HXT(IJ) + HXT(IJ+LW)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LW, K)
            DTDY(IJ, K, N) = (TX(IJ, K, N) - TX(IJ+LS, K, N)) *
     &                       RYM(IJ+LS) *
     &                       2.D0 / (HYT(IJ) + HYT(IJ+LS)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LS, K)
  400       CONTINUE
  410    CONTINUE
  420 CONTINUE

      DO 450 N = 1, NTDIM
         DO 440 K = KSTR+1, KEND
            DO 430 IJ = IJTSTR-NXDIM, IJTEND+NXDIM
               DTDZ(IJ, K, N) = (TX(IJ, K-1, N) - TX(IJ, K, N)) 
     &                          / DZM(IJ, K) * AMFTZ(IJ, K)
  430       CONTINUE
  440    CONTINUE
  450 CONTINUE

      DO 520 N = 1, NTDIM
         DO 510 K = KSTR, KEND
            DO 500 IJ = IJTSTR, IJTEND+NXDIM
               XDTDZ(IJ, K, N) =
     &         (  DTDZ(IJ, K, N)   + DTDZ(IJ+LW, K, N)
     &          + DTDZ(IJ, K+1, N) + DTDZ(IJ+LW, K+1, N)) *  0.25D0
               YDTDZ(IJ, K, N) =
     &         (  (DTDZ(IJ, K, N)    + DTDZ(IJ, K+1, N)) *
     &            DY(IJ)
     &          + (DTDZ(IJ+LS, K, N) + DTDZ(IJ+LS, K+1, N)) *
     &            DY(IJ+LS)) *
     &         RYM(IJ+LS) * 0.25D0
  500       CONTINUE
  510    CONTINUE
  520 CONTINUE

      DO 550 N = 1, NTDIM
         DO 540 K = KSTR+1, KEND
            DO 530 IJ = IJTSTR, IJTEND
               ZDTDX(IJ, K, N) = 
     &              (  (DTDX(IJ, K-1, N) + DTDX(IJ+LE, K-1, N)) *
     &                 DZ(IJ, K-1)
     &               + (DTDX(IJ, K, N) + DTDX(IJ+LE, K, N)) *
     &                 DZ(IJ, K))
     &              / DZM(IJ, K) * 0.25D0 * AMFTZ(IJ, K)
               ZDTDY(IJ, K, N) =
     &              (  (DTDY(IJ, K-1, N) + DTDY(IJ+LN, K-1, N)) *
     &                 DZ(IJ, K-1)
     &               + (DTDY(IJ, K, N) + DTDY(IJ+LN, K, N)) *
     &                 DZ(IJ, K))
     &              / DZM(IJ, K) * 0.25D0 * AMFTZ(IJ, K)
  530       CONTINUE
  540    CONTINUE
  550 CONTINUE

      RETURN
      END

#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE FLXTRB(
     O                ADT,  DIFFZ,
     I                 TX,
     I                 TY,
     I                 UY,     VY,      W,
     I                AHV)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      INTEGER     IDUMMY

!---- only a dummy routine to preserve a consistency
      IDUMMY = 1

      RETURN
      END
#endif
c *********************************************************************
      SUBROUTINE CHKFTX
      use zocdim
      implicit none

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      COMMON /WORKFT/ FTX, FTY, FTZ

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      CALL CHEKIN(   FTX,  'FTX',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTY,  'FTY',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTZ,  'FTZ',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTX(1, 1, 2),  'FSX',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTY(1, 1, 2),  'FSY',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTZ(1, 1, 2),  'FSZ',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')

      RETURN
      END 
