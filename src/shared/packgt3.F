************************************************************************
      SUBROUTINE GFWRT8
     I         ( GDATA ,
     I           JFILE , IJKDIM )
*
      IMPLICIT NONE
*
*   [PARAM]
      INTEGER    IJKDIM
*
*   [INPUT]
      REAL*8     GDATA ( IJKDIM )              !! data
      INTEGER    JFILE                         !! lun
*
#ifdef OPT_C_IO
*   [INTERNAL WORK]
      INTEGER    IOS
*
*   [EXTERNAL FUNC]
      INTEGER    XFWTR8
#endif
*
#ifdef OPT_C_IO
      IOS = XFWTR8( JFILE, GDATA, IJKDIM )
#else
      WRITE ( JFILE ) GDATA
#endif
*
      RETURN
      END
************************************************************************
      SUBROUTINE GFWRT4
     I         ( GDATA ,
     I           JFILE , IJKDIM )
*
      IMPLICIT NONE
*
*   [PARAM]
      INTEGER    IJKDIM
*
*   [INPUT]
      REAL*8     GDATA ( IJKDIM )              !! data
      INTEGER    JFILE                         !! lun
*
*   [INTERNAL WORK]
!      REAL*4     SDATA ( IJKDIM )
      REAL*4,  ALLOCATABLE  ::  SDATA(:)
*
#ifdef OPT_C_IO
      INTEGER    IOS
*
*   [EXTERNAL FUNC]
      INTEGER    XFWTR
#endif
*
      ALLOCATE ( SDATA(1:IJKDIM) )
      SDATA(1:IJKDIM) = SNGL( GDATA(1:IJKDIM ) )
*
#ifdef OPT_C_IO
      IOS = XFWTR( JFILE, SDATA, IJKDIM )
#else
      WRITE ( JFILE ) SDATA
#endif
*
      DEALLOCATE ( SDATA )
      RETURN
      END
************************************************************************
      SUBROUTINE GFWRM4         !! write MR4 (masked float)
     I         ( GDATA ,
     I           JFILE , IJKNUM ,
     I           VMISS             )
      IMPLICIT NONE
*
*   [INPUT]
      INTEGER    IJKNUM
      REAL*8     GDATA( IJKNUM )               !! Data
      REAL*8     VMISS
      INTEGER    JFILE                         !! unit number
*
*   [INTERNAL WORK]
      INTEGER*4  NNN                           !! Number Non NaN
!      REAL*4     GMSKD( IJKNUM )
!      INTEGER*4  IMSK ( IJKNUM )
!      INTEGER    IFLAG( IJKNUM )
      REAL(4),     ALLOCATABLE  ::  GMSKD(:)
      INTEGER(4),  ALLOCATABLE  ::   IMSK(:),  IFLAG(:)
      INTEGER    IJ, NMSK
*
#ifdef OPT_C_IO
      INTEGER IOS
      INTEGER XFWTR
#endif
*
      ALLOCATE ( IMSK (1:IJKNUM) )
      ALLOCATE ( GMSKD(1:IJKNUM) )
      ALLOCATE ( IFLAG(1:IJKNUM) )
*
      NNN = 0
      DO IJ = 1, IJKNUM
         IF ( GDATA( IJ ) .NE. VMISS ) THEN
            NNN = NNN + 1
            GMSKD( NNN ) = REAL( GDATA( IJ ) )
            IFLAG( IJ ) = 1
         ELSE
            IFLAG( IJ ) = 0
         ENDIF
      END DO
      CALL PACK_BITS_INTO32  !! 1-bit pack
     O     ( IMSK  , NMSK   ,
     I       IFLAG , IJKNUM , 1   )
*
#ifdef OPT_C_IO
      IOS = XFWTR(  JFILE, NNN,   1    )
      IOS = XFWTR(  JFILE, IMSK,  NMSK )
      IOS = XFWTR(  JFILE, GMSKD, NNN  )
#else
      WRITE ( JFILE ) NNN
      WRITE ( JFILE ) (IMSK(IJ),  IJ=1, NMSK)
      WRITE ( JFILE ) (GMSKD(IJ), IJ=1, NNN)
#endif
*
      DEALLOCATE ( IMSK, GMSKD, IFLAG )
      RETURN
      END
************************************************************************
      SUBROUTINE GFWRM8         !! write MR8 (masked double)
     I         ( GDATA ,
     I           JFILE , IJKNUM ,
     I           VMISS             )
      IMPLICIT NONE
*
*   [INPUT]
      INTEGER    IJKNUM
      REAL*8     GDATA( IJKNUM )               !! Data
      REAL*8     VMISS
      INTEGER    JFILE                         !! unit number
*
*   [INTERNAL WORK]
      INTEGER*4  NNN                           !! Number Non NaN
!      REAL*8     GMSKD( IJKNUM )
!      INTEGER*4  IMSK ( IJKNUM )
!      INTEGER    IFLAG( IJKNUM )
      REAL(8),     ALLOCATABLE  ::  GMSKD(:)
      INTEGER(4),  ALLOCATABLE  ::   IMSK(:),  IFLAG(:)
      INTEGER    IJ, NMSK
*
#ifdef OPT_C_IO
      INTEGER IOS
      INTEGER XFWTR, XFWTR8
#endif
*
      ALLOCATE ( IMSK (1:IJKNUM) )
      ALLOCATE ( GMSKD(1:IJKNUM) )
      ALLOCATE ( IFLAG(1:IJKNUM) )

      NNN = 0
      DO IJ = 1, IJKNUM
         IF ( GDATA( IJ ) .NE. VMISS ) THEN
            NNN = NNN + 1
            GMSKD( NNN ) = GDATA( IJ )
            IFLAG( IJ ) = 1
         ELSE
            IFLAG( IJ ) = 0
         ENDIF
      END DO
      CALL PACK_BITS_INTO32  !! 1-bit pack
     O   ( IMSK  , NMSK   ,
     I     IFLAG , IJKNUM , 1   )
*
#ifdef OPT_C_IO
      IOS = XFWTR(  JFILE, NNN,   1    )
      IOS = XFWTR(  JFILE, IMSK,  NMSK )
      IOS = XFWTR8( JFILE, GMSKD, NNN  )
#else
      WRITE ( JFILE ) NNN
      WRITE ( JFILE ) (IMSK(IJ),  IJ=1, NMSK)
      WRITE ( JFILE ) (GMSKD(IJ), IJ=1, NNN)
#endif
*
      DEALLOCATE ( IMSK, GMSKD, IFLAG )
      RETURN
      END
************************************************************************
      SUBROUTINE PACK_BITS_INTO32  !! pack N-bit integers
     O         ( PACKED , ILEN     ,
     I           IDATA  , LEN_ORIG , NBIT )
      IMPLICIT NONE
*
*   [OUTPUT]
      INTEGER*4 PACKED( 0:* )   !! packed array
      INTEGER   ILEN            !! length of packed array
*
*   [INPUT]
      INTEGER   IDATA( 0:* )    !! input(orignal) data
      INTEGER   LEN_ORIG        !! length of input data
      INTEGER   NBIT            !! the number of bits (between 1 and 31)

      INTEGER   BWIDTH
      PARAMETER (BWIDTH = 32)
*
#ifdef SYS_SX
      IF ( NBIT .GE. 8 ) THEN
         CALL PACKB_8           !! support 8 <= NBIT < 32
     O      ( PACKED , ILEN     ,
     I        IDATA  , LEN_ORIG , NBIT )
      ELSE
         CALL PACKB_1           !! support 1 <= NBIT < 32
     O      ( PACKED , ILEN     ,
     I        IDATA  , LEN_ORIG , NBIT )
      ENDIF
#else /* !SYS_SX */
      INTEGER MASK, I, IOFF, IP, IVAL
      INTEGER IPACK32LEN
*
      MASK = ISHFT( 1, NBIT ) - 1
      ILEN = IPACK32LEN( LEN_ORIG, NBIT )
*
      DO I = 0, ILEN - 1
         PACKED( I ) = 0
      END DO
*
      IP = 0
      IOFF = 0
      DO I = 0, LEN_ORIG - 1
         IF ( IOFF .GT. BWIDTH ) THEN
            IOFF = IOFF - BWIDTH
            IP = IP + 1
         ENDIF
*
         IVAL = IAND( IDATA( I ), MASK )
         PACKED( IP ) = IOR(
     &        PACKED( IP ),
     &        ISHFT( IVAL, BWIDTH - IOFF - NBIT ) )
         IF ( IOFF + NBIT .GT. BWIDTH ) THEN
            PACKED( IP+1 ) = IOR(
     &           PACKED( IP+1 ),
     &           ISHFT( IVAL, 2 * BWIDTH - IOFF - NBIT ) )
         ENDIF
*
         IOFF = IOFF + NBIT
      END DO
#endif /* !SYS_SX */
*
      RETURN
      END
************************************************************************
      FUNCTION IPACK32LEN( LEN_ORIG, NBIT )
      IMPLICIT NONE
      INTEGER IPACK32LEN
*
*   [INPUT]
      INTEGER LEN_ORIG          !! length of input data
      INTEGER NBIT              !! the number of bits (between 1 and 31)
*
      INTEGER BWIDTH
      PARAMETER (BWIDTH = 32 )
*
      !! (NBIT * LEN_ORIG + BWIDTH - 1) / BWIDTH
      IPACK32LEN = NBIT * (LEN_ORIG / BWIDTH)
     &     + ( NBIT * MOD(LEN_ORIG, BWIDTH) + BWIDTH - 1) / BWIDTH
      END

