      SUBROUTINE SRCVEL(
     O                 GX,     GY,     XX,     YY,
     I                 UX,     VX)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the bottom friction term of the equation of motion
c
c  HISTORY
c     '99.08.13  H.Hasumi: from CCSR2-MASK
c     '01.01.30  H.Hasumi: for partial step bottom topography
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)

      REAL*8     ABV
      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  BTMFRC
      SAVE BTMFRC
      NAMELIST /NMBTMF/ BTMFRC

      DATA OFIRST / .TRUE. /
      DATA BTMFRC / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBTMF, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMBTMF)
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            GX(IJ, K) = 0.D0
            GY(IJ, K) = 0.D0
            XX(IJ, K) = 0.D0
            YY(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE

      DO 210 K = KSTR, KEND
         DO 200 IJ = IJVSTR, IJVEND
            ABV = - BTMFRC / DZV(IJ, K) *
     &              SQRT(UX(IJ, K) * UX(IJ, K) + VX(IJ, K) * VX(IJ, K))
            GX(IJ, K) = ABV * UX(IJ, K) * AMSKB(IJ, K)
            GY(IJ, K) = ABV * VX(IJ, K) * AMSKB(IJ, K)
            XX(IJ, K) = GX(IJ, K)
            YY(IJ, K) = GY(IJ, K)
  200    CONTINUE
  210 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE SRCVLB(
     O                 GX,     GY,     XX,     YY,
     I                 UX,     VX)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Bottom friction term for the BBL momemtum equations. See Nakano and
c Suginohara (2002, JPO).
c
c  HISTORY
c     '01.02.08  H.Hasumi
c
c ---------------------------------------------------------------------

      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)

      REAL*8    RFRC(NXYDIM)
      SAVE RFRC

      REAL*8     ABV
      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  BTMFRC, RAYFRC
      INTEGER     MZ
      SAVE BTMFRC
      NAMELIST /NMBTMF/ BTMFRC
      NAMELIST /NMBBRF/ RAYFRC, MZ

      DATA OFIRST / .TRUE. /
      DATA BTMFRC / 0.D0 /
      DATA RAYFRC / 1.D0 /
      DATA MZ / NZ /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBTMF, END=99)
 99      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBBRF, END=98)
 98      CONTINUE
         WRITE(JFPAR, NMBBRF)

         DO IJ = 1, NXYDIM
            IF (NBOTV(IJ) .LE. MZ+KSTR-1) THEN
               RFRC(IJ) = RAYFRC * ABS(COR(IJ))
            ELSE
               RFRC(IJ) = 0.D0
            END IF
         END DO
      END IF

      DO IJ = IJVSTR, IJVEND
         ABV = - BTMFRC / DZV(IJ, KEND) *
     &           SQRT(  UX(IJ, KEND) * UX(IJ, KEND)
     &                + VX(IJ, KEND) * VX(IJ, KEND))
     &         - RFRC(IJ)
         GX(IJ, KEND) = ABV * UX(IJ, KEND) * AMSKVB(IJ)
         GY(IJ, KEND) = ABV * VX(IJ, KEND) * AMSKVB(IJ)
         XX(IJ, KEND) = GX(IJ, KEND)
         YY(IJ, KEND) = GY(IJ, KEND)
      END DO

      RETURN
      END
#endif
