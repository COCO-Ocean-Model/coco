      SUBROUTINE VDIFF(
     O               AMV,    AHV,
     I                UY,     VY,      R,   TAUX,   TAUY,
     I                TY,     HY)
      use zocdim
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Setting vertical diffusion and viscosity coefficients
c
c  HISTORY
c     '97.03.23  H.Hasumi
c     '98.05.08  H.Hasumi: revision of the consistency with the
c                          enhanced diffusivity version for the
c                          parameterization of the vertical convection
c     '98.05.11  H.Hasumi: namelist
c     '02.05.29  H.Nakano: tracer dimension
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       R(NXYDIM, NZDIM)
      REAL*8     AMV(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      HY(NXYDIM)

      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8    AMV0(NZ),   AHV0(NZ)
      SAVE      AMV0, AHV0
      NAMELIST /NMVISV/ AMV0
      NAMELIST /NMDIFV/ AHV0

      DATA OFIRST / .TRUE. /
      DATA AMV0 / NZ*0.D0 /
      DATA AHV0 / NZ*0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMVISV, END=99)
 99      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFV, END=98)
 98      CONTINUE
         WRITE(JFPAR, NMVISV)
         WRITE(JFPAR, NMDIFV)
      END IF

      DO 110 K = KSTR, KEND
         DO 100 IJ = 1, NXYDIM
            AMV(IJ, K) = AMV0(K-KSTR+1)
            AHV(IJ, K) = AHV0(K-KSTR+1)
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE VDIFFB(
     M                AMV,    AHV)
      use zocdim
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Vertical viscosity and diffusion coefficients for the bottom
c boundary layer.
c
c  HISTORY
c     '01.02.08  H.Hasumi
c
c ---------------------------------------------------------------------

      REAL*8     AMV(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      INTEGER     IJ,      K
      INTEGER     KT,     KV
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  AMVBBL, AHVBBL
      SAVE AMVBBL, AHVBBL
      NAMELIST /NMBBDV/ AMVBBL, AHVBBL

      DATA OFIRST / .TRUE. /
      DATA AMVBBL, AHVBBL / 0.D0, 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBBDV, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMBBDV)
      END IF

      DO IJ = IJSTR, IJEND
         KT = NBOT(IJ)
         KV = NBOTV(IJ)
         AMV(IJ, KV) = MAX(AMVBBL, AMV(IJ, KV)) * AMSKVB(IJ)
     &               + AMV(IJ, KV) * (1.D0 - AMSKVB(IJ))
         AHV(IJ, KT) = MAX(AHVBBL, AHV(IJ, KT)) * AMSKTB(IJ)
     &               + AHV(IJ, KT) * (1.D0 - AMSKTB(IJ))
      END DO

      DO IJ = IJSTR, IJEND
         KT = NBOT(IJ)
         KV = NBOTV(IJ)
         AMV(IJ, KEND) = AMV(IJ, KV)
         AHV(IJ, KEND) = AHV(IJ, KT)
      END DO

      RETURN
      END
#endif
