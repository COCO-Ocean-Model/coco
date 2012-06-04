      SUBROUTINE FORSTO(
     O                      UA,     VA,     TA,
     O                      HA,   UBTA,   VBTA,
     I                      UB,     VB,     TB,
     I                      HB,   UBTB,   VBTB)
      use zocdim
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '99.08.10  H.Hasumi
c     '02.05.29  H.Nakano: tracer dimension
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

      REAL*8      TA(NXYZDM, NTDIM),     TB(NXYZDM, NTDIM)
      REAL*8      UA(NXYZDM),     UB(NXYZDM)
      REAL*8      VA(NXYZDM),     VB(NXYZDM)
      REAL*8      HA(NXYDIM),     HB(NXYDIM)
      REAL*8    UBTA(NXYDIM),   UBTB(NXYDIM)
      REAL*8    VBTA(NXYDIM),   VBTB(NXYDIM)

      INTEGER     IJ,    IJK,      N

      DO 100 IJK = 1, NXYZDM
         UA(IJK) = UB(IJK)
         VA(IJK) = VB(IJK)
  100 CONTINUE
      DO 120 N = 1, NTDIM
         DO 110 IJK = 1, NXYZDM
            TA(IJK, N) = TB(IJK, N)
  110    CONTINUE
  120 CONTINUE
      DO 200 IJ  = 1, NXYDIM
         HA(IJ)   = HB(IJ)
         UBTA(IJ) = UBTB(IJ)
         VBTA(IJ) = VBTB(IJ)
  200 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE EXCNGO(
     M                      UA,     VA,     TA,
     M                      HA,   UBTA,   VBTA,
     M                      UB,     VB,     TB,
     M                      HB,   UBTB,   VBTB)
      use zocdim
      IMPLICIT NONE

      REAL*8      TA(NXYZDM, NTDIM),     TB(NXYZDM, NTDIM)
      REAL*8      UA(NXYZDM),     UB(NXYZDM)
      REAL*8      VA(NXYZDM),     VB(NXYZDM)
      REAL*8      HA(NXYDIM),     HB(NXYDIM)
      REAL*8    UBTA(NXYDIM),   UBTB(NXYDIM)
      REAL*8    VBTA(NXYDIM),   VBTB(NXYDIM)

      REAL*8       U,      V,      T
      REAL*8      SH,    UBT,    VBT
      INTEGER     IJ,    IJK,      N

      DO 100 IJK = 1, NXYZDM
         U       = UB(IJK)
         V       = VB(IJK)
         UB(IJK) = UA(IJK)
         VB(IJK) = VA(IJK)
         UA(IJK) = U
         VA(IJK) = V
  100 CONTINUE
      DO 120 N = 1, NTDIM
         DO 110 IJK = 1, NXYZDM
            T          = TB(IJK, N)
            TB(IJK, N) = TA(IJK, N)
            TA(IJK, N) = T
  110    CONTINUE
  120 CONTINUE      
      DO 200 IJ  = 1, NXYDIM
         SH       = HB(IJ)
         UBT      = UBTB(IJ)
         VBT      = VBTB(IJ)
         HB(IJ)   = HA(IJ)
         UBTB(IJ) = UBTA(IJ)
         VBTB(IJ) = VBTA(IJ)
         HA(IJ)   = SH
         UBTA(IJ) = UBT
         VBTA(IJ) = VBT
  200 CONTINUE

      RETURN
      END
