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
c     '02.10.10  H.Nakano: from COCO3.4
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

c *********************************************************************

      SUBROUTINE FORSTI(
     O                AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                AB,    HIB,    UIB,    VIB,    TIB,    HSB)
      use zocdim
      IMPLICIT NONE

      REAL*8     AA(NXYDIM, 0:NIC),     AB(NXYDIM, 0:NIC)
      REAL*8    HIA(NXYDIM, 0:NIC),    HIB(NXYDIM, 0:NIC)
      REAL*8    UIA(NXYDIM),    UIB(NXYDIM)
      REAL*8    VIA(NXYDIM),    VIB(NXYDIM)
      REAL*8    TIA(NXYDIM, 0:NIC),    TIB(NXYDIM, 0:NIC)
      REAL*8    HSA(NXYDIM, 0:NIC),    HSB(NXYDIM, 0:NIC)

      INTEGER     IJ,      K

      DO IJ = 1, NXYDIM
         UIA(IJ) = UIB(IJ)
         VIA(IJ) = VIB(IJ)
      END DO

      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            AA (IJ, K) = AB (IJ, K)
            HIA(IJ, K) = HIB(IJ, K)
            HSA(IJ, K) = HSB(IJ, K)
            TIA(IJ, K) = TIB(IJ, K)
         END DO
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE EXCNGI(
     M                AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     M                AB,    HIB,    UIB,    VIB,    TIB,    HSB)
      use zocdim
      IMPLICIT NONE
      
      REAL*8     AA(NXYDIM, 0:NIC),     AB(NXYDIM, 0:NIC)
      REAL*8    HIA(NXYDIM, 0:NIC),    HIB(NXYDIM, 0:NIC)
      REAL*8    UIA(NXYDIM),    UIB(NXYDIM)
      REAL*8    VIA(NXYDIM),    VIB(NXYDIM)
      REAL*8    TIA(NXYDIM, 0:NIC),    TIB(NXYDIM, 0:NIC)
      REAL*8    HSA(NXYDIM, 0:NIC),    HSB(NXYDIM, 0:NIC)

      INTEGER     IJ,      K
      REAL*8       A,     HI,     UI,     VI,     TI,     HS

      DO IJ = 1, NXYDIM
         UI      = UIB(IJ)
         VI      = VIB(IJ)
         UIB(IJ) = UIA(IJ)
         VIB(IJ) = VIA(IJ)
         UIA(IJ) = UI
         VIA(IJ) = VI
      END DO
      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            A          = AB (IJ, K)
            HI         = HIB(IJ, K)
            HS         = HSB(IJ, K)
            TI         = TIB(IJ, K)
            AB (IJ, K) = AA (IJ, K)
            HIB(IJ, K) = HIA(IJ, K)
            HSB(IJ, K) = HSA(IJ, K)
            TIB(IJ, K) = TIA(IJ, K)
            AA (IJ, K) = A
            HIA(IJ, K) = HI
            HSA(IJ, K) = HS
            TIA(IJ, K) = TI
         END DO
      END DO

      RETURN
      END
