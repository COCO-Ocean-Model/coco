      SUBROUTINE BTAVST(
     O                   UBTAV,  VBTAV)
      IMPLICIT NONE

c ---- information ----------------------------------------------------
c
c  HISTORY
c     '99.08.17  H.Hasumi
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"

      REAL*8   UBTAV(NXYDIM),  VBTAV(NXYDIM)
      REAL*8     UBT(NXYDIM),    VBT(NXYDIM)

      REAL*8    FAVR
      SAVE FAVR

      INTEGER     IJ

      DO 100 IJ = 1, NXYDIM
         UBTAV(IJ) = 0.D0
         VBTAV(IJ) = 0.D0
  100 CONTINUE

      FAVR = 1.D0 / DBLE(NTSS)

      RETURN

c =====================================================================

      ENTRY BTAVAD(
     M              UBTAV,  VBTAV,
     I                UBT,    VBT)

      DO 200 IJ = 1, NXYDIM
         UBTAV(IJ) = UBTAV(IJ) + FAVR * UBT(IJ)
         VBTAV(IJ) = VBTAV(IJ) + FAVR * VBT(IJ)
  200 CONTINUE

      RETURN
      END
