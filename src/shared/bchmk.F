#ifdef OPT_BBL
c --- information -----------------------------------------------------
c
c  Change the values of masking array for BBL grids.
c
c  HISTORY
c     '01.02.08  H.Hasumi
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

      SUBROUTINE RMMSKV
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKV(IJ, KEND) = 0.D0
         AMSKV(IJ, NBOTV(IJ)) = AMSKV0(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE ADMKV1
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKV(IJ, NBOTV(IJ)) = AMSKV1(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE ADMKVB
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKV(IJ, KEND) = AMSKVB(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE ADMSKV
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKV(IJ, KEND) = AMSKVB(IJ)
         AMSKV(IJ, NBOTV(IJ)) = AMSKV1(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE RMMSKT
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKT(IJ, KEND) = 0.D0
         AMSKT(IJ, NBOT(IJ)) = AMSKT0(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE ADMKTB
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKT(IJ, KEND) = AMSKTB(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE ADMKT1
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKT(IJ, NBOT(IJ)) = AMSKT1(IJ)
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE ADMSKT
      use zocdim
      use zocmsk

      INTEGER     IJ

      DO IJ = 1, NXYDIM
         AMSKT(IJ, KEND) = AMSKTB(IJ)
         AMSKT(IJ, NBOT(IJ)) = AMSKT1(IJ)
      END DO

      RETURN
      END
#else
c *********************************************************************
      SUBROUTINE ADMSKT
      RETURN
      END
#endif
