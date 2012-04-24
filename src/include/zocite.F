* Sea ice thermal energy
* You need to include zocphy.F before this
      REAL*8      EI
      REAL*8    TICE,   SICE

      REAL*8      TI
      REAL*8    EICE

      EI(TICE, SICE) = cpi * (dtds * SICE - TICE)
     &               + hfus * (1.D0 - dtds * SICE / TICE)

      TI(EICE, SICE) = (  hfus + cpi * dtds * SICE - EICE
     &                  - sqrt(  (hfus + cpi * dtds * SICE - EICE)**2
     &                         - 4.D0 * cpi * hfus * dtds * SICE))
     &                 / cpi * 0.5D0
