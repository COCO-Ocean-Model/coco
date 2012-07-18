! Sea ice thermal energy
! The following use statement is needed before this
! use zocphy, only:    cpi,   dtds,   hfus 

  real(8) ::     ei
  real(8) ::   tice,   sice

  real(8) ::     ti
  real(8) ::   eice

  ei(tice, sice) = cpi * (dtds * sice - tice) &
    &            + hfus * (1.d0 - dtds * sice / tice)

  ti(eice, sice) = (  hfus + cpi * dtds * sice - eice &
    &               - sqrt(  (hfus + cpi * dtds * sice - eice)**2 &
    &                      - 4.d0 * cpi * hfus * dtds * sice)) &
    &              / cpi * 0.5d0
