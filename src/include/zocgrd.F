      REAL*8      DY(NXYDIM),    DYM(NXYDIM)
      REAL*8      DZ(NXYDIM, NZDIM) ,    DZM(NXYDIM, NZDIM)
      REAL*8     DZV(NXYDIM, NZDIM)
      REAL*8     DZ0(NZDIM),      DS(NZDIM) ,    DSM(NZDIM)
      REAL*8      DT,     DX
      REAL*8      RX
      REAL*8      RY(NXYDIM),    RYM(NXYDIM)
      REAL*8      RS(NZDIM),    RSM(NZDIM)
      REAL*8     HIC(NIC+1)
      COMMON /GRID/ DT, DX, DY, DYM, DZ, DZM, DZV, DZ0, DS, DSM,
     &              RX, RY, RYM, RS, RSM, HIC

      REAL*8      TT,     TS,    TSS
      INTEGER     NT,    ITS,   ITST,   NTSS
      INTEGER IEULER
      COMMON /TIMSTP/ TT, TS, TSS, NT, ITS, ITST, NTSS, IEULER

      REAL*8     REA,   ZBOT
      REAL*8     COR(NXYDIM),   DEPT(NXYDIM),  RDEPV(NXYDIM)
      REAL*8     HXT(NXYDIM),    HXU(NXYDIM)
      REAL*8     HYT(NXYDIM),    HYU(NXYDIM)
      REAL*8    HXYT(NXYDIM),   HXYU(NXYDIM)
      REAL*8    HYXT(NXYDIM),   HYXU(NXYDIM)
      REAL*8     RXT(NXYDIM),    RXU(NXYDIM)
      REAL*8     RYT(NXYDIM),    RYU(NXYDIM)
      COMMON /MTRC/ REA, ZBOT,
     &              COR, DEPT, RDEPV,
     &              HXT, HXU, HYT, HYU, HXYT, HXYU, HYXT, HYXU,
     &              RXT, RXU, RYT, RYU
