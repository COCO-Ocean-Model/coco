c *** Write here the values for NX, NY, NZ, and KZ as you specified
c *** in the geometry file
c === the parameters for non-tripole (bipole) ===
      INTEGER   NX, NY, NZ
      PARAMETER(NX = 32, NY = 32, NZ = 28)
      INTEGER   KZ
      PARAMETER(KZ = 5)
c === end of the parameters for non-tripole (bipole) ===
c === the parameters for tripole ===
c      INTEGER   NX, NY, NY1, NY2,  NZ
c      PARAMETER(NX = 360, NY1 = 153, NY2=39,  NZ = 50)
c      PARAMETER(NY = NY1 + NY2)
c      INTEGER   KZ
c      PARAMETER(KZ = 8)
c === end of the parameters for tripole ===

c *** Write here the number of thickness category for sea ice.
c *** The minimum for NIC is 1, which corresponds to two-category.
      INTEGER   NIC
      PARAMETER(NIC = 5)

c *** Write here the number of tracer.
c *** The minimum number is 2, for temperature and salinity.
      INTEGER   NTDIM
      PARAMETER(NTDIM = 2)
