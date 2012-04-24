DSET ^ocn.grd
TITLE 2-deg interpolated ocean data
OPTIONS template big_endian
UNDEF -999.
XDEF 180 LINEAR 1. 2.
YDEF 90 LINEAR -89. 2.
ZDEF 40 LEVELS
2.5 7.5 14. 24. 40. 65. 100. 145. 200. 265.
340. 425. 520. 625. 740. 865. 1000. 1145. 1300. 1465.
1640. 1822.5 2010. 2202.5 2400. 2600. 2800. 3000. 3200. 3400.
3600. 3800. 4000. 4200. 4400. 4600. 4800. 5000. 5200. 5300.
TDEF 1 LINEAR  00:00Z1JAN0000 1yr
VARS 11
t 40 99 temperature
s 40 99 salinity
u 40 99 zonal velocity
v 40 99 meridional velocity
w 40 99 vertical velocity
c 40 99 convective adjustment index
sh 0 99 sea level
ai 0 99 ice concentration
hi 0 99 ice thickness
ui 0 99 ice zonal velocity
vi 0 99 ice meridional velocity
ENDVARS
