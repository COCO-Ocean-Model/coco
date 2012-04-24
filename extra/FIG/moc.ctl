DSET ^moc.grd
TITLE meridional overturning streamfunction
OPTIONS template big_endian
UNDEF -999.
XDEF 1 LINEAR 0. 1.
YDEF 91 LINEAR -90. 2.
ZDEF 41 LEVELS
0. 5. 10. 18. 30. 50. 80. 120. 170. 230. 300.
380. 470. 570. 680. 800. 930. 1070. 1220. 1380. 1550.
1730. 1915. 2105. 2300. 2500. 2700. 2900. 3100. 3300. 3500.
3700. 3900. 4100. 4300. 4500. 4700. 4900. 5100. 5300. 5400.
TDEF 1 LINEAR 0Z1JAN1000 1yr
VARS 5
mocg 41 99 global meridional overturn [Sv]
moca 41 99 Atlantic meridional overturn [Sv]
mocp 41 99 Pacific meridional overturn [Sv]
moci 41 99 Indian meridional overturn [Sv]
moco 41 99 Indo-Pacific meridional overturn [Sv]
ENDVARS
