include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJS.coco	= \
          ufile.o \
          ucaln.o \
          bgs2d.o \
          bgs3d.o \
          bshft.o \
          bshfi.o \
          bstbc.o \
          utrdg.o \
          utint.o \
          xprst.o \
          atmct.o \
          bchmk.o \
          brdge.o \
          brstt.o \
          cadvc.o \
          csrvl.o \
          ctnuv.o \
          cvisc.o \
          cvlrd.o \
          dvdif.o \
          dvlta.o \
          dvlva.o \
          dwdns.o \
          fbtav.o \
          fshdf.o \
          fshlw.o \
          qckot.o \
          tovtr.o \
          ucloc.o \
          bfrch.o \
          sfcng.o \
          tflxt.o \
          tslvt.o \
          binst.o \
          aprdc.o \
          aocea.o 

OBJS.icedcoco	= \
          ufile.o \
          ucaln.o \
          bgs2d.o \
          bgs3d.o \
          bshft.o \
          bshfi.o \
          utrdg.o \
	  utint.o \
          xprst.o \
          atmct.o \
          bchmk.o \
          bgsid.o \
          brdge.o \
          brstt.o \
          cadvc.o \
          csrvl.o \
          ctnuv.o \
          cvisc.o \
          cvlrd.o \
          dvdif.o \
          dvlta.o \
          dvlva.o \
          dwdns.o \
          fbtav.o \
          fshdf.o \
          fshlw.o \
          qckot.o \
          tovtr.o \
          ucloc.o \
          bfrch.o \
          bstbc.o \
          ictrn.o \
          ifhea.o \
          ifwtr.o \
          ipadv.o \
          ipmmt.o \
          iprdg.o \
          iptmp.o \
          ipthm.o \
          iprdc.o \
          sfcng.o \
          tflxt.o \
          tslvt.o \
          binst.o \
	  aprdc.o \
          aocea.o 

$(LIBOGCM):$(OBJS.$(MODEL))
	$(AR) $(LIBOGCM) $(OBJS.$(MODEL))
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.F90 *.F *.o *.mod *.trace *.cpp core
