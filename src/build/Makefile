include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJS.coco	= \
          utrdg.o \
          utint.o \
          xprst.o \
          atmct.o \
          bchmk.o \
          bgs2d.o \
          bgs3d.o \
          brdge.o \
          bshfi.o \
          bshft.o \
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
          ucaln.o \
          ucloc.o \
          ufile.o \
          bfrch.o \
          bstbc.o \
          sfcng.o \
          tflxt.o \
          tslvt.o \
          binst.o \
          aocea.o \
          aprdc.o

OBJS.icedcoco	= \
          utrdg.o \
	  utint.o \
          xprst.o \
          atmct.o \
          bchmk.o \
          bgs2d.o \
          bgs3d.o \
          brdge.o \
          bshfi.o \
          bshft.o \
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
          ucaln.o \
          ucloc.o \
          ufile.o \
          bfrch.o \
          bgsid.o \
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
