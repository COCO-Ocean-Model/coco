include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJS.coco	= \
          utrdg.o \
          xprst.o \
          atmct.o \
          bchmk.o \
          bgs2d.o \
          bgs3d.o \
          brdge.o \
          bshfi.o \
          bshft.o \
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
          tflxt.o \
          tovtr.o \
          ucaln.o \
          ucloc.o \
          ufile.o \
          bfrch.o \
          brstt.o \
          bstbc.o \
          sfcng.o \
          tslvt.o \
          utint.o \
          binst.o \
          aocea.o \
          aprdc.o

OBJS.icedcoco	= \
          utrdg.o \
          xprst.o \
          atmct.o \
          bchmk.o \
          bgs2d.o \
          bgs3d.o \
          brdge.o \
          bshfi.o \
          bshft.o \
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
          tflxt.o \
          tovtr.o \
          ucaln.o \
          ucloc.o \
          ufile.o \
          bfrch.o \
          bgsid.o \
          brstt.o \
          bstbc.o \
          ictrn.o \
          ifhea.o \
          ifwtr.o \
          ipadv.o \
          ipmmt.o \
          iprdc.o \
          iprdg.o \
          ipthm.o \
          iptmp.o \
          sfcng.o \
          tslvt.o \
          utint.o \
          binst.o \
          aocea.o \
          aprdc.o

$(LIBOGCM):$(OBJS.$(MODEL))
	$(AR) $(LIBOGCM) $(OBJS.$(MODEL))
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.F90 *.F *.o *.mod *.trace *.cpp core
