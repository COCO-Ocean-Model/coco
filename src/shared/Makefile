include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJSTD	= \
          atmct$(TSTAG).o \
          bchmk.o \
          bgs2d.o \
          bgs3d.o \
          brdge.o \
          bshfi.o \
          bshft.o \
          cadvc$(TSTAG).o \
          csrvl.o \
          ctnuv$(TSTAG).o \
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
          ucloc.$(SYSTEM).o \
          ufile.o \
          utrdg.o \
          xprst.o

$(LIBOGCM):$(OBJSTD)
	$(AR) $(LIBOGCM) $(OBJSTD)
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.o *.trace *.cpp core
