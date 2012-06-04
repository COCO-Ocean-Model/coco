include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJSTD1	= \
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

OBJSTD2	= 

OBJSTD	= $(OBJSTD1) \
          $(OBJSTD2)

$(LIBOGCM):$(OBJSTD)
	$(AR) $(LIBOGCM) $(OBJSTD)
	$(RANLIB) $(LIBOGCM)

lib1:$(OBJSTD1)
	$(AR) $(LIBOGCM) $(OBJSTD1)
	$(RANLIB) $(LIBOGCM)

lib2:$(OBJSTD2)
	$(AR) $(LIBOGCM) $(OBJSTD2)
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.o *.trace *.cpp core
