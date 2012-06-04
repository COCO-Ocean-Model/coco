include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJMOD	= \
          zocnod.o \
          zocdim.o \
          zocfil.o \
          zocgrd.o \
          zocmsk.o \
          zocout.o \
          zocphy.o \
          zocite.o \
          zqsat.o

$(LIBOGCM):$(OBJMOD)
	$(AR) $(LIBOGCM) $(OBJMOD)
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.o *.mod *.trace *.cpp core
