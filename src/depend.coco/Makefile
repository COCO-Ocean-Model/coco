include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJSTD	= \
          aocea.o \
          aprdc$(TSTAG).o \
          bfrch.o \
          binst.o \
          brstt.o \
          bstbc.o \
          sfcng.o \
          tslvt.o \
          utint.o

$(LIBOGCM):$(OBJSTD)
	$(AR) $(LIBOGCM) $(OBJSTD)
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.o *.trace *.cpp core
