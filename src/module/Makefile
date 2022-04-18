include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)

DIRETC = ../etc

DEPENDENCY = makefile.dep

OBJS	= \
	zocdim.o \
	zocnod.o \
	zocfil.o \
	zocgrd.o \
	zocmsk.o \
	zocout.o \
	zocphy.o \
	zocite.o \
	zqsat.o

SRCS	= $(OBJS:.o=.F90)

all lib: $(OBJS)
	$(AR) $(LIBOGCM) $(OBJS)
	$(RANLIB) $(LIBOGCM)

$(DEPENDENCY): $(SRCS)
	@python2 $(DIRETC)/fort_depend.py -w -o $@ -f $(SRCS)

clean:
	$(RM) *~ *.o *.mod *.trace *.cpp core $(DEPENDENCY)

##
-include $(DEPENDENCY)
