include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJSTD	= \
          aocea.o \
          aprdc$(TSTAG).o \
          bfrch.o \
          bgsid.o \
          binst.o \
          brstt.o \
          bstbc.o \
          ictrn.o \
          ifhea.o \
          ifwtr.o \
          ipadv.o \
          ipmmt.o \
          iprdc$(TSTAG).o \
          iprdg.o \
          ipthm.o \
          iptmp.o \
          sfcng.o \
          tslvt.o \
          utint.o

$(LIBOGCM):$(OBJSTD)
	$(AR) $(LIBOGCM) $(OBJSTD)
	$(RANLIB) $(LIBOGCM)

clean:
	$(RM) *~ *.o *.trace *.cpp core
