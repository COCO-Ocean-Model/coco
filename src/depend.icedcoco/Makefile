include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)

OBJSTD1	= \
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

OBJSTD2	= \
          aocea.o \
          aprdc$(TSTAG).o \

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
	$(RM) *~ *.o *.mod *.trace *.cpp core
