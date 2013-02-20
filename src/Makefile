include Mkinclude
include Sysdep.$(SYSTEM)
MAIN	= $(MODEL).o

SRCSHR	= \
          atmct$(TSTAG) \
          bchmk \
          bgs2d \
          bgs3d \
          brdge \
          bshfi \
          bshft \
          cadvc$(TSTAG) \
          csrvl \
          ctnuv$(TSTAG) \
          cvisc \
          cvlrd \
          dvdif \
          dvlta \
          dvlva \
          dwdns \
          fbtav \
          fshdf \
          fshlw \
          qckot \
          tflxt \
          tovtr \
          ucaln \
          ucloc.$(SYSTEM) \
          ufile \
          utrdg \
          xprst

SRCDPD.coco	= \
          aocea \
          aprdc$(TSTAG) \
          bfrch \
          binst \
          brstt \
          bstbc \
          sfcng \
          tslvt \
          utint

SRCDPD.icedcoco	= \
          aocea \
          aprdc$(TSTAG) \
          bfrch \
          binst \
          bgsid \
          brstt \
          bstbc \
          ictrn \
          ifhea \
          ifwtr \
          ipadv \
          ipmmt \
          iprdc$(TSTAG) \
          iprdg \
          ipthm \
          iptmp \
          sfcng \
          tslvt \
          utint

all:lib $(TARGET)

$(TARGET):$(LIBOGCM) $(MAIN)
	$(FC) $(FOPTION) -o $@ $(MAIN) $(LIBOGCM) $(LINKOPT)

lib:
	$(RM) build/*.F90 build/*.F
	for esrc in $(SRCSHR); do \
	  esrcb=`echo $$esrc | sed -e "s/\..*$///"` ; \
	  ( cd build ; \
	    if [ -f ../shared/$$esrc.F90 ]; then \
	     ln -s ../shared/$$esrc.F90 $$esrcb.F90 ; \
	    else \
	     ln -s ../shared/$$esrc.F $$esrcb.F ; \
	    fi ) ; \
	done
	for esrc in $(SRCDPD.$(MODEL)); do \
	  esrcb=`echo $$esrc | sed -e "s/\..*$///"` ; \
	  ( cd build ; \
	    if [ -f ../depend.$(MODEL)/$$esrc.F90 ]; then \
	     ln -s ../depend.$(MODEL)/$$esrc.F90 $$esrcb.F90 ; \
	    else \
	     ln -s ../depend.$(MODEL)/$$esrc.F $$esrcb.F ; \
	    fi ) ; \
	done
	for esrc in $(SRCENS); do \
	  esrcb=`echo $$esrc | sed -e "s/\..*$///"` ; \
	  (cd build ; $(RM) $$esrcb.F* ; \
	    if [ -f ../enstr/$$esrc.F90 ]; then \
	     ln -s ../enstr/$$esrc.F90 $$esrcb.F90 ; \
	    else \
	     ln -s ../enstr/$$esrc.F $$esrcb.F ; \
	    fi ) ; \
	done
	for esrc in $(SRCOPT); do \
	  esrcb=`echo $$esrc | sed -e "s/\..*$///"` ; \
	  (cd build ; $(RM) $$esrcb.F* ; \
	    if [ -f ../option/$$esrc.F90 ]; then \
	     ln -s ../option/$$esrc.F90 $$esrcb.F90 ; \
	    else \
	     ln -s ../option/$$esrc.F $$esrcb.F ; \
	    fi ) ; \
	done
	cd module; make
	cd build; make

clean:
	$(RM) *~ *.o *.mod *.trace core
	cd module; make clean
	cd build; make clean
	$(RM) ../lib/*
