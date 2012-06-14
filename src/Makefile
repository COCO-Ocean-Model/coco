include Mkinclude
include Sysdep.$(SYSTEM)$(PARALLEL)
MAIN	= $(MODEL).o

SRCSHR	= \
          atmct$(TSTAG).F \
          bchmk.F \
          bgs2d.F \
          bgs3d.F \
          brdge.F \
          bshfi.F \
          bshft.F \
          cadvc$(TSTAG).F \
          csrvl.F \
          ctnuv$(TSTAG).F90 \
          cvisc.F \
          cvlrd.F \
          dvdif.F \
          dvlta.F \
          dvlva.F \
          dwdns.F \
          fbtav.F \
          fshdf.F \
          fshlw.F \
          qckot.F \
          tflxt.F \
          tovtr.F90 \
          ucaln.F \
          ucloc.$(SYSTEM).F \
          ufile.F \
          utrdg.F90 \
          xprst.F

SRCDPD.coco	= \
          aocea.F \
          aprdc$(TSTAG).F \
          bfrch.F \
          binst.F \
          brstt.F \
          bstbc.F \
          sfcng.F \
          tslvt.F \
          utint.F

SRCDPD.icedcoco	= \
          aocea.F \
          aprdc$(TSTAG).F \
          bfrch.F \
          binst.F \
          bgsid.F \
          brstt.F \
          bstbc.F \
          ictrn.F \
          ifhea.F \
          ifwtr.F \
          ipadv.F \
          ipmmt.F \
          iprdc$(TSTAG).F \
          iprdg.F \
          ipthm.F \
          iptmp.F \
          sfcng.F \
          tslvt.F \
          utint.F

all:lib $(TARGET)

$(TARGET):$(LIBOGCM) $(MAIN)
	$(FC) $(FOPTION) -o $@ $(MAIN) $(LIBOGCM) $(LINKOPT)

lib:
	$(RM) build/*.F90 build/*.F
	for esrc in $(SRCSHR); do \
	  esrcm=`echo $$esrc | sed -e "s/\..*F/\.F/g"` ; \
	  (cd build ; ln -s ../shared/$$esrc $$esrcm) ; \
	done
	for esrc in $(SRCDPD.$(MODEL)); do \
	  esrcm=`echo $$esrc | sed -e "s/\..*F/\.F/g"` ; \
	  (cd build ; ln -s ../depend.$(MODEL)/$$esrc $$esrcm) ; \
	done
	for esrc in $(SRCENS); do \
	  esrcm=`echo $$esrc | sed -e "s/\..*F/\.F/g"` ; \
	  esrcb=`echo $$esrc | sed -e "s/\..*$///"` ; \
	  (cd build ; $(RM) $$esrcb.F* ; ln -s ../enstr/$$esrc $$esrcm) ; \
	done
	for esrc in $(SRCOPT); do \
	  esrcm=`echo $$esrc | sed -e "s/\..*F/\.F/g"` ; \
	  esrcb=`echo $$esrc | sed -e "s/\..*$///"` ; \
	  (cd build ; $(RM) $$esrcb.F* ; ln -s ../option/$$esrc $$esrcm) ; \
	done
	cd module; make
	cd build; make

clean:
	$(RM) *~ *.o *.mod *.trace core
	cd module; make clean
	cd build; make clean
	$(RM) ../lib/*
