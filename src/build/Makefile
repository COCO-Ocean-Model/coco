include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)

DIRETC = ../etc

SRCS_SHARED = \
	atmct.F90 \
	bchmk.F90 \
	bgs2d.F90 \
	bgs3d.F90 \
	brdge$(MPI).F90 \
	bshfi.F90 \
	bshft.F90 \
	bstbc.F90 \
	cadvc.F90 \
	csrvl.F90 \
	ctnuv.F90 \
	cvisc.F90 \
	cvlrd.F90 \
	dvdif.F90 \
	dvlta.F90 \
	dvlva.F90 \
	dwdns.F90 \
	fbtav.F90 \
	fshdf.F90 \
	fshlw.F90 \
	mpiio.F90 \
	packgt3.F \
	qckag.F90 \
	qckot$(MPI).F90 \
	tflxt.F90 \
	tovtr.F90 \
	ucaln.F90 \
	ucloc.F90 \
	ufile.F90 \
	utrdg.F90 \
	xprst.F90

SRCS_coco = \
	aocea.F90 \
	aprdc.F90 \
	bfrch.F90 \
	binst.F90 \
	brstt.F90 \
	sfcng.F90 \
	tslvt.F90 \
	utint.F90

SRCS_icedcoco = \
	aocea.F90 \
	aprdc.F90 \
	bfrch.F90 \
	bgsid$(MPI).F90 \
	binst.F90 \
	brstt$(MPI).F90 \
	ictrn.F90 \
	ifhea.F90 \
	ifwtr.F90 \
	ipadv.F90 \
	ipmmt.F90 \
	iprdc.F90 \
	iprdg.F90 \
	ipthm.F90 \
	iptmp.F90 \
	sfcng.F90 \
	tslvt.F90 \
	utint$(MPI).F90

SRCS_offline = \
	aocea.F90 \
	aprdc.F90 \
	brstt$(MPI).F90 \
	sfcng.F90 \
	tslvt.F90 \
	utint$(MPI).F90

SRCS_oeco2 =            \
	zoceco.F90      \
	dms_temp.F90    \
	n2osato.F90     \
	o2sato.F90      \
	uchem.F90       \
	tchem.F90       \
	utint$(MPI).F90 \
	sfcng.oeco.F90  \
	ecosys.F90      \
	brstt$(MPI).F90 \
	tslvt.F90       \
	aprdc.F90       \
	aocea.F90

## Options
SRCS_OPTION	= $(SRCOPT:%=%.F90)
SRCS_ENSTR	= $(SRCENS:%=%.F90)

## Sorted by lowest priority.
SRCS	= \
	$(SRCS_SHARED) \
	$(SRCS_$(MODEL)) \
	$(SRCS_ENSTR) \
	$(SRCS_OPTION)

DEPENDENCY	= makefile.dep

all: lib

lib: $(shell cat OBJLIST)
	$(AR) $(LIBOGCM) $(shell cat OBJLIST)
	$(RANLIB) $(LIBOGCM)

%.o : %.mod

.PHONY: FORCE
OBJLIST: FORCE
	-@$(RM) *.F90
	@sh $(DIRETC)/fetch.sh ../shared $(SRCS_SHARED)
	@sh $(DIRETC)/fetch.sh ../depend.$(MODEL) $(SRCS_$(MODEL))
	@sh $(DIRETC)/fetch.sh ../enstr $(SRCS_ENSTR)
	@sh $(DIRETC)/fetch.sh ../option $(SRCS_OPTION)
	@python3 $(DIRETC)/print_uniq_objs.py $(SRCS) > OBJLIST
	@python3 $(DIRETC)/fort_depend.py -w -o $(DEPENDENCY) \
		-f `cat OBJLIST | $(DIRETC)/srclist.sh` > /dev/null

clean:
	$(RM) *~ *.F90 *.F *.o *.mod *.trace *.L *.cpp core OBJLIST $(DEPENDENCY)

##
-include $(DEPENDENCY)
