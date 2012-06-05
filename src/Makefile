include Mkinclude
include Sysdep.$(SYSTEM)$(PARALLEL)
MAIN	= $(MODEL).o
OPTLST	= OPTION-LIST
ENSLST	= ENSTR-LIST
MODOPT	= $(shell echo $(OBJOPT) | sed -e "s/\.[^[:blank:]]*o/\.mod/g")
MODENS	= $(shell echo $(OBJENS) | sed -e "s/\.[^[:blank:]]*o/\.mod/g")

all:lib $(TARGET)

$(TARGET):$(LIBOGCM) $(MAIN)
	$(FC) $(FOPTION) -o $@ $(MAIN) $(OBJENS) $(OBJOPT) $(LIBOGCM) $(LINKOPT)
	$(RM) $(ENSLST) $(OBJENS) $(OPTLST) $(OBJOPT) $(MODOPT) $(MODENS)

$(ENSLST):
	cd enstr; make

$(OPTLST):
	cd option; make

lib:
	cd module; make
	cd depend.$(MODEL); make lib1
	cd shared; make lib1
	cd enstr; make lib1
	cd option; make lib1
	cd depend.$(MODEL); make lib2
	cd shared; make lib2
	cd enstr; make lib2
	cd option; make lib2
	cd enstr; make
	cd option; make

clean:
	$(RM) *~ *.o *.mod *.trace core
	cd module; make clean
	cd depend.$(MODEL); make clean
	cd shared; make clean
	cd enstr; make clean
	cd option; make clean
	$(RM) ../lib/*
