include Mkinclude
include Sysdep.$(SYSTEM)$(PARALLEL)
MAIN	= $(MODEL).o
OPTLST	= OPTION-LIST
ENSLST	= ENSTR-LIST

all:lib $(TARGET)

$(TARGET):$(LIBOGCM) $(MAIN) $(ENSLST) $(OPTLST)
	$(FC) $(FOPTION) -o $@ $(MAIN) $(OBJENS) $(OBJOPT) $(LIBOGCM) $(LINKOPT)
	$(RM) $(ENSLST) $(OBJENS) $(OPTLST) $(OBJOPT)

$(ENSLST):
	cd enstr; make

$(OPTLST):
	cd option; make

lib:
	cd depend.$(MODEL); make
	cd shared; make

clean:
	$(RM) *~ *.o *.trace core
	cd depend.$(MODEL); make clean
	cd shared; make clean
	cd enstr; make clean
	cd option; make clean
	$(RM) ../lib/*
