include ../src/Mkinclude
include Sysdep.$(SYSTEM)

GEO	= GEO
OUT	= zocdim.F90

.PHONY: all xmkmsk zocdim clean

all: xmkmsk

zocdim: FOPTION += -DOPT_ZOCDIM
zocdim:	xmkmsk

xmkmsk:
	$(FC) $(FOPTION) $(OPTEXMASK) $(OPTTRIPOLE) $(OPTIO) -o $@ xmkmsk.F
	./xmkmsk $(RUNOPT) < $(GEO) > $(OUT)

clean:
	$(RM) xmkmsk *.o *~ core *.trace *.cpp
