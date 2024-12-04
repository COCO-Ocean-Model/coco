include ../src/Mkinclude
include Sysdep.$(SYSTEM)

GEO	= GEO
OUT	= zocdim.F90

xmkmsk:
	$(FC) $(FOPTION) $(OPTEXMASK) $(OPTTRIPOLE) -o $@ xmkmsk.F
	./xmkmsk $(RUNOPT) < $(GEO) > $(OUT)

clean:
	$(RM) xmkmsk *.o *~ core *.trace *.cpp
