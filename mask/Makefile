include ../src/Mkinclude
include Sysdep.$(SYSTEM)

GEO	= GEO
OUT	= zocdim.F90

xmkmsk:
	$(FC) $(FOPTION) $(OPTEXMASK) $(OPTTRIPOLE) -o $@ xmkmsk.F
	ln -s MASK $(UNIT22)
	./xmkmsk $(RUNOPT) < $(GEO) > $(OUT)
	$(RM) $(UNIT22)

clean:
	$(RM) xmkmsk *.o *~ core *.trace *.cpp $(UNIT22)
