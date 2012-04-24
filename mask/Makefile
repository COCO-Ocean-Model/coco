include ../src/Mkinclude
include Sysdep.$(SYSTEM)$(PARALLEL)

GEO	= GEO
OUT	= zocdim.F

xmkmsk:
ifeq ($(OPTTRIPOLE), -DOPT_TRIPOLE)
	$(FC) $(FOPTION) -o $@ xmkmsk.tripole.F
else
	$(FC) $(FOPTION) -o $@ xmkmsk.F
endif
	ln -s MASK $(UNIT22)
	./xmkmsk $(RUNOPT) < $(GEO) > $(OUT)
	$(RM) $(UNIT22)

clean:
	$(RM) xmkmsk *.o *~ core *.trace *.cpp
