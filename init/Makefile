include ../src/Mkinclude
include Sysdep.$(SYSTEM)

initial$(PARALLEL):initial$(PARALLEL).o
	$(FC) $(FOPTION) -o $@ initial$(PARALLEL).o
	ln -s INITIAL $(UNIT41)
	initial$(PARALLEL) $(RUNOPT)
	$(RM) $(UNIT41)

clean:
	$(RM) *~ *.o *.trace *.cpp core
