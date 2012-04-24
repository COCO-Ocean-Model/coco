include ../src/Mkinclude
include Sysdep.$(SYSTEM)

bound$(PARALLEL):bound$(PARALLEL).o
	$(FC) $(FOPTION) -o $@ bound$(PARALLEL).o
	ln -s TSFC $(UNIT31)
	ln -s TDMP $(UNIT32)
	ln -s SSFC $(UNIT33)
	ln -s SDMP $(UNIT34)
	ln -s TAUX $(UNIT35)
	ln -s TAUY $(UNIT36)
	bound$(PARALLEL) $(RUNOPT)
	$(RM) $(UNIT31) $(UNIT32) $(UNIT33) $(UNIT34) $(UNIT35) $(UNIT36)

clean:
	$(RM) *~ *.o *.trace *.cpp core
