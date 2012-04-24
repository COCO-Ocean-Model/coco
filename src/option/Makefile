include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)
OPTLST	= OPTION-LIST

$(OPTLST):$(OBJOPT)
	echo $(OBJOPT) > $(OPTLST)
	cp $(OBJOPT) $(OPTLST) ..
	$(RM) $(OPTLST)

clean:
	$(RM) *~ *.o *.trace *.cpp core
