include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)
OPTLST	= OPTION-LIST

$(OPTLST):$(OBJOPT)
	echo $(OBJOPT) > $(OPTLST)
	cp $(OBJOPT) $(OPTLST) ..
	$(RM) $(OPTLST)

lib1:$(OBJOPT1)
	cp $(OBJOPT1) ..

lib2:$(OBJOPT2)
	cp $(OBJOPT2) ..

optlist:
	echo $(OBJOPT) > ../$(OPTLST)

clean:
	$(RM) *~ *.o *.trace *.cpp core
