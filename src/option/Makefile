include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)
OPTLST	= OPTION-LIST
#MODOPT	= $(OBJOPT:%.o=%.mod)
MODOPT1	= $(shell echo $(OBJOPT1) | sed -e "s/\.[^[:blank:]]*o/\.mod/g")
MODOPT2	= $(shell echo $(OBJOPT2) | sed -e "s/\.[^[:blank:]]*o/\.mod/g")

$(OPTLST):$(OBJOPT)
	echo $(OBJOPT) > $(OPTLST)
	cp $(OBJOPT) $(OPTLST) ..
#	@for emod in $(MODOPT); do \
#	 if [ -f $$emod ]; then \
#	  cp $$emod ..; \
#	 fi; \
#	done
	$(RM) $(OPTLST)

lib1:$(OBJOPT1)
	@for emod in $(MODOPT1); do \
	 if [ -f $$emod ]; then \
	  cp $$emod ..; \
	 fi; \
	done

lib2:$(OBJOPT2)
	@for emod in $(MODOPT2); do \
	 if [ -f $$emod ]; then \
	  cp $$emod ..; \
	 fi; \
	done

clean:
	$(RM) *~ *.o *.mod *.trace *.cpp core
