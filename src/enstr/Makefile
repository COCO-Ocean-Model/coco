include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)
ENSLST	= ENSTR-LIST
MODENS1	= $(shell echo $(OBJENS1) | sed -e "s/\.[^[:blank:]]*o/\.mod/g")
MODENS2	= $(shell echo $(OBJENS2) | sed -e "s/\.[^[:blank:]]*o/\.mod/g")

$(ENSLST):$(OBJENS)
	echo $(OBJENS) > $(ENSLST)
	cp $(OBJENS) $(ENSLST) ..
	$(RM) $(ENSLST)

lib1:$(OBJENS1)
	@for emod in $(MODENS1); do \
	 if [-f $$emod]; then \
	  cp $$emod ..; \
	 fi; \
	done

lib2:$(OBJENS2)
	@for emod in $(MODENS2); do \
	 if [-f $$emod]; then \
	  cp $$emod ..; \
	 fi; \
	done

clean:
	$(RM) *~ *.o *.mod *.trace *.cpp core
