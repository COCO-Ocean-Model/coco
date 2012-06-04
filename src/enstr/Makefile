include ../Mkinclude
include $(OGCMDIR)/src/Sysdep.$(SYSTEM)$(PARALLEL)
ENSLST	= ENSTR-LIST

$(ENSLST):$(OBJENS)
	echo $(OBJENS) > $(ENSLST)
	cp $(OBJENS) $(ENSLST) ..
	$(RM) $(ENSLST)

lib1:$(OBJENS1)
	cp $(OBJENS1) ..

lib2:$(OBJENS2)
	cp $(OBJENS2) ..

enslist:
	echo $(OBJENS) > ../$(ENSLST)

clean:
	$(RM) *~ *.o *.trace *.cpp core
