include Mkinclude
export OGCMDIR := $(shell dirname $(CURDIR))
include Sysdep.$(SYSTEM)
MAIN    = $(MODEL).o

all:lib $(TARGET)

$(TARGET):$(LIBOGCM) $(MAIN)
	$(FC) $(FOPTION) -o $@ $(MAIN) $(LIBOGCM)

lib: modules
	$(MAKE) -C build $@

modules:
	$(MAKE) -C module
ifeq "$(OPTTOUZA)" "-DOPT_TOUZA"
	$(MAKE) -C touza
endif
	@touch build/OBJLIST
	$(MAKE) -C build OBJLIST

clean:
	$(RM) *~ *.o *.mod *.trace *.L core ../lib/*
	$(MAKE) -C module $@
	$(MAKE) -C touza $@
	$(MAKE) -C build $@
