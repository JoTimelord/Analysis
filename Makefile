SUBDIRS := runlooper
LIBRARIES := $${HOME}/Utils/NanoTools/NanoCORE $${HOME}/Utils/rapido tools

all: $(LIBRARIES) $(SUBDIRS)

NanoCORE:
	$(MAKE) -C $${HOME}/Utils/NanoTools/NanoCORE

rapido:
	$(MAKE) -C $${HOME}/Utils/rapido

tools:
	$(MAKE) -C tools

$(SUBDIRS): NanoCORE rapido tools
	$(MAKE) -C $@

.PHONY: all $(LIBRARIES) $(SUBDIRS)

clean:
	cd runlooper/ && make clean;

cleanall:
	cd $${HOME}/Utils/rapido/ && make clean;
	cd tools/ && make clean;
	cd $${HOME}/Utils/NanoTools/NanoCORE/ && make clean;
	cd runlooper/ && make clean;
