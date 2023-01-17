SUBDIRS := runlooper
LIBRARIES := $${HOME}/Utils/NanoTools/NanoCORE $${HOME}/Utils/rooutil tools

all: $(LIBRARIES) $(SUBDIRS)

NanoCORE:
	$(MAKE) -C $${HOME}/Utils/NanoTools/NanoCORE

rooutil:
	$(MAKE) -C $${HOME}/Utils/rooutil

tools:
	$(MAKE) -C tools

$(SUBDIRS): NanoCORE rooutil tools
	$(MAKE) -C $@

.PHONY: all $(LIBRARIES) $(SUBDIRS)

clean:
	cd runlooper/ && make clean;

cleanall:
	cd $${HOME}/Utils/rooutil/ && make clean;
	cd tools/ && make clean;
	cd $${HOME}/Utils/NanoTools/NanoCORE/ && make clean;
	cd runlooper/ && make clean;
