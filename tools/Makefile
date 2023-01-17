include Makefile.arch

SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
LIB=ANALYSIS_TOOL.so

$(LIB): $(OBJECTS) 
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJECTS) $(ROOTLIBS) -lTMVA -lEG -lGenVector -lXMLIO -lMLP -lTreePlayer -o $@
	ln -sf $(LIB) lib$(LIB)

%.o:	%.cc
	$(CXX) $(CXXFLAGS) -I${CMSSW_BASE}/../../../external/boost/1.67.0/include -I${CMSSW_BASE}/src -c $< -o $@ -fno-var-tracking

all: $(LIB) 
clean:
	rm -f *.o \
	rm -f *.d \
	rm -f *.so \
