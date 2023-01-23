# Simple makefile
MAINDIR=../
EXE=$(MAINDIR)/runLooper

SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
HEADERS=$(SOURCES:.cc=.h)

CORRECTIONLIBDIR=${CMSSW_BASE}/../../../external/py3-correctionlib/2.0.0-0c4f44c8dd5561d8c0660135feeb81f4/lib/python3.9/site-packages/correctionlib

CC          = g++
CXX         = g++
CXXFLAGS    = -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual
LD          = g++
SOFLAGS     = -g -shared
CXXFLAGS    = -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual
LDFLAGS     = -g -O2 -L$${HOME}/Utils/rapido/src -L$${HOME}/Utils/NanoTools/NanoCORE -L$(MAINDIR)/tools -L${CMSSW_BASE}/lib/${SCRAM_ARCH}
LDFLAGS    += -L$(CORRECTIONLIBDIR)/lib
LDFLAGS    += $(shell root-config --glibs)
ROOTLIBS    = $(shell root-config --libs)
ROOTCFLAGS  = $(shell root-config --cflags)
CXXFLAGS   += $(ROOTCFLAGS)
CFLAGS      = $(ROOTCFLAGS) -Wall -Wno-unused-function -g -O2 -fPIC -fno-var-tracking
CFLAGS     += -I$${HOME}/Utils/NanoTools/NanoCORE -I$${HOME}/Utils/rapido/src -I$(MAINDIR)/tools -I${CMSSW_BASE}/src -DLorentzVectorPtEtaPhiM4D
CFLAGS     += -I${CMSSW_BASE}/../../../external/boost/1.67.0/include
CFLAGS     += -I$(CORRECTIONLIBDIR)/include
EXTRAFLAGS  = -fPIC -ITMultiDrawTreePlayer -Wunused-variable -lTMVA -lEG -lGenVector -lXMLIO -lMLP -lTreePlayer -lImt
EXTRAFLAGS += -lRAPIDO -lNANO_CORE -lCondFormatsJetMETObjects -lJetMETCorrectionsModules -lcorrectionlib

$(EXE): $(OBJECTS)
	$(LD) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(ROOTLIBS) $(EXTRAFLAGS) -o $@

%.o: %.cc
	$(CC) $(CFLAGS) $(EXTRACFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS) $(EXE)