#CGAL_MAKEFILE=/usr/local/cgal/CGAL-2.5-I-126/make/makefile_ppc_Linux-2.4.22-ben1_g++-3.2-3.2.3

include $(CGAL_MAKEFILE)

SOURCES := bench.C
OBJS := $(patsubst %.C,%.o,$(SOURCES))

CXXFLAGS := -I../../include $(CGAL_CXXFLAGS)
#CXXFLAGS += -O2
#CXXFLAGS += -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS -DCGAL_NO_ASSERTIONS -DCGAL_NO_WARNINGS
CXXFLAGS += -DCHR_NOTRACE
#CXXFLAGS += -DCHR_LEDA
CXXFLAGS += -DCHR_FILTERED_EXACT
#CXXFLAGS += -DCHR_GMP
#CXXFLAGS += -pg


LDFLAGS := $(CGAL_LDFLAGS)
#LDFLAGS += -pg -ldl

all: .depend bench

bench: $(OBJS)
	$(CGAL_CXX) -o $@ $(OBJS) $(LDFLAGS)

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f bench.o *~ bench

.depend: $(shell sed -e 's/[^:]*://' -e 's/\\$$//' <.depend)
	$(CGAL_CXX) $(CXXFLAGS) -MM $(SOURCES) >.depend

-include .depend
