# This is the makefile for compiling a CGAL application.

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# CGAL_CXX=$(HOME)/egcs/$(UNAME)/bin/g++
# Pour utiliser Purify:
# CGAL_CXX=purify g++
# CGAL_LDFLAGS= $(CGAL_LIB_DIR)/$(CGAL_OS_COMPILER)/libCGAL.a \
              # -L$(GMP_LIB_DIR) -lgmp -lm \

CXXFLAGS = \
        -I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS) \
	-O2 -DLOOPS=10000000 -W -Wall
	# -DCGAL_IA_NEW_FILTERS
	# -O2 -mcpu=pentiumpro -march=pentiumpro -W -Wall -DLOOPS=10000000
	# -g -Wall -W
	# -O2 -Wall -W
	# /Fatst1.s
	# /Op /Oi /Od /Za /O2 /Ox
	# -fpermissive
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS\
	# -DCGAL_NO_WARNINGS
	# -DCGAL_CHECK_EXPENSIVE \
	# -DCGAL_IA_DONT_STOP_CONSTANT_PROPAGATION
	# -DCGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT
	# -O2 -W -Winline -Wunused -malign-double
	# -mcpu=pentiumpro -march=pentiumpro -mieee-fp
	# -mcpu=ultrasparc -mtune=ultrasparc
	# -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	# -fomit-frame-pointer
	# -fprofile-arcs -ftest-coverage
	# -Wa,-xarch=v8plus   # g++ sans -m*=ultra*
	# -fhonor-std
	# -frepo
	# -pedantic
	# -DLOOPS=10000000

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

# LDFLAGS = Interval_arithmetic$(OBJ_EXT)
# LDFLAGS = ../../src/Interval_arithmetic/workaround_4_ms.obj

LDFLAGS = Interval_arithmetic$(OBJ_EXT) \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

.PRECIOUS: %$(OBJ_EXT)
.SUFFIXES: .test
.PHONY: %.test test

%$(EXE_EXT) : %.C

%$(OBJ_EXT) : %.C
	$(CGAL_CXX) $(CXXFLAGS) -c $<

%.s : %.C
	$(CGAL_CXX) $(CXXFLAGS) -S $<

%$(EXE_EXT) : %$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)$@ $+ $(LDFLAGS)
	./$@

%.test : %$(EXE_EXT)
	./$<

SOURCES  := $(wildcard  *.C)
OBJECTS  := $(SOURCES:.C=$(OBJ_EXT))
PROGRAMS := $(SOURCES:.C=$(EXE_EXT))
TEST_RUN := $(SOURCES:.C=.test)

all: ../../src/Interval_arithmetic$(OBJ_EXT) $(PROGRAMS)
# all: $(PROGRAMS)

test: $(TEST_RUN)

clean:
	rm -f *.o $(PROGRAMS) *.bb *.bbg *.da *.gcov core *.rpo .pure
	rm -f *.obj *.exe
