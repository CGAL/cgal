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
	-g -W # -Winline -O2 -mcpu=pentiumpro -march=pentiumpro -fstrict-aliasing # -pedantic
	# -g -W -Wall -Winline
	# -O2 -mcpu=pentiumpro -march=pentiumpro -W -Wall -DLOOPS=10000000
	# -O2 -W -Winline -Wunused -malign-double
	# -I$(CGAL)/include_advanced    # Advanced kernel, cf config.h aussi.
	# -mcpu=pentiumpro -march=pentiumpro -mieee-fp
	# -mcpu=ultrasparc -mtune=ultrasparc
	# -fomit-frame-pointer
	# -fprofile-arcs -ftest-coverage # for GCOV
	# -fpermissive
	# -fno-enforce-eh-specs # for GCC 3
	# -fhonor-std
	# -frepo
	# -pedantic
	# -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	# -Wa,-xarch=v8plus   # g++ sans -m*=ultra*
	# /Fatst1.s
	# /Op /Oi /Od /Za /O2 /Ox
	# -DCGAL_IA_NEW_FILTERS \
	# -DCGAL_CHECK_EXPENSIVE \
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS \
	# -DCGAL_NO_POSTCONDITIONS -DCGAL_NO_WARNINGS \
	# -DCGAL_IA_DONT_STOP_CONSTANT_PROPAGATION \
	# -DCGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT \
	# -DLOOPS=10000000

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

# LDFLAGS = Interval_arithmetic$(OBJ_EXT)
# LDFLAGS = ../../src/Interval_arithmetic/workaround_4_ms.obj

LDFLAGS = \
         $(TESTSUITE_LDFLAGS) \
         $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

.PRECIOUS: %$(OBJ_EXT) tst1 %$(EXE_EXT)
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
