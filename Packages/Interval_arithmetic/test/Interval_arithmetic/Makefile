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
# -I../../include -I$(HOME)/CGAL-2.1-I-9.modified/include \

CXXFLAGS = \
        -I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS) \
	# /Oi /Fatst2.s
	# /Fatst1.s
	# /Op /Oi /Od /Za /O2 /Ox
	# -g -Wall -W
	# -O2 -Wall -W
	# -fpermissive
	# -DCGAL_IA_STOP_CONSTANT_PROPAGATION \
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS\
	# -DCGAL_CHECK_EXPENSIVE \
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

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

# LDFLAGS = Interval_arithmetic$(OBJ_EXT)

LDFLAGS =  \
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
