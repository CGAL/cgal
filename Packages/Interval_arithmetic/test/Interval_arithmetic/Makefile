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
        -I../../../C2/include \
        -I../../../C3/include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS) \
	-g -Wall -W
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

LDFLAGS = Interval_arithmetic.o \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

.PRECIOUS: %.o
.SUFFIXES: .test
.PHONY: %.test test

% : %.C

%.o : %.C
	$(CGAL_CXX) $(CXXFLAGS) -c $<

%.s : %.C
	$(CGAL_CXX) $(CXXFLAGS) -S $<

% : %.o
	$(CGAL_CXX) $+ $(LDFLAGS) -o $@
	./$@

%.test : %
	./$<

EXE_EXT  =
OBJ_EXT  = .o
SOURCES  := $(wildcard  *.C)
OBJECTS  := $(SOURCES:.C=$(OBJ_EXT))
PROGRAMS := $(SOURCES:.C=$(EXE_EXT))
TEST_RUN := $(SOURCES:.C=.test)

all: ../../src/Interval_arithmetic.o $(PROGRAMS)
# all: $(PROGRAMS)

test: $(TEST_RUN)

clean:
	rm -f *$(OBJ_EXT) $(PROGRAMS) *.bb *.bbg *.da *.gcov core *.rpo .pure
