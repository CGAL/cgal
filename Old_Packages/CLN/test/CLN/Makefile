# This is the makefile for compiling a CGAL application.

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
        -I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS)
	# -g -Wall -W
	# /Op 
	# -O2 -Wall -W
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS\
	# -DCGAL_CHECK_EXPENSIVE \
	# -O2 -W -Winline -Wunused -malign-double
	# -mcpu=pentiumpro -march=pentiumpro -mieee-fp
	# -mcpu=ultrasparc -mtune=ultrasparc
	# -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	# -fomit-frame-pointer
	# -frepo
	# -pedantic

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

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

all: $(PROGRAMS)

test: $(TEST_RUN)

clean:
	rm -f *.o $(PROGRAMS) *.bb *.bbg *.da *.gcov core *.rpo .pure
	rm -f *.obj *.exe test_io.tmp
