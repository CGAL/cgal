# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) -O2

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

LDFLAGS = $(LONG_NAME_PROBLEM_LDFLAGS) $(CGAL_LDFLAGS)
ifeq (cygwin, $(findstring cygwin,$(OSTYPE)))
LDFLAGS+= ../../src/Benchmark/CGALBenchmark.lib
else
LIBPATH+= $(LIBPATH_OPT)../../src/Benchmark
LDFLAGS+= -lCGALBenchmark
endif

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: check_syntax$(EXE_EXT) simple$(EXE_EXT) leftturn$(EXE_EXT) parm_passing$(EXE_EXT) 

check_syntax$(EXE_EXT): check_syntax$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)check_syntax $? $(LDFLAGS)

simple$(EXE_EXT): simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)simple $? $(LDFLAGS)

leftturn$(EXE_EXT): leftturn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)leftturn $? $(LDFLAGS)

parm_passing$(EXE_EXT): parm_passing$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)parm_passing $? $(LDFLAGS)

clean: check_syntax.clean simple.clean leftturn.clean parm_passing.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# Until CGAL_MAKEFILE supports .cpp extension:

.cpp:
	$(MAKE) $*$(EXE_EXT)

.cpp.clean:
	rm -f $*$(EXE_EXT) $*$(OBJ_EXT) $*.moc
