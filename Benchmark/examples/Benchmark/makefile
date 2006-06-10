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
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH) -L../../src/Benchmark

LDFLAGS = $(LONG_NAME_PROBLEM_LDFLAGS) $(CGAL_LDFLAGS) -lCGALBenchmark

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: check_syntax$(EXE_EXT) simple$(EXE_EXT) leftturn$(EXE_EXT)

check_syntax$(EXE_EXT): check_syntax$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)check_syntax $? $(LDFLAGS)

simple$(EXE_EXT): simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)simple $? $(LDFLAGS)

leftturn$(EXE_EXT): leftturn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)leftturn $? $(LDFLAGS)

clean: check_syntax.clean simple.clean leftturn.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
