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
           -I../../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                simple$(EXE_EXT) 

Bench$(OBJ_EXT): Bench.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

Bench_option_parser$(OBJ_EXT): Bench_option_parser.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

simple$(EXE_EXT): simple$(OBJ_EXT) Bench_option_parser$(OBJ_EXT) Bench$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)simple $? $(LDFLAGS)

clean: \
                   simple.clean \
		   Bench_option_parser.clean \
		   Bench.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

vpath %.C ../../../include/CGAL
