# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-69/make/makefile_i686_CYGWINNT-5.0-1.3.12_g++-2.95
CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-69/make/makefile_i686_CYGWINNT-5.0-1.3.12_ICL.EXE-7
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../../include\
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            main$(EXE_EXT)

main$(EXE_EXT): main$(OBJ_EXT)  ../../../include/CGAL/Interval_skip_list.h
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)main main$(OBJ_EXT)  $(LDFLAGS)

clean: 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

