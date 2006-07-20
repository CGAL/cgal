# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = c:/cgal/CGAL-3.2-I-467/make/makefile_i686_CYGWINNT-5.1_CL.EXE-1310

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	   -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           -DSTRAIGHTSKELETONDLL_EXPORTS

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

all:            \
                StraightSkeleton$(EXE_EXT) \
	 	main$(EXE_EXT)

StraightSkeleton$(EXE_EXT): StraightSkeleton$(OBJ_EXT)
	$(CGAL_CXX)  -LD $(LIBPATH) $(EXE_OPT)StraightSkeleton StraightSkeleton$(OBJ_EXT) $(LDFLAGS)

main.exe : main.obj
	cl.exe   -Fe main main.obj StraightSkeleton.lib

clean: \
                   StraightSkeleton.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

