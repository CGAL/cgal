# Created by the script create_makefile
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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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
                barycenter$(EXE_EXT) \
                bounding_box$(EXE_EXT) \
                fitting_2$(EXE_EXT) \
                fitting_3$(EXE_EXT) \

barycenter$(EXE_EXT): barycenter$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)barycenter barycenter$(OBJ_EXT) $(LDFLAGS)

bounding_box$(EXE_EXT): bounding_box$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bounding_box bounding_box$(OBJ_EXT) $(LDFLAGS)

fitting_2$(EXE_EXT): fitting_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)fitting_2 fitting_2$(OBJ_EXT) $(LDFLAGS)

fitting_3$(EXE_EXT): fitting_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)fitting_3 fitting_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   barycenter.clean \
                   bounding_box.clean \
                   fitting_2.clean \
                   fitting_3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

