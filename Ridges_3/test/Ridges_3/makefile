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
                ridge_test$(EXE_EXT) 

PolyhedralSurf$(OBJ_EXT): PolyhedralSurf.cpp
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) PolyhedralSurf.cpp

ridge_test$(OBJ_EXT): ridge_test.cpp
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) ridge_test.cpp

ridge_test$(EXE_EXT): ridge_test$(OBJ_EXT) PolyhedralSurf$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ridge_test \
	ridge_test$(OBJ_EXT) PolyhedralSurf$(OBJ_EXT) $(LDFLAGS)

clean: \
                   PolyhedralSurf.clean \
                   ridge_test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

