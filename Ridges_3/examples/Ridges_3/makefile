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
                Compute_Ridges_Umbilics$(EXE_EXT) 

COMPUTE_RIDGES_UMBILICS_OBJS=PolyhedralSurf$(OBJ_EXT) Compute_Ridges_Umbilics$(OBJ_EXT)
Compute_Ridges_Umbilics$(EXE_EXT): $(COMPUTE_RIDGES_UMBILICS_OBJS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Compute_Ridges_Umbilics  $(COMPUTE_RIDGES_UMBILICS_OBJS) \
	$(LDFLAGS) 

clean: \
        Compute_Ridges_Umbilics.clean PolyhedralSurf.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

