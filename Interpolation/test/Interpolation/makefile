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
           -Iinclude \
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
                test_interpolation_functions_2$(EXE_EXT) \
                test_regular_neighbors_2$(EXE_EXT) \
                test_surface_neighbor_3$(EXE_EXT) 

test_interpolation_functions_2$(EXE_EXT): test_interpolation_functions_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_interpolation_functions_2 test_interpolation_functions_2$(OBJ_EXT) $(LDFLAGS)

test_regular_neighbors_2$(EXE_EXT): test_regular_neighbors_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_regular_neighbors_2 test_regular_neighbors_2$(OBJ_EXT) $(LDFLAGS)

test_surface_neighbor_3$(EXE_EXT): test_surface_neighbor_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_surface_neighbor_3 test_surface_neighbor_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_interpolation_functions_2.clean \
                   test_regular_neighbors_2.clean \
                   test_surface_neighbor_3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

