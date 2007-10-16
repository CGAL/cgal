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

CXXFLAGS = -I ../../../Mesh_2/include \
           -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(OTHER_CXXFLAGS)

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
                combined_spheres$(EXE_EXT) \
                implicit_surface_mesher_test$(EXE_EXT) \
                internal-test-surface-mesher$(EXE_EXT) \
                test_c2t3_iterators$(EXE_EXT) \
                test_robust_circumcenter$(EXE_EXT) 

combined_spheres$(EXE_EXT): combined_spheres$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)combined_spheres combined_spheres$(OBJ_EXT) $(LDFLAGS)

implicit_surface_mesher_test$(EXE_EXT): implicit_surface_mesher_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surface_mesher_test implicit_surface_mesher_test$(OBJ_EXT) $(LDFLAGS)

internal-test-surface-mesher$(EXE_EXT): internal-test-surface-mesher$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)internal-test-surface-mesher internal-test-surface-mesher$(OBJ_EXT) $(LDFLAGS)

test_c2t3_iterators$(EXE_EXT): test_c2t3_iterators$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_c2t3_iterators test_c2t3_iterators$(OBJ_EXT) $(LDFLAGS)

test_robust_circumcenter$(EXE_EXT): test_robust_circumcenter$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_robust_circumcenter test_robust_circumcenter$(OBJ_EXT) $(LDFLAGS)

run-test: all
	./implicit_surface_mesher_test$(EXE_EXT)
	./test_robust_circumcenter$(EXE_EXT)
	./test_c2t3_iterators$(EXE_EXT)

clean: \
                   combined_spheres.clean \
                   implicit_surface_mesher_test.clean \
                   internal-test-surface-mesher.clean \
                   test_c2t3_iterators.clean \
                   test_robust_circumcenter.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

#---------------------------------------------------------------------#
#
#                        dependencies
#
#      if you want deps, create a file my_makefile that contains
#      -include makefile.dependencies
#
#---------------------------------------------------------------------#

-include my_makefile
