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
           -I../../../H3/include \
           -I../../../H2/include \
           -I../../../Cartesian_kernel/include \
           -I../../../Number_types/include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
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

all:            \
                Cartesian$(EXE_EXT) \
                Homogeneous$(EXE_EXT) \
                Simple_cartesian$(EXE_EXT) \
                Simple_homogeneous$(EXE_EXT) \
                test_kernel__$(EXE_EXT) \
		test_kernel_archetype_2$(EXE_EXT) \
		test_kernel_archetype_3$(EXE_EXT)
		

Cartesian$(EXE_EXT): Cartesian$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Cartesian Cartesian$(OBJ_EXT) $(LDFLAGS)

Homogeneous$(EXE_EXT): Homogeneous$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Homogeneous Homogeneous$(OBJ_EXT) $(LDFLAGS)

Simple_cartesian$(EXE_EXT): Simple_cartesian$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Simple_cartesian Simple_cartesian$(OBJ_EXT) $(LDFLAGS)

Simple_homogeneous$(EXE_EXT): Simple_homogeneous$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Simple_homogeneous Simple_homogeneous$(OBJ_EXT) $(LDFLAGS)

test_kernel__$(EXE_EXT): test_kernel__$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernel__ test_kernel__$(OBJ_EXT) $(LDFLAGS)

test_kernel_archetype_2$(EXE_EXT): test_kernel_archetype_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernel_archetype_2 test_kernel_archetype_2$(OBJ_EXT) $(LDFLAGS)

test_kernel_archetype_3$(EXE_EXT): test_kernel_archetype_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernel_archetype_3 test_kernel_archetype_3$(OBJ_EXT) $(LDFLAGS)


clean: \
                   Cartesian.clean \
                   Homogeneous.clean \
                   Simple_cartesian.clean \
                   Simple_homogeneous.clean \
                   test_kernel__.clean \
		   test_kernel_archetype_2.clean \
		   test_kernel_archetype_3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

