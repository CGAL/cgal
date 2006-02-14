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
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                construction_test_suite_generator$(EXE_EXT) \
                ex_kernel_point$(EXE_EXT) \
                ex_kernel_segment$(EXE_EXT) \
                point_location$(EXE_EXT) \
                test_construction$(EXE_EXT) \
                test_insertion$(EXE_EXT) \
                test_iso_verts$(EXE_EXT) \
                test_removal$(EXE_EXT) \
                test_traits$(EXE_EXT) \
                test_vert_ray_shoot_vert_segments$(EXE_EXT) 

construction_test_suite_generator$(EXE_EXT): construction_test_suite_generator$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)construction_test_suite_generator construction_test_suite_generator$(OBJ_EXT) $(LDFLAGS)

ex_kernel_point$(EXE_EXT): ex_kernel_point$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ex_kernel_point ex_kernel_point$(OBJ_EXT) $(LDFLAGS)

ex_kernel_segment$(EXE_EXT): ex_kernel_segment$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ex_kernel_segment ex_kernel_segment$(OBJ_EXT) $(LDFLAGS)

point_location$(EXE_EXT): point_location$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)point_location point_location$(OBJ_EXT) $(LDFLAGS)

test_construction$(EXE_EXT): test_construction$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_construction test_construction$(OBJ_EXT) $(LDFLAGS)

test_insertion$(EXE_EXT): test_insertion$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_insertion test_insertion$(OBJ_EXT) $(LDFLAGS)

test_iso_verts$(EXE_EXT): test_iso_verts$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_iso_verts test_iso_verts$(OBJ_EXT) $(LDFLAGS)

test_removal$(EXE_EXT): test_removal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_removal test_removal$(OBJ_EXT) $(LDFLAGS)

test_traits$(EXE_EXT): test_traits$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_traits test_traits$(OBJ_EXT) $(LDFLAGS)

test_vert_ray_shoot_vert_segments$(EXE_EXT): test_vert_ray_shoot_vert_segments$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_vert_ray_shoot_vert_segments test_vert_ray_shoot_vert_segments$(OBJ_EXT) $(LDFLAGS)

clean: \
                   construction_test_suite_generator.clean \
                   ex_kernel_point.clean \
                   ex_kernel_segment.clean \
                   point_location.clean \
                   test_construction.clean \
                   test_insertion.clean \
                   test_iso_verts.clean \
                   test_removal.clean \
                   test_traits.clean \
                   test_vert_ray_shoot_vert_segments.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

