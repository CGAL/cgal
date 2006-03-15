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
	-I ../../include\
	-O2 -DNDEBUG \
	$(CGAL_CXXFLAGS) \
	$(LONG_NAME_PROBLEM_CXXFLAGS)

#	$(DEBUG_OPT) -DTRACE\
#
#	-O2 -DNDEBUG \
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
	test_triangulation_points_in_slices$(EXE_EXT) test_triangulation_points_in_polylines$(EXE_EXT) test_reconstruction_from_slices$(EXE_EXT)

test_triangulation_points_in_slices$(EXE_EXT): test_triangulation_points_in_slices$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) bin/test_triangulation_points_in_slices$(EXE_EXT) obj/test_triangulation_points_in_slices$(OBJ_EXT) $(LDFLAGS)

test_triangulation_points_in_polylines$(EXE_EXT): test_triangulation_points_in_polylines$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) bin/test_triangulation_points_in_polylines$(EXE_EXT) obj/test_triangulation_points_in_polylines$(OBJ_EXT) $(LDFLAGS)

test_reconstruction_from_slices$(EXE_EXT): test_reconstruction_from_slices$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) bin/test_reconstruction_from_slices$(EXE_EXT) obj/test_reconstruction_from_slices$(OBJ_EXT) $(LDFLAGS)

test_triangulation_points_in_slices$(OBJ_EXT): test_triangulation_points_in_slices.C ../../include/CGAL/Parser_SLC.h ../../include/CGAL/Triangulation_from_slices_3.h ../../include/CGAL/TFS_vertex_base_3.h

test_triangulation_points_in_polylines$(OBJ_EXT): test_triangulation_points_in_polylines.C ../../include/CGAL/Parser_CNT.h ../../include/CGAL/Parser_SLC.h ../../include/CGAL/Triangulation_from_slices_3.h ../../include/CGAL/TFS_polyline_vertex_base_3.h

test_reconstruction_from_slices$(OBJ_EXT): test_reconstruction_from_slices.C ../../include/CGAL/Parser_CNT.h  ../../include/CGAL/Parser_SLC.h ../../include/CGAL/Triangulation_from_slices_3.h ../../include/CGAL/TFS_polyline_vertex_base_3.h ../../include/CGAL/TFS_cell_base_3_for_reconstruction.h ../../include/CGAL/Reconstruction_from_slices_3.h

clean: \
	test_triangulation_points_in_slices.clean test_triangulation_points_in_polyline.clean test_reconstruction_from_slices.clean
	rm obj/* exe/*

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $< -o obj/$(@)

