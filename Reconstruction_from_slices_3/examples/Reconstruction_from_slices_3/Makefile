# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = c:/cgal/CGAL-3.2.1/make/makefile_i686_CYGWINNT-5.2-WOW64_CL.EXE-1400
#CGAL_MAKEFILE = /cygdrive/c/cgal/CGAL-3.3-I-120/make/makefile_i686_CYGWINNT-5.1_CL.EXE-1310
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
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_triangulation_points_in_slices$(EXE_EXT) test_triangulation_points_in_slices$(OBJ_EXT) $(LDFLAGS)

test_triangulation_points_in_polylines$(EXE_EXT): test_triangulation_points_in_polylines$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_triangulation_points_in_polylines$(EXE_EXT) test_triangulation_points_in_polylines$(OBJ_EXT) $(LDFLAGS)

test_reconstruction_from_slices$(EXE_EXT): test_reconstruction_from_slices$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_reconstruction_from_slices$(EXE_EXT) test_reconstruction_from_slices$(OBJ_EXT) $(LDFLAGS)

reconstruct$(EXE_EXT): reconstruct$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)reconstruct$(EXE_EXT) reconstruct$(OBJ_EXT) $(LDFLAGS)

test_triangulation_points_in_slices$(OBJ_EXT): test_triangulation_points_in_slices.cpp ../../include/CGAL/Parser_SLC.h ../../include/CGAL/Triangulation_from_slices_3.h ../../include/CGAL/TFS_vertex_base_3.h

test_triangulation_points_in_polylines$(OBJ_EXT): test_triangulation_points_in_polylines.cpp ../../include/CGAL/Parser_CNT.h ../../include/CGAL/Parser_SLC.h ../../include/CGAL/Triangulation_from_slices_3.h ../../include/CGAL/TFS_polyline_vertex_base_3.h

test_reconstruction_from_slices$(OBJ_EXT): test_reconstruction_from_slices.cpp ../../include/CGAL/Parser_CNT.h  ../../include/CGAL/Parser_SLC.h ../../include/CGAL/Triangulation_from_slices_3.h ../../include/CGAL/TFS_polyline_vertex_base_3.h ../../include/CGAL/TFS_cell_base_3.h ../../include/CGAL/Reconstruction_from_slices_3.h

clean: \
	test_triangulation_points_in_slices.clean test_triangulation_points_in_polyline.clean test_reconstruction_from_slices.clean
	rm *.o

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $< -o $(@)

