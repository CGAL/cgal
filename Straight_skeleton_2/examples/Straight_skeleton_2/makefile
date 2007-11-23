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
	Create_straight_skeleton_2$(EXE_EXT) \
	Create_straight_skeleton_from_polygon_with_holes_2$(EXE_EXT) \
	Create_offset_polygons_2$(EXE_EXT) \
	Create_skeleton_and_offset_polygons_2$(EXE_EXT) \
	Create_skeleton_and_offset_polygons_from_polygon_with_holes_2$(EXE_EXT) \
	Create_skeleton_and_offset_polygons_with_holes_2$(EXE_EXT) \
	Low_level_API$(EXE_EXT) 

Create_straight_skeleton_2$(EXE_EXT): Create_straight_skeleton_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Create_straight_skeleton_2$(EXE_EXT) Create_straight_skeleton_2$(OBJ_EXT) $(LDFLAGS)

Create_straight_skeleton_from_polygon_with_holes_2$(EXE_EXT): Create_straight_skeleton_from_polygon_with_holes_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Create_straight_skeleton_from_polygon_with_holes_2$(EXE_EXT) Create_straight_skeleton_from_polygon_with_holes_2$(OBJ_EXT) $(LDFLAGS)

Create_polygon_offset_2$(EXE_EXT): Create_polygon_offset_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Create_polygon_offset_2$(EXE_EXT) Create_polygon_offset_2$(OBJ_EXT) $(LDFLAGS)
	
Create_skeleton_and_offset_polygons_2$(EXE_EXT): Create_skeleton_and_offset_polygons_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Create_skeleton_and_offset_polygons_2$(EXE_EXT) Create_skeleton_and_offset_polygons_2$(OBJ_EXT) $(LDFLAGS)
	
Create_skeleton_and_offset_polygons_from_polygon_with_holes_2$(EXE_EXT): Create_skeleton_and_offset_polygons_from_polygon_with_holes_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Create_skeleton_and_offset_polygons_from_polygon_with_holes_2$(EXE_EXT) Create_skeleton_and_offset_polygons_from_polygon_with_holes_2$(OBJ_EXT) $(LDFLAGS)
	
Create_skeleton_and_offset_polygons_with_holes_2$(EXE_EXT): Create_skeleton_and_offset_polygons_with_holes_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Create_skeleton_and_offset_polygons_with_holes_2$(EXE_EXT) Create_skeleton_and_offset_polygons_with_holes_2$(OBJ_EXT) $(LDFLAGS)
	
Low_level_API$(EXE_EXT): Low_level_API$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Low_level_API$(EXE_EXT) Low_level_API$(OBJ_EXT) $(LDFLAGS)

clean: \
	Create_straight_skeleton_2.clean \
	Create_straight_skeleton_from_polygon_with_holes_2.clean \
	Create_offset_polygons_2.clean \
	Create_skeleton_and_offset_polygons_2.clean \
	Create_skeleton_and_offset_polygons_from_polygon_with_holes_2.clean \
	Create_skeleton_and_offset_polygons_with_holes_2.clean \
	Low_level_API.clean 
	
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

