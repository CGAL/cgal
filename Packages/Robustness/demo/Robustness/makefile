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
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                intersection_points_on_segment_2 \
                orientation_2 \
                orientation_IA_2 \
		hull_of_intersection_points_2 \
                triangulation_of_intersection_points_2 \
                cocircular_2 \
                rotation_2

intersection_points_on_segment_2$(EXE_EXT): intersection_points_on_segment_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)intersection_points_on_segment_2 intersection_points_on_segment_2$(OBJ_EXT) $(LDFLAGS)

orientation_2$(EXE_EXT): orientation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)orientation_2 orientation_2$(OBJ_EXT) $(LDFLAGS)

orientation_IA_2$(EXE_EXT): orientation_IA_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)orientation_IA_2 orientation_IA_2$(OBJ_EXT) $(LDFLAGS)

hull_of_intersection_points_2$(EXE_EXT): hull_of_intersection_points_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hull_of_intersection_points_2 hull_of_intersection_points_2$(OBJ_EXT) $(LDFLAGS)

triangulation_of_intersection_points_2$(EXE_EXT): triangulation_of_intersection_points_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_of_intersection_points_2 triangulation_of_intersection_points_2$(OBJ_EXT) $(LDFLAGS)

cocircular_2$(EXE_EXT): cocircular_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)cocircular_2 cocircular_2$(OBJ_EXT) $(LDFLAGS)

rotation_2$(EXE_EXT): rotation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rotation_2 rotation_2$(OBJ_EXT) $(LDFLAGS)


clean: \
                   intersection_points_on_segment_2.clean \
                   orientation_2.clean \
                   orientation_IA_2.clean \
                   hull_of_intersection_points_2.clean \
                   triangulation_of_intersection_points_2.clean \
                   cocircular_2.clean \
                   rotation_2.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

