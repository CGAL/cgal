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
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Circle_arr_from_file$(EXE_EXT) \
                Polyline_arr_from_file$(EXE_EXT) \
                Polyline_arr_from_mouse$(EXE_EXT) \
                Seg_circ_arr_from_file$(EXE_EXT) \
                Segment_arr_from_file$(EXE_EXT) \
                Segment_arr_from_mouse$(EXE_EXT) \
                Polygon_intersection_demo$(EXE_EXT) 

Circle_arr_from_file$(EXE_EXT): Circle_arr_from_file$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Circle_arr_from_file Circle_arr_from_file$(OBJ_EXT) $(LDFLAGS)

Polyline_arr_from_file$(EXE_EXT): Polyline_arr_from_file$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyline_arr_from_file Polyline_arr_from_file$(OBJ_EXT) $(LDFLAGS)

Polyline_arr_from_mouse$(EXE_EXT): Polyline_arr_from_mouse$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyline_arr_from_mouse Polyline_arr_from_mouse$(OBJ_EXT) $(LDFLAGS)

Seg_circ_arr_from_file$(EXE_EXT): Seg_circ_arr_from_file$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Seg_circ_arr_from_file Seg_circ_arr_from_file$(OBJ_EXT) $(LDFLAGS)

Segment_arr_from_file$(EXE_EXT): Segment_arr_from_file$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_arr_from_file Segment_arr_from_file$(OBJ_EXT) $(LDFLAGS)

Segment_arr_from_mouse$(EXE_EXT): Segment_arr_from_mouse$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_arr_from_mouse Segment_arr_from_mouse$(OBJ_EXT) $(LDFLAGS)

Polygon_intersection_demo$(EXE_EXT): Polygon_intersection_demo$(OBJ_EXT)  Polygon_intersection$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polygon_intersection_demo Polygon_intersection_demo$(OBJ_EXT)  Polygon_intersection$(OBJ_EXT) $(LDFLAGS)
   
clean: \
                   Circle_arr_from_file.clean \
                   Polyline_arr_from_file.clean \
                   Polyline_arr_from_mouse.clean \
                   Seg_circ_arr_from_file.clean \
                   Segment_arr_from_file.clean \
                   Segment_arr_from_mouse.clean \
                   Polygon_intersection.clean \
                   Polygon_intersection_demo.clean
  
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

