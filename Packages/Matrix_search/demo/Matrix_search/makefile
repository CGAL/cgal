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
                all_furthest_neighbors_2_demo$(EXE_EXT) \
                extremal_polygon_2_demo$(EXE_EXT) \
                rectangular_p_center_2_demo$(EXE_EXT) 

all_furthest_neighbors_2_demo$(EXE_EXT): all_furthest_neighbors_2_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)all_furthest_neighbors_2_demo all_furthest_neighbors_2_demo$(OBJ_EXT) $(LDFLAGS)

extremal_polygon_2_demo$(EXE_EXT): extremal_polygon_2_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)extremal_polygon_2_demo extremal_polygon_2_demo$(OBJ_EXT) $(LDFLAGS)

rectangular_p_center_2_demo$(EXE_EXT): rectangular_p_center_2_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rectangular_p_center_2_demo rectangular_p_center_2_demo$(OBJ_EXT) $(LDFLAGS)

clean: \
                   all_furthest_neighbors_2_demo.clean \
                   extremal_polygon_2_demo.clean \
                   rectangular_p_center_2_demo.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

