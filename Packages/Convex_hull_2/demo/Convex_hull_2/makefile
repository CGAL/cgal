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
                ch_of_polyline \
                convex_hull_2_demo \
                convex_hull_of_ipoints 

ch_of_polyline$(EXE_EXT): ch_of_polyline$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_of_polyline ch_of_polyline$(OBJ_EXT) $(LDFLAGS)

convex_hull_2_demo$(EXE_EXT): convex_hull_2_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex_hull_2_demo convex_hull_2_demo$(OBJ_EXT) $(LDFLAGS)

convex_hull_of_ipoints$(EXE_EXT): convex_hull_of_ipoints$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex_hull_of_ipoints convex_hull_of_ipoints$(OBJ_EXT) $(LDFLAGS)

clean: \
                   ch_of_polyline.clean \
                   convex_hull_2_demo.clean \
                   convex_hull_of_ipoints.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

