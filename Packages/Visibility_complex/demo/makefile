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

CXXFLAGS = -I../include \
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
                shortest_path$(EXE_EXT) \
                visibility_graph$(EXE_EXT) 

shortest_path$(EXE_EXT): shortest_path$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)shortest_path shortest_path$(OBJ_EXT) $(LDFLAGS)

visibility_graph$(EXE_EXT): visibility_graph$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visibility_graph visibility_graph$(OBJ_EXT) $(LDFLAGS)

clean: \
                   shortest_path.clean \
                   visibility_graph.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

