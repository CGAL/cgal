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
	   -I../../../Triangulation_3/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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
                alpha_shapes_3$(EXE_EXT) \
                alpha_shapes_with_hierarchy_3$(EXE_EXT) \
                weighted_alpha_shapes_3$(EXE_EXT) 

alpha_shapes_3$(EXE_EXT): alpha_shapes_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)alpha_shapes_3 alpha_shapes_3$(OBJ_EXT) $(LDFLAGS)

alpha_shapes_with_hierarchy_3$(EXE_EXT): alpha_shapes_with_hierarchy_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)alpha_shapes_with_hierarchy_3 alpha_shapes_with_hierarchy_3$(OBJ_EXT) $(LDFLAGS)

weighted_alpha_shapes_3$(EXE_EXT): weighted_alpha_shapes_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)weighted_alpha_shapes_3 weighted_alpha_shapes_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   alpha_shapes_3.clean \
                   alpha_shapes_with_hierarchy_3.clean \
                   weighted_alpha_shapes_3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

