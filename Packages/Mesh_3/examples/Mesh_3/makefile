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

CXXFLAGS = -I../../include \
           -I../../../Mesh_2/include\
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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
                essai$(EXE_EXT) 

HEADERS = ../../include/CGAL/Constrained_regular_triangulation_3.h \
	  ../../include/CGAL/Constrained_triangulation_vertex_base_3.h \
	  ../../include/CGAL/Constrained_triangulation_cell_base_3.h

essai$(OBJ_EXT): essai.C $(HEADERS)

essai$(EXE_EXT): essai$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)essai essai$(OBJ_EXT) $(LDFLAGS)

clean: \
                   essai.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

