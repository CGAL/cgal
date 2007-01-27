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
                edge_collapse_enriched_polyhedron$(EXE_EXT) \
                edge_collapse_polyhedron$(EXE_EXT) 

edge_collapse_enriched_polyhedron$(OBJ_EXT): edge_collapse_enriched_polyhedron.cpp
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) edge_collapse_enriched_polyhedron.cpp

edge_collapse_enriched_polyhedron$(EXE_EXT): edge_collapse_enriched_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)edge_collapse_enriched_polyhedron edge_collapse_enriched_polyhedron$(OBJ_EXT) $(LDFLAGS)

edge_collapse_polyhedron$(OBJ_EXT): edge_collapse_polyhedron.cpp
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) edge_collapse_polyhedron.cpp

edge_collapse_polyhedron$(EXE_EXT): edge_collapse_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)edge_collapse_polyhedron edge_collapse_polyhedron$(OBJ_EXT) $(LDFLAGS)

clean: \
                   edge_collapse_enriched_polyhedron.clean \
                   edge_collapse_polyhedron.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

