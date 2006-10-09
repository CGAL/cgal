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

CXXFLAGS = -I../../include/ \
           -I../../../BGL/include/ \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

LDFLAGS = $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                LT_edge_collapse_polyhedron \
                LT_edge_collapse_enriched_polyhedron \
                MP_edge_collapse_polyhedron \
                Edge_collapse_fully_cached_polyhedron \
                LT_edge_collapse_visited_polyhedron \
                LT_edge_collapse_with_fixed_vertices_polyhedron \
                LT_edge_collapse_with_fixed_vertices_enriched_polyhedron \

LT_edge_collapse_polyhedron$(EXE_EXT): LT_edge_collapse_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)LT_edge_collapse_polyhedron \
        LT_edge_collapse_polyhedron$(OBJ_EXT) $(LDFLAGS)
        
LT_edge_collapse_enriched_polyhedron$(EXE_EXT): LT_edge_collapse_enriched_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)LT_edge_collapse_enriched_polyhedron \
        LT_edge_collapse_enriched_polyhedron$(OBJ_EXT) $(LDFLAGS)
        
MP_edge_collapse_polyhedron$(EXE_EXT): MP_edge_collapse_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MP_edge_collapse_polyhedron \
        MP_edge_collapse_polyhedron$(OBJ_EXT) $(LDFLAGS)

Edge_collapse_fully_cached_polyhedron$(EXE_EXT): Edge_collapse_fully_cached_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Edge_collapse_fully_cached_polyhedron \
        Edge_collapse_fully_cached_polyhedron$(OBJ_EXT) $(LDFLAGS)

LT_edge_collapse_visited_polyhedron$(EXE_EXT): LT_edge_collapse_visited_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)LT_edge_collapse_visited_polyhedron \
        LT_edge_collapse_visited_polyhedron$(OBJ_EXT) $(LDFLAGS)
        
LT_edge_collapse_with_fixed_vertices_polyhedron$(EXE_EXT): LT_edge_collapse_with_fixed_vertices_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)LT_edge_collapse_with_fixed_vertices_polyhedron \
        LT_edge_collapse_with_fixed_vertices_polyhedron$(OBJ_EXT) $(LDFLAGS)

LT_edge_collapse_with_fixed_vertices_enriched_polyhedron$(EXE_EXT): LT_edge_collapse_with_fixed_vertices_enriched_polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)LT_edge_collapse_with_fixed_vertices_enriched_polyhedron \
        LT_edge_collapse_with_fixed_vertices_enriched_polyhedron$(OBJ_EXT) $(LDFLAGS)

clean: \
                   LT_edge_collapse_polyhedron.clean \
                   LT_edge_collapse_enriched_polyhedron.clean \
                   MP_edge_collapse_polyhedron.clean \
                   Edge_collapse_fully_cached_polyhedron.clean \
                   LT_edge_collapse_visited_polyhedron.clean \
                   LT_edge_collapse_with_fixed_vertices_polyhedron.clean \
                   LT_edge_collapse_with_fixed_vertices_enriched_polyhedron.clean \

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

