# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include \
	   -I../../../Data_structure_for_queries_3/include \
           $(OTHER_INCLUDES) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(OTHER_CXXFLAGS)

#           -I../../../Mesh_2/include \
 #          -I../../../Mesh_3/include \
  #         -I../../../Triangulation_3/include \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

OBJECTS = \
	implicit_function$(OBJ_EXT)

SURFACE_MESHER_POLYHEDRAL_OPT = -DSURFACE_MESHER_POLYHEDRAL

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                implicit_surface_mesher$(EXE_EXT) \
		polyhedral_surface_mesher$(EXE_EXT) \
                surfaces_union$(EXE_EXT)
#		test$(EXE_EXT)	
#                surface_mesher_ellipsoid$(EXE_EXT) \
#		surface_mesher_with_hierarchy$(EXE_EXT)

implicit_surface_mesher$(OBJ_EXT): surface_mesher.C
	$(CGAL_CXX) $(CXXFLAGS) $(EXE_OPT) implicit_surface_mesher$(OBJ_EXT) $(OBJ_OPT) surface_mesher.C

implicit_surface_mesher$(EXE_EXT): implicit_surface_mesher$(OBJ_EXT) $(OBJECTS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surface_mesher implicit_surface_mesher$(OBJ_EXT) $(OBJECTS) $(LDFLAGS)

polyhedral_surface_mesher$(OBJ_EXT): surface_mesher.C
	$(CGAL_CXX) $(CXXFLAGS) $(SURFACE_MESHER_POLYHEDRAL_OPT) $(EXE_OPT) polyhedral_surface_mesher$(OBJ_EXT) $(OBJ_OPT) surface_mesher.C

polyhedral_surface_mesher$(EXE_EXT): polyhedral_surface_mesher$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedral_surface_mesher polyhedral_surface_mesher$(OBJ_EXT) $(LDFLAGS)

surfaces_union$(EXE_EXT): surfaces_union$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)surfaces_union surfaces_union$(OBJ_EXT) $(LDFLAGS)

clean: implicit_function.clean surfaces_union.clean
	rm -f polyhedral_surface_mesher$(OBJ_EXT)
	rm -f polyhedral_surface_mesher$(EXE_EXT)
	rm -f implicit_surface_mesher$(OBJ_EXT)
	rm -f implicit_surface_mesher$(EXE_EXT)

run: implicit_surface_mesher$(EXE_EXT)
	./implicit_surface_mesher$(EXE_EXT)


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

