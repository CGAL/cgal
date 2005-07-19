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
	   -I../../../Data_structure_for_queries_3/include \
           -I../../../Mesh_2/include \
           -I../../../Mesh_3/include \
           -I../../../Triangulation_3/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(OTHER_CXXFLAGS)

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

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                surface_mesher$(EXE_EXT)
#		test$(EXE_EXT)	
#                surface_mesher_ellipsoid$(EXE_EXT) \
#		surface_mesher_with_hierarchy$(EXE_EXT)

surface_mesher$(EXE_EXT): surface_mesher$(OBJ_EXT) $(OBJECTS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)surface_mesher surface_mesher$(OBJ_EXT) $(LDFLAGS) $(OBJECTS)

clean: \
                   surface_mesher.clean \
		   implicit_function.clean

run: all
	./surface_mesher


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

