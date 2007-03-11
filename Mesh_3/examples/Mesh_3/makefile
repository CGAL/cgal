# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

DEFAULT_TARGETS = implicit_surfaces_mesher_3$(EXE_EXT)

all: $(DEFAULT_TARGETS)

.PHONY: default-rule

#---------------------------------------------------------------------#
#                    ImageIO files
#---------------------------------------------------------------------#

IMAGEIO_PATH = ../../../Surface_mesher/examples/Surface_mesher/ImageIO

include $(IMAGEIO_PATH)/imageio.mk

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

#OTHER_CXXFLAGS?= -D_GLIBCXX_CONCEPT_CHECKS 

CXXFLAGS = -I../../../Data_structure_for_queries_3/include \
           -I ../../../Surface_mesher/include \
           -I ../../../Mesh_2/include \
	   -I../../include \
           -I$(IMAGEIO_PATH) \
	   $(OTHER_INCLUDES) \
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

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

TOOLS=  cgal_to_medit \
        display_distribution \
        filter_remove_tets_from_medit \
        lanteri \
        lanteri_output_tet_mesh \
        off_to_ghs \
        read_mesh \
        slivers_exuder \
        stat_mesh \
	test_off

.PHONY: tools links

links:
	rm -f applications; ln -s ../../applications .
	for t in $(TOOLS); do rm -f $$t; ln -s applications/$$t .; done

tools:
	for t in $(TOOLS); do \
	  $(MAKE) -C ../../applications $$t; \
	done

implicit_surfaces_mesher_3$(OBJ_EXT): implicit_surfaces_mesher_3.cpp

implicit_surfaces_mesher_3$(EXE_EXT): implicit_surfaces_mesher_3$(OBJ_EXT) implicit_functions$(OBJ_EXT) parameters$(OBJ_EXT) $(IMAGEIO_OBJECTS) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surfaces_mesher_3 implicit_surfaces_mesher_3$(OBJ_EXT) implicit_functions$(OBJ_EXT) parameters$(OBJ_EXT) $(IMAGEIO_OBJECTS) $(LDFLAGS)

polyhedral_surface_mesher$(OBJ_EXT): polyhedral_surface_mesher.cpp 

polyhedral_surface_mesher$(EXE_EXT): polyhedral_surface_mesher$(OBJ_EXT) parameters$(OBJ_EXT) implicit_functions$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedral_surface_mesher polyhedral_surface_mesher$(OBJ_EXT) parameters$(OBJ_EXT)  implicit_functions$(OBJ_EXT) $(LDFLAGS)


clean: \
                   implicit_surfaces_mesher_3.clean \
                   polyhedral_surface_mesher.clean \
                   parameters.clean \
                   implicit_functions.clean \
		   imageio-clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
