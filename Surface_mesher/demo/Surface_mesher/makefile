# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

default-rule: all

.PHONY: default-rule

#---------------------------------------------------------------------#
#                    ImageIO files
#---------------------------------------------------------------------#

IMAGEIO_PATH = ../../examples/Surface_mesher/ImageIO/

include $(IMAGEIO_PATH)/imageio.mk

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
	   -I$(IMAGEIO_PATH) \
           -I../../../Data_structure_for_queries_3/include \
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

all:            \
                implicit_surface_mesher$(EXE_EXT) \
                polyhedron_remesher$(EXE_EXT)

implicit_surface_mesher$(EXE_EXT): implicit_surface_mesher$(OBJ_EXT) implicit_functions$(OBJ_EXT) parameters$(OBJ_EXT)  $(IMAGEIO_OBJECTS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surface_mesher implicit_surface_mesher$(OBJ_EXT) implicit_functions$(OBJ_EXT) parameters$(OBJ_EXT)  $(IMAGEIO_OBJECTS) $(LDFLAGS)

polyhedron_remesher$(EXE_EXT): polyhedron_remesher$(OBJ_EXT) implicit_functions$(OBJ_EXT) parameters$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_remesher implicit_functions$(OBJ_EXT) parameters$(OBJ_EXT) polyhedron_remesher$(OBJ_EXT) $(LDFLAGS)

clean: \
                   implicit_functions.clean \
                   implicit_surface_mesher.clean \
                   parameters.clean \
		   imageio-clean \
                   polyhedron_remesher.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

