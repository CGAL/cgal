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
	   -I./ImageIO \
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

.PHONY: imageio

OBJECTS = \
	implicit_functions$(OBJ_EXT)

all:            \
                3d_image_surface_mesher$(EXE_EXT) \
                implicit_surface_mesher$(EXE_EXT)

imageio:
	$(MAKE) -C ImageIO

3d_image_surface_mesher$(EXE_EXT): 3d_image_surface_mesher$(OBJ_EXT) imageio
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)3d_image_surface_mesher 3d_image_surface_mesher$(OBJ_EXT) $(LDFLAGS) -L./ImageIO -lImageIO

implicit_surface_mesher$(EXE_EXT): implicit_surface_mesher$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surface_mesher implicit_surface_mesher$(OBJ_EXT) $(LDFLAGS)

clean: \
                   3d_image_surface_mesher.clean \
                   implicit_functions.clean \
                   implicit_surface_mesher.clean
	$(MAKE) -C ImageIO clean



#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

