# Created by the script cgal_create_makefile
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

IMAGEIO_PATH = ImageIO/

include $(IMAGEIO_PATH)/imageio.mk

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
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

implicit_OBJECTS = implicit_functions$(OBJ_EXT)

TARGETS = 3d_image_surface_mesher$(EXE_EXT) \
          implicit_surface_mesher$(EXE_EXT)

all: $(TARGETS)

3d_image_surface_mesher$(EXE_EXT): 3d_image_surface_mesher$(OBJ_EXT) $(IMAGEIO_OBJECTS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)3d_image_surface_mesher 3d_image_surface_mesher$(OBJ_EXT) $(IMAGEIO_OBJECTS) $(LDFLAGS)

implicit_surface_mesher$(EXE_EXT): $(implicit_OBJECTS) implicit_surface_mesher$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surface_mesher implicit_surface_mesher$(OBJ_EXT) $(LDFLAGS)

clean: \
                   3d_image_surface_mesher.clean \
                   implicit_functions.clean \
                   implicit_surface_mesher.clean \
		   imageio-clean



#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

