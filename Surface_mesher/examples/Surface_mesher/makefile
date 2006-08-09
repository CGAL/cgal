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

TARGETS = mesh_a_3d_gray_image$(EXE_EXT) \
          mesh_an_implicit_function$(EXE_EXT)

all: $(TARGETS)

mesh_a_3d_gray_image$(EXE_EXT): mesh_a_3d_gray_image$(OBJ_EXT) $(IMAGEIO_OBJECTS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_a_3d_gray_image mesh_a_3d_gray_image$(OBJ_EXT) $(IMAGEIO_OBJECTS) $(LDFLAGS)

mesh_an_implicit_function$(EXE_EXT): mesh_an_implicit_function$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_an_implicit_function mesh_an_implicit_function$(OBJ_EXT) $(LDFLAGS)

clean: \
                   mesh_a_3d_gray_image.clean \
                   mesh_an_implicit_function.clean \
		   imageio-clean



#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

