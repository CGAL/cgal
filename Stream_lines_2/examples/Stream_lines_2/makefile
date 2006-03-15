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
                stl_regular_field$(EXE_EXT) \
                stl_triangular_field$(EXE_EXT) 

stl_regular_field$(EXE_EXT): stl_regular_field$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)stl_regular_field stl_regular_field$(OBJ_EXT) $(LDFLAGS)

stl_triangular_field$(EXE_EXT): stl_triangular_field$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)stl_triangular_field stl_triangular_field$(OBJ_EXT) $(LDFLAGS)

clean: \
                   stl_regular_field.clean \
                   stl_triangular_field.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

