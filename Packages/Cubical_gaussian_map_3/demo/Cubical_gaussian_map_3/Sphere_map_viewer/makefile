# author: Kapelushnik Lior, 040818304
#
# part of agc project about
# spherical arrangements of none intersecting arcs of great circles on a sphere
#
# a viewer of spherical arrangements makefile

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
       -I../../../include \
       -I../../../../Arrangement_2/include \
        $(CGAL_CXXFLAGS) \
        $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          -lglut

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:		viewer$(EXE_EXT)

viewer$(EXE_EXT): viewer$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)viewer viewer$(OBJ_EXT) $(LDFLAGS)

clean: \
		   viewer.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

viewer.clean:
	rm -rf ./viewer ./viewer.o
