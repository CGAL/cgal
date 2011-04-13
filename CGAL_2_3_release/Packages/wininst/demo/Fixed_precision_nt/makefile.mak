# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#		add your extra compilation flags here
#---------------------------------------------------------------------#
EXTRA_FLAGS=

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_WINDOW_LIBPATH)

LDFLAGS = $(CGAL_WINDOW_LDFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: delaunay$(EXE_EXT)

delaunay$(EXE_EXT): delaunay$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay delaunay$(OBJ_EXT) $(LDFLAGS)

clean: delaunay.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#


.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
