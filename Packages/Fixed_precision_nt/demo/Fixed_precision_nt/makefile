# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
	  $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
	  $(LONG_NAME_PROBLEM_LDFLAGS) \
	  $(CGAL_WINDOW_LDFLAGS)


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
