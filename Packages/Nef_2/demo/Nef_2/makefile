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

CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:    \
	Nef_polyhedron_2-demo Simple-demo

Nef_polyhedron_2-demo$(EXE_EXT): Nef_polyhedron_2-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Nef_polyhedron_2-demo Nef_polyhedron_2-demo$(OBJ_EXT) $(LDFLAGS)

Simple-demo$(EXE_EXT): Simple-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Simple-demo Simple-demo$(OBJ_EXT) $(LDFLAGS)

clean:  \
	Nef_polyhedron_2-demo.clean Simple-demo.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

