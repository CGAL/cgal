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
	chull_3d-windows-demo \
	chull_3d-geomview-demo \
	chull_dd-geomview-demo \
	chull_dd-runtime \
	delaunay_3d-windows-demo \
	delaunay_dd-runtime

chull_3d-windows-demo$(EXE_EXT): chull_3d-windows-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)chull_3d-windows-demo chull_3d-windows-demo$(OBJ_EXT) $(LDFLAGS)

chull_3d-geomview-demo$(EXE_EXT): chull_3d-geomview-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)chull_3d-geomview-demo chull_3d-geomview-demo$(OBJ_EXT) $(LDFLAGS)

chull_dd-geomview-demo$(EXE_EXT): chull_dd-geomview-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)chull_dd-geomview-demo chull_dd-geomview-demo$(OBJ_EXT) $(LDFLAGS)

chull_dd-runtime$(EXE_EXT): chull_dd-runtime$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)chull_dd-runtime chull_dd-runtime$(OBJ_EXT) $(LDFLAGS)

delaunay_3d-windows-demo$(EXE_EXT): delaunay_3d-windows-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay_3d-windows-demo delaunay_3d-windows-demo$(OBJ_EXT) $(LDFLAGS)

delaunay_dd-runtime$(EXE_EXT): delaunay_dd-runtime$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay_dd-runtime delaunay_dd-runtime$(OBJ_EXT) $(LDFLAGS)


clean:  \
	chull_3d-windows-demo.clean \
	chull_3d-geomview-demo.clean \
	chull_dd-geomview-demo.clean \
	chull_dd-runtime.clean \
	delaunay_3d-windows-demo.clean \
	delaunay_dd-runtime.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

