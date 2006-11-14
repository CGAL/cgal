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

CXXFLAGS = -I../../include/\
           -I../../../BGL/include\
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            test_edge_collapse_Polyhedron_3_LT \
                test_edge_collapse_Polyhedron_3_MP \

test_edge_collapse_Polyhedron_3_LT$(EXE_EXT): test_edge_collapse_Polyhedron_3_LT$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_LT test_edge_collapse_Polyhedron_3_LT$(OBJ_EXT) $(LDFLAGS)

test_edge_collapse_Polyhedron_3_MP$(EXE_EXT): test_edge_collapse_Polyhedron_3_MP$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_MP test_edge_collapse_Polyhedron_3_MP$(OBJ_EXT) $(LDFLAGS)


clean:          test_edge_collapse_Polyhedron_3_LT.clean \
                test_edge_collapse_Polyhedron_3_MP.clean \

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

