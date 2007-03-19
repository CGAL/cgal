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

all:            test_edge_collapse_Polyhedron_3 \
                test_MP_strategy \
                test_LT_strategy \

test_edge_collapse_Polyhedron_3$(EXE_EXT): test_edge_collapse_Polyhedron_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3 test_edge_collapse_Polyhedron_3$(OBJ_EXT) $(LDFLAGS)

test_MP_strategy$(EXE_EXT): test_MP_strategy$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_MP_strategy test_MP_strategy$(OBJ_EXT) $(LDFLAGS)

test_LT_strategy$(EXE_EXT): test_LT_strategy$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_LT_strategy test_LT_strategy$(OBJ_EXT) $(LDFLAGS)


clean:          test_edge_collapse_Polyhedron_3.clean \
                test_MP_strategy.clean \
                test_LT_strategy.clean \
  
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

